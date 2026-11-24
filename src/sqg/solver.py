"""Pseudo-spectral solver for active scalars  d_t th + u.grad th = 0,
u = grad-perp psi = (-d_y psi, d_x psi),  psi_hat = th_hat / |k|^alpha.

alpha = 1 : surface quasi-geostrophic (SQG), user's convention psi = (-Lap)^(-1/2) th
alpha = 2 : 2D Euler in vorticity form (th = omega, psi = (-Lap)^(-1) omega)

Domain [0, 2pi)^2, N^2 grid, real FFTs (scipy.fft, all cores), arrays indexed [iy, ix]
with x along the last (rfft) axis.  2/3 dealiasing: modes with |kx| >= N/3 or |ky| >= N/3
are zeroed in the state and in the nonlinear term.  Classical RK4 time stepping.
"""
from __future__ import annotations

import numpy as np
import numexpr as ne
import scipy.fft as sfft

WORKERS = -1


class Solver:
    def __init__(self, N: int, alpha: float = 1.0):
        self.N = N
        self.alpha = alpha
        self.dx = 2 * np.pi / N
        kx = np.fft.rfftfreq(N, 1.0 / N)            # 0..N/2
        ky = np.fft.fftfreq(N, 1.0 / N)             # 0..N/2-1, -N/2..-1
        self.KX = np.broadcast_to(kx[None, :], (N, N // 2 + 1))
        self.KY = np.broadcast_to(ky[:, None], (N, N // 2 + 1))
        self.kx1 = kx[None, :].astype(np.float64)
        self.ky1 = ky[:, None].astype(np.float64)
        kc = N / 3.0
        self.kcut = int(np.ceil(kc)) - 1           # largest retained |k_i|
        self.mask = ((np.abs(ky)[:, None] < kc) & (np.abs(kx)[None, :] < kc))
        K2 = kx[None, :] ** 2 + ky[:, None] ** 2
        with np.errstate(divide="ignore"):
            kinv = np.where(K2 > 0, K2 ** (-alpha / 2.0), 0.0)
        # fold the dealiasing mask into the inverse-Laplacian power
        self.kinv = (kinv * self.mask).astype(np.float64)
        self.maskf = self.mask.astype(np.float64)
        # spectrum bookkeeping: rfft half-plane, interior kx columns count twice
        w = np.full(N // 2 + 1, 2.0)
        w[0] = 1.0
        w[-1] = 1.0
        self.mult = np.broadcast_to(w[None, :], (N, N // 2 + 1))
        Kmag = np.sqrt(K2)
        self.shell = np.rint(Kmag).astype(np.int64).ravel()
        self.nshell = int(self.shell.max()) + 1
        with np.errstate(divide="ignore"):
            self.k_hinv = np.where(K2 > 0, 1.0 / Kmag, 0.0) ** 1  # |k|^-1 for H
        self.last_umax = None
        self.last_grad = None

    # ---------- transforms ----------
    def fwd(self, f):
        return sfft.rfft2(f, workers=WORKERS) / (self.N * self.N)

    def inv(self, fh):
        return sfft.irfft2(fh * (self.N * self.N), s=(self.N, self.N), workers=WORKERS,
                           overwrite_x=True)

    def to_spectral(self, f):
        return self.fwd(f) * self.maskf

    # ---------- right-hand side ----------
    def rhs(self, th, diag=False):
        kx, ky, kinv = self.kx1, self.ky1, self.kinv
        u = self.inv(ne.evaluate("-1j*ky*kinv*th"))
        v = self.inv(ne.evaluate("1j*kx*kinv*th"))
        tx = self.inv(ne.evaluate("1j*kx*th"))
        ty = self.inv(ne.evaluate("1j*ky*th"))
        if diag:
            self.last_umax = float(np.sqrt(ne.evaluate("max(u*u+v*v)")))
            g2 = ne.evaluate("tx*tx+ty*ty")
            self.last_grad = float(np.sqrt(g2.max()))
            flat = g2.ravel()
            cand = np.argpartition(flat, -16)[-16:]
            self.last_argmax = np.unravel_index(cand[np.argsort(flat[cand])[::-1]], g2.shape)
            del g2
        nl = ne.evaluate("u*tx+v*ty")
        del u, v, tx, ty
        nh = self.fwd(nl)
        m = self.maskf
        return ne.evaluate("-m*nh")

    def step(self, th, dt):
        """One classical RK4 step. Returns new state; computes umax/grad at the start."""
        k = self.rhs(th, diag=True)
        h2 = 0.5 * dt
        acc = k.copy()
        y = ne.evaluate("th + h2*k")
        k = self.rhs(y)
        ne.evaluate("acc + 2*k", out=acc)
        ne.evaluate("th + h2*k", out=y)
        k = self.rhs(y)
        ne.evaluate("acc + 2*k", out=acc)
        ne.evaluate("th + dt*k", out=y)
        k = self.rhs(y)
        h6 = dt / 6.0
        return ne.evaluate("th + h6*(acc + k)")

    # ---------- diagnostics ----------
    def spectrum(self, th):
        e = ne.evaluate("real(th*conj(th))*mult", local_dict={"th": th, "mult": self.mult})
        return 0.5 * np.bincount(self.shell, weights=e.ravel(), minlength=self.nshell)

    def invariants(self, th):
        e = ne.evaluate("real(th*conj(th))*mult", local_dict={"th": th, "mult": self.mult})
        L2 = 0.5 * float(e.sum())
        kin = self.k_hinv ** self.alpha
        H = 0.5 * float(ne.evaluate("sum(e*kin)"))
        return L2, H

    def grad_and_umax(self, th):
        self.rhs(th, diag=True)
        return self.last_grad, self.last_umax

    def _eval_grad(self, th, xs, ys):
        """|grad th| at the tensor points (ys[i], xs[j]) by direct Fourier summation over the
        retained modes (exact spectral interpolation)."""
        N = self.N
        kc = self.kcut
        rows = np.r_[0:kc + 1, N - kc:N]
        c = th[rows, : kc + 1]
        kx = np.arange(kc + 1, dtype=float)
        ky = self.ky1[rows, 0]
        m = np.full(kc + 1, 2.0)
        m[0] = 1.0
        ex = np.exp(1j * np.outer(kx, xs)) * m[:, None]          # (kx, j)
        ey = np.exp(1j * np.outer(ys, ky))                       # (i, ky)
        Ax = c @ (1j * kx[:, None] * ex)
        Ay = (c * (1j * ky[:, None])) @ ex
        fx = (ey @ Ax).real
        fy = (ey @ Ay).real
        return np.sqrt(fx * fx + fy * fy)

    def refined_gradmax(self, th, ncand=4):
        """Sub-grid maximum of |grad th| near the largest grid maxima: two levels of 9x9
        spectral-interpolation patches (spacing dx/4 over +-dx, then dx/16 over +-dx/4)."""
        iy, ix = self.last_argmax
        best = self.last_grad
        dx = self.dx
        used = []
        for a, b in zip(iy, ix):
            if any(abs(a - u) <= 2 and abs(b - v) <= 2 for u, v in used):
                continue
            used.append((a, b))
            x0, y0 = b * dx, a * dx
            for h in (dx / 4, dx / 16):
                off = h * np.arange(-4, 5)
                g = self._eval_grad(th, x0 + off, y0 + off)
                i, j = np.unravel_index(np.argmax(g), g.shape)
                x0, y0 = x0 + off[j], y0 + off[i]
                best = max(best, float(g[i, j]))
            if len(used) >= ncand:
                break
        return best

    def physical(self, th):
        return self.inv(th.copy())


def pad_state(th_c: np.ndarray, Nc: int, Nf: int) -> np.ndarray:
    """Zero-pad a (normalized, rfft) spectral state from Nc^2 to Nf^2 (exact interpolation)."""
    out = np.zeros((Nf, Nf // 2 + 1), dtype=np.complex128)
    h = Nc // 2
    out[:h, : h + 1] = th_c[:h, : h + 1]
    out[Nf - h:, : h + 1] = th_c[Nc - h:, : h + 1]
    # Nyquist rows/cols of the coarse grid are zero anyway after 2/3 dealiasing
    return out
