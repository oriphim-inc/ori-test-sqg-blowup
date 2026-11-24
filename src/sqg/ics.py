"""Initial-condition families (see PLAN.md). Every field is evaluated on an N^2 grid of
[0,2pi)^2, indexed [iy, ix]. All are trigonometric polynomials or entire periodic functions."""
from __future__ import annotations

import numpy as np


def grid(N):
    x = 2 * np.pi * np.arange(N) / N
    return np.meshgrid(x, x, indexing="xy")  # X[iy, ix] = x_ix, Y[iy, ix] = y_iy


def _pd(z, c):
    """Periodic, entire 'distance' 2 sin((z-c)/2): equals z-c near c."""
    return 2.0 * np.sin((z - c) / 2.0)


def _gauss(X, Y, xc, yc, a, b):
    return np.exp(-(_pd(X, xc) ** 2 / a ** 2 + _pd(Y, yc) ** 2 / b ** 2))


def cmt(X, Y):
    return np.sin(X) * np.sin(Y) + np.cos(Y)


def multisaddle(X, Y):
    return np.cos(2 * X) * np.cos(Y) + np.sin(X) * np.sin(Y) + np.cos(2 * X) * np.sin(3 * Y)


def ellipse(X, Y):
    return _gauss(X, Y, np.pi, np.pi, 0.8, 0.4)


def merger(X, Y):
    return _gauss(X, Y, np.pi - 0.7, np.pi, 0.5, 0.5) + _gauss(X, Y, np.pi + 0.7, np.pi, 0.5, 0.5)


def dipole(X, Y):
    return _gauss(X, Y, np.pi, np.pi - 0.7, 0.5, 0.5) - _gauss(X, Y, np.pi, np.pi + 0.7, 0.5, 0.5)


def random_smooth(X, Y, seed=1, k0=3.0, rms=0.5):
    """Random-phase field, |th_hat(k)| ~ k exp(-k^2/(2 k0^2)); defined on a fixed 64^2
    spectral set so that it is the same function at every N."""
    rng = np.random.default_rng(seed)
    M = 32
    th = np.zeros_like(X)
    ks = [(kx, ky) for kx in range(0, M + 1) for ky in range(-M, M + 1)
          if (kx > 0 or ky > 0)]
    amps = []
    for kx, ky in ks:
        k = np.hypot(kx, ky)
        amps.append(k * np.exp(-k * k / (2 * k0 * k0)))
    phases = rng.uniform(0, 2 * np.pi, len(ks))
    for (kx, ky), a, p in zip(ks, amps, phases):
        if a < 1e-18:
            continue
        th += a * np.cos(kx * X + ky * Y + p)
    # normalization constant computed analytically from amplitudes (grid-independent)
    norm = np.sqrt(0.5 * np.sum(np.array(amps) ** 2))
    return th * (rms / norm)


FAMILIES = {
    "A_cmt": cmt,
    "B_multisaddle": multisaddle,
    "C_ellipse": ellipse,
    "D_merger": merger,
    "E_random": random_smooth,
    "F_dipole": dipole,
}


def initial(name, N):
    X, Y = grid(N)
    return FAMILIES[name](X, Y)


def initial_spectral(name, solver):
    """Normalized, dealiased rfft coefficients of family `name` on solver's grid."""
    from .solver import pad_state
    N = solver.N
    if name == "E_random" and N > 128:
        f = initial(name, 128)
        c = np.fft.rfft2(f) / 128 ** 2
        return pad_state(c, 128, N) * solver.maskf
    return solver.to_spectral(initial(name, N))
