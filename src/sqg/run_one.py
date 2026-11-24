"""Integrate one run and store its time series (results/runs/<name>.npz)."""
from __future__ import annotations

import json
import os
import time

import numpy as np

from .analysis import strip_fit
from .ics import initial_spectral
from .solver import Solver, pad_state

HANDOFF_RATIO = 1e-20   # E(kcut)/max E below this -> state usable for zero-padded continuation


def run(name, family, N, alpha=1.0, tend=20.0, dt_out=0.02, cfl=0.4, dt_max=0.01,
        start=None, outdir="results/runs", statedir="results/states", log=print,
        stop_dx=1.0, max_wall=None, max_unresolved_time=2.0):
    """start: None (from t=0) or path to a handoff npz of a coarser run."""
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(statedir, exist_ok=True)
    S = Solver(N, alpha)
    if start is None:
        th = initial_spectral(family, S)
        t = 0.0
        parent = None
    else:
        d = np.load(start)
        th = pad_state(d["th"], int(d["N"]), N) * S.maskf
        t = float(d["t"])
        parent = os.path.basename(start)
    L20, H0 = S.invariants(th)
    rows, spectra = [], []
    nsteps = 0
    wall0 = time.perf_counter()
    handoff_t = None
    handoff = None
    t_lost = None
    reason = "tend"
    t_next = t
    while True:
        # ---- diagnostics at output time
        if t >= t_next - 1e-12:
            g, umax = S.grad_and_umax(th)
            gref = S.refined_gradmax(th)
            L2, H = S.invariants(th)
            E = S.spectrum(th)
            sf = strip_fit(E, S.kcut)
            Emax = E[1:].max()
            kc = S.kcut
            tail = E[int(0.9 * kc):kc + 1].sum() / E.sum()
            cut_ratio = E[kc] / Emax
            thmax = float(np.abs(S.physical(th)).max())
            rows.append([t, g, thmax, L2, H, (L2 - L20) / L20, (H - H0) / H0, sf["delta"],
                         sf["delta_err"], sf["delta_se"], sf["n"], sf["kb"], tail, cut_ratio,
                         umax, nsteps, gref])
            spectra.append(E[: kc + 1].astype(np.float64))
            if cut_ratio < HANDOFF_RATIO:
                handoff = th.copy()
                handoff_t = t
            dd = sf["delta"]
            log(f"[{name}] t={t:7.3f} |grad|={g:10.4g} ref={gref:10.4g} delta/dx={dd / S.dx:9.3g} "
                f"dL2={rows[-1][5]:+.2e} dH={rows[-1][6]:+.2e} steps={nsteps} "
                f"wall={time.perf_counter() - wall0:7.1f}s")
            t_next = t + dt_out
            if not np.isfinite(g) or not np.isfinite(L2):
                reason = "nonfinite"
                break
            if np.isfinite(dd) and dd < 2.0 * S.dx:
                if t_lost is None:
                    t_lost = t
            if np.isfinite(dd) and dd < stop_dx * S.dx:
                reason = "delta<dx"
                break
            if t_lost is not None and t - t_lost >= max_unresolved_time:
                reason = "delta<dx"   # under-resolved (delta<2dx) for max_unresolved_time
                break
            if t >= tend - 1e-12:
                reason = "tend"
                break
            if max_wall and time.perf_counter() - wall0 > max_wall:
                reason = "max_wall"
                break
        # ---- time step (umax from the last diagnostic or the previous step)
        umax = S.last_umax if S.last_umax else 1.0
        dt = min(cfl * S.dx / umax, dt_max, t_next - t)
        th = S.step(th, dt)
        t += dt
        if abs(t - t_next) < 1e-12:
            t = t_next
        nsteps += 1
    cols = ["t", "gradmax", "thmax", "L2", "H", "dL2", "dH", "delta", "delta_err", "delta_se",
            "n_exp", "kb", "tail", "cut_ratio", "umax", "steps", "gradmax_ref"]
    meta = dict(name=name, family=family, N=N, alpha=alpha, cfl=cfl, dt_out=dt_out,
                dt_max=dt_max, parent=parent, stop_reason=reason, handoff_t=handoff_t,
                wall_s=time.perf_counter() - wall0, steps=nsteps, dx=S.dx, kcut=S.kcut)
    if handoff is not None:
        np.savez(os.path.join(statedir, f"{name}_handoff.npz"), th=handoff, t=handoff_t, N=N)
    np.savez(os.path.join(statedir, f"{name}_final.npz"), th=th, t=t, N=N)
    np.savez(os.path.join(outdir, f"{name}.npz"), data=np.array(rows), cols=np.array(cols),
             spectra=np.array(spectra), meta=json.dumps(meta))
    return meta
