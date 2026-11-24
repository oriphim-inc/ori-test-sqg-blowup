"""Analyticity-strip fits and growth-law fits."""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

FLOOR = 1e-26      # shells with E(k) < FLOOR * max E are treated as round-off
MIN_PTS = 8


def strip_fit(E: np.ndarray, kcut: int):
    """Fit log E(k) = a - n log k - 2 delta k on shell spectrum E (index = |k|).

    Fit range: k in [kb*f0, kb*f1] where kb = min(kcut, last shell above FLOOR*max E).
    Main window (1/4, 1); alternatives (1/2, 1) and (1/4, 3/4).
    Returns dict(delta, delta_se, delta_err, n, kb) ; delta_err combines regression SE and
    window spread. NaN if too few points (spectrum still at round-off over most k).
    """
    E = np.asarray(E, float)
    kmax = min(kcut, len(E) - 1)
    Emax = E[1:].max()
    above = np.nonzero(E[1:kmax + 1] > FLOOR * Emax)[0]
    out = dict(delta=np.nan, delta_se=np.nan, delta_err=np.nan, n=np.nan, kb=np.nan)
    if above.size == 0:
        return out
    kb = int(above[-1] + 1)
    out["kb"] = kb
    res = []
    for f0, f1 in ((0.25, 1.0), (0.5, 1.0), (0.25, 0.75)):
        ka, kz = max(2, int(np.floor(kb * f0))), int(np.floor(kb * f1))
        k = np.arange(ka, kz + 1)
        k = k[E[k] > 0]
        if k.size < MIN_PTS:
            res.append(None)
            continue
        A = np.column_stack([np.ones(k.size), -np.log(k), -2.0 * k])
        y = np.log(E[k])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A @ coef
        dof = max(k.size - 3, 1)
        s2 = r @ r / dof
        cov = s2 * np.linalg.inv(A.T @ A)
        res.append((coef[2], np.sqrt(cov[2, 2]), coef[1]))
    if res[0] is None:
        return out
    d0, se0, n0 = res[0]
    alts = [r[0] for r in res if r is not None]
    spread = (max(alts) - min(alts)) if len(alts) > 1 else 0.0
    out.update(delta=d0, delta_se=se0, delta_err=float(np.hypot(se0, spread / 2)), n=n0)
    return out


# ------------------------------------------------------------------ growth-law fits
def _aic(rss, n, p):
    return n * np.log(max(rss, 1e-300) / n) + 2 * p


def _cov(res):
    J = res.jac
    n, p = J.shape
    rss = float(res.fun @ res.fun)
    try:
        cov = np.linalg.inv(J.T @ J) * rss / max(n - p, 1)
    except np.linalg.LinAlgError:
        cov = np.full((p, p), np.nan)
    return cov, rss


def fit_exponential(t, y):
    """log y = la + b t."""
    ly = np.log(y)
    A = np.column_stack([np.ones_like(t), t])
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    r = ly - A @ coef
    rss = float(r @ r)
    cov = np.linalg.inv(A.T @ A) * rss / max(len(t) - 2, 1)
    return dict(model="exponential", params={"a": float(np.exp(coef[0])), "b": float(coef[1])},
                errors={"a": float(np.exp(coef[0]) * np.sqrt(cov[0, 0])), "b": float(np.sqrt(cov[1, 1]))},
                rss=rss, aic=_aic(rss, len(t), 2), n=len(t))


def fit_double_exponential(t, y):
    """log y = la + b exp(c t)."""
    ly = np.log(y)
    best = None
    for c0 in (0.05, 0.1, 0.2, 0.4, 0.8, 1.5):
        e = np.exp(c0 * t)
        A = np.column_stack([np.ones_like(t), e])
        (la0, b0), *_ = np.linalg.lstsq(A, ly, rcond=None)
        f = lambda p: p[0] + p[1] * np.exp(p[2] * t) - ly
        try:
            r = least_squares(f, [la0, b0, c0], method="lm", max_nfev=20000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    cov, rss = _cov(best)
    la, b, c = best.x
    se = np.sqrt(np.abs(np.diag(cov)))
    return dict(model="double_exponential",
                params={"a": float(np.exp(la)), "b": float(b), "c": float(c)},
                errors={"a": float(np.exp(la) * se[0]), "b": float(se[1]), "c": float(se[2])},
                rss=rss, aic=_aic(rss, len(t), 3), n=len(t))


def fit_finite_time(t, y, increasing=True):
    """increasing:  y = B ((T - t)/(T - t0))^(-gamma)   (blow-up of y at T)
       otherwise :  y = B ((T - t)/(T - t0))^(+p)       (collapse of y to 0 at T)
    t0 = first time in the window, so B ~ y(t0) stays finite even when T -> infinity.
    T = t_last + exp(s), with T - t_last bounded to [1e-6, 1e3]; hitting the upper bound
    means 'no finite T*' (the law degenerates to an exponential)."""
    ly = np.log(y)
    tl, t0 = t[-1], t[0]
    sgn = -1.0 if increasing else 1.0

    def model(p):
        T = tl + np.exp(p[2])
        return p[0] + sgn * p[1] * np.log((T - t) / (T - t0))

    best = None
    for dT in (0.01, 0.1, 0.5, 1, 2, 5, 10, 30, 100):
        T0 = tl + dT
        X = np.column_stack([np.ones_like(t), sgn * np.log((T0 - t) / (T0 - t0))])
        (lB0, g0), *_ = np.linalg.lstsq(X, ly, rcond=None)
        try:
            r = least_squares(lambda p: model(p) - ly, [lB0, g0, np.log(dT)], method="trf",
                              max_nfev=20000, x_scale="jac",
                              bounds=([-np.inf, -np.inf, np.log(1e-6)], [np.inf, np.inf, np.log(1e3)]))
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    cov, rss = _cov(best)
    lB, g, s = best.x
    se = np.sqrt(np.abs(np.diag(cov)))
    T = tl + np.exp(s)
    name = "gamma" if increasing else "p"
    return dict(model="finite_time" if increasing else "finite_time_collapse",
                params={"B": float(np.exp(lB)), name: float(g), "T": float(T), "t0": float(t0)},
                errors={"B": float(np.exp(lB) * se[0]), name: float(se[1]),
                        "T": float(np.exp(s) * se[2]), "t0": 0.0},
                rss=rss, aic=_aic(rss, len(t), 3), n=len(t), T_minus_tlast=float(np.exp(s)),
                T_at_upper_bound=bool(np.exp(s) > 0.99e3))


def fit_growth_all(t, G):
    return [fit_exponential(t, G), fit_double_exponential(t, G), fit_finite_time(t, G, True)]


def fit_delta_all(t, d):
    """delta = d0 exp(-beta t)  versus  delta = c (T - t)^p."""
    e = fit_exponential(t, d)
    e["model"] = "exponential_decay"
    return [e, fit_finite_time(t, d, increasing=False)]
