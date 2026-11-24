"""Sceptic's pass/fail checks for the SQG blow-up campaign (written BEFORE the main run).

Run after `python run.py`:   .venv/bin/python -m pytest tests/checks -q
Everything skips until results/summary.json and results/runs/*.npz exist.

These checks never require a particular blow-up answer for SQG. They decide whether the
numbers that feed the answer can be trusted: no NaNs, invariants, L-infinity conservation,
initial values against analytic ones, coarse/fine agreement while the coarse run is
resolved, exact continuation hand-offs, the time-step check, the Euler sanity case (known
answer: globally regular) and the presence of every deliverable.

Environment overrides (only for smoke-testing these checks on a tiny run):
  SQG_RESULTS=<dir>       results directory (default <project>/results)
  SQG_REQUIRE_PLAN=0      skip the checks that demand the agreed grids 512..8192

The sampling-error tolerances on grid maxima assume N >= 512, i.e. dx = 0.012 well below
the smallest initial length scale of any family (~0.3: ellipse minor axis 0.4, Gaussian width
0.5, random field 1/k0 = 1/3). On tiny smoke grids (N = 32) they fail by design.

"Resolved" is defined independently of src/ here: a run is resolved at an output time if
delta >= 2 dx at every earlier output where delta was finite (once lost, stays lost), as in
PLAN.md. "Spectrally exact" means E(k_cut)/max E < 1e-20 (the plan's hand-off criterion).
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
RES = Path(os.environ.get("SQG_RESULTS", ROOT / "results"))
RUNS = RES / "runs"
REQUIRE_PLAN = os.environ.get("SQG_REQUIRE_PLAN", "1") != "0"

_have = (RES / "summary.json").exists() and any(RUNS.glob("*.npz")) if RUNS.exists() else False
pytestmark = pytest.mark.skipif(not _have, reason="results/ not produced yet (run `python run.py`)")

SQG_CASES = ["A_cmt", "B_multisaddle", "C_ellipse", "D_merger", "E_random", "F_dipole"]
EULER = "G_euler_cmt"
CASES = SQG_CASES + [EULER]
HANDOFF = 1e-20
MAIN_START = "0.5"


# ----------------------------------------------------------------------------- loading
def _case_of(name):
    return name.split("_N")[0]


def _N_of(name):
    return int(name.split("_N")[1].split("_")[0])


RUN_NAMES = sorted(p.stem for p in RUNS.glob("*.npz")) if RUNS.exists() else []
MAIN_NAMES = [n for n in RUN_NAMES if "_cfl" not in n]
PAIRS = []
for _c in CASES:
    _ns = sorted((n for n in MAIN_NAMES if _case_of(n) == _c), key=_N_of)
    PAIRS += [(a, b) for a, b in zip(_ns[:-1], _ns[1:])]


@pytest.fixture(scope="session")
def runs():
    out = {}
    for p in sorted(RUNS.glob("*.npz")):
        d = np.load(p, allow_pickle=False)
        meta = json.loads(str(d["meta"]))
        cols = [str(c) for c in d["cols"]]
        arr = np.atleast_2d(np.asarray(d["data"], float))
        data = {c: arr[:, i] for i, c in enumerate(cols)}
        out[p.stem] = dict(meta=meta, data=data, spectra=np.asarray(d["spectra"], float))
    return out


@pytest.fixture(scope="session")
def summary():
    with open(RES / "summary.json") as f:
        return json.load(f)


def resolved(r):
    dl = r["data"]["delta"]
    ok = np.ones(dl.shape, bool)
    bad = np.isfinite(dl) & (dl < 2.0 * r["meta"]["dx"])
    if bad.any():
        ok[np.argmax(bad):] = False
    return ok


def exact(r):
    return r["data"]["cut_ratio"] < HANDOFF


def common(ra, rb):
    ta, tb = np.round(ra["data"]["t"], 6), np.round(rb["data"]["t"], 6)
    _, ia, ib = np.intersect1d(ta, tb, return_indices=True)
    return ia, ib


def ref_thmax(runs, name):
    """max|theta| at t=0 for this case (continuation runs: from the finest base run)."""
    r = runs[name]
    if r["data"]["t"][0] == 0.0:
        return r["data"]["thmax"][0]
    base = [runs[n] for n in runs if _case_of(n) == _case_of(name) and "_cfl" not in n
            and runs[n]["data"]["t"][0] == 0.0]
    return max(base, key=lambda b: b["meta"]["N"])["data"]["thmax"][0]


# ============================================================================ per run
@pytest.mark.parametrize("name", RUN_NAMES)
def test_no_nonfinite_values(runs, name):
    """Every stored diagnostic except delta must be finite at every output, and the run must
    not have stopped on a non-finite state. A NaN/inf anywhere means the integration blew up
    numerically, which is a code failure, not a physical result."""
    r = runs[name]
    assert r["meta"]["stop_reason"] in ("tend", "delta<dx", "max_wall"), r["meta"]["stop_reason"]
    for c in ("t", "gradmax", "thmax", "L2", "H", "dL2", "dH", "umax", "cut_ratio", "tail"):
        assert np.all(np.isfinite(r["data"][c])), (name, c)
    assert np.all(np.isfinite(r["spectra"])), name


@pytest.mark.parametrize("name", RUN_NAMES)
def test_quadratic_invariants(runs, name):
    """The 2/3-truncated Galerkin system conserves L2 and H (SQG) / enstrophy and energy
    (Euler) exactly, so any drift is time-stepping error plus round-off.
    (a) while spectrally exact: |drift| <= 1e-9. Round-off over <=1e6 steps is ~sqrt(n)*1e-16
        ~ 1e-13, and the pilot saw <= 1e-10; 1e-9 leaves 10x margin.
    (b) over the whole resolved window: |drift| <= 1e-6. At C=0.8 explicit RK4 damps
        near-cutoff modes (|R| down to ~0.5 per step for corner modes), so drift measures the
        hidden numerical dissipation. 1e-6 of the total is still 1000x below the 1e-3 level at
        which coarse/fine gradient differences are judged; above it the scheme is not the
        'no dissipation' model the plan describes."""
    r = runs[name]
    d = r["data"]
    ex, rs = exact(r), resolved(r)
    for c in ("dL2", "dH"):
        if ex.any():
            assert np.max(np.abs(d[c][ex])) <= 1e-9, (name, c, np.max(np.abs(d[c][ex])))
        assert np.max(np.abs(d[c][rs])) <= 1e-6, (name, c, np.max(np.abs(d[c][rs])))


@pytest.mark.parametrize("name", RUN_NAMES)
def test_linf_conserved_while_resolved(runs, name):
    """max|theta| (max|omega| for Euler) is conserved by the PDE (transport by a
    divergence-free field). The stored value is a grid maximum, so it carries a sampling
    error ~ (dx/l)^2/8 per direction, and a truncation error ~ amplitude at the cutoff,
    ~exp(-delta*k_cut) with k_cut = N/3.
    (a) spectrally exact (delta >~ 7 dx): sampling <= 2/(8*49) = 5e-3 worst case -> tol 5e-3.
    (b) delta >= 4 dx: sampling <= 2/128 = 1.6e-2, truncation exp(-8.4) = 2e-4 -> tol 2e-2.
    (c) delta >= 2 dx (resolved): sampling <= 2/32, truncation exp(-4.2) = 1.5e-2 -> tol 8e-2.
    A larger change means the field is not being transported (sign/normalization error or
    Gibbs oscillations in a run labelled resolved)."""
    r = runs[name]
    d = r["data"]
    dev = np.abs(d["thmax"] / ref_thmax(runs, name) - 1.0)
    rs = resolved(r)
    dl = d["delta"]
    four = rs & ~(np.isfinite(dl) & (dl < 4 * r["meta"]["dx"]))
    for msk, tol, lab in ((exact(r), 5e-3, "exact"), (four, 2e-2, "delta>=4dx"), (rs, 8e-2, "resolved")):
        if msk.any():
            assert dev[msk].max() <= tol, (name, lab, dev[msk].max(), d["t"][msk][np.argmax(dev[msk])])


@pytest.mark.parametrize("name", RUN_NAMES)
def test_delta_well_defined_while_resolved(runs, name):
    """Inside the resolved window: delta may be NaN only before it first becomes finite
    (spectrum still at round-off over too few shells); after that it must be finite, strictly
    positive and carry a finite uncertainty. A negative or vanishing delta while the
    spectrum still decays at the cutoff is a fitting bug, not a singularity."""
    r = runs[name]
    d = r["data"]
    rs = resolved(r)
    dl = d["delta"][rs]
    fin = np.isfinite(dl)
    if not fin.any():
        pytest.fail(f"{name}: delta never finite in resolved window")
    first = np.argmax(fin)
    assert np.all(fin[first:]), (name, "NaN gap in delta", d["t"][rs][first:][~fin[first:]][:5])
    assert np.all(dl[first:] > 0), name
    assert np.all(np.isfinite(d["delta_err"][rs][first:])), name
    assert np.all(d["delta_err"][rs][first:] >= 0), name


@pytest.mark.parametrize("name", [n for n in RUN_NAMES if _case_of(n) in ("A_cmt", EULER)])
def test_initial_values_match_analytic(runs, name):
    """CMT data theta0 = sin x sin y + cos y is a trigonometric polynomial: on the same N^2
    grid, max|grad theta|, max|theta| and max|u| at t=0 are known in closed form
    (psi = sin x sin y / |k|^alpha + cos y with |k| = sqrt 2). Spectral derivatives are exact
    to round-off, so tol 1e-10 relative. Catches unit/normalization errors (FFT scaling,
    wavenumber 2pi factors, alpha) in the production runs themselves."""
    r = runs[name]
    if r["data"]["t"][0] != 0.0:
        pytest.skip("continuation run")
    N, alpha = r["meta"]["N"], r["meta"]["alpha"]
    x = 2 * np.pi * np.arange(N) / N
    X, Y = np.meshgrid(x, x, indexing="xy")
    s = 2.0 ** (alpha / 2.0)
    th = np.sin(X) * np.sin(Y) + np.cos(Y)
    g = np.sqrt((np.cos(X) * np.sin(Y)) ** 2 + (np.sin(X) * np.cos(Y) - np.sin(Y)) ** 2).max()
    u = -(np.sin(X) * np.cos(Y) / s - np.sin(Y))
    v = np.cos(X) * np.sin(Y) / s
    um = np.sqrt(u * u + v * v).max()
    d = r["data"]
    assert abs(d["gradmax"][0] / g - 1) < 1e-10, (d["gradmax"][0], g)
    assert abs(d["thmax"][0] / np.abs(th).max() - 1) < 1e-10
    assert abs(d["umax"][0] / um - 1) < 1e-10, (d["umax"][0], um)


@pytest.mark.parametrize("name", MAIN_NAMES)
def test_continuation_handoff_is_exact(runs, name):
    """A continuation run starts from the zero-padded parent state at the parent's last
    spectrally exact output. Padding is exact interpolation, so at the hand-off:
    - start time == parent handoff_t, and parent E(k_cut)/max E < 1e-20 there;
    - shell spectra agree for k <= parent k_cut to 1e-12 relative (same coefficients, only
      summation-order round-off) wherever E > 1e-24 max E;
    - L2 agrees to 1e-12;
    - child max|grad theta| >= parent's (the fine grid contains the coarse points) and exceeds
      it by at most the sampling error at delta >~ 7 dx_parent, 2/(8*49) ~ 5e-3."""
    r = runs[name]
    par = r["meta"].get("parent")
    if not par:
        pytest.skip("run starts at t=0")
    pname = par.replace("_handoff.npz", "")
    assert pname in runs, f"parent {pname} missing"
    p = runs[pname]
    assert p["meta"]["N"] < r["meta"]["N"]
    t0 = r["data"]["t"][0]
    assert p["meta"]["handoff_t"] is not None and abs(p["meta"]["handoff_t"] - t0) < 1e-9
    ip = int(np.argmin(np.abs(p["data"]["t"] - t0)))
    assert abs(p["data"]["t"][ip] - t0) < 1e-9
    assert p["data"]["cut_ratio"][ip] < HANDOFF
    Ep, Ec = p["spectra"][ip], r["spectra"][0][: len(p["spectra"][ip])]
    m = Ep > 1e-24 * Ep[1:].max()
    assert np.max(np.abs(Ec[m] / Ep[m] - 1)) < 1e-12
    assert abs(r["data"]["L2"][0] / p["data"]["L2"][ip] - 1) < 1e-12
    rg = r["data"]["gradmax"][0] / p["data"]["gradmax"][ip] - 1
    assert -1e-9 <= rg <= 5e-3, rg


# ============================================================================ coarse/fine
@pytest.mark.parametrize("pair", PAIRS, ids=[f"{a}-vs-N{_N_of(b)}" for a, b in PAIRS])
def test_spectra_agree_while_coarse_exact(runs, pair):
    """While the coarse run is spectrally exact, its truncation error is ~1e-10 of the peak
    amplitude, so a shell holding E >= 1e-4 max E differs by <~ 2e-10*sqrt(1e4) = 2e-8
    relative; time-stepping differences (dt halves with dx) are RK4 O(dt^4) and were 5e-10
    in the pilot. Tol 1e-6 on shells 1..10 leaves 50x margin. This is the decisive dynamical
    comparison: grid maxima carry sampling error, spectra do not."""
    a, b = (runs[n] for n in pair)
    ia, ib = common(a, b)
    sel = exact(a)[ia]
    if not sel.any():
        pytest.skip("no common spectrally exact times")
    worst = 0.0
    for i, j in zip(ia[sel], ib[sel]):
        Ea, Eb = a["spectra"][i][1:11], b["spectra"][j][1:11]
        m = Ea >= 1e-4 * a["spectra"][i][1:].max()
        worst = max(worst, float(np.max(np.abs(Ea[m] / Eb[m] - 1))))
    assert worst <= 1e-6, worst


@pytest.mark.parametrize("pair", PAIRS, ids=[f"{a}-vs-N{_N_of(b)}" for a, b in PAIRS])
def test_gradient_agrees_while_coarse_resolved(runs, pair):
    """max|grad theta| of N and 2N at common times, judged on the coarse run's state.
    Two-sided (sampling: fine grid finds a larger max between coarse points, ~2*(dx/delta)^2/8,
    plus truncation exp(-delta k_cut)):  exact 1e-2; delta>=4dx 2e-2; delta>=2dx 8e-2.
    One-sided (the fine grid contains every coarse point, so G_coarse <= G_fine up to the
    coarse truncation error only):  exact 1e-6; delta>=4dx 5e-3; delta>=2dx 5e-2.
    The one-sided bound is the sharp one: G_coarse noticeably above G_fine means the coarse
    run was not converged while labelled resolved, so its window cannot enter the fits."""
    a, b = (runs[n] for n in pair)
    ia, ib = common(a, b)
    if ia.size == 0:
        pytest.fail("no common output times")
    Ga, Gb = a["data"]["gradmax"][ia], b["data"]["gradmax"][ib]
    rel = np.abs(Ga - Gb) / Gb
    up = (Ga - Gb) / Gb
    rs = resolved(a)[ia]
    dl = a["data"]["delta"][ia]
    four = rs & ~(np.isfinite(dl) & (dl < 4 * a["meta"]["dx"]))
    t = a["data"]["t"][ia]
    for msk, two, one, lab in ((exact(a)[ia], 1e-2, 1e-6, "exact"), (four, 2e-2, 5e-3, ">=4dx"),
                               (rs, 8e-2, 5e-2, ">=2dx")):
        if msk.any():
            k2, k1 = np.argmax(rel[msk]), np.argmax(up[msk])
            assert rel[msk][k2] <= two, (lab, "two-sided", rel[msk][k2], t[msk][k2])
            assert up[msk][k1] <= one, (lab, "one-sided", up[msk][k1], t[msk][k1])


@pytest.mark.parametrize("pair", PAIRS, ids=[f"{a}-vs-N{_N_of(b)}" for a, b in PAIRS])
def test_delta_agrees_while_coarse_resolved(runs, pair):
    """delta is a property of the solution, not of the grid. Where the coarse run is
    resolved and both fits use >= 32 shells (kb >= 32), the two estimates must agree within
    max(3 x combined stated uncertainty, 10 % of delta). Failure means the stated delta
    uncertainty (which feeds the T* error bars) is too small, or the coarse delta is biased by
    the cutoff."""
    a, b = (runs[n] for n in pair)
    ia, ib = common(a, b)
    da, db = a["data"], b["data"]
    m = (resolved(a)[ia] & np.isfinite(da["delta"][ia]) & np.isfinite(db["delta"][ib])
         & (da["kb"][ia] >= 32) & (db["kb"][ib] >= 32))
    if not m.any():
        pytest.skip("no comparable times")
    dA, dB = da["delta"][ia][m], db["delta"][ib][m]
    tol = np.maximum(3 * np.hypot(da["delta_err"][ia][m], db["delta_err"][ib][m]), 0.1 * dB)
    bad = np.abs(dA - dB) > tol
    assert not bad.any(), list(zip(da["t"][ia][m][bad][:5], dA[bad][:5], dB[bad][:5]))


# ============================================================================ time step
def test_dt_check(runs, summary):
    """The C/2 run (same N, same grid, so no sampling difference) differs from the C run only
    by time-stepping error. Halving C cuts RK4 error 16x, so the difference ~ the error of the
    production run. Over the reference run's resolved window it must be <= 1e-4 relative in
    max|grad theta|: 10x below the 1e-3 convergence threshold used by the report, so time
    error cannot masquerade as resolution dependence. This also bounds the effect of RK4's
    damping of near-cutoff modes at C=0.8. Recomputed here from the npz files."""
    assert "dt_check" in summary
    cfl = [n for n in runs if "_cfl" in n]
    assert cfl, "no reduced-CFL run"
    for n in cfl:
        ref = runs[n.split("_cfl")[0]]
        r = runs[n]
        assert r["meta"]["cfl"] <= 0.5 * ref["meta"]["cfl"] + 1e-12
        ia, ib = common(r, ref)
        msk = resolved(ref)[ib]
        # Not vacuous: if dt_max=0.01 caps both runs they take identical steps and the check
        # measures nothing. Require >= 1.8x the steps over the reference's resolved window.
        j = np.nonzero(msk)[0][-1]
        ns_r, ns_ref = r["data"]["steps"][ia][j], ref["data"]["steps"][ib][j]
        assert ns_r >= 1.8 * ns_ref, (n, "dt_max cap makes the dt check vacuous", ns_r, ns_ref)
        rg = np.abs(r["data"]["gradmax"][ia] - ref["data"]["gradmax"][ib]) / ref["data"]["gradmax"][ib]
        assert rg[msk].max() <= 1e-4, (n, rg[msk].max())
        # the reduced-CFL run must reach as far as the reference's resolved window
        assert r["data"]["t"][ia][msk][-1] >= ref["data"]["t"][resolved(ref)][-1] - 1e-9


# ============================================================================ per case
@pytest.mark.parametrize("case", CASES)
def test_at_least_three_resolutions(runs, case):
    """The resolution question needs at least three grids per family (incl. Euler)."""
    Ns = {r["meta"]["N"] for n, r in runs.items() if _case_of(n) == case and "_cfl" not in n}
    assert len(Ns) >= 3, (case, sorted(Ns))


@pytest.mark.skipif(not REQUIRE_PLAN, reason="SQG_REQUIRE_PLAN=0")
@pytest.mark.parametrize("case", CASES)
def test_agreed_base_grids(runs, case):
    """PLAN.md: N = 512, 1024, 2048 run independently from t=0 for every family."""
    base = {r["meta"]["N"] for n, r in runs.items() if _case_of(n) == case and "_cfl" not in n
            and r["meta"].get("parent") is None and r["data"]["t"][0] == 0.0}
    assert {512, 1024, 2048} <= base, (case, sorted(base))


@pytest.mark.skipif(not REQUIRE_PLAN, reason="SQG_REQUIRE_PLAN=0")
def test_cmt_top_resolutions(runs):
    """PLAN.md: CMT is the headline case: 4096 must exist whenever 2048 lost resolution, and
    8192 whenever 4096 lost resolution (stop 'delta<dx'). A 4096 that stopped on max_wall is
    allowed but must be reported as such."""
    r2 = runs.get("A_cmt_N2048")
    assert r2 is not None
    if r2["meta"]["stop_reason"] == "delta<dx":
        assert "A_cmt_N4096" in runs
        if runs["A_cmt_N4096"]["meta"]["stop_reason"] == "delta<dx":
            assert "A_cmt_N8192" in runs


# ---------------------------------------------------------------- decision rule (fixed now)
def blowup_indicated(case_summary):
    """PLAN.md 'Blow-up indicated', operationalized BEFORE the run (main window start 0.5):
    at the top N with fits and the next lower N,
      (1) delta: AIC(c(T*-t)^p) < AIC(exp decay)            (delta extrapolates to 0),
      (2) gradient: AIC((T*-t)^-g) < min(AIC exp, AIC dexp) - 2   ('clearly better'),
      (3) T*(gradient) and T*(delta) agree within 2 combined standard errors,
      (4) T*(gradient) at top N and next lower N agree within 2 combined standard errors.
    Returns True only if all four hold."""
    pn = case_summary["fits"]["per_N"]
    Ns = sorted(int(n) for n in pn)
    if len(Ns) < 2:
        return False
    w = pn[str(Ns[-1])]["windows"][MAIN_START]
    w2 = pn[str(Ns[-2])]["windows"][MAIN_START]
    if not w["delta"]:
        return False
    g, d = w["gradient"], w["delta"]
    c1 = d["finite_time_collapse"]["aic"] < d["exponential_decay"]["aic"]
    c2 = g["finite_time"]["aic"] < min(g["exponential"]["aic"], g["double_exponential"]["aic"]) - 2
    Tg, eg = g["finite_time"]["params"]["T"], g["finite_time"]["errors"]["T"]
    Td, ed = d["finite_time_collapse"]["params"]["T"], d["finite_time_collapse"]["errors"]["T"]
    c3 = abs(Tg - Td) <= 2 * np.hypot(eg, ed)
    g2 = w2["gradient"]["finite_time"]
    c4 = abs(Tg - g2["params"]["T"]) <= 2 * np.hypot(eg, g2["errors"]["T"])
    return bool(c1 and c2 and c3 and c4)


def test_euler_sanity_not_flagged_as_blowup(summary):
    """2D Euler with smooth data is globally regular (Wolibner 1933, Kato 1967). If the same
    pipeline and decision rule call 'blow-up indicated' for Euler, the rule cannot
    discriminate and any SQG blow-up verdict from it is untrustworthy. (The pilot already
    showed short-window fits prefer collapse for Euler; the N-stability criterion must remove
    that.)"""
    assert EULER in summary["cases"]
    assert not blowup_indicated(summary["cases"][EULER])


@pytest.mark.parametrize("case", CASES)
def test_decision_rule_computable(summary, case):
    """Every family must have per-N fits at >= 2 resolutions with delta fits at the top one,
    so the agreed three-way verdict can be evaluated (no answer is required)."""
    pn = summary["cases"][case]["fits"]["per_N"]
    assert len(pn) >= 2, (case, list(pn))
    blowup_indicated(summary["cases"][case])


# ============================================================================ deliverables
@pytest.mark.parametrize("case", CASES)
def test_fits_present_with_uncertainties(summary, case):
    """Deliverable: gradient (exp, double-exp, finite-time) and strip-width (exp decay,
    collapse) fits with finite parameters and finite standard errors, for the composite and
    every per-N fit at the main window. The composite main window must hold >= 20 points
    (1 time unit at dt_out 0.05): fewer cannot discriminate 3-parameter laws. Fitted T* must
    lie after the data. Also: the composite resolved window cannot extend past the longest
    resolved window of any run (no under-resolved data in the fits)."""
    cs = summary["cases"][case]["fits"]
    sets = [("composite", cs["composite"])] + list(cs["per_N"].items())
    for lab, fs in sets:
        w = fs["windows"][MAIN_START]
        assert set(w["gradient"]) == {"exponential", "double_exponential", "finite_time"}, lab
        assert set(w["delta"]) == {"exponential_decay", "finite_time_collapse"}, (case, lab)
        for f in list(w["gradient"].values()) + list(w["delta"].values()):
            for k, v in f["params"].items():
                assert np.isfinite(v), (case, lab, f["model"], k)
                assert np.isfinite(f["errors"][k]), (case, lab, f["model"], k, "error")
            assert np.isfinite(f["aic"])
            if "T_minus_tlast" in f:
                assert f["T_minus_tlast"] > 0
    assert cs["composite"]["windows"][MAIN_START]["gradient"]["exponential"]["n"] >= 20, case


@pytest.mark.parametrize("case", CASES)
def test_summary_consistent_with_runs(runs, summary, case):
    """The summary's resolved windows must match an independent recomputation from the npz
    files (delta >= 2 dx, once lost stays lost), and the composite fit window may not end
    later than the latest resolved time of any run of the case."""
    names = [n for n in runs if _case_of(n) == case and "_cfl" not in n]
    tmax = 0.0
    for n in names:
        r = runs[n]
        mine = float(r["data"]["t"][resolved(r)][-1])
        tmax = max(tmax, mine)
        assert abs(summary["runs"][n]["t_end_resolved"] - mine) < 1e-9, n
    assert summary["cases"][case]["fits"]["composite"]["t_end_resolved"] <= tmax + 1e-9
    assert len(summary["cases"][case]["convergence"]) >= len(names) - 1


def test_figures_and_digest_present():
    """Deliverables: gradient curves, strip width, conservation, convergence, T* vs N, CMT
    spectra, per-case fit figures, and the human-readable digest."""
    fig = RES / "figures"
    need = ["gradient_growth.png", "strip_width.png", "conservation.png", "convergence.png",
            "Tstar_vs_N.png", "spectra_A_cmt.png"] + [f"fits_{c}.png" for c in CASES]
    missing = [f for f in need if not (fig / f).exists() or (fig / f).stat().st_size == 0]
    assert not missing, missing
    assert (RES / "summary.md").exists()


# ------------------------------------------------------------------ added at the review stage
# Found in review (2026-11-14): just after a zero-padded hand-off, the shells above the parent's
# cut-off are still filling in, kb < kcut, and the strip fit over [kb/4, kb] sees that cliff.
# delta in the continuation run is then 8-18 % ABOVE the coarse run's delta at the same time
# (CMT 2048->4096 and 4096->8192, ~0.5 time units), which is a diagnostic artefact, not dynamics.
# Those points must not enter the per-N delta fits that decide the verdict.
CONT_NAMES = [n for n in MAIN_NAMES if _case_of(n) in CASES]


@pytest.mark.parametrize("name", CONT_NAMES)
def test_continuation_delta_transient_excluded(runs, summary, name):
    r = runs[name]
    d = r["data"]
    if d["t"][0] == 0.0:
        pytest.skip("run starts at t=0")
    pn = summary["cases"][_case_of(name)]["fits"]["per_N"].get(str(r["meta"]["N"]))
    if not pn or not pn["windows"][MAIN_START]["delta"]:
        pytest.skip("no per-N delta fit for this run")
    w = pn["windows"][MAIN_START]
    kcut = r["meta"]["kcut"]
    m = (resolved(r) & np.isfinite(d["delta"]) & (d["t"] >= w["t0"] - 1e-9)
         & np.isfinite(d["kb"]) & (d["kb"] >= kcut))
    # Third review (2026-11-14): the analysis now also drops DELTA_SETTLE (disclosed, post hoc,
    # 0.2) time units after kb first reaches kcut. The original assertion (fit uses EXACTLY the
    # kb >= kcut points) failed by construction under that stricter rule. The intent is kept:
    # every fitted point must have kb >= kcut (subset), and the fit must equal an independent
    # refit on the settled points after the disclosed settling time.
    full = np.isfinite(d["kb"]) & (d["kb"] >= kcut)
    settle = 0.2
    m = m & (d["t"] >= d["t"][np.argmax(full)] + settle - 1e-9)
    if m.sum() < 6:
        pytest.skip("fewer than 6 settled delta points in the window")
    t, y = d["t"][m], np.log(d["delta"][m])
    b = np.polyfit(t, y, 1)[0]
    f = w["delta"]["exponential_decay"]
    assert f["n"] == int(m.sum()), (name, "delta fit uses", f["n"], "points; settled points:", int(m.sum()))
    assert abs(f["params"]["b"] - b) <= 1e-6 * abs(b), (name, f["params"]["b"], b)


# ------------------------------------------------------------------ added at the second review
# Found in the second review (2026-11-14): kb reaching kcut is NOT the end of the hand-off
# transient. For 0.02-0.1 time units afterwards the continuation delta is still far above the
# well-resolved parent (D_merger 4096: +72, +50, +30, +14 % at t = 3.46-3.52; B 4096: +18 %;
# E 4096: +118 %, +92 %, +50 %). A verdict that rests on the delta criterion must not change
# when a few more settling outputs are dropped. Uses the project's own fitting and verdict code
# (sqg.report._fit_set / verdict) on the stored runs, so it tests sensitivity, not correctness.
SETTLE = (0.06, 0.12, 0.2)


@pytest.mark.parametrize("case", SQG_CASES)
def test_verdict_robust_to_delta_settling_cut(runs, summary, case):
    import copy
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    import sqg.report as rp
    cs0 = summary["cases"][case]
    conts = [n for n in MAIN_NAMES if _case_of(n) == case and runs[n]["meta"].get("parent")
             and str(runs[n]["meta"]["N"]) in cs0["fits"]["per_N"]]
    if not conts:
        pytest.skip("no continuation runs with per-N fits")
    R = rp.load_runs(str(RUNS))
    v0 = cs0["verdict"]["verdict"]
    got = {}
    for extra in SETTLE:
        cs = copy.deepcopy(cs0)
        for n in conts:
            r = R[n]
            d, m = r["data"], rp.resolved_mask(r)
            df = d["delta_fit"].copy()
            fin = np.nonzero(np.isfinite(df))[0]
            if len(fin):
                df[d["t"] < d["t"][fin[0]] + extra - 1e-9] = np.nan
            fs = rp._fit_set(d["t"][m], d["G"][m], df[m])
            cs["fits"]["per_N"][str(r["meta"]["N"])] = json.loads(json.dumps(fs, default=float))
        got[extra] = rp.verdict(cs)["verdict"]
    assert all(v == v0 for v in got.values()), (case, "reported:", v0, "after extra settling cut:", got)
