> **Searching for blow-up in surface quasi-geostrophic flow.** A project Ori built and ran on a problem from Oriphim's held-out test set, which Ori was never tuned on. Short paper: [paper/paper.pdf](paper/paper.pdf). The Sceptic's review is in [notes/review.md](notes/review.md); limits are in [LIMITS.md](LIMITS.md). [oriphim.com](https://oriphim.com)
>
> How this repo was made: Ori built and ran the project and wrote the paper. Afterwards, Oriphim added the licence, edited the author line and wrote this note, with help from Claude Code.

# Inviscid 2D SQG: is there a numerical sign of finite-time blow-up?

By Oriphim (who posed the question and owns the project), with Ori.

## Answer
**No numerical sign of SQG blow-up within the resolved windows.** This covers the
Constantin–Majda–Tabak (CMT) data and five other families of initial data. The windows end
at t ≈ 7.9 for CMT, which is resolved at 8192² with max|∇θ| ≈ 46.5 there. Nothing is claimed
about later times, and none of this is a proof.

**Verdicts under the rule fixed in PLAN.md before the main run:**
- **"No blow-up indicated":** four families (A CMT, B multi-mode saddle, C elliptical vortex,
  F dipole).
- **"Inconclusive":** two families (D vortex merger, E random field). This is not a sign of
  blow-up; see "The δ test is weak" below.
- **"Blow-up indicated":** none. For every family, the gradient test never finds a
  finite-time law clearly better than exponential or double-exponential growth.

**CMT gradient.** max|∇θ| grows faster than exponentially. A double exponential has the best
AIC at every resolution, but its rate also shifts with the window, so no single simple law
fits the whole run. The (T*−t)^(−γ) fit puts T* later the further its window reaches:

| N | 512 | 1024 | 2048 | 4096 | 8192 |
|---|---|---|---|---|---|
| window ends at t | 5.2 | 5.9 | 6.6 | 7.3 | 7.9 |
| fitted T* | 6.4 | 7.3 | 8.1 | 10.5 | 13.7 |

A fixed singular time would show the opposite. max|∇θ| itself agrees across N to ≤ 5e-4, so
this drift reflects the window, not the resolution.

**CMT strip width.** δ(t) decays exponentially at every N. The decay rate β rises slowly
(0.82 → 0.99 from N = 512 to 8192), and the δ → 0 fit runs to its "no finite T*" bound.

**The δ test is weak.** D and E are "inconclusive" only because, on their 4096² runs, δ fits
better to a collapse at a finite T* than to exponential decay (D: T*≈6.6; E: T*≈2.2). The
same δ test applied to 2D Euler gives the same kind of collapse at every N (finite T*,
AIC margins 150–200), although 2D Euler is known to be regular. E's verdict is also fragile
in two ways:
- **It flips with the cut.** It reads "no blow-up" when 0.2 more time units after the
  hand-off are dropped.
- **Its top-grid δ fit can't be cross-checked.** That fit uses 13 points, and E's δ fails
  the agreement check between grids.

See LIMITS.md for the full list of caveats.

## Evidence
- **Gradient growth curves:** `results/figures/gradient_growth.png`. All families and
  resolutions; solid lines are resolved (δ ≥ 2Δx), dotted lines are under-resolved.
- **Strip width vs time:** `results/figures/strip_width.png`. δ(t) with bands, and 2Δx marked
  for each N.
- **Fits with uncertainties:**
  - `results/figures/fits_<family>.png` and `Tstar_vs_N.png`.
  - All parameters are in `results/summary.md` / `summary.json`. Each comes with a
    covariance standard error (SE) and a total error, which adds the spread over four
    fit-window starts. For the 4096/8192 runs the window is clipped at the hand-off, so the
    total equals the SE. Treat those ± as lower bounds. Totals in the hundreds mean "T* not
    determined".
- **Resolution convergence:** `results/figures/convergence.png`. max|∇θ| of N and 2N agree to
  ≤ 5e-4 (CMT; ≤ 6e-3 for the ellipse) while the coarser grid is resolved. The time at which
  they separate tracks δ ≈ 2Δx of the coarse grid.
- **Validation:**
  - L2 norm and Hamiltonian H = ½Σ|θ̂|²/|k| drift by ≤ 3e-12 while resolved and ≤ 8e-10 including the under-resolved tails (`conservation.png`). These
    are exact invariants of the truncated system, so this tests the code and time stepping,
    not spatial resolution.
  - max|θ| is conserved to ≤ 2.5e-4 while resolved (worst case: 512²).
  - Halving the time step changes max|∇θ| by ≤ 1e-11 (CMT at 512 and 2048).
  - 13 unit tests pass: steady states, the sign convention, RK4 order 4, exact interpolation,
    and fit recovery on synthetic data.
- **Checks:** the Sceptic's checks give 245 passed and 4 failed, all for the random family E:
  - three are the same issue: E's coarse-grid δ sits 17–35% below its fine-grid δ;
  - one is that E's verdict is not robust to the δ cut.

  Two of the checks were added or rewritten by the Sceptic during review, after results
  existed. LIMITS.md records which.
- **Post-hoc analysis choice (disclosed):** after the Sceptic's review, δ points of the
  4096/8192 continuation runs are dropped for 0.2 time units after the strip fit first
  reaches the cut-off. Before that, the high shells are still filling in after the
  hand-off. This choice was made after seeing results. It turned D and E from "no blow-up"
  into "inconclusive" and changed nothing else.

| family | top N | resolved to t | max\|∇θ\| there | best gradient law (AIC) | verdict |
|---|---|---|---|---|---|
| A CMT saddle | 8192 | 7.90 | 46.5 | double exp. | no blow-up indicated |
| B multi-mode saddle | 4096 | 2.26 | 18.8 | double exp. | no blow-up indicated |
| C elliptical vortex | 2048 | 20 (end) | 9.4 | double exp. (top N; finite-time at 1024) | no blow-up indicated |
| D vortex merger | 4096 | 4.64 | 8.9 | double exp. | inconclusive (δ collapse at top N only) |
| E random field | 4096 | 1.72 | 27.1 | double exp. | inconclusive, not robust (δ unreliable) |
| F dipole | 4096 | 11.9 | 4.3 | exponential (growth saturating) | no blow-up indicated |
| Euler, CMT data (control) | 4096 | 7.78 | 102 | (double) exp. | inconclusive (δ false positive) |

## Method in brief
- Pseudo-spectral scheme on [0,2π)²: real FFTs, 2/3 rule, classical RK4 with CFL 0.4,
  double precision, no dissipation. ψ̂ = θ̂/|k| and u = ∇⊥ψ.
- Grids 512, 1024 and 2048 run from t = 0. The 4096 and 8192 runs are zero-padded
  continuations from the coarser run, made while its spectrum was round-off-exact at the
  cut-off. This is the agreed scale-down in PLAN.md.
- δ(t) comes from fitting log E(k) = a − n log k − 2δk to the shell spectrum (Sulem, Sulem
  & Frisch 1983).
- Details, and the fixed verdict rule, are in PLAN.md.

## What this does not show
Read **LIMITS.md**. In short:
- The windows are finite, and CMT stops just short of t ≈ 8.
- The fit error bars are optimistic.
- δ depends on which shells are fitted, and it gives false positives on Euler.
- The initial-condition families are our choice.
- There is no 8192 run for Euler or the non-CMT families.

## Re-run
```
uv venv .venv && uv pip install --python .venv/bin/python numpy scipy matplotlib numexpr pytest
.venv/bin/python -m pytest tests/test_solver.py      # unit tests (~20 s)
.venv/bin/python run.py --nmax 512                   # quick version (~15 min)
.venv/bin/python run.py                              # full campaign (~14 h on an M5 Pro, 24 GB)
.venv/bin/python run.py --analyze-only               # refit and redraw from results/runs/*.npz
.venv/bin/python -m pytest tests/checks -q           # Sceptic's result checks
```
The full campaign needs about 8 GB of RAM for 8192². Finished runs in `results/runs/` are
reused; pass `--force` to recompute them. Logs are in `results/log_*.txt`.
