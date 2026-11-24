# LIMITS: what this campaign does NOT show

Written by the Sceptic at the review stage (2026-11-14), after the main campaign (commit 35d933f).
Updated at the second review (2026-11-14), after the analysis fix of commit 4fd362e (no
re-integration; only the fitting and reporting changed). Updated again at the third and final
review (2026-11-14), after the δ settling rule (DELTA_SETTLE = 0.2) was added to the analysis.

## Final verdicts (PLAN rule, fixed before the main run)
| Family | Verdict | What decides it |
|---|---|---|
| A CMT, B, C, F | no blow-up indicated | δ decays exponentially at least as well as a collapse, and the finite-time gradient law is not clearly best |
| D merger | inconclusive | δ collapse beats exponential decay at 4096 (T*δ ≈ 6.6), but the gradient criterion does not favour a finite-time law |
| E random | inconclusive, and not robust | as D at 4096 (T*δ ≈ 2.2), on a 13-point δ fit that cannot be checked against a resolved coarser grid |
| Euler control | inconclusive | δ collapse beats exponential decay at every N, even though 2D Euler is known to be regular |

No family reaches "blow-up indicated". "Inconclusive" for D and E does not point towards
blow-up: the Euler control shows the same δ behaviour (see below).

## What the answer is, and is not
- The result is a numerical indication over a finite resolved time window. It is not a proof
  either way. "No blow-up indicated" means that, up to the last resolved time of each family,
  the data do not look like a power-law singularity with a fixed T*. It does not mean that
  blow-up is ruled out at later times.
- Resolved windows (δ ≥ 2Δx at the finest grid): CMT up to t ≈ 7.9 (N = 8192), B up to t ≈ 2.3,
  C up to t = 20 (2048, never lost resolution), D up to t ≈ 4.6, E up to t ≈ 1.7,
  F up to t ≈ 11.9, Euler up to t ≈ 7.8 (N = 4096). Nothing is claimed beyond these times.
  Earlier studies proposed singular times for CMT near t ≈ 8 and later; t = 7.9 is just short
  of that region.

## Choices made by us, not by the user or the literature
- The six SQG initial-condition families were chosen by us. Other data (e.g. other saddle
  geometries, thinner initial fronts) are not tested.
- Family B is believed to be the data of Constantin, Lai, Sharma, Tseng & Wu (2012). This has
  NOT been checked against the paper. Do not cite B as a reproduction of that study.
- Scale-down and continuation: only 512, 1024 and 2048 run from t = 0. 4096 and 8192 are
  zero-padded continuations, started when the coarser run's spectrum at its cut-off fell below
  1e-20 of the peak. The hand-off itself is exact to round-off (checked). But see the
  δ transient below: it is only partly removed from the fits.
- Time budget (about 13–15 h on one laptop) set the grids: 8192 only for CMT; Euler only to
  4096 (no 8192); B, D, E, F to 4096; C only to 2048 (it stayed resolved to t = 20).
- Convergence in N is tested only up to these grids. "T* moves later with N" is not a
  resolution study in the usual sense (see below).

## How to read the fits
- max|∇θ| is converged across N to better than 5e-4 (relative) wherever the coarse run is
  resolved. So the per-N fits differ only because each N has a different fit window (the last
  half of its resolved window, or, for 4096 and 8192, the stretch from the hand-off to the end).
  The drift of T*(gradient) with N (CMT: 6.4, 7.3, 8.1, 10.5, 13.7) is a drift with the END
  TIME of the window, not a grid effect. It still means something: the converged gradient's
  local growth rate r = d log G/dt keeps rising but its rise slows down. A power law
  (T* − t)^(−γ) would make 1/r fall linearly to 0. Instead 1/r falls ever more slowly, and the
  local T* estimate, t + (1/r)/|d(1/r)/dt|, recedes from about 7.5 (t ≈ 4.6) to about 9.5
  (t ≈ 6.1) to about 16 (t ≈ 7.6). These are the Sceptic's rough estimates from 0.3-unit bins, not a fit.
- For 4096 and 8192 the "window-start sensitivity" is degenerate: the window is clipped at the
  hand-off, so the four starts give the same T*. Identical values there are not evidence of
  robustness.
- Each fitted parameter now has two errors. "SE" is the covariance standard error, which
  assumes independent residuals. The residuals are smooth and strongly correlated (outputs every
  0.02, model misfit, and step-like jumps in δ when the fit range kb changes). So SE alone is far
  too small: moving the window start moves T*(δ) by many SEs (Euler 2048: 15.4 ± 0.56 → 12.9 ±
  0.46 → 10.2 ± 0.28 → 8.2 ± 0.10). "Total" adds half the spread over the four window starts.
  This is only a partial fix:
  - For 4096 and 8192 the four windows are clipped at the hand-off and (nearly) coincide, so the
    spread is zero and "total" equals SE (e.g. CMT 8192 T*(grad) = 13.7 ± 0.08, Euler 4096
    T*(δ) = 8.8 ± 0.05). These ± are as optimistic as before.
  - Where some windows give a finite T* and others hit the "no finite T*" bound (≈ 1000), the
    spread is dominated by the bound and "total" is several hundred or thousand. That says
    "T* is not determined", not a real standard deviation.
  - The verdict rule (fixed in PLAN.md) still uses SE, not total. With total errors, criterion
    (3) "T*(grad) and T*(δ) agree" would hold more often (e.g. D 4096: 6.15 vs 6.59). It would
    not change any verdict, because the gradient criterion (2) fails for every family.
  Treat every ± as a lower bound on the uncertainty. The AIC margins are inflated for the same
  reason.
- "The finite-time gradient law is never clearly best" holds for the windows the verdict uses
  (each family's two finest grids) and for CMT at every N. It does NOT hold on every grid: on
  some coarser grids the finite-time law has the lowest AIC, by large margins (C at 1024, D at
  512 and 1024, E at 1024 for late window starts, F at 512; see the per-N fits in summary.md).
  For C, D and F the fitted T* falls inside the time range that a finer grid later resolved
  without blow-up, so they are not evidence of blow-up, but the blanket statement is wrong.
  E is the exception: at 1024 (window starts 0.6/0.7) T*(grad) ≈ 1.81–1.84, which lies just
  past the end of every finer E window (2048 to 1.42, 4096 to 1.72), so no finer grid has
  tested it. (Added at the draft review, 2026-11-23; E exception corrected at the re-review.)
- E_random's gradient T* does not recede: 2.33 at 2048 and 2.37 at 4096 (criterion (4) of the
  PLAN rule holds for E), less than 0.7 time units after E's resolved window ends at 1.72. E is
  kept from "blow-up indicated" by the gradient-AIC criterion (double exponential better by
  about 9) and by T*(grad) and T*(δ) disagreeing. "T* recedes as the window end moves later" is
  shown for CMT (and B, D), not for E.
- The double exponential has the lowest AIC for CMT at every N, but its parameters also drift
  with the window (c falls from 0.46 to 0.15). No single simple law describes the whole window.
  "Best AIC" ranks three curves; it does not identify a mechanism.

## The strip width δ is the weakest diagnostic
- δ comes from a least-squares fit to the shell spectrum over [kb/4, kb]. Its value depends on
  which shells are fitted. When kb differs between grids, δ differs by more than its stated
  error:
  - E_random: coarse δ is 17–35 % below fine δ even at δ ≈ 5Δx (the pre-registered check fails).
  - CMT: once the coarse run is below about 3.5Δx, coarse δ is 5–10 % above fine δ.
  - Euler: coarse δ is 10–14 % below fine δ at 4–5Δx, and up to about 37 % below at about 2Δx.
  The sign of the bias depends on the flow. So the end of every per-N δ window is biased, and
  the 2Δx "resolved" threshold is too loose for δ (it is fine for max|∇θ|).
- Continuation transient: for about 0.2–0.8 time units after each hand-off, the shells above
  the parent's cut-off are still filling in (kb < kcut). δ is then too high (CMT +8 to +18 %,
  B about +50 %, D up to +80 %, E about +115 %). And kb reaching kcut is not the end of the
  transient: for a few outputs afterwards δ is still well above the resolved parent (D 4096
  +72, +50, +30, +14 % at t = 3.46–3.52; B 4096 +18 %; E 4096 +118, +92, +50 %).
- Settling rule, chosen AFTER seeing results: the final analysis leaves out of the continuation
  δ fits every point with kb < kcut AND every point within 0.2 time units after kb first reaches
  kcut (DELTA_SETTLE = 0.2 in src/sqg/report.py). This value was picked after the second review
  had shown the verdicts changing with the cut. It is not part of the pre-registered plan. It
  moved D and E from "no blow-up indicated" to "inconclusive", i.e. towards the weaker claim,
  and left A, B, C, F unchanged.
- After the cut the D 4096 δ agrees with its resolved 2048 parent to within 1–2 % (t =
  3.66–4.06, Sceptic check). So D's δ collapse at 4096 is real behaviour of the diagnostic, not
  left over from the hand-off. It grows stronger, not weaker, with further cuts (ΔAIC 132 → 160,
  T*δ 6.6 → 6.1 for extra cuts of 0 → 0.2).
- E 4096: after the cut only 13 δ points remain (t = 1.48–1.72). The 2048 parent is resolved
  only to t = 1.42, so none of these points can be checked against a coarser grid, and E's
  coarse/fine δ check fails anyway (17–35 %). The δ collapse there is stable for extra cuts of
  0.06 and 0.12 (ΔAIC 64 → 57, T*δ 2.17 → 2.10). With an extra 0.2 too few points are left for
  any 4096 δ fit, the rule falls back to 2048 vs 1024, and the verdict flips to "no blow-up
  indicated" (T*δ moves later with N: 1.72 → 2.38). Check
  test_verdict_robust_to_delta_settling_cut[E_random] fails for this reason. E's verdict
  depends on which grid is the top one and is not robust. Report it as "inconclusive (not robust
  to the δ cut)".
- The check test_continuation_delta_transient_excluded originally asserted that the fit used
  exactly the kb ≥ kcut points. That failed by construction under the stricter rule. The
  Sceptic rewrote it at the third review to assert the original intent (no kb < kcut point in
  the fit) plus a match with an independent refit under the disclosed 0.2 settling. It passes
  for all continuation runs. This is a change to a check after results existed, and is recorded
  here for that reason.
- F 4096 shows a −19 to −32 % δ dip just after the cut that recovers by t ≈ 6.8. With the
  settling cut most of it is excluded; F's verdict does not depend on δ collapse in any case.
- The Euler control shows that the δ-collapse test cannot tell a regular flow from a singular
  one. Its δ-collapse fits beat exponential decay at 2048 (T* ≈ 12.9) and 4096 (T* ≈ 8.8), and even with
  δ ≥ 4Δx at 2048 they still prefer a collapse, with T* ≈ 76. Any slow rise in the δ decay rate
  is enough to trigger it, and regular flows show such a rise. The SQG verdicts therefore rest
  in practice on the gradient criterion. "δ decays exponentially" for SQG means only that SQG's
  δ looks no more singular than Euler's over the same window.
- The composite (mixed-N) δ fits should not be cited as evidence either way (summary.md now
  labels them "NOT citable"). With the kb < kcut points removed, the CMT composite still jumps
  with the window start (T* ≈ 29, 29, 35, then no finite T*), because each seam joins the end of
  a coarse window (biased δ) to the start of a fine one.
- The CMT δ decay rate in the main window rises with N: β = 0.82, 0.87, 0.90, 0.95, 0.99 for
  N = 512 … 8192. Since each N has a later window, this is the same slow steepening seen in the
  gradient, i.e. δ decays somewhat faster than one exponential over the whole run. Each per-N
  window still prefers exponential decay; with the earliest window start (0.4) at 1024 and 2048
  a δ collapse wins instead, with T* ≈ 44–46. Rising β is the kind of drift that produced the
  Euler false positive, so it is not evidence either way.

## Unmodelled / untested
- No dissipation, no forcing, a periodic box only, double precision only, and one dealiasing
  rule (2/3 sharp cut-off). Other truncations (e.g. Hou–Li smooth filter) were not compared.
- The time step: the C/2 runs (CMT at 512 and 2048) agree with C = 0.4 to about 1e-13 in
  max|∇θ| while resolved. RK4 still damps the near-cut-off modes a little at C = 0.4 (about
  1.6 % per step for the fastest-advected corner modes). This is checked to be harmless only for
  CMT, and only while resolved.
- Conservation of L2 and H (≤ 1e-11 here) tests the time stepping and the code, not the spatial
  resolution: the truncated system conserves both exactly even when under-resolved.
- One random seed for family E. Its δ diagnostic fails the coarse/fine check, and its top-grid
  δ fit spans only 0.24 time units.
- Under-resolved data (dotted in the figures) are shown but not used. Nothing should be read off
  them.
