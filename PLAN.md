# PLAN — inviscid 2D SQG: is there a numerical sign of finite-time blow-up?

Requested by: held-out runner. Planned by Ori (2026-11-13).

## Question
For the inviscid surface quasi-geostrophic (SQG) equation on the periodic box [0,2π)²,

    ∂t θ + u·∇θ = 0,   u = ∇⊥ψ = (−∂yψ, ∂xψ),   ψ = (−Δ)^(−1/2) θ   (ψ̂ = θ̂/|k|),

does max|∇θ|(t) grow like an exponential, a double exponential, or a finite-time singularity
A (T*−t)^(−γ)? Does the analyticity-strip width δ(t) reach zero at finite time? Do these
trends hold up when the resolution is increased?

Sign convention: the user's u = ∇⊥(−Δ)^(−1/2)θ. CMT (1994) write θ = −(−Δ)^(1/2)ψ, which is
the opposite sign. θ ↦ −θ maps one convention onto the other. For the CMT data, −θ0 is θ0
shifted by π in y, so the two give the same max|∇θ|(t) and δ(t).

## What counts as the answer
- **Blow-up indicated:** for the CMT data at the top resolution, δ(t) extrapolates to 0 at a
  finite T* while the run is still resolved (δ ≥ 2Δx). The (T*−t)^(−γ) fit to max|∇θ| is
  clearly better (AIC) than the exponential and double-exponential fits, and gives a
  consistent T* (overlapping uncertainties). T* is stable as N doubles.
- **No blow-up indicated:** δ(t) decays at most exponentially (δ = δ0 e^(−βt) fits at least as
  well as δ = c(T*−t)^p, or the fitted T* moves later as N doubles). The gradient is fit at
  least as well by (double-)exponential growth.
- **Inconclusive:** anything in between, reported as such. The answer is a numerical
  indication over the resolved window, not a proof.

## Method (agreed)
- Pseudo-spectral Fourier method with real FFTs, double precision. The 2/3 rule (sharp
  cut-off |kx|,|ky| < N/3) is applied to θ̂ and to the nonlinear term.
- Classical RK4 with an adaptive CFL step dt = C·Δx/max|u| (C = 0.4, capped at dt_max = 0.01).
  C = 0.8 was rejected: RK4 damps the near-cut-off modes by up to ~50% per step there, which
  acts as numerical dissipation. C/2 runs of CMT at 512 and 2048 check the time-stepping
  error.
- No dissipation, hyperviscosity or filter beyond the 2/3 truncation.
- The same code with ψ̂ = θ̂/|k|² is 2D Euler (θ = vorticity). This is the sanity case: Euler
  is globally regular, so δ must not hit 0.

## Initial-condition families (six for SQG)
A. **CMT** (Constantin, Majda & Tabak 1994): θ0 = sin x sin y + cos y. A hyperbolic saddle
   with a strong front.
B. **Multi-mode saddle:** θ0 = cos 2x cos y + sin x sin y + cos 2x sin 3y. I believe this is
   the data of Constantin, Lai, Sharma, Tseng & Wu (2012), but that is not verified.
C. **Elliptical Gaussian vortex** (aspect ratio 2), made periodic and entire:
   θ0 = exp(−[X²/0.8² + Y²/0.4²]) with X = 2 sin((x−π)/2), Y = 2 sin((y−π)/2).
D. **Two like-signed Gaussian vortices** (merger), with the same periodic distance, width 0.5,
   centres (π∓0.7, π).
E. **Random smooth field**, seed 1: random phases with |θ̂(k)| ∝ k·exp(−k²/(2·3²)),
   normalized to rms 0.5.
F. **Opposite-signed Gaussian dipole** (±1, width 0.5, centres (π, π∓0.7)), which propagates
   and compresses fronts.
G. **Euler sanity case:** family A with ψ̂ = θ̂/|k|².

All families are entire functions or trigonometric polynomials, so δ(0) is large and the
round-off floor is reached well inside the grid at t=0.

## Diagnostics (every Δt_out = 0.02)
- max|∇θ| on the N² grid, and a sub-grid maximum refined by exact spectral interpolation
  near the largest grid maxima (used for fits and convergence); max|θ|.
- L2 = ½Σ|θ̂|² and Hamiltonian H = ½Σ_{k≠0}|θ̂|²/|k| (for Euler: enstrophy and energy
  ½Σ|θ̂|²/|k|²), as relative drift from t=0.
  Caveat known in advance: the 2/3-truncated Galerkin system conserves both quadratic
  invariants exactly, even when under-resolved. Conservation therefore tests time stepping
  and code correctness, not spatial resolution.
- Shell spectrum E(k). The analyticity-strip width δ(t) comes from a least-squares fit
  log E(k) = a − n log k − 2δk (Sulem, Sulem & Frisch 1983). The fit window runs over k
  between the last point above the round-off floor (or k_max = N/3) and a quarter of that.
  The uncertainty combines the regression standard error with the spread over three fit
  windows.
- **Distance to the grid tolerance:** δ/Δx (resolved while δ ≥ 2Δx), plus the fraction of
  the spectrum in the last 10% of shells below the cut-off. A run is stopped once δ < Δx,
  after 2 time units with δ < 2Δx, or at T_end = 20.

## Fits (on the resolved window of each run)
- Gradient G = max|∇θ|: (i) G = a e^(bt); (ii) G = a exp(b e^(ct)); (iii) G = A(T*−t)^(−γ).
  Nonlinear least squares on log G, with covariance-based standard errors, AIC, and
  sensitivity to the start of the window.
- Strip width: δ = δ0 e^(−βt) versus δ = c(T*−t)^p, with T* and its standard error.

## Resolutions and scale-down (agreed in this plan)
The target is N = 512, 1024, 2048, 4096, 8192.
- N = 512, 1024 and 2048 are run independently from t=0 for every family.
- 4096 and 8192 are continuation runs. Each starts from the next coarser run, zero-padded in
  Fourier space, at the last output time when that run's spectrum at its cut-off is below
  1e−20 of the peak shell in energy (1e−10 in amplitude). The two trajectories then agree to
  about 1e−10, which costs no accuracy but saves the early-time steps.
- 8192² is run for family A (CMT), and for the Euler case when time allows. 4096² is run for
  the other families when the time budget allows. If the budget forces a smaller set, the
  report says which.
- Time budget for the whole campaign: about 13–15 h of wall time (4096 runs capped at
  1.5 h each, 8192 at 8 h) on this machine (Apple M5 Pro,
  15 cores, 24 GB).

## Verdict rule (fixed before the main run)
This is implemented in src/sqg/report.py `verdict()` and in the Sceptic's checks. It uses the
main fit window (the last half of each run's resolved window), at the top N and the next N.
- **Blow-up indicated:** the δ collapse fit has a lower AIC than exponential decay; the
  finite-time gradient fit beats both (double-)exponential fits by more than 2 in AIC; T*
  from the gradient and from δ agree within 2 combined standard errors; and T*(gradient)
  agrees between the two N within 2 combined standard errors.
- **No blow-up indicated:** either exponential decay of δ fits at least as well, or T*(δ)
  moves later by more than 2 SE as N doubles; and the finite-time gradient fit is not
  clearly better (by 2 in AIC).
- **Inconclusive:** otherwise.
The Euler sanity case must come out as not "blow-up indicated".

## Literature context (qualitative; no number from it is used as a target)
CMT 1994 (Nonlinearity 7:1495) suggested a singularity near t≈8 from moderate resolution.
Ohkitani & Yamada 1997 (Phys. Fluids 9:876) suggested double-exponential growth instead.
Constantin, Nie & Schörghofer 1998 (Phys. Lett. A 241:168) found no blow-up. Córdoba 1998
(Ann. Math. 148:1135) ruled out simple hyperbolic-saddle collapse. Constantin et al. 2012
(J. Sci. Comput. 50:1) used higher resolution and saw no blow-up up to later times.
