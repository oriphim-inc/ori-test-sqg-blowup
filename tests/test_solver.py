import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sqg.solver import Solver, pad_state  # noqa: E402
from sqg.ics import initial_spectral, grid, FAMILIES  # noqa: E402
from sqg.analysis import strip_fit, fit_growth_all, fit_delta_all  # noqa: E402


@pytest.mark.parametrize("alpha", [1.0, 2.0])
def test_single_shell_is_steady(alpha):
    # any combination of |k|=1 modes has psi proportional to theta -> u.grad theta = 0
    N = 64
    S = Solver(N, alpha)
    X, Y = grid(N)
    th = S.to_spectral(np.cos(X) + 0.7 * np.sin(Y) + 0.3 * np.cos(Y))
    th1 = th
    for _ in range(50):
        th1 = S.step(th1, 0.05)
    assert np.max(np.abs(th1 - th)) < 1e-13


def test_velocity_sign_convention():
    # theta = cos x -> psi = cos x (alpha=1) -> u = (-psi_y, psi_x) = (0, -sin x)
    N = 32
    S = Solver(N, 1.0)
    X, Y = grid(N)
    th = S.to_spectral(np.cos(X))
    v = S.inv(1j * S.kx1 * S.kinv * th)
    u = S.inv(-1j * S.ky1 * S.kinv * th)
    assert np.allclose(v, -np.sin(X), atol=1e-13)
    assert np.allclose(u, 0, atol=1e-13)


def test_rhs_matches_direct_evaluation():
    # cmt data: compute u.grad theta analytically for theta = sin x sin y + cos y (alpha=1)
    N = 64
    S = Solver(N, 1.0)
    X, Y = grid(N)
    th = S.to_spectral(np.sin(X) * np.sin(Y) + np.cos(Y))
    # psi = sin x sin y / sqrt2 + cos y
    px = np.cos(X) * np.sin(Y) / np.sqrt(2)
    py = np.sin(X) * np.cos(Y) / np.sqrt(2) - np.sin(Y)
    tx = np.cos(X) * np.sin(Y)
    ty = np.sin(X) * np.cos(Y) - np.sin(Y)
    adv = -py * tx + px * ty
    rhs = S.inv(S.rhs(th))
    assert np.allclose(rhs, -adv, atol=1e-12)


def test_rk4_fourth_order():
    N = 64
    S = Solver(N, 1.0)
    th0 = initial_spectral("A_cmt", S)
    T = 0.4

    def integ(n):
        th = th0
        for _ in range(n):
            th = S.step(th, T / n)
        return th

    ref = integ(320)
    e1 = np.abs(integ(10) - ref).max()
    e2 = np.abs(integ(20) - ref).max()
    order = np.log2(e1 / e2)
    assert 3.6 < order < 4.4, order


def test_invariants_conserved_short_run():
    N = 128
    for alpha in (1.0, 2.0):
        S = Solver(N, alpha)
        th = initial_spectral("D_merger", S)
        L0, H0 = S.invariants(th)
        for _ in range(100):
            th = S.step(th, 0.01)
        L1, H1 = S.invariants(th)
        assert abs(L1 / L0 - 1) < 1e-10 and abs(H1 / H0 - 1) < 1e-10


def test_pad_state_is_exact_interpolation():
    S1, S2 = Solver(64), Solver(128)
    th = initial_spectral("B_multisaddle", S1)
    f2 = S2.physical(pad_state(th, 64, 128))
    X, Y = grid(128)
    assert np.allclose(f2, FAMILIES["B_multisaddle"](X, Y), atol=1e-13)


def test_random_ic_is_grid_independent():
    S1, S2 = Solver(128), Solver(256)
    a = S1.physical(initial_spectral("E_random", S1))
    b = S2.physical(initial_spectral("E_random", S2))[::2, ::2]
    assert np.allclose(a, b, atol=1e-12)
    assert abs(np.sqrt(np.mean(a ** 2)) - 0.5) < 1e-3


def test_ics_are_resolved_at_512():
    for name in FAMILIES:
        S = Solver(512)
        E = S.spectrum(initial_spectral(name, S))
        assert E[S.kcut - 5:S.kcut + 1].max() < 1e-28 * E.max(), name


def test_strip_fit_recovers_delta():
    # synthetic shell spectrum E = k^-3 exp(-2 delta k) with mild noise
    rng = np.random.default_rng(0)
    k = np.arange(0, 400)
    for delta in (0.02, 0.05, 0.1):
        E = np.zeros(400)
        E[1:] = k[1:] ** -3.0 * np.exp(-2 * delta * k[1:]) * np.exp(0.05 * rng.standard_normal(399))
        r = strip_fit(E, 340)
        assert abs(r["delta"] / delta - 1) < 0.05, (delta, r)


def test_strip_fit_on_known_analytic_field():
    # f(x,y) = sum_k e^{-delta k} cos(kx): |f_hat| ~ e^{-delta k}, shell E ~ e^{-2 delta k}
    N = 256
    S = Solver(N)
    X, Y = grid(N)
    delta = 0.08
    r = np.exp(-delta)
    f = (1 - r * np.cos(X)) / (1 - 2 * r * np.cos(X) + r * r) - 1   # = sum_{k>=1} r^k cos kx
    out = strip_fit(S.spectrum(S.to_spectral(f)), S.kcut)
    assert abs(out["delta"] / delta - 1) < 0.02, out


def test_growth_fits_recover_parameters():
    t = np.linspace(0, 4, 200)
    G = 2.0 * (5.0 - t) ** -1.5
    fe, fd, ff = fit_growth_all(t, G)
    assert abs(ff["params"]["T"] - 5.0) < 1e-4 and abs(ff["params"]["gamma"] - 1.5) < 1e-4
    assert ff["aic"] < fe["aic"] and ff["aic"] < fd["aic"]
    G2 = 1.3 * np.exp(0.7 * np.exp(0.4 * t))
    fe, fd, ff = fit_growth_all(t, G2)
    assert abs(fd["params"]["c"] - 0.4) < 1e-4
    d = 0.3 * np.exp(-0.8 * t)
    fe, fc = fit_delta_all(t, d)
    assert abs(fe["params"]["b"] + 0.8) < 1e-8
    d2 = 0.4 * (6.0 - t) ** 1.0
    fe, fc = fit_delta_all(t, d2)
    assert abs(fc["params"]["T"] - 6.0) < 1e-4


def test_refined_gradmax_finds_offgrid_peak():
    # th = sin(x - x0) with x0 between grid points: true max|grad| = 1, grid max < 1
    N = 64
    S = Solver(N)
    X, Y = grid(N)
    x0 = 0.37 * S.dx
    th = S.to_spectral(0.5 * np.sin(X - x0) * np.cos(0.0 * Y) + 0.3 * np.cos(2 * (Y - 0.41 * S.dx)))
    g, _ = S.grad_and_umax(th)
    gr = S.refined_gradmax(th)
    # exact: |grad|^2 = 0.25 cos^2(x-x0) + 0.36 sin^2(2(y-y0)) -> max sqrt(0.61)
    assert gr >= g
    assert abs(gr - np.sqrt(0.61)) < 1e-4
