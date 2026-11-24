"""Analysis of finished runs: fits, convergence, figures, results/summary.json."""
from __future__ import annotations

import glob
import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .analysis import fit_growth_all, fit_delta_all  # noqa: E402
from .campaign import SQG, EULER  # noqa: E402

RESOLVED_DX = 2.0          # a time is "resolved" while delta >= 2 dx
WINDOW_STARTS = (0.4, 0.5, 0.6, 0.7)   # fit window = [s * t_end_resolved, t_end_resolved]
MAIN_START = 0.5
CONV_TOL = 1e-3
DELTA_SETTLE = 0.2         # extra settling time for delta after a hand-off


def load_runs(outdir="results/runs"):
    runs = {}
    for p in sorted(glob.glob(os.path.join(outdir, "*.npz"))):
        d = np.load(p, allow_pickle=False)
        meta = json.loads(str(d["meta"]))
        cols = [str(c) for c in d["cols"]]
        data = {c: d["data"][:, i] for i, c in enumerate(cols)}
        if "gradmax_ref" not in data:
            data["gradmax_ref"] = data["gradmax"]
        # G: sub-grid (spectrally interpolated) max|grad theta|, used for fits and convergence
        data["G"] = data["gradmax_ref"]
        # delta used in fits: for continuation runs drop the hand-off transient, i.e. points
        # whose strip-fit range has not yet reached the cut-off (high shells still filling)
        dfit = data["delta"].copy()
        if meta.get("parent"):
            full = np.isfinite(data["kb"]) & (data["kb"] >= meta["kcut"])
            dfit[~full] = np.nan
            # plus a settling time after kb first reaches kcut (post-hoc choice after the
            # Sceptic's review; D verdict is stable for 0.06-0.2, E_random is not)
            if full.any():
                t_full = data["t"][np.argmax(full)]
                dfit[data["t"] < t_full + DELTA_SETTLE - 1e-9] = np.nan
        data["delta_fit"] = dfit
        runs[meta["name"]] = dict(meta=meta, data=data, spectra=d["spectra"])
    return runs


def resolved_mask(r):
    d = r["data"]
    dl = d["delta"]
    ok = np.ones_like(dl, bool)
    bad = np.isfinite(dl) & (dl < RESOLVED_DX * r["meta"]["dx"])
    if bad.any():
        ok[np.argmax(bad):] = False   # once lost, stays lost
    return ok


def case_runs(runs, case):
    rs = [r for r in runs.values() if r["meta"]["name"].startswith(case + "_N")
          and r["meta"]["cfl"] == 0.4]
    return sorted(rs, key=lambda r: r["meta"]["N"])


def composite(rs):
    """Highest-resolution resolved value at each output time."""
    best = {}
    for r in rs:                       # ascending N: later overwrite earlier
        m = resolved_mask(r)
        d = r["data"]
        for i in np.nonzero(m)[0]:
            best[round(float(d["t"][i]), 6)] = (d["G"][i], d["delta_fit"][i], d["delta_err"][i],
                                                r["meta"]["N"])
    ts = np.array(sorted(best))
    arr = np.array([best[t] for t in ts]) if len(ts) else np.zeros((0, 4))
    return ts, arr


def _fit_set(t, G, dl):
    out = {"t_end_resolved": float(t[-1]), "windows": {}}
    for s in WINDOW_STARTS:
        sel = t >= s * t[-1]
        g = {f["model"]: f for f in fit_growth_all(t[sel], G[sel])}
        sd = sel & np.isfinite(dl)
        dfit = {f["model"]: f for f in fit_delta_all(t[sd], dl[sd])} if sd.sum() >= 6 else {}
        out["windows"][str(s)] = dict(t0=float(t[sel][0]), gradient=g, delta=dfit)
    # total error = hypot(covariance SE, half the spread over window starts); covariance SEs
    # assume independent residuals and understate the uncertainty
    ws = list(out["windows"].values())
    out["windows_clipped"] = len({round(w["t0"], 6) for w in ws}) < len(ws)
    for kind in ("gradient", "delta"):
        mdls = {m: w[kind][m]["params"] for w in ws if w[kind] for m in w[kind]}
        for mdl, pars in mdls.items():
            for par in pars:
                vals = [w[kind][mdl]["params"][par] for w in ws if w[kind] and mdl in w[kind]]
                spread = float(np.max(vals) - np.min(vals)) if vals else 0.0
                for w in ws:
                    if w[kind] and mdl in w[kind]:
                        f = w[kind][mdl]
                        f.setdefault("errors_total", {})[par] = float(np.hypot(f["errors"][par], spread / 2))
    return out


def fits_for_case(rs):
    res = {"per_N": {}}
    for r in rs:
        m = resolved_mask(r)
        d = r["data"]
        if m.sum() < 10:
            continue
        res["per_N"][str(r["meta"]["N"])] = _fit_set(d["t"][m], d["G"][m], d["delta_fit"][m])
    t, arr = composite(rs)
    res["composite"] = _fit_set(t, arr[:, 0], arr[:, 1])
    return res


def verdict(case_summary, start=MAIN_START):
    """Three-way verdict of PLAN.md, fixed before the main run (same rule as the Sceptic's
    tests/checks blowup_indicated). Uses the top N with fits and the next lower N.
    blow-up indicated: (1) delta collapse fit beats exp decay (AIC), (2) finite-time gradient
      fit beats exp and double-exp by > 2 in AIC, (3) T*(grad) and T*(delta) agree within 2
      combined SE, (4) T*(grad) at top N and next N agree within 2 combined SE.
    no blow-up indicated: [(a) exp decay of delta fits at least as well (AIC) OR (b) T*(delta)
      at top N is later than at the next N by > 2 combined SE] AND (c) the finite-time
      gradient fit is not clearly better (AIC_ft >= min(AIC_exp, AIC_dexp) - 2).
    otherwise inconclusive."""
    pn = case_summary["fits"]["per_N"]
    Ns = sorted((int(n) for n in pn if pn[n]["windows"][str(start)]["delta"]))
    if len(Ns) < 2:
        return dict(verdict="inconclusive", reason="fewer than two resolutions with delta fits")
    w = pn[str(Ns[-1])]["windows"][str(start)]
    w2 = pn[str(Ns[-2])]["windows"][str(start)]
    g, d = w["gradient"], w["delta"]
    aic_ft = g["finite_time"]["aic"]
    aic_ed = min(g["exponential"]["aic"], g["double_exponential"]["aic"])
    c1 = d["finite_time_collapse"]["aic"] < d["exponential_decay"]["aic"]
    c2 = aic_ft < aic_ed - 2
    Tg, eg = g["finite_time"]["params"]["T"], g["finite_time"]["errors"]["T"]
    Td, ed = d["finite_time_collapse"]["params"]["T"], d["finite_time_collapse"]["errors"]["T"]
    c3 = abs(Tg - Td) <= 2 * np.hypot(eg, ed)
    g2 = w2["gradient"]["finite_time"]
    c4 = abs(Tg - g2["params"]["T"]) <= 2 * np.hypot(eg, g2["errors"]["T"])
    d2 = w2["delta"]["finite_time_collapse"]
    a = not c1
    b = (Td - d2["params"]["T"]) > 2 * np.hypot(ed, d2["errors"]["T"])
    c = not c2
    if c1 and c2 and c3 and c4:
        v = "blow-up indicated"
    elif (a or b) and c:
        v = "no blow-up indicated"
    else:
        v = "inconclusive"
    return dict(verdict=v, N_top=Ns[-1], N_next=Ns[-2],
                criteria=dict(delta_collapse_beats_exp=bool(c1), grad_finite_time_clearly_best=bool(c2),
                              Tstar_grad_delta_agree=bool(c3), Tstar_grad_stable_in_N=bool(c4),
                              delta_exp_at_least_as_good=bool(a), Tstar_delta_moves_later_with_N=bool(b),
                              grad_finite_time_not_clearly_best=bool(c)),
                Tstar_grad=[Tg, eg], Tstar_delta=[Td, ed],
                Tstar_grad_nextN=[g2["params"]["T"], g2["errors"]["T"]],
                Tstar_delta_nextN=[d2["params"]["T"], d2["errors"]["T"]])


def convergence(rs):
    out = []
    for a, b in zip(rs[:-1], rs[1:]):
        da, db = a["data"], b["data"]
        ta = np.round(da["t"], 6)
        tb = np.round(db["t"], 6)
        common, ia, ib = np.intersect1d(ta, tb, return_indices=True)
        if common.size == 0:
            continue
        rg = np.abs(da["G"][ia] - db["G"][ib]) / db["G"][ib]
        rd = np.abs(da["delta"][ia] - db["delta"][ib]) / db["delta"][ib]
        div = np.nonzero(rg > CONV_TOL)[0]
        t_div = float(common[div[0]]) if div.size else None
        ma = resolved_mask(a)
        t_res_a = float(da["t"][ma][-1])
        out.append(dict(Ncoarse=a["meta"]["N"], Nfine=b["meta"]["N"],
                        t_overlap=[float(common[0]), float(common[-1])],
                        t_grad_differs_1e_3=t_div, t_coarse_resolved_until=t_res_a,
                        max_rel_grad_diff_while_coarse_resolved=float(
                            np.max(rg[common <= t_res_a])) if np.any(common <= t_res_a) else None,
                        max_rel_delta_diff_while_coarse_resolved=float(
                            np.nanmax(rd[common <= t_res_a])) if np.any(
                            (common <= t_res_a) & np.isfinite(rd)) else None,
                        t=common.tolist(), rel_grad=rg.tolist(), rel_delta=rd.tolist()))
    return out


def _model_table(fs):
    rows = []
    for name, f in fs.items():
        pe = ", ".join(f"{k}={v:.4g}±{f['errors'][k]:.2g}(SE)±{f.get('errors_total', {}).get(k, np.nan):.2g}(total)"
                       for k, v in f["params"].items() if k != "t0")
        if f.get("T_at_upper_bound"):
            pe += "  [T* at upper bound: no finite T*]"
        rows.append(f"{name}: {pe}  AIC={f['aic']:.1f}")
    return rows


def make_figures(runs, summary, figdir="results/figures"):
    os.makedirs(figdir, exist_ok=True)
    cases = [c for c in SQG + [EULER] if case_runs(runs, c)]
    colors = {512: "C0", 1024: "C1", 2048: "C2", 4096: "C3", 8192: "C4", 128: "C5", 256: "C6"}
    nc = len(cases)
    ncol = 4
    nrow = int(np.ceil(nc / ncol))

    for key, fname, ylabel in (("G", "gradient_growth.png", r"max $|\nabla\theta|$"),
                               ("delta", "strip_width.png", r"$\delta(t)$")):
        fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow), squeeze=False)
        for ax, c in zip(axs.flat, cases):
            for r in case_runs(runs, c):
                d, N = r["data"], r["meta"]["N"]
                m = resolved_mask(r)
                ax.semilogy(d["t"][m], d[key][m], "-", color=colors.get(N), label=f"N={N}")
                ax.semilogy(d["t"][~m], d[key][~m], ":", color=colors.get(N))
                if key == "delta":
                    ax.fill_between(d["t"][m], np.maximum(d["delta"][m] - d["delta_err"][m], 1e-6),
                                    d["delta"][m] + d["delta_err"][m], color=colors.get(N), alpha=.2)
                    ax.axhline(2 * r["meta"]["dx"], color=colors.get(N), lw=.6, ls="--")
            ax.set_title(c)
            ax.set_xlabel("t")
            ax.set_ylabel(ylabel)
            ax.legend(fontsize=7)
        for ax in list(axs.flat)[nc:]:
            ax.axis("off")
        fig.suptitle(("solid: resolved (δ ≥ 2Δx); dotted: under-resolved"
                      + ("; dashed lines: 2Δx" if key == "delta" else "")), fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, fname), dpi=130)
        plt.close(fig)

    # conservation
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow), squeeze=False)
    for ax, c in zip(axs.flat, cases):
        for r in case_runs(runs, c):
            d, N = r["data"], r["meta"]["N"]
            ax.semilogy(d["t"], np.abs(d["dL2"]) + 1e-17, "-", color=colors.get(N), label=f"L2 N={N}")
            ax.semilogy(d["t"], np.abs(d["dH"]) + 1e-17, "--", color=colors.get(N), label=f"H N={N}")
        ax.set_title(c + ": |relative drift|")
        ax.set_xlabel("t")
        ax.legend(fontsize=6, ncol=2)
    for ax in list(axs.flat)[nc:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "conservation.png"), dpi=130)
    plt.close(fig)

    # convergence
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.4 * nrow), squeeze=False)
    for ax, c in zip(axs.flat, cases):
        for cv in summary["cases"][c]["convergence"]:
            ax.semilogy(cv["t"], np.array(cv["rel_grad"]) + 1e-17,
                        label=f"{cv['Ncoarse']} vs {cv['Nfine']}", color=colors.get(cv["Ncoarse"]))
            if cv["t_coarse_resolved_until"] is not None:
                ax.axvline(cv["t_coarse_resolved_until"], color=colors.get(cv["Ncoarse"]), ls=":", lw=.8)
        ax.axhline(CONV_TOL, color="k", lw=.5)
        ax.set_title(c + ": |ΔG|/G between N and 2N")
        ax.set_xlabel("t")
        ax.legend(fontsize=7)
    for ax in list(axs.flat)[nc:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "convergence.png"), dpi=130)
    plt.close(fig)

    # fits on composite, per case
    for c in cases:
        rs = case_runs(runs, c)
        t, arr = composite(rs)
        fs = summary["cases"][c]["fits"]["composite"]["windows"][str(MAIN_START)]
        t0 = fs["t0"]
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
        axs[0].semilogy(t, arr[:, 0], "k.", ms=3, label="data (highest resolved N)")
        tt = np.linspace(t0, t[-1], 300)
        g = fs["gradient"]
        p = g["exponential"]["params"]
        axs[0].semilogy(tt, p["a"] * np.exp(p["b"] * tt), label="exp")
        p = g["double_exponential"]["params"]
        axs[0].semilogy(tt, p["a"] * np.exp(p["b"] * np.exp(p["c"] * tt)), label="double exp")
        p = g["finite_time"]["params"]
        axs[0].semilogy(tt, p["B"] * ((p["T"] - tt) / (p["T"] - p["t0"])) ** (-p["gamma"]), "--",
                        label=f"(T*-t)^-γ, T*={p['T']:.2f}")
        axs[0].axvline(t0, color="gray", lw=.5)
        axs[0].set_xlabel("t")
        axs[0].set_ylabel(r"max $|\nabla\theta|$")
        axs[0].legend(fontsize=8)
        axs[0].set_title(c + " gradient fits")
        ok = np.isfinite(arr[:, 1])
        axs[1].errorbar(t[ok], arr[ok, 1], yerr=arr[ok, 2], fmt="k.", ms=3, lw=.5)
        axs[1].set_yscale("log")
        if fs["delta"]:
            p = fs["delta"]["exponential_decay"]["params"]
            axs[1].semilogy(tt, p["a"] * np.exp(p["b"] * tt), label="δ0 e^{-βt}")
            p = fs["delta"]["finite_time_collapse"]["params"]
            axs[1].semilogy(tt, p["B"] * ((p["T"] - tt) / (p["T"] - p["t0"])) ** p["p"], "--",
                            label=f"c(T*-t)^p, T*={p['T']:.2f}")
        axs[1].set_xlabel("t")
        axs[1].set_ylabel("δ")
        axs[1].legend(fontsize=8)
        axs[1].set_title(c + " strip-width fits")
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, f"fits_{c}.png"), dpi=130)
        plt.close(fig)

    # T* vs N
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    for i, c in enumerate(cases):
        pn = summary["cases"][c]["fits"]["per_N"]
        Ns, Tg, Tge, Td, Tde = [], [], [], [], []
        for N, fs in pn.items():
            w = fs["windows"][str(MAIN_START)]
            Ns.append(int(N))
            f = w["gradient"]["finite_time"]
            Tg.append(f["params"]["T"])
            Tge.append(f["errors"]["T"])
            if w["delta"]:
                f = w["delta"]["finite_time_collapse"]
                Td.append(f["params"]["T"])
                Tde.append(f["errors"]["T"])
            else:
                Td.append(np.nan)
                Tde.append(np.nan)
        axs[0].errorbar(Ns, Tg, yerr=Tge, marker="o", capsize=3, label=c)
        axs[1].errorbar(Ns, Td, yerr=Tde, marker="o", capsize=3, label=c)
    for ax, ttl in zip(axs, ("T* from (T*-t)^-γ fit to max|∇θ|", "T* from c(T*-t)^p fit to δ")):
        ax.set_xscale("log", base=2)
        ax.set_xlabel("N")
        ax.set_ylabel("fitted T*")
        ax.set_yscale("log")
        ax.set_title(ttl + " (resolved window of each N)")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "Tstar_vs_N.png"), dpi=130)
    plt.close(fig)

    # spectra of CMT at top N
    rs = case_runs(runs, "A_cmt")
    if rs:
        fig, ax = plt.subplots(figsize=(6, 4.5))
        for r in rs[-2:]:
            d = r["data"]
            idx = np.linspace(0, len(d["t"]) - 1, 6).astype(int)
            for i in idx:
                E = r["spectra"][i]
                k = np.arange(len(E))
                ax.semilogy(k[1:], E[1:], lw=.8,
                            label=f"N={r['meta']['N']} t={d['t'][i]:.2f}")
        ax.set_xlabel("|k|")
        ax.set_ylabel("E(k)")
        ax.set_ylim(1e-34, 10)
        ax.legend(fontsize=6)
        ax.set_title("A_cmt shell spectra")
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, "spectra_A_cmt.png"), dpi=130)
        plt.close(fig)


def summarize(outdir="results/runs", resdir=None):
    resdir = resdir or os.path.dirname(os.path.abspath(outdir))
    runs = load_runs(outdir)
    summary = {"settings": dict(resolved_if_delta_ge_dx=RESOLVED_DX, window_starts=WINDOW_STARTS,
                                main_window_start=MAIN_START, conv_tol=CONV_TOL),
               "runs": {}, "cases": {}}
    for name, r in runs.items():
        d, m = r["data"], r["meta"]
        msk = resolved_mask(r)
        summary["runs"][name] = dict(
            meta=m, t_start=float(d["t"][0]), t_end=float(d["t"][-1]),
            t_end_resolved=float(d["t"][msk][-1]) if msk.any() else None,
            gradmax_end_resolved=float(d["gradmax"][msk][-1]) if msk.any() else None,
            max_abs_dL2=float(np.max(np.abs(d["dL2"]))), max_abs_dH=float(np.max(np.abs(d["dH"]))),
            max_abs_dL2_resolved=float(np.max(np.abs(d["dL2"][msk]))),
            max_abs_dH_resolved=float(np.max(np.abs(d["dH"][msk]))),
            thmax_drift_resolved=float(np.max(np.abs(d["thmax"][msk] / d["thmax"][0] - 1))))
    for c in SQG + [EULER]:
        rs = case_runs(runs, c)
        if not rs:
            continue
        summary["cases"][c] = dict(fits=fits_for_case(rs), convergence=convergence(rs))
        summary["cases"][c]["verdict"] = verdict(summary["cases"][c])
    # time-step check
    summary["dt_check"] = []
    for name, r in runs.items():
        if r["meta"]["cfl"] != 0.4:
            ref = runs.get(name.split("_cfl")[0])
            if ref is None:
                continue
            ta, tb = np.round(r["data"]["t"], 6), np.round(ref["data"]["t"], 6)
            common, ia, ib = np.intersect1d(ta, tb, return_indices=True)
            msk = resolved_mask(ref)[ib]
            rg = np.abs(r["data"]["G"][ia] - ref["data"]["G"][ib]) / ref["data"]["G"][ib]
            rd = np.abs(r["data"]["delta"][ia] - ref["data"]["delta"][ib]) / ref["data"]["delta"][ib]
            summary["dt_check"].append(dict(
                run=name, ref=ref["meta"]["name"], steps=int(r["meta"]["steps"]),
                ref_steps=int(ref["meta"]["steps"]),
                max_rel_grad_diff_resolved=float(rg[msk].max()),
                max_rel_delta_diff_resolved=float(np.nanmax(rd[msk]))))
    make_figures(runs, summary, figdir=os.path.join(resdir, "figures"))
    # strip bulky arrays from json
    slim = json.loads(json.dumps(summary))
    for c in slim["cases"].values():
        for cv in c["convergence"]:
            for k in ("t", "rel_grad", "rel_delta"):
                cv.pop(k)
    with open(os.path.join(resdir, "summary.json"), "w") as f:
        json.dump(slim, f, indent=1, default=float)
    # human-readable digest
    lines = []
    for c, cs in slim["cases"].items():
        lines.append(f"## {c}")
        vd = cs["verdict"]
        lines.append(f"- VERDICT (PLAN rule): {vd['verdict']}  {json.dumps(vd.get('criteria', vd.get('reason')))}")
        for N, fs in list(cs["fits"]["per_N"].items()) + [("composite", cs["fits"]["composite"])]:
            w = fs["windows"][str(MAIN_START)]
            lines.append(f"- N={N}: resolved to t={fs['t_end_resolved']:.2f}; window t>={w['t0']:.2f}"
                         + ("  [window starts clipped at hand-off: window spread not informative]"
                            if fs.get("windows_clipped") else "")
                         + ("  [composite: delta fit spans grid seams, NOT citable]" if N == "composite" else ""))
            lines += ["    grad  " + s for s in _model_table(w["gradient"])]
            lines += ["    delta " + s for s in _model_table(w["delta"])]
            Ts = [fs["windows"][k]["gradient"]["finite_time"]["params"]["T"] for k in fs["windows"]]
            Tds = [fs["windows"][k]["delta"]["finite_time_collapse"]["params"]["T"]
                   for k in fs["windows"] if fs["windows"][k]["delta"]]
            lines.append(f"    T*(grad) over window starts {WINDOW_STARTS}: "
                         + ", ".join(f"{x:.3g}" for x in Ts))
            lines.append("    T*(delta) over window starts: " + ", ".join(f"{x:.3g}" for x in Tds))
        for cv in cs["convergence"]:
            lines.append(f"- convergence {cv['Ncoarse']} vs {cv['Nfine']}: overlap {cv['t_overlap']}, "
                         f"coarse resolved until {cv['t_coarse_resolved_until']:.2f}, grad differs >1e-3 "
                         f"from t={cv["t_grad_differs_1e_3"]}, max rel grad diff while resolved "
                         f"{cv['max_rel_grad_diff_while_coarse_resolved']}")
    with open(os.path.join(resdir, "summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    return slim
