"""Which runs make up the campaign (PLAN.md 'Resolutions and scale-down')."""
from __future__ import annotations

import json
import os

import numpy as np

from .ics import FAMILIES
from .run_one import run

SQG = list(FAMILIES)
EULER = "G_euler_cmt"
TEND = 20.0
MAX_WALL = {4096: 1.5 * 3600, 8192: 8.0 * 3600}


def family_of(case):
    return "A_cmt" if case == EULER else case


def alpha_of(case):
    return 2.0 if case == EULER else 1.0


def run_name(case, N, cfl=0.4):
    s = f"{case}_N{N}"
    return s if cfl == 0.4 else s + f"_cfl{cfl}"


def load_meta(path):
    d = np.load(path, allow_pickle=False)
    return json.loads(str(d["meta"]))


def plan_runs(nmax):
    base_top = min(nmax, 2048)
    base = [base_top // 4, base_top // 2, base_top]
    cont = [N for N in (4096, 8192) if N <= nmax]
    return base, cont


def execute(nmax=8192, cases=None, force=False, outdir="results/runs", statedir="results/states",
            log=print):
    cases = cases or SQG + [EULER]
    base, cont = plan_runs(nmax)
    done = []

    def go(case, N, start=None, cfl=0.4):
        name = run_name(case, N, cfl)
        path = os.path.join(outdir, name + ".npz")
        if os.path.exists(path) and not force:
            log(f"skip {name} (exists)")
            return load_meta(path)
        m = run(name, family_of(case), N, alpha=alpha_of(case), tend=TEND, cfl=cfl,
                start=start, outdir=outdir, statedir=statedir, log=log,
                max_wall=MAX_WALL.get(N))
        done.append(m)
        log(f"done {name}: {m['stop_reason']} at {m['steps']} steps, {m['wall_s']:.0f}s")
        return m

    for case in cases:
        for N in base:
            m = go(case, N)
        # time-step check: half the CFL number (0.2 vs 0.4), CMT only, smallest and largest base grid
        if case == "A_cmt":
            go(case, base[0], cfl=0.2)
            go(case, base[-1], cfl=0.2)
        prev = m
        for N in cont:
            if N == 8192 and case not in ("A_cmt",):
                break
            if prev["stop_reason"] != "delta<dx" or prev["handoff_t"] is None:
                log(f"{case}: no continuation to {N} (parent stop: {prev['stop_reason']})")
                break
            start = os.path.join(statedir, f"{prev['name']}_handoff.npz")
            prev = go(case, N, start=start)
    return done
