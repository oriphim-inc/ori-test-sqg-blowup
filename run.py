"""Reproduce every result:  python run.py            (full campaign, N up to 8192)
                           python run.py --nmax 512  (quick: base grids 128/256/512)
Finished runs in results/runs/ are reused unless --force. --analyze-only skips integration."""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from sqg.campaign import execute, SQG, EULER  # noqa: E402
from sqg.report import summarize  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--nmax", type=int, default=8192)
ap.add_argument("--cases", nargs="*", default=None, help=f"subset of {SQG + [EULER]}")
ap.add_argument("--force", action="store_true")
ap.add_argument("--analyze-only", action="store_true")
ap.add_argument("--outdir", default="results/runs")
ap.add_argument("--statedir", default="results/states")
a = ap.parse_args()

if not a.analyze_only:
    execute(nmax=a.nmax, cases=a.cases, force=a.force, outdir=a.outdir, statedir=a.statedir,
            log=lambda s: print(s, flush=True))
summarize(outdir=a.outdir)
print(open(Path(a.outdir).resolve().parent / "summary.md").read())
