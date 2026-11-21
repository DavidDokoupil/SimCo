"""Timeouts / out-of-memory counts per tool and benchmark (Table "tab:to_mo" of the LPAR 2026 paper) and size comparisons.

usage (from experiments/):
  python3 summarize.py                       # the data as published (results/paper)
  python3 summarize.py --corrected           # with the corrected SimCo-Degen column on SYNTCOMP (results/corrected)
  python3 summarize.py --results results/run # a new run of run.sh

Needs pandas. Per input automaton one autcross CSV, <benchmark>/<idx>.csv.
"""
import argparse
import glob
import re

import pandas as pd

TOOLS = ["SimCo", "SimCo-Degen", "Autfilt", "Kofola", "Seminator"]


def load(folder, bench):
    rows = []
    for f in glob.glob(f"{folder}/{bench}/*.csv"):
        m = re.search(r"(\d+)\.csv$", f)
        if not m:
            continue
        d = pd.read_csv(f)
        d["idx"] = int(m.group(1))
        rows.append(d)
    d = pd.concat(rows, ignore_index=True)
    d["status"] = d["exit_status"].where(d["exit_status"].isin(["ok", "timeout"]), "MO")
    d.loc[d.status != "ok", "output.states"] = float("nan")
    return d.drop_duplicates(["idx", "tool"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results/paper")
    ap.add_argument("--corrected", action="store_true", help="replace SimCo-Degen on SYNTCOMP by results/corrected")
    a = ap.parse_args()
    table = {}
    for bench, name in (("random", "Random"), ("syntcomp", "SYNTCOMP")):
        d = load(a.results, bench)
        if a.corrected and bench == "syntcomp":
            c = load("results/corrected", bench)
            d = pd.concat([d[d.tool != "SimCo-Degen"], c[c.tool == "SimCo-Degen"]], ignore_index=True)
        for st, label in (("timeout", "Timeouts"), ("MO", "MOs")):
            table[f"{label} in {name}"] = {t: int(((d.tool == t) & (d.status == st)).sum()) for t in TOOLS}
        w = d.pivot(index="idx", columns="tool", values="output.states")
        ok = w[["SimCo", "SimCo-Degen"]].notna().all(axis=1)
        x = w[ok]
        print(f"{name}: SimCo vs SimCo-Degen on the {int(ok.sum())} inputs both finish: SimCo smaller "
              f"{int((x.SimCo < x['SimCo-Degen']).sum())}, larger {int((x.SimCo > x['SimCo-Degen']).sum())}, "
              f"states {int(x.SimCo.sum())} vs {int(x['SimCo-Degen'].sum())}")
    print()
    print(pd.DataFrame(table).T[TOOLS].to_string())


if __name__ == "__main__":
    main()
