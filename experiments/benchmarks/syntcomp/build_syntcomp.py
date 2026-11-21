"""Builds the SYNTCOMP benchmark (SYNTCOMPAutomata.hoa) from tlsf.csv, as for the LPAR 2026 paper.

usage (from experiments/):  python3 benchmarks/syntcomp/build_syntcomp.py

tlsf.csv: one line per formula, "<tlsf file>, <LTL formula>": the TLSF specifications of
https://github.com/SYNTCOMP/benchmarks/tree/master/tlsf converted to LTL with syfco (default instances of the
parametric specifications); 959 formulae. Each is translated with ltl2tgba (Spot) under a 60 s timeout (ltldo);
automata with an empty language and automata with a single acceptance set are dropped. With Spot 2.14.5 this gave the
209 TGBAs of the paper (SYNTCOMPAutomata.hoa.gz). The translation timeout makes the result depend on the machine, so
the file used in the paper is included.
"""
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
out_path = HERE / "SYNTCOMPAutomata.hoa"
count = 0
with open(HERE / "tlsf.csv") as f, open(out_path, "w") as out, tempfile.NamedTemporaryFile("w", suffix=".ltl") as tmp:
    rows = [line.strip().split(",", 1) for line in f if line.strip()]
    for i, (source, formula) in enumerate(rows, 1):
        print(f"\r{i}/{len(rows)}  kept: {count}", end="", flush=True)
        tmp.seek(0); tmp.truncate(); tmp.write(formula.strip()); tmp.flush()
        aut = subprocess.run(["ltldo", "--timeout=60", "--errors=ignore",
                              f"ltl2tgba --name='{source.strip()}' -F %F", "-F", tmp.name],
                             capture_output=True, text=True).stdout
        if not aut.strip():
            continue
        aut = subprocess.run(["autfilt", "-v", "--is-empty"], input=aut, capture_output=True, text=True).stdout
        aut = subprocess.run(["autfilt", "--acc-sets=2.."], input=aut, capture_output=True, text=True).stdout
        n = aut.count("--END--")
        if n:
            out.write(aut); out.flush(); count += n
print(f"\n{count} automata -> {out_path}")
