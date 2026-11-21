#!/bin/bash
# Runs the experiments of the LPAR 2026 paper on one benchmark.
#
# usage (from experiments/):  ./run.sh random|syntcomp [JOBS]
#   environment: PYTHON (default python3, must import Spot 2.14.5), KOFOLA, SEMINATOR (binaries), TOOLS (default:
#   "SimCo SimCo-Degen Kofola Seminator Autfilt"; e.g. TOOLS=SimCo-Degen runs only that configuration)
#   ONLY="1 5 7" runs only these inputs
#   output: results/run/<benchmark>/<idx>.csv (autcross CSV per input automaton, idx 1-based as autfilt -N)
#
# Settings as in the paper: 16 GB per tool, timeout 60 s (Random) / 1 h (SYNTCOMP). Outputs are postprocessed by
# `autfilt --small --low`, except SimCo's on SYNTCOMP (SimCo's own --postprocess only), as in the paper.
#
# SimCo-Degen degeneralizes the input with `autfilt -b` and complements the result with SimCo --merge --postprocess.
# The command of the paper's SYNTCOMP run passed the original file to SimCo as well (`... | simco ... %H`), so SimCo
# read the TGBA and the degeneralization was lost; see README.md, "Correction". This script uses the corrected pipeline.
set -e
bench=$1; jobs=${2:-4}
PYTHON=${PYTHON:-python3}; KOFOLA=${KOFOLA:-kofola}; SEMINATOR=${SEMINATOR:-seminator}
TOOLS=${TOOLS:-"SimCo SimCo-Degen Kofola Seminator Autfilt"}
simco="$PYTHON $(cd "$(dirname "$0")/.." && pwd)/SimCo/complement.py"
case $bench in
  random)   input=benchmarks/random/TGBA_new.hoa; timeout=60; spost="| autfilt --small --low" ;;
  syntcomp) input=$(mktemp --suffix=.hoa); gzip -dc benchmarks/syntcomp/SYNTCOMPAutomata.hoa.gz > "$input"
            timeout=3600; spost="" ;;
  *) echo "usage: $0 random|syntcomp [JOBS]"; exit 1 ;;
esac
declare -A cmd=(
  [SimCo]="$simco --jump --merge --postprocess %H 2>/dev/null $spost"
  [SimCo-Degen]="autfilt -b %H | $simco --merge --postprocess 2>/dev/null $spost"
  [Kofola]="$KOFOLA --complement %H 2>/dev/null | autfilt --small --low"
  [Seminator]="$SEMINATOR --complement %H 2>/dev/null | autfilt --small --low"
  [Autfilt]="autfilt --complement %H 2>/dev/null | autfilt --small --low"
)
tools=()
for t in $TOOLS; do tools+=("{$t} bash -c \"ulimit -v 16777216; ${cmd[$t]} >%O\""); done
out=results/run/$bench; mkdir -p "$out"
n=$(autfilt -c "$input")
export input timeout out
# GNU parallel cannot export arrays; the tool commands go through a file
printf '%s\n' "${tools[@]}" > "$out/tools.txt"
job_file() {  # IDX
  mapfile -t t < "$out/tools.txt"
  autfilt -N $1..$1 "$input" | autcross -F - --timeout=$timeout --ignore-execution-failures --no-check \
      --csv="$out/$1.csv" "${t[@]}" >/dev/null 2>&1
}
export -f job_file
{ if [ -n "$ONLY" ]; then printf "%s\n" $ONLY; else seq 1 "$n"; fi; } | parallel -j "$jobs" --joblog "$out/joblog.txt" job_file {}
[ "$bench" = syntcomp ] && rm -f "$input"
echo "$n automata -> $out"
