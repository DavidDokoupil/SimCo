#!/bin/bash
# The Random benchmark (TGBA.hoa) as described in the LPAR paper: randaut with -A2..5 acceptance sets and 2..5
# atomic propositions, 50 automata each (800), 2..7 states, density 0.5, default seed 0; automata with an empty
# language removed -> 491 TGBAs. Reproduces TGBA.hoa byte for byte with Spot 2.14.5 and 2.16.3 (checked 2026-11-21).
# usage (from the experiments root): benchmarks/random/build_random.sh
set -e
out=benchmarks/random/TGBA.hoa
tmp=$(mktemp)
for A in 2 3 4 5; do
  for ap in 2 3 4 5; do
    randaut -n50 -Q2..7 $ap -A$A -e0.5 --seed=0 >> "$tmp"
  done
done
autfilt "$tmp" -v --is-empty > "$out"
rm "$tmp"
echo "$(grep -c -- --END-- "$out") automata -> $out"
