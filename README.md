# SimCo

SimCo (**Sim**ple **Co**mplementation) is a prototype tool for the complementation of transition-based generalized
Büchi automata (TGBA), implementing the construction of

> David Dokoupil, Jan Strejček: *Simple Complementation of Generalized Büchi Automata*. LPAR 2026.

This repository contains the version of SimCo and the experiments used in that paper (tag `lpar-2026`), with one
correction of the published evaluation, see [Correction](#correction).

## Requirements

- Python 3.14 (the paper used 3.14.4)
- [Spot](https://spot.lre.epita.fr/) 2.14.5 with its Python bindings, built with `--enable-max-accsets=256` (needed
  for the SYNTCOMP benchmark; automata with many acceptance sets). The options used in the paper need an unmodified
  Spot. The experimental options `--sim`, `--sim_useful` and `--sim_input` (not used in the paper) need Spot with
  `patches/spot-2.14.5-simco.patch` (exports Spot's direct simulation relation to Python; `patch -p2` inside the
  unpacked Spot 2.14.5).

## Usage

```
python3 SimCo/complement.py [options] [file.hoa ...]   # reads standard input without a file
```

Every automaton of the input (HOA format) is complemented; the complement is printed in HOA format. The options of the
paper: `--jump` (level jumping), `--merge` (merging of tainted components), `--postprocess` (Spot's simplification
of the result, level Low). `--help` lists all options; the others are experimental.

## Experiments

Everything below runs from `experiments/`.

**Benchmarks** (`benchmarks/`):

- `random/TGBA_new.hoa` – **Random**, 491 TGBAs: Spot's `randaut` with 2–5 acceptance sets and 2–5 atomic
  propositions (50 automata each), 2–7 states, density 0.5, seed 0; automata with an empty language removed.
  `random/build_random.sh` regenerates it byte for byte (checked with Spot 2.14.5 and 2.16.3).
- `syntcomp/SYNTCOMPAutomata.hoa.gz` – **SYNTCOMP**, 209 TGBAs: the LTL formulae of the SYNTCOMP TLSF benchmarks
  (`syntcomp/tlsf.csv`, converted with syfco, default instances of parametric specifications; 959 formulae)
  translated with `ltl2tgba` under a 60 s timeout; empty automata and automata with a single acceptance set removed.
  `syntcomp/build_syntcomp.py` rebuilds it; the translation timeout makes the result machine-dependent, so the file of
  the paper is included.

**Running** (`run.sh random|syntcomp [JOBS]`, needs GNU parallel and Spot's `autcross`): SimCo, SimCo-Degen
(degeneralization with `autfilt -b`, then SimCo `--merge --postprocess`), Kofola, Seminator 2 and Spot's `autfilt
--complement`, 16 GB per tool, 60 s (Random) / 1 h (SYNTCOMP) per automaton. Results go to `results/run/<benchmark>/`.
Paths: `PYTHON` (a Python that imports Spot 2.14.5), `KOFOLA`, `SEMINATOR`. Versions used in the paper: Kofola commit
`f916ecc` ([VeriFIT/kofola](https://github.com/VeriFIT/kofola), branch `devel`), Seminator 2 commit `dbb71b3`
([adl/seminator](https://github.com/adl/seminator)), both built against Spot 2.14.5.

**Results** (`results/`):

- `paper/` – the data of the paper as published: one autcross CSV per input automaton, `<benchmark>/<idx>.csv`
  (idx 1-based as `autfilt -N`). Only the file names and one tool name were normalized (the Random runs called
  SimCo-Degen `SimCo-Classic`); all measured values are unchanged.
- `corrected/` – the corrected SimCo-Degen results on SYNTCOMP (see below).
- `python3 summarize.py [--corrected]` prints the table of timeouts and out-of-memory cases of the paper and compares
  SimCo with SimCo-Degen.

## Correction

In the paper's SYNTCOMP experiment, SimCo-Degen was run as

```
autfilt -b %H | python3 SimCo/complement.py --merge --postprocess %H
```

The trailing `%H` (the input file) made SimCo read the original TGBA instead of the degeneralized automaton coming
from `autfilt -b`. The published SimCo-Degen results on **SYNTCOMP** are therefore those of SimCo with `--merge` only
(without `--jump`) on the original TGBAs: they show the effect of level jumping, not of complementing TGBAs directly
instead of degeneralizing them first. This affects the SimCo-Degen column of the table of timeouts and out-of-memory cases
(SYNTCOMP rows), the SimCo-Degen plot of the SYNTCOMP benchmark in the figure comparing complement sizes, and the statement that only SimCo-Degen produces consistently larger
complements than SimCo on SYNTCOMP. The **Random** experiment is not affected: there the degeneralized automaton reached
SimCo (`autfilt -b %H | python3 SimCo/complement.py --jump --merge --postprocess`; `--jump` has no effect on
degeneralized automata, checked).

`run.sh` uses the intended pipeline. Corrected SimCo-Degen results on SYNTCOMP, with the paper's settings (this
version of SimCo, Spot 2.14.5, 1 h, 16 GB): **pending**.

## License

GNU General Public License v3.0 (`LICENSE`).
