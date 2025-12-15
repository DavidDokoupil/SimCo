# Instructions

## Shared Requirements

- Spot
  - Source: <https://spot.lre.epita.fr/>
  - Suggested version - 2.14.2, because the used version Kofola runs into segfaults with 2.14.2.
- Kofola
  - Version used in thesis - <https://github.com/VeriFIT/kofola/releases/tag/v1.0.0> (commit 7bbf4a1).
  - Install in the `tools\` directory or adjust the Jupyter Notebooks.
- Seminator 2
  - Version used in thesis - <https://github.com/adl/seminator> (commit dbb71b3).
  - **The standard Seminator 2 repository does not contain compatibility fixes for newer versions of Spot**.
  - Installation inside `tools\` directory is recommended, however not required.
- SimCo
  - Version used in thesis - <https://github.com/DavidDokoupil/SimCo> (commit 5e6116a)
  - Install in the `tools\` directory or adjust the Jupyter Notebooks.

## TGBA-Practice Dataset Requirements

- The Jupyter Notebook relies on the dataset - <https://github.com/ondrik/automata-benchmarks> (commit f375761), to be included in the `data\` directory. Since, the repository is quite large the user can also only download the folder - 
<https://github.com/ondrik/automata-benchmarks/tree/master/omega/pecan> and either change the Jupyter Notebook, or respect the path within the repository.

## Additional Information

- The notebook use the folder `data\` for the generation of inputs, `autcross` statistics, custom statistics, and plots. The user is free to use its own structure by modifying the Jupyter Notebooks
- By default:
  - Only plots of the most performant configuration of SimCo are saved. But all comparisons are shown in the notebooks.
  - The notebooks print some basic information during the creation of the custom statistics
