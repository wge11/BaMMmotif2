# Changelog

## Unreleased

### Changed results

These fixes change the output of some runs. Runs not listed here give the
same results as before, up to floating-point rounding (performance scores
such as AUSFC, pAUC and AUPRC agree to three decimals on the example and on
a 4,000-sequence test set).

- **Motif occurrences with the smallest p-values were missing.** The
  exponential tail fit in `BaMMScan` and `BaMMmotif --scoreSeqset` used the
  lowest instead of the highest background scores, so every hit scoring above
  all but 10 background scores got p = inf and was left out of the
  `.occurrence` file. On a 4,000-sequence JunD set the number of reported
  occurrences rises from 3,227 to 4,174.
- **`BaMMmotif --CGS` ignored its prior settings.** The `--verbose` flag was
  passed as the prior strength β, so β was 0 (1 with `--verbose`) instead of
  the documented 7, and `--noInitialZ`, `--noZSampling`, `--noAlphaOpti`,
  `--GibbsMH` and `--dissample` had no effect.
- Gibbs sampling now draws the motif fraction q from its own random number
  generator instead of the global `rand()`, which was not thread-safe. CGS
  results therefore differ from earlier versions (equally valid samples).
- Windows spanning the junction between the forward and reverse-complement
  strand are no longer reported as occurrences (they contained an `N`).
- p-values in `.zoops.stats` and `.occurrence` are capped at 1.
- `FDR --genericNeg` had no effect.

### Performance and reproducibility

- About 2.5x faster EM refinement with cross-validation (8.7 s → 3.4 s
  single-threaded on 4,000 sequences), mainly from replacing the atomic
  updates in the EM M-step with per-block buffers.
- Results are identical for any `--threads` value. Before, two runs with
  the same input could give different models and statistics.
- The default build is now optimised (`Release`); the documented `cmake ..`
  previously built without optimisation, which was about 2x slower.

### Fixed

- Heap-buffer-overflow in the cross-validation precision/recall calculation
  (found with AddressSanitizer).
- `BaMMScan`, `FDR` and `BaMMSimu` crashed at exit when `--alphabet` was given.
- Memory leaks of the sampled background sequences; `extractProbs` checked
  for the wrong number of arguments.
- `BaMMmotif --help` listed options that did not exist and wrong defaults;
  `-m`, `-n` and `-s` are now also accepted without `--FDR`.
- Python: `filterPWM.py` crashed on MEME files with two or more motifs, and
  every MEME model got uniform background frequencies; `occur2bed.py` did not
  run with pandas ≥ 1.0, dropped the first occurrence and produced shifted
  coordinates; `bamm2pwm.py` did not run with current NumPy; `BaMMmatch.py`
  failed with the `spawn`/`forkserver` start methods.
- R: `--SFC 1` style options were rejected by current versions of the
  `argparse` package.

### Build, tests and documentation

- CMake 3.10+, C++14, shared sources built once as libraries, options
  `BAMM_NATIVE` and `BAMM_WERROR`; no compiler warnings.
- `ctest` runs all tools on the example data and compares the results with
  `tests/reference`; Python unit tests and ruff checks.
- GitHub Actions replace Travis CI (travis-ci.org has shut down).
- Rewritten README: installation, quick start, complete option reference,
  output files, downstream analysis and citation.
- Removed `make.sh`, which copied an R script that no longer exists.
