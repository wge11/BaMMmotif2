# BaMM!motif 2

**Ba**yesian **M**arkov **M**odel **motif** discovery software, version 2.

[![CI](https://github.com/wge11/BaMMmotif2/actions/workflows/ci.yml/badge.svg)](https://github.com/wge11/BaMMmotif2/actions/workflows/ci.yml)

BaMM!motif learns higher-order Bayesian Markov models (BaMMs) of transcription
factor binding motifs from ChIP-seq, ATAC-seq or other sets of enriched
sequences. It starts from seed motifs (PWMs, BaMMs or binding sites), refines
them with expectation maximisation (EM) or collapsed Gibbs sampling, evaluates
them by cross-validation against sampled background sequences and reports
motif occurrences with p-values.

(C) Johannes Söding, Wanwan Ge, Anja Kiesel, Matthias Siebert

| Tool | Purpose |
|---|---|
| `BaMMmotif` | refine seed motifs into BaMMs, optionally evaluate them (`--FDR`) and scan the input (`--scoreSeqset`) |
| `BaMMScan` | scan sequences for occurrences of given motifs |
| `FDR` | evaluate given motifs by cross-validation (precision, recall, p-values) |
| `BaMMSimu` | sample background sequences, or embed/mask motifs in sequences |
| `extractProbs` | convert a BaMM (`.ihbcp`) into probabilities (`.ihbp`) |
| `R/*.R` | performance scores and plots (AUSFC, ROC, PR curves, logos, positional distribution) |
| `py/*.py` | motif comparison and format conversion |

For seed motifs we recommend our de novo motif discovery tool
[PEnG-motif](https://github.com/soedinglab/PEnG-motif).

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [BaMMmotif options](#bammmotif-options)
- [Other command-line tools](#other-command-line-tools)
- [Output files](#output-files)
- [Downstream analysis](#downstream-analysis)
- [BaMM file format](#bamm-file-format)
- [Reproducibility and performance](#reproducibility-and-performance)
- [Development](#development)
- [Citation](#citation)

## Installation

### Requirements

- a C++14 compiler: GCC 5 or later, or Clang 3.4 or later
- [CMake](https://cmake.org/) 3.10 or later
- [Boost](https://www.boost.org/) headers (Boost.Math)
- OpenMP (optional, for multi-threading; included with GCC)

On Debian/Ubuntu:

```bash
sudo apt-get install build-essential cmake libboost-dev
```

On macOS with [Homebrew](https://brew.sh/):

```bash
brew install cmake boost gcc
export CXX=g++-14   # use your Homebrew GCC version; Apple clang works but builds without OpenMP
```

### Build and install

```bash
git clone https://github.com/wge11/BaMMmotif2.git
cd BaMMmotif2
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$HOME/opt/BaMM"
cmake --build build --parallel
ctest --test-dir build          # optional: run all tools on the example data
cmake --install build
```

This installs the executables and the R scripts into `$HOME/opt/BaMM/bin` and
the Python scripts into `$HOME/opt/BaMM/share/bamm/py`. Add the `bin`
directory to your `PATH`, for example in `~/.bashrc`:

```bash
export PATH="$PATH:$HOME/opt/BaMM/bin"
```

Build options (pass with `-D<option>=ON`):

| Option | Effect |
|---|---|
| `CMAKE_BUILD_TYPE` | `Release` (default), `RelWithDebInfo`, `Debug` or `ASan` (AddressSanitizer) |
| `BAMM_NATIVE` | optimise for the CPU of the build machine (`-march=native`); the binaries may not run on other machines |
| `BAMM_WERROR` | treat compiler warnings as errors |

### Dependencies of the helper scripts

- R scripts: R with the packages `argparse`, `zoo`, `fdrtool`, `LSD` and
  `gdata` (`install.packages(c("argparse", "zoo", "fdrtool", "LSD", "gdata"))`;
  the `argparse` package also needs Python).
- Python scripts: Python 3.8 or later with `numpy`, `scipy` and, for
  `filterPWM.py`, `scikit-learn` (`pip install -r py/requirements.txt`).

## Quick start

The `example/` directory contains 300 JunD ChIP-seq peaks (`JunD.fasta`) and
seed PWMs from PEnG-motif (`PWM_peng10.meme`).

```bash
# 1. refine the first two seed PWMs into 2nd-order BaMMs with EM,
#    evaluate them by 4-fold cross-validation and scan the input sequences
BaMMmotif results example/JunD.fasta --PWMFile example/PWM_peng10.meme --maxPWM 2 \
    --EM --FDR --scoreSeqset

# 2. performance scores (AUSFC, pAUC, AUPRC) and evaluation plots
evaluateBaMM.R results JunD --SFC 1 --ROC5 1 --PRC 1

# 3. sequence logos of order 0, 1 and 2
for k in 0 1 2; do plotBaMMLogo.R results JunD_motif_1 $k; done

# 4. positional distribution of the motif occurrences
plotMotifDistribution.R results JunD
```

To scan other sequences with a learned model:

```bash
BaMMScan scan_out peaks.fasta --BaMMFile results/JunD_motif_1.ihbcp \
    --bgModelFile results/JunD.hbcp
```

## BaMMmotif options

```text
BaMMmotif OUTDIR SEQFILE (--PWMFile FILE | --BaMMFile FILE | --bindingSiteFile FILE) [OPTIONS]
```

`OUTDIR` is created if necessary; `SEQFILE` is a FASTA file with the positive
sequences. `BaMMmotif --help` prints the same list.

**Input**

| Option | Description | Default |
|---|---|---|
| `--alphabet STRING` | `STANDARD` (ACGT), `METHYLC` (ACGTM), `HYDROXYMETHYLC` (ACGTH) or `EXTENDED` (ACGTMH) | `STANDARD` |
| `--ss` | search the given strand only (not recommended for ChIP-seq data) | both strands |
| `--negSeqFile FILE` | FASTA file with background sequences; `BaMMmotif` only reports it in the summary and samples its own background set (`BaMMScan` learns its background model from it) | – |
| `--basename STRING` | prefix of all output files | basename of `SEQFILE` |

**Initial models** (exactly one is required)

| Option | Description |
|---|---|
| `--PWMFile FILE` | position weight matrices in [MEME format](https://meme-suite.org/meme/doc/meme-format.html) |
| `--BaMMFile FILE` | a BaMM (`.ihbcp`); needs `--bgModelFile` when scoring without optimisation |
| `--bindingSiteFile FILE` | binding sites of equal length, one per line |
| `--maxPWM INT` | number of motifs from `--PWMFile` to use (default: all) |

**Motif model**

| Option | Description | Default |
|---|---|---|
| `-k, --order INT` | model order | 2 |
| `-a, --alpha FLOAT...` | order-specific prior strengths; overrides `-b` and `-r` | 1 for k = 0, β·γ<sup>k</sup> for k > 0 |
| `-b, --beta FLOAT` | β in α<sub>k</sub> = β·γ<sup>k</sup> | 7 |
| `-r, --gamma FLOAT` | γ in α<sub>k</sub> = β·γ<sup>k</sup> | 3 |
| `--extend INT [INT]` | add uniform positions to both ends, or to the left and right end separately (`--extend 0 2`) | 0 |
| `-q FLOAT` | prior fraction of sequences that contain the motif | 0.3 |

**Background model**

| Option | Description | Default |
|---|---|---|
| `-K, --Order INT` | background model order | 2 |
| `-A, --Alpha FLOAT...` | prior strengths | 1 for k = 0, 10 for k > 0 |
| `--bgModelFile FILE` | read the background model from a `.hbcp` file | learned from `SEQFILE` |

**Optimisation** (without `--EM` or `--CGS` the initial model is used as it is)

| Option | Description |
|---|---|
| `--EM` | expectation maximisation |
| `--CGS` | collapsed Gibbs sampling (100 iterations) |
| `--noInitialZ` | CGS: start from random motif positions instead of one E-step |
| `--noZSampling` | CGS: do not sample motif positions |
| `--noQSampling` | CGS: do not sample the motif fraction q |
| `--noAlphaOpti` | CGS: do not optimise the prior strengths α |
| `--GibbsMH` | CGS: sample α with Metropolis–Hastings |
| `--dissample` | CGS: sample α from a discretised posterior |

**Evaluation**

| Option | Description | Default |
|---|---|---|
| `--FDR` | cross-validate the models and write precision/recall statistics (`.zoops.stats`) | off |
| `-n, --cvFold INT` | number of cross-validation folds | 4 |
| `-m, --mFold INT` | background sequences per positive sequence; raised automatically to give at least 5,000 | 1 |
| `-s, --sOrder INT` | k-mer order used to sample background sequences | 2 |
| `--mops` | also evaluate the multiple-occurrences-per-sequence model | off |
| `--zoops BOOL` | evaluate the zero-or-one-occurrence-per-sequence model | 1 |

**Motif occurrences**

| Option | Description | Default |
|---|---|---|
| `--scoreSeqset` | write motif occurrences with p- and E-values (`.occurrence`) | off |
| `--pvalCutoff FLOAT` | p-value cutoff for reported occurrences | 1e-4 |

**Output and performance**

| Option | Description | Default |
|---|---|---|
| `--saveBaMMs` | also write k-mer counts (`.counts`) and motif positions (`.positions`) | off |
| `--saveInitialBaMMs` | write the initial models (`_init_motif_<i>.ihbcp/.ihbp`) | off |
| `--savePvalues` | write p-values of the cross-validation scores (`.zoops.pvalues`) | off |
| `--saveLogOdds` | write log-odds scores of the positive and background sets | off |
| `--savePRs BOOL` | write `.zoops.stats` with `--FDR` | 1 |
| `--threads INT` | number of OpenMP threads | 4 |
| `--verbose` | print the progress of every iteration | off |
| `-h, --help` | print the help | |

## Other command-line tools

Each tool prints its options with `-h`.

```bash
# scan sequences with PWMs (writes <basename>_motif_<i>.occurrence) or a BaMM (<basename>.occurrence)
BaMMScan OUTDIR SEQFILE --PWMFile motifs.meme [--pvalCutoff 1e-4] [--negSeqFile bg.fasta]
BaMMScan OUTDIR SEQFILE --BaMMFile model.ihbcp --bgModelFile model.hbcp

# evaluate motifs by cross-validation (writes .zoops.stats, as BaMMmotif --FDR)
FDR OUTDIR SEQFILE --PWMFile motifs.meme [--EM | --CGS] [--cvFold 4] [--threads 4]

# sample background sequences, or embed / mask a motif in the input sequences
BaMMSimu OUTDIR SEQFILE --sampleBgset [-s 2] [-m 10]
BaMMSimu OUTDIR SEQFILE --BaMMFile model.ihbcp --embedSeqset [-q 0.5] [--at 50]

# probabilities (.ihbp, .hbp) from conditional probabilities (.ihbcp, .hbcp)
extractProbs OUTDIR model.ihbcp background.hbcp
```

## Output files

For an input `JunD.fasta` (or `--basename JunD`) BaMMmotif writes:

| File | Content | Written |
|---|---|---|
| `JunD.hbcp`, `JunD.hbp` | background model (conditional probabilities, probabilities) | always |
| `JunD_motif_<i>.ihbcp`, `.ihbp` | refined motif model *i* | always |
| `JunD_motif_<i>.zoops.stats` | TP, FP, FDR, recall and p-value per rank; header: background/positive ratio and motif occurrence fraction | `--FDR` |
| `JunD_motif_<i>.occurrence` | motif occurrences: sequence, length, strand, start..end, pattern, p-value, E-value | `--scoreSeqset` |
| `JunD_motif_<i>.counts`, `.positions` | k-mer counts; motif positions with responsibility ≥ 0.3 | `--saveBaMMs` |
| `JunD_motif_<i>.alphas` | learned prior strengths α (CGS) | `--saveBaMMs --CGS` |
| `JunD_motif_<i>.zoops.pvalues` | p-values of the cross-validation scores | `--savePvalues` |

Positions in `.occurrence` files are 1-based on the scanned sequence. With both
strands (the default), positions 1…L are the forward strand and positions
L+2…2L+1 the reverse complement; `strand` is `+` or `-` accordingly.
`py/occur2bed.py` converts the positions into genomic BED coordinates.

## Downstream analysis

The R scripts are installed next to the executables. On/off options take
`1`/`0` or `TRUE`/`FALSE`.

### Performance scores and evaluation plots

`evaluateBaMM.R` needs the `.zoops.stats` files written by `BaMMmotif --FDR`
or `FDR`. It calculates the area under the sensitivity–FDR curve (AUSFC), the
partial ROC AUC up to 5 % FPR (pAUC) and the area under the precision–recall
curve (AUPRC) and writes them to `<prefix>.bmscore`.

```bash
evaluateBaMM.R INPUT_DIR PREFIX [--SFC 1] [--ROC5 1] [--PRC 1]
```

`PREFIX` selects the files `INPUT_DIR/PREFIX*.zoops.stats`, e.g. `JunD` for
all motifs of a run. The options add plots:

| Option | Plot |
|---|---|
| `--SFC 1` | sensitivity–FDR curve |
| `--ROC5 1` | partial ROC curve up to 5 % false positive rate |
| `--PRC 1` | precision–recall curve |

![Sensitivity-FDR curve](example/images/JunD_motif_1_SFC.jpeg)
![Partial ROC curve](example/images/JunD_motif_1_pROC.jpeg)
![Precision-recall curve](example/images/JunD_motif_1_PRC.jpeg)

`plotSFC.R INPUT_DIR PREFIX` and `plotPvalStats.R INPUT_DIR PREFIX --plots 1`
produce further fdrtool-based evaluation plots.

### Sequence logos

```bash
plotBaMMLogo.R INPUT_DIR PREFIX ORDER [--revComp 1] [--stamp 1]
```

`PREFIX` selects `INPUT_DIR/PREFIX*.ihbcp` (with the matching `.ihbp`), `ORDER`
is the logo order (0, 1 or 2, at most the model order). `--revComp 1` plots the
reverse complement (order 0), `--stamp 1` omits the axes.

![Logo of order 0](example/images/JunD_motif_1-logo-order-0.png)
![Logo of order 1](example/images/JunD_motif_1-logo-order-1.png)
![Logo of order 2](example/images/JunD_motif_1-logo-order-2.png)

### Positional distribution of motif occurrences

```bash
plotMotifDistribution.R INPUT_DIR PREFIX
```

needs the `.occurrence` files from `BaMMmotif --scoreSeqset` or `BaMMScan`
and writes `<prefix>_motif_<i>_distribution.png`. It assumes that all input
sequences have the same length.

![Motif distribution](example/images/JunD_motif_1_ds_distribution.jpeg)

### Python scripts

Run from the `py/` directory (they import `utils.py`).

| Script | Purpose |
|---|---|
| `BaMMmatch.py QUERY.meme DB_DIR OUT.tsv` | compare motifs with a database of BaMMs (`DB_DIR/*/*.ihbcp`); reports p- and E-values of similar motifs |
| `filterPWM.py IN.meme OUT.meme` | remove redundant PWMs by clustering similar ones (affinity propagation) |
| `pwm2bamm_local.py IN.meme [-o DIR]` | convert PWMs into 0th-order BaMM files |
| `bamm2pwm.py MODEL.ihbcp OUT.meme` | convert the 0th order of a BaMM into a PWM with an IUPAC name |
| `occur2bed.py FILE.occurrence [-o DIR]` | convert occurrences into BED coordinates; needs FASTA headers of the form `chr1:1000-1205` (as written by `bedtools getfasta`) |

## BaMM file format

Each inhomogeneous (motif) BaMM of order *K* and length *W* is written to two
files with the same layout: `.ihbp` holds the probabilities and `.ihbcp` the
conditional probabilities. Blank lines separate motif positions; within a
position, line *k*+1 holds the values for order *k* = 0…*K*, with (*k*+1)-mers in
lexicographic order (A, C, G, T; AA, AC, …, TT; AAA, …).

`.ihbp` (order 2, one position *j*):

```text
Pj(A)   Pj(C)   Pj(G)   Pj(T)
Pj(AA)  Pj(AC)  Pj(AG)  ... Pj(TT)
Pj(AAA) Pj(AAC) Pj(AAG) ... Pj(TTT)
```

`.ihbcp` (order 2, one position *j*):

```text
Pj(A)     Pj(C)     Pj(G)     Pj(T)
Pj(A|A)   Pj(C|A)   Pj(G|A)   ... Pj(T|T)
Pj(A|AA)  Pj(C|AA)  Pj(G|AA)  ... Pj(T|TT)
```

The homogeneous background BaMM uses the same layout for a single position:
`.hbp` holds the probabilities and `.hbcp` the conditional probabilities. Its
first two lines are comments with the order and the prior strengths
(`# K = 2`, `# A = 1 10 10`).

## Reproducibility and performance

- Results do not depend on the number of threads: sums over sequences are
  accumulated in a fixed order and random numbers are drawn serially, so two
  runs with the same input and options give identical output.
- EM, Gibbs sampling, scoring and cross-validation use OpenMP; set the number
  of threads with `--threads`.
- Build in `Release` mode (the default) for production use; `-DBAMM_NATIVE=ON`
  adds CPU-specific optimisations.

## Development

```bash
cmake -S . -B build -DBAMM_WERROR=ON && cmake --build build --parallel
ctest --test-dir build --output-on-failure        # C++ tools vs tests/reference
python -m unittest discover -s tests              # Python scripts
ruff check py tests && ruff format --check py tests
```

`tests/run_examples.sh BIN_DIR OUT_DIR` runs every tool on the example data and
`tests/compare_outputs.py` compares two output directories with a numeric
tolerance, which is useful for checking that a change preserves results. After
an intended change of results, regenerate `tests/reference` from a new run
(gzip each output file). Continuous integration (GitHub Actions) builds with
GCC and Clang and runs both test suites. See [CHANGELOG.md](CHANGELOG.md) for
changes.

## Citation

If you use BaMM!motif, please cite:

- Ge W, Meier M, Roth C, Söding J. Bayesian Markov models improve the
  prediction of binding motifs beyond first order. *NAR Genomics and
  Bioinformatics* 3(2):lqab026 (2021).
  [doi:10.1093/nargab/lqab026](https://doi.org/10.1093/nargab/lqab026)
- Siebert M, Söding J. Bayesian Markov models consistently outperform PWMs at
  predicting motifs in nucleotide sequences. *Nucleic Acids Research*
  44(13):6055–6069 (2016).

## License

BaMM!motif is released under the GNU General Public License v3 or later; see
[LICENSE](LICENSE). It includes GetOpt_pp (GPLv3) by Daniel Gutson and
`cmake/FindASan.cmake` (MIT) by Matthew Arsenault.

## Contact

Bug reports and questions are welcome as
[GitHub issues](https://github.com/wge11/BaMMmotif2/issues).
