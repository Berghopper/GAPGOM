# Changelog of GAPGOM
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).\*

\*Changes of commit versions will not be kept track of, instead only major and minor versions will keep a changelog;

x.y.z (e.g. 1.2.5) --> changes of 1.2.z and 1.3.z will be kept track of but not the small "z" version changes, because this changes with every commit.
All z changes will be kept track of in the unreleased minor (y) version. There was an exception for the first version however, because this version didn't
use correct convention yet (0.0.1). Failed implementations might be re-evaluated in later updates.

## [0.4.0 - Unreleased]

---

### Planned changes/features (subject to changes) *=not yet implemented/finished

### Added

- Installation chapter to main vignette
- Extra tests
- Basic user input checks
- Tests for checking precalculated matrices versions.

### Changed

- `apply` functions
- How test results are stored -> they are now stored in a list instead of seperately.

### Removed

- All global assignments

### Fixed

- Precalculated topoicsim matrix
- Compliance with R CMD check
- Compliance with R CMD BiocCheck
- Dependency issues

## [0.3.0 - Vignette Benchy] - 2018-11-22

---


### Added

- Re-added data which lncRNApred's paper is based on.
- Support for all AnnotationDbi Key types.
- Vignette for GAPGOM.
- Vignette for benchmarks/performance measuring.
- Term algorithm as exported function.
- Set go data as exported function.
- Some frequently used params can now be set manually (go_data, All_Go_pairs).
- Seperate scoring function (If people want to do their own enrichment/only get semantic similarity scores).
- Some more arguments for the algorithms to make the user have more control.
- Custom gene interface for topoicsim.

### Changed

- `DESCRIPTION` file, now adheres to some smaller requirements Bioconductor enforces.
- Data preparation functions, they are now updated and more clear.
- Topoicsim for gene and geneset are now 1 generic function `topo_ic_sim_genes`.
- Documentation and comments (again).
- Argument parsing for topoicsim.
- Some tests and their results.
- Used arguments in both algorithms.
- Precalculated values are deprecated for this version and thus turned off, warning gets shown if it is tried to be used anyway. The values will be updated in upcoming releases and probably adhere to only certain `org.db` package versions (with a warning shown if incorrect versions are used).
- Some `lapply` loops, now for loops, as they were using global assignments and this is forbidden/unfavorable in conventions.
- Updated data-raw files, they now work properly for current version and are better documented.

### Removed

- All forms of parallelization support. It was buggy/inconsistent and where it could be implemented insignificant performance gains were made.
- Fantom5 example data
- Direct interface from fantom5 file to expset `fantom_load_expressionset()`.

### Fixed

- Inconsistent result of lncrnapred compared to original algorithm. Now is exactly the same.
- Documentation of internal functions is now hidden.
- `data-raw` missing in `.Rbuildignore`.
- `R CMD check` failing because of different lncRNApred results in factor indices. Fixed by converting unnecesary factors to chars.
- `R CMD check` failing because of `T` and `F` usage instead of `TRUE`/`FALSE`.
- `fantom_to_expset()` not always loading correctly; mouse and human had different metadata.

## [0.2.0 - Polished Potato] - 2018-10-20

---


### Added

- TopoICSim performance improvements.
- LncRNApred performance improvements.
- Parallelisation for combined method
- Some small extra options for similarity prediction.
- Precalculated similarity scores between frequent GO terms.
- Data preperation file, describing methods used for generating package data.
- Progress bars for TopoICSim.
- Basic tests (for main algorithms only).
- Some small tests for the main algorithms.
- LICENSE file.

### Changed

- Return value of main TopoICSim between gene sets to better match return value of TopoICSim between two genes.
- All functions are now pre-compiled (except for data_gen.R).
- Imports (per function).
- Comments, a lot of them.
- Implementation of entrez -> goid lookup. (some ids get looked up twice)
- `data.frame` subsetting/ddply functions --> replaced with way faster data.table alternatives.
- Vectorisation of parts that weren't yet vectorized.
- Changed some loops back to for loops instead of apply --> some weren't a good use case for apply.
- Some matrixes now use the Matrix library.

### Removed

- Nothing.

### Fixed

- R depends version (3.4.4 --> 3.5.1).
- Performance bug in TopoICSim (#4).
- Performance bug in entrez -> goid lookup (#5).
- `fdr` bug (#7).
- Bug where id_select_vector was actually selected for rather than filtered for (#8).
- Topoicsim being a total memory hog, it left a lot of unused items in ram --> using `gc()` now as an optional argument.
- Newline loading bar bug (#9).


## [0.1.0 - Technically Tidied] - 2018-10-04

---


### Added

- Changelog (NEWS.md).
- Better input interface and refactoring towards ExpressionSet based system.
- Replacement for current example data by ExpressionSet class of fantom5 annotated data (small subset of fantom5).
- Documentation file for previously mentioned new data.
- EntrezID based system for looking up GO ids instead of hardcoded translation table without source.

### Changed

- Small refactor to switch case in expression_prediction_function().
- A lot of refactoring in expression_prediction_function() and enrichment_analysis() to make it work with ExpressionSets.
- Documentation of lncRNApred related functions

### Removed

- Old example/expression data and their documentation, this will be re-added later once proper source is known.

### Fixed

- Bumped up y version number to correct Bioconductor version conventions.
- Small typo in zzz.R
- .gitignore mistake (still including .Rdata temp files... These got really big :( )
- File typo data_preperation.R --> data_preparation.R


## [0.0.1 - baby steps] - 2018-09-25

---


### Added

- Expression metric prediction algorithm made by Rezvan Ehsani and Finn Drablos.
- TopoICSim algorithm/measure made by Rezvan Ehsani and Finn Drablos.
- Documentation to the previously denoted algorithms and metrics.
- Default datasets to compute some examples with.