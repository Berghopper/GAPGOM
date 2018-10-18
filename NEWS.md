# Changelog of GAPGOM
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).\*

\*Changes of commit versions will not be kept track of, instead only major and minor versions will keep a changelog;

x.y.z (e.g. 1.2.5) --> changes of 1.2.z and 1.3.z will be kept track of but not the small "z" version changes, because this changes with every commit.
All z changes will be kept track of in the unreleased minor (y) version. There was an exception for the first version however, because this version didn't
use correct convention yet (0.0.1). Failed implementations might be re-evaluated in later updates.

## [0.2.0 - Unreleased]
### Planned changes/features (subject to changes) *=not yet implemented
### Added
- TopoICSim performance improvements
- lncRNApred performance improvements
- parallelisation for combined method
- some small extra options for similarity prediction
- precalculated similarity scores between frequent GO terms.
- data preperation file, describing methods used for generating package data.
- progress bars for TopoICSim.
- basic tests (for main algorithms only)
- Some small tests for the main algorithms
### Changed
- Return value of main TopoICSim between gene sets to better match return value of TopoICSim between two genes.
- All functions are now pre-compiled (except for data_gen.R)
- Imports (per function)
- Comments, a lot of them
- Implementation of entrez -> goid lookup. (some ids get looked up twice.) 
- data.frame subsetting/ddply functions --> replaced with way faster data.table alternatives
- Vectorisation of parts that aren't yet vectorized
- Changed some loops back to for loops instead of apply --> some weren't a good use case for apply.
- some matrixes now use the Matrix library.
### Removed
- Nothing.
### Fixed
- R depends version (3.4.4 --> 3.5.1)
- Performance bug in TopoICSim (#4)
- Performance bug in entrez -> goid lookup (#5)
- fdr bug (#7)
- bug where id_select_vector was actually selected for rather than filtered for. (#8)
- topoicsim being a total memory hog, it left a lot of unused items in ram --> using `gc()` now as an optional argument.


## [0.1.0 - Technically Tidied] - 2018-10-04
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
- small typo in zzz.R
- .gitignore mistake (still including .Rdata temp files... These got really big :( )
- file typo data_preperation.R --> data_preparation.R


## [0.0.1 - baby steps] - 2018-09-25
### Added
- Expression metric prediction algorithm made by Rezvan Ehsani and Finn Drablos.
- TopoICSim algorithm/measure made by Rezvan Ehsani and Finn Drablos.
- Documentation to the previously denoted algorithms and metrics.
- Default datasets to compute some examples with.