# Changelog of GAPGOM
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).\*

\*Changes of commit versions will not be kept track of, instead only major and minor versions will keep a changelog;

x.y.z (e.g. 1.2.5) --> changes of 1.2.z and 1.3.z will be kept track of but not the small "z" version changes, because this changes with every commit.
All z changes will be kept track of in the unreleased minor (y) version. There was an exception for the first version however, because this version didn't
use correct convention yet (0.0.1).

## [0.2.0 - Unreleased]
### Planned changes/features (subject to changes)
### Added
- TopoICSim performance improvements
- lncRNApred performance improvements
- parallelisation
- Rcpp performance tweaks
### Changed
- Implementations of measure calculations -> changed to Rcpp versions.
- Implementation of DAG functions - performance improvements (not sure how yet)
- Implementation of entrez -> goid lookup. (some ids get looked up twice.)
- data.frame subsetting with .subset2()?
### Removed
### Fixed


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