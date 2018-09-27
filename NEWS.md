# Changelog of GAPGOM
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).\*

\*Changes of commit versions will not be kept track of, instead only major and minor versions will keep a changelog;

x.y.z (e.g. 1.2.5) --> changes of 1.2.z and 1.3.z will be kept track of but not the small "z" version changes, because this changes with every commit.
All z changes will be kept track of in the unreleased minor version. There was an exception for the first version however, because this version didn't
use correct convention yet (0.0.1).

## [0.1.0 - Unreleased]
### Planned changes/features (subject to changes)
- Better input interface and refactoring towards ExpressionSet based system (this will be a lot of work...).
- replacement for current example data by ExpressionSet class of fantom5 annotated data.
### Added
- Changelog (NEWS.md).
### Changed
- Small refactor to switch case in expression_prediction_function().
- Bumped up y version number to correct Bioconductor convention.

## [0.0.1] - 2018-09-25
### Added
- Expression metric prediction algorithm made by Rezvan Ehsani and Finn Drablos.
- TopoICSim algorithm/measure made by Rezvan Ehsani and Finn Drablos.
- Documentation to the previously denoted algorithms and metrics.
- Default datasets to compute some examples with.