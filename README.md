# Authors

* Finn Drabløs
* Rezvan Ehsani
* Casper Peters @Berghopper

# LICENSE

Copyright (c) [2018] [Finn Drabløs, Rezvan Ehsani, Casper Peters] under the MIT license.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# ABOUT

GAPGOM (novel **G**ene **A**nnotation **P**rediction and other **GO** **M**etrics) is an R package with tools and algorithms for estimating correlation of gene expression enriched terms in gene sets, and semantic distance between sets of gene ontology (GO) terms.
This package has been made for predicting the annotation of un-annotated gene(s), in particular with respect to GO, and testing such predictions. 
The prediction is done by comparing expression patterns between a query gene and a library of annotated genes, and annotate the query gene by enriched terms from the set of genes with similar expression pattern (often described as "guilt by association").

For more info, read the package vignette after building and installing the package (This will later be available on Bioconductor).

# Installation

simply run `BiocManager::install("GAPGOM")` with Bioconductor version 3.9 or above.

# Issues

To let us help you better, please conform with the following rules so bug-squashing/issue tracking becomes easier and faster:

## Submitting an issue

- Make sure your issue is not a duplicate, and try search if your issue exists first.
- Please include a reproducable example, along with this, include your `sessionInfo()`.
- Feature requests should start with a header; `## FR ##` so it's easier to mark it as such.

## Commenting on existing issues

Please only comment if you have valuable information to add/can add to further discussion of the problem. Bumping issues with e.g. "+1" isn't allowed.

# Pull requests

- Only send PRs if it falls within the atomic scope of this project.

- Please keep them small and managable :). Also please make sure your PR version complies with all of Bioconductors package submission guidelines.

# Info on news file:

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).\*

\*Changes of commit versions will not be kept track of, instead only major and minor versions will keep a changelog;

x.y.z (e.g. 1.2.5) --> changes of 1.2.z and 1.3.z will be kept track of but not the small "z" version changes, because this changes with every commit.
All z changes will be kept track of in the unreleased minor (y) version. There was an exception for the first version however, because this version didn't
use correct convention yet (0.0.1). Failed implementations might be re-evaluated in later updates.

# Feature roadplan and versions

Here the roadplan of features is described with any upcoming/planned features. ETA's will not be provided, already implemented features will be ticked off. Only major features will be listed here, for further details on package development history, instead, refer to the [NEWS](NEWS) file.


Note: list is subject to changes and is not necessarily definetive.

- [x] Implement GO Annotation Prediction.
- [x] Implement the TopoICSim algorithm.
- [x] Optimization and refactors of original algorithms.
- [x] Adding benchmarks and discussion vignettes of both algorithms.
- [x] Adding documentation and mostly comply with package conventions.
- [x] Become part of the Bioconductor repository.
- [] Add geneset analysis for scoring expression similarities.
- [] Add chromosomal distance correlations on top of regular scoring. (THIS IS ALREADY ADDED AS CONCEPT SCRIPT)
- [] Along with the previously two mentioned GREAT-like features, include Enrichr into the pipeline. (THIS IS ALREADY ADDED AS CONCEPT SCRIPT)
- [] Added AnnotationHub intergration so more annotation data can be used.

# Building package from scratch

Clone the repository with `git clone https://github.com/Berghopper/GAPGOM.git`.
Make a new package project with an existing directory using the latest Rstudio.
After having the project loaded, first make sure you have all the packages from the DESCRIPTION file installed (both Imports and Suggests packages). `stats`, `utils` and `methods` are all standard R packages and packages that reside on cran are: `Matrix`, `plyr`, `magrittr`, `data.table`, `igraph`, `matrixStats`, `testthat`, `pryr`, `knitr`, `rmarkdown`, `prettydoc`, `ggplot2`, `kableExtra`, `profvis`, `reshape2`. The rest of the packages can be found on Bioconductor.

After installing all dependencies, the project can be built/compiled with Rstudio, click on the `build tab` > `more` and choose your option. Building might also require additional packages such as `devtools`. For building vignettes, you can use the `devtools::build_vignettes()` command. The built package should be able to be loaded with regular `install.packges()` using `type="source"`