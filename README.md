# Authors #

* Finn Drabløs
* Rezvan Ehsani
* Casper Peters @Berghopper

# LICENSE #

Copyright (c) [2018] [Finn Drabløs, Rezvan Ehsani, Casper Peters] under the MIT license.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# ABOUT #

GAPGOM (novel **G**ene **A**nnotation **P**rediction and other **GO** **M**etrics) is an R package with tools and algorithms for estimating correlation of gene expression enriched terms in gene sets, and semantic distance between sets of gene ontology (GO) terms.
This package has been made for predicting the annotation of un-annotated gene(s), in particular with respect to GO, and testing such predictions. 
The prediction is done by comparing expression patterns between a query gene and a library of annotated genes, and annotate the query gene by enriched terms from the set of genes with similar expression pattern (often described as "guilt by association").

For more info, read the package vignette after building and installing the package (This will later be available on Bioconductor).