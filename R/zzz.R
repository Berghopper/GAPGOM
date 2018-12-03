##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    pkgVersion <- packageDescription(pkgname, fields = "Version")
    msg <- paste0(pkgname, " v", pkgVersion, "\n",
                  "For help/issues, refer to the readme FAQ or report an ",
                  "issue on the issue page: ",
                  "https://github.com/Berghopper/GAPGOM/issues\n\n")

    citation <- paste0("If you use ", pkgname, " in any sort of publication, ",
                   "please cite:\n", "[1] Ehsani R, Drabl\u00f8s F: TopoICSim:", 
                   " a new semantic similarity measure based on gene ontology.",
                   " BMC Bioinformatics 2016, 17(1):296\n",
                   "[2] Ehsani R, Drabl\u00f8s F: ",
                   "Measures of co-expression for improved function prediction",
                   " of long non-coding RNAs, 2018. BMC Bioinformatics, 2018",
                   " Accepted."
                   )

    packageStartupMessage(paste0(msg, citation))
}
