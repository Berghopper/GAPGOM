#' GAPGOM internal - enrichment_analysis()
#'
#' This function is an internal function and should not be called by the user.
#' 
#' Enriches score results from multiple methods to give a better idea of
#' important similarities.
#' This function is specifically made for predicting lncRNA annotation by
#' assuming "guilt by association". For instance, the expression data in this
#' package is actually based on mRNA expression data, but correlated with
#' lncRNA. This expression data is the used in combination with mRNA GO
#' annotation to calculate similarity scores between GO terms,
#'
#' @section Notes:
#' Internal function used in expression_prediction_function().
#'
#' @param ordered_score_df the score dataframe see documentation on
#' GAPGOM::example_score_dataframe for formatting.
#' @param ontology desired ontology to use for prediction. One of three;
#' "BP" (Biological process), "MF" (Molecular function) or "CC"
#' (Cellular Component). Cellular Component is not included with the package's
#' standard data and will thus yield no results.
#' @param expression_data dataframe with at least 1 or more
#' rows. Some of the columns should be named certain ways and the dataframe
#' should contain certain information neccesary for calculation;
#' 1st neccesary col; name; "GeneID" (EnsemblIDs should go here).
#' 2nd neccesary col; name; "GeneType" (genetype as described in the Ensembl
#' database).
#' Next, a range of n-columns long, containing expression values, names do not
#' matter here. (You can use the existing dataset as an example)
#' @param id_translation_df dataframe with corresponding GOIDs,
#' EntrezIDs and General ID defined by ExpressionSet. (rownames of ).
#' @param enrichment_cutoff cutoff number for the amount of genes to be
#' enriched in the enrichment analysis. (default is 250)
#'
#' @return The resulting dataframe with prediction of similar GO terms.
#' These are ordered with respect to FDR values. The following columns will be
#' in the dataframe;
#' GOID - Gene Ontology ID,
#' Ontology - Ontology type (MF or BP),
#' FDR - False Positive Rate,
#' Term - description of GOID.
#' However, unlike in expression_prediction, this dataframe will have unsorted
#' row numbering. And it won't contain used method.
#'
#' @import AnnotationDbi
#' @importFrom plyr ddply .
enrichment_analysis <- compiler::cmpfun(function(
                                id_translation_df,
                                ordered_score_df,
                                organism,
                                ontology,
                                expression_set,
                                id_select_vector,
                                enrichment_cutoff = 250) {
  # extracted_genes -> Extracted genes with correct gene ontology.
  # now filter EG to also extract only genes that are present in user defined
  # expression data rows. 
  extracted_genes <- id_translation_df[(id_translation_df$ORIGID %in%
                                        id_select_vector), ]

  # List of top n (cutoff) genes (Ensembl ID)
  list_top_genes <- ordered_score_df[c(1:enrichment_cutoff), 1] ## FIX IN TOPLEVEL
  # List of gene ontologies given the Extracted genes that are in the top
  # 250 genes of the score dataframe.  for each ensemble ID there are more
  # gene ontologies.
  # list_of_gos, n genes but all unqiue corresponding GO IDs
  # REPLACE this also with help of a retrieval function via gosemsim
  list_of_gos <- extracted_genes[(extracted_genes$ORIGID %in%
                                    list_top_genes), 3]
  list_of_gos <- unique(list_of_gos)
  list_of_gos <- list_of_gos[which(!is.na(list_of_gos))]

  # Count the amount of genes with the same GO ID and quantify them.
  term_id_to_ext_id <- ddply(extracted_genes, .(GO), function(x) nrow(x))
  # set the column name to Freq Anno. (This was broken before.)
  colnames(term_id_to_ext_id)[2] <- "Freq_Anno"

  # Filter the quantification to only have the top genes where the go ID
  # corresponds
  qterm_id_to_ext_id <- term_id_to_ext_id[(term_id_to_ext_id$GO %in%
                                              list_of_gos), ]

  # ---

  # Quantify extracted_genes in List of top genes (grab every go_id for
  # corresponding ensembl IDs)
  quantified_ext_id_to_term_id <- ext_id_to_term_id(extracted_genes,
                                                    list_top_genes)

  # After this, filter it for existing Gene ontologies within the top GOs
  quantified_ext_id_to_term_id <- quantified_ext_id_to_term_id[(
      quantified_ext_id_to_term_id[, 1] %in% list_of_gos), 2]


  #			              GeneList | Genome
  #		                ------------------
  #	In Anno group 	  |   n1   |   n2  |
  #	------------------------------------
  #	Not in Anno group |   n3   |   n4  |
  #	------------------------------------
  #	Where, GeneList is the number of protein-coding genes that co-expressed
  # with lncRNA, Genome is the number of all protein coding-genes,
  #	In Anno group is the number of protein-coding genes that were both
  # co-expressed with lncRNA and annotated in the Term, and Not in Anno group
  #	is the number of protein-coding genes that were co-expressed with lncRNA
  # but were not annotated in the Term

  n1 <- quantified_ext_id_to_term_id
  n2 <- qterm_id_to_ext_id[, 2] - quantified_ext_id_to_term_id
  n3 <- length(unique(id_select_vector)) - enrichment_cutoff -
      qterm_id_to_ext_id[, 2] + quantified_ext_id_to_term_id
  n4 <- rep(enrichment_cutoff, nrow(qterm_id_to_ext_id)) # Issue #1 Bitbucket 

  # now bind into 1 df.
  qterm_id_to_ext_id <- cbind(qterm_id_to_ext_id, n1, n2, n3, n4)
  # select quantification values to at least be 5 for goID quantification.
  qterm_id_to_ext_id <- qterm_id_to_ext_id[(qterm_id_to_ext_id[, 2] >= 5), ]

  # select last 4 columns (n1,n2,n3,n4)
  args.df <- qterm_id_to_ext_id[, c(3:6)]
  # calculate p-values using the hypergeometrix distribution.
  pvalues <- apply(args.df, 1, function(n) min(phyper(0:n[1] - 1, n[2], n[3],
                                                      n[4],
                                                      lower.tail = FALSE)))
  # Grab corresponding go_ids
  go_id <- qterm_id_to_ext_id[, 1]
  # Replicate the ontology for the amount of rows
  ontology <- rep(ontology, nrow(args.df))
  # format the pvalues to have 3 digits at max.
  pvalues_formatted <- format(pvalues, scientific = TRUE, digits = 3)
  # Benjamini & Hochberg multiple comparisons adjustment
  fdr <- p.adjust(pvalues, method = "fdr", n = length(pvalues))
  # grab description of each gene ontology term using Term() from the 
  # annotationDbi package.

  term <- Term(qterm_id_to_ext_id[, 1])
  # Create a dataframe containing all results in a neat format.
  enrichment_dataframe <- data.frame(GOID = go_id, Ontology = ontology,
                                     Pvalue = pvalues_formatted,
                                     FDR = format(fdr, scientific = TRUE,
                                                  digits = 3),
                                     Term = term)

  # Omit NA's
  enrichment_dataframe <- na.omit(enrichment_dataframe)
  # Order the result by FDR (the adjusted p values)
  enrichment_dataframe <- enrichment_dataframe[(order(as.numeric(
      enrichment_dataframe[, 4]))), ]
  # Only keep every P-value that is significant
  enrichment_dataframe <- enrichment_dataframe[(as.numeric(
      enrichment_dataframe[, 4]) < 0.05), ]

  return(enrichment_dataframe)
})

generate_translation_df <- function(expression_set, organism, ontology) {
  entrezid_col <- resolve_entrezid_col(expression_set)
  go_data <- set_go_data(organism, ontology)
  rowtracker = 0
  entrez_go_dfs <- sapply(expression_set@featureData@data[,entrezid_col], function(entrezrawid) {
    rowtracker <<- rowtracker + 1
    #print(rowtracker)
    #print(rownames(expression_set@assayData$exprs)[rowtracker])
    entrez_id <- unlist(strsplit(entrezrawid, ":"))[2]
    goids <- go_data@geneAnno[go_data@geneAnno==entrez_id,]$GO
    if (length(goids) != 0){
      return(data.frame(ORIGID=rownames(expression_set@assayData$exprs)[rowtracker], ENTREZID=entrezrawid, GO=goids))
    }
  })
  entrez_go_df <- unique(as.data.frame(data.table::rbindlist(entrez_go_dfs)))
  return(entrez_go_df)
}



#' resolve function for entrez
resolve_entrezid_col <- function(expression_set) {
  colnames_vector <- colnames(expression_set@featureData@data)
  exp <- regexec(".*entrez.*", colnames_vector)
  regex_result <- unlist(regmatches(colnames_vector, exp))
  if (length(regex_result) < 1) {
    return(NULL)
  } else {
    return(regex_result)
  }
}
