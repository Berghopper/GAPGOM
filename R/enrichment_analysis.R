#' @export
enrichment_analysis <- function(ordered_score_df, ontology, ensembl_to_go_id_conversion_df, ensembl_id_pc, enrichment_cutoff) {
    # extracted_genes -> Extracted genes with correct gene ontology.
    if (ontology == "MF")
        extracted_genes <- ensembl_to_go_id_conversion_df[(ensembl_to_go_id_conversion_df[, 3] == "molecular_function"), ]
    if (ontology == "BP")
        extracted_genes <- ensembl_to_go_id_conversion_df[(ensembl_to_go_id_conversion_df[, 3] == "biological_process"), ]

    # now filter EG to also extract only genes that are present in protein coding expression data.
    extracted_genes <- extracted_genes[(extracted_genes$ensembl_gene_id %in% ensembl_id_pc), ]

    # List of top 250 genes (Ensembl ID)
    list_top_genes <- ordered_score_df[c(1:enrichment_cutoff), 1]
    # List of gene ontologylogies given the Extracted genes that are in the top 250 genes of the score dataframe.  for each ensemble ID there are more gene ontologylogies.
    # list_of_gos, 250 genes but all unqiue corresponding GO IDs
    list_of_gos <- extracted_genes[(extracted_genes$ensembl_gene_id %in% list_top_genes), 2]
    list_of_gos <- base::unique(list_of_gos)
    list_of_gos <- list_of_gos[base::which(!base::is.na(list_of_gos))]

    # Count the amount of genes with the same GO ID and quantify them.
    term_id_to_ext_id <- plyr::ddply(extracted_genes, plyr::.(go_id), function(x) base::nrow(x))
    # set the column name to Freq Anno. (This was broken before.)
    base::colnames(term_id_to_ext_id)[2] <- "Freq_Anno"

    # Filter the quantification to only have the top genes where the go ID corresponds
    qterm_id_to_ext_id <- term_id_to_ext_id[(term_id_to_ext_id$go_id %in% list_of_gos), ]

    # ---

    # Quantify extracted_genes in List of top genes (grab every go_id for corresponding ensembl IDs)
    quantified_ext_id_to_term_id <- UGenePred::ext_id_to_term_id(extracted_genes, list_top_genes)

    # After this, filter it for existing Gene ontologylogies within the top GOs
    quantified_ext_id_to_term_id <- quantified_ext_id_to_term_id[(quantified_ext_id_to_term_id[, 1] %in% list_of_gos), 2]


    #			  GeneList | Genome
    #		          ------------------
    #	In Anno group 	  |   n1   |   n2  |
    #	------------------------------------
    #	Not in Anno group |   n3   |   n4  |
    #	------------------------------------
    #	Where, GeneList is the number of protein-coding genes that co-expressed with lncRNA, Genome is the number of all protein coding-genes,
    #	In Anno group is the number of protein-coding genes that were both co-expressed with lncRNA and annotated in the term, and Not in Anno group
    #	is the number of protein-coding genes that were co-expressed with lncRNA but were not annotated in the term

    n1 = quantified_ext_id_to_term_id
    n2 = qterm_id_to_ext_id[, 2] - quantified_ext_id_to_term_id
    n3 = base::length(base::unique(ensembl_id_pc)) - enrichment_cutoff - qterm_id_to_ext_id[, 2] + quantified_ext_id_to_term_id
    n4 = base::rep(enrichment_cutoff, base::nrow(qterm_id_to_ext_id))


    # now bind into 1 df.
    qterm_id_to_ext_id <- base::cbind(qterm_id_to_ext_id, n1, n2, n3, n4)
    # select quantification values to at least be 5 for goID quantification.
    qterm_id_to_ext_id <- qterm_id_to_ext_id[(qterm_id_to_ext_id[, 2] >= 5), ]

    # select last 4 columns (n1,n2,n3,n4)
    args.df <- qterm_id_to_ext_id[, c(3:6)]
    # calculate p-values using the hypergeometrix distribution.
    pvalues <- base::apply(args.df, 1, function(n) base::min(stats::phyper(0:n[1] - 1, n[2], n[3], n[4], lower.tail = FALSE)))

    # Grab corresponding go_ids
    go_id <- qterm_id_to_ext_id[, 1]
    # Replicate the ontology for the amount of rows
    ontology <- base::rep(ontology, nrow(args.df))
    # format the pvalues to have 3 digits at max.
    pvalues_formatted <- base::format(pvalues, scientific = TRUE, digits = 3)
    # Benjamini & Hochberg multiple comparisons adjustment
    fdr <- stats::p.adjust(pvalues, method = "fdr", n = length(pvalues))
    # grab description of each gene ontology term using Term() from the annotationDbi package.
    term <- AnnotationDbi::Term(qterm_id_to_ext_id[, 1])
    # Create a dataframe containing all results in a neat format.
    enrichment_dataframe <- base::data.frame(GOID = go_id, Ontology = ontology, Pvalue = pvalues_formatted, FDR = base::format(fdr, scientific = TRUE, digits = 3), Term = term)

    # Omit NA's
    enrichment_dataframe <- stats::na.omit(enrichment_dataframe)
    # Order the result by FDR (the adjusted p values)
    enrichment_dataframe <- enrichment_dataframe[(base::order(base::as.numeric(enrichment_dataframe[, 4]))), ]
    # Only keep every P-value that is significant
    enrichment_dataframe <- enrichment_dataframe[(base::as.numeric(enrichment_dataframe[, 4]) < 0.05), ]
    return(enrichment_dataframe)
}
