#' @export
enrichment_analysis <- function(DF_, onto) {
  # EG_ -> Extracted genes with correct gene ontology.
  if(onto == "MF") EG_ <- EnsemblID2GOID[(EnsemblID2GOID[ ,3] == "molecular_function"), ]
  if(onto == "BP") EG_ <- EnsemblID2GOID[(EnsemblID2GOID[ ,3] == "biological_process"), ]

  # now filter EG to also extract only genes that are present in protein coding expression data.
  EG_ <- EG_[(EG_$ensembl_gene_id %in% EnsemblID_PC), ]

  # List of top 250 genes (Ensembl ID)
  List_Top_Genes <- DF_[c(1:CutOff), 1]
  # List of gene ontologies given the Extracted genes that are in the top 250 genes of the score dataframe.
  # for each ensemble ID there are more gene ontologies.
  # ListOfGos, 250 genes but all unqiue corresponding GO IDs
  ListOfGos <- EG_[(EG_$ensembl_gene_id %in% List_Top_Genes),2]
  ListOfGos <- unique(ListOfGos)
  ListOfGos <- ListOfGos[which(!is.na(ListOfGos))]

  # Count the amount of genes with the same GO ID and quantify them.
  TermID2ExtID <- ddply(EG_, .(go_id), function(x) nrow(x))
  # set the column name to Freq Anno. (This was broken before.)
  colnames(TermID2ExtID)[2] <- "Freq_Anno"

  # Filter the quantification to only have the top genes where the go ID corresponds
  qTermID2ExtID <- TermID2ExtID[(TermID2ExtID$go_id %in% ListOfGos), ]

  # ---

  # Quantify EG_ in List of top genes (grab every go_id for corresponding ensembl IDs)
  qExtID2TermID <- ExtID2TermID(EG_, List_Top_Genes)

  # After this, filter it for existing Gene ontologies within the top GOs
  qExtID2TermID <- qExtID2TermID[(qExtID2TermID[ ,1] %in% ListOfGos),2]


  #			  GeneList | Genome
  #		          ------------------
  #	In Anno group 	  |   n1   |   n2  |
  #	------------------------------------
  #	Not in Anno group |   n3   |   n4  |
  #	------------------------------------
  #	Where, GeneList is the number of protein-coding genes that co-expressed with lncRNA, Genome is the number of all protein coding-genes,
  #	In Anno group is the number of protein-coding genes that were both co-expressed with lncRNA and annotated in the Trem, and Not in Anno group
  #	is the number of protein-coding genes that were co-expressed with lncRNA but were not annotated in the Trem

  n1 = qExtID2TermID
  n2 = qTermID2ExtID[ ,2]-qExtID2TermID
  n3 = length(unique(EnsemblID_PC)) - CutOff - qTermID2ExtID[ ,2] + qExtID2TermID
  n4 = rep(CutOff, nrow(qTermID2ExtID))


  # now bind into 1 df.
  qTermID2ExtID <- cbind(qTermID2ExtID, n1, n2, n3, n4)
  # select quantification values to at least be 5 for goID quantification.
  qTermID2ExtID <- qTermID2ExtID[(qTermID2ExtID[ ,2]>=5),]

  # select last 4 columns (n1,n2,n3,n4)
  args.df<-qTermID2ExtID[,c(3:6)]
  # calculate p-values using the hypergeometrix distribution.
  pvalues <- apply(args.df, 1, function(n)
    min(phyper(0:n[1]-1,n[2], n[3], n[4], lower.tail=FALSE)))

  # Grab corresponding GOIDs
  GOID <- qTermID2ExtID[ ,1]
  # Replicate the ontology for the amount of rows
  Ontology <- rep(onto, nrow(args.df))
  # format the pvalues to have 3 digits at max.
  Pvalue <- format(pvalues, scientific=TRUE, digits = 3)
  # Benjamini & Hochberg multiple comparisons adjustment
  fdr  <- p.adjust(pvalues, method = "fdr", n = length(pvalues))
  # grab description of each gene ontology term using Term() from the annotationDbi package.
  TERM <- Term(qTermID2ExtID[ ,1])
  # Create a dataframe containing all results in a neat format.
  D_EN <- data.frame(GOID=GOID, Ontology=Ontology, Pvalue=Pvalue, FDR=format(fdr, scientific=TRUE, digits = 3),Term=TERM)

  # Omit NA's
  D_EN <- na.omit(D_EN)
  # Order the result by FDR (the adjusted p values)
  D_EN <- D_EN[(order(as.numeric(D_EN[ ,4]))), ]
  # Only keep every P-value that is significant
  D_EN <- D_EN[(as.numeric(D_EN[ ,4])<.05),]
  return(D_EN)
}
