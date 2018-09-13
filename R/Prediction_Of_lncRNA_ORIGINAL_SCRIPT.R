######################################################################################################################
#
#
#   - This code predicts functions of un-annotated genes based on Gene Ontology and correlated expression patterns.
#     - Rezvan Ehsani <rezvanehsani74@gmail.com>
#     - Finn Drablos <finn.drablos@ntnu.no>
#
#   - You need to enter the path of two data files in the code:
#     - ExpressionData - Gene_ID GO_ID GO_type
#       - ("ENSG00000275323	GO:0003723	molecular_function"
#     - EnsemblID2GOID - Gene_ID Gene_type Gene_name FPKM FPKM ...
#       - ("ENSG00000185097	protein_coding	OR4F16	9.305e-17	3.320169e-15 ...")
#
#   - Main function is Prediction Function;
#
#                       Input: 1- GeneID - Ensembl ID of the gene that you want to predict, e.g. ENSG00000228630
#                              2- Onto - Ontology type, currently two options; MF(Molecular Function) or BP(Biological Process)
#                              3- Method - Currently five options; Pearson, Spearman, Fisher, Sobolev, combine
#
#                       Output: the output is a list of ontology terms, ordered with respect to FDR values
#                              1- GOID - Gene Ontology ID
#                              2- Ontology - Ontology type (MF or BP)
#                              3- FDR - False Positive Rate
#                              4- Term - description of GOID
#
#   - There is an example for HOTAIR lncRNA in last line of code.
#
######################################################################################################################

# load libs
library(plyr)
library(GOSim)

# set stringAsFactors to false so
options(stringsAsFactors = FALSE)
# read data from files
# expression data containing fpkm expression values. fpkm is fragements per kilobase million. fragments means that is is for paired-end data.
ExpressionData = read.table('/media/casper/USB_ccpeters/workspace internship project/pred_lncRNA/lncRNA2function_data.txt', sep='\t', head=TRUE)
# Gene ID's and annotation term type.
EnsemblID2GOID = read.table('/media/casper/USB_ccpeters/workspace internship project/pred_lncRNA/EG2GO.txt', sep='\t',head=TRUE)
# expressiondata_pc => protein coding only. So known functional genes.
ExpressionData_PC <- ExpressionData[(ExpressionData$GeneType == "protein_coding"), ]

# grab ensembl IDs of protein coding expression data.
EnsemblID_PC <- ExpressionData_PC$GeneID
# set cutoff. 250 has the best performance for GO term prediction.
CutOff <- 250
# l = Total amount of collumn in expressiondata.
l <- ncol(ExpressionData)

# lncRNA similarity metrics.
SobolevMetric <- function(x, y) {
            x1 <- x**2 / (sum(x**2))
            y1 <- y**2 / (sum(y**2))
	    z1 <- x1 - y1
	    FT <- fft(z1)
	    w <- 2*pi*(1:length(FT))/(length(FT))
	    s <- sum((1+w)*abs(FT)**2)**(1/2)
	    return(s)
	    }

FisherMetric <- function(x, y) {
	    x1 <- x**2 / (sum(x**2))
	    y1 <- y**2 / (sum(y**2))
	    t <- x1 * y1
	    s <- acos(sum(sqrt(t)))
	    return(s)
            }

# convert an ext ID to a term ID
# dat = EG_ (Extracted genes with correct gene ontology and narrowed selection of top genes.)
# List_Top_Genes (List of top 250 genes)
ExtID2TermID <- function(dat, List_Top_Genes){
  # add 1 for each gene that is in list top genes.
	return( ddply(dat, .(go_id), function(x) nrow(x[(x$ensembl_gene_id %in% List_Top_Genes),])))
	}

# function for enrichment
# We use an enrichment analysis to identify enriched Gene Ontology terms in the co-expressed gene sets
# These are then used to predict GO terms for the un-annotated gene.
# DF_ is the scoring/metrix dataframe from novel/chosen gene compared to all other genes.
Enrichment_func <- function(DF_, onto) {
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


# Main function is Prediction Function;
#
#                       Input: 1- GeneID - Ensembl ID of the gene that you want to predict, e.g. ENSG00000228630
#                              2- Onto - Ontology type, currently two options; MF(Molecular Function) or BP(Biological Process)
#                              3- Method - Currently five options; Pearson, Spearman, Fisher, Sobolev, combine
#
#                       Output: the output is a list of ontology terms, ordered with respect to FDR values
#                              1- GOID - Gene Ontology ID
#                              2- Ontology - Ontology type (MF or BP)
#                              3- FDR - False Positive Rate
#                              4- Term - description of GOID
#
Prediction_Function<-function(GeneID, Onto, Method){

	# Target expression data
  Target_EX <- ExpressionData[( ExpressionData[,1] == GeneID ),]
	Tareget_EX <- as.numeric( Target_EX[1,c(4:l)] ) # spelling error tarEget.

	# these functions find score between target expression of target gene vs the rest of the desired protein coding genes.
	# this is correlation for spearman and pearson, but the fisher and sobolev metrics are metrics of their own.
	###
	if( Method == "Pearson" )       SCORE <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) abs(cor(as.numeric(x), Tareget_EX)))
        if( Method == "Spearman" )      SCORE <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) abs(cor(as.numeric(x), Tareget_EX, method = "spearman")))
	if( Method == "FisherMetric" )  SCORE <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) SobolevMetric(as.numeric(x), Tareget_EX))
	if( Method == "SobolevMetric" ) SCORE <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) FisherMetric(as.numeric(x), Tareget_EX))
	if( Method == "combine" ) {
				  SCORE_Pearson  <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) abs(cor(as.numeric(x), Tareget_EX)))
				  SCORE_Spearman <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) abs(cor(as.numeric(x), Tareget_EX, method = "spearman")))
                                  SCORE_Fisher   <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) FisherMetric(as.numeric(x), Tareget_EX))
				  SCORE_Sobolev  <- apply( ExpressionData_PC[, c(4:l)], 1, function(x) SobolevMetric(as.numeric(x), Tareget_EX))
				}
  ###

	# DFScore is the "score" (correlation/metric) dataframe against target gene.
	if( Method == "combine" ) DFScore <- data.frame( EnsemblID_PC, SCORE_Pearson, SCORE_Spearman, SCORE_Sobolev, SCORE_Fisher )
	else DFScore <- data.frame( EnsemblID_PC, SCORE )

	# filter NA's
	DFScore <- na.omit(DFScore)
	# filter gene ID's (Select everything except the chosen gene).
	DFScore <- DFScore[(DFScore[,1]!=GeneID),]

	# run the enrichments analysis'.
	if( Method == "Pearson" | Method == "Spearman")  EnrichResult<-Enrichment_func( DFScore[ rev(order(DFScore[,2])), ], Onto)
	if( Method == "SobolevMetric" | Method == "FisherMetric") EnrichResult<-Enrichment_func( DFScore[order(DFScore[,2]), ], Onto)
	if( Method == "combine" ) {
	      # Run enrichment for each method
				EnrichPearson  <- Enrichment_func( DFScore[ rev(order(DFScore$SCORE_Pearson)), ], Onto)
				EnrichSpearman <- Enrichment_func( DFScore[ rev(order(DFScore$SCORE_Spearman)), ], Onto)
				EnrichSobolev  <- Enrichment_func( DFScore[     order(DFScore$SCORE_Sobolev) , ], Onto)
				EnrichFisher   <- Enrichment_func( DFScore[     order(DFScore$SCORE_Fisher) , ], Onto)
				# bind all results by row
				EnrichCombine  <- rbind( EnrichPearson, EnrichSpearman, EnrichSobolev,  EnrichFisher )
				# and keep the rows with the lowest corrected P-values.
				EnrichCombine  <- ddply(EnrichCombine, .(GOID), function(x) x[which.min(x$FDR),])
				# Order p-values on lowest > bigest
				EnrichResult   <- EnrichCombine[ order(as.numeric(EnrichCombine$FDR)), ]
				}
	if(nrow(EnrichResult)>0){
	    # number the rownames and return the enrichment results.
			rownames(EnrichResult) <- c(1:nrow(EnrichResult))
			return(EnrichResult)
			}
	else{print("Could not find anything!")
	}
}

# Example for prediction of HOTAIR lncRNA
start.time <- Sys.time()
result <- Prediction_Function( GeneID="ENSG00000228630", Onto="BP", Method="combine" )
end.time <- Sys.time()
time.taken <- end.time - start.time

time.taken

