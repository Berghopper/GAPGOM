# this code serves as a test to see if the package's algorithm matches up with the original (prediction_of_lncRNA_ORIGINALSCRIPT.R).
# For this it also uses the original lncRNA2Function dataset. While the original source is unavaible, there's still an intermediary version available.

## ORIGINAL DATA

# link: https://tare.medisin.ntnu.no/pred_lncRNA/ based on lncRNA2Function: http://bio-annotation.cn/lncrna2function/

# expression data containing fpkm expression values. fpkm is fragements per kilobase million. fragments means that is is for paired-end data.
options(stringsAsFactors = FALSE)
ExpressionData <- read.table('/media/casper/USB_ccpeters/internship_thesis/data/rezvan_lncrna2function/lncRNA2function_data.txt', sep='\t', head=TRUE)

# Gene ID's and annotation term type.
EnsemblID2GOID <- read.table('/media/casper/USB_ccpeters/internship_thesis/data/rezvan_lncrna2function/EG2GO.txt', sep='\t', head=TRUE)

## Conversion to expressionset

expression_matrix <- as.matrix(ExpressionData[,4:ncol(ExpressionData)])
rownames(expression_matrix) <- ExpressionData$GeneID
featuredat <- as.data.frame(ExpressionData[,1:3])
rownames(featuredat) <- ExpressionData$GeneID
expset <- ExpressionSet(expression_matrix, 
                        featureData = new("AnnotatedDataFrame", 
                                          data=featuredat))

## selection of filter

# keep everything that is a protein coding gene
filter_vector <- pData(featureData(expset))[(pData(featureData(expset))$GeneType=="protein_coding"),]$GeneID

## make translational data.frame (BP only)
bp_only <- EnsemblID2GOID[(EnsemblID2GOID[ ,3] == "biological_process"), ]
id_translation_df <- bp_only[,1:2]
colnames(id_translation_df) <- c("ORIGID","GO")

## algorithm run

# set gid and run.
gid <- "ENSG00000228630"

result <- GAPGOM::expression_prediction(gid, 
                                        expset,
                                        "human", 
                                        "BP",
                                        id_select_vector = filter_vector,
                                        id_translation_df = id_translation_df,
                                        method = "pearson", verbose = TRUE, 
                                        filter_pvals = TRUE)
