# proof-of-concept script for implementing great-like features 
# like taking into account genomic positions of transcripts and trying to 
# annotate them

# this script is just meant to show off how a potential pipeline would work,
# it might not be 100% accurate.

rm(list=ls())
options(stringsAsFactors = FALSE)

library(GAPGOM)
library(Biobase)
library(biomaRt)
library(HelloRanges)
library(GenomicRanges)

# STEP 1 LOAD DATA

# Load table with basic info about genes etc (http://pheweb.sph.umich.edu:5000/pheno/20001_1044)
# (Only first 50 used.)

geneset_table <- read.csv("~/Downloads/Prostate_UKBiobank_PheWeb_hg19.txt", sep="\t", strip.white = TRUE)
# preloaded here:
geneset_table <- data.frame(list(Position = c("17:46,805,705", "8:128,077,146", 
"10:51,549,496", "17:36,101,156", "11:69,010,651", "11:2,229,782", 
"17:46,114,750", "19:51,361,757", "17:69,110,923", "2:173,309,402", 
"17:45,441,520", "22:43,500,212", "10:74,268,244", "3:170,083,629", 
"2:171,370,372", "12:53,303,331", "4:27,631,917", "5:2,853,074", 
"12:121,665,478", "8:143,451,323", "17:47,460,588", "14:53,380,014", 
"3:25,316,576", "11:65,791,581", "8:54,058,632", "2:234,567,995", 
"2:63,111,004", "7:116,894,926", "10:122,798,413", "1:76,825,837", 
"7:97,774,859", "14:67,367,121", "7:27,976,563", "4:7,134,928", 
"12:97,880,059", "2:122,043,913", "16:81,739,478", "8:82,706,381", 
"17:64,758,694", "17:12,579,132", "6:92,112,016", "8:43,423,953", 
"4:114,175,084", "5:126,859,798", "6:33,937,525", "10:125,687,637", 
"4:106,084,643", "4:115,009,369", "10:64,528,351", "9:2,410,927"
), SNP = c("C / T", "G / A", "T / C", "T / C", "T / G", "C / T", 
"G / A", "T / C", "A / G", "T / C", "C / T", "G / T", "T / G", 
"C / G", "A / G", "T / C", "G / A", "C / T", "C / T", "C / T", 
"C / T", "A / G", "C / G", "C / T", "A / G", "C / T", "G / A", 
"A / G", "A / G", "G / A", "A / G", "G / A", "G / A", "T / G", 
"G / A", "G / A", "C / T", "T / C", "G / T", "C / T", "C / T", 
"C / T", "A / G", "T / A", "T / C", "G / A", "G / A", "C / T", 
"G / A", "C / T"), Variant = c("rs138213197", "rs77541621", "rs10993994", 
"rs7501939", "rs12799883", "rs10743182", "rs554574584", "rs17632542", 
"rs7217073", "rs80353656", "rs188447985", "rs5759167", "rs139990447", 
"rs61436251", "rs934859", "rs17120257", "rs143137339", "rs187177973", 
"rs116932569", "rs186745647", "rs186051044", "rs4901309", "rs189178590", 
"rs139150042", "rs547142665", "rs115416905", "rs2539982", "rs38908", 
"rs2420906", "rs185254575", "rs6977321", "rs181925254", "rs10486567", 
"rs545155977", "rs552249643", "rs556907436", "rs150716819", "rs368791485", 
"rs184470928", "rs148907705", "rs184594722", "rs180742137", "rs187068793", 
"rs186476772", "rs114158129", "rs182276278", "rs6825684", "rs142678070", 
"rs182966144", "rs149368004"), Gene = c("HOXB13", "POU5F1B", 
"MSMB", "HNF1B", "MYEOV", "TH", "COPZ2", "KLK3", "AC007461.1", 
"ITGA6", "EFCAB13", "BIK", "MICU1", "SKIL", "MYO3B", "KRT8", 
"STIM2", "C5orf38", "P2RX4", "TSNARE1", "ZNF652", "FERMT2", "RARB", 
"CATSPER1", "OPRK1", "UGT1A10, UGT1A8", "EHBP1", "WNT2", "WDR11", 
"ST6GALNAC3", "LMTK2", "GPHN", "JAZF1", "SORCS2", "NEDD1", "TFCP2L1", 
"CMIP", "SNX16", "PRKCA", "MYOCD", "MAP3K7", "HGSNAT", "ANK2", 
"PRRC1", "GRM4", "CPXM2", "TET2", "ARSJ", "ADO", "VLDLR"), MAF = c(0.0018, 
0.029, 0.61, 0.6, 0.51, 0.84, 0.002, 0.073, 0.53, 0.062, 0.0028, 
0.5, 0.0013, 0.2, 0.5, 0.11, 0.0012, 0.0012, 0.0072, 0.0013, 
0.0038, 0.81, 0.0014, 0.0012, 0.0038, 0.0013, 0.54, 0.47, 0.59, 
0.0013, 0.54, 0.0014, 0.23, 0.0076, 0.0012, 0.0013, 0.0022, 0.0011, 
0.0022, 0.092, 0.0017, 0.0055, 0.001, 0.0092, 0.0021, 0.003, 
0.13, 0.00093, 0.0054, 0.0033), P = c(4.8e-41, 7.9e-20, 1.3e-16, 
1.1e-14, 1.1e-14, 3e-14, 3.4e-13, 2.9e-12, 3e-11, 2e-09, 2.2e-09, 
3.8e-09, 3.8e-09, 1.8e-08, 2e-08, 2.6e-08, 4.3e-08, 5e-08, 5.6e-08, 
5.9e-08, 7.9e-08, 8.4e-08, 9e-08, 1.1e-07, 1.4e-07, 1.6e-07, 
1.7e-07, 1.9e-07, 2e-07, 2.1e-07, 2.2e-07, 2.3e-07, 2.3e-07, 
2.4e-07, 2.5e-07, 2.6e-07, 2.7e-07, 2.7e-07, 2.9e-07, 3e-07, 
3.3e-07, 3.6e-07, 3.9e-07, 4.3e-07, 4.5e-07, 4.8e-07, 4.8e-07, 
4.9e-07, 5.7e-07, 6.4e-07)))

## ORIGINAL lncRNA2GOA DATA

# link: https://tare.medisin.ntnu.no/pred_lncRNA/

expression_vals <- '/run/media/casper/USB_ccpeters/internship_thesis/data/rezvan_lncrna2function/lncRNA2function_data.txt'
ensembl_to_go <- '/run/media/casper/USB_ccpeters/internship_thesis/data/rezvan_lncrna2function/EG2GO.txt'

# expression data containing fpkm expression values. fpkm is fragements per kilobase million. fragments means that is is for paired-end data.
options(stringsAsFactors = FALSE)
ExpressionData <- read.table(expression_vals, sep = '\t', head = TRUE)

# Gene ID's and annotation term type.
EnsemblID2GOID <- read.table(ensembl_to_go, sep = '\t', head = TRUE)

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

# CHROMOSOME POSITION FILES

# bed files / other tables are obtained from the human genome tables browser: 
# https://genome.ucsc.edu/cgi-bin/hgTables
# GENOME BED FILE: (hg19, Ensembl genes track, BED format)
# ENSEMBL TRANSCRIPTS TO GENES: (hg19, ensemblToGeneName table, Ensembl genes track)

# commented code to generate snp bed.
# tmp_var <- lapply(
#             sapply(
#               as.character(geneset_table$Position), function(x) {strsplit(x, ":")} # split chromosome and position
#               ), function(x) {c(paste0("chr",x[1]),rep(gsub(",", "", unlist(x[2], FALSE, FALSE)),2))}) # format chromosome number and replicate position 
# for (i in seq_len(length(tmp_var))) {
#   string <- tmp_var[[i]]
#   # edit string to be seperated by tabs and write to file "snp_hg19.bed"
#   write(paste0(string, collapse = "\t"), "blah.txt", sep = "\t", append = TRUE)
# }

genome_bed_file <- "/run/media/casper/USB_ccpeters/internship_thesis/data/chromosome positions/hg19.bed"
snp_bed_file <- "/run/media/casper/USB_ccpeters/internship_thesis/data/chromosome positions/snp_hg19.bed"
transcript_to_gene_file <- "/run/media/casper/USB_ccpeters/internship_thesis/data/chromosome positions/hg19_transcript_to_genename.txt"

genome_bed <- import(genome_bed_file)
snp_bed <- import(snp_bed_file)
transcript_to_gene <- read.csv(transcript_to_gene_file, sep="\t")

# FUNCTION DEFINTIONS FOR CALCULATIONS

symbol_to_ensembl <- function(geneset) {
  # convert the genenames to ensembl genes
  gene_list <- as.character(geneset)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  convertion_table <- getBM(attributes <- c("hgnc_symbol", "ensembl_gene_id"), 
                            filters = "hgnc_symbol", values = gene_list, 
                            bmHeader = TRUE, mart = mart)
  colnames(convertion_table) <- c("hgnc_symbol","ensembl_gene_id")
  convertion_table$ensembl_gene_id <- as.character(
    convertion_table$ensembl_gene_id)
  return(convertion_table)
}

parse_comma_seperated_vector <- function(gene_list_unfiltered) {
  # parse gene list by further seperating comma seperates items
  gene_list <- c()
  for (gene in gene_list_unfiltered) {
    gene_list <- c(gene_list, stringr::str_trim(unlist(strsplit(gene, ","))))
  }
  return (gene_list)
}

# enrichr and lncrna2goa comparison (all biological process)
# lncrna2goa_gos = the amount of gos to be taken from the lncrna2goa algoritm.
# 0 = ALL. The reason for this delimiter is that the lncrna2goa can generate
# a lot of gos.
enrichr_and_lncrna2goa <- function(geneset, convertion_table, 
                                   lncrna2goa_gos = 0) {
  # first run enrichr
  enrichr_result <- enrichR::enrichr(geneset, "GO_Biological_Process_2018")
  # then run lncrna2goa (single gene expression correlative enrichment analysis)
  # create empty columns
  convertion_table[["EnrichR Annotation"]] <- rep(NA, nrow(convertion_table))
  convertion_table[["lncRNA2GOA Annotation"]] <- rep(NA, nrow(convertion_table))
  
  ens_ids <- convertion_table$ensembl_gene_id
  
  gene_enrichments <- list()
  # run enrichments (for each gene)
  for (gene in ens_ids[ens_ids %in% ExpressionData$GeneID]) {
    print(gene)
    enrichment <- expression_prediction(gene, 
                                        expset, 
                                        "human", 
                                        "BP", 
                                        idtype = "ENSEMBL", 
                                        verbose = TRUE, 
                                        id_select_vector = filter_vector, 
                                        id_translation_df = id_translation_df)
    gene_enrichments[[gene]] <- enrichment
    # add the top x gos to the conversion table
    if (lncrna2goa_gos == 0) {
      convertion_table[ens_ids==gene,]$`lncRNA2GOA Annotation` <- 
        paste0(enrichment[, "GOID"], collapse=";")      
    } else {
      if (nrow(enrichment) >= lncrna2goa_gos) {
        convertion_table[ens_ids==gene,]$`lncRNA2GOA Annotation` <- 
          paste0(enrichment[1:lncrna2goa_gos, "GOID"], collapse=";")
      } else { # add everything else if the length is less than 10
        convertion_table[ens_ids==gene,]$`lncRNA2GOA Annotation` <- 
          paste0(enrichment[, "GOID"], collapse=";")
      }  
    }
    
  }
  
  # now the results need to be combined somehow. --> GO annotations for each
  # gene. This will be added to the convertion table.
  
  # transform EnrichR result first.
  
  regex_str <- "GO:\\d*"
  for (i in seq_len(nrow(enrichr_result$GO_Biological_Process_2018))) {
    row <- unlist(enrichr_result$GO_Biological_Process_2018[i,], FALSE, FALSE)
    exp <- regexec(regex_str, row[1])
    # GO needs to be seperated from description string with regex
    regex_result <- unlist(regmatches(row[1], exp), FALSE, FALSE)
    hgnc_symbols <- convertion_table$hgnc_symbol
    for (gene in unlist(strsplit(row[9], ";"), FALSE, FALSE)) {
      existing_gos <- convertion_table[hgnc_symbols==gene,]$`EnrichR Annotation`
      if (!is.na(regex_result)) {
        if (all(is.na(existing_gos))) {
          convertion_table[hgnc_symbols==gene,]$`EnrichR Annotation` <- 
            regex_result
        } else {
          convertion_table[hgnc_symbols==gene,]$`EnrichR Annotation` <- 
            paste(convertion_table[hgnc_symbols==gene,]$
                    `EnrichR Annotation`, regex_result, sep=";")
        }
      }
    }
  }
  
  # now add GO ratios
  
  # create empty column
  convertion_table[["GO ratio"]] <- rep(NA, nrow(convertion_table))
  for (i in seq_len(nrow(convertion_table))) {
    # count for each gopair how many overlap between the two analysis'
    boolvec <- unlist(strsplit(convertion_table$`EnrichR Annotation`, ";")[i], 
                    FALSE, FALSE) %in% 
      unlist(strsplit(convertion_table$`lncRNA2GOA Annotation`, ";")[i], 
             FALSE, FALSE)
    convertion_table[i, "GO ratio"] <- length(boolvec[boolvec==TRUE])
  }
  
  return(convertion_table)
}

find_closest_genes <- function(snp_bed, genome_bed, transcript_to_gene) {
  nearest_ensembl_transcripts <- as.data.frame(distanceToNearest(snp_bed, 
                                                                 genome_bed))
  # now translate all the transcripts to genes
  nearest_genes <- apply(nearest_ensembl_transcripts, 1, function(row) {
    hitno <- row[2]
    transcript <- genome_bed[hitno]$name
    gene <- transcript_to_gene[transcript_to_gene$X.name==transcript,2]
    return(gene)
    })
  return(nearest_genes)
}

# STEP 2: PREPARSE SOME OF THE DATA

geneset <- parse_comma_seperated_vector(geneset_table$Gene)
conversion_table <- symbol_to_ensembl(geneset)

# STEP 3: RUN ALREADY PREDICTED GENES

result1 <- enrichr_and_lncrna2goa(geneset, conversion_table, lncrna2goa_gos = 0)

# STEP 4: CALCULATE WITH CURRENTLY ANNOTATED CLOSESTS GENES AND COPMARE

nearest_genes <- find_closest_genes(snp_bed, genome_bed, transcript_to_gene)
View(cbind(geneset_table, nearest_genes))
# Annotated genes are different! Let's compare them now with Enrichr again.
result2 <- enrichr_and_lncrna2goa(nearest_genes, conversion_table, lncrna2goa_gos = 0)
