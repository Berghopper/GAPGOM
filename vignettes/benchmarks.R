## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_libs, warning = F, message = F, include = T--------------------
library(GAPGOM)
library(profvis)
library(GO.db)
library(graph)
library(ggplot2)
library(reshape2)

## ---- eval=FALSE---------------------------------------------------------
#  # Example with default dataset, take a look at the data documentation
#  # to fully grasp what's going on with making of the filter etc. (Biobase
#  # ExpressionSet)
#  
#  # keep everything that is a protein coding gene
#  filter_vector <- expset@featureData@data[(expset@featureData@data$GeneType=="protein_coding"),]$GeneID
#  # set gid and run.
#  gid <- "ENSG00000228630"
#  
#  p <- profvis::profvis({GAPGOM::expression_prediction(gid,
#                                          GAPGOM::expset,
#                                          "human",
#                                          "BP",
#                                          id_translation_df = GAPGOM::id_translation_df,
#                                          id_select_vector = filter_vector,
#                                          method = "combine", verbose = T, filter_pvals = T
#                         )})
#  time <- max(p$x$message$prof$time)*10
#  mem <- max(p$x$message$prof$memalloc)

## ------------------------------------------------------------------------
# laptop
310
# time
124.91
# mem
# ---
# server
# time
560
# mem
72.47

## ---- eval=FALSE---------------------------------------------------------
#  # prepare the godata for mouse and some other calculations later needed in benchmarking
#  organism <- "human"
#  ontology <- "BP"
#  go_data <- GAPGOM::set_go_data(organism, ontology)

## ----topotitj, eval=FALSE------------------------------------------------
#  # grab 15 random GOs (for term algorithm)
#  ## sample(unique(go_data@geneAnno$GO), 15)
#  random_gos <- c("GO:0030177", "GO:0001771", "GO:0045715", "GO:0044330", "GO:0098780",
#  "GO:1901006", "GO:0061143", "GO:0060025", "GO:0015695", "GO:0090074",
#  "GO:0035445", "GO:0008595", "GO:1903634", "GO:1903826", "GO:0048389"
#  )
#  # print them for reproducability
#  ## dput(random_gos)
#  # now compare all unique random GO pairs. (105 uniques).
#  unique_pairs <- GAPGOM:::.unique_combos(random_gos, random_gos)
#  
#  times <- c()
#  mem_usages <- c()
#  for (i in seq_len(nrow(unique_pairs))) {
#    prof_toptitj <- profvis({
#      pair <- unique_pairs[i]
#      go1 <- pair[[1]]
#      go2 <- pair[[2]]
#      GAPGOM::topo_ic_sim_term(organism, ontology, go1, go2, go_data = go_data)
#    })
#    time <- max(prof_toptitj$x$message$prof$time)*10
#    mem <- max(prof_toptitj$x$message$prof$memalloc)
#    mem_usages <- c(mem_usages, mem)
#    times <- c(times, time)
#    gc()
#  }
#  times_term <- times
#  mems_term <- mem_usages

## ---- fig.width=3, fig.height=3------------------------------------------
times_term_df <- data.frame(GAPGOM:::benchmarks$server_times_term, 
                          GAPGOM:::benchmarks$laptop_times_term, 
                          seq_len(105)
                          )
colnames(times_term_df) <- c("server", "laptop", "n")
times_term_df_melted <- melt(times_term_df, id="n")
colnames(times_term_df_melted) <- c("n", "machine", "milliseconds")
ggplot(times_term_df_melted, aes(x=machine, y=milliseconds, colour=machine)) + 
  geom_boxplot(notch=F) + 
  scale_y_continuous(breaks=pretty(times_term_df_melted$milliseconds, n = 5)) + 
  labs(title = paste(strwrap("Speed of term algorithm", width = 20), collapse = "\n")) + 
  theme(plot.title = element_text(hjust = 0.5))

## ---- fig.width=3, fig.height=3------------------------------------------
mems_term_df <- data.frame(GAPGOM:::benchmarks$server_mems_term, 
                          GAPGOM:::benchmarks$laptop_mems_term, 
                          seq_len(105)
                          )
colnames(mems_term_df) <- c("server", "laptop", "n")
mems_term_df_melted <- melt(times_term_df, id="n")
colnames(mems_term_df_melted) <- c("n", "machine", "RAM")
ggplot(mems_term_df_melted, aes(x=machine, y=RAM, colour=machine)) + 
  geom_boxplot(notch=F) + 
  scale_y_continuous(breaks=pretty(mems_term_df_melted$RAM, n = 5)) + 
  labs(title = paste(strwrap("RAM usage for gene algorithm", width = 20), collapse = "\n")) + 
  theme(plot.title = element_text(hjust = 0.5))

## ---- eval=FALSE---------------------------------------------------------
#  ## dput(sample(unique(go_data@geneAnno$ENTREZID), 5))
#  random_genes <- c("3848", "2824", "65108", "3988", "10800")
#  
#  unique_pairs <- GAPGOM:::.unique_combos(random_genes, random_genes)
#  
#  times <- c()
#  mem_usages <- c()
#  for (i in seq_len(nrow(unique_pairs))) {
#    prof_topg1g2 <- profvis({
#      pair <- unique_pairs[i]
#      gene1 <- pair[[1]]
#      gene2 <- pair[[2]]
#      GAPGOM::topo_ic_sim_genes(organism, ontology, gene1, gene2, go_data=go_data)
#    })
#    time <- max(prof_topg1g2$x$message$prof$time)*10
#    mem <- max(prof_topg1g2$x$message$prof$memalloc)
#    mem_usages <- c(mem_usages, mem)
#    times <- c(times, time)
#    gc()
#  }
#  times
#  mem_usages
#  times_gen <- times
#  mems_gen <- mem_usages

## ---- fig.width=3, fig.height=3------------------------------------------
times_gene_df <- data.frame(GAPGOM:::benchmarks$server_times_gen, 
                          GAPGOM:::benchmarks$laptop_times_gen, 
                          seq_len(10)
                          )
colnames(times_gene_df) <- c("server", "laptop", "n")
times_gene_df_melted <- melt(times_gene_df, id="n")
times_gene_df_melted$value <- times_gene_df_melted$value/1000 
colnames(times_gene_df_melted) <- c("n", "machine", "seconds")
ggplot(times_gene_df_melted, aes(x=machine, y=seconds, colour=machine)) + 
  geom_boxplot(notch=F) + 
  labs(title = paste(strwrap("Speed of gene algorithm", width = 20), collapse = "\n")) + 
  theme(plot.title = element_text(hjust = 0.5))

## ---- fig.width=3, fig.height=3------------------------------------------
mems_gene_df <- data.frame(GAPGOM:::benchmarks$server_mems_gen, 
                          GAPGOM:::benchmarks$laptop_mems_gen, 
                          seq_len(10)
                          )
colnames(mems_gene_df) <- c("server", "laptop", "n")
mems_gene_df_melted <- melt(mems_gene_df, id="n")
colnames(mems_gene_df_melted) <- c("n", "machine", "RAM")
ggplot(mems_gene_df_melted, aes(x=machine, y=RAM, colour=machine)) + 
  geom_boxplot(notch=F) + 
  scale_y_continuous(breaks=pretty(mems_gene_df_melted$RAM, n = 5)) + 
  labs(title = paste(strwrap("RAM usage for gene algorithm", width = 20), collapse = "\n")) + 
  theme(plot.title = element_text(hjust = 0.5))

## ---- eval=FALSE---------------------------------------------------------
#  list1=c("126133","221","218","216","8854","220","219","160428","224","222","8659","501","64577","223","217","4329","10840","7915","5832")
#  times <- c()
#  mem_usages <- c()
#  for (i in seq(length(list1)-1)) {
#    sampled_list <- list1[1:(i+1)]
#    print(sampled_list)
#    p <- profvis({
#       GAPGOM::topo_ic_sim_genes(organism, ontology, sampled_list, sampled_list, drop=NULL, go_data=go_data)
#     })
#    time <- max(p$x$message$prof$time)*10
#    mem <- max(p$x$message$prof$memalloc)
#    mem_usages <- c(mem_usages, mem)
#    times <- c(times, time)
#    gc()
#  }
#  times
#  mem_usages
#  times_genset <- times
#  mems_genset <- mem_usages

## ---- fig.width=6, fig.height=3------------------------------------------
final_times <- data.frame(GAPGOM:::benchmarks$server_times_genset, 
                          GAPGOM:::benchmarks$laptop_times_genset, 
                          seq_along(GAPGOM:::benchmarks$laptop_times_genset)
                          )
colnames(final_times) <- c("server", "laptop", "n")
final_times_melted <- melt(final_times, id="n")
final_times_melted$value <- final_times_melted$value/1000/60
colnames(final_times_melted) <- c("Genes in geneset", "Machine", "Minutes")
p <- ggplot(data=final_times_melted, aes(x=`Genes in geneset`, y=Minutes, colour=Machine)) + 
  geom_point() + 
  labs(title = paste(strwrap("Speed of TopoICSim geneset algorithm", width = 20), collapse = "\n")) + 
  theme(plot.title = element_text(hjust = 0.5))
p
p + stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

## ---- fig.width=3, fig.height=3------------------------------------------
mems_geneset_df <- data.frame(GAPGOM:::benchmarks$server_mems_genset, 
                          GAPGOM:::benchmarks$laptop_mems_genset, 
                          seq_len(18)
                          )
colnames(mems_geneset_df) <- c("server", "laptop", "n")
mems_geneset_df_melted <- melt(mems_geneset_df, id="n")
colnames(mems_geneset_df_melted) <- c("Genes in geneset", "machine", "RAM")
ggplot(mems_geneset_df_melted, aes(x=machine, y=RAM, colour=machine)) + 
  geom_boxplot(notch=F) + 
  scale_y_continuous(breaks=pretty(mems_geneset_df_melted$RAM, n = 5)) + 
  labs(title = paste(strwrap("RAM usage for geneset algorithm", width = 20), collapse = "\n")) + 
  theme(plot.title = element_text(hjust = 0.5))

## ---- fig.width=6, fig.height=3------------------------------------------
p <- ggplot(mems_geneset_df_melted, aes(x=`Genes in geneset`, y=RAM, colour=machine)) + 
  geom_point()  + 
  scale_y_continuous(breaks=pretty(mems_geneset_df_melted$RAM, n = 5)) + 
  labs(title = paste(strwrap("RAM usage for geneset algorithm", width = 20), collapse = "\n")) + 
  theme(plot.title = element_text(hjust = 0.5))
p
p + stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

## ---- eval=F-------------------------------------------------------------
#  ## set stringsasfactors to false to ensure data is loaded properly.
#  options(stringsAsFactors = F)
#  
#  ## load data
#  #
#  expdata <- read.table("~/Downloads/GSE63733_m.txt")
#  # http://software.broadinstitute.org/gsea/msigdb/cards/HALLMARK_GLYCOLYSIS.html
#  # GLYCOLOSIS HALLMARK
#  geneset <- read.table("~/Downloads/geneset_glyc.txt", sep="\n")
#  
#  colnames(expdata)[1:2] <- c("ENSEMBL", "SYMBOL") # Set important colnames
#  
#  ## make an expressionset
#  
#  expression_matrix <- as.matrix(expdata[,3:ncol(expdata)])
#  rownames(expression_matrix) <- expdata[[1]]
#  featuredat <- as.data.frame(expdata[,1:2])
#  rownames(featuredat) <- expdata[[1]]
#  expset <- ExpressionSet(expression_matrix,
#                          featureData = new("AnnotatedDataFrame",
#                          data=featuredat))
#  
#  ## make a selection on the geneset
#  geneset <- as.vector(geneset[3:nrow(geneset),])
#  geneset_ensembl <- expdata[expdata[[2]] %in% geneset,][[1]]
#  # select on the 10 genes with the most variance
#  geneset_ensembl_topvar <- names(rev(sort(apply(
#    expression_matrix[geneset_ensembl,], 1, var)))[1:10])
#  # convert back to symbol
#  geneset_selection <- expdata[expdata$ENSEMBL %in% geneset_ensembl_topvar,][[2]]
#  
#  ## predict annotation of the geneset selection (per ontology).
#  
#  bench_results <- list()
#  for (ontology in c("MF", "BP", "CC")) {
#    predicted_annotations <- list()
#    for (gene in geneset_ensembl_topvar) {
#      symbol <- expdata[expdata$ENSEMBL==gene,]$SYMBOL
#      new_gene_name <- paste0(symbol, "_pred")
#  
#      result <- GAPGOM::expression_prediction(gene,
#                                            expset,
#                                            "human",
#                                            ontology,
#                                            idtype = "ENSEMBL")
#      go_annotation_prediction <- unique(as.vector(result$GOID))
#      predicted_annotations[[new_gene_name]] <- go_annotation_prediction
#    }
#  
#    ## calculate original set
#    orig <- GAPGOM::topo_ic_sim_genes("human", ontology, geneset_selection,
#                                      geneset_selection, idtype = "SYMBOL")
#    ## compare with custom now
#    cus <- GAPGOM::topo_ic_sim_genes("human", ontology, c(), geneset_selection,
#                                     custom_genes1 = predicted_annotations,
#                                     idtype = "SYMBOL",
#                                     all_go_pairs = orig$AllGoPairs) # use original
#                                     # GO matrix, as similar GOs are likely
#                                     # to occur.
#    ## append results to list
#    bench_results[[ontology]] <- list(orig=orig, cus=cus)
#  }

## ------------------------------------------------------------------------
sessionInfo()

