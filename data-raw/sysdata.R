# sysdata.rda generation file (included params for tests)

# variable requests from package for resaving/updating
# Part of actual data generation are commented
freq_go_pairs <- GAPGOM:::freq_go_pairs # see freq_go_pairs.R
lncrnapred_baseresult <- GAPGOM:::lncrnapred_baseresult # documented example was used for the result
maintopo_baseresult <- GAPGOM:::maintopo_baseresult # documented example was used for the result
genetopo_baseresult <- GAPGOM:::genetopo_baseresult # documented example was used for the result
termtopo_baseresult <- GAPGOM:::termtopo_baseresult # documented example was used for the result
benchmarks <- GAPGOM:::benchmarks

#####
#' Benchmark prep, this is ran on multiple machine of choice and concatted into a list at the en

#' prepare

library(GAPGOM)
library(profvis)
library(GO.db)
library(graph)

# prepare the godata for mouse and some other calculations later needed in benchmarking
organism <- "human"
ontology <- "BP"
go_data <- GAPGOM::set_go_data(organism, ontology)

#' term
# grab 15 random GOs (for term algorithm)
## sample(unique(go_data@geneAnno$GO), 15)
random_gos <- c("GO:0030177", "GO:0001771", "GO:0045715", "GO:0044330", "GO:0098780",
                "GO:1901006", "GO:0061143", "GO:0060025", "GO:0015695", "GO:0090074",
                "GO:0035445", "GO:0008595", "GO:1903634", "GO:1903826", "GO:0048389"
)
# print them for reproducability
## dput(random_gos)
# now compare all unique random GO pairs. (105 uniques).
unique_pairs <- GAPGOM:::.unique_combos(random_gos, random_gos)

times <- c()
mem_usages <- c()
for (i in seq_len(nrow(unique_pairs))) {
  prof_toptitj <- profvis({
    pair <- unique_pairs[i]
    go1 <- pair[[1]]
    go2 <- pair[[2]]
    GAPGOM::topo_ic_sim_term(ontology, organism, go1, go2, go_data = go_data)
  })
  time <- max(prof_toptitj$x$message$prof$time)*10
  mem <- max(prof_toptitj$x$message$prof$memalloc)
  mem_usages <- c(mem_usages, mem)
  times <- c(times, time)
  gc()
}
times_term <- times
mems_term <- mem_usages

#' gene

## dput(sample(unique(go_data@geneAnno$ENTREZID), 5))
random_genes <- c("3848", "2824", "65108", "3988", "10800")

unique_pairs <- GAPGOM:::.unique_combos(random_genes, random_genes)

times <- c()
mem_usages <- c()
for (i in seq_len(nrow(unique_pairs))) {
  prof_topg1g2 <- profvis({
    pair <- unique_pairs[i]
    gene1 <- pair[[1]]
    gene2 <- pair[[2]]
    GAPGOM::topo_ic_sim_genes(ontology, organism, gene1, gene2, go_data=go_data)
  })
  time <- max(prof_topg1g2$x$message$prof$time)*10
  mem <- max(prof_topg1g2$x$message$prof$memalloc)
  mem_usages <- c(mem_usages, mem)
  times <- c(times, time)
  gc()
}
times
mem_usages
times_gen <- times
mems_gen <- mem_usages

#' geneset

list1=c("126133","221","218","216","8854","220","219","160428","224","222","8659","501","64577","223","217","4329","10840","7915", "5832")
times <- c()
mem_usages <- c()
for (i in seq(length(list1)-1)) {
  sampled_list <- list1[1:(i+1)]
  print(sampled_list)
  p <- profvis({
    GAPGOM::topo_ic_sim_genes(ontology, organism, sampled_list, sampled_list, drop=NULL, go_data=go_data)
  })
  time <- max(p$x$message$prof$time)*10
  mem <- max(p$x$message$prof$memalloc)
  mem_usages <- c(mem_usages, mem)
  times <- c(times, time)
  gc()
}
times
mem_usages
times_genset <- times
mems_genset <- mem_usages

#' combine to list (seperateley done per machine)

benchmarks <- list()

benchmarks$server_times_term <- times_term
benchmarks$server_times_gen <- times_gen
benchmarks$server_times_genset <- times_genset
benchmarks$server_mems_term <- mems_term
benchmarks$server_mems_gen <- mems_gen
benchmarks$server_mems_genset <- mems_genset
##
benchmarks$laptop_times_term <- times_term
benchmarks$laptop_times_gen <- times_gen
benchmarks$laptop_times_genset <- times_genset
benchmarks$laptop_mems_term <- mems_term
benchmarks$laptop_mems_gen <- mems_gen
benchmarks$laptop_mems_genset <- mems_genset

c(bench1, bench2) #... combines all machine results

#####

save(freq_go_pairs, lncrnapred_baseresult, maintopo_baseresult, genetopo_baseresult, termtopo_baseresult, benchmarks, file = "/media/casper/USB_ccpeters/Repositories/Personal/gapgom/R/sysdata.rda", compress = "xz", compression_level = 9)
