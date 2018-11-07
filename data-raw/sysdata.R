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






# tmp
a<-list(laptop_times_term = c(1020, 1020, 1510, 1100, 1110, 1160, 
                              1000, 1070, 960, 930, 1050, 1090, 940, 970, 940, 980, 960, 970, 
                              1020, 1070, 960, 1140, 960, 950, 960, 960, 1010, 1020, 1220, 
                              1190, 1160, 1100, 1190, 1480, 960, 940, 960, 1060, 1340, 1230, 
                              1720, 1370, 990, 1580, 1040, 1140, 1260, 1140, 1140, 1030, 1140, 
                              1150, 1000, 980, 980, 950, 1060, 1000, 980, 1300, 1330, 1140, 
                              1480, 1280, 1480, 1500, 1620, 1590, 1030, 1040, 1350, 1070, 1180, 
                              1470, 1470, 1180, 1120, 980, 1160, 1040, 940, 1000, 1020, 1010, 
                              1040, 1090, 1090, 1160, 1130, 1000, 1020, 990, 1020, 1000, 1010, 
                              1060, 990, 1140, 990, 1110, 960, 980, 1000, 1070, 1000), laptop_times_gen = c(43960, 
                                                                                                            55740, 31120, 3460, 43640, 28020, 3250, 32750, 3720, 2340), laptop_times_genset = c(2550, 
                                                                                                                                                                                                8620, 34400, 225930, 526920, 533160, 643340, 793870, 769890, 
                                                                                                                                                                                                965980, 1106340, 1273660, 1199040, 1163100, 1431080, 1545920, 
                                                                                                                                                                                                2274760, 2748600), laptop_mems_term = c(158.06298828125, 158.032371520996, 
                                                                                                                                                                                                                                        158.033241271973, 158.068199157715, 158.039283752441, 158.043952941895, 
                                                                                                                                                                                                                                        158.042121887207, 158.032447814941, 158.038566589355, 158.029670715332, 
                                                                                                                                                                                                                                        158.027336120605, 158.039558410645, 158.040382385254, 158.028800964355, 
                                                                                                                                                                                                                                        157.374839782715, 157.376609802246, 157.377494812012, 157.375633239746, 
                                                                                                                                                                                                                                        158.160850524902, 157.377342224121, 157.387321472168, 157.376853942871, 
                                                                                                                                                                                                                                        157.39054107666, 157.375450134277, 157.37621307373, 157.375953674316, 
                                                                                                                                                                                                                                        157.376243591309, 158.013893127441, 158.015266418457, 158.034675598145, 
                                                                                                                                                                                                                                        158.025215148926, 158.019554138184, 158.022621154785, 158.025276184082, 
                                                                                                                                                                                                                                        158.044746398926, 158.010520935059, 158.009208679199, 158.010688781738, 
                                                                                                                                                                                                                                        158.018852233887, 157.966667175293, 158.103630065918, 158.162582397461, 
                                                                                                                                                                                                                                        157.96019744873, 158.164497375488, 157.988624572754, 158.165061950684, 
                                                                                                                                                                                                                                        158.172355651855, 158.151550292969, 157.949424743652, 157.95189666748, 
                                                                                                                                                                                                                                        157.739776611328, 158.012924194336, 157.743774414062, 157.739013671875, 
                                                                                                                                                                                                                                        157.736129760742, 157.73616027832, 157.732711791992, 157.742889404297, 
                                                                                                                                                                                                                                        157.737548828125, 157.737152099609, 158.129653930664, 158.074546813965, 
                                                                                                                                                                                                                                        158.172897338867, 158.08911895752, 158.068656921387, 158.086372375488, 
                                                                                                                                                                                                                                        158.091957092285, 158.09188079834, 158.091346740723, 158.161552429199, 
                                                                                                                                                                                                                                        158.163734436035, 158.092735290527, 158.162986755371, 158.173652648926, 
                                                                                                                                                                                                                                        157.983642578125, 158.057540893555, 157.981201171875, 157.953758239746, 
                                                                                                                                                                                                                                        157.946907043457, 157.963592529297, 157.948348999023, 157.945610046387, 
                                                                                                                                                                                                                                        157.948272705078, 157.950103759766, 157.737846374512, 157.739723205566, 
                                                                                                                                                                                                                                        157.744514465332, 157.747383117676, 157.751533508301, 157.749015808105, 
                                                                                                                                                                                                                                        157.37914276123, 157.380012512207, 157.380378723145, 157.380653381348, 
                                                                                                                                                                                                                                        157.381828308105, 158.047439575195, 158.054092407227, 158.049194335938, 
                                                                                                                                                                                                                                        158.063262939453, 158.163208007812, 158.079330444336, 158.164794921875, 
                                                                                                                                                                                                                                        158.164573669434, 158.164909362793, 158.087265014648), laptop_mems_gen = c(158.173652648926, 
                                                                                                                                                                                                                                                                                                                   158.17366027832, 158.17366027832, 158.108894348145, 158.17366027832, 
                                                                                                                                                                                                                                                                                                                   158.17366027832, 158.172836303711, 158.173652648926, 158.160057067871, 
                                                                                                                                                                                                                                                                                                                   158.173652648926), laptop_mems_genset = c(158.171363830566, 158.17366027832, 
                                                                                                                                                                                                                                                                                                                                                             158.17366027832, 158.17366027832, 190.418739318848, 333.050483703613, 
                                                                                                                                                                                                                                                                                                                                                             400.270927429199, 400.270927429199, 480.935447692871, 462.186309814453, 
                                                                                                                                                                                                                                                                                                                                                             539.966064453125, 555.201210021973, 574.89380645752, 579.264678955078, 
                                                                                                                                                                                                                                                                                                                                                             572.894660949707, 580.388862609863, 600.780853271484, 669.066604614258
                                                                                                                                                                                                                                                                                                                   ))

