# sysdata.rda generation file

freq_go_pairs <- GAPGOM:::freq_go_pairs
lncrnapred_baseresult <- GAPGOM:::lncrnapred_baseresult
maintopo_baseresult <- GAPGOM:::maintopo_baseresult
genetopo_baseresult <- GAPGOM:::genetopo_baseresult
benchmarks <- GAPGOM:::benchmarks
save(freq_go_pairs, lncrnapred_baseresult, maintopo_baseresult, genetopo_baseresult, benchmarks, file = "/media/casper/USB_ccpeters/Repositories/Personal/gapgom/R/sysdata.rda", compress = "xz", compression_level = 9)
