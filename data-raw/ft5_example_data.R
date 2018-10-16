# In this file internal procedures will be defined for generating data, this is for "developers only".
# ft5_example_data data generation script

library(GAPGOM)
#ft5_raw <- fantom_load_raw("/media/casper/USB_ccpeters/internship_thesis/data/f5/mouse/mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt", verbose = T)
ft5_raw <- fantom_load_raw(fantom_download("~/", organism = "mouse"), verbose = T)

ft5_raw$df <- ft5_raw$df[1:10000, 1:15]
ft5_example_data <- fantom_to_expset(fanraw = ft5_raw, verbose = T) 
