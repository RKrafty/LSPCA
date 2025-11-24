## code to prepare `mydataset` dataset goes here
HC <- load("data-raw/HC.rda")
FEP <- load("data-raw/FEP.rda")
usethis::use_data(HC, FEP, overwrite = TRUE)
