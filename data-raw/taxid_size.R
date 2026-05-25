## code to prepare `taxid_size` dataset goes here
taxid_size <- read_tsv('taxid_genomelen_lineage.txt')
usethis::use_data(taxid_size, overwrite = TRUE)
