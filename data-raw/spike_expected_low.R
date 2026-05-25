## code to prepare `spike_expected_low` dataset goes here
spike_expected_low <- readxl::read_xlsx('spike.xlsx', sheet = 'LowLoad')
usethis::use_data(spike_expected_low, overwrite = TRUE)
