## code to prepare `spike_expected_high` dataset goes here
spike_expected_high <- readxl::read_xlsx('spike.xlsx', sheet = 'HighLoad')
usethis::use_data(spike_expected_high, overwrite = TRUE)
