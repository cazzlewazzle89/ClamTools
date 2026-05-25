# ClamTools
### Miscellaneous functions useful for analysing microbiome abundance data

#### Installation
```R
install.packages("pak")
pak::pkg_install("cazzlewazzle89/ClamTools")
```

#### Example workflow
```R
# create list of spike count files
sample_list <- list.files(
  path = "SPIKE_OUT/",
  pattern = '.*_spike.txt$',
  full.names = T
)

# import and process spike counts
spike_counts <- map(sample_list, import_spike_counts) %>%
  list_rbind() %>%
  mutate(Species = str_replace(Species, '_', ' '))

# calculate the scaling factor for each sample
scaling_factors <- calculate_scaling_factors(spike_counts, spike_expected_low)

# import bracken abundance table
bracken <- read_tsv('bracken_combined.tsv') %>%
  dplyr::select(name, taxonomy_id, ends_with('_num'))

# combine with bracken microbiome profiles
bracken_long_size <- bracken %>%
  pivot_longer(
    cols = -c(name, taxonomy_id),
    names_to = 'Sample',
    values_to = 'Count'
  ) %>%
  mutate(Sample = str_remove(Sample, '_out.txt_num')) %>%
  left_join(taxid_size, by = c('taxonomy_id' = 'species_taxid'))

# combine scaling factor data with bracken abundance table
# and calculate absolute abundance
bracken_abundance <- bracken_long_size %>%
  left_join(scaling_factors, by = 'Sample') %>%
  mutate(ReadsPerBp = Count / `genome_size:median`, ##normalising read counts by species genome size
         AbsoluteCells = ReadsPerBp * ScalingFactor) %>% ## calaulating the absolute number of species genome copies using ReadsPerBp ReadsPerBp* ScalingFactor
  dplyr::select(Sample, name, Count, ScalingFactor, AbsoluteCells) %>%
  filter(!is.na(AbsoluteCells)) # remove species with no genome size in database
```