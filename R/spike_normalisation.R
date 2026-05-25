
# import counts for each spiked species in the form of idxstats output
import_spike_counts <- function(x) {

  sampleid <- str_remove_all(x, '.*\\/')
  sampleid <- str_remove(sampleid, '_spike.txt')

  read_tsv(x, col_names = c('Species', 'Length', 'Count', 'Unmapped')) %>%
    filter(! str_starts(Species, '\\*')) %>%
    mutate(Sample = sampleid)

}

# make like of spike count files
sample_list <- list.files(
  path = '~/Desktop/TEMP/ClamToolsTest/', 
  pattern = '.*_spike.txt$',
  full.names = T
)

# import and process spike counts
spike_counts <- map_dfr(sample_list, import_spike_counts) %>%
  mutate(Species = str_replace(Species, '_', ' '))

# import expected cell count for spike
spike_expected_low <- readxl::read_xlsx('spike.xlsx', sheet = 'LowLoad')
spike_expected_high <- readxl::read_xlsx('spike.xlsx', sheet = 'HighLoad')

# import bracken abundance table
bracken <- read_tsv('~/Desktop/TEMP/ClamToolsTest/bracken_combined.tsv') %>%
  dplyr::select(name, taxonomy_id, ends_with('_num'))

# import summary stats of genome length for each species in database
taxid_size <- read_tsv('taxid_genomelen_lineage.txt')

# combine taxid with species genome size summary stats
name_taxid_size <- bracken %>% 
  dplyr::select(name, taxonomy_id) %>%
  left_join(taxid_size, by = c('taxonomy_id' = 'species_taxid')) %>%
  dplyr::select(-name) %>%
  rename(Lineage = 5)

# combine with bracken microbiome profiles
bracken_long_size <- bracken %>%
  pivot_longer(
    cols = -c(name, taxonomy_id),
    names_to = 'Sample',
    values_to = 'Count'
  ) %>%
  mutate(Sample = str_remove(Sample, '_brackenout.tsv_num')) %>%
  left_join(name_taxid_size, by = 'taxonomy_id')

# calculate scaling factor for abundance normalisation
calculate_scaling_factor <- function(sample, spike_expected) {
  
  # subset relevant sample and normalised count by genome size
  df <- spike_counts %>% filter(Sample == sample) %>%
    mutate(ReadsNormalised = Count / Length) %>%
    left_join(spike_expected, by = 'Species')
  
  # Logic to handle missing/low spike-ins
  if (nrow(df) < 1) {
    return(tibble(Sample = sample,
                  ScalingFactor = 0, 
                  Method = "Failed: No spikes detected"))
  }
  
  if (nrow(df) == 1) {
    # Simple ratio if only one spike exists (cannot do regression)
    temp_factor <- df$Expected_Cells / df$ReadsNormalised
    return(tibble(Sample = sample,
                  ScalingFactor = temp_factor, 
                  Method = "Manual: Single spike"))
  }
  
  # ROBUST LINEAR REGRESSION
  # Model: Known_Cells ~ (Reads / Genome_Size)
  # + 0 forces the intercept through the origin (0 reads must = 0 cells)
  tryCatch({
    model <- rlm(Expected_Cells ~ ReadsNormalised + 0, data = df)
    
    return(tibble(Sample = sample,
                  ScalingFactor = coef(model)[["ReadsNormalised"]], 
                  Method = "RLM: Robust Regression"))
  }, error = function(e) {
    return(tibble(Sample = sample, 
                  ScalingFactor = 0, 
                  Method = "Failed: RLM Error"))
  })
}

# calc scaling factor for each sample
scaling_factors <- map(
  unique(spike_counts$Sample),
  ~ calculate_scaling_factor(.x, spike_expected_low)
) %>%
  list_rbind()

# combine scaling factor data with bracken abundance table
# and calculate absolute abundance
bracken_abundance_temp <- bracken_long_size %>%
    left_join(scaling_factors, by = 'Sample') %>%
    mutate(ReadsPerBp = Count / `genome_size:median`,
           AbsoluteCells = ReadsPerBp * ScalingFactor) %>%
    dplyr::select(Sample, name, Count, ScalingFactor, AbsoluteCells) %>%
    filter(!is.na(AbsoluteCells)) # remove species with no genome size in database

