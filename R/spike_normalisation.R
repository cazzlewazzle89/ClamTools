#' Import Spike Counts from idxstats Output
#'
#' Extracts the sample identifier from an idxstats file path, reads the tab-delimited 
#' count data, filters out unmapped or unplaced reads, and appends the sample ID.
#'
#' @param x Character string. The file path to a spike count text file (e.g., idxstats output).
#'
#' @return A tibble with five columns:
#' \describe{
#'   \item{Species}{Character. Name of the spiked-in reference species.}
#'   \item{Length}{Numeric. Sequence length of the reference genome.}
#'   \item{Count}{Numeric. Number of mapped reads.}
#'   \item{Unmapped}{Numeric. Number of unmapped reads.}
#'   \item{Sample}{Character. Extracted sample identifier.}
#' }
#' 
#' @importFrom stringr str_remove_all str_remove str_starts
#' @importFrom readr read_tsv
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
#' @export
import_spike_counts <- function(x) {

  sampleid <- stringr::str_remove_all(x, '.*\\/')
  sampleid <- stringr::str_remove(sampleid, '_spike.txt')

  readr::read_tsv(x, col_names = c('Species', 'Length', 'Count', 'Unmapped'), show_col_types = FALSE) %>%
    dplyr::filter(!stringr::str_starts(Species, '\\*')) %>%
    dplyr::mutate(Sample = sampleid)

}

#' Calculate Abundance Scaling Factors for All Samples
#'
#' Takes a complete data frame of spike counts, splits it by sample, and applies
#' robust linear regression to calculate scaling factors for every sample automatically.
#'
#' @param spike_counts A data frame containing columns \code{Sample}, \code{Count}, \code{Length}, and \code{Species}.
#' @param spike_expected A data frame containing expected spike loads (\code{Species} and \code{Expected_Cells}).
#'
#' @return A bundled tibble with columns \code{Sample}, \code{ScalingFactor}, and \code{Method} for all samples.
#' @importFrom dplyr left_join mutate filter tibble
#' @importFrom purrr map list_rbind
#' @importFrom MASS rlm
#' @importFrom stats coef
#' @export
calculate_scaling_factors <- function(spike_counts, spike_expected) {
  
  # split the data frame into a list of smaller data frames, one per Sample
  sample_groups <- split(spike_counts, spike_counts$Sample)
  
  # map over each sub-data frame and calculate the scaling factor
  purrr::map(sample_groups, \(sub_df) {
    
    current_sample <- unique(sub_df$Sample)
    
    # join and normalise data for this sample
    df <- sub_df %>% 
      dplyr::mutate(ReadsNormalised = Count / Length) %>%
      dplyr::left_join(spike_expected, by = 'Species')
    
    # conditional logic pathway
    if (nrow(df) < 1) {
      return(dplyr::tibble(Sample = current_sample, ScalingFactor = 0, Method = "Failed: No spikes detected"))
    }
    
    if (nrow(df) == 1) {
      temp_factor <- df$Expected_Cells / df$ReadsNormalised
      return(dplyr::tibble(Sample = current_sample, ScalingFactor = temp_factor, Method = "Manual: Single spike"))
    }
    
    # robust linear regression
    tryCatch({
      model <- MASS::rlm(Expected_Cells ~ ReadsNormalised + 0, data = df)
      return(dplyr::tibble(
        Sample = current_sample,
        ScalingFactor = stats::coef(model)[["ReadsNormalised"]], 
        Method = "RLM: Robust Regression"
      ))
    }, error = function(e) {
      return(dplyr::tibble(Sample = current_sample, ScalingFactor = 0, Method = "Failed: RLM Error"))
    })
  }) %>% 
    # bind the list of results back into a single clean data frame
    purrr::list_rbind()
}

