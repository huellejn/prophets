#' Analysis of PFSr using count-based method
#'
#' @param data a data frame containing columns 'PFS1' (numeric), 'PFS2' (numeric), and 'status' (0/1).
#' @param delta a numeric value above 0. Desired PFSr (PFS2/PFS1).
#' @param prob a numeric value. Desired probability.
#' @param modified a boolean value indicating if the PFS2 values should be modified.
#' @param min_pfs2 a numeric value. The minimal PFS2 that should be modified measured in months.
#' @param tidy a boolean value indicating if the output should be formatted as a tibble.
#'
#' @return
#' @export
#' @importFrom dplyr filter mutate
#' @importFrom broom tidy
#'
#' @examples
count_PFSr <- function(data,
                       delta = 1,
                       prob = 0.5,
                       modified = FALSE, min_pfs2 = 6,
                       tidy = FALSE
                       ) {

  # Check that all required columns are present
  required_columns <- c("PFS1", "PFS2", "status")

  if(! all(required_columns %in% colnames(data)) ) {
    missing_columns <- required_columns[!required_columns %in% colnames(data)]
    stop(paste0("Error: The following columns are required, but missing: ", paste(missing_columns, collapse = ", ")))
  }

  if(modified == T){
    data <- data %>%
      dplyr::mutate(
        PFS1  = pmax(PFS1, 2.0),
        PFS2  = ifelse(PFS2 > min_pfs2 & ratio < delta, PFS1*delta+0.25, PFS2)
      )
  }

  data <- data %>%
    dplyr::mutate(
      ratio = round(PFS2/PFS1, 3)
    ) %>%
    dplyr::filter(status == 1 | (status == 0 & ratio >= delta))

  pos <- sum(data$ratio >= delta)
  neg <- sum(data$ratio < delta)
  res <- binom.test(pos, (neg+pos), p = prob)

  if(tidy == T) {
    res <- broom::tidy(stats) %>%
      mutate(
        PFSr_method = "count-based",
        delta = delta
        )
  }

  return(res)

}
