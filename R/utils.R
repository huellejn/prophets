#' @title Utility functions
#' @section Functions:
#' \describe{
#'     \item{check_columns}{Checks if all required columns are present in the input data.}
#'     \item{melt_data}{Transfers the data to a long format}
#' }
#' @importFrom dplyr mutate select
#' @importFrom tidyr pivot_longer


# Check data

## Check if required columns are present.
check_columns <- function(data, required_columns) {
  present_columns <- names(data)
  missing_columns <- required_columns[! required_columns %in% present_columns]
  if(length(missing_columns) > 0) stop( paste("The following required columns are missing from the data:", missing_columns))
}

# Melt data

melt_data <- function(data) {
  
  data <- data %>%
    dplyr::mutate(
      status1 = 1,
      idx = 1:nrow(.)
    ) %>%
    dplyr::select(idx, PFS1, status1, PFS2, status2 = status) %>%
    tidyr::pivot_longer(cols = c("PFS1", "PFS2", "status1", "status2"), 
                 names_to = c(".value", "ct.line"), 
                 names_pattern = "^([A-Za-z]+)(\\d+)"
    ) 
  
  return(data)
  
}