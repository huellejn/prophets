#' Check if all required columns are necessary
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )
#' @param required_columns vector of expected column names
#' @noRd
check_columns <- function(data, required_columns) {
  present_columns <- names(data)
  missing_columns <- required_columns[! required_columns %in% present_columns]
  if(length(missing_columns) > 0) stop( paste("The following required columns are missing from the data:", missing_columns))
}


#' Check if data has the appropiate format
#' @param 
#' @noRd
check_data <- function(value, type_expected, range) {
  
  # expected types
  expected_type <- list(
    alpha = "numeric",
    delta = "numeric",
    prob = "numeric",
    plot = "logical",
    min_pfs2 = "numeric",
    null_HR = "numeric",
    alt_HR = "numeric",
    rho = "numeric",
    model = "vector",
    selected_PFS = "vector",
    verbose = "logical",
    ges = "numeric",
    sample_size = "integer",
    power = "numeric",
    lost = "numeric",
    pfsratio = "numeric",
    p0 = "numeric",
    p1 = "numeric",
    beta = "numeric",
    k = "numeric",
    tidy = "logical"
  )
  
  # expected ranges
  expected_type <- list(
    alpha = c(min = 0, max = 1),
    delta = c(min = 0, max = Inf),
    prob = c(min = 0, max = 1),
    min_pfs2 = c(min = 0, max = 1),
    null_HR = c(min = 0, max = 1),
    alt_HR = c(min = 0, max = 1),
    rho = c(min = 0, max = 1),
    model = "vector",
    selected_PFS = "vector",
    ges = c(min = 0, max = 1),
    sample_size = c(min = 0, max = Inf),
    power = c(min = 0, max = 1),
    lost = c(min = 0, max = 1),
    pfsratio = c(min = 0, max = 1),
    p0 = c(min = 0, max = 1),
    p1 = c(min = 0, max = 1),
    beta = c(min = 0, max = 1),
    k = c(min = 0, max = 1)
  )
  
}

#' Melt data
#' @param data Data frame with columns PFS1, PFS2 and status.
#' @noRd
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