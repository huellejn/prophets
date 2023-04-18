#' Analysis of progression free survival ratio (PFSr) using the midrank method
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric)
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#'
#' @return midrank-based statistics for PFSr analysis
#' @export
#' @importFrom tibble tibble
#'
#' @examples
#' data(input)
#' midrank_PFSr(data = input, delta = 1.3)
midrank_PFSr <- function(data,  
                         delta = 1 ) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("PFS1", "PFS2"))
  
  rankmax <- Vectorize(function(a, vec)  rev(which(vec <= a))[1], "a")
  rankmin <- Vectorize(function(a, vec)  which(vec >= a)[1], "a")
  
  data$PFS1m <- delta*data$PFS1
  dataMOD <- data[, c('PFS1m', 'PFS2')]
  dataMOD$PFS2[data$status == 0] <- Inf
  
  lsort = sort(unlist(data[, c('PFS1m', 'PFS2')]))
  rsort = sort(unlist(dataMOD))
  mid1 = (rankmin(dataMOD$PFS1m, rsort) + rankmax(dataMOD$PFS1m, lsort)) /2
  mid2 = (rankmin(data$PFS2, rsort) + rankmax(dataMOD$PFS2, lsort)) /2
  prop = as.numeric((mid2 - mid1) >= 0)
  estim <- mean(prop)
  sd <- sqrt((estim*(1-estim))/nrow(data))
  
  res <- tibble::tibble(
    method = "midrank",
    delta = delta,
    estimate = estim,
    conf.low = estim + qnorm(.025) * sd,
    conf.high = estim + qnorm(.975) * sd
  )
  
  return(res)
  
  }