#' Analysis of PFSr using the non-parametric (midrank) method
#'
#' @param data
#' @param delta
#' @param tidy
#' @param modified
#' @param min_pfs2
#'
#' @return
#' @export
#'
#' @examples
midrank_PFSr <- function(data,  delta=1, tidy=T,
                         modified=F, min_pfs2=6
){

  x <- c("PFS1", "PFS2", "status")
  if(all(x %in% colnames(data))==F)
    stop("Error: some required column is missing (PFS1, PFS2 or status)")

  if(modified==T){
    data = data %>%
      mutate(
        PFS1  = pmax(PFS1,2.0),
        ratio = round(PFS2/PFS1,3),
        PFS2  = ifelse(PFS2>min_pfs2 & ratio < delta, PFS1*delta+0.25, PFS2),
        ratio = round(PFS2/PFS1,3)
      )} else{
        data = data
      }
  rankmax <- Vectorize(function(a, vec)  rev(which(vec <= a))[1], "a")
  rankmin <- Vectorize(function(a, vec)  which(vec >= a)[1], "a")
  #data <- ana_s1
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
  res <- tibble("Method" = "Midrank",
                "estimate"=estim,
                "delta"=delta,
                "conf.low" = estim + qnorm(.025) * sd,
                "conf.high" = estim + qnorm(.975) *sd
  )

  res <- res %>%
    #mutate(delta=delta) %>%
    relocate(delta, .before = estimate) %>%
    relocate(Method, .before = delta)
  if(tidy==T){
    return(res)
  }else{
    res <- res %>%
      gt(caption = "Midrank-based method PFSratio result") %>%
      fmt_number(columns = 2:5,decimals = 2)
    return(res)
  }
}
