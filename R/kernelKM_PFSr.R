#' Kernel conditional Kaplan Meier method
#'
#' @param data a data frame with values for PFS1/PFS2 ratio, status (1/0) and PFS1 (time0)
#' @param ratio.name name of the column that contains the PFS1/PFS2 ratio ('numeric')
#' @param status.name name of the column that contains the status ('boolean')
#' @param time0.name name of the column that contains the PFS1 ('numeric')
#' @param points 
#' @param conf.int a boolean indicating if confidende intervals should be calculated
#' @param n.boot number of permutations
#'
#' @return a Kaplan-Meier based statistics according to the kernel conditional KM method
#' @export
#'
#' @examples
#' # Example with default values:
#' kernelKM_PFSr(data, "ratio", "status", "PFS1")
#' # Example with modified values
#' kernelKM_PFSr(data, "ratio", "status", "PFS1", points = 2, conf.int = TRUE)
kernelKM_PFSr <- function(data, ratio.name = "ratio", status.name = "status", time0.name = "PFS1", points = NULL, conf.int = FALSE, n.boot = 2000) {
  
  result = estimator(data, ratio.name, status.name, time0.name)
  ratio = result$ratio
  surv = result$surv
  
  if( min(surv) > 0.5 ) {
    medsurv = NA
  } else {
    medsurv = ratio[sum(surv>0.5)+1]
  }
  
  res = list(ratio=ratio, surv=surv, med=medsurv)
  
  # If conf.int is TRUE, calculate confidence intervals
  if (conf.int) {
    
    bootresult = boot(n.boot, data, ratio.name, status.name, time0.name, ratio)
    se = bootresult$se
    low = c(1, exp(-exp(log(-log(surv[-1]))-1.96*se[-1]/surv[-1]/log(surv[-1]))))
    upp = c(1, exp(-exp(log(-log(surv[-1]))+1.96*se[-1]/surv[-1]/log(surv[-1]))))
    low.med = bootresult$low.med
    upp.med= bootresult$upp.med
    
    res = c(res, list(se=se, low=low, upp=upp, low.med=low.med, upp.med=upp.med))
    
  }
  
  # If points is not null
  if( !is.null(points) ) {
    
    index = apply(outer(ratio, points, "<="), 2, sum)
    surv.points = surv[index]
    ind = which(points > max(data[[ratio.name]]))
    
    if( length(ind) > 0 ) {
      
      surv.points[ind] = NA
      print(paste0("The GMI survival probability is not estimable at points greater than ", 
                   max(data[[ratio.name]]), 
                   ", the largest GMI value in the data"))
      
    }
    
    res = c(res, list(points=points, surv.points=surv.points))
    
    if(conf.int) {
      
      se.points = se[index]
      low.points = low[index]
      upp.points = upp[index]
      
      if (length(ind)>0){
        se.points[ind] = NA; low.points[ind] = NA; upp.points[ind] = NA
      }
      
      res = c(res, list( se.points=se.points, low.points=low.points, upp.points=upp.points))
      
    }
    
  }
  
  return(res)
  
}
