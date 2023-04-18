#' Study power calculation
#'
#' @param sample_size a numeric value. The sample size of the study.
#' @param null_HR a numeric value. H~0~ (null hypothesis, unsatisfying PFSr).
#' @param alt_HR a numeric value. H~1~ (satisfying PFSr).
#' @param rho a numeric value. the correlation between PFS1 and PFS2.
#' @param alpha a numeric value. Probability (range 0.01-1.00)
#' @param k a numeric value. Common shape parameter k (k </=/> 1 for decreasing, constant, and increasing hazard functions).
#' @param model a factor with levels c("GBVE", "Weibull"). Indicates the model used for the calculation.
#' @param verbose logical value. If TRUE, results output is in text format.
#'
#' @return the calculated required sample size for the defined parameters
#' @export
#'
#' @examples
#' PFSr_power_calculate(sample_size = 65, null_HR = 1, alt_HR = 1.3, rho = 0.1, model = "GBVE")
PFSr_power_calculate <- function(sample_size = NULL, 
                                 null_HR = 1,
                                 alt_HR = NULL, 
                                 rho = NULL, 
                                 alpha = 0.05,
                                 k = 1, 
                                 model = c("GBVE", "Weibull"),
                                 verbose = FALSE) {
  
  x <- c("GBVE", "Weibull")
  if( ! model  %in% x)
    stop("Error: specified model is not GBVE or Weibull")
  
  if(alpha < 0.01 | alpha > 1.00)
    stop("Error: alpha is out of range")
  if(sample_size < 1)
    stop("Error: sample size needs to be > 0!")
  
  if(rho < -1 | rho > 1)
    stop("Error: rho is out of range")
  
  if(alt_HR <= null_HR)
    stop("Error: alt_HR must be higher than null_HR")
  
  v.fnct <- function(v) {(2*gamma(v+1)^2/gamma(2*v+1)-1) - rho}
  v <- uniroot(v.fnct, c(0,1))$root
  HR <- alt_HR/null_HR
  
  if(model=="GBVE") {
    
    ges <- (1+HR^(-1/v))^-1
    
  } else {
    
    ges <- 1/(1+HR^-k)
    
  }
  
  delta <- 4*sample_size*(ges-0.5)^2
  power <- 1-pchisq(qchisq(alpha, 1,0, lower.tail = FALSE), 1, delta)
  
  if(verbose == FALSE) {
    res <- data.frame(
      model = model,
      alpha = alpha,
      sample_size = sample_size,
      alt_HR = alt_HR,
      null_HR = null_HR,
      rho = rho,
      power = power
    )
  } else {
    cat(paste("For a clinical study adopting the PFSratio as primary endpoint, with:",
              paste("- a ", model,"-based model,", sep = ""),
              paste("- alpha = ", alpha,",", sep = ""),
              paste("- a known or expected sample size of ", sample_size, " patients,", sep = ""),
              paste("- assuming an alternative PFSratio of ", alt_HR, ",", " and null PFSratio of ", null_HR, ",", sep = ""),
              paste("- and a ", ifelse(rho<0.3, "weak ", ifelse(rho<0.7, "moderate ", "strong ")), "(", rho, ") ", "correlation between PFS1 and PFS2,", sep = ""),
              paste("the study power for the PFSratio-based analysis is of: ", round(power*100,0), "%", sep = ""),
              sep="\n")
    )
    
  }
}