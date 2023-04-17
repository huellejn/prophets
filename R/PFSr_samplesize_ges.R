#' Study sample size calculation using the generalized effect size method
#'
#' @param ges a numeric value. Generalized effect size.
#' @param alpha  numeric value. Probability (range 0.01-1.00)
#' @param power numeric value. Statistical power
#' @param lost numeric value. Expected proportion of non-informative pairs.
#' @param verbose logical value. If TRUE, results output is in text format.
#'
#' @return
#' @export
#'
#' @examples
ges_PFSr_samplesize <- function( ges = 0.3, 
                                 alpha = 0.05, 
                                 power = 0.80, 
                                 lost = 0.1, 
                                 verbose = FALSE) {
  
  delta.fnct <- function(delta) { 1-pchisq(qchisq(alpha, 1,0, lower.tail = FALSE), 1, delta) - power }
  delta <- uniroot(delta.fnct, c(0, 15))$root
  np <- delta/4/(ges-0.5)^2
  N = np/(1-lost)
  
  if(verbose == FALSE) {
    res <- data.frame(
      alpha = alpha,
      power = power,
      ges = ges,
      lost = lost,
      n = N
    )
    return(res)
  } else {
    
    cat(paste("For a clinical trial adopting the PFSratio as primary endpoint, with:",
              paste("- alpha = ", alpha,",", sep = ""),
              paste("- power = ", round(power*100,0), "%,", sep = ""),
              paste("- assuming a probability of PFSratio being greater than ∂ of ", ges, ",", sep = ""),
              paste("- with an expected proportion of non-informative pairs of ", round(lost*100,0), "%,", sep = ""),
              paste("the required sample size is of: ", round(N,0), " patients", sep = ""),
              sep="\n")
    )
    
  }
  
}