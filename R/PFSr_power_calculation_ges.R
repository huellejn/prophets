#' Study power calculation using the GES method
#'
#' @param sample_size a numeric value. The sample size of the study.
#' @param ges a numeric value. Generalized effect size.
#' @param alpha a numeric value. Probability (range 0.01-1.00)
#' @param verbose logical value. If TRUE, results output is in text format.
#'
#' @return the calculated required sample size for the defined parameters
#' @export
#'
#' @examples
PFSr_power_calculation_ges <- function( sample_size = NULL, 
                                        ges = 0.3, 
                                        alpha = 0.05,
                                        verbose = FALSE) {
  
  delta <- 4*sample_size*(ges-0.5)^2
  power <- 1-pchisq(qchisq(alpha, 1,0, lower.tail = FALSE), 1, delta)
  
  if(verbose == FALSE) {
    res <- data.frame(
      alpha = alpha,
      sample_size = sample_size,
      ges = ges,
      power = power
    )
    return(res)
    
  } else {
    
    cat(paste("For a clinical study adopting the PFSratio as primary endpoint, with:",
              paste("- alpha = ", alpha,",", sep = ""),
              paste("- a known or expected sample size of ", sample_size, " patients,", sep = ""),
              paste("- assuming a probability of PFSratio being greater than ∂ of ", ges, ",", sep = ""),
              paste("the study power for the PFSratio-based analysis is of: ", round(power*100,0), "%", sep = ""),
              sep="\n")
    )
    
  }
}