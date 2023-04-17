#' Study sample size calculation according to a pre-specified proportion of patients with positive PFSr
#'
#' @param p0 a numeric value. Proportion of patients with a PFSratio > pfsratio according to the null hypothesis
#' @param p1 a numeric value. True proportion of patients with a PFSratio > pfsratio
#' @param alpha a numeric value. Probability (range 0.01-1.00)
#' @param beta a numeric value.
#' @param pfsratio a numeric value. Desired PFSr that the study should achieve.
#' @param lost a numeric value. Expected proportion of non-informative pairs.
#' @param verbose a logical value. If TRUE, results output is in text format.
#'
#' @return
#' @export
#' @importFrom gsDesign nBinomial1Sample
#'
#' @examples
PFSr_samplesize_proportion <- function( p0, 
                                        p1, 
                                        alpha, 
                                        beta, 
                                        pfsratio, 
                                        lost = 0.1,
                                        verbose = FALSE) {
  
  np <- gsDesign::nBinomial1Sample(p0, p1, alpha, beta, n = 1:1000, outtype = 1, conservative = FALSE)
  r <- qbinom(p = 1 - alpha, size = np, prob = p0) + 1
  N = np/(1-lost)
  
  if(verbose == FALSE) {
    
    res <- data.frame(
      alpha = alpha,
      beta = beta,
      power = 1-beta,
      p0 = p0,
      p1 = p1,
      PFSr = pfsratio,
      lost = lost,
      r = r,
      n = np
    )
    
    return(res)
    
  } else {
    cat(paste(paste("For a clinical trial adopting the PFSratio threshold of ", pfsratio, " as primary endpoint, with:", sep = ""),
              paste("- alpha = ", alpha,",", sep = ""),
              paste("- power = ", round((1-beta)*100,0), "%,", sep = ""),
              paste("- assuming a true proportion of patients with a PFSratio >", pfsratio,  " as equal to ", round(p1*100,1),"%", sep = ""),
              paste("- a null hypothesis that the PFSratio is > ", pfsratio, " in ", round(p0*100,1),"% of patients", sep = ""),
              paste("- with an expected proportion of non-informative pairs of ", round(lost*100,0), "%,", sep = ""),
              paste("the required sample size is of: ", round(np,0), " patients, with a critical value for H0 rejection: r < ",r, sep = ""),
              sep="\n")
    )
    
  }
} 