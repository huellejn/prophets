#' Study sample size calculation
#'
#' @param alpha a numeric value. between 0.01-1.00. Probability
#' @param power a numeric value. statistical power
#' @param rho a numeric value. the correlation between PFS1 and PFS2.
#' @param alt_HR a numeric value. H~1~ (satisfying PFSr).
#' @param null_HR a numeric value. H~0~ (null hypothesis, unsatisfying PFSr).
#' @param k a numeric value. Common shape parameter k (k </=/> 1 for decreasing, constant, and increasing hazard functions).
#' @param lost a numeric value. expected proportion of non-informative pairs
#' @param model a factor with levels c("GBVE", "Weibull"). Indicates the model used for the calculation.
#' @param verbose if TRUE, results output is in text format
#'
#' @return the calculated required sample size for the defined parameters
#' @export
#'
#' @examples
PFSr_samplesize <- function(alpha = 0.05, 
                            power = 0.8,
                            rho = 0.5,
                            alt_HR = 1.33,
                            null_HR = 1,
                            k = 1,
                            lost = 0.1,
                            model = c("GBVE", "Weibull"),
                            verbose = FALSE) {
  
  model_available <- c("GBVE", "Weibull")
  if( ! model  %in% model_available)
    stop("Error: specified model is not GBVE or Weibull")
  
  if(alpha < 0.01 | alpha > 1.00)
    stop("Error: alpha is out of range.")
  if(power < 0.1 | power > 1.00)
    stop("Error: power is out of range")
  
  if(rho < -1 | rho > 1)
    stop("Error: rho is out of range")
  
  if(alt_HR <= null_HR)
    stop("Error: alt_HR must be higher than null_HR")
  
  delta.fnct <- function(delta, alph = alpha, pow = power) { 
    1-pchisq(qchisq(alph, 1,0, lower.tail = FALSE), 1, delta) - pow } 
  v.fnct <- function(v, r = rho) {(2*gamma(v+1)^2/gamma(2*v+1)-1) - r}
  
  # Calculations 
  HR <- alt_HR/null_HR
  
  # Non-central parameter
  delta <- uniroot(delta.fnct, c(0, 15))$root 
  
  # Dependence parameter v (v=1 corresponds to independence)
  v <- uniroot(v.fnct, c(0,1))$root
  
  if(model == "GBVE") {
    
    # GBVE model
    ## Dependence parameter v (v=1 corresponds to independence)
    v <- uniroot(v.fnct, c(0,1))$root
    ## Calculate the generalized treatment effect size (p_ges): prob(GMI)>1
    p_ges <- (1+HR^(-1/v))^-1 
    
  } else {
    
    # Weibull model
    ## Calculate the generalized treatment effect size (p_ges): prob(GMI)>1
    p_ges <- 1/(1+HR^-k)
    
  }
  
  # Calculate the number of paired events (np)
  np <- delta/4/(p_ges-0.5)^2 
  
  N = np/(1-lost)
  
  if(verbose == FALSE) {
    res <- data.frame(
      model = model,
      alpha = alpha,
      power = power,
      PFSratio_alternative = alt_HR,
      PFSratio_null = null_HR,
      rho = rho,
      correlation = ifelse(rho < 0.3, "weak", ifelse(rho < 0.7, "moderate", "strong")),
      lost = lost,
      n = N
    )
    return(res)
      
  } else {
    cat(paste("For a clinical trial adopting the PFSratio as primary endpoint, with:",
              paste("- a ", model,"-based model,", sep = ""),
              paste("- alpha = ", alpha,",", sep = ""),
              paste("- power = ", round(power*100,0), "%,", sep = ""),
              paste("- assuming an alternative PFSratio of ", alt_HR, ",", " and null PFSratio of ", null_HR, ",", sep = ""),
              paste("- and a ", ifelse(rho<0.3, "weak ", ifelse(rho<0.7, "moderate ", "strong ")), "(", rho, ") ", "correlation between PFS1 and PFS2,", sep = ""),
              paste("- with an expected proportion of non-informative pairs of ", round(lost*100,0), "%,", sep = ""),
              paste("the required sample size is of: ", round(N,0), " patients", sep = ""),
              sep="\n")
    )
    
  }
}