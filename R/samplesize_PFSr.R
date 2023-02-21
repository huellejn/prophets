#' Calculation of sample size for the design of a clinical trial with PFSr as an initial endpoint
#'
#' @param alpha a numeric value between 0.01 and 1. Desired level of significance.
#' @param power  a numeric value between 0.1 and 1. Desired power (1-beta) value of the test.
#' @param rho a numeric value between -1 and 1. Pearson correlation coefficient
#' @param alt_HR  a numeric value that is <= null_HR. Desired minimal hazard ratio of the alternative hypothesis (PFS2 is larger than PFS1).
#' @param null_HR  a numeric value that is > alt_HR. Hazard ratio of the null hypothesis (PFS2 is equal to PFS1).
#' @param k  a numeric value. Common shape parameter for decreasing (<1), constant (=1), and increasing (>1) hazard functions.
#' @param lost a numeric value
#' @param model a character value that is either 'GBVE' or 'Weibull'. Sets the model for the sample size calculation, either Gumbel's type B bivariate extreme-value (GBVE) or the Weibull model.
#'
#' @return A named vector with the input parameters and the calculated required number of samples (N).
#' @export
#'
#' @importFrom stats pchisq qchisq uniroot
#'
#' @examples
#' samplesize_PFSr(model = "Weibull")
samplesize_PFSr <- function(alpha = 0.05,
                            power = 0.8,
                            rho = 0.5,
                            alt_HR = 1.33,
                            null_HR = 1,
                            k = 1,
                            lost = 0.1,
                            model = c("GBVE", "Weibull")
                            ) {

  # Check data validity
  x <- c("GBVE", "Weibull")

  if( ! model  %in% x | ! length(model)  == 1)
    stop("Error: Specified model is not 'GBVE' or 'Weibull'.")

  if(alpha < 0.01 | alpha > 1)
    stop("Error: Alpha is out of range. The valid range is 0.01 - 1.")

  if(power < 0.1 | power > 1)
    stop("Error: Power is out of range. The valid range is 0.1 - 1.")

  if(rho < -1 | rho > 1)
    stop("Error: rho is out of range. The valid range is -1 - 1.")

  if(alt_HR <= null_HR)
    stop("Error: alt_HR must be higher than null_HR.")

  delta.fnct <- function(delta, alph = alpha, pow = power) {
    1-pchisq(qchisq(alph, 1, 0, lower.tail = FALSE), 1, delta) - pow
  }

  v.fnct <- function(v, r = rho) {
    (2*gamma(v+1)^2/gamma(2*v+1)-1) - r
    }

  # Calculations
  HR <- alt_HR/null_HR

  ## Non-central parameter
  delta <- uniroot(delta.fnct, c(0, 15))$root

  ## Dependence parameter v (v=1 corresponds to independence)
  ## Generalized treatment effect size p_ges : prob(GMI)>1

  if(model == "GBVE"){

    # GBVE model
    v <- uniroot(v.fnct, c(0,1))$root
    p_ges <- (1+HR^(-1/v))^-1

  } else {

    # Weibull model
    v <- uniroot(v.fnct, c(0,1))$root
    p_ges <- 1/(1+HR^-k)

  }

  ## Number of paired events np
  np <- delta/4/(p_ges-0.5)^2
  N = np/(1-lost)

  # Export results
  res <- list(
    model = model,
    alpha = alpha,
    power = power,
    alt_HR = alt_HR,
    null_HR = null_HR,
    rho = rho,
    lost = lost,
    N = N)

  return(res)

}

