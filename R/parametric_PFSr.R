#' Analysis of progression free survival ratio (PFSr) using the parametric method
#'
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) ) and 'ratio' (numeric, ratio PFS2/PFS1)
#' @param delta a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.
#'
#' @return Parametric-based statistics for PFSr analysis
#' @export
#' @importFrom survival survreg Surv
#' @importFrom msm deltamethod
#' @importFrom tibble tibble
#'
#' @examples
#' data(input)
#' input$ratio <- input$PFS2/input$PFS1
#' parametric_PFSr(data = input, delta = 1.3)
parametric_PFSr <- function(data, 
                            delta = 1) {
  
  # Check if required columns are present
  check_columns(data = data, required_columns = c("status", "ratio"))
  
  obj_survreg = survival::survreg(survival::Surv(ratio, status) ~ 1, data = data, dist = "loglogistic")
  
  shape = exp(-obj_survreg$coef)
  scale = 1/(obj_survreg$scale)
  estmean = c(obj_survreg$coef, obj_survreg$scale)
  estim = 1/(1+(shape*delta)^(scale))
  form <- sprintf("~1/(1+(exp(-x1 * %f ))^(1/x2))", delta)
  sd = msm::deltamethod(as.formula(form), estmean, obj_survreg$var)
  mean = (pi/(shape*scale*sin(pi/scale)))
  
  res = tibble::tibble(method = "parametric",
               delta = delta,
               estimate = estim,
               conf.low = estim + qnorm(.025) * sd,
               conf.high = estim + qnorm(.975) * sd
  )
  
  return(res)
  
  }