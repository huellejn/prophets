#' Analysis of PFSr using the parametric method
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
parametric_PFSr <- function(data, delta = 1, tidy=T,
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
  reg = survreg(Surv(ratio, status)~1, data = data, dist="loglogistic")
  shape = exp(-reg$coef)
  scale = 1/(reg$scale) # I modified the function to have delta here, is this correct?
  estmean = c(reg$coef,reg$scale)
  estim = 1/(1+(shape*delta)^(scale))
  form <- sprintf("~1/(1+(exp(-x1 * %f ))^(1/x2))", delta)
  sd =deltamethod(as.formula(form), estmean, reg$var)
  mean = (pi/(shape*scale*sin(pi/scale)))
  res = tibble("Method" = "Parametric",
               "estimate"=estim,
               "delta" = delta,
               # "SD"=sd,
               "conf.low" = estim + qnorm(.025) * sd,
               "conf.high" = estim + qnorm(.975) * sd
  )
  res <- res %>%
    # mutate(delta=delta) %>%
    relocate(delta, .before = estimate) %>%
    relocate(Method, .before = delta)

  if(tidy==T){
    return(res)
  }else{
    res <- res %>%
      gt(caption = "Parametric-based method PFSratio result") %>%
      fmt_number(columns = 2:5,decimals = 2)
    return(res)
  }
}
