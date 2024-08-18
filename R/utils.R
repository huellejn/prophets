#' Check if all required columns are necessary
#' @param data a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )
#' @param required_columns vector of expected column names
#' @noRd
check_columns <- function(data, required_columns) {
  present_columns <- names(data)
  missing_columns <- required_columns[! required_columns %in% present_columns]
  if(length(missing_columns) > 0) stop( paste("The following required columns are missing from the data:", missing_columns))
}


#' Check if data has the appropiate format
#' @param 
#' @noRd
check_data <- function(value, type_expected, range) {
  
  # expected types
  expected_type <- list(
    alpha = "numeric",
    delta = "numeric",
    prob = "numeric",
    plot = "logical",
    min_pfs2 = "numeric",
    null_HR = "numeric",
    alt_HR = "numeric",
    rho = "numeric",
    model = "vector",
    selected_PFS = "vector",
    verbose = "logical",
    ges = "numeric",
    sample_size = "integer",
    power = "numeric",
    lost = "numeric",
    pfsratio = "numeric",
    p0 = "numeric",
    p1 = "numeric",
    beta = "numeric",
    k = "numeric",
    tidy = "logical"
  )
  
  # expected ranges
  expected_type <- list(
    alpha = c(min = 0, max = 1),
    delta = c(min = 0, max = Inf),
    prob = c(min = 0, max = 1),
    min_pfs2 = c(min = 0, max = 1),
    null_HR = c(min = 0, max = 1),
    alt_HR = c(min = 0, max = 1),
    rho = c(min = 0, max = 1),
    model = "vector",
    selected_PFS = "vector",
    ges = c(min = 0, max = 1),
    sample_size = c(min = 0, max = Inf),
    power = c(min = 0, max = 1),
    lost = c(min = 0, max = 1),
    pfsratio = c(min = 0, max = 1),
    p0 = c(min = 0, max = 1),
    p1 = c(min = 0, max = 1),
    beta = c(min = 0, max = 1),
    k = c(min = 0, max = 1)
  )
  
}

#' Melt data
#' @param data Data frame with columns PFS1, PFS2 and status.
#' @noRd
melt_data <- function(data) {
  
  data <- data %>%
    dplyr::mutate(
      status1 = 1,
      idx = 1:nrow(.)
    ) %>%
    dplyr::select(idx, PFS1, status1, PFS2, status2 = status) %>%
    tidyr::pivot_longer(cols = c("PFS1", "PFS2", "status1", "status2"), 
                 names_to = c(".value", "ct.line"), 
                 names_pattern = "^([A-Za-z]+)(\\d+)"
    ) 
  
  return(data)
  
}

#' Silverman: Utils function for the kernelKM_PFSr function
#' @param x a matrix of PFS1 values
#' @noRd
silverman<-function(x)
{
  abs(0.5*exp(-abs(x)/sqrt(2))*sin(abs(x)/sqrt(2)+pi/4))
}

#' calc.W: Utils function for the kernelKM_PFSr function
#' @param T0used a matrix of PFS1 values
#' @noRd
calc.W <- function(T0used){
  n = length(T0used)
  T0matrix = matrix(rep(T0used, n), ncol=n)
  sigmaterm = sd(T0used)
  nterm = n^(-2/5)
  conterm = 1
  a_n = conterm * nterm * sigmaterm
  W = silverman((t(T0matrix) - T0matrix)/a_n)
  return(W)
}

#' estimator: Utils function for the kernelKM_PFSr function
#' @param data a data frame with columns ratio, status, and PFS1 (time0)
#' @param ratio.name name of the column that contains the PFSr
#' @param status.name name of the column that contains the status
#' @param time0.name name of the column that contains PFS1
#' @param points number of points
#' @noRd
estimator <- function(data, ratio.name, status.name, time0.name, points = NULL){
  R = data[[ratio.name]]
  Delta = data[[status.name]]
  T0 = data[[time0.name]]
  n = length(T0)
  T0used = log(T0) 
  W = calc.W(T0used)
  tt = sort(unique(R[Delta==1]))
  D = length(tt)
  Rmatrix = matrix(rep(R, D), ncol=D)
  ttmatrix = matrix(rep(tt, n), nrow=n, byrow = T)
  Rgeqtt = (Rmatrix >= ttmatrix)
  ReqttDelta = Delta * (Rmatrix == ttmatrix)
  tmpu = W %*% ReqttDelta
  tmpd = W %*% Rgeqtt
  tmp = tmpu / tmpd
  LP = 1 - tmp
  CKM = t(apply(LP, 1, cumprod))
  surv = colMeans(CKM)
  
  ## add time 0 in case index=0 
  tt = c(0, tt)
  surv = c(1, surv)
  
  if(!is.null(points)){
    index =apply(outer(tt, points, "<="), 2, sum)
    return(surv[index]) ## this output is needed in bootstrap 
  } else{
    return(list(ratio=tt, surv=surv)) ## ratio: the vector of 0 and unique ratio points 
  }
}

#' getbootestimate: Utils function for the kernelKM_PFSr function
#' @param boot.index 
#' @param data a data frame with columns ratio, status, and PFS1 (time0)
#' @param ratio.name name of the column that contains the PFSr
#' @param status.name name of the column that contains the status
#' @param time0.name name of the column that contains PFS1
#' @param points number of points
#' @noRd
getbootestimate <-function(boot.index, data, ratio.name, status.name, time0.name, points){
  return(estimator(data=data[boot.index,], 
                   ratio.name=ratio.name, status.name=status.name,
                   time0.name=time0.name, points=points))
}


#' med_func: Utils function for the kernelKM_PFSr function
#' @param surv a vector of column means 
#' @param ratio vector of PFSr
#' @noRd
med_func = function(surv, ratio){
  if (min(surv)>0.5){
    medsurv=NA
  } else {
    medsurv = ratio[sum(surv>0.5)+1]
  }
}

#' boot: Utils function for the kernelKM_PFSr function
#' @param n.boot 
#' @param data a data frame with columns ratio, status, and PFS1 (time0)
#' @param ratio.name name of the column that contains the PFSr
#' @param status.name name of the column that contains the status
#' @param time0.name name of the column that contains PFS1
#' @param points number of points
#' @noRd
boot = function(n.boot, data, ratio.name, status.name, time0.name, points){
  n = nrow(data)
  bootsample = matrix(sample(1:n,n*n.boot,replace=T), nrow=n)
  est.boot = t(apply(bootsample, 2, getbootestimate, data=data, 
                     ratio.name=ratio.name, status.name=status.name,
                     time0.name=time0.name, points=points))
  se = apply(est.boot,2,sd)
  med.boot = apply(est.boot, 1, med_func, ratio=points)
  
  if(mean(as.numeric(is.na(med.boot)))>0.1){
    print("The confidence interval for median is set to NA because the median is NA for more than 10% bootstrap samples")
    low.med = NA
    upp.med = NA
  }
  else{
    low.med = quantile(med.boot, probs=0.025, na.rm=T) 
    upp.med = quantile(med.boot, probs=0.975, na.rm=T)
  }
  return(list(se=se, low.med=low.med, upp.med=upp.med))
}
