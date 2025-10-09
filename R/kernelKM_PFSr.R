#' Kernel conditional Kaplan Meier method
#'
#' @param data a data frame with values for PFS1/PFS2 ratio, status (1/0) and PFS1 (time0)
#' @param ratio.name name of the column that contains the PFS1/PFS2 ratio ('numeric')
#' @param status.name name of the column that contains the status ('boolean')
#' @param time0.name name of the column that contains the PFS1 ('numeric')
#' @param delta a numeric value indicating the desired difference between PFS2 and PFS1 that is seen as a success.
#' @param conf.int a boolean indicating if confidende intervals should be calculated
#' @param n.boot an integer value indicating the number of permutations to obtain confidence intervals for the survival estimates
#'
#' @return a Kaplan-Meier based statistics according to the kernel conditional KM method
#' @export
#'
#' @examples
#' # Example with default values:
#' kernelKM_PFSr(data, "ratio", "status", "PFS1")
#' # Example with modified values
#' kernelKM_PFSr(data, "ratio", "status", "PFS1", delta = 2, conf.int = TRUE)
kernelKM_PFSr <- function(
  data,
  ratio.name,
  status.name,
  time0.name,
  delta = NULL,
  conf.int = FALSE,
  n.boot = 2000
) {
  # Get base estimates
  result <- estimator(data, ratio.name, status.name, time0.name)
  ratio <- result$ratio
  surv <- result$surv

  # Calculate median survival
  medsurv <- ifelse(min(surv) > 0.5, NA, ratio[sum(surv > 0.5) + 1])

  # Helper function to get survival at delta points
  get_surv_at_delta <- function(surv_vec, ratio, delta, data, ratio.name) {
    index <- apply(outer(ratio, delta, "<="), 2, sum)
    surv_points <- surv_vec[index]

    # Handle values beyond data range
    ind <- which(delta > max(data[[ratio.name]]))
    if (length(ind) > 0) {
      surv_points[ind] <- NA
      message(
        "PFSratio survival probability is not estimable at delta > ",
        max(data[[ratio.name]]),
        " (max PFSratio in data)"
      )
    }

    surv_points
  }

  # Without confidence intervals
  if (!conf.int) {
    if (!is.null(delta)) {
      surv_points <- get_surv_at_delta(surv, ratio, delta, data, ratio.name)
      return(list(
        ratio = ratio,
        surv = surv,
        `median PFSratio` = medsurv,
        delta = delta,
        `PFSr_estimator` = surv_points
      ))
    }
    return(list(ratio = ratio, surv = surv, `median PFSratio` = medsurv))
  }

  # With confidence intervals
  bootresult <- boot(n.boot, data, ratio.name, status.name, time0.name, ratio)
  se <- bootresult$se

  # Calculate confidence bounds
  log_transform <- log(-log(surv[-1]))
  se_ratio <- se[-1] / (surv[-1] * log(surv[-1]))

  low <- c(1, exp(-exp(log_transform - 1.96 * se_ratio)))
  upp <- c(1, exp(-exp(log_transform + 1.96 * se_ratio)))

  if (!is.null(delta)) {
    # Get point estimates at delta
    surv_points <- get_surv_at_delta(surv, ratio, delta, data, ratio.name)
    index <- apply(outer(ratio, delta, "<="), 2, sum)

    out <- tibble(
      delta = delta,
      estimate = surv_points,
      conf.low = low[index],
      conf.high = upp[index]
    )

    survtab <- tibble(
      `median PFSratio` = medsurv,
      conf.low = bootresult$low.med,
      conf.high = bootresult$upp.med
    )

    return(list(`median PFSratio` = survtab, `PFSr_estimator` = out))
  }

  list(
    ratio = ratio,
    surv = surv,
    se = se,
    low = low,
    upp = upp,
    `median PFSratio` = medsurv
  )
}
