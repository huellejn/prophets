% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PFSr_power_calculate.R
\name{PFSr_power_calculate}
\alias{PFSr_power_calculate}
\title{Study power calculation}
\usage{
PFSr_power_calculate(
  sample_size = NULL,
  null_HR = 1,
  alt_HR = NULL,
  rho = NULL,
  alpha = 0.05,
  k = 1,
  model = c("GBVE", "Weibull"),
  verbose = FALSE
)
}
\arguments{
\item{sample_size}{a numeric value. The sample size of the study.}

\item{null_HR}{a numeric value. H~0~ (null hypothesis, unsatisfying PFSr).}

\item{alt_HR}{a numeric value. H~1~ (satisfying PFSr).}

\item{rho}{a numeric value. the correlation between PFS1 and PFS2.}

\item{alpha}{a numeric value. Probability (range 0.01-1.00)}

\item{k}{a numeric value. Common shape parameter k (k </=/> 1 for decreasing, constant, and increasing hazard functions).}

\item{model}{a factor with levels c("GBVE", "Weibull"). Indicates the model used for the calculation.}

\item{verbose}{logical value. If TRUE, results output is in text format.}
}
\value{
the calculated required sample size for the defined parameters
}
\description{
Study power calculation
}
\examples{
PFSr_power_calculate(sample_size = 65, null_HR = 1, alt_HR = 1.3, rho = 0.1, model = "GBVE")
}
