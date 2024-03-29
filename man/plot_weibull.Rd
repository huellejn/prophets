% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_weibull.R
\name{plot_weibull}
\alias{plot_weibull}
\title{Plot Weibull distribution for PFS values}
\usage{
plot_weibull(data)
}
\arguments{
\item{data}{a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric), 'status' (logical or integer with values c(0,1) )}
}
\value{
A plot of the Weibull distribution for PFS values
}
\description{
Plot Weibull distribution for PFS values
}
\examples{
data(input)
plot_weibull(data = input)
}
