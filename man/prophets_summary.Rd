% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prophets_summary.R
\name{prophets_summary}
\alias{prophets_summary}
\title{Title}
\usage{
prophets_summary(data, delta = 1)
}
\arguments{
\item{data}{a dataframe with columns 'PFS1' (numeric), 'PF2' (numeric)}

\item{delta}{a numeric value. The value 'delta' is the desired difference between PFS2 and PFS1 that is seen as a success.}
}
\value{
a table with PFSr analysis results using four different methods
}
\description{
Title
}
\examples{
data(input)
input$ratio <- input$PFS2/input$PFS1
prophets_summary(data = input, delta = 1.3)
}
