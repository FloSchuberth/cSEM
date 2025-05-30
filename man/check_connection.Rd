% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_plot.cSEMResults.R
\name{check_connection}
\alias{check_connection}
\title{Internal: Check whether two indicators belong to the same construct.}
\usage{
check_connection(
  .indicator1,
  .indicator2,
  .model_measurement,
  .model_error_cor
)
}
\arguments{
\item{.indicator1}{Character string. The name of the indicator 1.}

\item{.indicator2}{Character string. The name of the indicator 1.}

\item{.model_measurement}{Matrix. The measurement matrix indicating the relationship
between constructs and indicators.}

\item{.model_error_cor}{Matrix. The matrix indicates the error correlation structure.}
}
\value{
TRUE if both indicators belong to the same construct, FALSE otherwise.
}
\description{
Checks whether two indicators belong to the same construct.
}
\keyword{internal}
