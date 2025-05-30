% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_plot.cSEMResults.R
\name{buildDotCode}
\alias{buildDotCode}
\title{Internal: Build DOT code for the SEM plot, including construct correlations.}
\usage{
buildDotCode(
  title,
  graph_attrs,
  constructs,
  r2_values,
  measurement_edge_fun,
  path_coefficients,
  path_p_values,
  correlations,
  plot_significances,
  plot_correlations,
  plot_structural_model_only,
  plot_labels,
  is_second_order = FALSE,
  model_measurement = NULL,
  model_error_cor = NULL,
  construct_correlations = NULL,
  indicator_correlations = NULL
)
}
\arguments{
\item{title}{The title of the plot.}

\item{graph_attrs}{Optional graph attributes.}

\item{constructs}{A vector of constructs.}

\item{r2_values}{Named vector of R2 values.}

\item{measurement_edge_fun}{Function to generate measurement edge code.}

\item{path_coefficients}{Matrix/data frame of path coefficients.}

\item{path_p_values}{Named vector of path p-values. Used for construct correlations too.}

\item{correlations}{List containing correlations (exogenous and indicator).}

\item{plot_significances}{Logical. Whether to display significance levels.}

\item{plot_correlations}{Option for indicator correlations ("none", "exo", or "all").}

\item{plot_structural_model_only}{Logical. Whether to display only the structural model.}

\item{is_second_order}{Logical. Whether the model is second-order.}

\item{model_measurement}{a matrix. The measurement matrix.}

\item{model_error_cor}{a matrix.}

\item{construct_correlations}{A matrix. The construct correlation matrix.}

\item{indicator_correlations}{A matrix. The indicator correlation matrix.}
}
\value{
A character string containing the complete DOT code.
}
\description{
Constructs the DOT script for the SEM path diagram, now including
correlations between constructs (not just exogenous ones). Correctly
handles drawing only one edge per correlation.
}
\keyword{internal}
