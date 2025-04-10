% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputefinder.R
\name{classify_missingness}
\alias{classify_missingness}
\title{Classify Missingness as MAR or MNAR Based on Intensity Distribution}
\usage{
classify_missingness(se, threshold = 0.35, min_non_na = 5, return_plot = FALSE)
}
\arguments{
\item{se}{A \code{SummarizedExperiment} object.}

\item{threshold}{A numeric value for the second derivative threshold (default = 0.35).}

\item{min_non_na}{Integer. Minimum number of non-missing values required for MAR retention (default = 5).}

\item{return_plot}{Logical. If TRUE, a ggplot object is returned alongside results.}
}
\value{
A list with:
\item{data}{Filtered \code{SummarizedExperiment} object.}
\item{MAR}{Character vector of MAR feature names.}
\item{MNAR}{Character vector of MNAR feature names.}
\item{cutoffs}{Named numeric vector of intensity cutoffs per condition.}
\item{plot}{List of ggplot objects, if \code{return_plot = TRUE}.}
}
\description{
This function identifies MNAR and MAR features using intensity distributions and missingness patterns.
}
