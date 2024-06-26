% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_functions.R
\name{get_isodists}
\alias{get_isodists}
\title{Get isotopic distributions}
\usage{
get_isodists(seqs, ndeam, nhyds, long_format = F)
}
\arguments{
\item{nhyds}{Integer, vector}
}
\value{
Array of intensity distributions, without the masses
}
\description{
Wrapper for bacollite ms_iso.
Calculates theoretical isotopic deistributions for a list of
sequences
}
