% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screening-tests.R
\name{screening-tests}
\alias{screening-tests}
\alias{screening_globalTest}
\alias{screening_hotellingT}
\alias{screening_invnormT}
\alias{screening_nettleton}
\title{Global test procedures}
\usage{
screening_hotellingT(set, perm, d, l)

screening_globalTest(set, perm, d, l)

screening_invnormT(set, perm, d, l)

screening_nettleton(set, perm, d, l)
}
\arguments{
\item{set}{identifiers of probes included in the gene set}

\item{perm}{number of permutations to be applied in resampling tests}

\item{d}{numeric matrix containing expression values}

\item{l}{vector or factor containing treatment/phenotype labels}
}
\value{
p-value of the test applied to the data
}
\description{
Various global test procedures that can be applied in the screening stage of twoStageGSA
}
\author{
Florian Klinglmueller
}

