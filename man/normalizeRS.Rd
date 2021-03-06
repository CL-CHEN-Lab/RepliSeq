% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalizeRS.R
\name{normalizeRS}
\alias{normalizeRS}
\title{Normalize Repli-seq assay}
\usage{
normalizeRS(rs_assay, norm_ratios, round_digits = 3)
}
\arguments{
\item{rs_assay}{a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx}

\item{norm_ratios}{an array with the ratios in the order of rs_assay fractions}

\item{round_digits}{amount of digits after comma to conserve in the output when rounding (default = 3)}
}
\value{
a Repli-seq assay after normalizing (dividing fractions by ratios) as a data.frame and formatted as chr,start,stop,S1,...,Sx
}
\description{
Calculate Repli-seq assay count matrices after normalizing (dividing counts by ratios)
}
