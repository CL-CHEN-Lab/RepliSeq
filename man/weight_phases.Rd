% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weight_phases.R
\name{weight_phases}
\alias{weight_phases}
\title{Weighting phases}
\usage{
weight_phases(
  rs_assay,
  G1_peak,
  G2_peak,
  gates,
  SimSPhase = NA,
  TotalReads = NA
)
}
\arguments{
\item{rs_assay}{a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx}

\item{G1_peak}{aproximate position of the G1 peak}

\item{G2_peak}{aproximate position of the G2 peak}

\item{gates}{a data.frame containing start and end of each gate ordered as in rs_assay}

\item{SimSPhase}{a user constumise simulation of the S-phase progression (data.frame) formatted as Percentage of replicated DNA, Number of active forks}

\item{TotalReads}{Totala ammount of reads for normalization, by defaul 15000000 per phase}
}
\value{
A  named vector containing the normalization factors for each phase to use with normalizeRS
}
\description{
Using a simulation of the S-phase (active forks in function of the percentage replicated genome), this function calculates a normalization factor for each phase based on the FACS gates and G1/G2 peaks
}
