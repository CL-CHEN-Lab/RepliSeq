# RepliSeq

Analysis of Repli-Seq data to study DNA replication timing program in R.

## Description :

An R package that features a set of functions to conduct Repli-seq data analysis in R.  
We propose with this package to analyze Repli-seq data from within data.frames which lets you easily 
complete your analysis with dplyr, calculate intersections with tidygenomics and plot your results with ggplot vizualizations.   
RepliSeq functions includes **loading** multi-fractions Repli-seq assays data as count matrices 
(from 1 fraction for controls to x defined by you hardware capabilities) ; **normalization** ; **rescaling** ; 
calculting metrics such as **S50** (replication timing) and **URI** (Under replication index got 
from two repliseq assays comparison).


## Installation :

You can install this package by entering the following from within R :

```{r}

devtools::install_github("CL-CHEN-Lab/RepliSeq")

```
