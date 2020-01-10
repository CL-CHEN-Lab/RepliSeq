# RepliSeq   
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Analysis of Repli-Seq data to study DNA replication timing program in R.  


## Description:

An R package that features a set of functions to conduct Repli-seq data analysis.  

We propose this package to analyze Repli-seq data within data.frames, which lets you easily complete your analysis with dplyr, calculate intersections with tidygenomics and plot your results with ggplot vizualizations.

RepliSeq functions include **loading** multi-fractions (from 2 to N fractions defined by your experiment dessign and your hardware capabilities) Repli-seq assay data as count matrices; **rescaling** profiles; **smoothing** profiles; calculting metrics such as **Replication timing** (calculated as the S50, on a scale from 0, early, to 1, late, which is the moment in S phase when a sequence has been replicated in 50% of cells replication timing) and **URI** (Under replication index got from two repliseq assays comparison).


## Installation:

You can install this package by entering the following within R:

```{r}

devtools::install_github("CL-CHEN-Lab/RepliSeq")

```

## Requirements:

This package depends on:

* R version >= 3.4.4 (2018-03-15)

As mentionned in the DESCRIPTION, this packages imports: 

* dplyr (>= 0.8.3)  
* magrittr (>= 1.5)

In addition, the function **writeBigWig()** requires UCSC's **wigToBigWig** application to be installed on the computer. It can be found at [encodeproject](https://www.encodeproject.org/software/wigtobigwig/) 



## Authors:

Sami EL HILALI : elhilali.sami@gmail.com 

Chunlong CHEN : chunlong.chen@curie.fr

Don't hesitate to contact the authors or open an issue for a question or if you wish to see new features to be added to this package.



## References: 

Brison, O., El-Hilali, S., Azar, D. et al. [Transcription-Mediated Organization of the Replication Initiation Program Across Large Genes Sets Up Common Fragile Sites Genome-Wide.](https://doi.org/10.1038/s41467-019-13674-5) Nat Commun 10, 5693 (2019) doi:10.1038/s41467-019-13674-5

Chen C.L., Rappailles A., Duquenne L., Huvet M., Guilbaud G., Farinelli L, Audit B, d'Aubenton-Carafa Y., Arneodo A., Hyrien O. and Thermes C. (2010) [Impact of replication timing on non-CpG and CpG substitution rates in mammalian genomes](https://genome.cshlp.org/content/20/4/447.long). *Genome. Res.* 20, 447-457. 



## Usage examples: 

We propose an overview of some function usage. For extended documentation, please refer to the Vignette *how-to-use*.

#### readRS(path_data,fractions):

This function reads Repli-seq assays from multiple files (one file for one fraction) and outputs a dataframe from it.   
It requires bedgraph inputs [(see bedgraph spec)](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) with a one line header but **no other comments** such as:  

track 	type=bedGraph 	name=NT_chr22-s1	description=50kbprofile   
chr22	0	50000	0   
chr22	50000	100000	0   

```{r}

### args :
temp_paths <- c("../inst/extdata/NT_chr22-s1.bdg","../inst/extdata/NT_chr22-s2.bdg",
                "../inst/extdata/NT_chr22-s3.bdg","../inst/extdata/NT_chr22-s4.bdg",
                "../inst/extdata/NT_chr22-s5.bdg","../inst/extdata/NT_chr22-s6.bdg")
temp_fractions <- c("S1","S2","S3","S4","S5","S6")

### 2 fractions RepliSeq
# apply function :
RS_early <- readRS(temp_paths[1:2],temp_fractions[1:2])

### 6 fractions RepliSeq
# apply function :
RS_all <- readRS(temp_paths,temp_fractions)

### 1 fraction RepliSeq ( for S0 controls )
# apply function : 
RS_S0 <- readRS("../inst/extdata/NT_chr22-s0.bdg","S0")

### Result :

tail(RS_early)

```

| chr    | start    | stop     | S1     | S2    |
|--------|----------|----------|--------|-------|
| <fctr> | <int>    | <int>    | <dbl>  | <dbl> |
| chr22  | 51000000 | 51050000 | 12.392 | 4.929 |
| chr22  | 51050000 | 51100000 | 11.604 | 5.887 |
| chr22  | 51100000 | 51150000 | 12.568 | 7.941 |
| chr22  | 51150000 | 51200000 | 9.853  | 5.887 |
| chr22  | 51200000 | 51250000 | 2.584  | 1.711 |
| chr22  | 51250000 | 51300000 | 0.000  | 0.000 |


#### calculateS50(rs_assay):

This function returns a dataframe composed of genomic coordinates associated with replication timing as an S50 value comprised within 0 (early replicating) and 1 (late replicating).

```{r}

temp_rs <- data.frame(chr = rep("chr1",7),
                      start = seq(0,6000,1000),
                      stop = seq(1000,7000,1000),
                      S1 = c(0,0,0,1,1,1,1),
                      S2 = c(0,0,1,1,1,1,0),
                      S3 = c(0,1,1,1,1,0,0),
                      S4 = c(1,1,1,1,0,0,0))


temp_S50 <- RepliSeq::calculateS50(temp_rs)

# Result :

print(temp_S50)

```

| chr    | start | stop  | S50   |
|--------|-------|-------|-------|
| <fctr> | <dbl> | <dbl> | <dbl> |
| chr1   | 0     | 1000  | 0.875 |
| chr1   | 1000  | 2000  | 0.750 |
| chr1   | 2000  | 3000  | 0.625 |
| chr1   | 3000  | 4000  | 0.500 |
| chr1   | 4000  | 5000  | 0.375 |
| chr1   | 5000  | 6000  | 0.250 |
| chr1   | 6000  | 7000  | 0.125 |


#### calculateURI(rs_x, rs_y):

This function calculates URI between two Repli-seq assays. It returns a dataframe with the following columns:   
chr,start,stop,sum_x,sum_y,mean_xy,URI


```{r}
####### load second Repli-seq assay for comparison 
####### 6 fractions RepliSeq

# args :

aph_paths <- c("../inst/extdata/Aph_chr22-s1.bdg","../inst/extdata/Aph_chr22-s2.bdg",
               "../inst/extdata/Aph_chr22-s3.bdg","../inst/extdata/Aph_chr22-s4.bdg",
               "../inst/extdata/Aph_chr22-s5.bdg","../inst/extdata/Aph_chr22-s6.bdg")
aph_fractions <- temp_fractions

# read :

RS_aph_all <- readRS(aph_paths,aph_fractions)

# apply function :

aph_nt_uri <- calculateURI(RS_aph_all,RS_all)

# Result :

tail(aph_nt_uri)


```

| chr    | start    | stop     | sum_x  | sum_y  | mean_xy | URI         |
|--------|----------|----------|--------|--------|---------|-------------|
| <fctr> | <int>    | <int>    | <dbl>  | <dbl>  | <dbl>   | <dbl>       |
| chr22  | 51000000 | 51050000 | 27.581 | 30.869 | 29.225  | -1.37048107 |
| chr22  | 51050000 | 51100000 | 30.556 | 31.274 | 30.915  | -0.66372243 |
| chr22  | 51100000 | 51150000 | 38.770 | 36.226 | 37.498  | 0.05718338  |
| chr22  | 51150000 | 51200000 | 32.394 | 26.028 | 29.211  | 1.24529116  |
| chr22  | 51200000 | 51250000 | 10.063 | 8.533  | 9.298   | 0.82273039  |
| chr22  | 51250000 | 51300000 | 0.000  | 0.000  | 0.000   | NaN         |
  
