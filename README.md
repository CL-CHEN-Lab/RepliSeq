# RepliSeq <a><img src='https://github.com/SamiLhll/RepliSeq/blob/7f8770efaa51c2e0ac824576196ae4fd5f58b837/inst/img/Repliseq_logo.png' align="right" height="230" /></a>

This is an R package which aims at processing Repli-seq data. It takes raw counts (Bedgraph file format) as input and makes it then quick and easy to further analyze the DNA replication timing with a set of functions to manipulate and vizualize the data.   
RepliSeq functions include **loading** multi-fraction Repli-seq assay data as count matrices (from 2 to N fractions depending on the experimental design) but also **rescaling** profiles to any resolution and calculting the **Replication timing** as the S50 (moment in the S-phase where a loci reaches 50% of its total replication on a scale from 0, early, to 1, late)

<!-- badges: start -->
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->


-----------------------------------------------------------------------  

## Installation:

You can install this package by entering the following within R:

```{r}
# get the development version from GitHub using devtools :
# install.packages("devtools")
devtools::install_github("SamiLhll/RepliSeq",build_vignettes = TRUE)
# building the vignette makes the installation a bit longer but its mandatory so ou can access it by doing :   
vignette("RepliSeq")

```

## Requirements:

The function **writeBigWig()** requires UCSC's **wigToBigWig** application to be installed on the computer.   
It can be found at [encodeproject](https://www.encodeproject.org/software/wigtobigwig/) 


## Usage : 

Check out the vignette for extended documentation.

#### Loading repliseq data :

The function *readRS(path_data,fractions)* reads Repli-seq assays from multiple files (one file per fraction) and returns a dataframe.   
It requires bedgraph inputs [(see bedgraph specifications)](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) with a one line header but no other comments such as: 

track 	type=bedGraph 	name=NT_chr22-s1	description=50kb   
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


#### Compute the replication timing (S50) :

The function *calculateS50(rs_assay)* returns a dataframe composed of genomic coordinates associated with replication timing as an S50 [(Chen et al. (2010))](https://doi.org/10.1101/gr.098947.109) value comprised within 0 (early replicating) and 1 (late replicating).

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

# Compare the total replication among Repli-Seq assays :

As introduced in [Brison,.O, El-Hilali,S. et al. (2019)](https://doi.org/10.1038/s41467-019-13674-5), Repli-Seq assays could be compared to quantitatively assess which parts of DNA were the most affected by Aphidicolin. The function *calculateURI()* calculates this Under Replication Index (URI) from two Repli-Seq assays loaded with *readRS()*.

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

  
## Getting help :
  
Need help, Identified a bug, or want to see other features implemented ?   
Feel free to open an issue here or send an email to the authors :   
  
[Sami EL HILALI](mailto:elhilali.sami@gmail.com) and [Chunlong CHEN](mailto:chunlong.chen@curie.fr)


## References: 

Brison O., El-Hilali S., Azar, D., Koundrioukoff1 S., Schmidt M., Naehse-Kumpf V., Jaszczyszyn Y., Lachages A.M., Dutrillaux B., Thermes C., Debatisse M. and Chen C.L. (2019) [Transcription-Mediated Organization of the Replication Initiation Program Across Large Genes Sets Up Common Fragile Sites Genome-Wide.](https://doi.org/10.1038/s41467-019-13674-5) ***Nat. Commun.*** 10, 5693

Chen C.L., Rappailles A., Duquenne L., Huvet M., Guilbaud G., Farinelli L, Audit B, d'Aubenton-Carafa Y., Arneodo A., Hyrien O. and Thermes C. (2010) [Impact of replication timing on non-CpG and CpG substitution rates in mammalian genomes](https://genome.cshlp.org/content/20/4/447.long). ***Genome. Res.*** 20, 447-457. 



  
