#' Read Repli-seq assay data (bedgraph format).
#' Bedgraphs with a one line header and with similar organization (same length and sorted)
#' @param paths_data list of paths to repliseq fractions files
#' @param fractions list of fractions names
#' @return a dataframe composed of genomic coordinates plus all the specified fractions as for example
#'     chr,start,stop,S1,S2,S3,S4,S5,S6

readRS <- function(paths_data,fractions) {
  #### Control that path data and fractions are the same length
  if (length(paths_data) != length(fractions)) {
    stop("Error : paths_data and fractions argument don't have the same sizes")
  }
  ### read genomic coordinates from first fraction :
  rs_temp <- read.table(paths_data[1],header=FALSE,skip=1,sep="") %>%
    select(V1,V2,V3) %>%
    rename(chr = V1, start = V2, stop = V3)
  ### read other fractions and create the corresponding columns in rs_temp
  for (i in c(1:length(path_data))) {
    temp_fraction <- fractions[i]
    temp_path <- paths_data[i]
    temp_fraction_df <- read.table(i,header=FALSE,skip=1,sep="")
    rs_temp[temp_fraction] <- as.list(temp_fraction_df$V4)
  }
  return(rs_temp)
}
