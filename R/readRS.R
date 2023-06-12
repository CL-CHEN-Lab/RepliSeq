#' @title Read Repli-Seq assay data (Bedgraph format).
#' @description Read and integrate all Repli-Seq fractions in a single dataframe.
#'  Fractions must be formatted as Bedgraph with a one line header and with similar resolutions.
#'
#' @param paths_data Array. List of paths to Repli-Seq fractions files
#' @param fraction_names Array. List of fractions names (default = S1, ..., Sn)
#'
#' @return a dataframe composed of genomic coordinates plus all the specified fractions as for example : chr,start,stop,S1,S2,S3,S4,S5,S6
#' @import magrittr
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom utils read.table
#' 
#' @examples 
#' # basic usage of read_RS :
#' 
#' S1_path <- system.file("extdata", "NT_chr22-s1.bdg", package="RepliSeq")
#' S2_path <- system.file("extdata", "NT_chr22-s2.bdg", package="RepliSeq")
#' 
#' my_rs_assay <- readRS(paths_data = c(S1_path , S2_path),
#'                       fraction_names = c("S1","S2"))
#' 
#' @export


readRS <- function(paths_data = NULL,
                   fraction_names = NULL) {
  
  V1 <- V2 <- V3 <- NULL
  
  ### Create default faction_names
  if (is.null(fraction_names)) {
    amount_of_fractions <- length(paths_data)
    fractions_names <- paste0("S",c(1:amount_of_fractions))
  }
  
  #### Control that path data and fractions are the same length
  if (length(paths_data) != length(fraction_names)) {
    stop("arguments paths_data and fractions must have the same sizes")
  }
  ### read genomic coordinates from first fraction :
  rs_temp <- utils::read.table(paths_data[1] , header = FALSE , skip=1 , sep="") %>%
    select(V1,V2,V3) %>%
    rename(chr = V1, start = V2, stop = V3)
  
  ### read other fractions and create the corresponding columns in rs_temp
  for (i in c(1:length(paths_data))) {
    temp_fraction <- fraction_names[i]
    temp_path <- paths_data[i]
    temp_fraction_df <- utils::read.table(temp_path,header=FALSE,skip=1,sep="")
    rs_temp[temp_fraction] <- temp_fraction_df$V4
  }
  
  return(rs_temp)
  
}
