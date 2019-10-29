#' @title smooth Repli-seq assay
#' @description Calculate the smoothed count matrices of a Repli-seq assay
#' 
#' @param rs_assay a dataframe for a Repli-seq assay loaded with readRS()
#' @param smooth_factor the factor to apply to the scale (going from 1kb to 50kb gives smooth_factor = 50 ; from 50kb to 100kb : smooth_factor = 2)
#'
#' @return a dataframe composed of genomic coordinates plus all the fractions from rs_assay as for example : chr,start,stop,S1,S2,S3,S4,S5,S6
#' @import magrittr
#' @importFrom dplyr mutate_each
#' @importFrom dplyr funs
#' @export
#'
#' 


smoothRS <- function(rs_assay,smooth_factor) {
  # calculate_inital scale :
  initial_scale <- (rs_assay$stop - rs_assay$start)[1]
  # get names
  rs_names <- names(rs_assay)
  rs_fractions <- rs_names[rs_names != "chr" & rs_names != "start" & rs_names != "stop"]
  # initialize to_return variable
  to_return <- NULL
  rs_chroms <- levels(unique(rs_assay$chr))
  # iterate on chroms :
  for (i in rs_chroms) {
    # select data :
    rs_chr <- subset(rs_assay, chr == i, select = rs_fractions)
    # calculate rs_rescaled :
    rs_chr_cumsums <- rs_chr %>% mutate_each(funs(cumsum))
    rs_chr_smoothed <- as.data.frame(apply(rs_chr_cumsums,2,diff,lag=smooth_factor))
    # select coordinates :
    temp_coords_df <- subset(rs_assay,chr == i, select = c("chr","start","stop"))
    ## need to remove last values from coordinates as diff functions did
    temp_coords_df <- temp_coords_df[(1:length(rs_chr_smoothed[,1])),]
    # combine rs with coordinates :
    chr_all <- cbind.data.frame(temp_coords_df,rs_chr_smoothed)
    colnames(chr_all) <- rs_names
    # append data :
    to_return <- rbind(to_return,chr_all)
  }
  return(to_return)
}