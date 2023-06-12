#' @title rescale Repli-seq assay
#' @description Calculate the count matrices for a lower resolution ( larger genomic windows )
#' 
#' @param rs_assay a dataframe for a Repli-seq assay loaded with readRS()
#' @param scale_factor the factor to apply to the scale (going from 1kb to 50kb gives scale_factor = 50 ; from 50kb to 100kb : scale_factor = 2)
#'
#' @return a dataframe composed of genomic coordinates plus all the fractions from rs_assay as for example : chr,start,stop,S1,S2,S3,S4,S5,S6
#' @import magrittr
#' @importFrom dplyr mutate_each
#' @importFrom dplyr funs
#' @export
#'
#' 


rescaleRS <- function(rs_assay,
                      scale_factor) {
  
  chr <- NULL
  # calculate_inital scale :
  initial_scale <- (rs_assay$stop - rs_assay$start)[1]
  # get names
  rs_names <- names(rs_assay)
  rs_fractions <- rs_names[rs_names != "chr" & rs_names != "start" & rs_names != "stop"]
  # initialize to_return variable
  to_return <- NULL
  rs_chroms <- unique(rs_assay$chr)
  # iterate on chroms :
  for (i in rs_chroms) {
    # select data :
    rs_chr <- subset(rs_assay, chr == i, select = rs_fractions)
    # calculate rs_rescaled :
    rs_chr_cumsums <- rs_chr %>% mutate_each(funs(cumsum))
    rs_chr_smoothed <- as.data.frame(apply(rs_chr_cumsums,2,diff,lag=scale_factor))
    rs_chr_rescaled <- as.data.frame(rs_chr_smoothed[seq(1,length(rs_chr_smoothed[,1,drop=TRUE]),by = scale_factor),])
    # calculate corresponding coordinates :
    chr_chr <- rep(i,times=length(rs_chr_rescaled[,1,drop=TRUE]))
    new_scale = scale_factor * initial_scale
    chr_start <- seq(0,(length(rs_chr_rescaled[,1,drop=TRUE]) * new_scale) - new_scale,by = new_scale)
    chr_stop <- chr_start + new_scale
    # combine rs with coordinates :
    temp_coords_df <- data.frame(chr = chr_chr, start = chr_start, stop = chr_stop)
    chr_all <- cbind.data.frame(temp_coords_df,rs_chr_rescaled)
    colnames(chr_all) <- rs_names
    # append data :
    to_return <- rbind(to_return,chr_all)
  }
  return(to_return)
}