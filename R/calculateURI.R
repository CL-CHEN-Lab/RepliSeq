#' @title calculate URI
#' @description Calculate the under replication index of two compared Repli-seq assays
#'
#' @param rs_x a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param rs_y a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#'
#' @return a dataframe composed of genomic coordinates plus sum_x, sum_y, mean_xy and URI : chr,start,stop,sum_x,sum_y,mean_xy,URI
#' @export
#'

calculateURI <- function(rs_x,rs_y) {
  # Control the resolution :
  scale_x <- (rs_x$stop - rs_x$start)[1]
  scale_y <- (rs_y$stop - rs_y$start)[1]
  if (scale_x != scale_y) {
    stop("the Repliseq assays don't have the same resolution. use rescaleRS() before")
  }
  # select fractions :
  temp_x <- rs_x[,!(names(rs_x) %in% c("chr","start","stop"))]
  temp_y <- rs_y[,!(names(rs_y) %in% c("chr","start","stop"))]
  # Calculate mean replication for each condition :
  sum_x <- rowSums(temp_x, na.rm = FALSE)
  sum_y <- rowSums(temp_y, na.rm = FALSE)
  # Calculate difference of means
  temp_diff <- sum_x - sum_y
  # calculate mean of two conditions
  temp_mean <- (sum_x + sum_y) / 2
  # Normalize Difference
  temp_diffNorm <- as.numeric(temp_diff / temp_mean)
  # Standardize using scale function
  temp_diffNormStd <- scale(temp_diffNorm,center=TRUE,scale=TRUE)
  # put results in dataframe and return
  to_return <- data.frame(chr=rs_x$chr,start=rs_x$start,stop=rs_x$stop,
                          sum_x=sum_x,
                          sum_y=sum_y,
                          mean_xy=temp_mean,
                          URI=temp_diffNormStd)
  
  return(to_return)
}