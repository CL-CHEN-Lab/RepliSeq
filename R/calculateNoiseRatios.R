#' @title Calculate noise ratios
#' @description Calculate noise ratios of a Repli-seq assay versus a Repli-seq control
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param rs_control a Repli-seq control assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S0
#' @param noise_method either "max_density"(default) or "min_downslope"
#' @return an array with the ratios in the order of rs_assay fractions
#' 
#' @importFrom stats density
#' 
#' @export
#'


calculateNoiseRatios <- function(rs_assay, rs_control, noise_method = "max_density") {
  # select control and fractions :
  temp_fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))]
  temp_control <- rs_control[,!(names(rs_control) %in% c("chr","start","stop"))]
  # calculate ratios :
  to_return <- NULL
  for (fraction in names(temp_fractions)) {
    temp_fraction_vs_control <- round(rs_assay[,fraction] / temp_control, digits = 3)
    temp_density <- stats::density(temp_fraction_vs_control,na.rm = TRUE)
    if (noise_method == "min_downslope") {
      temp_downslopes <- which(diff(temp_density$y) < 0)
      temp_ratio <- temp_density$x[min(temp_downslopes)]
    }
    else if (noise_method == "max_density") {
      temp_max <- which.max(temp_density$y)
      temp_ratio <- round(temp_density$x[temp_max],digits = 3)
    }
    to_return <- c(to_return,temp_ratio)
  }
  return(to_return)
}