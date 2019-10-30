#' @title remove Noise
#' @description Calculate Repli-seq assay count matrices after substracting noise
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param rs_control a Repli-seq control assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S0
#' @param noise_ratios an array with the ratios in the order of rs_assay fractions calculated with calculateNoiseRatios()
#'
#' @return a Repli-seq assay after removing noise (data.frame) and formatted as chr,start,stop,S1,...,Sx
#' @export
#'

removeNoise <- function(rs_assay, rs_control,noise_ratios) {
  # select control and fractions :
  temp_fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))]
  temp_control <- rs_control[,!(names(rs_control) %in% c("chr","start","stop"))]
  # calculate ratios :
  to_return <- rs_assay[,(names(rs_assay) %in% c("chr","start","stop"))]
  for (i in c(1:length(names(temp_fractions)))) {
    temp_noise_removed <- as.numeric(round((rs_assay[,names(temp_fractions)[i]] / noise_ratios[i]) - temp_control, digits  = 3))
    temp_noise_removed[temp_noise_removed < 0] <- 0
    to_return[,names(temp_fractions)[i]] <- temp_noise_removed
  }
  return(to_return)
}