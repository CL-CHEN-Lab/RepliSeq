#' @title Normalize Repli-seq assay
#' @description Calculate Repli-seq assay count matrices after normalizing (dividing counts by ratios)
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param norm_ratios an array with the ratios in the order of rs_assay fractions
#'
#' @return a Repli-seq assay after normalizing (dividing fractions by ratios) as a data.frame and formatted as chr,start,stop,S1,...,Sx
#' @export
#'

normalizeRS <- function(rs_assay, norm_ratios) {
  # select fractions :
  temp_fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))]
  to_return <- rs_assay[,(names(rs_assay) %in% c("chr","start","stop"))]
  # proceed to normalization :
  for (i in c(1:length(names(temp_fractions)))) {
    temp_normalized <- as.numeric(round((rs_assay[,names(temp_fractions)[i]] / norm_ratios[i]),digits = 3))
    to_return[,names(temp_fractions)[i]] <- temp_normalized
  }
  return(to_return)
}