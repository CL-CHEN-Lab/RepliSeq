#' @title Normalize Repli-seq assay
#' @description Calculate Repli-seq assay count matrices after normalizing (dividing counts by ratios)
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param norm_ratios an array with the ratios in the order of rs_assay fractions
#' @param round_digits amount of digits after comma to conserve in the output when rounding (default = 3)
#'
#' @return a Repli-seq assay after normalizing (dividing fractions by ratios) as a data.frame and formatted as chr,start,stop,S1,...,Sx
#' @export
#'

normalizeRS <- function(rs_assay, norm_ratios,round_digits = 3) {
  # select fractions :
  temp_fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))]
  to_return <- rs_assay[,(names(rs_assay) %in% c("chr","start","stop"))]
  # proceed to normalization :
  for (i in c(1:length(names(temp_fractions)))) {
    ratio_to_divide <- as.double(norm_ratios[i][[1]])
    list_to_divide <- rs_assay[,names(temp_fractions)[i]]
    temp_normalized <- list_to_divide / ratio_to_divide
    temp_normalized <- round(temp_normalized,digits = round_digits)
    to_return[,names(temp_fractions)[i]] <- temp_normalized
  }
  return(to_return)
}