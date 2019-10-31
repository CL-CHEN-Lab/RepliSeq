#' @title double the values of X chr (for human male samples)
#' @description Calculate Repli-seq assay count matrices after doubling the values of chrX windows
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx (or S0)
#' @param chr level of chr to double the value (default is set to "chrX")
#'
#' @return a Repli-seq assay doubling values of chrX windows (data.frame) and formatted as chr,start,stop,S1,...,Sx
#' @export
#' 

doubleXchr <- function(rs_assay,chr = "chrX") {
  # select fractions and genomic coords :
  fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop")),drop = FALSE]
  to_return <- rs_assay[,names(rs_assay) %in% c("chr","start","stop")]
  for (i in names(fractions)) {
    # select fraction :
    temp_fraction <- fractions[,i]
    # selct indexes of chr :
    temp_fraction[rs_assay$chr == chr] <- 2 * temp_fraction[rs_assay$chr == chr]
    # add fraction to to_return :
    to_return[i] <- temp_fraction
  }
  return(to_return)
}