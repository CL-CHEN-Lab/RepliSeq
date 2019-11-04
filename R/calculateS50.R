#' @title Calculate S50.
#' @description Calculate replication timing S50 values from Repli-seq assay (data.frame)
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#'
#' @return a dataframe composed of genomic coordinates plus S50 : chr,start,stop,S50
#' @export
#'


calculateS50 <- function(rs_assay) {
  fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))]
  # extract fractions as matrix :
  assay <- as.matrix(fractions)
  # calculate S50 with function s50
  temp_s50 <- apply(assay,1,s50)
  # concatenate genomic coordinates and calculated s50 in data frame to return
  temp <- data.frame(chr = rs_assay$chr, start = rs_assay$start, stop = rs_assay$stop,S50 = round(temp_s50,digits = 3))
}

###########

s50 <- function(rs_line,scale_value = 1000) {
  if (sum(rs_line) == 0) {
    ## No replication at all, then should be NA ?
    to_return <- 0
  } else {
    ## There's replication so calculate the S50 value.
    scaled_line <- NULL
    for (i in c(1:length(rs_line))) {
      scaled_line <- c(scaled_line,rep(rs_line[i],times = scale_value))
    }
    rs_sum <- cumsum(scaled_line)
    rs_pos <- min(which(rs_sum>=(sum(scaled_line)/2)))
    to_return  <- rs_pos/length(scaled_line)
  }
}