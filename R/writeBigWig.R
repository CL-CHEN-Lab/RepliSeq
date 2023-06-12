#' @title Write BigWig files of Repli-seq assays
#' @description writes one BigWig file per fraction in the provided Repli-seq assay
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx (or S0)
#' @param path_file a path to the files to write
#' @param sample_name a string for the file names
#' @param chromsizes path to the chromsizes file to pass to WigToBigWig as argument
#' @param wiggle_start a numeric value for wiggle header
#' @param wiggle_step a numeric value for wiggle header
#' @param wiggle_span a numeric value for wiggle header
#' @param system_separator default is "/" as dedicated for linux file system
#'
#' @importFrom utils write.table
#'
#' @return NULL
#' @export
#' 

writeBigwig <- function(rs_assay,path_file,sample_name,chromsizes,wiggle_start,wiggle_step,wiggle_span,system_separator = "/") {
  # select fractions :
  fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop")),drop = FALSE]
  # iterate over fractions :
  for (i in names(fractions)) {
    temp_fraction <- rs_assay[,i]
    to_write <- NULL
    # iterate over chromosomes :
    for (j in levels(rs_assay$chr)) {
      temp_header <- paste0("fixedStep chrom=",j," start=",wiggle_start," step=",wiggle_step," span=",wiggle_span)
      temp_fraction_chr <- temp_fraction[rs_assay$chr == j]
      to_write <- c(to_write,temp_header,temp_fraction_chr)
    }
    # write wiggle ; convert to BigWig and rm wiggle :
    file_name <- paste0(path_file,system_separator,sample_name,"-",i,".wig")
    utils::write.table(to_write,file = file_name,quote = FALSE, col.names = temp_header,row.names = FALSE)
    system(paste0("wigToBigWig -clip ",file_name," ",chromsizes," ",file_name,".bw"))
    system(paste0("rm ",file_name))
  }
}

