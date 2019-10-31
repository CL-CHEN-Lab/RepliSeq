#' @title Write Bedgraph files of Repli-seq assays
#' @description writes one bedgraph file per fraction in the provided Repli-seq assay
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx (or S0)
#' @param path_file a path to the files to write
#' @param sample_name a string for the file names
#' @param system_separator default is "/" as dedicated for linux file system
#'
#' @return NULL
#' @export
#'

writeBedgraph <- function(rs_assay,path_file,sample_name,system_separator = "/") {
  # select fractions and genomic coords :
  fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop")),drop = FALSE]
  temp_coords <- rs_assay[,names(rs_assay) %in% c("chr","start","stop")]
  
  for (i in names(fractions)) {
    header <- c("track ","type=bedGraph ","name=RepliSeq",paste0("description=",sample_name))
    temp_fraction <- rs_assay[,i,drop=FALSE]
    to_write <- temp_coords
    to_write[i] <- temp_fraction
    file_name <- paste0(path_file,system_separator,sample_name,"-",i,".bdg")
    write.table(to_write,file = file_name,quote = FALSE, col.names = header,row.names = FALSE)
  }
}