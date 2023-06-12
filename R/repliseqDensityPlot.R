#' @title Calculate and plot density plots.
#' @description Calculate and draw a density plot from Repli-seq assay (data.frame)
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param grided boolean (default = F) if True then each fraction s displayed on a different panel, otherwise all the fractions are colored in a single panel
#' @param stat "counts" (default) or "density" to use in the stat argument
#' @param ncols (default = 3) amount of columns when grided == T (see ggplot2 facet_wrap ncol parameter) 
#' @param rs_resolution resolution of the rs_assay to display in the labels in the plots. Default equals "50kb"
#'
#' @return a ggplot object displaying density for the provided repliseq assay
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#'


repliseqDensityPlot <- function (rs_assay,
                                 stat = "counts",
                                 grided = F,
                                 ncols = 3,
                                 rs_resolution = "50kb") {
  
  value <- variable <- NULL
  # extract and melt the fractions :
  temp_fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))]
  temp_fractions_melted <- melt(temp_fractions)
  
  
  label_x <- paste0("Counts of reads / ",rs_resolution, "windows")
  label_y <- paste0(stat," of windows")
  
  # grided plot :
  if (grided) {
    p <- ggplot(temp_fractions_melted,aes(x = value)) +
      geom_line(stat=stat) +
      facet_wrap(~ variable,ncol=ncols) +
      theme_bw(20) +
      labs(x = label_x, y = label_y)
  }
  
  # not grided plot :
  else {
    p <- ggplot(temp_fractions_melted,aes(x = value,color=variable)) +
      geom_line(stat=stat) +
      theme_bw(20) +
      labs(x = label_x, y = label_y ,color="sample")
  }
  
  return(p)
}