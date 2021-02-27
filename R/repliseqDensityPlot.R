#' @title Calculate and plot density plots.
#' @description Calculate and draw a density plot from Repli-seq assay (data.frame)
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param grided boolean (default = F) if True then each fraction s displayed on a different panel, otherwise all the fractions are colored in a single panel
#' @param ncols (default = 3) amount of columns when grided == T (see ggplot2 facet_wrap ncol parameter) 
#'
#' @return a ggplot object displaying density for the provided repliseq assay
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#'


repliseqDensityPlot <- function (rs_assay,grided = F,ncols = 3) {
  
  # extract and melt the fractions :
  temp_fractions <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))]
  temp_fractions_melted <- melt(temp_fractions)
  
  # grided plot :
  if (grided) {
    p <- ggplot(temp_fractions_melted,aes(x = value)) +
      geom_line(stat="density") +
      facet_wrap(~ variable,ncol=ncols) +
      theme_bw() +
      labs(x = "")
  }
  
  # not grided plot :
  else {
    p <- ggplot(temp_fractions_melted,aes(x = value,color=variable)) +
      geom_line(stat="density") +
      theme_bw() +
      labs(x = "",color="sample")
  }
  
  return(p)
}