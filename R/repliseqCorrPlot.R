#' @title Calculate and plot correlation plots.
#' @description Calculate and draw a correlation plot from Repli-seq assay (data.frame) see http://www.sthda.com/french/wiki/visualiser-une-matrice-de-correlation-par-un-correlogramme
#'
#' @param rs_assay a Repli-seq assay (data.frame) loaded with readRS() and formatted as chr,start,stop,S1,...,Sx
#' @param corr_method character in either "pearson","spearman" or "kendall" to use as a parameter for cor.test function
#' @param ordering a character in either "original","AOE", "FPC", "hclust", "alphabet" to reorder the compared variables (to fill the corrplot() order parameter )
#' @param plottype from corrplot type argument. either "upper", "lower" or "mixed". default = "mixed".
#' 
#' @return a mixed correlation plot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom corrplot corrplot.mixed
#' @export
#'


repliseqCorrPlot <- function(rs_assay,corr_method = "pearson",ordering = "original",plottype="mixed") {
  temp_rs_assay <- rs_assay[,!(names(rs_assay) %in% c("chr","start","stop"))]
  # ressource : http://www.sthda.com/french/wiki/visualiser-une-matrice-de-correlation-par-un-correlogramme
  # Calculate correlation matrix
  temp.cor <- cor(temp_rs_assay,method = corr_method)
  # Calculate p-values associated with each correlation
  p.mat <- cor.mtest(temp.cor)
  # plot the correlation matrix
  if (type == "mixed") {
    corrplot.mixed(temp.cor,is.corr = FALSE,
                   order=ordering, tl.col="black",addrect = 2,
                   upper.col = brewer.pal(n = 10, name = "PuOr"),
                   lower.col= brewer.pal(n = 10, name = "PuOr"),
                   tl.cex=0.8,number.cex = 0.8,cl.cex = 0.8)
  }
  else {
    corrplot(temp.cor,is.corr = FALSE, type = plottype,
                   order=ordering, tl.col="black",addrect = 2,
                   upper.col = brewer.pal(n = 10, name = "PuOr"),
                   lower.col= brewer.pal(n = 10, name = "PuOr"),
                   tl.cex=0.8,number.cex = 0.8,cl.cex = 0.8)
  }
  
}

###########

# Calculation of the p-values associated with a correlation matrix
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}