# ------------------------------------------------------------------------------ #
#
# ------------------------------------------------------------------------------ #
cor_wrapper <- function(data, features) {
  library(lattice)

  # Get correlations
  #correlations <- 0
  #for (ind in 1:20) {
  #  correlations <- correlations + 1/20 * cor( complete(data$imp,ind)[cor_items] )
  #}
  correlations <- cor(data$raw[features], use = 'pairwise.complete.obs')
  #cor.test(data$raw[['tci_gen_imp']], data$raw[['tci_fin_imp']])$p.value*100

  # This fn will put correlations on the panels of the matrix heatmap
  customPanel <- function(x, y, z, ...) {
    panel.levelplot(x, y, z, ...)
    panel.text(x, y, round(z, 1))
  }

  # Actual Plot
  levelplot(round(correlations,4), main = 'Pairwise Correlations', scales = list(x = list(draw = FALSE)),
            xlab = '', ylab = '', at = seq(-1, 1, length.out = 100), panel = customPanel)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Wrapper for density plots
# ------------------------------------------------------------------------------ #
density_wrapper <- function(data, names, title) {
  library(ggplot2)
  library(reshape2)

  ggplot( melt(data.frame(data$raw[names]), id.vars = NULL), aes(x = value)) +
    facet_wrap(~variable, scales = "free") +
    geom_histogram(na.rm = TRUE) +
    ggtitle(title) # + theme_bw()
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Wrapper for ggpairs
# ------------------------------------------------------------------------------ #
ggpairs_wrap <- function(data) {
  library(ggplot2)
  library(GGally)
  ggpairs(data, lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.75)))
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#
# ------------------------------------------------------------------------------ #
cca_wrapper <- function(dset1, dset2){

    # Setup:
    dset1_names <- colnames(dset1)
    dset2_names <- colnames(dset2)

    library(CCA)
    cca_res <- cancor(dset1, dset2)
    nxcoef  <- length(cca_res$xcoef[,1])
    nycoef  <- length(cca_res$ycoef[,1])

    # Set plot options
    op <- par(mar = c(7.5, 4.1, 4.1, 2.1))
    #options(repr.plot.width = 9, repr.plot.height = 9)

    # One figure in row 1 and two figures in row 2
    layout(matrix(c(1, 2, 2, 2, 1, 3, 3, 3), 2, 4, byrow = TRUE))

    # Significance testing:
    library(CCP)
    n_obs         <- dim(dset1)[1]
    n_task_vars   <- dim(dset1)[2]
    n_survey_vars <- dim(dset2)[2]

    ps <- p.asym(rho = cca_res$cor, n_obs, n_task_vars, n_survey_vars, tstat = "Wilks")
    nlines <- max(sum(ps$p.value < 0.05), 1)
    if (nlines < 1) stop('nlines < 1 (no significant correlations)')
    colors <- rainbow(nlines)


    # U-V Correlation Plot:
    print(cca_res$cor)
    #plot(cca_res$cor,type="b", main='Task-Survey Cannonical Correlation Results',
    #     xlab='Canonical Covariate Index', ylab='U-V Correlation')
    print('P-values')
    print(ps$p.value[1:nlines])
    barplot(cca_res$cor[1:nlines], col = colors[1:nlines], main = 'Canonical Correlations')
    #abline(h = 0)

    # U Coefficients Plot:
    #plotchar <- seq(18,18+nlines,1)

    #xrange <- range(1:length(dset1_names))
    #yrange <- range(cca_res$xcoef[,1:nlines])

    #op <- par(mar=c(12.1, 4.1, 4.1, 2.1))

    #plot(xrange, yrange, type="n", xlab="", ylab="Canonical Coefficient", xaxt="n")

    scale_fun <- function(x) {x*1/apply(apply(cca_res$xcoef, 2, abs ), 2, max)}
    scaled_xc <- apply(cca_res$xcoef, 1, scale_fun)
    scaled_xc <- cca_res$xcoef

    print(scaled_xc)
    barplot(t(scaled_xc[,1:nlines]), beside = TRUE, col = colors[1:nlines], scales = list(x = list('',rot = 90)), xaxt = 'n')
    abline(h = 0)
    linetype <- c(1:nlines)

    #for (i in 1:nlines) {
    #  lines(1:nxcoef, cca_res$xcoef[,i], type="b", lwd=1.5, lty=linetype[1], col=colors[i], pch=plotchar[i])
    #}

    # add a title and subtitle
    title("CCA Task Coefficients")

    # add a legend
    legend(1, 50, 1:nlines, cex = 0.8, col = colors, lty = linetype)
    axis(1, labels = dset1_names, at = seq(2, (nlines+1)*length(dset1_names)+1, (nlines+1)), las = 2)
    grid()

    ## V Coefficients Plot:
    plotchar <- seq(18,18 + nlines,1)

    xrange <- range(1:length(dset2_names))
    yrange <- range(cca_res$ycoef[,1:nlines])

    #op <- par(mar=c(10.1, 4.1, 4.1, 2.1))
    #plot(xrange, yrange, type="n", xlab="", ylab="Canonical Coefficient", xaxt="n")
    #abline(h = 0)

    scale_fun <- function(x) {x*1/apply(apply(cca_res$ycoef, 2, abs ), 2, max)}
    scaled_yc <- apply(cca_res$ycoef, 1, scale_fun)
    scaled_yc <- cca_res$ycoef

    print(scaled_yc)
    barplot(t(scaled_yc[,1:nlines]), beside = TRUE, col = colors[1:nlines], scales = list(x = list('',rot = 90)), xaxt = 'n')
    axis(1, labels = dset2_names, at = seq(2, (nlines+1)*length(dset2_names) + 1,(nlines+1)), las = 2)
    abline(h = 0)
    #for (i in 1:nlines) {
    #  lines(1:nycoef, cca_res$ycoef[,i], type="b", lwd=1.5, lty=linetype[1], col=colors[i], pch=plotchar[i])
    #}

    # add a title and subtitle
    title("CCA Survey Coefficients")

    # add a legend
    legend(xrange[1], yrange[2], 1:nlines, cex=0.8, col=colors, lty=linetype)
    #axis(1, labels = dset2_names, at = 1:nycoef, las=2)
    grid()

    return(cca_res)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#
# ------------------------------------------------------------------------------ #
pca_wrapper <- function(data) {
  pca_res  <- prcomp(data, center = TRUE, scale. = TRUE)
  pca_pvar <- pca_res$sdev^2/sum(pca_res$sdev^2)

  ### PLOT THE VARIATION CAPTURED BY PCAs:
  op <- par(mar = c(4, 4, 4, 1))
  layout(1)
  #layout(matrix(c(1, 2), 2, 1, byrow = TRUE))
  #par(mfrow=c(1,3), mai = c(0.1, 0.5, 0.1, 0.3))
  #options(repr.plot.width = 9, repr.plot.height = 3)

  ylab_top <- "Proportion of Variance [%]"
  main <- 'PCA Variance Explained, for fMRI Stop - Go Contrast'
  # Task PCA Plot:
  barplot(pca_pvar,  col = c('red'), las = 2, cex.axis = 0.7, xlab = c('Principle Components'), ylab = ylab_top, main = main)

  ## PLOT THE LOADINGS FOR EACH PCA:
  #par(mfrow=c(1,3), mai = c(0.1, 0.5, 0.1, 0.3))
  options(repr.plot.width = 5, repr.plot.height = 5)

  n_comp <- 3
  n_cols <- dim(data)[2]
  barplot(t(pca_res$rotation[, 1:n_comp]), main = "PCA Weights, for fMRI Stop - Go Contrast", horiz = FALSE, beside = TRUE, col = rainbow(3),  scales = list(x = list('', rot = 90)))
  axis(2, labels = colnames(data), at = seq(1 + n_comp/2, n_cols*(1 + n_comp), by = (n_comp + 1)), las = 2)
  legend("topright", legend = c('PC1', 'PC2', 'PC3'), fill = rainbow(3))

  # project new data onto the PCA space
  new_df <- scale(data, pca_res$center, pca_res$scale) %*% pca_res$rotation
}
# ------------------------------------------------------------------------------ #
