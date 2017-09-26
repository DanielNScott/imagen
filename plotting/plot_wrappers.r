library(viridis)

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
ggpairs_wrap <- function(data, title, names = NULL) {
  library(ggplot2)
  library(GGally)
  if (!is.null(names)) {
    colnames(data) <- names
  }
  ggpairs(data, lower = list(continuous = wrap("smooth", alpha = 0.3, size = 0.75))) +
    ggtitle(title)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
cca_wrapper <- function(dset1, dset2, title1, title2, title3){

  # Setup:
  dset1_names <- colnames(dset1)
  dset2_names <- colnames(dset2)

  library(CCA)
  #PMA::CCA(dset1, dset2)
  #browser()
  cca_res <- cancor(dset1, dset2)

  library(CCP)
  n_obs         <- dim(dset1)[1]
  n_task_vars   <- dim(dset1)[2]
  n_survey_vars <- dim(dset2)[2]

  ps <- p.asym(rho = cca_res$cor, n_obs, n_task_vars, n_survey_vars, tstat = "Wilks")
  nlines <- sum(ps$p.value < 0.05)

  # Task PCA Plot:
  plot <- ggplot(data.frame('index' = 1:length(cca_res$cor), 'value' = cca_res$cor), aes(x = index, y = value)) +
    geom_bar(stat = 'identity') + xlab('Canonical Covariate Index') + ylab('Correlation') +
    ggtitle(title1)
  print(plot)

  scale_fun <- function(x) {x*1/apply(apply(cca_res$xcoef, 2, abs ), 2, max)}
  scaled_xc <- apply(cca_res$xcoef, 1, scale_fun)
  scaled_xc <- cca_res$xcoef

  browser()
  colnames(scaled_xc) <- c('Cov. 1', 'Cov. 2', 'Cov. 3')
  plot <- ggplot(data = melt(t(scaled_xc)), aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat = "identity", position = 'dodge') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('Feature') + ylab('Weight') + guides(fill = guide_legend(title = NULL)) +
    ggtitle(title2)
  print(plot)

  scale_fun <- function(x) {x*1/apply(apply(cca_res$ycoef, 2, abs ), 2, max)}
  scaled_yc <- apply(cca_res$ycoef, 1, scale_fun)
  scaled_yc <- cca_res$ycoef
  colnames(scaled_yc) <- c('Cov. 1', 'Cov. 2', 'Cov. 3', 'Cov. 4')

  plot <- ggplot(data = melt(t(scaled_yc[,1:3])), aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat = "identity", position = 'dodge') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('Feature') + ylab('Weight') + guides(fill = guide_legend(title = NULL)) +
    ggtitle(title3)
  print(plot)

  print('P-values')
  print(ps$p.value[1:nlines])

  return(cca_res)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
pca_wrapper <- function(data, title, max_comp = NULL, subset = NULL) {
  cmpl <- complete.cases(data)
  data <- data[cmpl,]

  pca_res  <- prcomp(data, center = TRUE, scale. = TRUE)
  pca_pvar <- pca_res$sdev^2/sum(pca_res$sdev^2)

  if (is.null(max_comp)) {max_comp <- 3}
  if (is.null(subset  )) {subset   <- 1:max_comp}

  # Task PCA Plot:
  plot <- ggplot(data.frame('index' = subset, 'value' = pca_pvar[subset]), aes(x = index, y = value)) +
    geom_bar(stat = 'identity') + xlab('PCA Index') + ylab('Proportion of Variance') +
    ggtitle(title)
  print(plot)

  plot <- ggplot(data = melt(t(pca_res$rotation[,subset])), aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('Feature') + ylab('Weight') + guides(fill = guide_legend(title = NULL)) +
    ggtitle(title)
  print(plot)

  # project new data onto the PCA space
  proj <- scale(data, pca_res$center, pca_res$scale) %*% pca_res$rotation
  return(list(proj = proj, mask = cmpl))
}
# ------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------ #
tsne_wrapper <- function(data, color_by, title){

  # Get basic result
  res_full <- Rtsne::Rtsne(data)
  res      <- res_full$Y

  # Append a color column
  res  <- cbind(res, data[color_by])
  colnames(res) <- c('x','y', 'z')

  ggplot(data.frame(res), aes(x = x, y = y, colour = z)) +
    geom_point() + scale_color_viridis(option = 'plasma') +
    xlab('Dimension 1 [A.U.]') + ylab('Dimension 2 [A.U.]') +
    labs(colour = color_by) + ggtitle(title)

}
# ------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------ #
color_point_wrap <- function(data, title, xlabel = NULL, ylabel = NULL, zlabel = NULL){

  var_names <- colnames(data)
  if (is.null(xlabel)) {xlabel <- var_names[1]}
  if (is.null(ylabel)) {ylabel <- var_names[2]}
  if (is.null(zlabel)) {zlabel <- var_names[3]}

  colnames(data) <- c('x','y', 'z')
  ggplot(data, aes(x = x, y = y, colour = z)) +
    geom_point() + scale_color_viridis(option = 'plasma') +
    xlab(xlabel) + ylab(ylabel) +
    labs(colour = zlabel) + ggtitle(title)
}
# ------------------------------------------------------------------------------ #



