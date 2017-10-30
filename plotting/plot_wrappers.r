library(viridis)
library(ggplot2)
library(reshape2)
library(GGally)
library(ggthemr)
ggthemr('fresh')

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
level_wrapper <- function(data, title, ylabel = '', legend = 'value') {

  ggplot(data, aes(x = Var1, y = Var2, z = value)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(label = round(value, 1)), color = 'orange', size = 7) +
    scale_fill_gradientn(colours = viridis(10, option = 'viridis'), limits = c(-1,1)) +
    guides(fill = guide_legend(title = legend)) +
    ggtitle(title) + ylab(ylabel) +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())

}
# ------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------ #
# Wrapper for density plots
# ------------------------------------------------------------------------------ #
density_wrapper <- function(data, names, title, xlabel = NULL, ylabel = NULL, facet_labels = NULL) {
  library(ggplot2)
  library(reshape2)

  facet_labeller <- function(variable, value){
    return(facet_labels[value])
  }

  plot <- ggplot( melt(data[names], id.vars = NULL), aes(x = value)) +
    geom_density(na.rm = TRUE) +
    xlab(xlabel) + ylab(ylabel) +
    theme( strip.text = element_text(margin = margin(b = 20))) +
    ggtitle(title)

  if (!is.null(facet_labels)) {
    plot <- plot + facet_wrap(~variable, scales = "free", labeller = facet_labeller)
  } else {
    plot <- plot + facet_wrap(~variable, scales = 'free')
  }

  print(plot)
    # + theme_bw()
  #theme(plot.title = element_text(vjust=1.5, face="bold", size = 20),
    #      axis.title.x = element_blank(), axis.title.y = element_blank())
}
# ------------------------------------------------------------------------------ #

#----
bar_wrap <- function(data, colnames = NULL, title = NULL ){
  #data <- data.frame(data, colnames)
  #colnames(data) <- c('value', 'names')

  plot <- ggplot(data, aes(x = names, y = mean)) +
    geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab('MID Model Term') + ylab('Regression Coefficient') +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), colour = 'black', width = 0.2) +
    scale_x_discrete(labels = c('(Intercept)', 'Location','Baseline RT.','Var in Duration\nModulation','Var in Reward\nModulation')) +
    ggtitle(title)
  print(plot)
}
#----


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
cca_wrapper <- function(dset1, dset2, title1, title2, title3, subset = NULL,
                        max_comp = NULL, yscale = NULL, flip = FALSE, xticks1 = NULL,
                        xticks2 = NULL){

  # Setup:
  dset1_names <- colnames(dset1)
  dset2_names <- colnames(dset2)

  if (is.null(max_comp)) {max_comp <- 3}
  if (is.null(subset  )) {subset   <- 1:max_comp}

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

  # Correlation Plot
  plot <- ggplot(data.frame('index' = subset, 'value' = cca_res$cor[subset]), aes(x = index, y = value)) +
    geom_bar(stat = 'identity') + xlab('Canonical Covariate Index') + ylab('Correlation') +
    ggtitle(title1)

  if (!is.null(yscale)) {plot <- plot + ylim(yscale)}
  print(plot)

  # Setup for Canonical Covar U plot
  scale_fun <- function(x) {x*1/apply(apply(cca_res$xcoef, 2, abs ), 2, max)}
  scaled_xc <- apply(cca_res$xcoef, 1, scale_fun)
  scaled_xc <- cca_res$xcoef
  if (flip) {scaled_xc <- -scaled_xc}

  # Plot Canonical Covar U
  if (is.null(xticks1)) {xticks1 <- colnames(dset1)}
  colnames(scaled_xc) <- paste('Cov.', 1:dim(scaled_xc)[2], sep = '')
  scaled_xc <- scaled_xc[,subset]
  if (length(subset) > 1) {
    plot <- ggplot(data = melt(t(scaled_xc)), aes(x = Var2, y = value, fill = Var1))
  } else {
    plot <- ggplot(data = melt(t(scaled_xc)), aes(x = Var2, y = value))
  }
  plot <- plot + geom_bar(stat = "identity", position = 'dodge') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab(NULL) + ylab('Weight') + guides(fill = guide_legend(title = NULL)) +
    geom_hline(yintercept = 0) +
    scale_x_discrete(labels = xticks1) +
    ggtitle(title2)

  if (length(subset == 1)) {plot <- plot + theme(legend.position = "none")}
  print(plot)

  # Setup for Canonical Covariate V
  scale_fun <- function(x) {x*1/apply(apply(cca_res$ycoef, 2, abs ), 2, max)}
  scaled_yc <- apply(cca_res$ycoef, 1, scale_fun)
  scaled_yc <- cca_res$ycoef
  if (flip) {scaled_yc <- -scaled_yc}
  colnames(scaled_yc) <- paste('Cov.', 1:dim(scaled_yc)[2], sep = '')

  # Plot Canonical Covariate V
  if (is.null(xticks2)) {xticks2 <- colnames(dset2)}
  scaled_yc <- scaled_yc[,subset]
  if (length(subset) > 1) {
    plot <- ggplot(data = melt(t(scaled_yc)), aes(x = Var2, y = value, fill = Var1))
  } else {
    plot <- ggplot(data = melt(t(scaled_yc)), aes(x = Var2, y = value))
  }
  plot <- plot + geom_bar(stat = "identity", position = 'dodge') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab(NULL) + ylab('Weight') + guides(fill = guide_legend(title = NULL)) +
    geom_hline(yintercept = 0) +
    scale_x_discrete(labels = xticks2) +
    ggtitle(title3)

  if (length(subset == 1)) {plot <- plot + theme(legend.position = "none")}
  print(plot)

  print('P-values')
  print(ps$p.value[1:nlines])

  return(cca_res)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
pca_wrapper <- function(data, title, max_comp = NULL, subset = NULL, xticks = NULL) {
  cmpl <- complete.cases(data)
  data <- data[cmpl,]

  pca_res  <- prcomp(data, center = TRUE, scale. = TRUE)
  pca_pvar <- pca_res$sdev^2/sum(pca_res$sdev^2)

  if (is.null(max_comp)) {max_comp <- 3}
  if (is.null(subset  )) {subset   <- 1:max_comp}
  if (is.null(xticks  )) {xticks   <- colnames(data)}

  # Task PCA Plot:
  plot <- ggplot(data.frame('index' = subset, 'value' = pca_pvar[subset]), aes(x = index, y = value)) +
    geom_bar(stat = 'identity') + xlab('PCA Index') + ylab('Proportion of Variance') +
    ggtitle(title)
  print(plot)

  plot <- ggplot(data = melt(t(pca_res$rotation[,subset])), aes(x = Var2, y = value, fill = Var1)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = xticks) +
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



