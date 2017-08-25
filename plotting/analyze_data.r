# ------------------------------------------------------------------------------ #
plot_cors <- function(data, features) {
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
    #panel.text(x, y, round(z, 1))
  }

  # Actual Plot
  levelplot(round(correlations,4), main = 'Pairwise Correlations', scales = list(x = list(draw = FALSE)),
            xlab = '', ylab = '', at = seq(-1, 1, length.out = 100), panel = customPanel)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Wrapper for density plots
# ------------------------------------------------------------------------------ #
density_plot <- function(data, names, title) {
  library(ggplot2)
  library(reshape2)

  ggplot( melt(data.frame(data$raw[names]), id.vars = NULL), aes(x = value)) +
  facet_grid(~variable, scales = "free") + geom_histogram(na.rm = TRUE) +
  theme_bw() + ggtitle(title)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Function to call all the plotting
# ------------------------------------------------------------------------------ #
analysis <- function(data) {

  # First do diagnostics
  # --- some code will go here ---
  # Messy, comprehensive histograms, sploms, correlations, etc of everything
  # ------------------------------

  # Training vs test set parameter recovery
  d <- melt( data$raw[c(data$names$sst, 'set')], id.vars = 'set')
  ggplot(d, aes(x = value, fill = set)) +
    facet_wrap(~variable, scales = 'free') +
    geom_histogram(alpha = 0.2, position = 'identity') +
    ggtitle('SST Params by Train / Test Set')

  # Task - FMRI CCA
  feats_left <- c('mid_rew', 'mid_high_rew', 'SSRTGo', 'SSRTStop', 'shift_go', data$names$CANTAB, setdiff(data$names$Age, 'bmi'))
  std_fmri_feats <- as.vector(outer(rois, sfx, function(x,y) {paste(x,y, sep = '')}))
  cmpl <- complete.cases(data$scores[c(feats_left, std_fmri_feats)])
  cca_wrapper(data$scores[feats_left][cmpl,], data$scores[std_fmri_feats][cmpl,])

  ### Start plotting things that qualify as 'results'...

  ###------------------###
  ### HISTOGRAMS OF... ###
  ###------------------###
  # ... impulsivity stuff
  feats <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  density_plot(data, feats, 'Impulsivity Surveys')

  # ... mid params
  feats <- c('mid_rew', 'mid_high_rew', 'mid_rew_var', 'mid_high_rew_var', 'shift_go', 'SSRTGo', 'SSRTStop')
  density_plot(data, feats, 'Model Parameters')

  # ... sst params
  d <- melt( data$raw[c(setdiff(data$names$sst, 'pf_stop'), 'set')], id.vars = 'set')
  ggplot(d, aes(x = value, fill = set)) +
    facet_wrap(~variable, scales = 'free') +
    geom_histogram(alpha = 0.2, position = 'identity') +
    ggtitle('SST Params by Train / Test Set')

  # ... standard fmri weights
  rois <- c('rPreSMA', 'rIFG', 'rCaudate', 'rSTN', 'rGPe', 'rGPi', 'rThalamus')
  titles <- c('SST Go Weights', 'SST Stop Weights', 'SST Stop Respond Weights',
              'SST Stop > Go Contrasts', 'SST Stop > Stop Respond Contrasts')
  sfx <- c('_go', '_st', '_sr', '_st_go', '_st_sr')
  for (i in 1:length(sfx)) {
    std_fmri_feats <- paste(rois, sfx[i], sep = '')
    print(density_plot(data, std_fmri_feats, titles[i]))
  }

  # ... ancestral graph fmri weights
  rois <- c('presma.ifg', 'presma..stn', 'ifg..stn', 'presma..caud', 'ifg..caud',
            'caud..gpe', 'gpe..gpi', 'stn..gpi', 'gpi..thalamus')
  titles <- c('SST Stop Conn.', 'SST Stop Respond Conn', 'SST Stop, Stop Respond Contrasts')
  sfx <- c('_st', '_sr', '_st_sr')
  for (i in 1:length(sfx)) {
    feats <- paste(rois, sfx[i], sep = '')
    print(density_plot(data, feats, titles[i]))
  }


  ### TASK - SURVEY ###
  # Pairwise correlations
  feats <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount',
              'adhd_teacher', 'adhd_parent','adhd_child', 'nic_use', 'alc_use',
              'alc_regret', 'mid_rew', 'mid_high_rew', 'SSRTGo', 'SSRTStop', 'shift_go',
             data$names$CANTAB, data$names$Age)
  plot_cors(data, feats)

  # SPLOMS (Scatter Plot Matrices)
  feats <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount',
             'mid_rew', 'mid_high_rew', 'SSRTGo', 'SSRTStop', 'shift_go')
  splom(data$raw[feats], pscales = 0)

  # Task - Survey CCA
  feats_left <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount',
                  'adhd_teacher', 'adhd_parent','adhd_child', 'nic_use', 'alc_use', 'alc_regret', 'bmi')
  feats_right <- c('mid_rew', 'mid_high_rew', 'SSRTGo', 'SSRTStop', 'shift_go', data$names$CANTAB, setdiff(data$names$Age, 'bmi'))
  cmpl <- complete.cases(data$scores[c(feats_left, feats_right)])
  cca_wrapper(data$scores[feats_left][cmpl,], data$scores[feats_right][cmpl,])



  ### TASK - FMRI ###
  # Pairwise correlations
  feats <- c('mid_rew', 'mid_high_rew', 'SSRTGo', 'SSRTStop', 'shift_go', data$names$CANTAB, setdiff(data$names$Age, 'bmi'),
             std_fmri_feats)
  plot_cors(data, feats)

  # Task - FMRI CCA
  feats_left <- c('mid_rew', 'mid_high_rew', 'SSRTGo', 'SSRTStop', 'shift_go', data$names$CANTAB, setdiff(data$names$Age, 'bmi'))
  std_fmri_feats <- as.vector(outer(rois, sfx, function(x,y) {paste(x,y, sep = '')}))
  cmpl <- complete.cases(data$scores[c(feats_left, 'rSTN_st_go', 'rCaudate_st_go')])
  cca_wrapper(data$scores[feats_left][cmpl,], data$scores[c('rSTN_st_go', 'rCaudate_st_go')][cmpl,])

  fm <- as.formula(paste('rSTN_st_go', paste(feats_left, collapse = ' + '), sep = ' ~ ') )
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









# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
#task_imp_cc <- function(data, feats_left, feats_right){
# What items do we want pairwise correlations for?
#  cor_items <- c('tci_gen_imp','tci_fin_imp','surps_imp','log10.k.','bmi',
#                 'sj1a', 'sj1b','sj1c', 'espad_6_life_nic', 'espad_8a_alc_life')

#  task_names <- c(data$names$sst, data$names$mid, data$names$genes,
#                  data$names$CANTAB)
#  cca_wrapper( complete(data$imp,1)[task_names][data$train_inds,], complete(data$imp,1)[cor_items][data$train_inds,])
#}