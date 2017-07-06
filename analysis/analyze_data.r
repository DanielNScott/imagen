compute_correlations <- function(data) {

  # What items do we want pairwise correlations for?
  cor_items <- c('tci_gen_imp','tci_fin_imp','surps_imp','log10.k.','bmi',
                  'sj1a', 'sj1b','sj1c', 'espad_6_life_nic', 'espad_8a_alc_life', 'binge')

  # Get correlations
  correlations <- 0
  for (ind in 1:20) {
    correlations <- correlations + 1/20 * cor( complete(data$imp,ind)[cor_items] )
  }

  # This fn will put correlations on the panels of the matrix heatmap
  myPanel <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    panel.text(x, y, round(z,1))
  }

  # Plotting setup
  scale <- list(x = list(rot = 90))
  #options(repr.plot.width = 9, repr.plot.height = 7.5)

  # Plot...
  library(lattice)
  levelplot(round(correlations,4), scales = scale, main = 'Pairwise Correlations', xlab = '', ylab = '',
            at = seq(-1, 1, length.out = 100), panel = myPanel)
}


task_imp_cc <- function(data){
  # What items do we want pairwise correlations for?
  cor_items <- c('tci_gen_imp','tci_fin_imp','surps_imp','log10.k.','bmi',
                 'sj1a', 'sj1b','sj1c', 'espad_6_life_nic', 'espad_8a_alc_life')

  task_names <- c(data$names$sst, data$names$mid, data$names$genes,
                  data$names$CANTAB)
  cca_wrapper( complete(data$imp,1)[task_names][data$train_inds,], complete(data$imp,1)[cor_items][data$train_inds,])
}


function() {
  library(CCA)
  library(mice)
  source('cca_wrapper.r')
  #with(imp, cancor())
  #with(imp, cancor(imp['IQ_PR_14','IQ_VC_14'],betas_df))

  tasks <- c('IQ_PR_14', 'IQ_VC_14', 'cgt_delay_avers_14', 'cgt_delib_14', 'cgt_quality_14', 'cgt_prop_bet_14',
             'cgt_risk_adjust_14',  'log10.k._14', 'mu_go_14', 'mu_stop_14', 'sigma_go_14',
             'sigma_stop_14')

  #tasks <- c('IQ_PR_14', 'IQ_VC_14', 'mu_go_14', 'mu_stop_14', 'sigma_go_14', 'sigma_stop_14')

  #betas <- c('stn_go','stn_stop')


  # Perform All Data PCA:
  pca_res  <- prcomp(betas_df, center = TRUE, scale. = TRUE)
  pca_pvar <- pca_res$sdev^2/sum(pca_res$sdev^2)

  ### PLOT THE VARIATION CAPTURED BY PCAs:
  #par(mfrow=c(1,3), mai = c(0.1, 0.5, 0.1, 0.3))
  options(repr.plot.width=9, repr.plot.height=3)

  ylab_top <- "Proportion of Variance [%]"

  # Task PCA Plot:
  barplot(pca_pvar,  col=c('red'), las=2, cex.axis=0.7, xlab=c('Principle Components'), ylab=ylab_top)

  ## PLOT THE LOADINGS FOR EACH PCA:
  #par(mfrow=c(1,3), mai = c(0.1, 0.5, 0.1, 0.3))
  options(repr.plot.width=5, repr.plot.height=5)

  barplot(t(pca_res$rotation[,1:3]), main="", horiz=TRUE, beside=TRUE, col=rainbow(3), yaxt='n')
  axis(2, labels=colnames(betas_df), at=seq(3,6*5,by=5), las=2)
  legend("topright", legend=c('PC1', 'PC2', 'PC3'), fill = rainbow(3))


  # project new data onto the PCA space
  new_df <- scale(betas_df, pca_res$center, pca_res$scale) %*% pca_res$rotation

  dim(complete(imp,)[tasks])
  for (i in 4:5) {
    cc <- cca_wrapper(complete(imp,1)[tasks], new_df[,1:3])
  }
}