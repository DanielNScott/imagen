library(ggplot2)
library(reshape2)
library(GGally)
library(ggthemr)
ggthemr('fresh')

# ------------------------------------------------------------------------------ #
# Function to call all the plotting
# ------------------------------------------------------------------------------ #
analysis <- function(data) {


  source('./plot_wrappers.r')

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



  # ------------------------------------
  ### Densities
  # ------------------------------------

  # SST 1
  title <- 'Stop Signal Task Parameters'
  feats <- c("tau_stop", "tau_go", "mu_go", "mu_stop")
  ggpairs_wrap(data$scores[feats], title)

  # SST 2
  title <- 'Stop Signal Task Parameters'
  feats_sst <- c("shift_go", "shift_stop", "GoRT", "SSRT")
  ggpairs_wrap(data$scores[feats], title)

  # CANTAB Params
  cgt_names <- c('delay avers', 'delib time', 'quality', 'prop. bet', 'risk adj.', 'risk take')
  feats_b <- data$names$CANTAB[5:10]
  ggpairs_wrap(data$raw[feats_b], title = 'CANTAB Statistics', cgt_names)

  # ROI GLM
  title <- 'ROI GLM Coefficients'
  rob_feats <- c('SS_STN', 'SS_pSMA','SS_rVS')
  ggpairs_wrap(data$raw[rob_feats], title)

  # Impulsivity Metrics
  feats_imp <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  ggpairs_wrap(data$scores[feats_a], 'Impulsivity Metrics')

  # ADHD
  feats_b <- c('adhd_teacher', 'adhd_parent', 'adhd_child')
  ggpairs_wrap(data$scores[feats_b], 'Hyperactivity Metrics')

  # ------------------------------------
  ### PAIR-WISE
  # ------------------------------------

  # Imp vs ADHD vs SU
  feats_a <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  feats_b <- c('adhd_teacher', 'adhd_parent', 'adhd_child', 'nic_use', 'alc_use', 'alc_regret', 'bmi')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  #------

  # Imp vs ADHD
  feats_a <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  feats_b <- c('adhd_teacher', 'adhd_parent', 'adhd_child')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)], 'Impulsivity Metrics')

  # Imp vs SST
  feats_a <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  feats_b <- c('GoRT', 'SSRT', 'shift_go')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # Imp vs MID
  feats_a <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  feats_b <- c('mid_rew', 'mid_high_rew', 'mid_rew_var', 'mid_high_rew_var')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # IMP vs CANTAB
  feats_a <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  feats_b <- data$names$CANTAB[5:11]
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # IMP vs AGN
  feats_a <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  feats_b <- data$names$CANTAB[1:4]
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  #------

  # ADHD vs SST
  feats_a <- c('adhd_teacher', 'adhd_parent', 'adhd_child')
  feats_b <- c('GoRT', 'SSRT', 'shift_go')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # ADHD vs MID
  feats_a <- c('adhd_teacher', 'adhd_parent', 'adhd_child')
  feats_b <- c('mid_rew', 'mid_high_rew', 'mid_rew_var', 'mid_high_rew_var')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # ADHD vs CANTAB
  feats_a <- c('adhd_teacher', 'adhd_parent', 'adhd_child')
  feats_b <- data$names$CANTAB[5:11]
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # ADHD vs AGN
  feats_a <- c('adhd_teacher', 'adhd_parent', 'adhd_child')
  feats_b <- data$names$CANTAB[1:4]
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  #------

  # SU vs SST
  feats_a <- c('nic_use', 'alc_use', 'alc_regret', 'bmi')
  feats_b <- c('GoRT', 'SSRT', 'shift_go')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # SU vs MID
  feats_a <- c('nic_use', 'alc_use', 'alc_regret', 'bmi')
  feats_b <- c('mid_rew', 'mid_high_rew', 'mid_rew_var', 'mid_high_rew_var')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # SU vs CANTAB
  feats_a <- c('nic_use', 'alc_use', 'alc_regret', 'bmi')
  feats_b <- data$names$CANTAB[5:11]
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # SU vs AGN
  feats_a <- c('nic_use', 'alc_use', 'alc_regret', 'bmi')
  feats_b <- data$names$CANTAB[1:4]
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  #------

  # CANTAB VS AGN
  feats_a <- data$names$CANTAB[1:4]
  feats_b <- data$names$CANTAB[6:11]
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # CANTAB VS SST
  feats_a <- data$names$CANTAB[6:11]
  feats_b <- c('GoRT', 'SSRT', 'shift_go')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  # CANTAB VS MID
  feats_a <- data$names$CANTAB[6:11]
  feats_b <- c('GoRT', 'SSRT', 'shift_go')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  #------

  # SST vs MID
  feats_a <- c('GoRT', 'SSRT', 'shift_go')
  feats_b <- c('mid_rew', 'mid_high_rew', 'mid_rew_var', 'mid_high_rew_var')
  ggpairs_wrap(data$scores[c(feats_a, feats_b)])

  #------

  # CANTAB
  feats_a <- data$names$CANTAB[5:11]
  ggpairs_wrap(data$scores[ feats_a])

  # Std fMRI
  rois <- c('rPreSMA', 'rIFG', 'rCaudate', 'rSTN', 'rGPe', 'rGPi', 'rThalamus')
  titles <- c('SST Go Weights', 'SST Stop Weights', 'SST Stop Respond Weights',
              'SST Stop > Go Contrasts', 'SST Stop > Stop Respond Contrasts')
  sfx <- c('_st_go','_go','_st')
  for (i in 1:length(sfx)) {
    std_fmri_feats <- paste(rois, sfx[i], sep = '')
    print(ggpairs_wrap(data$raw[std_fmri_feats], title = 'fMRI GLM Weights'))
  }

  # --------------------------------------------------------
  ###                          PCAs
  # --------------------------------------------------------
  feats_a <- c("SS_rVS", "SS_STN", "SS_pSMA", "SS_BA10")
  feats_b <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount', 'adhd_teacher', 'adhd_parent', 'adhd_child', 'bmi')

  title <- 'PCA: fMRI & Impulsivity Measures'
  pca_wrapper(data$scores[c(feats_a, feats_b)][cmpl,], title)

  feats_a <- c("SS_rVS", "SS_SN_STN", "SS_pSMA_PCG", "SS_lIFG_BA10")
  feats_b <- c('nic_use', 'alc_use', 'alc_regret', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  title <- 'PCA: fMRI & Substance Use'
  pca_wrapper(data$scores[c(feats_a, feats_b)][cmpl,], title)

  ###
  cmpl <- complete.cases(data$scores[c(feats_sst, rob_feats)])
  title <- 'PCA: SST Parameters & fMRI'
  pca_wrapper(data$scores[c(feats_sst, rob_feats)][cmpl,], title)
  title <- 'CCA: SST Params & fMRI'
  cca_wrapper(data$scores[feats_sst][cmpl,],data$scores[rob_feats][cmpl,], title, title, title)

  cmpl <- complete.cases(data$scores[c(rob_feats, feats_imp)])
  title <- 'PCA: fMRI & Impulsivity'
  pca_wrapper(data$scores[c(rob_feats, feats_imp)][cmpl,], title)
  title <- 'CCA: fMRI & Impulsivity'
  cca_wrapper(data$scores[rob_feats][cmpl,],data$scores[feats_imp][cmpl,], title)
  cca_wrapper(data$scores[rob_feats][cmpl,],data$scores[feats_imp][cmpl,], title, title, title)

  # --------------------------------------------------------
  ###                          CCAs
  # --------------------------------------------------------
  feats_a <- data$names$sst
  feats_b <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount', 'adhd_teacher', 'adhd_parent', 'adhd_child', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  cca_wrapper(data$scores[feats_a][cmpl,], data$scores[feats_b][cmpl,])

  feats_a <- c("SS_rVS", "SS_SN_STN", "SS_pSMA_PCG", "SS_lIFG_BA10")
  feats_b <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount', 'adhd_teacher', 'adhd_parent', 'adhd_child', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  cca_wrapper(data$scores[feats_a][cmpl,], data$scores[feats_b][cmpl,])

  # --------------------------------------------------------
  ###                          t-SNE
  # --------------------------------------------------------
  feats_a <- c("SS_rVS", "SS_SN_STN", "SS_pSMA_PCG", "SS_lIFG_BA10")
  feats_b <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount', 'adhd_teacher', 'adhd_parent', 'adhd_child', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])

  title <- 'Nonlinear Dim. Reduction by TCI Impulsivity'
  tsne_wrapper(data$scores[c(feats_a, feats_b)][cmpl,], 'tci_gen_imp', title)

  title <- 'Nonlinear Dim. Reduction by SURPS Impulsivity'
  tsne_wrapper(data$scores[c(feats_a, feats_b)][cmpl,], 'surps_imp', title)

  # --------------------------------------------------------
  ###                          CCAs
  # --------------------------------------------------------

  # Task - FMRI CCA
  feats_left <- c('mid_rew', 'mid_high_rew', 'GoRT', 'SSRT', 'shift_go', data$names$CANTAB, setdiff(data$names$Age, 'bmi'))
  std_fmri_feats <- as.vector(outer(rois, sfx, function(x,y) {paste(x,y, sep = '')}))
  cmpl <- complete.cases(data$scores[c(feats_left, 'rSTN_st_go', 'rCaudate_st_go')])
  cca_wrapper(data$scores[feats_left][cmpl,], data$scores[c('rSTN_st_go', 'rCaudate_st_go')][cmpl,])

  # Task - FMRI CCA
  feats_left <- c('mid_rew', 'mid_high_rew', 'GoRT', 'SSRT', 'shift_go', data$names$CANTAB, setdiff(data$names$Age, 'bmi'))
  std_fmri_feats <- as.vector(outer(rois, sfx, function(x,y) {paste(x,y, sep = '')}))
  cmpl <- complete.cases(data$scores[c(feats_left, std_fmri_feats)])
  cca_wrapper(data$scores[feats_left][cmpl,], data$scores[std_fmri_feats][cmpl,])

  #  act. by features as linear model
  fm <- as.formula(paste('rSTN_st_go', paste(feats_left, collapse = ' + '), sep = ' ~ ') )

  # Scatterplots of STN connectivity
  plotdata <- data.frame(x = imputed$ifg..stn, y = imputed$presma..stn, z = imputed$stn..gpi)
  ggplot(plotdata, aes(x, y, colour = z)) + scale_color_viridis(option = 'plasma') +
    geom_point() + theme_bw() +
    labs(x = 'IFG to STN', y = 'preSMA to STN', colour = 'STN to GPi', title = 'STN Connectivity') +
    theme(legend.position = c(0.9, 0.75)) +
    geom_abline(slope = -0.399, intercept = -0.385, lty = 2) +
    annotate("text", x = 0.1, y = -0.6, label = "Slope: -0.4 \n R: 0.32")

  # Repeat of prev plot but with thresholding at midpoint of STN - GPi connectivty
  z <- imputed$stn..gpi
  npts <- 2
  colors <- viridis(npts)
  zcolor <- colors[(z - min(z))/diff(range(z))*npts + 1]
  plotdata <- data.frame(x = imputed$ifg..stn, y = imputed$presma..stn)
  ggplot(plotdata, aes(x, y, colour = zcolor)) + geom_point() + theme_bw() +
    labs(x = 'IFG to STN', y = 'preSMA to STN', colour = 'STN to GPi', title = 'STN Connectivity') +     theme(legend.position = c(0.85, 0.8)) +
    scale_color_identity('STN to GPi', labels = c('-0.64 < val < 0.03', ' 0.03 < val < 0.58'), breaks = colors, guide = "legend")


  ##########################3
  ## Completely different plots
  ######

  # Impulsivity as a function of SS_STN and MID_VS
  color_point_wrap(data$scores[c('SS_STN','MID_VS','tci_gen_imp')], 'Impulsivity by VS & STN Activity')
  color_point_wrap(data$scores[c('SS_STN','MID_VS',  'surps_imp')], 'Impulsivity by VS & STN Activity')
  color_point_wrap(data$scores[c('SS_STN','MID_VS',   'discount')], 'Impulsivity by VS & STN Activity')
  color_point_wrap(data$scores[c('SS_STN','MID_VS','tci_fin_imp')], 'Impulsivity by VS & STN Activity')

  color_point_wrap(data$scores['tci_fin_imp'], 'Impulsivity by Network PCAs')

}
# ------------------------------------------------------------------------------ #

#{
#  # ... ancestral graph fmri weights
#  rois <- c('presma.ifg', 'presma..stn', 'ifg..stn', 'presma..caud', 'ifg..caud',
#            'caud..gpe', 'gpe..gpi', 'stn..gpi', 'gpi..thalamus')
#  titles <- c('SST Stop Conn.', 'SST Stop Respond Conn', 'SST Stop, Stop Respond Contrasts')
#  sfx <- c('_st', '_sr', '_st_sr')
#  for (i in 1:length(sfx)) {
#    feats <- paste(rois, sfx[i], sep = '')
#    print(density_plot(data, feats, titles[i]))
#  }
#}