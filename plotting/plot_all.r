# ------------------------------------------------------------------------------ #
# Function to call all the plotting
# ------------------------------------------------------------------------------ #
analysis <- function(data) {
  source('./plot_wrappers.r')
  ################################################################################
  #                                                                              #
  #                               MODEL PARAMETERS                               #
  #                                                                              #
  ################################################################################
  # Subset Histograms
  ggthemr('flat')
  d <- melt( data$raw[c(data$names$sst, 'set')], id.vars = 'set')
  ggplot(d, aes(x = value, fill = set)) +
    facet_wrap(~variable, scales = 'free') +
    geom_histogram(alpha = 0.4, position = 'identity') +
    ggtitle('SST Parameter Histograms by Subset') +
    xlab('time [ms]')

  d <- melt( data$raw[c(data$names$mid, 'set')], id.vars = 'set')
  ggplot(d, aes(x = value, fill = set)) +
    facet_wrap(~variable, scales = 'free') +
    geom_histogram(alpha = 0.4, position = 'identity') +
    ggtitle('MID Parameter Histograms by Subset') +
    xlab('time [ms]')

  # SPLOMS
  title <- 'Group Distributions of (Some) Stop Signal Task Parameters'
  feats_sst <- c('GoRT', 'SSRT', 'GoRTSd', 'shift_go')
  feats_rename <- c('SSRT', 'Sequential Change in GoRT', 'Sequential Change in StopRT')
  ggpairs_wrap(data$raw[feats_sst], title)#, 'time [ms]', 'density [ ]', xticks = feats_rename )

  title <- 'Group Distributions of (Some) Stop Signal Task Parameters'
  ggpairs_wrap(data$raw[data$names$mid], title)#, 'time [ms]', 'density [ ]', xticks = feats_rename )

  ggthemr('fresh')
  ################################################################################
  #                                                                              #
  #                                LINEAR MODELS                                 #
  #                                                                              #
  ################################################################################
  # Function for extracting coeffiient values etc, to be used in each plot below.
  task_lm <- function(data, inds, dep, names) {
    task_vars <- as.matrix(data$scores[inds, names])
    colnames(task_vars) <- names
    dep_var <- as.vector(data$scores[dep][inds,])
    model   <- lm(dep_var ~ task_vars)
    smry    <- summary(model)
    rownames(smry$coefficients) <- c('Intercept', colnames(task_vars))
    coef    <- smry$coefficients
    print(smry)
    return(coef)
  }
  # --------------------------- PRIMARY COMPARISONS -----------------------------#
  ###### MID #####
  roi    <- 'MID_VS'
  names  <- c('mid_dur', 'mid_loc', 'mid_rew', 'mid_high_rew', 'mid_int')
  #names  <- c('mid_RT_ar', 'mid_RT_hr', 'mid_RT_var')
  xlabel <- 'Model Parameter'
  ylabel <- 'Coefficient Value'

  coef   <- task_lm(data, data$train_inds, roi, names)
  title  <- 'Train Subset'
  plt1   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  coef   <- task_lm(data, data$test_inds, roi, names)
  title  <- 'Test Subset'
  plt2   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  coef   <- task_lm(data, c(data$train_inds, data$test_inds), roi, names)
  title  <- paste(roi, 'Contrast ~ MID Model:', 'All')
  plt3   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  print(ggarrange(plt3, ggarrange(plt1, plt2, ncol = 2, labels = c('B', 'C')), nrow = 2, labels = 'A'))

  ###### SST #####
  #ctrst_rois <- c('rIFG_st_go', 'rPreSMA_st_go')#, 'rSTN_st_go', 'rGPe_st_go', 'rGPi_st_go')
  ctrst_rois <- c('rGPe_st_go')
  reliable_sst <- c('GoRT', 'SSRT', 'GoRTSd', 'shift_go')
  xlabel <- 'Model Parameter'
  ylabel <- 'Coefficient Value'
  title_core <- 'Contrast ~ SST Model:'

  for (roi in ctrst_rois) {
    print(paste('--------------------- ROI:', roi, '----------------------------'))
    coef   <- task_lm(data, data$train_inds, roi, reliable_sst)
    title  <- 'Train Subset'
    plt1   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

    coef   <- task_lm(data, data$test_inds, roi, reliable_sst)
    title  <- 'Test Subset'
    plt2   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

    coef   <- task_lm(data, c(data$train_inds, data$test_inds), roi, reliable_sst)
    title  <- paste(roi, title_core, 'All')
    plt3   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

    print(ggarrange(plt3, ggarrange(plt1, plt2, ncol = 2, labels = c('B', 'C')), nrow = 2, labels = 'A'))
  }

  # --------------------------- Secondary Comparisons ---------------------------#
  ###### MID_VS by SST Model #####
  roi    <- 'MID_VS'
  names  <- c('GoRT', 'SSRT', 'GoRTSd', 'shift_stop')
  xlabel <- 'Model Parameter'
  ylabel <- 'Coefficient Value'

  coef   <- task_lm(data, data$train_inds, roi, names)
  title  <- 'Train Subset'
  plt1   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  coef   <- task_lm(data, data$test_inds, roi, names)
  title  <- 'Test Subset'
  plt2   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  coef   <- task_lm(data, c(data$train_inds, data$test_inds), roi, names)
  title  <- paste(roi, 'Contrast ~ SST Model:', 'All')
  plt3   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  print(ggarrange(plt3, ggarrange(plt1, plt2, ncol = 2, labels = c('B', 'C')), nrow = 2, labels = 'A'))

  ###### MID_VS by SST ROIS #####
  roi    <- 'MID_VS'
  names  <- c('rIFG_st_go', 'rPreSMA_st_go', 'rCaudate_st_go', 'rSTN_st_go', 'rGPe_st_go', 'rGPi_st_go')
  xlabel <- 'SST ROI Contrast'
  ylabel <- 'Coefficient Value'

  coef   <- task_lm(data, data$train_inds, roi, names)
  title  <- 'Train Subset'
  plt1   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  coef   <- task_lm(data, data$test_inds, roi, names)
  title  <- 'Test Subset'
  plt2   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  coef   <- task_lm(data, c(data$train_inds, data$test_inds), roi, names)
  title  <- paste(roi, 'Contrast ~ SST Model:', 'All')
  plt3   <- bar_wrap(coef[,1:2], title = title, xlabel = xlabel, ylabel = ylabel)

  print(ggarrange(plt3, ggarrange(plt1, plt2, ncol = 2, labels = c('B', 'C')), nrow = 2, labels = 'A'))


  ################################################################################
  #                                                                              #
  #                          CANONICAL CORRELATIONS                              #
  #                                                                              #
  ################################################################################

  # ---------------------------- Primary Comparisons ----------------------------#
  # SST Params ~ Stop Network ---- Train Subset
  feats_task <- c('GoRT', 'GoRTSd', 'SSRT', 'shift_go')
  feats_fmri <- c('rIFG_st_go', 'rPreSMA_st_go', 'rSTN_st_go', 'rCaudate_st_go', 'rGPi_st_go')
  scores <- data$scores[c(feats_task, feats_fmri)]
  cmpl   <- complete.cases(scores[data$train_inds,])
  cca_wrapper(scores[feats_task][data$train_inds,][cmpl,],
              scores[feats_fmri][data$train_inds,][cmpl,],
              'Correlations', 'SST Model Parameters', 'fMRI ROI Activity', subset = 1)

  # SST Params ~ Stop Network ---- Test Subset
  feats_task <- c('GoRT', 'GoRTSd', 'SSRT', 'shift_go')
  feats_fmri <- c('rIFG_st_go', 'rPreSMA_st_go', 'rSTN_st_go', 'rCaudate_st_go', 'rGPi_st_go')
  scores <- data$scores[c(feats_task, feats_fmri)]
  cmpl   <- complete.cases(scores[data$test_inds,])
  cca_wrapper(scores[feats_task][data$test_inds,][cmpl,],
              scores[feats_fmri][data$test_inds,][cmpl,],
              'Correlations', 'SST Model Parameters', 'fMRI ROI Activity')

  # SST Params ~ Stop Network ---- All
  feats_task <- c('GoRT', 'GoRTSd', 'SSRT', 'shift_go')
  feats_fmri <- c('rIFG_st_go', 'rPreSMA_st_go', 'rSTN_st_go', 'rCaudate_st_go', 'rGPi_st_go')
  scores <- data$scores[c(feats_task, feats_fmri)]
  inds   <- c(data$train_inds, data$test_inds)
  cmpl   <- complete.cases(scores[inds,])
  cca_wrapper(scores[feats_task][inds,][cmpl,],
              scores[feats_fmri][inds,][cmpl,],
              'Correlations', 'SST Model Parameters', 'fMRI ROI Activity')

  ####

  # SST Params ~ Stop Network ---- Training Subset
  feats_task <- c('GoRT', 'GoRTSd', 'SSRT', 'shift_go')
  feats_fmri <- c('rIFG_st_go', 'rPreSMA_st_go', 'rSTN_st_go', 'rCaudate_st_go', 'rGPi_st_go', 'MID_VS')
  scores <- data$scores[c(feats_task, feats_fmri)]
  cmpl   <- complete.cases(scores[data$train_inds,])
  cca_wrapper(scores[feats_task][data$train_inds,][cmpl,],
              scores[feats_fmri][data$train_inds,][cmpl,],
              '', 'Stop Signal Model Parameters', 'Region of Interest fMRI Activity')

  # MID + SST ~ MID_VS + STOP_NETWORK
  feats_task <- c('GoRT', 'GoRTSd', 'SSRT', 'shift_go', data$names$mid)
  feats_fmri <- c('rIFG_st_go', 'rPreSMA_st_go', 'rSTN_st_go', 'rCaudate_st_go', 'rGPi_st_go', 'MID_VS')
  scores <- data$scores[c(feats_task, feats_fmri)]
  cmpl   <- complete.cases(scores[data$train_inds,])
  cca_wrapper(scores[feats_task][data$train_inds,][cmpl,],
              scores[feats_fmri][data$train_inds,][cmpl,],
              '', 'Stop Signal Model Parameters', 'Region of Interest fMRI Activity', subset = 1:2)


  feats_left <- c('mid_rew', 'mid_high_rew', 'GoRT', 'SSRT', 'shift_go', data$names$CANTAB, setdiff(data$names$Age, 'bmi'))
  std_fmri_feats <- as.vector(outer(rois, sfx, function(x,y) {paste(x,y, sep = '')}))
  cmpl <- complete.cases(data$scores[c(feats_left, 'rSTN_st_go', 'rCaudate_st_go')])
  cca_wrapper(data$scores[feats_left][cmpl,], data$scores[c('rSTN_st_go', 'rCaudate_st_go')][cmpl,])

  # Task - FMRI CCA
  feats_left <- c('mid_rew', 'mid_high_rew', 'GoRT', 'SSRT', 'shift_go', data$names$CANTAB, setdiff(data$names$Age, 'bmi'))
  std_fmri_feats <- as.vector(outer(rois, sfx, function(x,y) {paste(x,y, sep = '')}))
  cmpl <- complete.cases(data$scores[c(feats_left, std_fmri_feats)])
  cca_wrapper(data$scores[feats_left][cmpl,], data$scores[std_fmri_feats][cmpl,])


  # CCA of IMP and FMRI
  feats_imp <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  cmpl <- complete.cases(data$scores[c(feats_imp, feats_core_fmri)])
  cca_wrapper(data$scores[feats_imp][cmpl,], data$scores[feats_core_fmri][cmpl,],
              '', '', '')

  # --------------------------- Secondary Comparisons ---------------------------#
  # SST -- IMP
  feats_a <- data$names$sst
  feats_b <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount', 'adhd_teacher', 'adhd_parent', 'adhd_child', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  cca_wrapper(data$scores[feats_a][cmpl,], data$scores[feats_b][cmpl,], '','','')

  # CANTAB -- IMP
  feats_a <- data$names$CANTAB
  feats_b <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount', 'adhd_teacher', 'adhd_parent', 'adhd_child', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  cca_wrapper(data$scores[feats_a][cmpl,], data$scores[feats_b][cmpl,],'','','')

  # ROB -- IMP
  feats_a <- c("SS_rVS", "SS_SN_STN", "SS_pSMA_PCG", "SS_lIFG_BA10")
  feats_b <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount', 'adhd_teacher', 'adhd_parent', 'adhd_child', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  cca_wrapper(data$scores[feats_a][cmpl,], data$scores[feats_b][cmpl,])

  ################################################################################
  #                                                                              #
  #                                   SPLOMS                                     #
  #                                                                              #
  ################################################################################
  # SST 1
  title <- 'Stop Signal Task Parameters'
  feats <- c("tau_stop", "tau_go", "mu_go", "mu_stop")
  ggpairs_wrap(data$scores[feats], title)

  title <- 'Group Distributions of (Some) MID Task Parameters'
  feats_mid <- c('mid_rew', 'mid_high_rew', 'mid_rew_var', 'mid_high_rew_var')
  density_wrapper(data$raw, feats_mid, title)

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

  # Std fMRI
  #rois <- c('rPreSMA', 'rIFG', 'rCaudate', 'rSTN', 'rGPe', 'rGPi', 'rThalamus')
  rois <- c('rIFG', 'rPreSMA', 'rCaudate', 'rSTN', 'rGPe', 'rGPi')
  titles <- c('SST Go Weights', 'SST Stop Weights', 'SST Stop Respond Weights',
              'SST Stop > Go Contrasts', 'SST Stop > Stop Respond Contrasts')
  sfx <- c('_st_go') #,'_go','_st')
  for (i in 1:length(sfx)) {
    std_fmri_feats <- paste(rois, sfx[i], sep = '')
    print(ggpairs_wrap(data$raw[std_fmri_feats], title = 'fMRI GLM Weights'))
  }

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

  ################################################################################
  #                                                                              #
  #                         PRINCIPAL COMPONENTS ANALYSES                        #
  #                                                                              #
  ################################################################################
  feats_a <- c("SS_rVS", "SS_STN", "SS_pSMA", "SS_BA10")
  feats_b <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount', 'adhd_teacher', 'adhd_parent', 'adhd_child', 'bmi')
  title <- 'PCA: fMRI & Impulsivity Measures'
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  pca_wrapper(data$scores[c(feats_a, feats_b)][cmpl,], title)

  feats_a <- c("SS_rVS", "SS_SN_STN", "SS_pSMA_PCG", "SS_lIFG_BA10")
  feats_b <- c('nic_use', 'alc_use', 'alc_regret', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  title <- 'PCA: fMRI & Substance Use'
  pca_wrapper(data$scores[c(feats_a, feats_b)][cmpl,], title)

  #
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


  ################################################################################
  #                                                                              #
  #                                    MISC.                                     #
  #                                                                              #
  ################################################################################

  #---------------------------- Plots in FYP pres. ------------------------------#

  # SST Param Density
  title <- 'Group Distributions of (Some) Stop Signal Task Parameters'
  feats_sst <- c("SSRT", "shift_go", "shift_stop")
  feats_rename <- c('SSRT', 'Sequential Change in GoRT', 'Sequential Change in StopRT')
  density_wrapper(data$raw, feats_sst, title, 'time [ms]', 'density [ ]', xticks = feats_rename )

  # PCA of ROI fMRIs
  title <- 'PCA of fMRI Activity for Stop-Go Contrast'
  feats_core_fmri_plain <- c('rIFG', 'rPreSMA', 'rCaudate', 'rSTN', 'rGPe', 'rGPi')
  feats_core_fmri <- paste(feats_core_fmri_plain, '_st_go', sep = '')
  pca.res <- pca_wrapper(data$scores[feats_core_fmri][cmpl,], title, subset = 1:3, xticks = feats_core_fmri_plain)
  proj <- pca.res$proj
  mask <- pca.res$mask

  # Impulsivity ~ PCA projections
  feats_imp <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  cmpl2 <- complete.cases(cbind(proj, data$scores[feats_imp][cmpl,]))
  cca_wrapper( proj[cmpl2,], data$scores[feats_imp][cmpl,][cmpl2,], '', '', '')

  # CANTAB -- ESPAD+
  feats_a <- data$names$CANTAB[5:10]
  feats_b <- c('nic_use', 'alc_use', 'alc_regret', 'bmi')
  cmpl <- complete.cases(data$scores[c(feats_a, feats_b)])
  cca_wrapper(data$scores[feats_a][cmpl,], data$scores[feats_b][cmpl,],'','','', subset = 1)

  # CCA of MID and SST on std fMRI & VS
  feats_sst_less <- c('SSRT', 'shift_go')
  feats_core_fmri_plain <- c('rIFG', 'rPreSMA', 'rSTN', 'rCaudate', 'rGPi')
  feats_core_fmri <- paste(feats_core_fmri_plain, '_st_go', sep = '')
  cmpl <- complete.cases(data$scores[c(feats_sst_less, c(feats_core_fmri, 'MID_VS'))])
  cca_wrapper(data$scores[feats_sst_less][cmpl,], data$scores[c(feats_core_fmri, 'MID_VS')][cmpl,], '',
              'Stop Signal Model Parameters', 'Region of Interest fMRI Activity', subset = 1)

  # CCA of MID
  feats_imp <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  feats_core_fmri_plain <- c('rIFG', 'rPreSMA', 'rSTN', 'rCaudate', 'rGPi')
  feats_core_fmri <- paste(feats_core_fmri_plain, '_st_go', sep = '')
  cmpl <- complete.cases(data$scores[c(feats_imp, c(feats_core_fmri, 'MID_VS'))])
  cca_wrapper(data$scores[feats_imp][cmpl,], data$scores[c(feats_core_fmri, 'MID_VS')][cmpl,], '',
              'Stop Signal Model Parameters', 'Region of Interest fMRI Activity', subset = 1)

  # Impulsivity Metrics
  feats_imp <- c('tci_gen_imp', 'tci_fin_imp', 'surps_imp', 'discount')
  correlations <- cor(data$raw[feats_imp], use = 'pairwise.complete.obs')
  correlations <- matrix(unlist(correlations), ncol = 4, byrow = TRUE)
  colnames(correlations) <- feats_imp
  yticklabel <- c('General (TCI)', 'Financial (TCI)', 'General (SURPS)', 'Discounting')
  level_wrapper(melt(correlations), 'Impulsivity Metrics', ylabel = 'Assessment')

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

  # Impulsivity as a function of SS_STN and MID_VS
  color_point_wrap(data$scores[c('SS_STN','MID_VS','tci_gen_imp')], 'Impulsivity by VS & STN Activity')
  color_point_wrap(data$scores[c('SS_STN','MID_VS',  'surps_imp')], 'Impulsivity by VS & STN Activity')
  color_point_wrap(data$scores[c('SS_STN','MID_VS',   'discount')], 'Impulsivity by VS & STN Activity')
  color_point_wrap(data$scores[c('SS_STN','MID_VS','tci_fin_imp')], 'Impulsivity by VS & STN Activity')
  color_point_wrap(data$scores[c('rSTN_st_go','MID_VS','tci_gen_imp')], 'Impulsivity by VS & STN Activity')
  color_point_wrap(data$scores['tci_fin_imp'], 'Impulsivity by Network PCAs')

  # Covariance matrices
  avg <- 0; n <- 0;
  for (subj_num in 1:10) {
    roi_nums <- c(2,3,5,6)

    subj_cov <- rebuild_covar(fmri_cov, subj_num, roi_nums)
    melt_cov <- melt(cov2cor(subj_cov))

    title <- paste('Correlations in ROI Activity:', toString(subj_id))
    plot  <- level_wrapper(melt_cov, title, ylabel = NULL)
    print(plot)

    # Also, compute average as we're doing this,
    # masking out the dubious ones.
    if (sum(cov2cor(subj_cov) < 8)) {avg <- avg + subj_cov; n <- n + 1; print(subj_cov)}
  }
  avg <- avg/n

  #------------------------------------------------------------------------------#
  #          Newman-Girvan community detection on covariance matrices            #
  #------------------------------------------------------------------------------#
  # source(paste(dir_analy, '/community_detection.r'   , sep = ''))
  #dist <- matrix(nrow = 190, ncol = 190)
  #for (i in 1:190) {
  #  for (j in i:190) {
  #    dist[i,j] <- cov_dist(ag$fitST[[i]]$Shat, ag$fitST[[j]]$Shat)
  #  }
  #}
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