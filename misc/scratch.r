options(repr.plot.width=8, repr.plot.height=3)

nplts <- 3
for (subj_num in 1:nplts) {
    dmelt <- melt(cbind(lm_list[[subj_num]]$conv_regs[1:40,c(1,3,8)], stn_series[1:40,subj_num]/12.5), id.vars = 'TR', variable.name = 'Condition')
    plt <- ggplot(dmelt, aes(TR, value)) + geom_line(aes(colour=Condition))
    print(plt)
}

##### DESIGN MATRIX w/ DRIFTS
#dmatdf <- data.frame(design_mat)
#colnames(dmatdf) <- c(labels[setdiff(1:7, bad_cond)], 'Intercept', 'Linear Drift', 'Sq. Drift')
#dmatdf['TR'] <- 1:444

#dmatdf.m <- melt(dmatdf, id.vars = 'TR')
#dmatdf.m <- ddply(dmatdf.m, .(variable), transform, rescale = scale(value))

#options(repr.plot.width=6, repr.plot.height=5)

#p <- ggplot(dmatdf.m, aes(variable, TR)) +
#      geom_tile(aes(fill = value)) +
#      scale_fill_gradient(low = "white", high = "steelblue") +
#      labs(x = "") + scale_x_discrete(expand = c(0, 0)) + #scale_y_discrete(expand = c(0, 0)) +
#      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#              panel.background = element_blank(),
#              legend.position = "none",
#              axis.ticks = element_blank(),
#              axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))

#print(p)

options(repr.plot.width=6, repr.plot.height=3)

dmelt <- melt(betas[c('STOP_SUCCESS', 'STOP_FAILURE')], id.vars = 'STOP_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(STOP_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_y_continuous('Stop Failure Beta') + scale_x_continuous('Stop Success Beta') + theme(aspect.ratio = 1)

dmelt <- melt(betas[c('GO_SUCCESS', 'STOP_SUCCESS')], id.vars = 'GO_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(GO_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_x_continuous('GO Success Beta') + scale_y_continuous('Stop Success Beta') + theme(aspect.ratio = 1)

dmelt <- melt(betas[c('GO_SUCCESS', 'STOP_FAILURE')], id.vars = 'GO_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(GO_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_x_continuous('Go Success Beta') + scale_y_continuous('Stop Failure Beta') + theme(aspect.ratio = 1)

#ggplot(melt(data.frame(contrast=contrasts), id.vars=NULL), aes(x=value)) + geom_histogram(bins=50)

if (FU) {
    betas_BLFU <- array(dim=c(198,8))
    contrasts_BLFU <- array(dim=c(198,2))

    cols <- c('design_matGO_SUCCESS','design_matSTOP_SUCCESS','design_matSTOP_FAILURE','design_matGO_TOO_LATE')
    fus_lost <- 0
    for (subj_num in 1:50) {
        if (subj_num %in% subs_lost) {
            fus_lost <- fus_lost + 1
            next
        }
        betas_BLFU[subj_num,1:4] <- lm_list[[subj_num]]$lm$coefficients[cols]
        betas_BLFU[subj_num,5:8] <- lm_list_fu2[[subj_num-fus_lost]]$lm$coefficients[cols]

        # Contrast go conds with stop conds ignoring 'failures'
        #contr <- as.vector(c(-1, 0.5, 0.5))
        #contrasts[subj_num,1] <- contr %*% as.vector(betas[subj_num,1:3])
    }
    betas_BLFU <- data.frame(betas_BLFU)
    colnames(betas_BLFU) <- c('GO_SUCCESS','STOP_SUCCESS','STOP_FAILURE','GO_TOO_LATE','GO_SUCCESS_FU2','STOP_SUCCESS_FU2','STOP_FAILURE_FU2','GO_TOO_LATE_FU2')
}

###

dmelt <- melt(betas_BLFU[c('STOP_SUCCESS', 'STOP_SUCCESS_FU2')], id.vars = 'STOP_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(STOP_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_y_continuous('Stop Success Beta FU2') + scale_x_continuous('Stop Success Beta') + theme(aspect.ratio = 1)

dmelt <- melt(betas_BLFU[c('GO_SUCCESS', 'GO_SUCCESS_FU2')], id.vars = 'GO_SUCCESS', na.rm=TRUE)
ggplot(dmelt, aes(GO_SUCCESS, value)) + geom_point(colour='dark cyan') + scale_x_continuous('GO Success Beta') + scale_y_continuous('Go Success Beta FU2') + theme(aspect.ratio = 1)

dmelt <- melt(betas_BLFU[c('STOP_FAILURE', 'STOP_FAILURE_FU2')], id.vars = 'STOP_FAILURE', na.rm=TRUE)
ggplot(dmelt, aes(STOP_FAILURE, value)) + geom_point(colour='dark cyan') + scale_x_continuous('Stop Failure Beta') + scale_y_continuous('Stop Failure Beta FU2') + theme(aspect.ratio = 1)








# This is ok for now, but some information about the standard deviatins should also be encoded.
# I'm uncertain what at this point, but distributional overlap seems important.
imp_ind <- (data$train['mu_go_14'] - data$train['mu_stop_14']) / data$train['sigma_go_14']
colnames(imp_ind) <- c('imp_ind')

options(repr.plot.width=4, repr.plot.height=2)
dmelt  <- melt(imp_ind, id.vars = NULL)

probs  <- c(0.25, 0.5, 0.75)
quants <- quantile(imp_ind, prob=probs, na.rm=TRUE)

ggplot(dmelt, aes(x = value)) +
   facet_wrap(~variable, scales = "free") +
   geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +
   geom_histogram(aes(y = ..density..), alpha = 0.4) +
   scale_colour_manual(name = 'Density', values = c('red', 'blue')) +
   theme(legend.position = c(0.85, 0.85)) +
   geom_vline(data=dmelt, aes(xintercept=quants[1],), linetype="dashed", size=0.5) +
   geom_vline(data=dmelt, aes(xintercept=quants[2],), linetype="dashed", size=0.5) +
   geom_vline(data=dmelt, aes(xintercept=quants[3],), linetype="dashed", size=0.5)

dmelt  <- melt(data$raw['stn_std'], id.vars = NULL)

probs  <- c(0.25, 0.5, 0.75)
stn_quants <- quantile(data$raw['stn_std'], prob=probs, na.rm=TRUE)

options(repr.plot.width=4, repr.plot.height=2)
ggplot(dmelt, aes(x = value)) +
   facet_wrap(~variable, scales = "free") +
   geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density', na.rm=TRUE) +
   geom_histogram(aes(y = ..density..), alpha = 0.4, bins=30, na.rm=TRUE) +
   scale_colour_manual(name = 'Density', values = c('red', 'blue')) +
   theme(legend.position = c(0.85, 0.85)) +
   geom_vline(data=dmelt, aes(xintercept=stn_quants[1],), linetype="dashed", size=0.5) +
   geom_vline(data=dmelt, aes(xintercept=stn_quants[2],), linetype="dashed", size=0.5) +
   geom_vline(data=dmelt, aes(xintercept=stn_quants[3],), linetype="dashed", size=0.5) +
   ggtitle("STN Activation StdDevs & Quartiles")





  ##################################################
  # fmri stuff following what's in fit_fmri_glm
  ##################################################

  # ------------------------- #
  # Actual processing is done #
  # ------------------------- #
  # Now stuff is just getting printed
  lm_desc <- summary(linear_model)
  coeff   <- lm_desc$coefficients

  # STOP_SUCCESS is coeff. 4 originally (cond 3)
  stop_ind <- 4
  if (2 %in% bad_cond) {
    stop_ind <- stop_ind - 1
  }

  # Information
  print(paste('Subject ', subj_num,
            ' GO_SUCCESS beta: ', round(coeff[2,'Estimate'], digits=2),
            '. Sig: ', coeff[2,'Pr(>|t|)'] < 0.05,
            ' STOP_SUCCESS beta: ', round(coeff[stop_ind,'Estimate'], digits=2),
            '. Sig: ', coeff[stop_ind,'Pr(>|t|)'] < 0.05, sep=''))

  if (coeff[2,'Pr(>|t|)'] < 0.05) {
    frac_go_sig <- frac_go_sig + 1
  }
  if (coeff[stop_ind,'Pr(>|t|)'] < 0.05) {
    frac_stop_sig <- frac_stop_sig + 1
  }

  # Report on betas:
  #frac_go_sig   <- frac_go_sig  /n_to_fit
  #frac_stop_sig <- frac_stop_sig/n_to_fit

  print('')
  print(paste('Fraction of GO_SUCCESS   betas which are significant:', frac_go_sig))
  print(paste('Fraction of STOP_SUCCESS betas which are significant:', frac_stop_sig))

  return (lm_list)







  replace_bad_names <- function(data){

  dict <- list(
    c('', ''),
    c('', ''),
    c('IQ_PR_14', ''),
    c('IQ_VC_14', ''),
    c('GB_SSRT_14', ''),
    c('agn_mean_correct_latency_negative_14', ''),
    c('agn_mean_correct_latency_positive_14', ''),
    c('agn_total_omissions_negative_14', ''),
    c('agn_total_omissions_positive_14', ''),
    c('cgt_delay_aversion_14', ''),
    c('cgt_deliberation_time_14', ''),
    c('cgt_quality_of_decision_making_14', ''),
    c('cgt_overall_proportion_bet_14', ''),
    c('cgt_risk_adjustment_14', ''),
    c('cgt_risk_taking_14', ''),
    c('prm_percent_correct_14', ''),
    c('rvp_a_14', ''),
    c('swm_between_errors_14', ''),
    c('swm_strategy_14', ''),
    c('log10.k._14', ''),
    c('mu_targ_dur_co_14', ''),  # MID task parameters
    c('mu_targ_left_co_14', ''),
    c('mu_rewarded_co_14', ''),
    c('mu_high_rewarded_co_14', ''),
    c('mu_std_inv_rt_14', ''),
    c('mu_int_14', ''),
    c('sig_targ_dur_co_14', ''),
    c('sig_targ_left_co_14', ''),
    c('sig_rewarded_co_14', ''),
    c('sig_high_rewarded_co_14', ''),
    c('sig_std_inv_rt_14', ''),
    c('sig_int_14' '')
    )
  # Task data at age 14:
  features_14_task_raw <- c('IQ_PR_14', 'IQ_VC_14', 'GB_SSRT_14',
      'agn_mean_correct_latency_negative_14',
      'agn_mean_correct_latency_positive_14',
      'agn_total_omissions_negative_14',
      'agn_total_omissions_positive_14',
      'cgt_delay_aversion_14',
      'cgt_deliberation_time_14',
      'cgt_quality_of_decision_making_14',
      'cgt_overall_proportion_bet_14',
      'cgt_risk_adjustment_14',
      'cgt_risk_taking_14',
      'prm_percent_correct_14',
      'rvp_a_14',
      'swm_between_errors_14',
      'swm_strategy_14',
      'log10.k._14', setdiff(colnames(sst_params),'Subject'), setdiff(colnames(mid_params),'Subject'))


   new_sst_names <- c('mu_go_14', 'mu_stop_14', 'sigma_go_14','sigma_stop_14','tau_go_14','tau_stop_14','p_tf_14')
   new_mid_names <- c('mu_targ_dur_co_14', 'mu_targ_left_co_14',  'mu_rewarded_co_14',
                      'mu_high_rewarded_co_14', 'mu_std_inv_rt_14', 'mu_int_14',
                      'sig_targ_dur_co_14','sig_targ_left_co_14','sig_rewarded_co_14',
                      'sig_high_rewarded_co_14','sig_std_inv_rt_14','sig_int_14')

   new_agn_names <- c('agn_mean_corr_lat_neg_14', 'agn_mean_corr_lat_pos_14',
                      'agn_num_omis_neg_14', 'agn_num_omis_pos_14')

   new_cgt_names <- c('cgt_delay_avers_14','cgt_delib_14','cgt_quality_14','cgt_prop_bet_14',
                      'cgt_risk_adjust_14','cgt_risk_taking_14')

   new_espad_names <- c('alc_14', 'nic_14', 'amphet_14', 'coke_14', 'crack_14', 'ghb_14',
                        'glue_14', 'hash_14', 'ketamine_14', 'lsd_14', 'mushrooms_14', 'narc_14')

   features_14_task <- c('IQ_PR_14', 'IQ_VC_14', 'GB_SSRT_14', new_agn_names, new_cgt_names,
                         'prm_perc_corr_14', 'rvp_a_14', 'swm_btwn_errs_14', 'swm_strategy_14',
                         'log10.k._14', new_sst_names, new_mid_names)

   misc_task_names <- c('IQ_PR_14', 'IQ_VC_14', 'GB_SSRT_14', 'prm_perc_corr_14', 'rvp_a_14',
                        'swm_btwn_errs_14', 'swm_strategy_14', 'log10.k._14')

    features_14_survey_raw <- c('Sex_best_M0_14', 'PDS_14','All_Alc_14', 'All_Nic_14',
           'espad_life_amphet_14',
           'espad_life_coke_14',
           'espad_life_crack_14',
           'espad_life_ghb_14',
           'espad_life_glue_14',
           'espad_life_hash_14',
           'espad_life_ketamine_14',
           'espad_life_lsd_14',
           'espad_life_mushrooms_14',
           'espad_life_narcotic_14')

   # Features left out:
   #       'espad_life_anabolic_14',
   #       'espad_life_heroin_14',
   #       'espad_life_mdma_14',
   #       'espad_life_tranq_14'

   features_14_survey <- c('bio_sex_14', 'pds_14', new_espad_names)

   features_18_task <- c('log10.k._18', 'agn_mean_correct_latency_negative_18', 'agn_mean_correct_latency_neutral_18',
     'agn_mean_correct_latency_positive_18', 'agn_total_omissions_negative_18', 'agn_total_omissions_neutral_18',
     'agn_total_omissions_positive_18', 'cgt_delay_aversion_18', 'cgt_deliberation_time_18',
     'cgt_overall_proportion_bet_18', 'cgt_quality_of_decision_making_18', 'cgt_risk_adjustment_18',
     'cgt_risk_taking_18', 'prm_percent_correct_18', 'rvp_a_18', 'swm_between_errors_18', 'swm_strategy_18')

    #features_18_survey <- c('age_18', 'X6.life.nic_18', 'X8a.life.alc_18', 'Life.amph_18', 'Life.anab_18',
    #    'Life.coke_18', 'Life.crack_18', 'Life.hash.thc_18', 'Life.heroin_18', 'Life.GHB_18', 'Life.glue_18',
    #    'Life.ketamine_18', 'Life.lsd_18', 'Life.MDMA_18', 'Life.mushrooms_18', 'Lif.narc_18', 'Life.tranq_18')

    features_18_survey <- c('X6.life.nic_18', 'X8a.life.alc_18', 'Life.amph_18', 'Life.anab_18',
        'Life.coke_18', 'Life.crack_18', 'Life.hash.thc_18', 'Life.heroin_18',
        'Life.ketamine_18', 'Life.lsd_18', 'Life.MDMA_18', 'Life.mushrooms_18', 'Lif.narc_18')


    n_feat_task   <- length(features_14_task)
    n_feat_survey <- length(features_14_survey)

    for (i in 1:n_feat_task){
        old_name <- features_14_task_raw[[i]]
        new_name <- features_14_task[[i]]
        colnames(raw_df)[colnames(raw_df) == old_name] <- new_name
    }
    for (i in 1:n_feat_survey){
        old_name <- features_14_survey_raw[[i]]
        new_name <- features_14_survey[[i]]
        names(raw_df)[names(raw_df) == old_name] <- new_name
    }


    output <- list('raw' = raw_df, 'survey_names_14' = features_14_survey,
                   'task_names_14' = features_14_task, 'sst_names' = new_sst_names,
                   'mid_names' = new_mid_names, 'espad_names' = new_espad_names,
                   'cgt_names' = new_cgt_names, 'agn_names' = new_agn_names,
                   'msc_names' = misc_task_names,
                   'task_names_18' = features_18_task,
                   'survey_names_18' = features_18_survey
                   )
    return(output)
}
# ------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------ #
#
# ------------------------------------------------------------------------------ #
somethingorother <- function () {
  # Remove the wierd negative one values from the continuous task variables.
  data$raw[data$task_names_14][data$raw[data$task_names_14] == -1] <- NA
  data$raw[data$task_names_18][data$raw[data$task_names_18] == -1] <- NA


  cat('Hence, the following features are retained:\n')
  print(colnames(data$raw))

}
# ------------------------------------------------------------------------------ #


group_params <- c('mu_go', 'mu_go_var', 'mu_stop', 'mu_stop_var', 'tau_go', 'tau_stop', 'tau_go_var',
                  'tau_stop_var', 'sigma_go', 'sigma_stop', 'sigma_go_var', 'sigma_stop_var', 'p_tf')

dir_lists<- c( paste('/home/dan/documents/lncc/',
                    c('sst_data_train_020/', 'sst_data_train_050/',
                      'sst_data_train_100/', 'sst_data_train_150/',
                      'sst_data_train/'), sep=''))

param_file_names <- c( paste( paste('parameters',c('1','2','3','4'), sep='') , '.csv', sep=''))


stats <- array(NA,c(13,4,5))

ndirs   <- length(dir_lists)
npfiles <- length(param_file_names)

for (dir_num in 1:ndirs) {
    dir <- dir_lists[dir_num]

    for (file_num in 1:npfiles) {
        file <- param_file_names[file_num]

        data <- read.csv(paste(dir, file, sep=''), sep=';') #, row.names=group_params)

        data  <- data[group_params]
        means <- as.vector(colMeans(data))

        stats[ , file_num, dir_num] <- means
    }
}

summary <- array(NA, c(5,13))
summary[1,] <- apply(stats[,,1], 1, mean)
summary[2,] <- apply(stats[,,2], 1, mean)
summary[3,] <- apply(stats[,,3], 1, mean)
summary[4,] <- apply(stats[,,4], 1, mean)
summary[5,] <- apply(stats[,,5], 1, mean)

library(ggplot2)
library(reshape2)

## V Coefficients Plot:
nlines = 5
plotchar <- seq(18,18+nlines,1)

xrange <- range(1:5)
yrange <- range(summary)

op <- par(mar=c(10.1, 4.1, 4.1, 2.1))
#plot(xrange, yrange, type="n", xlab="", ylab="", xaxt="n")
