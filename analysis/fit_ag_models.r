fit_ag_models <- function(data,  models, cond, rois) {
  # Fits the ancestral graph models
  #
  # Args:
  #   data:
  #   models:
  #   cond:
  #   rois:
  #
  # Returns:
  #   Fit statistics and model edge weights.

  n_models <- length(models)
  n_subj   <- length(data)
  n_conds  <- length(cond)
  subj_ids <- names(data)

  # Most of the matrices below have the same attributes
  setup_matrix <- function(name, suffixes = 1:n_models, n_rows = n_subj, n_cols = n_models, other_names = c()) {
    matrix <- round(matrix(, n_rows, n_cols), digits=2)

    colnames(matrix) <- c(other_names, paste(name, suffixes, sep=''))
    rownames(matrix) <- subj_ids

    return(matrix)
  }

  # prepare matrix to save AIC values for individual subject
  AIC_direct <- setup_matrix('aic')

  # prepare matrix to save chi (yuan bentler chi test) for individuals
  Fit_direct <- setup_matrix('chi')

  # save P-val chi test,  if p<0.05 model is no fit for subject compared to saturated model
  P_direct   <- setup_matrix('chi-p')

  # make matrix to save log-likelihood per models en per ppn
  logL_direct <- setup_matrix('logL')

  # save number of observations for each subject to compute BIC
  nObs <- setup_matrix('nobs')

  # save individual AIC values
  #AMT <- setup_matrix('nfit-', suffixes = names(models), n_rows = n_conds, n_cols = n_models*2, other_names = names(models))
  AMT <- round(matrix(,n_conds, n_models*2), digits=2)
  colnames(AMT) <- c(names(models), paste('nfit-', names(models), sep=''))
  rownames(AMT) <- names(cond)

  # save individual BIC values
  #BMT <- setup_matrix('nfit-', suffixes = names(models), n_rows = n_conds, n_cols = n_models*2, other_names = names(models))
  BMT <- round(matrix(,n_conds, n_models*2), digits=2)
  colnames(AMT) <- c(names(models), paste('nfit-', names(models), sep=''))
  rownames(AMT) <- names(cond)


  for (cond_num in 1:n_conds) {

    for (subj_num in 1:n_subj) {
        Aic   <- c()
        Fit   <- c()
        P     <- c()
        nPars <- c() # vecotr met nr par voor ieder models
        nobs  <- c() # vector with  number of observations for each model

        # var models contains the models
        # data contains the data for all subjects in a list
        # rois indicates the Region of input for the models
        # cond contains the index numbers for each Condition,  each subject in data

        for (model_num in 1:n_models) {
          model     <- models[[model_num]]
          subj_data <- data[subj_num]
          condition <- cond[[cond_num]][[subj_num]]

          fitted_ag <- fitAG(subj_data,  model, rois, condition, detail = 'both')

          Aic[  model_num] <- fitted_ag$aic
          nPars[model_num] <- 0 #fitted_ag$npars
          Fit[  model_num] <- fitted_ag$fit$chi2
          P[    model_num] <- fitted_ag$fit$p

          nobs[model_num]  <- length(cond[[cond_num]][[subj_num]])
        }

    AIC_direct[subj_num, ]  <- round(Aic, digits=2)
    Fit_direct[subj_num, ]  <- round(Fit, digits=2)
    P_direct[   subj_num, ] <- round(P, digits=2)
    logL_direct[subj_num, ] <- Aic - 2*nPars

    nObs[subj_num, ] <- round(nobs[subj_num], digits=2)

    }

    A = cbind(AIC_direct, P_direct)

    # Radom effects AIC with pooled log likelihood over subjects
    AIC.direct <- apply(logL_direct, 2, sum)+2*nPars

    # Random effects BIC with pooled log likelihood over subjects
    BIC.totaal <- apply(logL_direct, 2, sum)+(nPars*log(apply(nObs, 2, sum)))

    # Get the number of subjects where the model fits
    nfit.rfx   <- apply(ifelse(A[, (n_models+1):(2*n_models)]>0.049, 1, 0), 2, sum)

    # Save the random effects AIC and BIC but also individual subjects
    # (from chi test) where the model did not fit (key!)
    AMT[cond_num, ] <- c(AIC.direct, nfit.rfx)
    BMT[cond_num, ] <- c(BIC.totaal, nfit.rfx)

  } # close cond loop

  return(list(AMT = AMT, BMT = BMT))
}
