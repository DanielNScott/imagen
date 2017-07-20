ag_analysis <- function(dir_analy, data) {
  # Ancestral graph analysis for the stop signal task.
  #
  # Args:
  #   dir_analy: Directory with the analysis source files
  #   data: Array-like activation data with trials as rows, ROIs as cols.
  #
  # Returns:
  #
  print('Sourcing code and setting up models...')

  # Where are the ancestral graph scripts
  dir_src <- '/users/dscott3/projects/ancestral-graphs/AG_codes/src'

  # Source the general AG codes
  file.sources <- list.files(dir_src, pattern = "*.R", full.names = TRUE)
  sapply(file.sources, source, .GlobalEnv)

  source(paste(dir_analy, '/sst_ag_index_data.r'   , sep = ''))
  source(paste(dir_analy, '/define_ag_models.r'    , sep = ''))
  source(paste(dir_analy, '/sst_ag_model_fitting.r', sep = ''))

  # Define indices for the coefficient data
  indices <- sst_ag_index_data('/users/dscott3/projects/imagen/data')

  # Define conditions to be evalauted seperatly in sst_ag_model_fitting.r
  cond <- list(indices$ST, indices$SR, indices$Go)
  names(cond) = c('ST', 'SR', 'Go')

  # Select the graph nodes for the connectivity networks to examine and source models of them
  #M0roi=c("CaudateR40exc" ,"PreSMARsmall","IFGR", "maxSTNR25exc","maxGPeR30exc","maxGPiR30exc","ThalamusR40exc")
  rois   <- c("rCaudate", "rPreSMA", "rIFG", "rSTN", "rGPe", "rGPi", "rThalamus")
  models <- define_ag_models()


  # Now do the connectivity evaluation for each subject, each model/condition, to compute aic/bic and n-fits
  print('Fitting models...')
  model_fits <- fit_ag_models(data, models, cond, rois)

  # Print results based on random effects aic and bic, + number of fits based on yuan
  # and bentler chi (n = number of subject where model fits)
  model_fits$AMT
  model_fits$BMT


  #------------------------------------------------------------------------------#
  #         Compute (indvidual) connectivity strengths of winning model          #
  #------------------------------------------------------------------------------#

  # ag0 was the hindi model in StopModel.R
  hindi <- models$ag0

  #conditions to get beta for [these are the conditions where the model fitted all ppn]
  cond2 <- list(indices$ST, indices$SR)
  names(cond2) <- c('ST', 'SR')

  # Make succesfull stop and failed stop lists using data
  ST <- list()
  SR <- list()

  n_subj <- length(data)

  for (subj_num in 1:n_subj){
    ST[[subj_num]] <- data[[subj_num]][cond2[[1]][[subj_num]],rois]
    SR[[subj_num]] <- data[[subj_num]][cond2[[2]][[subj_num]],rois]
  }

  # compute covarience matrix
  covList.ST <- lapply(ST, cov)
  covList.SR <- lapply(SR, cov)

  fitST = list()
  fitSR = list()

  for (subj_num in 1:n_subj){
    fitST[[subj_num]] <- fitAncestralGraph(hindi, covList.ST[[subj_num]], dim(ST[[subj_num]])[1])
    fitSR[[subj_num]] <- fitAncestralGraph(hindi, covList.SR[[subj_num]], dim(SR[[subj_num]])[1])
    }

  # compute the variance of beta per subject
  var.g.ST <- ag.var.group(fitST, covList.ST, ST)
  var.g.SR <- ag.var.group(fitSR, covList.SR, SR)

  # save matrix for var beta and beta itself (zijn nog leeg)
  varbeta.ST <- matrix(,dim(var.g.ST)[3], dim(var.g.ST)[1])
  varbeta.SR <- matrix(,dim(var.g.SR)[3], dim(var.g.SR)[1])

  beta.ST    <- matrix(,n_subj, length(ag.theta(fitST[[1]])))
  beta.SR    <- matrix(,n_subj, length(ag.theta(fitSR[[1]])))

  for (subj_num in 1:n_subj){
    varbeta.ST[subj_num,] <- diag(var.g.ST[,,subj_num])
    varbeta.SR[subj_num,] <- diag(var.g.SR[,,subj_num])

    beta.ST[subj_num,]    <- ag.theta(fitST[[subj_num]])
    beta.SR[subj_num,]    <- ag.theta(fitSR[[subj_num]])
    }

  # these are the standarized connections (to use for group comparisons)
  sbeta.ST <- round(beta.ST / varbeta.ST, dig = 8)
  sbeta.SR <- round(beta.SR / varbeta.SR, dig = 8)

  # note here beta and sbeta contain the estimated for all the defined connections in the model

  #------------------------------------------------------------------------------#
  #           give names to the rows and columns of the beta matrix              #
  #------------------------------------------------------------------------------#

  # select the defined directed and undirected connections
  beta.ST  <- beta.ST[,  c(1:8,10)]
  beta.SR  <- beta.SR[,  c(1:8,10)]
  sbeta.ST <- sbeta.ST[, c(1:8,10)]
  sbeta.SR <- sbeta.SR[, c(1:8,10)]


  colnames(beta.ST) <- paste('ST-',c('caud->gpe','presma->caud','presma->stn','ifg->caud',
                                   'ifg->stn','stn->gpi','gpe->gpi','gpi->thalamus',
                                   'presma-ifg'),sep='')

  colnames(beta.SR) <- paste('SR-',c('caud->gpe','presma->caud','presma->stn','ifg->caud',
                                   'ifg->stn','stn->gpi','gpe->gpi','gpi->thalamus',
                                   'presma-ifg'),sep='')

  colnames(sbeta.ST) <- colnames(beta.ST)
  colnames(sbeta.SR) <- colnames(beta.SR)

  Beta_output <- list(Connectivity_ST = beta.ST, Connectivity_SR = beta.SR,
                      standarized_con_ST = sbeta.ST, standarized_con_SR = sbeta.SR)

  for (rn in 1:length(Beta_output)) {
    rownames(Beta_output[[rn]]) = names(data)
  }
}
