# Ancestral graph analysis for the stop signal task.
sst_ag_analysis <- function(ag_analysis_dir, data) {

  # Where are the ancestral graph scripts
  ag_source_dir <- '/home/dan/projects/ancestral-graphs/AG_codes'

  # source AG_codes
  file.sources <- list.files(paste(ag_source_dir, '/rFunctions2', sep = ''), pattern = "*.R", full.names = TRUE, recursive = TRUE)
  sapply(file.sources, source, .GlobalEnv)

  # Define indices for the coefficient data
  source(paste(ag_analysis_dir, '/sst_ag_index_data.r', sep = ''))
  indices <- sst_ag_index_data('/home/dan/projects/imagen/data')

  # these are the conditions that are evalauted seperatly in sst_ag_model_fitting.r
  cond <- list(indices$ST, indices$SR, indices$Go)
  names(cond) = c('ST', 'SR', 'Go')

  # select the graph nodes for connectivity network
  #M0roi=c("CaudateR40exc" ,"PreSMARsmall","IFGR", "maxSTNR25exc","maxGPeR30exc","maxGPiR30exc","ThalamusR40exc")
  #M0roi <- c("rCaudate", "rPreSMA", "rIFG", "rSTN", "rGPe", "rGPi", "rThalamus")
  M0roi <- c("rCaudate", "rPreSMA", "rIFG", "rGPe", "rGPi", "rThalamus")
  C.Lab <- M0roi

  # source the models to use, these are the defined PFC-BG models to evaluate for connectivity/fits
  source(paste(ag_analysis_dir, '/define_ag_models.r', sep = ''))
  models <- define_ag_models()

  # Now do the connectivity evaluation for each subject, each model/condition, to compute aic/bic and n-fits
  source(paste(ag_analysis_dir, '/sst_ag_model_fitting.r', sep = ''), chdir = F)
  model_fits <- fit_ag_models(data, models)

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
  cond2 <- list(ST, SR)
  names(cond2) <- c('ST', 'SR')

  # Make succesfull stop and failed stop lists using data
  ST <- list()
  SR <- list()

  n_subj <- length(data)

  for (subj_num in 1:n_subj){
    ST[[subj_num]] <- data[[subj_num]][cond2[[1]][[subj_num]],C.Lab]
    SR[[subj_num]] <- data[[subj_num]][cond2[[2]][[subj_num]],C.Lab]
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
