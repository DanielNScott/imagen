index_ag_data <- function(data) {

  # Initialize lists of trial numbers for diff conds.
  conds <- list( list(), list(), list(),
                 list(), list(), list(),
                 list(), list(), list())

  names(conds) <- c('Go', 'GoL', 'GoR', 'St', 'StL', 'StR',
                    'Sr', 'SrL', 'SrR')

  # Populate lists
  for (subj in 1:length(data)) {
    trial_names <- rownames(data[[subj]])

    conds$Go[[subj]]    <- grep('GO_SUCCESS'  , trial_names)
    conds$St[[subj]]    <- grep('STOP_SUCCESS', trial_names)
    conds$Sr[[subj]]    <- grep('STOP_FAILURE', trial_names)

    conds$GoR[[subj]]   <- grep('GO_SUCCESS_RIGHT'  , trial_names)
    conds$StR[[subj]]   <- grep('STOP_SUCCESS_RIGHT', trial_names)
    conds$SrR[[subj]]   <- grep('STOP_FAILURE_RIGHt', trial_names)

    conds$GoL[[subj]]   <- grep('GO_SUCCESS_LEFT'  , trial_names)
    conds$StL[[subj]]   <- grep('STOP_SUCCESS_LEFT', trial_names)
    conds$SrL[[subj]]   <- grep('STOP_FAILURE_LEFT', trial_names)
  }

  # Copy subject names
  names(Go)    <- names(data)
  names(ST)    <- names(data)
  names(SR)    <- names(data)

  names(GoR)    <- names(data)
  names(STR)    <- names(data)
  names(SRR)    <- names(data)

  names(GoL)    <- names(data)
  names(STL)    <- names(data)
  names(SRL)    <- names(data)
  #names(error) <- names(data)
  #names(miss)  <- names(data)

  #return(list(Go = Go, ST = ST, SR = SR))
  return(conds)
}
