index_ag_data <- function(data) {

  # Initialize lists of trial numbers for diff conds.
  conds <- list( list(), list(), list(),
                 list(), list(), list(),
                 list(), list(), list())

  names(conds) <- c('Go', 'GoL', 'GoR', 'St', 'StL', 'StR',
                    'Sr', 'SrL', 'SrR')

  # TODO: lapply-ify this loop.
  #patterns <- c('GO_SUCCESS', 'STOP_SUCCESS', 'STOP_FAILURE')
  # For the example data ...
  #patterns <- names(conds)

  # Populate lists
  for (subj in 1:length(data)) {
    trial_names <- rownames(data[[subj]])

    conds$Go[[subj]]    <- grep('GO_SUCCESS'  , trial_names)
    conds$St[[subj]]    <- grep('STOP_SUCCESS', trial_names)
    conds$Sr[[subj]]    <- grep('STOP_FAILURE', trial_names)

    conds$GoR[[subj]]   <- grep('GO_SUCCESS_RIGHT'  , trial_names)
    conds$StR[[subj]]   <- grep('STOP_SUCCESS_RIGHT', trial_names)
    conds$SrR[[subj]]   <- grep('STOP_FAILURE_RIGHT', trial_names)

    conds$GoL[[subj]]   <- grep('GO_SUCCESS_LEFT'  , trial_names)
    conds$StL[[subj]]   <- grep('STOP_SUCCESS_LEFT', trial_names)
    conds$SrL[[subj]]   <- grep('STOP_FAILURE_LEFT', trial_names)

    # Code for the example data...
    #conds$Go[[subj]]    <- grep('Go', trial_names)
    #conds$St[[subj]]    <- grep('St', trial_names)
    #conds$Sr[[subj]]    <- grep('Sr', trial_names)

    #conds$GoR[[subj]]   <- grep('GoR', trial_names)
    #conds$StR[[subj]]   <- grep('StR', trial_names)
    #conds$SrR[[subj]]   <- grep('SrR', trial_names)

    #conds$GoL[[subj]]   <- grep('GoL', trial_names)
    #conds$StL[[subj]]   <- grep('StL', trial_names)
    #conds$SrL[[subj]]   <- grep('SrL', trial_names)
  }

  # Index by subject id.
  conds <- lapply(conds, function(x){names(x) <- names(data); return(x)} )

  return(conds)
}
