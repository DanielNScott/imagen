

write_stim_1D <- function(task_data, dir) {

  # Possible conditions are outer(-, -) of these.
  events <- c('GO_SUCCESS', 'STOP_SUCCESS', 'STOP_FAILURE', 'GO_TOO_LATE',
             'GO_WRONG_KEY_RESPONSE', 'GO_FAILURE', 'STOP_TOO_EARLY_RESPONSE')
  hands  <- c('LeftArrow', 'RightArrow')

  # Define which conditions to collapse across hands
  collapse = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

  # Determine which conditions get files.
  conditions <- c(outer(events[!collapse], hands, paste), events[collapse])

  # Write a file for everything in outer(), modulo any collapsing
  for (cond in conditions) {

    # Figure out if we have a handedness split for this event and set mask.
    split      <- strsplit(cond, split = ' ')[[1]]
    stim       <- split[1]
    task_stims <- task_data['Stimulus.Presented']

    if (length(split) == 2) {
      hand      <- split[2]
      hand_msk  <- task_stims == hand
      hand_abrv <- ifelse(hand == 'LeftArrow', '_l', '_r')
    } else {
      # Hand mask should always be TRUE here.
      hand_msk  <- task_stims == hands[1] || task_stims == hands[2]
      hand_abrv <- ''
      assertthat::assert_that(hand_msk)
    }

    # Select based on mask, etc.
    mask  <- (task_data['Response.Outcome'] == stim) & (hand_msk)
    times <- t(task_data[mask, 'Trial.Start.Time..Onset.']) / 1000
    fname <- paste(dir, '/1Ds/', tolower(stim), hand_abrv, '.1D', sep = '')

    if (length(times) != 0) {
      # Write table doesn't work here becauase it inserts row indexing...
      write(times, file = fname, ncolumns = 100000)

    } else {
      # This is the convention in afni...
      write(x = '*', file = fname)
    }
  }
}


# Sets up the */afni/ file structure so fits can be created
setup_afni <- function(loc, subj_ids) {

  # Different base directories for cluster or local use
  if (loc == 'local') { base_dir <- '/home/dan/projects/imagen/'
  } else {              base_dir <- '/users/dscott3/projects/imagen/'}

  visit <- 'BL'
  task  <- 'SST'

  if (loc == 'local') {
    # For cluster use:
    taskdata_dir   <- paste(base_dir, '/data/fmri/', sep = '')
    timeseries_dir <- paste(base_dir, '/data/fmri/', sep = '')
    movement_dir   <- paste(base_dir, '/data/fmri/', sep = '')
  } else {
    # For local testing:
    taskdata_dir   <- paste(base_dir, '/data/BL_SST_task/', sep = '')
    timeseries_dir <- paste(base_dir, '/data/BL_SST_AFNI/', sep = '')
    movement_dir   <- paste(base_dir, '/data/BL_SST_move/', sep = '')
  }

  # Setup subject everything
  for (subj_id in subj_ids[[1]]) {
    # (subj_ids is a data frame so needs conversion)

    # Convert subject id to appropriate string
    subj_id_str <- formatC(subj_id, width = 12, format = 'd', flag = '0')

    # Set filename prefixes and suffixes
    task_pfx <- 'ss_'
    fmri_sfx <- '_BL_SST.nii.gz'
    move_pfx <- 'EPI_stop_'

    # Filename beurocracy...
    fmri_fname <- paste(subj_id_str, fmri_sfx, sep = '')
    task_fname <- paste(task_pfx, subj_id_str, '.csv', sep = '')
    move_fname <- paste(move_pfx, subj_id_str, '.txt', sep = '')

    fmri_full_name <- paste(timeseries_dir, fmri_fname, sep = '')
    task_full_name <- paste(taskdata_dir  , task_fname, sep = '')
    move_full_name <- paste(movement_dir  , move_fname, sep = '')

    # Read the SST series
    if (file.exists(task_full_name)) {
      # Only want to read 2 columns, time and outcome
      cols <- c('NULL','NULL','numeric','NULL','NULL',NA    ,'NULL','NULL',
                'NULL','NULL','NULL'   , NA   ,'NULL','NULL','NULL','numeric')
      task_data <- read.csv(file = task_full_name, header = TRUE, sep = '\t', colClasses = cols, skip = 1)
      print(paste('Successful read of ', task_fname))
    } else {
      print(paste('Failure to read ', task_fname))
      next
    }

    # Setup names of stuff and softlink command
    afni_dir   <- paste(base_dir, '/afni/', sep = '')
    subj_dir   <- paste(afni_dir, subj_id_str, '_DATA', sep = '')
    subj_1Ds   <- paste(subj_dir, '/1Ds', sep = '')

    # Name subject directories, and write symlink commands
    proc_fname <- 'process_subj.tcsh'
    link_cmd1  <- paste('ln -s ', afni_dir, proc_fname, ' ', subj_dir, '/', proc_fname, sep = '')
    link_cmd2  <- paste('ln -s ', base_dir, '/data/ROI_masks', ' ', subj_dir, '/', 'ROI_masks', sep = '')
    link_cmd3  <- paste('ln -s ', move_full_name, ' ', subj_dir, '/1Ds/motion_demean.1D', sep = '')
    link_cmd4  <- paste('ln -s ', fmri_full_name, ' ', subj_dir, '/', sep = '')

    # Execute all that
    dir.create(subj_dir)
    dir.create(subj_1Ds)
    system(link_cmd1)
    system(link_cmd2)
    system(link_cmd3)
    system(link_cmd4)

    # Write regression files
    write_stim_1D(task_data, subj_dir)
  }
}


