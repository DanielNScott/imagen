#-----------------------------------------------------------------------#
write_stim_1D <- function(task_data, dir, ag = FALSE) {
  # Writes a set of AFNI .1D stimulus files for the SST.
  #
  # Args:
  #   task_data: A data frame with task timing and outcome info.
  #   dir: The directory in which to create and fill a subdir '1Ds'.
  #
  # Returns:
  #   Nothing

  # Possible conditions are outer(-, -) of these.
  events <- c('GO_SUCCESS', 'STOP_SUCCESS', 'STOP_FAILURE', 'GO_TOO_LATE',
             'GO_WRONG_KEY_RESPONSE', 'GO_FAILURE', 'STOP_TOO_EARLY_RESPONSE')
  hands  <- c('LeftArrow', 'RightArrow')

  # Define which conditions to make minimal regressors for.
  # 'collapse' collapses across hands AND turns off by-trial regressor file-writing.
  if (ag) {
    collapse <- c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)
  } else {
    collapse <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  }

  # Determine which conditions get files.
  conditions <- c(outer(events[!collapse], hands, paste), events[collapse])

  # Which establish which conditions get seperation of trials
  n_handed_conds  <- 2*length(events[!collapse])
  n_conds         <- length(conditions)
  seperate_trials <- c(rep_len(TRUE, n_handed_conds), rep_len(FALSE, n_conds - n_handed_conds))

  stim_text <- ''
  global_instance_num <- 1

  # Write a file for everything in outer(), modulo any collapsing
  for (cond in conditions) {

    # Figure out if we have a handedness split for this event and set mask.
    split      <- strsplit(cond, split = ' ')[[1]]
    stim       <- split[1]
    task_stims <- task_data['Stimulus.Presented']

    if (length(split) == 2) {
      hand      <- split[2]
      hand_msk  <- task_stims == hand
      hand_abrv <- ifelse(hand == 'LeftArrow', '_left', '_right')
    } else {
      hand_msk  <- task_stims == hands[1] || task_stims == hands[2]
      hand_abrv <- ''

      # Hand mask should always be TRUE here.
      assertthat::assert_that(hand_msk)
    }

    # Select based on mask, etc.
    mask  <- (task_data['Response.Outcome'] == stim) & (hand_msk)
    times <- t(task_data[mask, 'Trial.Start.Time..Onset.']) / 1000
    fname <- paste(dir, '/1Ds/', tolower(stim), hand_abrv, '.1D', sep = '')

    if (length(times) != 0) {
      # Write table doesn't work here becauase it inserts row indexing...
      write(times, file = fname, ncolumns = 100000)

      cond_index <- grep(cond, conditions)
      if (seperate_trials[cond_index]) {

        local_stim_num <- 0
        rep_set <- 1:length(times)
        for (instance in rep_set) {
          afni_stim_name  <- paste(tolower(stim), hand_abrv, sep = '')

          fname_inst <- paste(dir, '/1Ds/', afni_stim_name, '_', toString(instance), '.1D', sep = '')
          fname_rest <- paste(dir, '/1Ds/', afni_stim_name, '_', toString(instance), '_rest.1D', sep = '')

          this_time   <- times[instance]
          other_times <- times[setdiff(rep_set, instance)]
          write(this_time  , file = fname_inst, ncolumns = 1)
          write(other_times, file = fname_rest, ncolumns = 100000)

          if (cond == 'STOP_SUCCESS LeftArrow' || cond == 'STOP_FAILURE LeftArrow') {
            local_stim_num <- local_stim_num + 1
            stim_num   <- toString(13 + global_instance_num)
            stim_text  <- paste(stim_text,  ' -stim_file ', stim_num, ' ${dir1D}', afni_stim_name, '_',
                                local_stim_num, '.1D', ' \'SPMG1\' ', '-stim_label ', stim_num, ' ',
                                afni_stim_name, '_', local_stim_num ,' \\\\\n', sep = '')

            global_instance_num <- global_instance_num + 1
          }
        }
      }
    } else {
      # This is the convention in afni...
      write(x = '*', file = fname)
    }
  }

  # Fix up AG deconvolution file.
  deconv_file <- 'deconvolve.sh'
  deconv_text <- readChar(deconv_file, file.info(deconv_file)$size)

  deconv_text <- sub('NSTIM_ANCHOR\n', paste(' -num_stimts', toString(global_instance_num - 1), ' \\\\\n'), deconv_text)
  deconv_text <- sub('STIMS_ANCHOR\n', stim_text, deconv_text)

  write(x = deconv_text, file = deconv_file)

}
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
setup_afni <- function(loc, subj_ids) {
  # Sets up the */afni/ file structure so fits can be created

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

  # Begin writing shell script to submit everything
  write('#!/bin/bash \n', file = 'submitter.sh')

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
    afni_dir   <- paste(base_dir, '/data/afni/', sep = '')
    subj_dir   <- paste(afni_dir, subj_id_str, '_DATA', sep = '')
    subj_1Ds   <- paste(subj_dir, '/1Ds', sep = '')

    # Name subject directories, and write symlink commands
    proc_fname <- 'process_subj.sh'
    wipe_fname <- 'wipe_results.sh'
    link_cmd1  <- paste('ln -s ', afni_dir, proc_fname, ' ', subj_dir, '/', proc_fname, sep = '')
    link_cmd2  <- paste('ln -s ', base_dir, '/data/ROI_masks', ' ', subj_dir, '/', 'ROI_masks', sep = '')
    link_cmd3  <- paste('ln -s ', move_full_name, ' ', subj_dir, '/1Ds/motion_demean.1D', sep = '')
    link_cmd4  <- paste('ln -s ', fmri_full_name, ' ', subj_dir, '/BL_SST.nii.gz', sep = '')
    link_cmd5  <- paste('cp    ', afni_dir, 'deconvolve.sh', ' ', subj_dir, '/', deconvolve.sh, sep = '')
    link_cmd6  <- paste('ln -s ', afni_dir, 'extract.sh', ' ', subj_dir, '/', 'extract.sh', sep = '')

    # Execute all that
    dir.create(subj_dir)
    dir.create(subj_1Ds)
    system(link_cmd1)
    system(link_cmd2)
    system(link_cmd3)
    system(link_cmd4)
    #system(link_cmd5)
    system(link_cmd6)

    # Write regression files
    write_stim_1D(task_data, subj_dir, ag = TRUE)

    jobname  <- paste(subj_id_str, 'glm', sep = '_')
    cdcmd1   <- paste('cd ', subj_id_str, '_DATA;', sep = '')
    cdcmd2   <- 'cd ../; sleep 1;'
    srunpfx  <- 'sbatch -t 0:05:00 -n 1 --nodes 1 --cpus-per-task 1 -J'
    srun_cmd <- paste(cdcmd1, srunpfx, jobname, 'process_subj.sh ;', cdcmd2, sep = ' ')

    write(srun_cmd, file = 'submitter.sh', append = TRUE)
  }
}
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
read_stats_dump <- function(fname) {
  roi_stats <- read.table(fname, sep = ' ', header = FALSE, blank.lines.skip = FALSE)
  paste_fun  <- function(x,y) {paste(x, y, sep = '_')}
  str_outer  <- function(x,y) {as.vector(t(outer(x, y, paste_fun)))}

  types  <- c('hrf', 'hrf_dt', 'hrf_dd')
  stats  <- c('beta', 'tval')
  conds  <- c('gs', 'ss', 'sf', 'gf', 'gtl', 'gwk', 'gte')
  ctrsts <- c('ss_go', 'ss_sf')

  stat_names <- c('all_Fval')
  for (cond in conds) {
    # Setup "cond_type_stat" fields
    block <- str_outer(cond, types)
    block <- str_outer(block, stats)

    # Append the condition level F-value and save
    block <- c(block, paste_fun(cond, 'Fval'))
    stat_names <- c(stat_names, block)
  }
  for (ctrst in ctrsts) {
    block <- str_outer(ctrst, stats)
    stat_names <- c(stat_names, block, paste_fun(ctrst, 'Fval'))
  }

  colnames(roi_stats) <- stat_names
  roi_stats <- t(roi_stats)
  return(roi_stats)
}
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
read_all_stats <- function(subj_ids, afni_dir) {

  #roi_mask_names <- c('L_AAL_ACC', 'R_AAL_ACC', 'R_Caudate_AAL',
  #                     'R_GPe', 'R_GPi', 'R_IFG', 'R_NAcc', 'R_preSMA',
  #                     'R_STN', 'R_Thalamus_AAL')
  rois <- c('rIFG', 'rPreSMA', 'rCaudate', 'rGPe', 'rGPi',
            'rSTN', 'rThalamus', 'rNAcc', 'lACC', 'rACC')
  cois <- c('gs_hrf_beta', 'ss_hrf_beta', 'sf_hrf_beta', 'ss_go_beta', 'ss_sf_beta')
  cond <- c('go', 'st', 'sr', 'st_go', 'st_sr')

  paste_fun   <- function(x,y) {paste(x, y, sep = '_')}
  roi_by_cond <- as.vector(outer(rois, cond, paste_fun))

  results <- data.frame('Subject' = subj_ids)

  std_fmri_feats <- as.vector(outer(rois, cond, function(x,y) {paste(x,y, sep = '')}))

  # Cycle through subjects directories, doing function "fn"
  merge_in <- TRUE
  for (subj_id in subj_ids[[1]]) {
    subj_id_str <- formatC(subj_id, width = 12, format = 'd', flag = '0')

    # Setup names of stuff
    subj_dir <- paste(afni_dir, subj_id_str, '_DATA/', sep = '')

    stats <- tryCatch({
      stats <- read_stats_dump(paste(subj_dir, '/stats_masked.out', sep=''))
    }, error = function(e){
      print(paste('Unable to read file stats for', subj_id_str))
      return(TRUE)
    })

    if (is.logical(stats)) {next}

    mean_betas <- stats[cois,]
    colnames(mean_betas) <- rois
    for (cnum in 1:length(cois)) {
      condi <- cois[cnum]
      condo <- cond[cnum]

      cond_frame <- data.frame('Subject' = subj_id, matrix(mean_betas[condi,], ncol = length(rois)))
      cond_names <- as.vector(outer(rois, condo, paste_fun))
      colnames(cond_frame) <- c('Subject', cond_names)

      if (merge_in) {
        results <- merge(results, cond_frame, by = 'Subject', all = TRUE)
      } else {
        #
        row_indices <- match(subj_id, results[['Subject']])
        results[cond_names][row_indices, ] <- cond_frame[cond_names]
      }
    }
    merge_in <- FALSE
  }

  # Convert zeros to NAs
  results[results == 0] <- NA
  return(results)
}
# ------------------------------------------------------------------------------ #


#-----------------------------------------------------------------------#
do_for_subj_dirs <- function(loc, subj_ids, fn) {

  ###################################################################
  #        This does nothing right now. Deal with it later...       #
  ###################################################################

  # Different base directories for cluster or local use
  if (loc == 'local') { base_dir <- '/home/dan/projects/imagen/'
  } else {              base_dir <- '/users/dscott3/projects/imagen/'}

  afni_dir <- paste(base_dir, '/afni/', sep = '')
  results  <- list()

  # Cycle through subjects directories, doing function "fn"
  for (subj_id in subj_ids[[1]]) {
    subj_id_str <- formatC(subj_id, width = 12, format = 'd', flag = '0')

    # Setup names of stuff
    subj_dir   <- paste(afni_dir, subj_id_str, '_DATA/', sep = '')
    subj_1Ds   <- paste(subj_dir, '/1Ds/', sep = '')

    results <- append(results, list())
  }
}
#-----------------------------------------------------------------------#


