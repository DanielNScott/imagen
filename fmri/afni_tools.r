## Contains:
# - write_stim_1D
# - setup_afni
# - read_stats_dump
# - read_all_stats

#-----------------------------------------------------------------------#
write_stim_1D <- function(task_data, subj_dir, visit, ag = FALSE) {
  # Writes a set of AFNI .1D stimulus files for the SST.
  #
  # Args:
  #   task_data: A data frame with trial condition timing and outcome info.
  #   subj_dir: The directory in which to create and '/1Ds'.
  #
  # Returns:
  #   Nothing

  # Possible conditions are events x hands.
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
  # ToDo:
  # - stop_too_early_response into go_success
  # - what is go_failure? as opposed to go_too_late?

  # Determine which conditions get files.
  conditions <- c(outer(events[!collapse], hands, paste), events[collapse])

  # Establish which conditions get seperation of trials by hand
  n_handed_conds <- 2*length(events[!collapse])
  n_conds        <- length(conditions)
  handed_trials  <- c(rep_len(TRUE, n_handed_conds), rep_len(FALSE, n_conds - n_handed_conds))

  # stim_text gets written into the afni processing command
  # as the set of lines describing the stimuli
  stim_text <- ''
  global_instance_num <- 1

  # Write a file for every condition
  for (cond in conditions) {

    # Figure out if we have a handedness split for this event
    split      <- strsplit(cond, split = ' ')[[1]]
    stim       <- split[1]
    task_stims <- task_data['Stimulus.Presented']

    # Determine which hand is being used
    if (length(split) == 2) {
      hand      <- split[2]
      hand_msk  <- task_stims == hand
      hand_abrv <- ifelse(hand == 'LeftArrow', '_left', '_right')

    } else {
      hand_msk  <- task_stims == hands[1] || task_stims == hands[2]
      hand_abrv <- ''

      # Hand mask should always be TRUE here, i.e. the record indicates
      # either the left or right hand is in use.
      assertthat::assert_that(hand_msk)
    }

    # Either way, write out standard (i.e. non-collapsed) .1D file
    # There is some duplication of work (and overwriting) here, but this makes it
    # easier to swap between analysis conditions.

    # Determine the set of times this event ocurred
    event_mask  <- (task_data['Response.Outcome'] == stim)
    event_times <- t(task_data[event_mask, 'Trial.Start.Time..Onset.']) / 1000
    event_fname <- paste(subj_dir, '/1Ds_', visit, '/', tolower(stim), '.1D', sep = '')
    #print(times[length(times)])

    if (length(event_times) != 0) {
      # Write table doesn't work here becauase it inserts row indexing...
      write(event_times, file = event_fname, ncolumns = 100000)

    } else {
      # This is the convention in afni...
      write(x = '*', file = event_fname)
    }

    # Get list of condition times and a filename
    cond_mask  <- (task_data['Response.Outcome'] == stim) & (hand_msk)
    cond_times <- t(task_data[cond_mask, 'Trial.Start.Time..Onset.']) / 1000
    cond_fname <- paste(subj_dir, '/1Ds_', visit, '/', tolower(stim), hand_abrv, '.1D', sep = '')

    if (length(cond_times) != 0) {
      # Write table doesn't work here becauase it inserts row indexing...
      write(cond_times, file = cond_fname, ncolumns = 100000)

      cond_index <- grep(cond, conditions)
      if (handed_trials[cond_index]) {

        local_stim_num <- 0
        rep_set <- 1:length(cond_times)
        for (instance in rep_set) {
          afni_stim_name  <- paste(tolower(stim), hand_abrv, sep = '')

          fname_inst <- paste(subj_dir, '/1Ds_', visit, '/', afni_stim_name, '_', toString(instance), '.1D', sep = '')
          fname_rest <- paste(subj_dir, '/1Ds_', visit, '/', afni_stim_name, '_', toString(instance), '_rest.1D', sep = '')

          this_time   <- cond_times[instance]
          other_times <- cond_times[setdiff(rep_set, instance)]
          #write(this_time  , file = fname_inst, ncolumns = 1)
          #write(other_times, file = fname_rest, ncolumns = 100000)

          if (cond == 'STOP_SUCCESS LeftArrow' || cond == 'STOP_FAILURE LeftArrow') {
            local_stim_num <- local_stim_num + 1
            stim_num   <- toString(13 + global_instance_num)
            stim_text  <- paste(stim_text,  ' -stim_times ', stim_num, ' ${dir1D}', afni_stim_name, '_',
                                local_stim_num, '.1D', ' \'SPMG1(0)\' ', '-stim_label ', stim_num, ' ',
                                afni_stim_name, '_', local_stim_num ,' \\\\\n', sep = '')

            global_instance_num <- global_instance_num + 1
          }
        }
      }

    } else {
      # By AFNI convention, conditions with no instances use files with '*' in them.
      write(x = '*', file = cond_fname)
    }
  }

  # Fix up AG deconvolution file.
  #deconv_file <- paste(dir, '/deconvolve.sh', sep = '')
  #deconv_text <- readChar(deconv_file, file.info(deconv_file)$size)

  #deconv_text <- sub('NSTIM_ANCHOR\n', paste(' -num_stimts', toString(13 + global_instance_num - 1), ' \\\\\n'), deconv_text)
  #deconv_text <- sub('STIMS_ANCHOR\n', stim_text, deconv_text)

  #write(x = deconv_text, file = deconv_file)

}
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
setup_afni <- function(loc, subj_ids, task = 'SST') {
  # Sets up the ./data/afni/ file structure so fits can be created.
  #
  # Args:
  #   loc: 'local' or 'cluster' for determining the base-path
  #   subj_ids: the subjects to set up the analysis for
  #   task: 'SST' or 'MID' - currently a dummy variable.
  #
  # Returns:
  #   Nothing
  #
  # Details:
  #   Each subject needs a folder with .1D files for motion (softlinked),
  #   and task condition onset times, as well as further links to ROI
  #   masks, fMRI data, and shell scripts for running afni and extracting
  #   data from the afni output.

  # Assign paths for behavioral, fMRI, and fMRI movement data
  home_dir <- ifelse (loc == 'local', '/home/dan/', '/users/dscott3/')
  base_dir <- paste(home_dir, 'projects/imagen/', sep='')

  # We'll put afni files under ./data/afni/
  afni_dir   <- paste(base_dir, '/data/afni/', sep = '')

  for (visit in c('BL', 'FU2')) {
    taskdata_dir   <- paste(base_dir, '/data/fmri/', visit, '_SST_task/', sep = '')
    timeseries_dir <- paste(base_dir, '/data/fmri/', visit, '_SST_AFNI/', sep = '')
    movement_dir   <- paste(base_dir, '/data/fmri/', visit, '_SST_move/', sep = '')

    # Set filename prefixes and suffixes
    task_pfx <- 'ss_'
    fmri_sfx <- paste('_', visit, '_SST.nii.gz', sep = '')
    move_pfx <- 'EPI_stop_'

    # Begin writing shell script (which gets filled out below)
    # to submit everything to SLURM
    submitter_fname <- paste(afni_dir, 'submitter_', visit, '.sh', sep='')
    write('#!/bin/bash \n', file = submitter_fname)

    # Setup specific subjects' directories
    for (subj_id in subj_ids) {
      # ( subj_ids is a data frame so needs conversion using [[ ]] )

      # Convert subject id to appropriate string
      subj_id_str <- formatC(subj_id, width = 12, format = 'd', flag = '0')

      # Filename beurocracy...
      fmri_fname <- paste(subj_id_str, fmri_sfx, sep = '')
      task_fname <- paste(task_pfx, subj_id_str, '.csv', sep = '')
      move_fname <- paste(move_pfx, subj_id_str, '.txt', sep = '')

      fmri_full_name <- paste(timeseries_dir, fmri_fname, sep = '')
      task_full_name <- paste(taskdata_dir  , task_fname, sep = '')
      move_full_name <- paste(movement_dir  , move_fname, sep = '')

      # Read the SST series
      if (file.exists(task_full_name)) {
        # Only want to read 2 columns, trial time and outcome
        # ToDo: Figure out why there are 4 here.
        # ToDo: Collapse stop-too-early into go condition
        cols <- c('NULL','NULL','numeric','NULL','NULL',NA    ,'NULL','NULL',
                  'NULL','NULL','NULL'   , NA   ,'NULL','NULL','NULL','numeric')

        task_data <- read.csv(file = task_full_name, header = TRUE, sep = '\t',
         colClasses = cols, skip = 1)

        print(paste('Successful read of', task_fname))
      } else {
        print(paste('File', task_fname, 'does not exist. Aborting setup for this subject.'))
        next
      }

      # Create subjedt directory and 1D directory
      subj_dir   <- paste(afni_dir, subj_id_str, '_DATA', sep = '')
      subj_1Ds   <- paste(subj_dir, '/1Ds_', visit, sep = '')
      dir.create(subj_dir)
      dir.create(subj_1Ds)

      # Soft link commands
      proc_fname <- 'process_subj.sh'
      #wipe_fname <- 'wipe_results.sh'
      link_cmd1  <- paste('ln -s ', base_dir, '/fmri/', proc_fname, ' ', subj_dir, '/', proc_fname, sep = '')
      link_cmd2  <- paste('ln -s ', base_dir, '/data/ROI_masks', ' ', subj_dir, '/', 'ROI_masks', sep = '')
      link_cmd3  <- paste('cp ', move_full_name, ' ', subj_dir, '/motion_demean_', visit, '.1D', sep = '')
      link_cmd4  <- paste('ln -s ', fmri_full_name, ' ', subj_dir, '/', visit, '_SST.nii.gz', sep = '')
      #link_cmd5  <- paste('cp    ', afni_dir, 'deconvolve.sh', ' ', subj_dir, '/deconvolve.sh', sep = '')
      #link_cmd6  <- paste('ln -s ', afni_dir, 'extract.sh', ' ', subj_dir, '/', 'extract.sh', sep = '')

      # Create links
      if (visit == 'BL') {
        system(link_cmd1)
        system(link_cmd2)
      }
      system(link_cmd3)
      system(link_cmd4)
      #system(link_cmd5)
      #system(link_cmd6)

      # Write the behavioral 1D files
      write_stim_1D(task_data, subj_dir, visit, ag = TRUE)

      # Create the set of commands associated with entering this directory
      # and submitting processing for this subject to SLURM
      jobname  <- paste(subj_id_str, 'glm', sep = '_')
      cdcmd1   <- paste('cd ', subj_id_str, '_DATA;', sep = '')
      cdcmd2   <- 'cd ../; sleep 1;'
      srunpfx  <- 'sbatch -t 0:05:00 -n 1 --nodes 1 --cpus-per-task 1 -J'
      srun_cmd <- paste(cdcmd1, srunpfx, jobname, 'process_subj.sh ;', cdcmd2, sep = ' ')

      # Append the commands to the mass-submission file
      write(srun_cmd, file = submitter_fname, append = TRUE)
    }
  }
}
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
read_stats_dump <- function(subj_dir, ag, hrf) {
  # Reads the file of glm regression weights "stats_masked.out" for the
  # subjects' rois, names the statistics, and (possibly) removes some.
  #
  # Args:
  #   subj_dir: (string) the subject fmri processing folder
  #   ag: (logical) produce output for the ancestral graphs processing
  #   hrf: 'SPMG1' or 'SPMG3' for the simple and multi-basis standard HRFs
  #
  # Output:
  #   An array of roi regression weights from the GLM fit.
  #
  # Notes:
  #   Currently assumes statsTags is set to '' in process_subj.sh.
  #   This matters for the indexing in the ancestral graphs processing.

  # Read the results of running AFNI:
  # roi_stats should be a matrix with ROIs as rows and statistics as columns.

  fname     <- paste(subj_dir, 'stats_masked.out', sep = '')
  roi_stats <- read.table(fname, sep = ' ', header = FALSE, blank.lines.skip = FALSE)

  # Two helper functions
  paste_fun <- function(x,y) {paste(x, y, sep = '_')}
  str_outer <- function(x,y) {as.vector(t(outer(x, y, paste_fun)))}

  # Throw away some statistics we don't want
  if (ag) {
    # If using AFNI outlier removal tool on individual-modulation betas
    #censor_file <- paste(subj_dir, 'stats.rm.out.cen.1D', sep = '')
    #censor_list <- read.csv(censor_file)

    # Trim off non-beta values:
    # The first 5 are Full_Fstat and go_success statistics
    # The s 9 from the end and 5 from the front
    # We find the F-stat in the middle on the basis of the trial-data later
    n_stats     <- length(roi_stats)
    remove_inds <- c(1:5, seq(n_stats, n_stats-8, -1))
    keep_inds   <- setdiff(1:n_stats,  c(1:5, seq(n_stats, n_stats-8, -1)))
    roi_stats   <- roi_stats[, keep_inds]
    #censor_list <- censor_list[keep_inds, ]

    colnames(roi_stats) <- paste('stop_signal:', 1:dim(roi_stats)[2])
    return(roi_stats)
  }

  # Below we name all the statistics being read
  #
  # Regressor types
  # For SPMG3(0) regressors have HRF, derivative and dispersion deriv.
  if (hrf == 'SPMG3') {
    types  <- c('hrf', 'hrf_dt', 'hrf_dd')
  } else {
    # Don't need to label it when there's only one component.
    types  <- c('')
  }

  # We assume the file contains betas and T-values for each of these
  #stats  <- c('beta', 'tval')

  # Not anymore! Don't need to label this if there's only one either.
  stats <- c('')

  # Conditions: Go success, stop success, etc.
  conds  <- c('gs', 'ss', 'sf', 'gf', 'gtl', 'gwk', 'ste')

  # The contrasts being computed
  ctrsts <- c('ss-sf', 'ss-go')

  # These two loops build the long list of names
  #stat_names <- c('all_Fval')
  stat_names <- c()
  for (cond in conds) {
    # Setup "cond_type_stat" fields
    block <- str_outer(cond, types)
    block <- str_outer(block, stats)

    # Append the condition level F-value and save
    if (hrf == 'SPMG3') { # No F-val if only one basis fn
      block <- c(block, paste_fun(cond, 'Fval'))
    }
    stat_names <- c(stat_names, block)
  }
  for (ctrst in ctrsts) {
    block <- str_outer(ctrst, stats)
    if (hrf == 'SPMG3') {
      stat_names <- c(stat_names, block, paste_fun(ctrst, 'Fval'))
    }
      stat_names <- c(stat_names, block)
  }

  # Assign the names and return the stats
  colnames(roi_stats) <- stat_names
  roi_stats <- t(roi_stats)

  fn   <- function(x,y) { paste(x, y, sep='_') }
  rois  <- c('rIFG', 'rPreSMA', 'rCaudate', 'rGPe', 'rGPi', 'rSTN', 'rThalamus', 'rNAcc', 'rACC', 'lACC')
  names     <- as.vector(outer(stat_names, rois, fn))
  roi_stats <- as.vector(roi_stats)
  names(roi_stats) <- names

  return(roi_stats)
}
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
read_all_stats <- function(subj_ids, afni_dir, ag = FALSE, hrf) {

  # Regions of interest
  rois <- c('rIFG', 'rPreSMA', 'rCaudate', 'rGPe', 'rGPi',
    'rSTN', 'rThalamus', 'rNAcc', 'rACC', 'lACC')

  # Conditions of interest:
  # Go Success, Stop Success, Stop Failure, contrasts
  cois <- c('gs_hrf_beta', 'ss_hrf_beta', 'sf_hrf_beta', 'ss_go_beta', 'ss_sf_beta')
  cond <- c('go', 'st', 'sr', 'st_go', 'st_sr')

  # Names for each roi, e.g. rIFG_go
  paste_fun   <- function(x,y) {paste(x, y, sep = '_')}
  roi_by_cond <- as.vector(outer(rois, cond, paste_fun))

  # Intialize output frame everything goes in
  results <- data.frame('Subject' = subj_ids)

  # Name the standard fMRI analysis results
  std_fmri_feats <- as.vector(outer(rois, cond, function(x,y) {paste(x,y, sep = '')}))

  # This is opaque...
  merge_in <- TRUE

  data <- c()
  # Process each subject's data.
  for (subj_id in subj_ids) {
    subj_id_str <- formatC(subj_id, width = 12, format = 'd', flag = '0')

    # Get the name of the subject's data folder
    subj_dir <- paste(afni_dir, subj_id_str, '_DATA/', sep = '')
    # Attempt to read the stats dump file
    tryCatch({
      stats <- read_stats_dump(subj_dir, ag, hrf)

      if (ag) {
        # ToDo: I should probably save all the stats not just the covar. matrices...
        #
        stop('AG/covariance processing needs to be fixed before being used.')

        # The stats consist of a regression weight for each trial,
        # which we'll turn into a covariance matrix.
        stats <- t(stats)

        # Standard covariance matrix estimation isn't robust, so
        # we'll perform outlier removal based on Median-Absolute-Deviation
        med     <- apply(stats, 2, median, na.rm = T)
        mad_est <- apply(stats, 2, mad,    na.rm = T)
        mad_sd  <- 1.48*mad_est
        mad_bnd <- rbind(med + 2*mad_sd, med - 2*mad_sd)
        out_fn  <- function(x) {x > mad_bnd[1,] | x < mad_bnd[2,]}
        out_msk <- t(apply(stats, 1, out_fn))

        stats[out_msk] <- NA

        # An alternative is to use something like the LW estimator
        # ToDo: Robust estimator

        # Convert betas to covariance matrices
        colnames(stats) <- rois
        stats <- var(stats, na.rm = T)
        dim(stats) <- NULL
        stats <- matrix(stats, ncol = length(rois)^2)

        # Sort out proper names
        cond_names <- outer(rois, 1:length(rois), paste)
        colnames(stats) <- cond_names

        # Agglomerate stuff
        cond_frame <- data.frame('Subject' = subj_id, stats)
        if (merge_in) {
          results <- merge(results, cond_frame, by = 'Subject')
        } else {
          results <- rbind(results, cond_frame)
        }

      } else {

        new_data <- data.frame('Subject' = subj_id, t(stats))
        data <- rbind(data, new_data)

        # Collect subject data
        #mean_betas <- stats[cois,]
        #colnames(mean_betas) <- rois
        #for (cnum in 1:length(cois)) {
        #   condi <- cois[cnum]
        #   condo <- cond[cnum]

        #   cond_frame <- data.frame('Subject' = subj_id, matrix(mean_betas[condi,], ncol = length(rois)))
        #   cond_names <- as.vector(outer(rois, condo, paste_fun))
        #   colnames(cond_frame) <- c('Subject', cond_names)

        #   if (merge_in) {
        #     results <- merge(results, cond_frame, by = 'Subject', all = TRUE)
        #   } else {
           #
        #     row_indices <- match(subj_id, results[['Subject']])
        #     results[cond_names][row_indices, ] <- cond_frame[cond_names]
        #   }
        #}
      }
      merge_in <- FALSE

    }, error = function(e){
      print(paste('Unable to read file stats for', subj_id_str))
    })
  }

  # Convert zeros to NAs
  data[data == 0] <- NA

  return(data)
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


