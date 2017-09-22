

write_stim_1D <- function(task_data, subj_id_str) {

  # Create a directory to put files in
  dir1D <- paste(subj_id_str, '_1Ds', sep = '')
  dir.create(dir1D)

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
    fname <- paste('./', dir1D, '/', tolower(stim), hand_abrv, '.1D', sep = '')

    if (length(times) != 0) {
      # Write table doesn't work here becauase it inserts row indexing...
      write(times, file = fname, ncolumns = 100000)

    } else {
      # This is the convention in afni...
      write(x = '*', file = fname)
    }
  }
}
