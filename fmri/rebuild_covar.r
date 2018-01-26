rebuild_covar <- function(fmri_cov, subj_num, roi_nums){
  # Covariance matrices
  roi_nums <- c(2,3,5,6)
  subj_id  <- fmri_cov['Subject'][subj_num,]
  flat_cov <- fmri_cov[subj_num,2:101]

  subj_cov <- matrix(flat_cov, ncol = 10)[roi_nums, roi_nums]
  subj_cov <- matrix(unlist(subj_cov), ncol = length(roi_nums), byrow = TRUE)
  colnames(subj_cov) <- c('rPreSMA',  'rCaudate', 'rGPi', 'rSTN')

  return(subj_cov)
}