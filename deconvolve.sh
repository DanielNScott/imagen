#!/bin/bash

mask=${1}
dir1D=${2}

# ANCHOR will be replaced depending on the number of
# stop success conditions each subject has.
3dDeconvolve -input prct_change+tlrc.HEAD           \
             -force_TR 2.2                          \
             -TR_times 2.2                          \
             -polort 7                              \
             -censor censor__combined_2.1D          \
             -x1D_uncensored X.nocensor.xmat.1D     \
             -mask ${mask}                          \
             -local_times                           \
             -allzero_OK                            \
             -GOFORIT 3                             \
             -num_stimts 13                         \
             -stim_times 1 ${dir1D}go_success.1D              'SPMG3'  -stim_label 1 go_success    \
             -stim_times 2 ${dir1D}stop_success.1D            'SPMG3'  -stim_label 2 stop_success  \
             -stim_times 3 ${dir1D}stop_failure.1D            'SPMG3'  -stim_label 3 stop_failure  \
             -stim_times 4 ${dir1D}go_failure.1D              'SPMG3'  -stim_label 4 go_failure    \
             -stim_times 5 ${dir1D}go_too_late.1D             'SPMG3'  -stim_label 5 go_too_late   \
             -stim_times 6 ${dir1D}go_wrong_key_response.1D   'SPMG3'  -stim_label 6 go_wrong_key  \
             -stim_times 7 ${dir1D}stop_too_early_response.1D 'SPMG3'  -stim_label 7 go_too_early  \
             -stim_file 8 ${dir1D}motion_demean.1D'[0]'  -stim_base 8  -stim_label 8  roll  \
             -stim_file 9 ${dir1D}motion_demean.1D'[1]'  -stim_base 9  -stim_label 9  pitch \
             -stim_file 10 ${dir1D}motion_demean.1D'[2]' -stim_base 10 -stim_label 10 yaw   \
             -stim_file 11 ${dir1D}motion_demean.1D'[3]' -stim_base 11 -stim_label 11 dS    \
             -stim_file 12 ${dir1D}motion_demean.1D'[4]' -stim_base 12 -stim_label 12 dL    \
             -stim_file 13 ${dir1D}motion_demean.1D'[5]' -stim_base 13 -stim_label 13 dP    \
             -fitts fitts.                          \
             -errts errts.                          \
             -x1D X.xmat.1D -xjpeg X.jpg            \
             -bucket stats.                         \
             -x1D_stop




