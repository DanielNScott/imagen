#!/bin/tcsh -xef

# =========================== auto block: setup ============================
# script setup
# take note of the AFNI version
afni -ver

# check that the current AFNI version is recent enough
afni_history -check_date 28 Oct 2015
if ( $status ) then
    echo "** this script requires newer AFNI binaries (than 28 Oct 2015)"
    echo "   (consider: @update.afni.binaries -defaults)"
    exit
endif

set dir1D = 1Ds/
set mask = ROI_masks/All_ROIs.nii.gz
set dataFile = BL_SST.nii.gz

# assign output directory name
set output_dir = results

# verify that the results directory does not yet exist
if ( -d $output_dir ) then
    echo output dir "results" already exists
    exit
endif


# create results and stimuli directories
mkdir $output_dir
cd $output_dir
ln -s ../$dir1D
ln -s ../ROI_masks
ln -s ../$dataFile
ln -s ../extract.sh

# ============================ auto block: tcat ============================
# apply 3dTcat to copy input dsets to results dir, while
# removing the first 0 TRs
3dTcat -prefix tcat BL_SST.nii.gz'[0..$]'

# -------------------------------------------------------
# enter the results directory (can begin processing data)
#cd $output_dir


# ========================== auto block: outcount ==========================
# data check: compute outlier fraction for each volume
touch out.pre_ss_warn.txt

3dToutcount -automask -fraction -polort 7 -legendre                     \
             tcat+tlrc > outcount.1D

# censor outlier TRs per run, ignoring the first 0 TRs
# - censor when more than 0.1 of automask voxels are outliers
# - step() defines which TRs to remove via censoring
1deval -a outcount.1D -expr "1-step(a-0.1)" > rm.out.cen.1D

# outliers at TR 0 might suggest pre-steady state TRs
if ( `1deval -a outcount.1D"{0}" -expr "step(a-0.4)"` ) then
    echo "** TR #0 outliers: possible pre-steady state TRs in run $run" \
        >> out.pre_ss_warn.txt
endif

# catenate outlier counts into a single time series
cat outcount.1D > outcount_rall.1D

# catenate outlier censor files into a single time series
cat rm.out.cen.1D > outcount__censor.1D

# ================================ despike =================================
# apply 3dDespike to each run
3dDespike -NEW -nomask -prefix despike tcat+tlrc

# ================================ regress =================================
# compute de-meaned motion parameters (for use in regression)
1d_tool.py -infile ${dir1D}motion_demean.1D -set_nruns 1                \
           -demean -write motion_demean.1D

# compute motion parameter derivatives (just to have)
1d_tool.py -infile motion_demean.1D -set_nruns 1                        \
           -derivative -demean -write motion_deriv.1D

# create censor file motion__censor.1D, for censoring motion
1d_tool.py -infile motion_demean.1D -set_nruns 1                        \
    -show_censor_count -censor_prev_TR                                  \
    -censor_motion 0.3 motion_

# combine multiple censor files
1deval -a motion__censor.1D -b outcount__censor.1D        \
       -expr "a*b" > censor__combined_2.1D

# create bandpass regressors (instead of using 3dBandpass, say)
#1dBport -nodata 444 2.2 -band 0.01 99999 -invert -nozero > bandpass_rall.1D

# note TRs that were not censored
set ktrs = `1d_tool.py -infile censor__combined_2.1D             \
                       -show_trs_uncensored encoded`

# Calculate mean and transform to % signal change
3dTstat -mask ${mask} -prefix voxel_mean despike+tlrc.HEAD

3dcalc -a 'despike+tlrc.HEAD' -b 'voxel_mean+tlrc.BRIK'    \
       -expr '100*a/b*ispositive(b-150)' -prefix prct_change

#3dmema - put roi mask instead of whole brain mask, but generally treat the roi as if it's the whole brain
# ------------------------------
# run the regression analysis
3dDeconvolve -input prct_change+tlrc.HEAD           \
             -polort 7                              \
             -censor censor__combined_2.1D          \
             -x1D_uncensored X.nocensor.xmat.1D     \
             -mask ${mask}                          \
             -force_TR 2.2                          \
             -TR_times 2.2                          \
             -local_times                           \
             -allzero_OK                            \
             -GOFORIT 3                             \
             -num_stimts 13                         \
             -stim_times 1 ${dir1D}go_success.1D              'SPMG3'   -stim_label 1 go_success    \
             -stim_times 2 ${dir1D}stop_success.1D            'SPMG3'   -stim_label 2 stop_success  \
             -stim_times 3 ${dir1D}stop_failure.1D            'SPMG3'   -stim_label 3 stop_failure  \
             -stim_times 4 ${dir1D}go_failure.1D              'SPMG3'   -stim_label 4 go_failure    \
             -stim_times 5 ${dir1D}go_too_late.1D             'SPMG3'   -stim_label 5 go_too_late   \
             -stim_times 6 ${dir1D}go_wrong_key_response.1D   'SPMG3'   -stim_label 6 go_wrong_key  \
             -stim_times 7 ${dir1D}stop_too_early_response.1D 'SPMG3'   -stim_label 7 go_too_early  \
             -stim_file  8 ${dir1D}motion_demean.1D'[0]' -stim_base 8  -stim_label 8  roll  \
             -stim_file  9 ${dir1D}motion_demean.1D'[1]' -stim_base 9  -stim_label 9  pitch \
             -stim_file 10 ${dir1D}motion_demean.1D'[2]' -stim_base 10 -stim_label 10 yaw   \
             -stim_file 11 ${dir1D}motion_demean.1D'[3]' -stim_base 11 -stim_label 11 dS    \
             -stim_file 12 ${dir1D}motion_demean.1D'[4]' -stim_base 12 -stim_label 12 dL    \
             -stim_file 13 ${dir1D}motion_demean.1D'[5]' -stim_base 13 -stim_label 13 dP    \
             -num_glt 2                                 \
             -glt_label 1 ss_go_diff                    \
             -gltsym 'SYM: +stop_success -go_success'   \
             -glt_label 2 ss_sf_diff                    \
             -gltsym 'SYM: +stop_success -stop_failure' \
             -fitts fitts.                         \
             -errts errts.                       \
             -fout -tout -x1D X.xmat.1D -xjpeg X.jpg    \
             -tout -fout                                \
             -bucket stats.                        \
             -x1D_stop

# Invoke REML fit
3dREMLfit -matrix X.xmat.1D                      \
          -GOFORIT -input prct_change+tlrc.HEAD  \
          -mask ${mask}        \
          -fout -tout -Rbuck stats.subj_REML     \
          -Rvar stats.subj_REMLvar               \
          -Rfitts fitts.subj_REML                \
          -Rerrts errts.subj_REML -verb $*

sh extract.sh

exit
# -- use 3dTproject to project out regression matrix --
3dTproject -polort 0 -input prct_change+tlrc.HEAD \
           -ort X.xmat.1D -prefix errts..tproject

# if 3dDeconvolve fails, terminate the script
if ( $status != 0 ) then
    echo '---------------------------------------'
    echo '** 3dDeconvolve error, failing...'
    echo '   (consider the file 3dDeconvolve.err)'
    exit
endif

# display any large pairwise correlations from the X-matrix
1d_tool.py -show_cormat_warnings -infile X.xmat.1D |& tee out.cormat_warn.txt

# create an all_runs dataset to match the fitts, errts, etc.
3dTcat -prefix all_runs. prct_change+tlrc.HEAD

# --------------------------------------------------
# create a temporal signal to noise ratio dataset
#    signal: if 'scale' block, mean should be 100
#    noise : compute standard deviation of errts
3dTstat -mean -prefix rm.signal.all all_runs.+tlrc
3dTstat -stdev -prefix rm.noise.all errts..tproject+tlrc
3dcalc -a rm.signal.all+tlrc                                \
       -b rm.noise.all+tlrc                                 \
       -expr 'a/b' -prefix TSNR.

# --------------------------------------------------------
# compute sum of non-baseline regressors from the X-matrix
# (use 1d_tool.py to get list of regressor colums)
set reg_cols = `1d_tool.py -infile X.xmat.1D -show_indices_interest`
3dTstat -sum -prefix sum_ideal.1D X.xmat.1D"[$reg_cols]"

# also, create a stimulus-only X-matrix, for easy review
1dcat X.xmat.1D"[$reg_cols]" > X.stim.xmat.1D

# ================== auto block: generate review scripts ===================

# generate a review script for the unprocessed EPI data
gen_epi_review.py -script @epi_review. \
    -dsets tcat+tlrc.HEAD

# generate scripts to review single subject results
# (try with defaults, but do not allow bad exit status)
gen_ss_review_scripts.py                    \
    -errts_dset errts..tproject+tlrc.HEAD -exit0

# ========================== auto block: finalize ==========================

# remove temporary files
\rm -f rm.*

# if the basic subject review script is here, run it
# (want this to be the last text output)
if ( -e @ss_review_basic ) ./@ss_review_basic |& tee out.ss_review..txt

# return to parent directory
#cd ..


#echo "execution finished: `date`"

