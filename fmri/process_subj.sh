#!/bin/tcsh -xef

#
set IMTag = '_IM'              # '' or '_IM' for individual modulation
set IMHRF = 'SPMG1(0)' #'SPMG1(0)'         # 'SPMG3' or 'SPMG1(0)'
set HRF = 'SPMG1(0)' #'SPMG1(0)'         # 'SPMG3' or 'SPMG1(0)'
set statsTags = '-fout -tout'             # '-fout -tout'
set dataFile = BL_SST.nii.gz

# Since TCSH is a headache to do math in, these need to be set consistently
# These don't work since substitution syntax isn't working out for me atm...
#set start_tr = 4      # E.g. 4 if skipping 0,1,2,3
#set trs_left = 440    # total num TRs - start_tr
set secs_rmed = 4.4   # 2.2*(start_tr)

# execute via :
#   tcsh -xef proc.s01 |& tee output.proc.s01

# =========================== auto block: setup ============================
# script setup

# take note of the AFNI version
afni -ver

# check that the current AFNI version is recent enough
afni_history -check_date 23 Sep 2016
if ( $status ) then
    echo "** this script requires newer AFNI binaries (than 23 Sep 2016)"
    echo "   (consider: @update.afni.binaries -defaults)"
    exit
endif

# the user may specify a single subject to run with
if ( $#argv > 0 ) then
    set subj = $argv[1]
else
    set subj = s01
endif

# assign output directory name
set output_dir = results

# verify that the results directory does not yet exist
if ( -d $output_dir ) then
    echo output dir "results" already exists
    exit
endif

mkdir $output_dir

# set list of runs
set runs = (`count -digits 2 1 1`)


# ============================ auto block: tcat ============================
# apply 3dTcat to copy input dsets to results dir, while
# removing the first 2 TRs
3dTcat -prefix $output_dir/pb00.$subj.r01.tcat \
    000020320615_BL_SST.nii.gz'[2..$]'

# and make note of repetitions (TRs) per run
set tr_counts = ( 442 )

# -------------------------------------------------------

# enter the results directory (can begin processing data)
cd $output_dir

set mask = ROI_masks/All_ROIs.nii.gz
set dir1D = 1Ds/

# Link stuff
cp -r ../$dir1D ./
ln -s ../ROI_masks
ln -s ../$dataFile
ln -s ../extract.sh

foreach f ( ./$dir1D/*.1D )
    timing_tool.py -timing $f -add_offset -$secs_rmed -write_timing $f
end

# ========================== auto block: outcount ==========================
# data check: compute outlier fraction for each volume
touch out.pre_ss_warn.txt
foreach run ( $runs )
    3dToutcount -automask -fraction -polort 7 -legendre                     \
                pb00.$subj.r$run.tcat+tlrc > outcount.r$run.1D

    # censor outlier TRs per run, ignoring the first 0 TRs
    # - censor when more than 0.07 of automask voxels are outliers
    # - step() defines which TRs to remove via censoring
    1deval -a outcount.r$run.1D -expr "1-step(a-0.07)" > rm.out.cen.r$run.1D

    # outliers at TR 0 might suggest pre-steady state TRs
    if ( `1deval -a outcount.r$run.1D"{0}" -expr "step(a-0.4)"` ) then
        echo "** TR #0 outliers: possible pre-steady state TRs in run $run" \
            >> out.pre_ss_warn.txt
    endif
end

# catenate outlier counts into a single time series
cat outcount.r*.1D > outcount_rall.1D

# catenate outlier censor files into a single time series
cat rm.out.cen.r*.1D > outcount_${subj}_censor.1D

# Plot the outlier counts
1dplot -jpeg outliers -one '1D: 442@0.07' outcount_rall.1D

# ================================= tshift =================================
# time shift data so all slice timing is the same
#foreach run ( $runs )
#    3dTshift -tzero 0 -quintic -prefix pb01.$subj.r$run.tshift \
#             pb00.$subj.r$run.tcat+tlrc
#end

# --------------------------------
# extract volreg registration base
#3dbucket -prefix vr_base pb01.$subj.r01.tshift+tlrc"[2]"

# ================================= volreg =================================
# align each dset to base volume
#foreach run ( $runs )
#    # register each volume to the base
#    3dvolreg -verbose -zpad 1 -base vr_base+tlrc                    \
#             -1Dfile dfile.r$run.1D -prefix pb02.$subj.r$run.volreg \
#             -cubic                                                 \
#             pb01.$subj.r$run.tshift+tlrc
#end

# make a single file of registration params
#cat dfile.r*.1D > dfile_rall.1D

# ================================== blur ==================================
# blur each volume of each run
#foreach run ( $runs )
#    3dmerge -1blur_fwhm 4.0 -doall -prefix pb03.$subj.r$run.blur \
#            pb02.$subj.r$run.volreg+tlrc
#end

# ================================== mask ==================================
# create 'full_mask' dataset (union mask)
foreach run ( $runs )
    3dAutomask -dilate 1 -prefix rm.mask_r$run pb00.$subj.r$run.tcat+tlrc
end

# create union of inputs, output type is byte
3dmask_tool -inputs rm.mask_r*+tlrc.HEAD -union -prefix full_mask.$subj

# ================================= scale ==================================
# scale each voxel time series to have a mean of 100
# (be sure no negatives creep in)
# (subject to a range of [0,200])
foreach run ( $runs )
    3dTstat -prefix rm.mean_r$run pb00.$subj.r$run.tcat+tlrc
    3dcalc -a pb00.$subj.r$run.tcat+tlrc -b rm.mean_r$run+tlrc \
           -expr 'min(200, a/b*100)*step(a)*step(b)'           \
           -prefix pb04.$subj.r$run.scale
end

# ================================ regress =================================

# compute de-meaned motion parameters (for use in regression)
#1d_tool.py -infile dfile_rall.1D -set_nruns 1                            \
#           -demean -write motion_demean.1D

# compute motion parameter derivatives (just to have)
#1d_tool.py -infile dfile_rall.1D -set_nruns 1                            \
#           -derivative -demean -write motion_deriv.1D

# Plot motion parameters
1dplot -jpeg motion -volreg ../motion_demean.1D

# Plot motion euclidean norm
#1dplot motion_FT_enorm.1D

# create censor file motion_${subj}_censor.1D, for censoring motion
1d_tool.py -infile '../motion_demean.1D{2..$}' -set_nruns 1 \
           -show_censor_count -censor_prev_TR -censor_motion 0.3 motion_${subj}

# combine multiple censor files
1deval -a motion_${subj}_censor.1D -b outcount_${subj}_censor.1D         \
       -expr "a*b" > censor_${subj}_combined_2.1D

# note TRs that were not censored
set ktrs = `1d_tool.py -infile censor_${subj}_combined_2.1D              \
                       -show_trs_uncensored encoded`

# ------------------------------
# run the regression analysis
3dDeconvolve \
    -force_TR 2.2 \
    -mask ${mask}                          \
    -input pb04.$subj.r*.scale+tlrc.HEAD                                 \
    -censor censor_${subj}_combined_2.1D                                 \
    -TR_times 2.2                                                        \
    -GOFORIT 5 \
    -polort 7                                                            \
    -num_stimts 13                                                       \
    -stim_times 1 ${dir1D}go_success.1D              "${HRF}" -stim_label 1 go_success     \
    -stim_times${IMTag} 2 ${dir1D}stop_success.1D  "${IMHRF}"   -stim_label 2 stop_success   \
    -stim_times${IMTag} 3 ${dir1D}stop_failure.1D  "${IMHRF}"   -stim_label 3 stop_failure   \
    -stim_times 4 ${dir1D}go_failure.1D              "${HRF}" -stim_label 4 go_failure     \
    -stim_times 5 ${dir1D}go_too_late.1D             "${HRF}" -stim_label 5 go_too_late    \
    -stim_times 6 ${dir1D}go_wrong_key_response.1D   "${HRF}" -stim_label 6 go_wrong_key   \
    -stim_times 7 ${dir1D}stop_too_early_response.1D "${HRF}" -stim_label 7 stop_too_early \
    -stim_file  8 ../motion_demean.1D'[0]' -stim_base 8  -stim_label 8  roll  \
    -stim_file  9 ../motion_demean.1D'[1]' -stim_base 9  -stim_label 9  pitch \
    -stim_file 10 ../motion_demean.1D'[2]' -stim_base 10 -stim_label 10 yaw   \
    -stim_file 11 ../motion_demean.1D'[3]' -stim_base 11 -stim_label 11 dS    \
    -stim_file 12 ../motion_demean.1D'[4]' -stim_base 12 -stim_label 12 dL    \
    -stim_file 13 ../motion_demean.1D'[5]' -stim_base 13 -stim_label 13 dP    \
    -gltsym 'SYM: stop_success - stop_failure'                           \
    -gltsym 'SYM: stop_success - go_success'                             \
    -glt_label 1 ss-sf                                                   \
    -glt_label 2 ss-go                                                   \
    -jobs 3                                                              \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                              \
    -x1D_uncensored X.nocensor.xmat.1D                                   \
    -fitts fitts.$subj                                                   \
    -errts errts.${subj}                                                 \
    -bucket stats.$subj

# Plot X.mat
#1dplot X.xmat.1D
1dplot -jpeg xmat_stims  X.xmat.1D'[8..13$]'

# if 3dDeconvolve fails, terminate the script
if ( $status != 0 ) then
    echo '---------------------------------------'
    echo '** 3dDeconvolve error, failing...'
    echo '   (consider the file 3dDeconvolve.err)'
    exit
endif


# display any large pairwise correlations from the X-matrix
1d_tool.py -show_cormat_warnings -infile X.xmat.1D |& tee out.cormat_warn.txt

# -- execute the 3dREMLfit script, written by 3dDeconvolve --
#tcsh -x stats.REML_cmd
3dREMLfit -matrix X.xmat.1D -input pb04.s01.r01.scale+tlrc.HEAD \
 -fout -tout -Rbuck stats.s01_REML -Rvar stats.s01_REMLvar \
 -mask ${mask} -GOFORIT 5 -Rfitts fitts.s01_REML -Rerrts errts.s01_REML -verb $*

bash extract.sh

# if 3dREMLfit fails, terminate the script
if ( $status != 0 ) then
    echo '---------------------------------------'
    echo '** 3dREMLfit error, failing...'
    exit
endif


# create an all_runs dataset to match the fitts, errts, etc.
3dTcat -prefix all_runs.$subj pb04.$subj.r*.scale+tlrc.HEAD

# --------------------------------------------------
# create a temporal signal to noise ratio dataset
#    signal: if 'scale' block, mean should be 100
#    noise : compute standard deviation of errts
3dTstat -mean -prefix rm.signal.all all_runs.$subj+tlrc"[$ktrs]"
3dTstat -stdev -prefix rm.noise.all errts.${subj}_REML+tlrc"[$ktrs]"
3dcalc -a rm.signal.all+tlrc                                             \
       -b rm.noise.all+tlrc                                              \
       -c full_mask.$subj+tlrc                                           \
       -expr 'a/b' -prefix TSNR.$subj

# ---------------------------------------------------
# compute and store GCOR (global correlation average)
# (sum of squares of global mean of unit errts)
3dTnorm -norm2 -prefix rm.errts.unit errts.${subj}_REML+tlrc
3dmaskave -quiet -mask full_mask.$subj+tlrc rm.errts.unit+tlrc           \
          > gmean.errts.unit.1D
3dmaskave -quiet rm.errts.unit+tlrc           \
          > gmean.errts.unit.1D
3dTstat -sos -prefix - gmean.errts.unit.1D\' > out.gcor.1D
echo "-- GCOR = `cat out.gcor.1D`"

# ---------------------------------------------------
# compute correlation volume
# (per voxel: average correlation across masked brain)
# (now just dot product with average unit time series)
3dcalc -a rm.errts.unit+tlrc -b gmean.errts.unit.1D -expr 'a*b' -prefix rm.DP
3dTstat -sum -prefix corr_brain rm.DP+tlrc

# create ideal files for fixed response stim types
1dcat X.nocensor.xmat.1D'[4]' > ideal_go_success.1D
1dcat X.nocensor.xmat.1D'[5]' > ideal_go_failure.1D
1dcat X.nocensor.xmat.1D'[6]' > ideal_stop_success.1D
1dcat X.nocensor.xmat.1D'[7]' > ideal_stop_failure.1D
1dcat X.nocensor.xmat.1D'[8]' > ideal_go_too_late.1D
1dcat X.nocensor.xmat.1D'[9]' > ideal_go_wrong_key.1D
1dcat X.nocensor.xmat.1D'[10]' > ideal_stop_too_early.1D

# --------------------------------------------------------
# compute sum of non-baseline regressors from the X-matrix
# (use 1d_tool.py to get list of regressor colums)
set reg_cols = `1d_tool.py -infile X.nocensor.xmat.1D -show_indices_interest`
3dTstat -sum -prefix sum_ideal.1D X.nocensor.xmat.1D"[$reg_cols]"

# also, create a stimulus-only X-matrix, for easy review
1dcat X.nocensor.xmat.1D"[$reg_cols]" > X.stim.xmat.1D

# ================== auto block: generate review scripts ===================

# generate a review script for the unprocessed EPI data
gen_epi_review.py -script @epi_review.$subj \
    -dsets pb00.$subj.r*.tcat+tlrc.HEAD

# generate scripts to review single subject results
# (try with defaults, but do not allow bad exit status)
gen_ss_review_scripts.py -mot_limit 0.3 -out_limit 0.07 -exit0

# ========================== auto block: finalize ==========================

# remove temporary files
\rm -f rm.*

# if the basic subject review script is here, run it
# (want this to be the last text output)
if ( -e @ss_review_basic ) ./@ss_review_basic |& tee out.ss_review.$subj.txt

# return to parent directory
cd ..

echo "execution finished: `date`"




# ==========================================================================
# script generated by the command:
#
# afni_proc.py -subj_id s01 -dsets 000020320615_BL_SST.nii.gz                \
#     -tcat_remove_first_trs 2 -regress_stim_times 1Ds/stim1_go_success.1D   \
#     1Ds/stim2_go_failure.1D 1Ds/stim3_stop_success.1D                      \
#     1Ds/stim4_stop_failure.1D 1Ds/stim5_go_too_late.1D                     \
#     1Ds/stim6_go_wrong_key.1D 1Ds/stim7_stop_too_early.1D                  \
#     -regress_stim_labels go_success go_failure stop_success stop_failure   \
#     go_too_late go_wrong_key stop_too_early -regress_basis 'SPMG1(0)'      \
#     -regress_opts_3dD -gltsym 'SYM: stop_success - stop_failure' -gltsym   \
#     'SYM: stop_success - go_success' -glt_label 1 ss-sf -glt_label 2 ss-go \
#     -jobs 3 -regress_reml_exec -regress_make_ideal_sum sum_ideal.1D        \
#     -regress_censor_outliers 0.07 -regress_censor_motion 0.3
