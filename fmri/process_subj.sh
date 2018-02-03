#!/bin/tcsh -xef

#----------------------------- Settings ------------------------------------#
set apt = 'BL'                 # Appointment
set IMTag = ''                 # '' or '_IM' for by-trial regressors
set IMHRF = 'SPMG1(0)'         # HRF for any by-trial regressors
set HRF = 'SPMG1(0)'           # HRF for pooled trial regressors
set statsFlags = ''            # Stats from deconvolves: -fout -tout -rout
                               # NOTE: needs to be off for extract.sh to work
                               #       properly
set do_setup = 1
set do_preproc = 1
set do_regress = 1
set do_postproc = 1

# Since TCSH is a headache to do math in, these need to be set consistently...
# ... but don't work since afni uses literals to escape regexes (stupid!!!)
#
#set start_tr = 4     # E.g. 4 if skipping 0,1,2,3
#set trs_left = 440   # total num TRs - start_tr

# This one works:
set secs_rmed = 4.4   # 2.2*(start_tr)
#---------------------------------------------------------------------------#

set dataFile = ${apt}_SST.nii.gz
set output_dir = results

set subj = 's01'
set run = '01'


if ($do_setup == 1 ) then
   # =========================== auto block: setup ============================
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

   # verify that the results directory does not yet exist
   if ( -d $output_dir ) then
       echo output dir "results" already exists
       exit
   endif

   mkdir $output_dir
   mkdir $output_dir/plots
   mkdir $output_dir/regress
   mkdir $output_dir/preproc
   mkdir $output_dir/mask_avgs
   mkdir $output_dir/ideals

   # set list of runs
   set runs = (`count -digits 2 1 1`)


   # ============================ auto block: tcat ============================
   # apply 3dTcat to copy input dsets to results dir, while
   # removing the first 2 TRs
   3dTcat -prefix $output_dir/pb00.$subj.r01.tcat BL_SST.nii.gz'[2..400]'

   # and make note of repetitions (TRs) per run
   set tr_counts = ( 442 )

   # ======================= Manual block: fix timing  ========================
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
endif

if ($do_preproc == 1) then
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

   # ========================== auto block: outcount ==========================
   # data check: compute outlier fraction for each volume
   touch out.pre_ss_warn.txt
   foreach run ( $runs )
       3dToutcount -automask -fraction -polort 7 -legendre                     \
                   pb00.$subj.r$run.tcat+tlrc > outcount.r$run.1D

       # censor outlier TRs per run, ignoring the first 0 TRs
       # - censor when more than 0.05 of automask voxels are outliers
       # - step() defines which TRs to remove via censoring
       1deval -a outcount.r$run.1D -expr "1-step(a-0.05)" > rm.out.cen.r$run.1D

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

   # ================== mask out very high variance voxels ===================
   # old: 'min(200, a/b*100)*step(a)*step(b)'
   # new: z-score
   #foreach run ( $runs )
   #    3dTstat -prefix rm.mean_r$run pb00.$subj.r$run.tcat+tlrc
   #    3dcalc -a pb00.$subj.r$run.tcat+tlrc -b rm.mean_r$run+tlrc \
   #           -expr '(a-b)/stdev(b))'                             \
   #           -prefix pb04.$subj.r$run.scale
   #end

   # ================================ censor =================================
   # Would compute motion params here, but we already have them.

   # create censor file motion_${subj}_censor.1D, for censoring motion
   1d_tool.py -infile '../motion_demean.1D{2..400}' -set_nruns 1 \
              -show_censor_count -censor_prev_TR -censor_motion 0.2 motion_${subj}

   # combine multiple censor files
   1deval -a motion_${subj}_censor.1D -b outcount_${subj}_censor.1D         \
          -expr "a*b" > censor_${subj}_combined_2.1D

   # note TRs that were not censored
   set ktrs = `1d_tool.py -infile censor_${subj}_combined_2.1D              \
                          -show_trs_uncensored encoded`

   # ======================== detrend, bandpass  ==============================
   # 3dTproject dumps time points even when told not to!
   #3dTproject -polort 7 -prefix detrended \ #-cenmode NTRP -censor censor_${subj}_combined_2.1D \
   #  -passband 0.01 0.1 -TR 2.2 -input pb00.$subj.r$run.tcat+tlrc

   3dDespike -prefix despike pb00.$subj.r$run.tcat+tlrc
   3dBandpass -dt 2.2 -nodetrend -band 0.01 9999 -prefix bandpass despike+tlrc
   3dDetrend -polort 7 -prefix detrend bandpass+tlrc

   # ================================= scale ==================================
   # old: 'min(200, a/b*100)*step(a)*step(b)'
   # new: z-score
   foreach run ( $runs )
       3dTstat -prefix mean detrend+tlrc
       3dTstat -prefix stdev -stdev detrend+tlrc
       3dcalc -a detrend+tlrc -b mean+tlrc -c stdev+tlrc -expr '(a-b)/c' -prefix zscored
   end
endif

if ($do_regress == 1) then
   # ================================ regress =================================
   3dDeconvolve \
       -force_TR 2.2 \
       -mask ${mask}                          \
       -input zscored+tlrc \  #pb04.$subj.r*.scale+tlrc.HEAD                                 \
       -censor censor_${subj}_combined_2.1D                                 \
       -TR_times 2.2                                                        \
       -GOFORIT 20 \
       -polort 0                                                            \
       -num_stimts 13                                                       \
       -stim_times 1 ${dir1D}go_success.1D              "${HRF}" -stim_label 1 go_success     \
       -stim_times${IMTag} 2 ${dir1D}stop_success.1D  "${IMHRF}" -stim_label 2 stop_success   \
       -stim_times${IMTag} 3 ${dir1D}stop_failure.1D  "${IMHRF}" -stim_label 3 stop_failure   \
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
       ${statsFlags} -x1D X.xmat.1D -xjpeg X.jpg                              \
       -x1D_uncensored X.nocensor.xmat.1D                                   \
       -fitts fitts.$subj                                                   \
       -errts errts.${subj}                                                 \
       -bucket stats.$subj


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
   3dREMLfit -matrix X.xmat.1D -input zscored+tlrc \
             ${statsFlags} -Rbuck stats.s01_REML -Rvar stats.s01_REMLvar \
             -mask ${mask} -GOFORIT 20 -Rfitts fitts.s01_REML -Rerrts errts.s01_REML -verb $*
endif

if ($do_postproc == 1) then
   # ============================= dump betas ===============================
   # Remove any old stats dump file.
   rm stats_masked.out L_* R_*

   # Set ROIs
   rois=(R_IFG R_preSMA R_Caudate_AAL R_GPe R_GPi R_STN R_Thalamus_AAL R_NAcc R_AAL_ACC L_AAL_ACC)

   # Cycle through ROIs and extract betas
   foreach mask ${rois[@]}
      echo "#--------------------------------------------------------------------#"
      echo "  Dumping stats for voxels from ${mask} mask."
      echo "#--------------------------------------------------------------------#"
      rm zyxt+tlrc.*
      rm fstat_mask_${mask}*

      # Compute mask with full F-statistic q < 0.05
      3dFDR -input stats.s01_REML+tlrc.HEAD'[0]' -qval  -mask_file ROI_masks/${mask}.nii.gz
      3dcalc -a zyxt+tlrc.BRIK'[0]' -expr 'step(0.1-a)' -prefix fstat_mask_${mask}

      # Get average stats over the mask, but only for beta values (sub-briks 1,3,...,17)
      3dmaskave -quiet -mask fstat_mask_${mask}+tlrc.BRIK.gz stats.s01_REML+tlrc.BRIK.gz'[1..$(2)]' \
      | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked.out

      # Insert a newline after each ROI
      printf '\n' >> stats_masked.out

      # Create roi response and fit
      3dmaskave -quiet -mask fstat_mask_${mask}+tlrc.BRIK.gz zscored+tlrc.BRIK > ${mask}_ave.1D
      3dmaskave -quiet -mask fstat_mask_${mask}+tlrc.BRIK.gz fitts.s01_REML+tlrc.BRIK.gz > ${mask}_fitts.1D
   end

   # Create copies
   mv stats_masked.out stats_masked.out.bckp
   cp stats_masked.out.bckp ../stats_masked.out

   # Old snippets:

   #
   # Remove the pesky F-statistic in the middle of the stop-condition block...
   # Applies when using SPMGX with X > 1
   #
   #3dinfo -verb stats.subj_REML+tlrc.HEAD > 3dinfo.txt
   #badCol=`grep 'stop_success_Fstat' 3dinfo.txt | cut -d ' ' -f 6 | tr '#' ' '`
   #let "badCol=${badCol}+1"
   #cut -d ' ' -f ${badCol} stats_masked.out.bckp --complement > stats_masked.out
   #cp stats_masked.out ../

   #
   # Cluster tool for selecting functional ROIs
   #
   #3dclust -savemask clust_mask -1clip 0.95 -NN3 1 FFDRq.s01_REML+tlrc.BRIK.gz'[0]'
   #3dmask_tool -inputs cluster_mask+tlrc.HEAD ROI_masks/${mask}.nii.gz -intersection -prefix ${mask}_cluster_mask

   # ======================= auto block: post proc ==========================
   #
   # if 3dREMLfit fails, terminate the script
   if ( $status != 0 ) then
       echo '---------------------------------------'
       echo '** 3dREMLfit error, failing...'
       exit
   endif

   # create an all_runs dataset to match the fitts, errts, etc.
   3dTcat -prefix all_runs.$subj zscored+tlrc.HEAD #pb04.$subj.r*.scale+tlrc.HEAD

   # --------------------------------------------------
   # create a temporal signal to noise ratio dataset
   #    signal: if 'scale' block, mean should be 100
   #    noise : compute standard deviation of errts
   3dTstat -mean -prefix rm.signal.all all_runs.$subj+tlrc"[$ktrs]"
   3dTstat -stdev -prefix rm.noise.all errts.${subj}_REML+tlrc"[$ktrs]"
   3dcalc -a rm.signal.all+tlrc                                             \
          -b rm.noise.all+tlrc                                              \
          -c full_mask.$subj+tlrc                                           \
          -expr 'c*a/b' -prefix TSNR.$subj

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

   # more and better plots...
   matlab -nosplash -nodisplay -r "addpath(genpath('../../../fmri')); plot_diagnostics; exit"

   echo "execution finished: `date`"
endif

