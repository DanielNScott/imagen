#!/bin/bash

#----------------------------- Settings ------------------------------------#
apt='BL'                 # Appointment
IMTagSS=''               # '' or '_IM' for by-trial regressors, Stop Success
IMTagSF=''               # '' or '_IM' for by-trial regressors, Stop Failure
IMHRF='SPMG1(0)'         # HRF for any by-trial regressors
HRF='SPMG1(0)'           # HRF for pooled trial regressors
njobs=1                  #
REML='_REML'             # Use REML or deconvolve output? Set to: '' or '_REML'
pval=0.5                 # Threshold for roi voxel extraction
statsFlags=''            # Stats from deconvolves: -fout -tout -rout
bayes=0                  # Use precision-weighted averaging?

# ROI mask names (.nii.gz)
rois=(R_IFG R_preSMA R_Caudate R_GPe R_GPi R_STN R_Thal R_NAcc R_ACC L_ACC)

# Control blocks
do_setup=0
do_preproc=0
do_regress=0
do_postproc=1

# Since TCSH is a headache to do math in, these need to be consistently...
# ... but don't work since afni uses literals to escape regexes (stupid!!!)
#
#start_tr = 4     # E.g. 4 if skipping 0,1,2,3
#trs_left = 440   # total num TRs - start_tr

# This one works:
shift_stims=-4.4   # 2.2*(start_tr)

# Location of ./imagen/
imagen_path='/users/dscott3/projects/imagen/'
#---------------------------------------------------------------------------#

dataFile=${apt}_SST.nii.gz
output_dir=results_${apt}${IMTagSS}

mask=ROI_masks/All_ROIs.nii.gz
dir1D=1Ds_${apt}/

# If we're skipping setup, jump to the results dir.
if [ $do_setup = 0 ] && [ $do_preproc = 0 ]; then
   echo "Setup and preproc off, entering $output_dir"
   cd $output_dir
fi


if [ $do_setup = 1 ]; then
   rm -rf results_${apt}${IMTagSS}

   # =========================== auto block: setup ============================
   # take note of the AFNI version
   afni -ver

   # check that the current AFNI version is recent enough
   afni_history -check_date 23 Sep 2016
   if [ $status ]; then
       echo "** this script requires newer AFNI binaries (than 23 Sep 2016)"
       echo "   (consider: @update.afni.binaries -defaults)"
       exit
   fi

   mkdir $output_dir
   mkdir $output_dir/plots
   mkdir $output_dir/ideals

   # ============================ auto block: tcat ============================
   # apply 3dTcat to copy input dsets to results dir, while
   # removing the first 2 TRs
   3dTcat -prefix $output_dir/tcat ${apt}_SST.nii.gz'[2..$]'


   # ======================= Manual block: fix timing  ========================
   # enter the results directory (can begin processing data)
   cd $output_dir

   # Link stuff
   cp -r ../${dir1D} ./
   ln -s ../ROI_masks
   ln -s ../$dataFile

   # 1dcat doesn't work consistently here for some reason... (so need to copy... wtf?)
   cp ../motion_demean_${apt}.1D ./
   1dcat motion_demean_${apt}.1D'{2..$}' > motion_demean.1D

   for f in ./$dir1D/*.1D
   do
      timing_tool.py -timing ${f} -add_offset ${shift_stims} -write_timing ${f}
   done
fi

if [ $do_preproc = 1 ]; then
   rm outcount* motion_c* motion_e* motion_C* rm* bandpass*
   rm censored* despike* detrend* full* out.pre* final* zscore*

   # ================================== mask ==================================
   # create 'full_mask' data(union mask)
   3dAutomask -dilate 1 -prefix rm.mask tcat+tlrc

   # create union of inputs, output type is byte
   3dmask_tool -inputs rm.mask+tlrc.HEAD -union -prefix full_mask

   # ========================== auto block: outcount ==========================
   # data check: compute outlier fraction for each volume
   touch out.pre_ss_warn.txt
   3dToutcount -automask -fraction -polort 7 -legendre tcat+tlrc > outcount.1D

   # censor outlier TRs per run, ignoring the first 0 TRs
   # - censor when more than 0.05 of automask voxels are outliers
   # - step() defines which TRs to remove via censoring
   1deval -a outcount.1D -expr "1-step(a-0.05)" > rm.out.cen.1D

   # outliers at TR 0 might suggest pre-steady state TRs
   if [ `1deval -a outcount.1D"{0}" -expr "step(a-0.4)"` ]; then
     echo "** TR #0 outliers: possible pre-steady state TRs" \
         >> out.pre_ss_warn.txt
   fi

   # catenate outlier counts into a single time series
   cat outcount.1D > outcount_rall.1D

   # catenate outlier censor files into a single time series
   cat rm.out.cen.1D > outcount_censor.1D

   # ================== mask out very high variance voxels ===================

   # ================================ censor =================================
   # create censor file motion_${subj}_censor.1D, for censoring motion
   1d_tool.py -infile motion_demean.1D -set_nruns 1 \
              -show_censor_count -censor_prev_TR -censor_next_TR \
              -censor_motion 0.2 motion

   # combine multiple censor files
   1deval -a motion_censor.1D -b outcount_censor.1D -expr "a*b" > censor_combined.1D

   # note TRs that were not censored
   ktrs = `1d_tool.py -infile censor_combined.1D -show_trs_uncensored encoded`

   # ======================== detrend, bandpass  ==============================
   # 3dTproject dumps time points even when told not to!
   #3dTproject -polort 7 -prefix detrended \ #-cenmode NTRP -censor censor_${subj}_combined_2.1D \
   #  -passband 0.01 0.1 -TR 2.2 -input pb00.$subj.r$run.tcat+tlrc

   3dDespike -NEW25 -prefix despike tcat+tlrc
   3dBandpass -dt 2.2 -nodetrend -band 0.01 9999 -prefix bandpass despike+tlrc
   3dDetrend -polort 7 -prefix detrend bandpass+tlrc

   # ================================= scale ==================================
   # Z-Score
   3dTstat -prefix detrend_mean detrend+tlrc.
   3dTstat -prefix detrend_stdev -stdev detrend+tlrc.
   3dcalc -a detrend+tlrc -b detrend_mean+tlrc -c detrend_stdev+tlrc -expr '(a-b)/c' -prefix detrend_zscored

   # Censor
   3dcalc -a detrend_zscored+tlrc -b censor_combined.1D -expr 'a*b' -prefix censored_detrend

   # Z-Score
   3dTstat -prefix censored_mean censored_detrend+tlrc.
   3dTstat -prefix censored_stdev -stdev censored_detrend+tlrc.
   3dcalc -a censored_detrend+tlrc -b censored_mean+tlrc -c censored_stdev+tlrc -expr '(a-b)/c' -prefix zscored

   # Censor
   3dcalc -a zscored+tlrc -b censor_combined.1D -expr 'a*b' -prefix final
fi

if [ $do_regress = 1 ]; then
   # ================================ regress =================================
   3dDeconvolve \
       -force_TR 2.2 \
       -mask ${mask}                          \
       -input final+tlrc.                               \
       -TR_times 2.2                                                        \
       -GOFORIT 20 \
       -polort 0                                                            \
       -num_stimts 14                                                       \
       -censor censor_combined.1D                                 \
       -stim_times 1 ${dir1D}go_success.1D              "${HRF}" -stim_label 1 go_success     \
       -stim_times${IMTagSS} 2 ${dir1D}stop_success.1D  "${IMHRF}" -stim_label 2 stop_success   \
       -stim_times${IMTagSF} 3 ${dir1D}stop_failure.1D  "${IMHRF}" -stim_label 3 stop_failure   \
       -stim_times 4 ${dir1D}go_failure.1D              "${HRF}" -stim_label 4 go_failure     \
       -stim_times 5 ${dir1D}go_too_late.1D             "${HRF}" -stim_label 5 go_too_late    \
       -stim_times 6 ${dir1D}go_wrong_key_response.1D   "${HRF}" -stim_label 6 go_wrong_key   \
       -stim_times 7 ${dir1D}stop_too_early_response.1D "${HRF}" -stim_label 7 stop_too_early \
       -stim_file  8 motion_demean.1D'[0]' -stim_base 8  -stim_label 8  roll  \
       -stim_file  9 motion_demean.1D'[1]' -stim_base 9  -stim_label 9  pitch \
       -stim_file 10 motion_demean.1D'[2]' -stim_base 10 -stim_label 10 yaw   \
       -stim_file 11 motion_demean.1D'[3]' -stim_base 11 -stim_label 11 dS    \
       -stim_file 12 motion_demean.1D'[4]' -stim_base 12 -stim_label 12 dL    \
       -stim_file 13 motion_demean.1D'[5]' -stim_base 13 -stim_label 13 dP    \
       -stim_file 14 motion_enorm.1D       -stim_base 14 -stim_label 14 norm  \
       -gltsym 'SYM: stop_success - stop_failure'                           \
       -gltsym 'SYM: stop_success - go_success'                             \
       -glt_label 1 ss-sf                                                   \
       -glt_label 2 ss-go                                                   \
       -jobs $njobs                                                         \
       ${statsFlags} -x1D X.xmat.1D -xjpeg X.jpg                            \
       -x1D_uncensored X.nocensor.xmat.1D                                   \
       -fitts fitts                                                   \
       -errts errts                                                 \
       -bucket stats


   # if 3dDeconvolve fails, terminate the script
   if [ $status != 0 ]; then
       echo '---------------------------------------'
       echo '** 3dDeconvolve error, failing...'
       echo '   (consider the file 3dDeconvolve.err)'
       exit
   fi

   # display any large pairwise correlations from the X-matrix
   1d_tool.py -show_cormat_warnings -infile X.xmat.1D |& tee out.cormat_warn.txt

   # -- execute the 3dREMLfit script, written by 3dDeconvolve --
   #tcsh -x stats.REML_cmd
   3dREMLfit -matrix X.xmat.1D -input zscored+tlrc \
             ${statsFlags} -Rbuck stats_REML -Rvar stats_REMLvar \
             -mask ${mask} -GOFORIT 20 -Rfitts fitts_REML -Rerrts errts_REML -verb $*
fi

if [ $do_postproc = 1 ]; then
   # ============================= dump betas ===============================
   # 3dmaskave returns an error when there are no voxels, and csh has no error
   # handling, so use bash to circumvent. tcsh is such a poor choice...
   # Remove any old stats dump file.
   rm stats_masked*.out L_* R_*

   # Cycle through ROIs and extract betas using F-statistic threshold.
   for roi in ${rois[@]}
   do
      # Cycle over F-stats.
      #for fBrik in $(seq 2 2 4)
      #do
      # betabrik
      #let "betaBrik = ${fBrik} - 1"
      echo "#--------------------------------------------------------------------#"
      echo "  Dumping stats for voxels from ${roi} mask."
      echo "#--------------------------------------------------------------------#"
      rm zyxt+tlrc.*
      rm thresh_${roi}*
      rm pos_thresh_${roi}*
      rm neg_thresh_${roi}*

      # Magic numbers, sorry, see commands.
      # Go-success will have same dof as every other F-stat assuming same regressors.
      dof=`3dAttribute BRICK_STATAUX stats_REML+tlrc"[2]" | awk {'print $5'}`
      fval=`cdf -p2t fift ${pval} 1 ${dof} | awk {'print $3'}`

      # Compute mask with full F-statistic q < 0.05
      #3dFDR -input stats${REML}+tlrc.HEAD'[0]' -qval  -mask_file ROI_masks/${roi}.nii.gz

      # Create thresholded f-statistic mask
      3dcalc -a stats${REML}+tlrc.BRIK'[2..$(2)]'     \
         -b ROI_masks/${roi}.nii.gz                \
         -expr "step(a-${fval})*b" -prefix thresh_${roi}${pval}

      # Create positive-beta mask
      3dcalc -a stats${REML}+tlrc.BRIK'[1..$(2)]'        \
         -b thresh_${roi}${pval}+tlrc                    \
         -expr "step(a)*b" -prefix pos_thresh_${roi}${pval}

      # Create negative-beta mask (could have inverted pos. but this is cleaner.)
      3dcalc -a stats${REML}+tlrc.BRIK'[1..$(2)]'        \
         -b thresh_${roi}${pval}+tlrc                    \
         -expr "step(-a)*b" -prefix neg_thresh_${roi}${pval}

      if [ ${bayes} = 1 ]; then
         rm prec_* masked_prec_* wgted_pos_betas* scaled_pos_betas*

         # Precision-weighted averaging:
         # Stderr = beta/T  ->  (1/stderr)^2 = (T/beta)^2  ->  prec. = F/(beta^2)

         # Create weights mask
         3dcalc -a stats${REML}+tlrc.BRIK'[2..$(2)]'     \
            -b stats${REML}+tlrc.BRIK'[1..$(2)]'         \
            -c ROI_masks/${roi}.nii.gz                   \
            -expr "(a/(b^2))*c" -prefix prec_${roi}${pval}

         # Mask by positive values
         3dcalc -a pos_thresh_${roi}${pval}+tlrc               \
            -b prec_${roi}${pval}+tlrc                    \
            -expr "a*b" -prefix masked_prec_${roi}${pval}

         # Get the normalization factor (sum of weights)
         sums=`3dmaskave -quiet -mask ROI_masks/${roi}.nii.gz -sum masked_prec_${roi}${pval}+tlrc.BRIK.gz | tr -s '\n' ' '`

         # Generate weighted betas
         3dcalc -a stats${REML}+tlrc.BRIK'[1..$(2)]'              \
            -b masked_prec_${roi}${pval}+tlrc           \
            -expr "a*b" -prefix wgted_pos_betas${roi}${pval}

         for regNum in $(seq 1 9)
         do
            rm scaled_pos*
            let "betaNum = regNum*2 - 1"
            let "regID = regNum - 1"

            sum=`echo ${sums} | awk {'print $'${regNum}}`

            3dcalc -a wgted_pos_betas${roi}${pval}+tlrc \
               -expr "a/${sum}" -prefix scaled_pos_betas${roi}${pval}

            # Get average stats over the mask, but only for beta values
            3dmaskave -quiet -sum scaled_pos_betas${roi}${pval}+tlrc.BRIK.gz"[${regID}]" \
            | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked_pos.out

            if [ $regNum -lt 9 ]; then
               echo -n ',' >> stats_masked_pos.out
               echo -n ',' >> stats_masked_neg.out
            fi
         done
      else
         for regNum in $(seq 1 9)
         do
            let "betaNum = regNum*2 - 1"
            let "regID = regNum - 1"
            # Get average stats over the mask, but only for beta values
            3dmaskave -quiet -mask pos_thresh_${roi}${pval}+tlrc.BRIK.gz"[${regID}]" stats${REML}+tlrc.BRIK.gz"[${betaNum}]" \
            | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked_pos.out

            # Get average stats over the mask, but only for beta values
            3dmaskave -quiet -mask neg_thresh_${roi}${pval}+tlrc.BRIK.gz"[${regID}]" stats${REML}+tlrc.BRIK.gz"[${betaNum}]" \
            | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked_neg.out

            if [ $regNum -lt 9 ]; then
               echo -n ',' >> stats_masked_pos.out
               echo -n ',' >> stats_masked_neg.out
            fi
         done
      fi

      # Insert a newline after each ROI
      printf '\n' >> stats_masked_pos.out
      printf '\n' >> stats_masked_neg.out

      # Create roi response and fit
      # These need updating to play nicely with the precision weighting.
      #3dmaskave -quiet -mask pos_thresh_${roi}${pval}+tlrc.BRIK.gz final+tlrc.BRIK           > ${roi}_pos_ave.1D
      #3dmaskave -quiet -mask neg_thresh_${roi}${pval}+tlrc.BRIK.gz fitts${REML}+tlrc.BRIK.gz > ${roi}_pos_fitts.1D

      # Create copies
      cp stats_masked_pos.out ../stats_masked_pos_${apt}.out
      cp stats_masked_neg.out ../stats_masked_neg_${apt}.out

      # create ideal files for fixed response stim types
      1dcat X.nocensor.xmat.1D'[1]' > ideals/ideal_go_success.1D
      1dcat X.nocensor.xmat.1D'[2]' > ideals/ideal_stop_success.1D
      1dcat X.nocensor.xmat.1D'[3]' > ideals/ideal_stop_failure.1D
      1dcat X.nocensor.xmat.1D'[4]' > ideals/ideal_go_failure.1D
      1dcat X.nocensor.xmat.1D'[5]' > ideals/ideal_go_too_late.1D
      1dcat X.nocensor.xmat.1D'[6]' > ideals/ideal_go_wrong_key.1D
      1dcat X.nocensor.xmat.1D'[7]' > ideals/ideal_stop_too_early.1D

      # more and better plots...
      #matlab -nosplash -nodisplay -r "addpath(genpath('${imagen_path}')); plot_diagnostics; exit"
   done

#if [ $do_postproc = 1 ] & [ $IMTag = '_IM' ]; then
#      for regNum in $(seq 1 9)
#      do
#         let "betaNum = regNum*2 - 1"
#         let "regID = regNum - 1"
#         # Get average stats over the mask, but only for beta values
#         3dmaskave -quiet -mask pos_thresh_${roi}${pval}+tlrc.BRIK.gz"[${regID}]" stats${REML}+tlrc.BRIK.
#         | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked_pos.out
#
#         # Get average stats over the mask, but only for beta values
#         3dmaskave -quiet -mask neg_thresh_${roi}${pval}+tlrc.BRIK.gz"[${regID}]" stats${REML}+tlrc.BRIK.
#         | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked_neg.out
#
#         if [ $regNum -lt 9 ]; then
#            echo -n ',' >> stats_masked_pos.out
#            echo -n ',' >> stats_masked_neg.out
#         fi
#      done
#fi


   # rem ve temporary files
   \rm -f rm.*

   # return to parent directory
   cd ..

   echo "execution finished: `date`"
fi

