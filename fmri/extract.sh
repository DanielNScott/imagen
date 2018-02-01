#!/bin/bash

# Remove any old stats dump file.
rm stats_masked.out L_* R_*

# Cycle through ROIs and fit GLM
#rois=(All_ROIs)
rois=(R_IFG R_preSMA R_Caudate_AAL R_GPe R_GPi R_STN R_Thalamus_AAL R_NAcc R_AAL_ACC L_AAL_ACC)
for mask in "${rois[@]}"
do
  echo "#--------------------------------------------------------------------#"
  echo "  Dumping stats for voxels from ${mask} mask."
  echo "#--------------------------------------------------------------------#"
  rm zyxt+tlrc.*

  3dFDR -input stats.s01_REML+tlrc.HEAD'[0]' -qval -mask_file ROI_masks/${mask}.nii.gz
  3dcalc -a zyxt+tlrc.BRIK'[0]' -expr 'step(0.05-a)' -prefix cluster_mask_${mask}
  3dmaskave -quiet -mask cluster_mask_${mask}+tlrc.BRIK stats.s01_REML+tlrc.BRIK.gz \
    | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked.out

  # Insert a newline after each ROI
  printf '\n' >> stats_masked.out

  #3dclust -savemask clust_mask -1clip 0.95 -NN3 1 FFDRq.s01_REML+tlrc.BRIK.gz'[0]'
  #3dmask_tool -inputs cluster_mask+tlrc.HEAD ROI_masks/${mask}.nii.gz -intersection -prefix ${mask}_cluster_mask
  # Average over each voxel & append to file, transposing and stripping trailing whitespace
  #3dmaskave -quiet -mask ROI_masks/${mask}.nii.gz stats.s01_REML+tlrc.BRIK.gz \
  #  | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked.out

done

mv stats_masked.out stats_masked.out.bckp

# Remove the pesky F-statistic in the middle of the stop-condition block...
#3dinfo -verb stats.subj_REML+tlrc.HEAD > 3dinfo.txt
#badCol=`grep 'stop_success_Fstat' 3dinfo.txt | cut -d ' ' -f 6 | tr '#' ' '`
#let "badCol=${badCol}+1"
#cut -d ' ' -f ${badCol} stats_masked.out.bckp --complement > stats_masked.out
#cp stats_masked.out ../

cp stats_masked.out.bckp ../stats_masked.out
#EOF
