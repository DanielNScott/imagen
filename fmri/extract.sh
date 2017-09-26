#!/bin/bash

# Remove any old stats dump file.
rm stats_masked.out L_* R_*

# Cycle through ROIs and fit GLM
rois=( R_IFG R_preSMA R_Caudate_AAL R_GPe R_GPi R_STN R_Thalamus_AAL R_NAcc R_AAL_ACC L_AAL_ACC )
for mask in "${rois[@]}"
do
  echo "#--------------------------------------------------------------------#"
  echo "  Dumping stats for voxels from ${mask} mask."
  echo "#--------------------------------------------------------------------#"

  # Average over each voxel & append to file, transposing and stripping trailing whitespace
  3dmaskave -quiet -enorm -mask ROI_masks/${mask}.nii.gz stats.subj_REML+tlrc.BRIK.gz \
    | tr -s '\n' ' ' | sed -e 's/[[:space:]]*$//' >> stats_masked.out

  # Insert a newline after each ROI
  printf '\n' >> stats_masked.out

done
cp stats_masked.out ../
#EOF
