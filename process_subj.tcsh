#!/bin/tcsh -xef

set subj = ''
set dir1D = 1Ds/
set mask = ./ROI_masks/All_ROIs.nii.gz
set dataFile = BL_SST.nii.gz
set outDir = ./Results/

#set subj = 000051571318
#set dir1D = ${subj}_1Ds/
#set dataFile = ${subj}_BL_SST.nii.gz
#set outDir = ./${subj}_Results/
#rm -r
#mkdir

3dDeconvolve -input ${dataFile}                     \
             -jobs 16                               \
             -polort 3                              \
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
             -stim_times 4 ${dir1D}go_failure.1D              'SPMG3'   -stim_label 4 go_failure      \
             -stim_times 5 ${dir1D}go_too_late.1D             'SPMG3'   -stim_label 5 go_too_late     \
             -stim_times 6 ${dir1D}go_wrong_key_response.1D   'SPMG3'   -stim_label 6 go_wrong_key    \
             -stim_times 7 ${dir1D}stop_too_early_response.1D 'SPMG3'   -stim_label 7 go_too_early   \
             -stim_file  8 ${dir1D}motion_demean.1D'[0]' -stim_base 8  -stim_label 8  roll  \
             -stim_file  9 ${dir1D}motion_demean.1D'[1]' -stim_base 9  -stim_label 9  pitch \
             -stim_file 10 ${dir1D}motion_demean.1D'[2]' -stim_base 10 -stim_label 10 yaw   \
             -stim_file 11 ${dir1D}motion_demean.1D'[3]' -stim_base 11 -stim_label 11 dS    \
             -stim_file 12 ${dir1D}motion_demean.1D'[4]' -stim_base 12 -stim_label 12 dL    \
             -stim_file 13 ${dir1D}motion_demean.1D'[5]' -stim_base 13 -stim_label 13 dP    \
             -num_glt 2                                \
             -glt_label 1 ss_go_diff                    \
             -gltsym 'SYM: +stop_success -go_success'   \
             -glt_label 2 ss_sf_diff                      \
             -gltsym 'SYM: +stop_success -stop_failure' \
             -fitts fitts.$subj                         \
             -errts errts.${subj}                       \
             -fout -tout -x1D X.xmat.1D -xjpeg X.jpg    \
             -tout -fout                                \
             -bucket stats.$subj                        \
             -x1D_stop

# Invoke REML fit
3dREMLfit -matrix X.xmat.1D -input ${dataFile} \
 -mask ${mask} -fout -tout -GOFORIT 3 \
 -Rbuck stats.REML -Rvar stats.REMLvar \
 -Rfitts fitts.REML -Rerrts errts.REML -verb $*


# -- use 3dTproject to project out regression matrix --
3dTproject -polort 0 -input ${dataFile} \
           -ort X.xmat.1D -prefix errts.${subj}.tproject

# -- Send voxelwise statistics off into the world -- #

set rois = (R_IFG R_preSMA R_Caudate_AAL R_GPe R_GPi R_STN R_Thalamus_AAL R_NAcc R_AAL_ACC L_AAL_ACC)
foreach mask (${rois})
   3dmaskdump -noijk -mask ROI_masks/$mask.nii.gz stats.REML+tlrc.HEAD > stats_masked_$mask.out
end

#EOF
