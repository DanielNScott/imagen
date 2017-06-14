# Imagen ML Analysis #

This repository contains files for pre-processing (QC, aggregating, outlier handling etc.), visualization, and analysis of the subset of the IMAGEN Consortium data that the LNCC has requested access to.

## Overview ##
There is a great deal of task, behavioral, clinical, and genetic data which needs...
 - cleaning
 - imputation as appropriate
 - a missingness model
 - a-priori theoretical feature selection (i.e. pruning to useful stuff)
 - modelling (e.g. BEESTS, HDDM, LATER)

Additionally, there is a great deal of fMRI data, which needs...
 - cleaning / preprocessing
 - standard regressor & statistical parametric map analysis
 - conversion to a form suitable for ancestral graph analysis

These pieces will ideally be used to...
 - determine neurological endophenotypes of impulsivity linked to fMRI
 - link such endophenotypes predictively to clinical/behavioral outcomes
 - determine the extent to which the STN and striatal dopamine seperately impact impulsivity

## References ##
TBC.

## Directories ##
preproc  - code for aggregating data from files, etc.
analysis - code for running e.g. CCAs
plotting - code for plotting

## Other Resources ##
Rik Henson fMRI Notes:
http://imaging.mrc-cbu.cam.ac.uk/imaging/DesignEfficiency
