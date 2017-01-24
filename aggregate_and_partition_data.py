"""
Aggregates the IMAGEN subject data for analysis.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
import itertools
import json
import matplotlib
from subprocess import call
import sklearn.utils as skutils
#from sklearn.model_selection import train_test_split

import ipdb
ipdb.set_trace()

# Main routine
def main():

  # Files to read and their location:
  path = '/home/dan/documents/lncc/From Catherine/Round 2/Exported CSVs/'
  file_list = ['SST_Groups_BSL_14.csv',
               'CANTAB_14.csv',
               'CANTAB_18.csv',
               'DelayDiscounting_K_14.csv',
               'DelayDiscounting_K_16.csv',
               'DelayDiscounting_K_18.csv',
               'ESPAD_Life_14.csv',
               'ESPAD_Life_16.csv',
               'ESPAD_Life_18.csv',
               ]

  # Initialize data frame we'll merge all data into.
  data = pd.DataFrame([], index=[], columns=['Subject'])
  
  # Cycle through files to read
  for file in file_list:
    file_data = pd.read_csv(path + file, sep=',')
    age = file[-6:-4]

    # Masks for completely duplicated rows and rows with same subj ids, which
    #   shouldn't exist.
    full_duplicate_msk = file_data.duplicated()
    subj_duplicate_msk = file_data.duplicated(['Subject'])

    # Alert User to duplicate data being dropped
    if not file_data[file_data.duplicated(['Subject'])].empty:
      print('File ' + file + ' contains ' + str(sum(subj_duplicate_msk)) + ' subject duplicate fields') 
      print('File ' + file + ' contains ' + str(sum(full_duplicate_msk)) + ' fully   duplicate fields')
      print('These are: ')
      print(file_data[subj_duplicate_msk]['Subject'])

      file_data = file_data.drop_duplicates(['Subject'])

    # Change keys to reflect subject age
    for key in set(file_data.keys()) - set(['Subject']):
      file_data.rename(columns = {key: key + '_' + age}, inplace = True)

    # Change subject id numbers to integers, which they should be
    # It may be that merging on a float dimension is bad as well, in which case
    #   in which case using ints here is better too.
    file_data['Subject'] = file_data['Subject'].apply(lambda x: int(x))
    
    # Merge and report
    data = pd.merge(data, file_data, how='outer', on='Subject')
    print('Finished file ' + file)

  # Clean and write big table.
  data.drop_duplicates()
  data.to_csv('all_data.csv', sep=',')

  # Write a sub-table for just participants with valid SSRT data
  no_bsl_ssrt_msk = pd.isnull( data['GB_SSRT_14'] )

  valid_ssrt_data = data[ ~no_bsl_ssrt_msk].copy()
  valid_ssrt_data = skutils.shuffle(valid_ssrt_data, random_state=0)

  # Partition this into training, test, and safe-sub-training categories.
  n_rows   = valid_ssrt_data.shape[0]
  n_train  = round(n_rows/2)
  n_safety = round(n_train/2)

  safe_train_data = valid_ssrt_data.iloc[0:n_safety]
  safe_test_data  = valid_ssrt_data.iloc[n_safety:n_train]

  # Write
  safe_test_data.to_csv('safe_test_data.csv', sep=',')
  safe_train_data.to_csv('safe_train_data.csv', sep=',')


#
# Send system arguments to main:
if __name__ == "__main__":
    import sys
#    import time
#
#    # Check args
#    if len(sys.argv) != 1:
#      print("Incorrect number of arguments.")
#      print("Proper usage: python partition_subjects.py filename")
#      exit
#
    # Run the actual data processing
    main()
