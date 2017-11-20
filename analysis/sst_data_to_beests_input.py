"""
Process IMAGEN stop signal task data into BEESTS-friendly form.
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

import ipdb
ipdb.set_trace()

# Main routine to import the primary csv and do some preliminary parsing
def main():

  # Where is the list of subjects' data?
  train_file = '/home/dan/projects/imagen/data/archived/synthetic data/safe_train_data.csv'
  test_file  = '/home/dan/projects/imagen/data/archived/synthetic data/safe_test_data.csv'

  # Read the list of subjects who will be in 'training' and 'test' sets
  subj_list_train = pd.read_csv(train_file, sep=',', usecols=['Subject'])
  subj_list_test  = pd.read_csv(test_file , sep=',', usecols=['Subject'])

  # Format to 12 digit leading-zero-padded integer strings.
  subj_list_train = subj_list_train['Subject'].apply(lambda x: str(int(x)).zfill(12))
  subj_list_test  = subj_list_test['Subject'].apply(lambda x: str(int(x)).zfill(12))

  # Data frame in which to accumulate data from subjects' files.
  df_test  = pd.DataFrame()
  df_train = pd.DataFrame()

  # Fieldnames to read.
  # Trailing white space on item 1 is not an error here.
  fields = ["Go Stimulus Presentation Time " ,
            "Stop Stimulus Presentation Time",
            "Stimulus Presented"             ,
            "Response made by subject"       ,
            "Relative Response Time"         ,
            "Response Outcome"]

  # Location of the csv files to be read:
  data_dir = '/home/dan/projects/imagen/data/archived/From Catherine/SST/SST_BL_behavioraldata'

  # For the iterator below...
  test_tuple  = (df_test , subj_list_test , 'test' )
  train_tuple = (df_train, subj_list_train, 'train')

  # Two passes, one for train and one for test
  for df, subj_list, file_label in [test_tuple, train_tuple]:

    f_num = 0
    # Read subject files:
    for padded_subj_str in subj_list:

      # Some files may be corrupted, so try-catch:
      try:
        f_num = f_num + 1
        file  = data_dir + '/' + 'ss_' + padded_subj_str + '.csv'
        frame = pd.read_csv(file, header=1, sep='\t', usecols=fields)
        
        sid = pd.Series(np.ones([frame.shape[0]]) * f_num, index=frame.index)
        sid = sid.astype('int')

        ssd = frame["Stop Stimulus Presentation Time"] - frame["Go Stimulus Presentation Time " ]
        rt  = frame["Relative Response Time"].copy()

        # Create boolean mask and indicator column for stop signals
        ss_bool = ssd.copy()
        ss_msk  = ~ss_bool.isnull().copy()

        ss_bool[~ss_msk] = 0
        ss_bool[ ss_msk] = 1
        ss_bool = ss_bool.astype('int')

        # Boolean inhibition success column
        # 
        # -999: No data since no stop signal
        #    0: Failure to inhibit
        #    1: Successful inhibition 
        inhib = ss_bool.copy()
        inhib[~ss_msk] = -999
        inhib[ss_msk & ~rt.isnull()] = 0
        inhib = inhib.astype('int')

        # Format for output
        rt[rt.isnull()]   = -999
        ssd[ssd.isnull()] = -999
        
        rt  = rt.astype('int')
        ssd = ssd.astype('int')

        col_names = ['subj_idx', 'ss_presented', 'inhibited', 'ssd', 'rt']
        tmp = pd.concat([sid, ss_bool, inhib, ssd, rt], axis=1)
        tmp.columns = col_names

        # SSD has some bad rows so they should be removed, according to the
        # following three issues:
        #
        # 1) The ssd is < 0 which makes no sense
        # 2) The RT is unrealistically fast
        # 3) BEESTS can't handle no-response data when no signal was present.
        min_rt_allowed = 150

        bad_inds = ssd[ss_msk][ssd[ss_msk] < 0].index
        bad_inds = bad_inds.append(tmp[rt != -999][rt < min_rt_allowed].index)
        bad_inds = bad_inds.append(tmp[ss_bool == 0][rt == -999].index)

        if len(bad_inds) > 0:
          #print(' ')
          #print('The following bad data will be dropped from ' + file + ' data:')
          #print(tmp.iloc[bad_inds])

          loss_prct = round(100*len(bad_inds)/len(tmp))
          #print('This represents ' + str(loss_prct) + '% data loss.')

          tmp = tmp.drop(bad_inds)

          #print(' ')
          #print('RT quantiles:')
        #else:
          #print(' ')
          #print('RT quantiles for ' + file + ' data:')

        # Report max, min, and quartiles for each dataset.
        #print( rt[rt != -999].quantile([0, 0.25, 0.5, 0.75, 1]) )

        # Save the current data into the output dataframe
        # tmp.sort_values(col_names, inplace=True)
        df = pd.concat([df, tmp])

        print('Finished file: ' + file.split('/')[-1] + ' (file num ' + str(f_num) + ')')

      except:
        print('Failure to read or process file: ' + file.split('/')[-1])
        continue

    # Once everything is loaded into the output dataframe, save it as a csv.
    print('Saving sst_data_' + file_label + '.csv')
    df.to_csv('sst_data_' + file_label + '.csv', index=False)


# Send system arguments to main:
if __name__ == "__main__":
  import sys
  import time

  # Check args
  #if len(sys.argv) != 3:
  # print("Incorrect number of arguments.")
  # print("Proper usage: python process_sst_data data_dir file_limit")
  # exit

  # Run the actual data processing
  main()
