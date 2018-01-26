"""
Process IMAGEN stop signal task data into BEESTS-friendly form.
"""
import pandas as pd
import numpy as np
from copy import copy
import itertools
import ipdb

# Main routine to import the primary csv and do some preliminary parsing
def main():

  # Training and test subject lists, as files
  train_file = '/home/dan/projects/imagen/data/protected/sandbox_subject_train.csv'
  test_file  = '/home/dan/projects/imagen/data/protected/sandbox_subject_test.csv'

  # Read the list of subjects who will be in 'training' and 'test' sets
  subj_list_train = pd.read_csv(train_file, sep=',', usecols=['Subject'])
  subj_list_test  = pd.read_csv(test_file , sep=',', usecols=['Subject'])

  # Format to 12 digit leading-zero-padded integer strings.
  subj_list_train = subj_list_train['Subject'].apply(lambda x: str(int(x)).zfill(12))
  subj_list_test  = subj_list_test['Subject'].apply(lambda x: str(int(x)).zfill(12))

  # Data frame in which to accumulate data from subjects' files.
  df_test  = pd.DataFrame()
  df_train = pd.DataFrame()

  misc_metrics_test  = pd.DataFrame()
  misc_metrics_train = pd.DataFrame()

  # Fieldnames to read.
  # Trailing white space on item 1 is not an error here.
  fields = ["Go Stimulus Presentation Time " , "Stop Stimulus Presentation Time",
    "Stimulus Presented", "Response made by subject", "Relative Response Time",
    "Response Outcome", "Delay"]

  # Location of the csv files to be read:
  data_dir = '/home/dan/projects/imagen/data/archived/recieved/SST/SST_BL_behavioraldata'

  # For the iterator below...
  test_tuple  = (df_test , misc_metrics_test , subj_list_test , 'test' )
  train_tuple = (df_train, misc_metrics_train, subj_list_train, 'train')

  # Two passes, one for train and one for test
  for df, misc, subj_list, file_label in [test_tuple, train_tuple]:

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
      res = frame["Response Outcome"]
      dly = frame["Delay"]

      # Create boolean mask and indicator column for stop signals
      ss_bool = ssd.copy()
      ss_msk  = ~ss_bool.isnull().copy()

      ss_bool[~ss_msk] = 0
      ss_bool[ ss_msk] = 1
      ss_bool = ss_bool.astype('int')

      # Some ss_presented indicators will be wrong for our purposes:
      # Responding during the SSD period should be considered "success."
      too_early = res == "STOP_TOO_EARLY_RESPONSE"
      ssd[too_early] = np.nan
      ss_bool[too_early] = 0

      # Some SSDs are also computed wrong when calculated above, for example from row 34
      # of the file /SST_FU2_behavioraldata/ss_000036529694.csv, which gives a negative 
      # SSD, with value -9178. Hence, we mask and check against the 'Delay' column
      bad_ssd_msk = ssd < 0
      bad_ssd = ssd[bad_ssd_msk].copy()
      ssd[bad_ssd_msk] = dly[ss_msk]

      # There is also an issue with subjects having no recorded responses,
      # as in ss_000024585027.csv.
      # ToDo: Deal with this case

      # The output of this file is not in a useful format for creating other
      # diagnostics such as no-response rates
      # ToDo: determine what to do about no-response rate quantification
      # ToDo: determine what about min-RT diagnostics, cleaning, etc... 
      #       ... such as calculating a fraction of total trials in which 
      #       ... subjects are responding too quickly to be really doing the task

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

      # BEESTS can't handle no-response data when no signal was present.
      drop_inds = tmp[ss_bool == 0][rt == -999].index

      min_rt   = 100
      way_fast = tmp[rt != -999][rt < min_rt]

      no_response_frac = 0
      if len(drop_inds) > 0:
        print(' ')
        print('The following bad data will be dropped from ' + file + ' data:')
        print( pd.concat([tmp.iloc[drop_inds], res.iloc[drop_inds]],1) )

        no_response_frac = len(drop_inds)/len(tmp[ss_bool == 0 ])
        no_response_prct = round(150*no_response_frac)
        print('This represents ' + str(no_response_prct) + '% data loss from go trials.')

      low_rt_frac = 0
      if len(way_fast > 0):
        print(' ')
        low_rt_frac = len(way_fast)/len(tmp[ss_bool == 0 ])
        low_rt_prc  = round(100*low_rt_frac)
        print('This participant also exhibits ' + str(low_rt_prc) + '% trials with RTs < 100')

      pre_signal_frac = 0
      if len(ssd[too_early]) > 0:
        print(' ')
        pre_signal_frac = len(ssd[too_early])/(len(tmp[ss_bool == 1]) + len(ssd[too_early]))
        pre_signal_prc  = round(100*pre_signal_frac)
        print('This participant also exhibits ' + str(pre_signal_prc) + '% trials pre-signal responding')

      # Bad SSD computations s
      if len(ssd[bad_ssd_msk]) > 0:
        print(' ')
        print('Additionally, the following bad SSDs were discovered:')
        print(pd.concat( [bad_ssd, ssd[bad_ssd_msk]], 1))

      # Drop fields and append data
      tmp = tmp.drop(drop_inds)
      df = pd.concat([df, tmp])

      # Save additional subject diagnostics
      cols = ('Subject','no_response_frac', 'pre_signal_frac', 'low_rt_frac')
      row  = [[int(padded_subj_str),no_response_frac,pre_signal_frac,low_rt_frac]]
      misc = pd.concat([misc, pd.DataFrame(data = row, columns = cols)])

      print('Finished file: ' + file.split('/')[-1] + ' (file num ' + str(f_num) + ')')

    except:
      print('Failure to read or process file: ' + file.split('/')[-1])
      continue

   # Once everything is loaded into the output dataframe, save it as a csv.
   print('Saving sst_data_' + file_label + '.csv')
   df.to_csv('sst_data_' + file_label + '.csv', index=False)

   print('Saving misc_sst_data_' + file_label + '.csv')
   misc.to_csv('misc_sst_data_' + file_label + '.csv', index=False)

# Send system arguments to main:
if __name__ == "__main__":

  # Run the actual data processing
  main()
