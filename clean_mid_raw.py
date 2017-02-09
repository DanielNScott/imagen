import os
import pandas as pd
from glob import glob

try:
    os.makedirs('./preprocessed/MID_FU2_behavioraldata/')
except:
    pass
    
# remove any quotes from the file
files = glob('./raw/MID/MID_FU2_BehavData/mid_*.csv')
for ii, fname in enumerate(files):
    with open(fname) as f:
        content = f.readlines()
        _content = ''
        for line in content:
            line = line.replace('"', '')
            
            _content += line
    idx = fname.find('/mid_') + 1
    with open('./preprocessed/MID_FU2_behavioraldata/' + fname[idx:], 'w') as f:
        f.write(_content)
        
files = glob('./raw/MID_FU2_behavioraldata/mid_*.csv')
for ii, fname in enumerate(files):
    with open(fname) as f:
        content = f.readlines()
        _content = ''
        for line in content:
            if line.find('"') > -1:
                print ii, line, fname
                raise Exception
