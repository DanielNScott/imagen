import os
import pandas as pd
from glob import glob

try:
    os.makedirs('/home/dan/projects/imagen/data/archived/Recieved/MID/preprocessed/MID_BL_behavioraldata/')
except:
    pass

# remove any quotes from the file
files = glob('/home/dan/projects/imagen/data/archived/Recieved/MID/MID_BL_behavioraldata/mid_*.csv')
for ii, fname in enumerate(files):
    print('Preprocessing file: ' + fname)
    with open(fname) as f:
        content = f.readlines()
        _content = ''
        for line in content:
            line = line.replace('"', '')
            line = line.replace("'", '')

            _content += line
    idx = fname.find('/mid_') + 1
    with open('/home/dan/projects/imagen/data/archived/Recieved/MID/preprocessed/MID_BL_behavioraldata/' + fname[idx:], 'w') as f:
        f.write(_content)

files = glob('/home/dan/projects/imagen/data/archived/Recieved/MID/raw/MID_BL_behavioraldata/mid_*.csv')
for ii, fname in enumerate(files):
    with open(fname) as f:
        content = f.readlines()
        _content = ''
        for line in content:
            if line.find('"') > -1:
                print(ii, line, fname)
                raise Exception
