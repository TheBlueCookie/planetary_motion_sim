import os
import pandas as pd

def read_sim_data(subdir_path):
    files = [f for f in os.listdir(subdir_path)]
    data = []

    for fname in files:
        fname = os.path.join(subdir_path, fname)
        data.append(pd.read_csv(fname, sep='\\s+', header=None,
                                names=['posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'accx', 'accy', 'accz', 'forx',
                                       'fory', 'forz']))
    return data