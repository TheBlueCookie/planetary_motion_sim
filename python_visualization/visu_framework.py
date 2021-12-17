import os

import matplotlib.pyplot as plt
import pandas as pd


def read_sim_data(subdir_path: str):
    files = [f for f in os.listdir(subdir_path)]
    data = []

    for fname in files:
        fname = os.path.join(subdir_path, fname)
        data.append(pd.read_csv(fname, sep='\\s+', header=None,
                                names=['posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'accx', 'accy', 'accz', 'forx',
                                       'fory', 'forz']))
    return data


def plot_3d_trajectories(planet_data: pd.DataFrame, ax):
    ax.plot(planet_data.posx, planet_data.posy, planet_data.posz)


def plot_2d_trajectories(planet_data: pd.DataFrame, ax):
    ax.plot(planet_data.posx, planet_data.posy)
