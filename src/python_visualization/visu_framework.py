import os
import pandas as pd
import numpy as np
from scipy.signal import find_peaks, argrelmin
from scipy.optimize import curve_fit


def read_sim_data(subdir_path: str):
    files = [f for f in os.listdir(subdir_path)]
    data = []

    for fname in files[:-1]:
        fname = os.path.join(subdir_path, fname)
        data.append(pd.read_csv(fname, sep='\\s+', header=None,
                                names=['posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'accx', 'accy', 'accz', 'forx',
                                       'fory', 'forz', 'time']))

    energies = pd.read_csv(os.path.join(subdir_path, files[-1]), sep='\\s+', header=None,
                           names=['epot', 'ekin', 'time'])
    return data, energies


def plot_3d_trajectories(planet_data: pd.DataFrame, ax):
    ax.plot(planet_data.posx, planet_data.posy, planet_data.posz)


def plot_2d_trajectories(planet_data: pd.DataFrame, ax):
    ax.plot(planet_data.posx, planet_data.posy)


def planet_motion_with_trace(i, scat, lines, data, trace_lens, dim, fig, time):
    print(f'\r{i}', end='')
    fig.suptitle(f'time = {(time[i] / (3600 * 24)):.2f} earth days')
    if dim == 2:
        scat.set_offsets(data[0:8, :, i])
    elif dim == 3:
        for ind, dat in enumerate(data):
            update_traced_line(i, dat, scat[ind], ind, [1] * len(trace_lens), dim)
    for ind, dat in enumerate(data):
        update_traced_line(i, dat, lines[ind], ind, trace_lens, dim)
    return lines


def update_traced_line(i, data, line, ind, trace_lens, dim):
    trace_len = trace_lens[ind]
    if dim == 2:
        if i <= trace_len:
            line.set_data(data[0, 0:i], data[1, 0:i])
        else:
            line.set_data(data[0, i - trace_len:i], data[1, i - trace_len:i])

    elif dim == 3:
        if i <= trace_len:
            line.set_data(data[0, 0:i], data[1, 0:i])
            line.set_3d_properties(data[2, 0:i])
        else:
            line.set_data(data[0, i - trace_len:i], data[1, i - trace_len:i])
            line.set_3d_properties(data[2, i - trace_len:i])


def init_lines(lines):
    for line in lines:
        line.set_data([], [])


def get_orbital_periods_fft(r_data, time):
    orb_periods = np.zeros(len(r_data))
    for i, dat in enumerate(r_data):
        r_fft = np.abs(np.fft.fft(dat))
        peaks = find_peaks(r_fft)
        if not peaks[0].size == 0:
            orb_periods[i] = peaks[0][0]
        else:
            orb_periods[i] = np.nan

    orb_periods = 1 / (orb_periods * 24 * 3600) * np.max(time)

    return orb_periods


def get_orbital_periods_avg_mins(r_data, time):
    orb_periods = np.zeros(len(r_data))
    min_inds = []
    for i, dat in enumerate(r_data):
        extrema = argrelmin(dat)
        if not extrema[0].size == 0:
            ex_days = list(time[list(extrema[0])] / (24 * 3600))
            min_inds.append(ex_days)
            periods = [np.abs(p - ex_days[i + 1]) for i, p in enumerate(ex_days[:-1])]
            orb_periods[i] = np.mean(periods)
        else:
            min_inds.append([])
            orb_periods[i] = np.nan

    return orb_periods, min_inds


def get_linear_drift(r_data, time):
    popt, pcov = curve_fit(lambda x, a, b: a * x + b, time, r_data)
    return popt, pcov
