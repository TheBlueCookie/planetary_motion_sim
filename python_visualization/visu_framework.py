import os
import pandas as pd


def read_sim_data(subdir_path: str):
    files = [f for f in os.listdir(subdir_path)]
    data = []

    for fname in files[:-1]:
        fname = os.path.join(subdir_path, fname)
        data.append(pd.read_csv(fname, sep='\\s+', header=None,
                                names=['posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'accx', 'accy', 'accz', 'forx',
                                       'fory', 'forz']))

    energies = pd.read_csv(os.path.join(subdir_path, files[-1]), sep='\\s+', header=None, names=['epot', 'ekin'])
    return data, energies


def plot_3d_trajectories(planet_data: pd.DataFrame, ax):
    ax.plot(planet_data.posx, planet_data.posy, planet_data.posz)


def plot_2d_trajectories(planet_data: pd.DataFrame, ax):
    ax.plot(planet_data.posx, planet_data.posy)


def planet_motion_with_trace(i, scat, lines, data, trace_lens, dim):
    print(f'\r{i}', end='')
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


def init_2dlines(lines):
    for line in lines:
        line.set_data([], [])


def init_3dlines(lines):
    for line in lines:
        line.set_data([], [], [])
