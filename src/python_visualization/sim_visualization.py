import numpy

import matplotlib.pyplot as plt
import numpy as np
from visu_framework import *
import scipy.signal
from matplotlib.animation import FuncAnimation

dir_name = os.path.dirname(os.path.realpath(__file__))

data_dir = os.path.join('..', '..', 'run', 'data')

data_subdirs = [f[0] for f in os.walk(data_dir)][1:]

masses = np.array([1988500e24, 4.87e24, 3.30e24, 5.972e24, 0.642e24, 1898e24, 568e24, 86.8e24, 102e24])

max_marker_size = 15
normalized_marker_sizes = masses * max_marker_size / masses[5]

plt.rcParams.update({'font.size': 23})

# distance from origin plot
fig_r, ax_r = plt.subplots(figsize=(13, 10))
ax_r2 = ax_r.twinx()

# fft plot
fig_fft, ax_fft = plt.subplots()

# energy plot
fig_en, ax_en = plt.subplots()
ax_en2 = ax_en.twinx()

# 2d animation
fig_2d, ax_2d = plt.subplots()

# 3d animation
fig_3d, ax_3d = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(12, 12))

moon = False
save_anim = False
dat_ind = 10
sel_inds = [0, 5] #[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
spacing_fac = 1.03
rv_ind = 1
frame_number = 250

if moon:
    names = ['sun', 'mercury', 'venus', 'earth', 'moon', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
    colors = ['gold', 'tan', 'lightcoral', 'royalblue', 'darkgray', 'orangered', 'wheat', 'khaki', 'dodgerblue', 'darkblue']
    size = [30, 10, 10, 10, 5, 10, 15, 15, 10, 10]
else:
    names = ['sun', 'mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
    colors = ['gold', 'tan', 'lightcoral', 'royalblue', 'orangered', 'wheat', 'khaki', 'dodgerblue', 'darkblue']
    size = [30, 10, 10, 10, 10, 15, 15, 10, 10]

print(data_subdirs[dat_ind])

sim_data, energies = read_sim_data(os.path.join(data_subdirs[dat_ind]))
sim_steps = len(sim_data[0].posx)
print(sim_steps)
bodies = len(sim_data)

r_planets = np.zeros((bodies, sim_steps))
v_planets = np.zeros((bodies, sim_steps))

x_points = np.arange(0, sim_steps)
# ax_en.plot(x_points, energies.epot, label='potential energy')
# ax_en2.plot(x_points, energies.ekin, label='kinetic energy')
ax_en.plot(x_points, energies.epot + energies.ekin, label='total energy')
ax_en.legend()

xy_coords = np.zeros((bodies, 2, sim_steps))
xyz_coords = np.zeros((bodies, 3, sim_steps))

lines2d = []
lines3d = []
lines3d_scat = []

for i, dat in enumerate(sim_data):
    xy_coords[i, 0, :] = dat.posx
    xy_coords[i, 1, :] = dat.posy
    xyz_coords[i, 0, :] = dat.posx
    xyz_coords[i, 1, :] = dat.posy
    xyz_coords[i, 2, :] = dat.posz
    print(names[sel_inds[i]], colors[sel_inds[i]])
    line2d, = ax_2d.plot([], [], label=names[sel_inds[i]])
    line3d_1, = ax_3d.plot([], [], [], linewidth=2, color=colors[sel_inds[i]], label=names[sel_inds[i]])
    line3d_2, = ax_3d.plot([], [], [], marker='o', markersize=size[sel_inds[i]], color=colors[sel_inds[i]])
    lines2d.append(line2d)
    lines3d.append(line3d_1)
    lines3d_scat.append(line3d_2)
    if i != 0:
        r_planets[i, :] = [np.linalg.norm(
            [dat.posx[k] - sim_data[0].posx[k], dat.posy[k] - sim_data[0].posy[k], dat.posz[k] - sim_data[0].posz[k]])
            for k in range(sim_steps)]
    else:
        r_planets[i, :] = [np.linalg.norm([dat.posx[k], dat.posy[k], dat.posz[k]]) for k in range(sim_steps)]
    v_planets[i, :] = [np.linalg.norm([dat.velx[k], dat.vely[k], dat.velz[k]]) for k in range(sim_steps)]

days_time = np.copy(energies.time / (3600 * 24))

print(get_orbital_periods_fft(r_planets, energies.time))

if rv_ind == 0:
    popt, pcov = get_linear_drift(r_planets[rv_ind], days_time)
    correction = np.copy(days_time * popt[0] + popt[1])
    # ax_r.plot(days_time, correction, color='black', linestyle='dashed', label='linear drift')
    r_planets[rv_ind] -= correction

orbital_periods, min_indices = get_orbital_periods_avg_mins(r_planets, energies.time)
print(orbital_periods, min_indices)


ax_r_min = np.min(r_planets[rv_ind])
ax_r_max = np.max(r_planets[rv_ind])

ax_r2_min = np.min(v_planets[rv_ind])
ax_r2_max = np.max(v_planets[rv_ind])

ax_r.set_ylim(ax_r_min, ax_r_max * spacing_fac)
ax_r2.set_ylim(ax_r2_min, ax_r2_max * spacing_fac)

if len(min_indices[rv_ind]) == 1:
    ax_r.vlines([0] + min_indices[rv_ind], ax_r_min, ax_r_max, colors='gray', linestyles='dashed',
                label=f'mean orbital period: {min_indices[rv_ind][0]:.2f}')
else:
    pass
    ax_r.vlines(min_indices[rv_ind], ax_r_min, ax_r_max, colors='gray', linestyles='dashed',
                label=f'mean orbital period: {orbital_periods[rv_ind]:.2f}')

ax_r.plot(energies.time / (24 * 3600), r_planets[rv_ind], color='tab:blue', label='distance')
ax_r2.plot(energies.time / (24 * 3600), v_planets[rv_ind], color='tab:orange', label='velocity')

ax_r.legend(loc=1)
ax_r2.legend(loc=2)
ax_r.set_xlabel('earth days')
ax_r.set_ylabel('$[m]$')
ax_r2.set_ylabel('$[m/s]$')
fig_r.suptitle(f'{names[sel_inds[rv_ind]]}: scalar distance to sun and absolute velocity')

scat_2d = ax_2d.scatter([], [])

ax_2d.set_xlim(np.min(xy_coords[:, 0, :]) * spacing_fac, np.max(xy_coords[:, 0, :]) * spacing_fac)
ax_2d.set_ylim(np.min(xy_coords[:, 1, :]) * spacing_fac, np.max(xy_coords[:, 1, :]) * spacing_fac)

ax_3d.set_xlim(np.min(xyz_coords[:, 0, :]) * spacing_fac, np.max(xyz_coords[:, 0, :]) * spacing_fac)
ax_3d.set_ylim(np.min(xyz_coords[:, 1, :]) * spacing_fac, np.max(xyz_coords[:, 1, :]) * spacing_fac)
ax_3d.set_zlim(np.min(xyz_coords[:, 2, :]) * spacing_fac, np.max(xyz_coords[:, 2, :]) * spacing_fac)

trace_lens = [50] * bodies
# trace_lens = [500, 50000000, 50, 500000, 10, 50] * bodies


def init2d():
    scat_2d.set_offsets([])
    init_lines(lines2d)
    return lines2d


def init3d():
    init_lines(lines3d)
    init_lines(lines3d_scat)
    return lines3d

frames = np.arange(0, sim_steps, int(sim_steps / frame_number))

# planet_anim_2d = FuncAnimation(fig_2d, planet_motion_with_trace,
#                                fargs=(scat_2d, lines2d, xy_coords, trace_lens, 2, fig_2d, energies.time),
#                                init_func=init2d, frames=frames,
#                                interval=1,
#                                blit=True, repeat=False)
#
planet_anim_3d = FuncAnimation(fig_3d, planet_motion_with_trace,
                               fargs=(lines3d_scat, lines3d, xyz_coords, trace_lens, 3, fig_3d, energies.time),
                               init_func=init3d, frames=frames,
                               interval=1,
                               blit=True, repeat=False)

ax_2d.legend()
ax_3d.legend()
#
if save_anim:
    # planet_anim_2d.save('2d_anim.gif', writer='pillow')
    planet_anim_3d.save('3d_anim.gif', writer='pillow')
plt.show()
