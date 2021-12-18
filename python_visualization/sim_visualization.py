import matplotlib.pyplot as plt
import numpy as np
from visu_framework import *
from matplotlib.animation import FuncAnimation

dir_name = os.path.dirname(os.path.realpath(__file__))

data_dir = os.path.join(os.path.dirname(dir_name), 'fortran_sim', 'data')

data_subdirs = [f[0] for f in os.walk(data_dir)][1:]

masses = np.array([1988500e24, 4.87e24, 3.30e24, 5.972e24, 0.642e24, 1898e24, 568e24, 86.8e24, 102e24])

max_marker_size = 20
normalized_marker_sizes = masses * max_marker_size / masses[5]

# energy plot
fig_en, ax_en = plt.subplots()
ax_en2 = ax_en.twinx()

# 2d animation
fig_2d, ax_2d = plt.subplots()

# 3d animation
fig_3d, ax_3d = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(12, 12))
names = ['sun', 'mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptun']
colors = ['gold', 'tan', 'lightcoral', 'royalblue', 'orangered', 'wheat', 'khaki', 'dodgerblue', 'darkblue']

dat_ind = 2

sim_data, energies = read_sim_data(os.path.join(data_subdirs[dat_ind]))
sim_steps = len(sim_data[0].posx)
bodies = len(sim_data)

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
    line2d, = ax_2d.plot([], [])
    line3d_1, = ax_3d.plot([], [], [], linewidth=3, color=colors[i])
    line3d_2, = ax_3d.plot([], [], [], marker='o', markersize=10, color=colors[i])
    lines2d.append(line2d)
    lines3d.append(line3d_1)
    lines3d_scat.append(line3d_2)

scat_2d = ax_2d.scatter([], [])

ax_2d.set_xlim(-4e12, 4e12)
ax_2d.set_ylim(-4e12, 4e12)

ax_3d.set_xlim(-4e12, 4e12)
ax_3d.set_ylim(-4e12, 4e12)
ax_3d.set_zlim(-4e12, 4e12)

trace_lens = [50000] * bodies


def init2d():
    scat_2d.set_offsets([])
    init_2dlines(lines2d)
    return lines2d


def init3d():
    init_3dlines(lines3d)
    init_3dlines(lines3d_scat)
    return lines3d


frames = 500

planet_anim_2d = FuncAnimation(fig_2d, planet_motion_with_trace, fargs=(scat_2d, lines2d, xy_coords, trace_lens, 2),
                               init_func=init2d, frames=100,
                               interval=10,
                               blit=True)

planet_anim_3d = FuncAnimation(fig_3d, planet_motion_with_trace, fargs=(lines3d_scat, lines3d, xyz_coords, trace_lens, 3),
                               init_func=init2d, frames=frames,
                               interval=int(sim_steps/frames),
                               blit=True)

# planet_anim_2d.save('2d_anim.gif', writer='pillow')
planet_anim_3d.save('3d_anim.gif', writer='pillow')
plt.show()