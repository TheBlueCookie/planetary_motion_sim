import os
import matplotlib.pyplot as plt
import numpy as np
from visu_framework import read_sim_data, plot_3d_trajectories, plot_2d_trajectories
from mpl_toolkits import mplot3d

dir_name = os.path.dirname(os.path.realpath(__file__))

data_dir = os.path.join(os.path.dirname(dir_name), 'fortran_sim', 'data')

data_subdirs = [f[0] for f in os.walk(data_dir)][1:]

fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
names = ['sun', 'mercury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptun']

# init_pos = np.array([[0, 0, 0],
#                      [64316320834.295082, 8797906245.0492535, -7637385988.6148157],
#                      [58228360555.309898, 20703430197.104145, -25891345677.120247],
#                      [211292623544.15613, -29626845299.383781, 92565325445.287781],
#                      [364209770599.89368, 344838817540.59583, -554059027197.60303]], dtype=float)
#
# for i, pos in enumerate(init_pos):
#     ax2.scatter(pos[0], pos[1], pos[2], label=names[i])

# ax2.set_xlim(min(init_pos[:, 0]), max(init_pos[:, 0]))
# ax2.set_ylim(min(init_pos[:, 1]), max(init_pos[:, 1]))
# ax2.set_zlim(min(init_pos[:, 2]), max(init_pos[:, 2]))
#

dat_ind = 0

test_data = read_sim_data(os.path.join(data_subdirs[dat_ind]))
fig1, ax1 = plt.subplots(figsize=(8, 8))

for i, dat in enumerate(test_data):
    ax1.scatter(dat.posx[0], dat.posy[0], label=names[i])
    ax2.scatter(dat.posx[0], dat.posy[0], dat.posz[0], label=names[i])
    plot_3d_trajectories(dat, ax2)
    plot_2d_trajectories(dat, ax1)

# fig2, ax2 = plt.subplots()
# ax2.scatter(test_data[1].forx, test_data[1].fory, s=2, color='tab:blue')
ax2.legend()
plt.show()
