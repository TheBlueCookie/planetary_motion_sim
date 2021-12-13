import os
import matplotlib.pyplot as plt
from visu_framework import read_sim_data

dir_name = os.path.dirname(os.path.realpath(__file__))

data_dir = os.path.join(os.path.dirname(dir_name), 'fortran_sim', 'data')

data_subdirs = [f[0] for f in os.walk(data_dir)][1:]



test_data = read_sim_data(os.path.join(data_subdirs[7]))

fig1, ax1 = plt.subplots()
ax1.scatter(test_data[0].posx, test_data[0].posy, s=1, color='tab:orange')
ax1.scatter(test_data[1].posx, test_data[1].posy, s=1, color='tab:blue')
ax1.scatter(test_data[1].posx, test_data[2].posy, s=1, color='tab:green')

# fig2, ax2 = plt.subplots()
# ax2.scatter(test_data[1].forx, test_data[1].fory, s=2, color='tab:blue')
plt.show()
