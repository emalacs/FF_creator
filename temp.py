import pandas as pd
import numpy as np

pep_atoms = pd.read_csv('input/pep_atoms', sep = '\s+', header = None)

dat = pd.read_csv('analysis/histo-time.dat', sep = '\s+')

print(dat.to_string())

x = list(range(0, (len(dat.columns))))
print(x)
dat.columns = [x]
dat.rename(columns = {0:'frame'}, inplace = True)
print(dat)
# print(x)
y = dat[['frame']]

print(y)

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
#ax = plt.add_subplot(111, projection = '3d')

hist, xedges, yedges = np.histogram2d(x = dat, y = dat, bins = 4, range = [[0, 27], [0, 117]])

num_bars = 3248
x_pos = len(x)
y_pos = len(x)
z_pos = len(x)
x_size = x
y_size = y
z_size = dat[1:]

ax.bar3d(x_pos, y_pos, z_pos, x_size, y_size, z_size, color = 'aqua')
plt.show()

exit()
