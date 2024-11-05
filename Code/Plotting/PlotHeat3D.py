# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# # Create a grid of data points
# x = np.linspace(-5, 5, 50)
# y = np.linspace(-5, 5, 50)
# x, y = np.meshgrid(x, y)

# # Z values as a function of X and Y
# z = np.sin(np.sqrt(x**2 + y**2))

# # Create figure and 3D axis
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Plot 3D surface with a heatmap
# heatmap = ax.plot_surface(x, y, z, cmap='viridis')

# # Add color bar for the heatmap
# fig.colorbar(heatmap)

# plt.show()


import matplotlib.pyplot as plt
import numpy as np

# different from yours, see below
x = y = z = np.linspace(-2, 2, 41)
X, Y, Z = np.meshgrid(x, y, z)
values = 2*X*X - Y*Y + 1/(Z*Z+1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(X, Y, Z, c=values, cmap='PRGn')
fig.colorbar(scatter, ax=ax)

plt.show()