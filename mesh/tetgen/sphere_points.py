import random 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

num_p = 100

u_list = [ random.random() for dummy in range(num_p)]
v_list = [ random.random() for dummy in range(num_p)]
phi_list = [2*np.pi*u for u in u_list]
theta_list = [np.arccos(2*v - 1) for v in v_list]

sphere_points = [(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)) for theta in theta_list for phi in phi_list]

x_list = []
y_list = []
z_list = []
for i in range(num_p):
    theta = theta_list[i]
    phi = phi_list[i]
    x_list.append(np.sin(theta)*np.cos(phi))
    y_list.append(np.sin(theta)*np.sin(phi))
    z_list.append(np.cos(theta))

pl = plt.subplot(111, projection='3d')
pl.scatter(x_list, y_list, z_list)

plt.show()
