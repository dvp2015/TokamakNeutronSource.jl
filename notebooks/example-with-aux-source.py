# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

import sys
sys.path.insert(0, r'C:\Users\Roma\OneDrive\Projects\VNC\profile')
sys.path.insert(0, r'C:\Users\Roma\Projects\aux-source')
from plasma_scenario import *
from aux_source import ptrac


# %%


data = np.loadtxt('trt.txt', delimiter='\t') * 100


# %%


fitter = Fitter(a_max = 60, shift_max=10, triang_max=0.2, ellip_max=1.9, R0 = x0+10, Z0=y0)


# %%


scenario = fitter.fit(data[:, 0].reshape(1, -1), data[:, 1].reshape(1, -1), np.array([50]), np.linspace(0, 2 * np.pi, 100))


# %%


scenario = ParametricMagneticSurface(theta=0, shift0=5, shift_slope=0, triang_x=30, triang_y=0.05, ellip0=1, ellipx=30,
                                    ellipy=1.1, shift_max=10, triang_max=0.2, ellip_max=1.9, a_max=60, R0=x0+5, Z0=y0)


# %%


xx, yy = scenario.rz_eval(50, phi)
plt.plot(xx, yy)
plt.axes().set_aspect('equal')


# %%


x0, y0


# %%


a = np.linspace(0, 60, 100)
sh = scenario.shift_eval(a)
plt.plot(a, sh)


# %%


x0, y0 = np.mean(data, axis=0)
plt.plot(data[:, 0], data[:, 1])
a_max, b_max = np.max(data, axis=0) - np.min(data, axis=0)
plt.plot(x0, y0, 'k.')
Nphi = 200
phi = np.linspace(0, 2 * np.pi, Nphi)
d = b_max / a_max - 1
npts = 1000
minor_pts = np.zeros(npts)

Nx = 25
Ny = 50
x_bins = np.linspace(np.min(data[:, 0]), np.max(data[:, 0]), Nx + 1) 
y_bins = np.linspace(np.min(data[:, 1]), np.max(data[:, 1]), Ny + 1)
x_mid = 0.5 * (x_bins[1:] + x_bins[:-1])
y_mid = 0.5 * (y_bins[1:] + y_bins[:-1])
intensity = np.zeros((Nx, Ny))
dx = np.mean(np.diff(x_bins))
dy = np.mean(np.diff(y_bins))
tol = np.sqrt(dx**2 + dy**2) / 100

x = np.zeros((npts, phi.size))
y = np.zeros((npts, phi.size))

for i, k in enumerate(np.linspace(0, 1, npts)):
    a = a_max * k / 2
    minor_pts[i] = a
    b = (b_max - (1 - k) * d * a_max) * k / 2
    xx, yy = scenario.rz_eval(a, phi)
    x[i, :] = xx #x0 + a * np.cos(phi) + 5
    y[i, :] = yy #y0 + b * np.sin(phi)
    #print(a)
    if i % 60 == 0:
        xx, yy = scenario.rz_eval(a, phi)
        plt.plot(xx, yy, 'r')
plt.axes().set_aspect('equal')

for i in range(Nx):
    print(i)
    for j in range(Ny):
        distance = np.sqrt((x - x_mid[i])**2 + (y - y_mid[j])**2)
        ind = np.argmin(distance) // Nphi
        #print(ind)
        intensity[i, j] = 1 / np.sqrt(2 * np.pi * sig**2) * np.exp(-0.5 * (minor_pts[ind] / sig)**2)


# %%


plt.pcolor(x_bins, y_bins, intensity.transpose(), cmap='jet')
plt.plot(data[:, 0], data[:, 1], 'r')
plt.axes().set_aspect('equal')


# %%


sig = 17
intensity_pts = 1 / np.sqrt(2 * np.pi * sig**2) * np.exp(-0.5 * (minor_pts / sig)**2)
plt.plot(minor_pts, intensity_pts)
plt.xlabel('a, cm')
plt.savefig('trt_profile.png', dpi=600)


# %%


from mckit import source


# %%


def create_bin_distributions(bins, start_name):
    """Creates individual distributions for every bin.

    Parameters
    ----------
    bins : array_like 
        Bin boundaries.
    start_name : int
        Starting name of the distributions.

    Returns
    -------
    free_name : int
        Distribution name, that can be used for new distributions.
    distributions : list
        A list of created distributions.
    """
    distributions = []
    for i in range(len(bins) - 1):
        distributions.append(
            source.Distribution(start_name, bins[i:i+2], [1])
        )
        start_name += 1
    return start_name, distributions


# %%


start_name=1
aux_name = start_name + 10
aux_name, x_distr = create_bin_distributions(x_bins, aux_name)
# ybins
aux_name, y_distr = create_bin_distributions(y_bins, aux_name)
probs = []
r_indices = []
z_indices = []
for i in range(Nx):
    for j in range(Ny):
        if intensity[i, j] > 0:
            r_indices.append(x_distr[i])
            z_indices.append(y_distr[j])
            probs.append(intensity[i, j])
r_dist = source.Distribution(start_name + 1, r_indices, probs, 'RAD')
z_dist = source.Distribution(start_name + 2, z_indices, r_dist, 'EXT')
src = source.Source(RAD=r_dist, EXT=z_dist)


# %% jupyter={"outputs_hidden": true}


print(src.mcnp_repr())


# %%


datap = ptrac.read_ptrac('src_test.p')
pos = np.zeros((len(datap), 3))
for i, p in enumerate(datap):
    pos[i, :] = p.position


# %%


rad = np.linalg.norm(pos[:, :2], axis=1)
ext = pos[:, 2]
index_x = np.searchsorted(x_bins, rad) - 1
index_y = np.searchsorted(y_bins, ext) - 1
probs = np.zeros_like(intensity)
out_num = 0
for i, j in zip(index_x, index_y):
    if i < x_bins.size - 1 and j < y_bins.size - 1:
        probs[i, j] += 1
    else:
        out_num += 1
print(out_num)


# %%


plt.pcolor(x_bins, y_bins, probs.transpose() / np.sum(probs), cmap='jet')
plt.plot(data[:, 0], data[:, 1], 'r')
plt.axes().set_aspect('equal')
plt.colorbar()
plt.xlabel('r, cm')
plt.ylabel('z, cm')
plt.savefig('trt_source.png', dpi=600)


# %%


plt.pcolor(x_bins, y_bins, intensity.transpose(), cmap='jet')
plt.plot(data[:, 0], data[:, 1], 'r')
plt.axes().set_aspect('equal')


# %%


distrs = source.expand_matrix_distribution(intensity, x_bins, y_bins)


# %%


data2 = np.random.permutation(data)


# %%


plt.plot(data2[:, 0], data2[:, 1])
plt.plot(data2[:, 0], data2[:, 1], 'r.')
plt.axes().set_aspect('equal')


# %%


def get_right_order(x):
    N = x.shape[0] 
    mid_point = np.mean(x, axis=0)
    deltas = x - mid_point
    phi = np.arctan2(deltas[:, 0], deltas[:, 1])
    indices = np.argsort(phi)
    data = np.zeros_like(x)
    data[:, 0] = x[indices, 0]
    data[:, 1] = x[indices, 1]
    return data


# %%


data3 = get_right_order(data2)


# %%


plt.plot(data3[:, 0], data3[:, 1], 'r.')
plt.plot(data3[:, 0], data3[:, 1], 'k')
plt.axes().set_aspect('equal')


# %%




