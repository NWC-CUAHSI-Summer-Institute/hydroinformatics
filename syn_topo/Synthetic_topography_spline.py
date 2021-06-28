# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 10:36:17 2021

@author: Research_Lab
"""

import numpy as np
from scipy.interpolate import CubicSpline
from matplotlib import pyplot as plt
import seaborn as sns; sns.set_theme()
#%%
import pandas
import os
datadir = os.path.join('C:\\Users\\jbarrett.carter\\OneDrive\\CUAHSI-SI\\Topography') # directory for some sample data files
#filename = 'Elevation_profile2.csv'
filename = "site1.csv"
filepath = os.path.join(datadir, filename)
df = pandas.read_csv(filepath)

elev = df.iloc[:,1]
#%%
dist = df.iloc[:,0]
#%%
#Setting a cutoff threshold at the first 

cs = CubicSpline(dist,elev)

syn_elev = cs(dist)

#%%
# Estimating the standard deviation for the noise

SD_sample = np.std(elev-syn_elev)
noise = np.random.normal(0, .25, len(elev))
syn_elev_rough = syn_elev + noise

#%%

#Generating Plots

plt.figure(figsize=(10, 10))
plt.subplot(2, 1, 1)
plt.plot(dist, elev, linewidth=2, label='Elev Prof Site 1')
plt.plot(dist, syn_elev, linewidth=1, label='Syn topo smooth')
plt.xlabel('Distance (degrees)')
plt.ylabel('Elevation (m)')
plt.legend(loc='best')

plt.subplot(2, 1, 2)
plt.plot(dist, elev, linewidth=2, label='Elev Prof Site 1')
plt.plot(dist, syn_elev_rough, linewidth=1, label='Syn topo rough')
plt.xlabel('Distance (degrees)')
plt.ylabel('Elevation (m)')
plt.legend(loc='best')

#%%
# Making it 3D

y = dist.to_numpy()
x = np.linspace(0,5,num=len(y))
z = np.zeros((len(x),len(y)))
xyz = np.zeros((1,3))
col=1
# to make 2d grid of z values
for col in range(z.shape[0]):
    # z[:,col]=syn_elev*abs(np.sin(10*x[col]))
    # z[:,col]=syn_elev*x[col]
    z[:,col]=syn_elev+10*(np.sin(x[col]))
    # z[:,col]=syn_elev


# to make xyz matrix
# for col in range(len(y)):
#     zs=syn_elev*abs(np.sin(10*y[col]))
#     ys = np.ones((len(dist)))
#     ys = ys*y[col]
#     new_vals = np.stack((x,ys,zs),axis = 1)
#     topo = np.append(topo,new_vals,axis=0)

#%%

#visualize cross sections

cs = 1000

plt.plot(y,z[:,cs])

#%%

#heatmap
xis = np.linspace(40,4000,num = 100).astype(int)
yis = np.linspace(40,4000,num = 100).astype(int)
zis = np.meshgrid(xis,yis,indexing = 'ij')
zsub = z[zis]
ax = sns.heatmap(zsub)

#%%
# contour plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(x,y,z)

#%%

a = np.array([1,2,3])
b = np.array([4,5,6])
c = np.array([7,8,9])
d = a+b
e=a*b
