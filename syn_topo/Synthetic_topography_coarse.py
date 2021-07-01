# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 10:36:17 2021

@author: Research_Lab
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as plt
import seaborn as sns; sns.set_theme()
from scipy.interpolate import interp2d

#%%
import pandas as pd
import os
datadir = os.path.join('C:\\Users\\jbarrett.carter\\OneDrive\\CUAHSI-SI\\Topography') # directory for some sample data files
#filename = 'Elevation_profile2.csv'
filename = "Shallotte_profile.csv"
file2 = "OceanDepth.csv"
filepath = os.path.join(datadir, filename)
filepath2 = os.path.join(datadir, file2)
df = pd.read_csv(filepath)
df2 = pd.read_csv(filepath2)
df2 = df2.loc[0:16,:]

#%%
elev = df.loc[:,'Elevation']
dist = (df.loc[:,'X_axis']-1)*3
dist = (dist-max(dist))*(-1)
elev2 = df2.loc[:,'dpeth']
dist2 = df2.loc[:,'X-axis']

elev = elev.to_numpy()
elev2 = elev2.to_numpy()
dist = dist.to_numpy()
dist2=dist2.to_numpy()

#%%
plt.plot(dist,elev)
plt.figure()
plt.plot(dist2,elev2)

#%%

#making the two match up

dist2 = (dist2-max(dist2))*(-1)
elev2 = elev2+elev[0]-elev2[0]
dist=dist+max(dist2)+1

elev = np.flip(elev)
dist = np.flip(dist)


elev = np.concatenate([elev,elev2])
dist = np.concatenate([dist,dist2])
plt.plot(dist,elev)

#%%

# make length = 1000 km

dist = dist + 10**6- max(dist)
dist =np.append(dist,0)
elev = np.append(elev,min(elev))
plt.plot(dist,elev)
#%%
#Create spline models
dist = np.flip(dist)/1000 #changed to meters
elev = np.flip(elev)
cs = CubicSpline(dist,elev)
# us = UnivariateSpline(dist, elev)

y = np.linspace(1,1000,num=1000)
syn_elev = cs(y)
# syn_elev2 = us(y)


#%%

#Generating Plots

plt.figure(figsize=(10,10))
plt.subplot(3, 1, 1)
plt.plot(dist, elev, 'bo-', linewidth=2, label='Elev Prof Site 1')
plt.plot(y, syn_elev,color = 'orange', linewidth=1, label='Syn topo smooth')
plt.xlabel('Distance (km)')
plt.ylabel('Elevation (m)')
plt.legend(loc='best')

plt.subplot(3, 1, 2)
plt.plot(dist[dist>920], elev[dist>920],'bo-', linewidth=2, label='Elev Prof Site 1')
plt.plot(y[y>920], syn_elev[y>920], color = 'orange',linewidth=1, label='Syn topo rough')
plt.xlabel('Distance (km)')
plt.ylabel('Elevation (m)')
plt.legend(loc='best')

plt.subplot(3, 1, 3)
plt.plot(dist[dist>980], elev[dist>980], 'bo-',linewidth=2, label='Elev Prof Site 1')
plt.plot(y[y>980], syn_elev[y>980], color = 'orange',linewidth=1, label='Syn topo rough')
plt.xlabel('Distance (km)')
plt.ylabel('Elevation (m)')
plt.legend(loc='best')

#%%

#Generating Plots (Univariate spline)

# plt.figure()
# plt.subplot(3, 1, 1)
# plt.plot(dist, elev, 'bo-', linewidth=2, label='Elev Prof Site 1')
# plt.plot(y, syn_elev2,color = 'orange', linewidth=1, label='Syn topo smooth')
# plt.xlabel('Distance (km)')
# plt.ylabel('Elevation (m)')
# plt.legend(loc='best')

# plt.subplot(3, 1, 2)
# plt.plot(dist[dist>900], elev[dist>900],'bo-', linewidth=2, label='Elev Prof Site 1')
# plt.plot(y[y>900], syn_elev2[y>900], color = 'orange',linewidth=1, label='Syn topo rough')
# plt.xlabel('Distance (km)')
# plt.ylabel('Elevation (m)')
# plt.legend(loc='best')

# plt.subplot(3, 1, 3)
# plt.plot(dist[dist>980], elev[dist>980], 'bo-',linewidth=2, label='Elev Prof Site 1')
# plt.plot(y[y>980], syn_elev2[y>980], color = 'orange',linewidth=1, label='Syn topo rough')
# plt.xlabel('Distance (km)')
# plt.ylabel('Elevation (m)')
# plt.legend(loc='best')

#%%

# Need to fix model in range of y = 930 to 970

cs2 = CubicSpline(dist[dist>925],elev[dist>925],bc_type='natural')
y_upper = y[y>925]
syn_elev_upper = cs2(y_upper)
plt.plot(dist[dist>980], elev[dist>980], 'bo-',linewidth=2)
plt.plot(y_upper[y_upper>980],syn_elev_upper[y_upper>980],color='orange')

syn_elev[y>925]=syn_elev_upper

#%%
# Making it 3D

x = y
z = np.zeros((len(y),len(x)))
# xyz = np.zeros((1,3))

#%%
# Defining the shorline shoreline (z=0 contour)
below_sl = syn_elev<0
base_shore = int(max(np.argwhere(below_sl)))
y_bs = y[base_shore]
x_max = max(x)
y_max = max(y)
x_max_fine = max(df['X_axis'])*3/1000

# s = np.ones(len(y))*y_bs # flat shore
# s = y_bs+2*abs(np.sin((x-500+x_max_fine/2)*np.pi/x_max_fine)) # one big curve
s = y_bs+2*abs(np.sin((x-500+x_max_fine/4)*2*np.pi/x_max_fine)) #curved with points

# add in triangular bay centered at x = x_bay
sb = np.empty((len(y)))
sb[:]=np.NaN
x_bay = x_max/2
bay_width = 4
bay_height = 8
out_bay = np.logical_or(x<(x_bay-bay_width/2),x>(x_bay+bay_width/2))
in_bay = out_bay==False
sb[out_bay]= np.NaN
sb[in_bay] = s[int(len(s)/2+bay_width/2)]+bay_height-abs(x[in_bay]-x_bay)*\
    bay_height/(bay_width/2)

# add in transition

trans = np.empty((len(y)))
trans[:]=np.NaN
trans_width = 1
out_trans = np.logical_or(x<(x_bay-bay_width/2-trans_width),x>(x_bay+bay_width/2+trans_width))
in_trans = out_trans==False
trans[out_trans]= np.NaN
trans[in_trans] = s[int(len(s)/2+(bay_width/2+trans_width))]+bay_height-\
    abs(x[in_trans]-x_bay)*bay_height/(bay_width/2)+trans_width*bay_height*2/bay_width

# add in river

river_width = bay_width/10
x_riv = np.logical_and(x > x_bay-river_width/2, x < x_bay+river_width/2)

plt.plot(x,s, label = 'shorline')
plt.plot(x,sb, label = 'bay')
plt.plot(x,trans, label = 'transition')
plt.legend()

sub = np.logical_and(x>450,x<550)
plt.figure()
plt.plot(x[sub],s[sub], label = 'shorline')
plt.plot(x[sub],sb[sub], label = 'bay')
plt.plot(x[sub],trans[sub], label = 'transition')
plt.legend()
# plt.ylim([0,30000])
#%%

col=1
# to make 2d grid of z values
for col in range(z.shape[1]):
    # z[:,col]=syn_elev*abs(np.sin(10*x[col]))
    # z[:,col]=syn_elev*x[col]
    # z[:,col]=syn_elev+5*(np.sin(x[col]*np.pi/x_max))
    # z[:,col]=syn_elev
    ynew = y-s[col]+y_bs
    ynew[ynew<0]=0
    ynew[ynew>y_max]=y_max
    z[:,col]=cs(ynew)
    if in_bay[col]==True:
        y_bay = np.logical_and(y>s[col],y<sb[col])
        z[y_bay,col]=-5
        if x_riv[col]==True:
            y_riv = y>sb[col]
            z[y_riv,col]=-5
    # if in_trans[col]==True:
    #     y_trans = np.logical_and(y>s[col],y<trans[col])
    #     y_trans2 = np.logical_and(y>s[col],y<trans[col])
    # z[:,col]=cs(y)


# to make xyz matrix
# for col in range(len(y)):
#     zs=syn_elev*abs(np.sin(10*y[col]))
#     ys = np.ones((len(dist)))
#     ys = ys*y[col]
#     new_vals = np.stack((x,ys,zs),axis = 1)
#     topo = np.append(topo,new_vals,axis=0)


#%%

#visualize cross sections

# xind = 0

# plt.plot(y[6000:7000],z[6000:7000,xind])
# plt.ylim((-10,10))

#%%

#heatmap
xis = np.linspace(450,550,num = 100,endpoint=False).astype(int)
yis = np.linspace(980,1000,num = 20,endpoint=False).astype(int)
# xis = np.linspace(0,len(x),num = 100,endpoint=False).astype(int)
# yis = np.linspace(0,len(y),num = 100,endpoint=False).astype(int)
zis = np.meshgrid(yis,xis,indexing = 'ij')
zsub = z[zis]
zsub = np.flip(zsub)
# sns.color_palette("crest", as_cmap=True)
# ax = sns.heatmap(np.flip(z),cmap='seismic',center=0)
ax = sns.heatmap(zsub,cmap='seismic',center=0)
