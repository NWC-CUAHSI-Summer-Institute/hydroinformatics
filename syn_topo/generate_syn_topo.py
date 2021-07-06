# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:01:00 2021

@author: J. Barrett Carter

This is a function that will be used to automatically generate synthetic topo-
graphy from a specified elevation profile.

The 'site_file_name' variable must be a csv file containing a column for distance
labeled 'Distance' and a column for elevation labeled 'Elevation'. The distance
should be measured from a point on land to a point at sea, and the 'Distance' values
should be orded from smallest to largest.

The options for the shape parameter are either 'straight', 'curved', or 'points'. The options
for the bay parameter are True or False. The options for the plot parameter are
0, 1, 2, or 3; 0 meaning no plot, 1 meaning just profile plot, 2 meaning just
heatmap, 3 meaning both.

"""
def generate_topo(site_file_name,directory,shape,bay = False,plot = 3):

    import numpy as np
    from scipy.interpolate import CubicSpline
    from scipy.interpolate import interp1d
    from matplotlib import pyplot as plt
    import seaborn as sns; sns.set_theme()
    
    #%%
    import pandas as pd
    import os
    datadir = os.path.join(directory) # directory for elevation profile data files
    filename = site_file_name
    filepath = os.path.join(datadir, filename)
    df = pd.read_csv(filepath)
    
    #%%
    elev = df.loc[:,'Elevation']
    dist = df.loc[:,'Distance']
    elev = elev[0:10001]
    dist = dist[0:10001]
    elev = elev.to_numpy()
    elev = np.flip(elev)
    dist = dist.to_numpy()
    #%%
    #Create spline model
    
    cs = CubicSpline(dist,elev)
    
    syn_elev = cs(dist)
    
    #%%
    
    #Generating Plots
    
    if (plot == 1) or (plot == 3):
    
        plt.figure()
        plt.plot(dist, elev, linewidth=2, label='Elevation Profile')
        plt.plot(dist, syn_elev, linewidth=1, label='Spline Model')
        plt.xlabel('Distance (degrees)')
        plt.ylabel('Elevation (m)')
        plt.title(site_file_name)
        plt.legend(loc='best')
    
    #%%
    # Making it 3D
    
    y = dist
    x = y
    z = np.zeros((len(y),len(x)))
    
    #%%
    # Defining the shorline shoreline (z=0 contour)
    above_sl = elev>-5
    # above_sl = above_sl.to_numpy()
    base_shore = int(min(np.argwhere(above_sl)))
    y_bs = y[base_shore]
    x_max = max(x)
    y_max = max(y)
    
    # give shoreline shape
    if shape == 'straight':
        s = np.ones(len(y))*y_bs # flat shoreline
    if shape == 'curved':
        s = y_bs+2000*abs(np.sin((x)*np.pi/x_max)) # one big curve
    if shape == 'points':
        s = y_bs+2000*abs(np.sin((x+x_max/4)*2*np.pi/x_max)) #curved with points
    
    # add in triangular bay centered at x = x_bay
    if bay == True:
        sb = np.empty((len(y)))
        sb[:]=np.NaN
        x_bay = x_max/2
        bay_width = 4000
        bay_height = 8000
        out_bay = np.logical_or(x<(x_bay-bay_width/2),x>(x_bay+bay_width/2))
        in_bay = out_bay==False
        sb[out_bay]= np.NaN
        sb[in_bay] = s[int(len(s)/2+bay_width/2/3)]+bay_height-abs(x[in_bay]-x_bay)*\
            bay_height/(bay_width/2)
        
        # add in transition
        
        trans = np.empty((len(y)))
        trans[:]=np.NaN
        trans_width = 1000
        out_trans = np.logical_or(x<(x_bay-bay_width/2-trans_width),x>(x_bay+bay_width/2+trans_width))
        in_trans = out_trans==False
        trans[out_trans]= np.NaN
        trans[in_trans] = s[int(len(s)/2+(bay_width/2+trans_width)/3)]+bay_height-\
            abs(x[in_trans]-x_bay)*bay_height/(bay_width/2)+trans_width*bay_height*2/bay_width
        
        # add in river
        
        river_width = bay_width/10
        x_riv = np.logical_and(x > x_bay-river_width/2, x < x_bay+river_width/2)
    #%%

    # to make 2d grid of z values
    for col in range(z.shape[1]):

        ynew = y-s[col]+y_bs
        ynew[ynew<0]=0
        ynew[ynew>y_max]=y_max
        z[:,col]=cs(ynew)
        if (bay == True) and (in_bay[col]==True):
            y_bay = np.logical_and(y>s[col],y<sb[col])
            z[y_bay,col]=-5
            if x_riv[col]==True:
                y_riv = y>sb[col]
                z[y_riv,col]=-5
    
    #%%
    
#Add in transition
    if bay == True:
        for yind in range(z.shape[0]):
            if y[yind]>=min(trans[np.isnan(trans)==False]):
                if y[yind]<max(sb[np.isnan(sb)==False]):
                    
                    zcross = z[yind,:]
                    
                    left_trans = np.logical_and(x<x_bay,trans>=y[yind])
                    right_trans = np.logical_and(x>x_bay,trans>=y[yind])
                    
                    zcross[left_trans]=np.nan
                    zcross[right_trans]=np.nan
                    
                    xsub = x[np.isnan(zcross)==False]
                    zcross = zcross[np.isnan(zcross)==False]

                    cs2 = interp1d(xsub,zcross)
                    
                    zcross_syn = cs2(x)
                    
                    z[yind,:]=zcross_syn
                    
                if y[yind]>=max(sb[np.isnan(sb)==False]):
                    zcross = z[yind,:]
                    
                    left_trans = np.logical_and(x<x_bay,x>x_bay-river_width/2-trans_width)
                    right_trans = np.logical_and(x>x_bay,x<x_bay+river_width/2+trans_width)
                    
                    zcross[left_trans]=np.nan
                    zcross[right_trans]=np.nan
                    
                    xsub = x[np.isnan(zcross)==False]
                    zcross = zcross[np.isnan(zcross)==False]
                    
                    cs2 = interp1d(xsub,zcross)
                    
                    zcross_syn = cs2(x)
                    
                    z[yind,:]=zcross_syn
    
    #%%
    
    #heatmap
    if (plot == 2) or (plot == 3):
        xis = np.linspace(0,len(x),num = 1000,endpoint=False).astype(int)
        yis = np.linspace(0,len(y),num = 1000,endpoint=False).astype(int)
        zis = np.meshgrid(yis,xis,indexing = 'ij')
        zsub = z[zis]
        zsub = np.flip(zsub)
        plt.figure()
        ax = sns.heatmap(zsub,cmap='seismic',center=0)
    
    
    #%%

my_dir = 'C:\\Users\\jbarrett.carter\\OneDrive\\CUAHSI-SI\\Topography\\Profiles'
profiles = ['Atlantic_city_NJ.csv','Marley_beach_SC.csv','Melbourne_FL.csv',
            'Savannah_GA.csv','Shallotte_NC.csv']
generate_topo(site_file_name='Melbourne_FL.csv', directory = my_dir,shape = 'points', bay = True, plot = 3)
# generate_topo(site_file_name='Shallotte_profile.csv',shape = 'curved', bay = True, plot = 3)
for profile in profiles:
    generate_topo(site_file_name=profile, directory = my_dir,shape = 'straight', bay = True, plot = 3)
