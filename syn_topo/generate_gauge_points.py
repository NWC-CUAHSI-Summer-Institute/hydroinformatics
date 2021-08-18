# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:01:00 2021

@author: J. Barrett Carter

This is a function that will be used to automatically generate gauge locations
for the synthetic topography from a specified elevation profile.

The 'site_file_name' variable must be a csv file containing a column for distance
labeled 'Distance' and a column for elevation labeled 'Elevation'. The distance
should be measured from a point on land to a point at sea, and the 'Distance' 
values should be orded from smallest to largest.

The options for the shape parameter are either 'straight', 'curved', or 'points'. 
The options for the bay parameter are True or False. The options for the plot 
parameter are True or False, which determines if plot are generated.

return_gauge_points is either True or False and determines if the xy coordinates
for the gauges are returned by the function.

dim is the dimension of the coarse z-array in km. This must match that used in
the 'generate_syn_topo' function.

"""
def generate_gauge_points(site_file_name,directory,shape = 'straight',bay = False,
                          plot = True,return_gauge_points = True, dim = 5000):

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
    file2 = "OceanDepth.csv"
    filepath2 = os.path.join(datadir, file2)
    filepath = os.path.join(datadir, filename)
    df = pd.read_csv(filepath)
    df2 = pd.read_csv(filepath2)
    df2 = df2.loc[0:16,:]
    
    #%%
    ### Generate fine topography
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
    
    # #Generating Plots
    
    # # if (plot == 1) or (plot == 3):
    # if (plot == True):
    
    #     plt.figure()
    #     plt.plot(dist, elev, linewidth=2, label='Elevation Profile')
    #     plt.plot(dist, syn_elev, linewidth=1, label='Spline Model')
    #     plt.xlabel('Distance (degrees)')
    #     plt.ylabel('Elevation (m)')
    #     plt.title(site_file_name)
    #     plt.legend(loc='best')
    
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
    # if (plot == 2) or (plot == 3):
    if (plot == True):
        xis = np.linspace(0,len(x),num = 1000,endpoint=False).astype(int)
        yis = np.linspace(0,len(y),num = 1000,endpoint=False).astype(int)
        zis = np.meshgrid(yis,xis,indexing = 'ij')
        zsub = z[zis]
        zsub = np.flip(zsub)
        plt.figure()
        sns.heatmap(zsub,cmap='seismic',center=0)
        
    #%%

    # define gauge points
    
    x_max_ggs = x_max/2-3000
    x_min_ggs = x_max_ggs-12*1000
    x_ggs = np.linspace(x_min_ggs,x_max_ggs,num = 10)
    
    gauge_points = np.zeros((1,2))
    x_gg = x_ggs[1]
    for x_gg in x_ggs:
        x_ind = int(x_gg/3)
        zs = z[:,x_ind]
        y_shore_ggs = y[np.logical_and(zs>-0.5,zs<0)]
        # y_shore_ggs = y[np.round(zs)==0]
        y_sg_ind = y_shore_ggs/3
        y_sg_ind = y_sg_ind.astype(int)
        y_rem = []
        for i in range(len(y_sg_ind)-1):
            if y_sg_ind[i+1]==y_sg_ind[i]+1:
                y_rem=np.append(y_rem,i)
        y_sg_ind = np.delete(y_sg_ind,y_rem.astype(int))
        y_shore_ggs = y[y_sg_ind]
        y_ind_new = []
        for i in range(len(y_sg_ind)-1):
            y_ind_new = np.append(y_ind_new,int((y_sg_ind[i]+y_sg_ind[i+1])/2))
        y_ind = np.append(y_sg_ind,y_ind_new).astype(int)
        y_ggs = (dim-30)*1000+y[y_ind]
        for y_gg in y_ggs:    
            gauge_points = np.append(gauge_points,[[x_gg,y_gg]],axis = 0)
    
    gauge_points = np.delete(gauge_points,[0],axis = 0)
    
    
    if bay == True:
        x_bay_ggs = x_bay    
        y_bay_ggs = np.linspace((dim-30)*1000+y_bs,(dim-30)*1000+y_max,10)
        
        for y_b_gg in y_bay_ggs:
            gauge_points = np.append(gauge_points,[[x_bay_ggs,y_b_gg]],axis = 0)
    
    if (plot == True):
        plt.figure()
        plt.plot((dim-30)*1000+y,zs)
        # plt.plot(y_shore_ggs,zs[y_sg_ind],'r.')
        plt.title(site_file_name)
        plt.plot(y_ggs,zs[y_ind],'r.')
    
    
    #%%
    
    # plot gauge locations
    if (plot == True):
        plt.figure()
        plt.plot(gauge_points[:,0],gauge_points[:,1],'g.')
    
    #%%
    
    if return_gauge_points == True:
        return(gauge_points)
    
    #%%

my_dir = 'C:\\Users\\jbarrett.carter\\OneDrive\\CUAHSI-SI\\Topography\\Profiles'

### to generate one set of gauges at a time
# ggs = generate_gauge_points(site_file_name='Melbourne_FL.csv', directory = my_dir,
#                             shape = 'straight', bay = False, plot = True)

### to generate all sets of gauge locations

shapes = ['straight','points','curved']

bays = [True, False]

profiles = ['Atlantic_city_NJ.csv','Marley_beach_SC.csv','Melbourne_FL.csv',
            'Savannah_GA.csv','Shallotte_NC.csv']

gauge_point_sets = {}

for p in profiles:
    name0 = p.split(sep = '.')[0]
    for b in bays:
        name1 = str(b)
        for s in shapes:
            name2 = s
            name = name0+'-'+name1+'-'+name2
            gauge_point_sets[name]=\
                generate_gauge_points(site_file_name=p,
                                      directory = my_dir,shape = s, 
                                      bay = b, plot = True)
            
            


# generate_topo(site_file_name='Shallotte_profile.csv',shape = 'curved', bay = True, plot = 3)
# for profile in profiles:
#     generate_topo(site_file_name=profile, directory = my_dir,shape = 'straight', 
#                   bay = True, plot = 3)
# topo_fine = zs[0]
# topo_coarse = zs[1]

### experimenting

# import numpy as np
# my_dict = {}
# my_array = np.array([[1,2,3],[4,5,6]])
# my_array2 = np.array([[10,11,12,13],[14,15,16,17]])
# my_var = 'C'
# my_array3 = np.array([[20,30],[40,50]])
# my_dict['A']=my_array
# my_dict['B']=my_array2
# my_dict[my_var]=my_array3
