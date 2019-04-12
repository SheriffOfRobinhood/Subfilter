#!/usr/bin/env python

import matplotlib
matplotlib.use('agg') 
# This means that can run this in script, but stops 
#    outputting plots through plt.show().
import matplotlib.pyplot as plt
from turbulent_fluxes import filter_2d, filtered_variable, \
     explicit_TKE, vertical_flux, delta_horizontal, level_heights,\
     filtered_strain_rate
from matplotlib import cm as cm
import numpy as np
import iris
import iris.quickplot as qplt	   #iris plotting module
#from mayavi import mlab

# Global constants
#===============================================================================
# Directory to which plots will be output.
plot_dir = '' 

plot_type = '.eps'
data_dir = '' # Directory containing data
figshow = True

def plot_q(q_r, q_s,  ilev, ix1, ix2, iy1, iy2, iy, j):
    meanfield= np.mean(q_r.data,axis=(1,2),keepdims=1)
    pltdat = (q_r.data-meanfield)*1000.

    lev1 = np.arange(-1,1.1,0.1)
    lev2 = np.arange(-1,1.1,0.1)

    plt.figure(1)
    plt.clf
    plt.subplot(2, 2, 1)
    Cs1 = plt.contourf(pltdat[ilev,iy1:iy2,ix1:ix2],levels=lev1)
    plt.title(r'$q^r$ pert level '+'%03d'%ilev)
#    plt.xlabel('word length anomaly')
#    plt.ylabel('sentence length anomaly')

# Make a colorbar for the ContourSet returned by the contourf call.
    cbar1 = plt.colorbar(Cs1)
    cbar1.ax.set_ylabel('q')
# Add the contour line levels to the colorbar
#  cbar.add_lines(CS2)

    plt.subplot(2, 2, 2)
    Cs2 = plt.contourf(q_s.data[ilev,iy1:iy2,ix1:ix2]*1000.,levels=lev2)
    plt.title(r'$q^s$ level '+'%03d'%ilev)
    cbar2 = plt.colorbar(Cs2)
    cbar2.ax.set_ylabel('q')


    plt.subplot(2, 2, 3)
    Cs3 = plt.contourf(q_r.coord('longitude').points[ix1:ix2], \
      q_r.coord('level_height').points, \
      pltdat[:,iy,ix1:ix2],levels=lev1)

    plt.ylim([0,5000])
    plt.title(r'$q^r$ pert at iy '+'%03d'%iy)
# Make a colorbar for the ContourSet returned by the contourf call.
    cbar3 = plt.colorbar(Cs3)
    cbar3.ax.set_ylabel('q')

    plt.subplot(2, 2, 4)
    Cs4 = plt.contourf(q_s.coord('longitude').points[ix1:ix2], \
      q_s.coord('level_height').points, \
      q_s.data[:,iy,ix1:ix2]*1000.,levels=lev2)

    plt.ylim([0,5000])
    plt.title('$q^s$ x-z at iy '+'%03d'%iy)
# Make a colorbar for the ContourSet returned by the contourf call.
    cbar4 = plt.colorbar(Cs4)
    cbar4.ax.set_ylabel('q')
    plt.savefig('q_lev_'+'%03d'%ilev+'_x_z'+'%03d'%iy+'_%02d'%j+'_20km.png')

    plt.close()
#    plt.show()

    return

def test_filter(pp_file1, pp_file2, pp_file3, twod_filter, contents=0):
#   Work out how many times in data
    n_times = len(pp_file1)
#   Set up list of times
    times = []
   
    for j in range(n_times):
#       Calculate time for use in title

        if contents == 1 :
            cube = iris.load(data_dir+pp_file2[j])
            for c in cube :
                print c.attributes['STASH'], c.coord('time')
            cube=cube[0]
            time = cube.coord('time')
            time = str(time.units.num2date(time.points[0]).time())[0:5]
            print time
            times = times + [time]

        q = filtered_variable('specific humidity', pp_file2[j], twod_filter,\
            mean_profile=0)

        print 'q[0]\n', q[0]
        print 'q[1]\n', q[1]
        q_r = q[0]
        q_s = q[1]
        ilev = 10

        iy = 512

        ix1 = 400
        ix2 = 600
        iy1 = 400
        iy2 = 600

        plot_q(q_r, q_s,  ilev, ix1, ix2, iy1, iy2, iy, j)    

    return q

def main():
    '''
    Top level code, a bit of a mess.
    '''
#   Non-global variables that are set once
#    sigma_list = [0.2,0.25,0.3,0.4,0.8,1.0,1.5,2.0] # Sigma used in Gaussian filter function
#    sigma_list = [2.0] # Sigma used in Gaussian#
    sigma_list = [0.5]
    filter_name = 'gaussian'
    n_hours = 12
#    n_hours = 1
    

    pp_file1a = ['anyha1_' + str(i).zfill(2) + '.pp' for i in range(n_hours)]
    pp_file2a = ['anyha2_' + str(i).zfill(2) + '.pp' for i in range(n_hours)]
    pp_file3a = ['anyha3_' + str(i).zfill(2) + '.pp' for i in range(n_hours)]

    pp_file1b = ['anyhb1_' + str(i).zfill(2) + '.pp' for i in range(n_hours)]
    pp_file2b = ['anyhb2_' + str(i).zfill(2) + '.pp' for i in range(n_hours)]
    pp_file3b = ['anyhb3_' + str(i).zfill(2) + '.pp' for i in range(n_hours)]

# Just looking at anyhb at present

    for sigma in sigma_list:
        twod_filter = filter_2d(filter_name, sigma=sigma, delta_x=0.1)
        res = test_filter(pp_file1a[0:], pp_file2a[0:], pp_file3a[0:], twod_filter)

if __name__ == "__main__":
    main()
