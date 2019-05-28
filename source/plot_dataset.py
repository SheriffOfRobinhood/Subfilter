# -*- coding: utf-8 -*-
"""

  plot_dataset.py

Created on Tue Oct 23 09:52:58 2018

@author: Peter Clark
"""
import os 
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from netCDF4 import Dataset

import subfilter

def print_ncattr(nc_fid,key):
    """
    Prints the NetCDF file attributes for a given key

    Parameters
    ----------
    key : unicode
        a valid netCDF4.Dataset.variables key
    """
    try:
        print "\t\ttype:", repr(nc_fid.variables[key].dtype)
        for ncattr in nc_fid.variables[key].ncattrs():
            print '\t\t%s:' % ncattr,\
                  repr(nc_fid.variables[key].getncattr(ncattr))
    except KeyError:
        print "\t\tWARNING: %s does not contain variable attributes" % key

def ncattrprint(nc_fid) :
    # Attribute information
    nc_attrs = nc_fid.ncattrs()
    print "NetCDF Global Attributes:"
    for nc_attr in nc_attrs:
        print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))


def ncdimprint(nc_fid) :
   # Dimension shape information.
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    print "NetCDF dimension information:"
    for dim in nc_dims:
        print "\tName:", dim 
        print "\t\tsize:", len(nc_fid.dimensions[dim])
        print_ncattr(nc_fid,dim)

def ncvarprint(nc_fid) :
    # Variable information.
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    print "NetCDF variable information:"
    for var in nc_vars:
        if var not in nc_dims :
 #       if (var not in nc_dims) and  (len(nc_fid.variables[var].dimensions) == 4) :
            print '\tName:', var
            print "\t\tdimensions:", nc_fid.variables[var].dimensions
            print "\t\tsize:", nc_fid.variables[var].size
            print_ncattr(nc_fid,var)


def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        ncattrprint(nc_fid)
        ncdimprint(nc_fid)
        ncvarprint(nc_fid)        
    return nc_attrs, nc_dims, nc_vars

#dir = 'C:/Users/paclk/OneDrive - University of Reading/python/SG/'
#file = 'diagnostics_ts_14400.0.nc'
#dir = 'C:/Users/paclk/OneDrive - University of Reading/python/subfilter/test_data/BOMEX/'
#file = 'diagnostics_ts_18000.0.nc'
#file = 'diagnostics_ts_18000.0_filter_00.nc'

dir = '/gws/nopw/j04/paracon_rdg/users/toddj/MONC_RCE_1.5/MONC_RCE_1.5_m1600_g0084/diagnostic_files/'
file = 'diagnostics_ts_4784400.0.nc'

dataset = Dataset(dir+file, 'r') # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class

nc_attrs, nc_dims, nc_vars = ncdump(dataset,verb=True)

#ncdimprint(dataset)
#ncvarprint(dataset)

#w = dataset.variables['w_r']
w = dataset.variables['w']

w_tvar = w.dimensions[0]
w_time = dataset.variables[w_tvar]
print w.shape
print w_tvar
print w_time.dtype

#iz = 20
#for i in range(0,w.shape[0]):
#    fig = plt.figure(i)
#    plt.contourf(w[i,:,:,iz],20)
#    
#plt.show()


    

#for variable in dataset.variables :
#    print(variable, variable.dimension())
    #vars(variable)
 #   print variable
#    print dataset.variables[variable].dimensions
