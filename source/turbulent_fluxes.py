'''
Calculate TKE, stress and buoyancy flux terms from UM output

This includes all the code I've written so far as part of my secondment, looking at turbulent mixing in idealised UM simulations. At some point I shall probably split this into different files and try to generalise some of the code. This includes:
- Explicit TKE calculation and calculation using the strain rate and buoyancy gradient. 
- Three different 2D filters (and 1D versions). Now moved to filters.py
- Save calculated TKEs to netCDF files, and check whether TKE is already in pickle file before recalculating (Explicit TKE calculations at high res take a while)
- Plotting programs to show calculated TKEs.
'''

import numpy as np
import iris
from scipy.signal import fftconvolve
from sys import float_info
from time import strftime
from os import listdir
from copy import copy
import filters # user written module in /home/mm0100/phill/python/
import velocity_gradients
#import scalar_gradients

# Global constants
#===============================================================================
eps = float_info.min # smallest possible float
pi = 3.141592653589793
earth_radius = 6371229.0 # from earth_constants_mod
data_dir = '' # Directory containing data
filter_fileroot = 'filter_fns'
explicit_TKE_fileroot = 'explicit_TKE'
Smag_TKE_fileroot = 'Smag_TKE'
Smag_buoyancy_flux_fileroot = 'Smag_buoyancy_flux'
diagnostic_buoyancy_flux_fileroot = 'diag_buoyancy_flux'
Smag_humidity_flux_fileroot = 'Smag_humidity_flux'
Smag_dev_stress_fileroot = 'Smag_dev_stress'
dev_stress_fileroot = 'dev_stress'
vertical_flux_fileroot = 'vertical_flux'
hist2d_fileroot = 'hist2d'
filtered_vars_fileroot = 'filtered_vars'
vertical_gradient_fileroot = 'vertical_grad'
filtered_strain_rate_fileroot = 'filtered_strain_rate'
derived_data_dir = data_dir + '.'
# dictionary links my variable names to stash codes
stash_dict = {'strain rate':'13192', \
              'buoyancy gradient':'3469', \
              'liquid water':'254', \
              'frozen water':'12', \
              'specific humidity':'10', \
              'buoyancy':'-9999', \
              'potential temperature':'4', \
              'liquid water potential temperature':'theta_L', \
              'vgrad liquid water potential temperature':'dtheta_Ldz', \
              'virtual potential temperature':'-8888', \
              'vertical velocity':'150', \
              'u wind':'2', \
              'v wind':'2', \
              'visc_h':'13191', \
              'visc_m':'13190'} 
cp = 1005.
lc = 2.501E6

#===============================================================================


def read_variable(pp_file, code, time=None):
    '''
    Reads variable identified by the given stash code from the given pp file.

    Args:
      pp_file (str): Name of the pp file that will be read
      code (str): Up to five digit stash code of the variable requested.
      time (str, default=None): point in time required. Functionality not
                                coded yet

    Returns:
      cubes (list): list of iris cubes contained in the pp file that match the
                    criteria.
    '''
    if code == 'theta_L':
        cubes = calc_theta_l(pp_file)
    elif code == 'dtheta_Ldz':
        cubes = calc_theta_l(pp_file)
        heights = cubes[0].coords('level_height')[0].points
        for cube in cubes :
            d_theta_l_dz = velocity_gradients.dw_dz(heights, heights, cube.data, surface_value='same')
            print np.shape( cube.data), np.shape(d_theta_l_dz)
            cube.data = d_theta_l_dz
#            cube.data[-2:-1,:,:] = 0.0
            cube.attributes['var_name'] = 'vgrad liquid water potential temperature'
        
    elif code == '-8888':
        cubes = calc_theta_v(pp_file)
    elif code == '-9999':
# This code corresponds to the buoyancy, which must be calculated from other variables rather than read directly.
        cubes = calc_theta_v(pp_file)[0]
        cubes_mean = cubes.data.mean(axis=2,dtype=np.float64).mean(axis=1,dtype=np.float64)
        cubes.data = cubes.data - cubes_mean[:, np.newaxis, np.newaxis]
        cubes = [cubes]
    else:
        stash_code = iris_stash_code(code)
        if time == None: # Read all times in the file
            cubes = iris.load(data_dir + pp_file,
                              iris.AttributeConstraint(STASH = stash_code))
        else: # Read particular time requested
            print 'This has not been coded yet'
            cubes = -1
    return cubes


def iris_stash_code(code):
    '''
    Converts stash code to format required by iris

    Args:
      code (str): Stash code string of up to 5 digits

    Returns:
      iris_stash_code (str): Stash code in format required by iris.
    '''
    temp = "%05d" % int(code)
    iris_stash_code = 'm01s'+temp[0:2]+'i'+temp[2:]
    return iris_stash_code


def Smag_TKE(pp_file1, pp_file2, mean_profile=1):
    '''
    Cube contains TKE calculated from strain rate and buoyancy gradient.

    Uses turblent kinetic energy = lambda^2 * (strain rate^2
                                               - buoyancy gradient).
    If buoyancy gradient > strain rate ^2, defaults to zero.
    Before it is caculated, a netCDF file is read. If the TKE has already been
    calculated, it will be stored in the netCDF file and will be read instead
    of recalculated. Once calculated, the TKE is stored in a netCDF file.

    Args:
      pp_file1 (str): Name of pp file containing subgrid turbulence scheme
                      diagnostics (i.e. strain rate and lambda)
      pp_file2 (str): Name of pp_file containing boundary layer diagnostics
                     (i.e. buoyancy gradient).
      mean_profile (bool, default=0): Only store the mean profile in the
                                          netCDF file (saves a lot of disk
                                          space).

    Attributes:
      pp_file1 (str): As for arg above.
      pp_file2 (str): As for arg above.
      mean_profile (bool): As for arg above.
      level_heights (nd_array): Array of floats, describes height of model 
                                levels in metres
      data (nd_array): Array of floats, which give the magnitude of the TKE. If
                       mean_profile=1, this is the mean profile, if not it is
                       the full 3D array.
    Dimensions:
 
    '''
    attributes = {'pp_file1':pp_file1, 
                  'pp_file2':pp_file2,
                  'mean_profile':mean_profile}
    new_calc, result = read_cubes(Smag_TKE_fileroot, attributes) # Read existing TKEs (if any exist)
    if new_calc == 1: # If data does not already exists, calculate now.
        heights = level_heights(pp_file2, '3469')[0:-1] #NB This needs checking
        buoyancy_grad = read_variable(pp_file2, '3469')[0].data[1:,:,:] 
# Based on their level_height coordinates, the buoyancy_grad and strain rate 
# variables are defined on the same levels. 
# However, comparing the BL shear diagnostic to the strain rate suggests that 
# they are offset by a level, but I'm not certain about this. 
# Moreover this offset gives me more realistic results. 
# Good place to check if I have problems though.
        my_lambda = read_variable(pp_file1, '13193')[0].data[:-1,:,:]#NB This needs checking
        strain_rate = read_variable(pp_file1, '13192')[0].data[:-1,:,:]/np.sqrt(2.0)#NB This needs checking
        Pr = 1.0
        TKE = (my_lambda ** 2) * ( (strain_rate ** 2) - (Pr * buoyancy_grad) )
        TKE[ TKE < 0.0 ] = 0.0 # Lots of negative TKE due to positive buoyancy.
        if mean_profile == 1:
            TKE = calc_mean_profile(TKE)
# Write this new filter to the file
        result = write_cubes(Smag_TKE_fileroot, TKE, attributes, heights=heights)
    return result


def explicit_TKE(pp_file, twod_filter, mean_profile=1):
    '''
    Class contains resolved TKE calculated from filtered velocity fields

    Uses turblent kinetic energy = 0.5 * (u'**2 + v'**2 + w'**2)

    Before it is caculated, a netCDF file is read. If the TKE has already been
    calculated, it will be stored in the netCDF file and will be read instead
    of recalculated. Once calculated, the TKE is stored in a netCDF file.

    Args:
      pp_file (str): Name of pp file containing velocity diagnostics.
      twod_filter(class, filter_2d): Class containing data and metadata for the 
                                     2D filter that is used to calculate the
                                     velocity perturbations used in the TKE
                                     calculation
      mean_profile (bool, default=0): Only store the mean profile in the
                                          netCDF file (saves a lot of disk
                                          space).

    Attributes:
      filter_name (str): Name of filter used. Either Gaussian, wave-cutoff or 
                         running-mean.
      wavenumber (float): If a wave-cutoff filter is used, contains the cutoff 
                          wavenumber.
      delta_x (float): Distance between points in the horizontal, 
                       used to caculate the filter
      width (int): If set, controls the width of the filter. Must be set for
                   running-mean filter.
      cutoff (float): If float is not set, this controls the width of the
                      filter. The width of the filter is extended until the
                      minimum value in the filter is less than this cutoff
                      value.
      high_pass (bool): If a wave-cutoff filter is used, this determines whether
                        it is high or low pass (note high pass hasn't actually
                        been coded yet!)
      sigma (float): If a Gaussian filter is used, this is the lengthscale of
                     the filter.
      pp_file (str): As for arg above.
      mean_profile (bool): As for arg above.
      level_heights (nd_array): Array of floats, describes height of model 
                                levels in metres
      data (nd_array): Array of floats, which give the magnitude of the TKE. If
                       mean_profile=1, this is the mean profile, if not it is
                       the full 3D array. 
    '''
    attributes = {'var_name' : '', \
                  'filter_name' : twod_filter.attributes['filter_name'],
                  'wavenumber' : twod_filter.attributes['wavenumber'],
                  'delta_x' : twod_filter.attributes['delta_x'],
                  'width' : twod_filter.attributes['width'],
                  'cutoff' : twod_filter.attributes['cutoff'],
                  'high_pass' : twod_filter.attributes['high_pass'],
                  'sigma' : twod_filter.attributes['sigma'],
                  'pp_file' : pp_file,
                  'mean_profile' : mean_profile}
# Read existing explicit TKEs (if any exist)

    new_calc, result = read_cubes(explicit_TKE_fileroot, attributes, any_var = 1)
    if type(result) is int :
        print 'result is an integer'
	new_calc = 1
    elif len(result) != 4 : 
        print 'result is not a list of length 4'
	new_calc = 1

    if new_calc == 1:# If data does not already exists, calculate now.
        heights = level_heights(pp_file, '150')
        regrid_cube = read_variable(pp_file, '150')[0]
        print 'Computing u_part of TKE'
        u_part = quadratic_subfilter(pp_file, pp_file, '2', '2', \
                                     twod_filter.data, mean_profile, \
                                     levels=heights, regrid_cube=regrid_cube)
        print 'Computing v_part of TKE'
        v_part = quadratic_subfilter(pp_file, pp_file, '3', '3',\
                                     twod_filter.data, mean_profile, \
                                     levels=heights, regrid_cube=regrid_cube)
        print 'Computing w_part of TKE'
        w_part = quadratic_subfilter(pp_file, pp_file, '150', '150', \
                                     twod_filter.data, mean_profile)
        print 'Computing TKE'
        TKE = 0.5 * (u_part[0] + v_part[0] + w_part[0])
# Write this variable to file
        result = write_cubes(explicit_TKE_fileroot, \
                 [TKE, u_part[0], v_part[0], w_part[0]], attributes, \
                 heights=heights, set_var_name = ['TKE', 'TKE_u', 'TKE_v', 'TKE_w'])

    return result


def total_TKE(Smag_TKE, explicit_TKE):
    '''
    Cube containing the total (Smagorinsky plus explicit) TKE
    '''
    new_explicit_TKE = explicit_TKE.interpolate([('level_height', Smag_TKE.coord('level_height').points)], iris.analysis.Linear(extrapolation_mode='nan'))
    new_explicit_TKE.coord('model_level_number').points = copy(
        Smag_TKE.coord('model_level_number').points)
    total_TKE = iris.analysis.maths.add(Smag_TKE, new_explicit_TKE)
    return total_TKE
    


def filtered_variable(var_name, pp_file, twod_filter, mean_profile=1):
    '''
    Class contains resolved variables calculated by filtering the full field

    Before it is caculated, a netCDF file is read. If the filtered variable has
    already been calculated, it will be stored in the netCDF file and will be 
    read instead of recalculated. Once calculated, the filtered variable is 
    stored in a netCDF file.

    Args:
      var_name (str): Name of variable to be filtered. In order for this code
                      to work, this must be linked to a stash code in the 
                      stash_dict dictionary.
      pp_file (str): Name of pp file containing variable.
      twod_filter(class, filter_2d): Class containing data and metadata for the 
                                     2D filter that is used to calculate the
                                     velocity perturbations used in the TKE
                                     calculation
      mean_profile (bool, default=0): Only store the mean profile in the
                                          netCDF file (saves a lot of disk
                                          space).

    Attributes:
      var_name (str): As for arg above.
      filter_name (str): Name of filter used. Either Gaussian, wave-cutoff or 
                         running-mean.
      wavenumber (float): If a wave-cutoff filter is used, contains the cutoff 
                          wavenumber.
      delta_x (float): Distance between points in the horizontal, 
                       used to caculate the filter
      width (int): If set, controls the width of the filter. Must be set for
                   running-mean filter.
      cutoff (float): If float is not set, this controls the width of the
                      filter. The width of the filter is extended until the
                      minimum value in the filter is less than this cutoff
                      value.
      high_pass (bool): If a wave-cutoff filter is used, this determines whether
                        it is high or low pass (note high pass hasn't actually
                        been coded yet!)
      sigma (float): If a Gaussian filter is used, this is the lengthscale of
                     the filter.
      pp_file (str): As for arg above.
      mean_profile (bool): As for arg above.
      level_heights (nd_array): Array of floats, describes height of model 
                                levels in metres
      data (nd_array): Array of floats, which give the magnitude of the TKE. If
                       mean_profile=1, this is the mean profile, if not it is
                       the full 3D array. 
    '''
    print 'var_name =', var_name

    print 'twod_filter\n',twod_filter

    attributes = {'var_name' : var_name, \
                  'filter_name' : twod_filter.attributes['filter_name'], \
                  'wavenumber' : twod_filter.attributes['wavenumber'], \
                  'delta_x' : twod_filter.attributes['delta_x'], \
                  'width' : twod_filter.attributes['width'], \
                  'cutoff' : twod_filter.attributes['cutoff'], \
                  'high_pass' : twod_filter.attributes['high_pass'], \
                  'sigma' : twod_filter.attributes['sigma'], \
                  'pp_file' : pp_file, \
                  'sub_filter' : 0, \
                  'mean_profile' : mean_profile}

    print 'attributes\n' , attributes

# Read existing filtered variables
    new_calc, result = read_cubes(filtered_vars_fileroot, attributes)

    print 'result\n', result

    if type(result) is int :
        print 'result is an integer'
	new_calc = 1
    elif len(result) != 2 : 
        print 'result is not a list of length 2'
	new_calc = 1

    if new_calc == 1:# If data does not already exists, calculate now.
        heights = level_heights(pp_file, stash_dict[var_name])
        data = calc_filtered_field(pp_file, stash_dict[var_name],
                                   twod_filter.data)#[0]

        print 'filtered_variable data type:',type(data)
# Write this new filter to the file
        result = write_cubes(filtered_vars_fileroot, data, attributes, \
                             heights=heights, set_sub_filter=[0,1])
        print 'result\n', result

    
    return result
        

def calc_filtered_field(pp_file, stash_code, twod_filter, levels=None, regrid_cube=None):
    '''
    Filters a field and returns the resolved and subfilter parts

    Args:
      pp_file (str): Name of pp_file containing variable.
      stash_code (str): Stash code for variable, in 5 digit format.
      twod_filter (nd_array):: Kernel for the filtering.

    Returns:
      field_r (nd_array): resolved part of the filtered field.
      field_s (nd_array): Subfilter part of the filtered field.    
    '''
    if levels == None:
        field = read_variable(pp_file, stash_code)[0]
    else :
        field = read_variable(pp_file, stash_code)[0].interpolate([('level_height', levels)], iris.analysis.Linear(extrapolation_mode='nan'))
    if regrid_cube != None:
        field = periodic_regrid(field, regrid_cube)
    field = field.data
    field_r, field_s = filtered_field_calc(field, twod_filter)
    print 'Result from calc_filtered_field', type(field_r), type(field_s)
    return [field_r, field_s]
    

def filtered_field_calc(field, twod_filter):
    field_r = np.zeros(field.size).reshape(field.shape)
    for i in range(field.shape[0]):
        field_r[i,:,:] = convolve(field[i,:,:], twod_filter)
    field_s = field - field_r
    return [field_r, field_s]


def filter_2d(filter_name, wavenumber=-1, delta_x=1.0, width=-1,
                 cutoff=0.0001, high_pass=0, sigma=-1):
    '''
    Class to link a filter to it's properties

    Args:
      filter_name (str): Name of filter used. Either Gaussian, wave-cutoff or 
                         running-mean.
      wavenumber (float): If a wave-cutoff filter is used, contains the cutoff 
                          wavenumber.
      delta_x (float): Distance between points in the horizontal, 
                       used to caculate the filter
      width (int): If set, controls the width of the filter. Must be set for
                   running-mean filter.
      cutoff (float): If float is not set, this controls the width of the
                      filter. The width of the filter is extended until the
                      minimum value in the filter is less than this cutoff
                      value.
      high_pass (bool): If a wave-cutoff filter is used, this determines whether
                        it is high or low pass (note high pass hasn't actually
                        been coded yet!)
      sigma (float): If a Gaussian filter is used, this is the lengthscale of
                     the filter.

    Args:
      filter_name (str): As arg above.
      wavenumber (float): As arg above.
      delta_x (float): As arg above.
      width (int): As arg above.
      cutoff (float): As arg above.
      high_pass (bool): As arg above.
      sigma (float): As arg above.
      data (nd_array): Array of filter values.
    '''
    attributes = {'filter_name' : filter_name,
                  'wavenumber' : wavenumber,
                  'delta_x' : delta_x,
                  'width' : width,
                  'cutoff' : cutoff,
                  'high_pass' : high_pass,
                  'sigma' : sigma}
# Read existing filters (if any exist)
    new_calc, result = read_cubes(filter_fileroot, attributes)
# If data does not already exists, calculate now.
    if new_calc:
        print 'No existing filter, so calculating now'
        if (filter_name == 'gaussian'):
            if (sigma == -1):
                data = filter_2d_error(filter_name, 'sigma')
            else:
                 data = filters.gaussian_filter_2d(sigma, delta_x, cutoff, width)
        elif (filter_name == 'running_mean'):
            if (width == -1):
                data = filter_2d_error(filter_name, 'width')
            else:
                data = filters.running_mean_filter_2d(width)
        elif (filter_name == 'wave_cutoff'):
            if (wavenumber == -1):
                data = filter_2d_error(filter_name, 'wavenumber')
            else:
                data = filters.wave_cutoff_filter_2d(wavenumber, delta_x, width,
                                             cutoff, high_pass)     
        else:
            print 'This filter type is not available.'
            print 'Available filters are:'
            print 'gaussian, running_mean & wave_cutoff'
            data = -9999
# Write this new filter to the file
        result = write_cubes(filter_fileroot, data, attributes)

    return result[0]


def read_cubes(fileroot, attributes, any_var=None):
    '''
    Reads objects from a netCDF file to a list.

    Filtering variables for high resolution simulations is very slow, so once 
    they're calculated they're netCDFd. This function reads the netCDFd files.
    It may be better to write to pp_files, but have yet to succesfully read a 
    cube that I've written to a pp_file.

    Args:
      fileroot (str): Name of netCDF file

    Returns:
      recalculate (int) :
      objs (list): List of objects contained in the netCDF file.
    '''
    recalculate = 1
    objs = - 9999
    constraint = iris.AttributeConstraint()
    for attribute in attributes:
        if attribute != 'var_name' or any_var != 1  :
            print 'Reading attribute', attribute
            if attribute != 'sub_filter' :
                constraint = constraint & eval('iris.AttributeConstraint('\
                 +attribute+' = attributes["'+attribute+'"])')
    for file in listdir(derived_data_dir):
        if file.startswith(fileroot) and file.endswith('.nc'):
            objs = iris.load(derived_data_dir + file, constraint)
            if len(objs) > 0:
#                objs = objs[0]
                recalculate = 0
                print '********************************************************'
                if 'var_name' in attributes.keys():
                    print 'Variable:', attributes['var_name']
                    print 'There are ', len(objs), ' cubes in the list.' 
                print 'Cube was read from ', file
                print '********************************************************'
                break
    print 'recalculate=', recalculate
    print "'mean_profile' in attributes.keys() =", 'mean_profile' in attributes.keys()
    if (recalculate == 1 and 'mean_profile' in attributes.keys()):
        if (attributes['mean_profile'] == 1):
            attributes_temp = copy(attributes)
            attributes_temp['mean_profile'] = 0
            recalculate, objs_temp = read_cubes(fileroot, attributes_temp, any_var=any_var)
            if recalculate == 0:
                heights = objs_temp.coords('level_height')[0].points
                data = objs_temp.data.mean(axis=2,dtype=np.float64).mean(axis=1,dtype=np.float64)
                data = data.reshape(len(data), 1, 1)
                objs = write_cubes(fileroot, data, attributes, heights=heights)
    return recalculate, objs


def write_cubes(fileroot, data_list, attributes, heights=None, \
                set_sub_filter = None, set_var_name = None):
    '''
    Writes list of objects to a netCDF file

    Filtering variables for high resolution simulations is very slow, so once 
    they're calculated they're netCDFd. This function reads the netCDFd files.
    It may be better to write to pp_files, but have yet to succesfully read a 
    cube that I've written to a pp_file.

    Args:
      filename (str): Name of netCDF file
      data (ndarray): data to be written
      heights (ndarray): height of levels, if not None, arbitrary values
                         are used for other dimensions.
    '''

    filename = derived_data_dir + fileroot+'_'+strftime("%d%m%Y_%H%M%S")+'.nc'

    objs = list([])
    if not isinstance(data_list, list) : 
        data_list = [data_list]

    obj_counter = 0
    for data in data_list : 
        attributes_temp = copy(attributes)
    
        if isinstance(set_sub_filter, list) :
            attributes_temp['sub_filter'] = set_sub_filter[obj_counter]

        if isinstance(set_var_name, list) :
            attributes_temp['var_name'] = set_var_name[obj_counter]

        if heights == None:
            obj = iris.cube.Cube(data, attributes=attributes_temp)
        else:
            heights_dim = iris.coords.AuxCoord(heights, long_name='level_height')
            model_level = iris.coords.DimCoord(np.arange(data.shape[0]), 
                            standard_name='model_level_number', units='1')
            latitude = iris.coords.DimCoord(np.arange(data.shape[1]), 
                            standard_name='latitude', units='degrees')
            longitude = iris.coords.DimCoord(np.arange(data.shape[2]), 
                                standard_name='longitude', units='degrees')
            obj = iris.cube.Cube(data, dim_coords_and_dims=[(model_level,0),
                                   (latitude,1), (longitude,2)], 
                                    aux_coords_and_dims=[(heights_dim,0)], 
                                    attributes=attributes_temp)
        objs.append(obj) 
        obj_counter += 1 

    print 'Saving to ',filename,  obj      
    iris.save(objs, filename)

    cubes = iris.load(filename)
    print 'Cubes in ' , filename
    for c in cubes :
        print c
    return objs


def filter_2d_error(filter_name, problem):
    '''
    Prints error when parameter required by filter does not exist.

    Args:
      filter_name (str): Name of filter
      problem (str): Name of parameter that has not been set

    Returns:
      filter_2d (-9999): Error code for filter.
    '''
    print 'A ' + filter_name + ' filter was selcted, but a suitable value'
    print 'for the ' + problem + ' was not chosen'
    filter_2d = -9999
    return filter_2d


def level_heights(pp_file, stash_code):
    '''
    Returns the values of the vertical axis (height in metres)

    Args:
      pp_file (str): Name of pp_file
      stash_code (str): Stash code for variable

    Returns:
      heights (nd_array): Height of model levels on which variable is defined.
    '''
    cube = read_variable(pp_file, stash_code)
    heights = cube[0].coords('level_height')[0].points
    return heights


def get_level_height(cube):
    level_height  = [coord.points for coord in cube.aux_coords if 
                     coord.long_name == 'level_height'][0]
    return level_height


def binned_var_comparison(var1, var2, pp_list1, pp_list2, twod_filter, n_bins=512, extremes=[None, None, None, None]):
#       Calculate bins
    if extremes[0] != None:
        var1_min = extremes[0]
        var1_max = extremes[1]
        var2_min = extremes[2]
        var2_max = extremes[3]
        var1_bin_edges = calc_hist_bins(var1_min, var1_max, n_bins, 'normal')
        var2_bin_edges = calc_hist_bins(var2_min, var2_max, n_bins, 'normal')
        bins = [var1_bin_edges, var2_bin_edges]
        extremes_str = ','.join(['%.5e' % i for i in extremes])
    else:
        extremes_str = str(extremes)

#   Attributes used to identify unique calculations
    attributes = {'filter_name' : twod_filter.attributes['filter_name'],
                  'wavenumber' : twod_filter.attributes['wavenumber'],
                  'delta_x' : twod_filter.attributes['delta_x'],
                  'width' : twod_filter.attributes['width'],
                  'cutoff' : twod_filter.attributes['cutoff'],
                  'high_pass' : twod_filter.attributes['high_pass'],
                  'sigma' : twod_filter.attributes['sigma'],
                  'var1' : var1,
                  'var2' : var2,
                  'pp_list1' : ','.join(pp_list1),
                  'pp_list2' : ','.join(pp_list2),
                  'extremes' : extremes_str,
                  'n_bins' : n_bins}

    new_calc, hist2d = read_cubes(hist2d_fileroot, attributes) # Read existing TKEs (if any exist)
    print 'pp_list1=', pp_list1
    print 'var1=', var1
    print 'twod_filter=', twod_filter
    if new_calc == 1: # If data does not already exists, calculate now.
        var2_var1_hist = np.zeros(n_bins*n_bins).reshape(n_bins, n_bins)
        if var2 == 'explicit TKE':
            for i in range(len(pp_list1)):
                data1 = filtered_variable(var1, pp_list1[i], twod_filter,
                  mean_profile = 0)
                level_heights = data1[0:,:,:].coord('level_height').points
                if (var1 == 'strain rate'):
                    data1 = ((data1[:-1,:,:])/np.sqrt(2.0))**2
                elif (var1 == 'buoyancy gradient'):
               	    data1 = data1[1:,:,:]
                data2 = explicit_TKE(pp_list2[i], twod_filter, mean_profile = 0)[0:-3,:,:]
                print 'data1.data.min()=', data1.data.min()
                print 'data1.data.max()=', data1.data.max()
                print 'data2.data.min()=', data2.data.min()
                print 'data2.data.max()=', data2.data.max()
                if extremes[0] != None:
                    check_bin_limits(data1, var1_min, var1_max)
                    check_bin_limits(data2, var2_min, var2_max)
                    heatmap, var1_edges, var2_edges = np.histogram2d(data1.data.flatten(), data2.data.flatten(), bins=bins)
                else:
                    heatmap, var1_edges, var2_edges = np.histogram2d(data1.data.flatten(), data2.data.flatten(), bins=n_bins)
                print 'heatmpa.shape=', heatmap.shape
                print 'var2_var1_hist.shape=', var2_var1_hist.shape
                var2_var1_hist += heatmap
                attributes['var1_bin_edges'] = var1_edges
                attributes['var2_bin_edges'] = var2_edges
            hist2d = write_cubes(hist2d_fileroot, var2_var1_hist, attributes)
            result = hist2d
        else:
            print "No code exists to do things other than explicit TKE"
    else:
        result = hist2d
    return result


def calc_hist_bins(data_min, data_max, n_bins, bin_type):
    if bin_type == 'normal':
        bin_edges = data_min + (data_max - data_min) * (np.arange(n_bins+1)/float(n_bins))
    elif bin_type == 'log':
        if data_min > 0:
            bin_edges = data_min + (data_max- data_min) * np.exp(np.arange(n_bins+1))/float(np.exp(n_bins))
        else:
            n_bins = n_bins/2
            temp = np.max([abs(data_min), data_max]) * np.exp(np.arange(n_bins+1))/float(np.exp(n_bins))
            bin_edges = np.concatenate((-temp[::-1], temp))
    elif bin_type == 'square':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**2.)/float(n_bins**2.)
    elif bin_type == 'cubic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**3.)/float(n_bins**3.)
    elif bin_type == 'quartic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**4.)/float(n_bins**4.)
    elif bin_type == 'quintic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**5.)/float(n_bins**5.)
    elif bin_type == 'sextic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**6.)/float(n_bins**6.)
    elif bin_type == 'septic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**7.)/float(n_bins**7.)
    elif bin_type == 'octic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**8.)/float(n_bins**8.)
    elif bin_type == 'nonic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**9.)/float(n_bins**9.)
    elif bin_type == 'decic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**10.)/float(n_bins**10.)
    elif bin_type == 'inverse_decic':
        bin_edges = data_min + (data_max- data_min) * (np.arange(n_bins+1)**0.1)/float(n_bins**0.1)
    else:
        print bin_type, ' is not an available bin type'
    return bin_edges


def check_bin_limits(data, min_val, max_val):
    if data.data.min() < min_val:
        print 'Min value in data ', data.data.min(), 'is less than minimum bin size', min_val
        stop
    if data.data.max() > max_val:
        print 'Max value in data ', data.data.max(), 'is greater than maximum bin size', max_val
        stop


def dev_stress_tensor(pp_file, stash1, stash2, twod_filter, mean_profile=1):
#   To prevent duplicate calculations of the same thing, force stash2 >= stash1
    if (int(stash1) > int(stash2)):
        stash1_old = stash1
        stash1 = stash2
        stash2 = stash1_old
    attributes = {'filter_name' : twod_filter.attributes['filter_name'],
                  'wavenumber' : twod_filter.attributes['wavenumber'],
                  'delta_x' : twod_filter.attributes['delta_x'],
                  'width' : twod_filter.attributes['width'],
                  'cutoff' : twod_filter.attributes['cutoff'],
                  'high_pass' : twod_filter.attributes['high_pass'],
                  'sigma' : twod_filter.attributes['sigma'],
                  'pp_file' : pp_file,
                  'mean_profile' : mean_profile,
                  'stash1' : stash1,
                  'stash2' : stash2}
    new_calc, result = read_cubes(dev_stress_fileroot, attributes)# Read existing stress_tensors (if any exist)
    if new_calc == 1:# If data does not already exists, calculate now.
        heights = level_heights(pp_file, '150')
        regrid_cube = read_variable(pp_file, '150')[0]
        stress_data = -quadratic_subfilter(pp_file, pp_file
                        , stash1, stash2, twod_filter.data
                        , mean_profile, levels=heights, regrid_cube=regrid_cube)
        if (stash1 == stash2):
            stress_data += (explicit_TKE(pp_file, twod_filter, mean_profile=mean_profile).data)/3.0
# Write this new filter to the file
        result = write_cubes(dev_stress_fileroot, stress_data, attributes, heights=heights)
    return result


def vertical_flux(var_name, pp_file, twod_filter, mean_profile=1):
    attributes = {'var_name' : var_name,
                  'filter_name' : twod_filter.attributes['filter_name'],
                  'wavenumber' : twod_filter.attributes['wavenumber'],
                  'delta_x' : twod_filter.attributes['delta_x'],
                  'width' : twod_filter.attributes['width'],
                  'cutoff' : twod_filter.attributes['cutoff'],
                  'high_pass' : twod_filter.attributes['high_pass'],
                  'sigma' : twod_filter.attributes['sigma'],
                  'pp_file' : pp_file,
                  'mean_profile' : mean_profile}
    new_calc, result = read_cubes(vertical_flux_fileroot, attributes)# Read existing stress_tensors (if any exist)
    if new_calc == 1:# If data does not already exists, calculate now.
        heights = level_heights(pp_file, '150')
        regrid_cube = read_variable(pp_file, '150')[0]
        data = quadratic_subfilter(pp_file, pp_file, '150', 
                                   stash_dict[var_name], twod_filter.data,
                                   mean_profile, levels=heights, regrid_cube=regrid_cube)
# Write this new filter to the file
        result = write_cubes(vertical_flux_fileroot, data, attributes, heights=heights)
    return result


def quadratic_subfilter(pp_file1, pp_file2, stash_code1, stash_code2, \
                        twod_filter, mean_profile, levels=None, \
                        regrid_cube=None) :
    '''
    Calculates s(a,b)=(a*b)^r-a^r*b^r
    a = stash_code1 in pp_file1
    b = stash_code2 in pp_file2

    ()^r is twod_filter

    If needed, will interpolate to levels or regrid to regrid_cube

    Modified to return [s(a,b), a^r, b^r] as list of arrays
    '''
    if levels == None:
        field1 = read_variable(pp_file1, stash_code1)[0]
    else :
        field1 = read_variable(pp_file1, stash_code1)[0]
        field1 = field1.interpolate([('level_height', levels)] \
                                      , iris.analysis.Linear(extrapolation_mode='nan'))
    if regrid_cube != None:
        field1 = periodic_regrid(field1, regrid_cube)
    field1 = field1.data
    if stash_code2 == stash_code1:
        field2 = field1
    else:
        if levels == None:
            field2 = read_variable(pp_file2, stash_code2)[0]
        else:
            field2 = read_variable(pp_file2, stash_code2)[0]
            field2 = field2.interpolate([('level_height', levels)] \
                                          , iris.analysis.Linear(extrapolation_mode='nan'))
        if regrid_cube != None:
            field2 = periodic_regrid(field2, regrid_cube)
        field2 = field2.data
    field1_r = np.zeros(field1.size).reshape(field1.shape)
    field2_r = np.zeros(field2.size).reshape(field2.shape)
    field1_field2_r = np.zeros(field1.size).reshape(field1.shape)
    for i in range(field1.shape[0]):
        field1_r[i,:,:] = convolve(field1[i,:,:], twod_filter)
        field1_field2_r[i,:,:] = convolve( (field1[i,:,:] * field2[i,:,:]), 
                                           twod_filter)
        if stash_code2 != stash_code1:
            field2_r[i,:,:] = convolve(field2[i,:,:], twod_filter)
        else:
            field2_r[i,:,:] = field1_r[i,:,:]
    result = field1_field2_r -(field1_r * field2_r)
    if mean_profile == 1:
        result = [ calc_mean_profile(result), \
                   calc_mean_profile(field1_r), \
                   calc_mean_profile(field1_r) ]
    else :
        result = [result, field1_r, field2_r]
    return result


def convolve(field, twod_filter):
    pad_len = int(np.ceil(len(twod_filter)))
    field = np.pad(field, pad_len, mode='wrap')
    result = fftconvolve(field, twod_filter, mode='same')
    return result[pad_len:-pad_len, pad_len:-pad_len]


def delta_horizontal(cube, dimension='latitude'):
    dim_vals = cube.coords('grid_'+dimension)[0].points
    return ((dim_vals[1]-dim_vals[0])*pi)/180.


def periodic_regrid(cube1, cube2):
    cube1_lat = cube1.coords('grid_latitude')[0].points
    cube1_lon = cube1.coords('grid_longitude')[0].points
    cube2_lat = cube2.coords('grid_latitude')[0].points
    cube2_lon = cube2.coords('grid_longitude')[0].points
    cube1_new = cube1.regrid(cube2, iris.analysis.Linear(extrapolation_mode='nan'))
    if cube1_lat[0] > cube2_lat[0]:
        wt_mn = (cube1_lat[0] - cube2_lat[0]) / (cube2_lat[1] - cube2_lat[0])
        wt_pl = 1.0 - wt_mn
        cube1_new.data[:,0,:] = (wt_pl * cube1.data[:,0,:]) + (wt_mn * cube1.data[:,-1,:])
    if cube1_lat[0] < cube2_lat[0]:
        wt_mn = (cube1_lat[0] - cube2_lat[0]) / (cube2_lat[1] - cube2_lat[0])
        wt_pl = 1.0 - wt_mn
        cube1_new.data[:,-1,:] = (wt_pl * cube1.data[:,0,:]) + (wt_mn * cube1.data[:,-1,:])
    if cube1_lon[0] > cube2_lon[0]:
        wt_mn = (cube1_lon[0] - cube2_lon[0]) / (cube2_lon[1] - cube2_lon[0])
        wt_pl = 1.0 - wt_mn
        cube1_new.data[:,:,0] = (wt_pl * cube1.data[:,:,0]) + (wt_mn * cube1.data[:,:,-1])
    if cube1_lon[0] < cube2_lon[0]:
        wt_mn = (cube1_lon[0] - cube2_lon[0]) / (cube2_lon[1] - cube2_lon[0])
        wt_pl = 1.0 - wt_mn
        cube1_new.data[:,:,-1] = (wt_pl * cube1.data[:,:,0]) + (wt_mn * cube1.data[:,:,-1])
    return cube1_new


def Smag_buoyancy_flux(pp_file1, pp_file2, mean_profile=1):
    attributes = {'pp_file1':pp_file1,
                  'pp_file2':pp_file2, 
                  'mean_profile':mean_profile}
    new_calc, result = read_cubes(Smag_buoyancy_flux_fileroot, attributes) # Read existing buoyancy fluxess (if any exist)
    if new_calc == 1: # If data does not already exists, calculate now.
        heights = level_heights(pp_file1, '150') # defined on theta levels and points
        theta_v = read_variable(pp_file1, '-8888')
        k_h = -1.0 * read_variable(pp_file2, stash_dict['visc_h'])[0].data[:-1,:,:]
        d_theta_v_dz = velocity_gradients.dw_dz(heights, heights, theta_v[0].data, surface_value='same')[0:-3,:,:]
        heights = heights[0:-3]
        buoyancy_flux = k_h * d_theta_v_dz
        if mean_profile == 1:
            buoynacy_flux = calc_mean_profile(buoyancy_flux)
        result = write_cubes(Smag_buoyancy_flux_fileroot, buoyancy_flux, attributes, heights=heights)
    return result


def Smag_humidity_flux(pp_file1, pp_file2, mean_profile=1):
    attributes = {'pp_file1':pp_file1, 
                  'pp_file2':pp_file2, 
                  'mean_profile':mean_profile}
    new_calc, result = read_cubes(Smag_humidity_flux_fileroot, attributes) # Read existing humidity fluxess (if any exist)
    if new_calc == 1: # If data does not already exists, calculate now.
        heights = level_heights(pp_file1, '150') # defined on theta levels and points
        q = read_variable(pp_file1, '10')
        k_h = -1.0 * read_variable(pp_file2, stash_dict['visc_h'])[0].data[:-1,:,:]
        dq_dz = velocity_gradients.dw_dz(heights, heights, q[0].data, surface_value='same')[0:-3,:,:]
        heights = heights[0:-3]
        humidity_flux = k_h * dq_dz
        if mean_profile == 1:
            humidity_flux = calc_mean_profile(humidity_flux)
        result = write_cubes(Smag_humidity_flux_fileroot, humidity_flux, attributes, heights=heights)
    return result


def Smag_dev_stress(pp_file, stash1, stash2, mean_profile=1):
#   To prevent duplicate calculations of the same thing, force stash2 >= stash1
    if (int(stash1) > int(stash2)):
        stash1_old = stash1
        stash1 = stash2
        stash2 = stash1_old
    attributes = {'pp_file':pp_file,
                  'stash1':stash1, 
                  'stash2':stash2, 
                  'mean_profile':mean_profile}
    wind = {'2':'u', '3':'v', '150':'w'}
    dimension = {'2':'x', '3':'y', '150':'z'}
    new_calc, result = read_cubes(Smag_dev_stress_fileroot, attributes) # Read existing humidity fluxess (if any exist)
    if new_calc == 1: # If data does not already exists, calculate now.
        theta_levels = level_heights(pp_file, '150') # defined on theta levels and points
        rho_levels = level_heights(pp_file, '2')
        velocity1 = read_variable(pp_file, stash1)[0]
        velocity2 = read_variable(pp_file, stash2)[0]
        delta_lambda = delta_horizontal(velocity1)
        s_ij = eval('0.5*(velocity_gradients.d'+wind[stash1]+'_d'+dimension[stash2]+'(theta_levels, rho_levels, velocity1.data, delta_lambda=delta_lambda, surface_value="zero")+velocity_gradients.d'+wind[stash2]+'_d'+dimension[stash1]+'(theta_levels, rho_levels, velocity2.data, delta_lambda=delta_lambda, surface_value="zero"))')
        k_m = 1.0
        dev_stress = k_m * s_ij
        if mean_profile == 1:
            dev_stress = calc_mean_profile(dev_stress)
        result = write_cubes(Smag_dev_stress_fileroot, dev_stress, attributes, heights=theta_levels)
    return result


def calc_theta_v(pp_file):
    theta = iris.load(data_dir + pp_file,
                      iris.AttributeConstraint(STASH = 'm01s00i004'))[0]
    q_v = iris.load(data_dir + pp_file,
                    iris.AttributeConstraint(STASH = 'm01s00i010'))[0]
    q_l = iris.load(data_dir + pp_file,
                    iris.AttributeConstraint(STASH = 'm01s00i254'))[0]
    cubes = theta*(1.0+0.61*q_v-q_l)
    return [cubes]
      

def calc_theta_l(pp_file):
    theta = iris.load(data_dir + pp_file,
                      iris.AttributeConstraint(STASH = 'm01s00i004'))[0]
    q_l = iris.load(data_dir + pp_file,
                    iris.AttributeConstraint(STASH = 'm01s00i254'))[0]
    lc_over_cp = theta.copy()
    lc_over_cp.data[:,:,:] = lc/cp
    cubes = theta - (lc_over_cp*q_l)
    return [cubes]
      

def calc_mean_profile(variable):
    variable = variable.mean(axis=2,dtype=np.float64).mean(axis=1,dtype=np.float64)
    variable = variable.reshape(len(variable), 1, 1)
    return variable

#def dev_stress_comparison():
    
    


def main():
    print 'Have moved all code to separate files to make it easier to maintain'

if __name__ == "__main__":
    main()


def diagnostic_buoyancy_flux(pp_file1, pp_file2, mean_profile=1):
    print 'not ready yet'
    attributes = {'pp_file1':pp_file1,
                  'pp_file2':pp_file2, 
                  'mean_profile':mean_profile}
    new_calc, result = read_cubes(diagnostic_buoyancy_flux_fileroot, attributes) # Read existing buoyancy fluxess (if any exist)
    if new_calc == 1: # If data does not already exists, calculate now.
        heights = level_heights(pp_file1, '150') # defined on theta levels and points
        theta_v = read_variable(pp_file1, '-8888')
        rho = read_variable(pp_file1, '253')[0].interpolate([('level_height', heights[0:-3])], iris.analysis.Linear(extrapolation_mode='nan')).data/(earth_radius**2)
        heat_flux_diag = read_variable(pp_file2, '3216')[0].interpolate([('level_height', heights[0:-3])], iris.analysis.Linear(extrapolation_mode='nan')).data
        heat_flux_diag = heat_flux_diag/(rho*cp)
        theta_l = read_variable(pp_file1, '-7777')
        d_theta_v_dz = velocity_gradients.dw_dz(heights, heights, theta_v[0].data, surface_value='same')[0:-3,:,:]
        d_theta_l_dz = velocity_gradients.dw_dz(heights, heights, theta_l[0].data, surface_value='same')[0:-3,:,:]
        print 'heat_flux_diag=', heat_flux_diag.mean(axis=2,dtype=np.float64).mean(axis=1,dtype=np.float64)
        print 'd_theta_l_dz =', d_theta_l_dz.mean(axis=2,dtype=np.float64).mean(axis=1,dtype=np.float64)
        k_h = heat_flux_diag/(d_theta_l_dz+eps) # This assumes that the heat flux is still using the k approach. Need to check this is indeed the case.
        ind = np.where(np.isinf(k_h.flatten()))[0]
        if len(ind) > 1:
            ind = ind[0]
            print 'k_h.flatten()[ind]=', k_h.flatten()[ind]
            print 'heat_flux_diag.flatten()[ind]=',heat_flux_diag.flatten()[ind]
            print 'd_theta_l_dz.flatten()[ind]=',d_theta_l_dz.flatten()[ind]
        heights = heights[0:-3]
        print 'k_h=', k_h.mean(axis=2,dtype=np.float64).mean(axis=1,dtype=np.float64)
        print 'd_theta_v_dz=', d_theta_v_dz.mean(axis=2,dtype=np.float64).mean(axis=1,dtype=np.float64)
        buoyancy_flux = k_h * d_theta_v_dz
        if mean_profile == 1:
            buoyancy_flux = calc_mean_profile(buoyancy_flux)
        result = write_cubes(Smag_buoyancy_flux_fileroot, buoyancy_flux, attributes, heights=heights)
    return result


def filtered_strain_rate(pp_file, twod_filter, mean_profile=1):
    attributes = {'filter_name' : twod_filter.attributes['filter_name'],
                  'wavenumber' : twod_filter.attributes['wavenumber'],
                  'delta_x' : twod_filter.attributes['delta_x'],
                  'width' : twod_filter.attributes['width'],
                  'cutoff' : twod_filter.attributes['cutoff'],
                  'high_pass' : twod_filter.attributes['high_pass'],
                  'sigma' : twod_filter.attributes['sigma'],
                  'pp_file' : pp_file,
                  'mean_profile' : mean_profile}
    new_calc, result = read_cubes(filtered_strain_rate_fileroot, attributes) # Read existing data (if exists)
    if new_calc == 1: # If data does not already exists, calculate now.
        s_11_r = filtered_field_calc(Smag_dev_stress(pp_file, '2', '2', mean_profile).data, twod_filter.data)[0]
        s_12_r = filtered_field_calc(Smag_dev_stress(pp_file, '2', '3', mean_profile).data, twod_filter.data)[0]
        s_13_r = filtered_field_calc(Smag_dev_stress(pp_file, '2', '150', mean_profile).data, twod_filter.data)[0]
        s_22_r = filtered_field_calc(Smag_dev_stress(pp_file, '3', '3', mean_profile).data, twod_filter.data)[0]
        s_23_r = filtered_field_calc(Smag_dev_stress(pp_file, '3', '150', mean_profile).data, twod_filter.data)[0]
        s_33_r = filtered_field_calc(Smag_dev_stress(pp_file, '150', '150', mean_profile).data, twod_filter.data)[0]
        strain_rate = np.sqrt(2.0 * ( (s_11_r * s_11_r) +
                                    (2.0 * (s_12_r * s_12_r)) +
                                    (2.0 * (s_13_r * s_13_r)) +
                                    (2.0 * (s_23_r * s_23_r)) +
                                    (s_22_r * s_22_r) + 
                                    (s_33_r * s_33_r) ))
        heights = level_heights(pp_file, '150') # defined on theta levels and points
        result = write_cubes(filtered_strain_rate_fileroot, strain_rate, attributes, heights=heights)
    return result
    
    
