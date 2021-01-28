import os
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import subfilter as sf
import filters as filt
import difference_ops as do

options = {
#        'FFT_type': 'FFTconvolve',
#        'FFT_type': 'FFT',
        'FFT_type': 'RFFT',
        'save_all': 'Yes',
        'th_ref': 300.0,
          }


dir = '/storage/silver/scenario/si818415/phd/first_data/'
#dir = '/storage/silver/wxproc/xm904103/'
#odir = '/storage/silver/scenario/si818415/phd/first_data/'
odir = '/storage/silver/wxproc/xm904103/'
odir = odir + 'test_op_' + options['FFT_type']+'/'

os.makedirs(odir, exist_ok = True)

file = 'diagnostics_3d_ts_13200.nc'
#ref_file = 'diagnostics_3d_ts_13200.nc'
ref_file = None
#w = dataset.variables['w']

plot_dir = odir + 'plots/'
os.makedirs(plot_dir, exist_ok = True)

plot_type = '.png'
data_dir = '' # Directory containing data
figshow = True

def main():
    '''
    Peter's top level code, altered slightly.
    '''
#   Non-global variables that are set once
    width_list = [3, 5, 20, 40, 80, 200]
#    width_list = [3]
    filter_name = 'running_mean'
#    filter_name = 'wave_cutoff'

    dataset = Dataset(dir+file, 'r') # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
                                 
    od = sf.options_database(dataset)
    if od is None:
        dx = 5.0
        dy = 5.0
    else :
        dx = float(od['dxx'])
        dy = float(od['dyy'])
    
    if ref_file is not None:
        ref_path = dir+ref_file
        ref_dataset = Dataset(dir+ref_file, 'r')
    else:
        ref_path = ref_file
        ref_dataset = None
    ilev = 15
    iy = 40

    opgrid = 'w'
    fname = 'test_plot'

    derived_dataset_name, derived_data, exists = \
        sf.setup_derived_data_file( dir+file, odir, ref_path, fname,
                                   options, override=True)
        
    print(f"Derived dataset name:{derived_dataset_name:s}")
    print("Variables in derived dataset.")
    print(derived_data.variables)

    filter_list = list([])

    for i,width in enumerate(width_list):
        if filter_name == 'gaussian':
            filter_id = 'filter_ga{:02d}'.format(i)
            twod_filter = filt.filter_2d(filter_id,
                                       filter_name,
                                       sigma=sigma, width=width,
                                       delta_x=dx)
        elif filter_name == 'wave_cutoff':
            filter_id = 'filter_wc{:02d}'.format(i)
            twod_filter = filt.filter_2d(filter_id,
                                       filter_name, wavenumber=np.pi/(2*sigma),
                                       width=width,
                                       delta_x=dx)
        elif filter_name == 'running_mean':
            filter_id = 'filter_rm{:02d}'.format(i)
            twod_filter = filt.filter_2d(filter_id,
                                       filter_name,
                                       width=width,
                                       delta_x=dx)

        print(twod_filter)
        filter_list.append(twod_filter)

# Add whole domain filter
    filter_name = 'domain'
    # filter_id = 'filter_do{:02d}'.format(len(filter_list))
    twod_filter = filt.filter_2d(filter_id, filter_name, delta_x=dx)
    filter_list.append(twod_filter)

    print(filter_list)

    z = do.last_dim(dataset["z"])
    zn = do.last_dim(dataset["zn"])

    for twod_filter in filter_list:

        print(twod_filter)

        filtered_dataset_name, filtered_data, exists = \
            sf.setup_filtered_data_file( dir+file, odir, ref_path, fname,
                                       options, twod_filter, override=True)
        print("Variables in filtered dataset.")
        print(filtered_data.variables)
#        exists = False
        if exists :
            print('Derived data file exists' )
        else :

            field_list = sf.filter_variable_list(dataset, ref_dataset,
                                                 derived_data, filtered_data,
                                                 options,
                                                 twod_filter, var_list=None,
                                                 grid = opgrid)

            quad_field_list = sf.filter_variable_pair_list(dataset, 
                                                  ref_dataset,
                                                  derived_data, filtered_data,
                                                  options,
                                                  twod_filter, var_list=None,
                                                  grid = opgrid)


        z = derived_data["z"]
        zn = derived_data["zn"]

        if twod_filter.attributes['filter_type']!='domain' :
            fig1 = plt.figure(1)
            plt.contourf(twod_filter.data,20)
            plt.savefig(plot_dir+'Filter_'+\
                        twod_filter.id+plot_type)
            plt.close()

            fig2 = plt.figure(2)
            plt.plot(twod_filter.data[np.shape(twod_filter.data)[0]//2,:])
            plt.savefig(plot_dir+'Filter_y_xsect_'+\
                        twod_filter.id+plot_type)
            plt.close()
            
        filtered_data.close()
    derived_data.close()
    dataset.close()

if __name__ == "__main__":
    main()