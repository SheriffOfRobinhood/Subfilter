"""
Created on Tue Feb  1 10:40:07 2022

@author: paclk
"""
#%%
import os

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import subfilter
import subfilter.subfilter as sf
import subfilter.utils.cloud_monc as cldm
from subfilter.io.datain import (configure_model_resolution,
                                 get_data_on_grid,
                                 )
from subfilter.utils.string_utils import get_string_index

test_case = 0
plot_type = '.png'
figshow = True

if test_case == 0:
    config_file = 'config_test_case_0.yaml'
    indir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
    odir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
    file = 'diagnostics_3d_ts_21600.nc'
    ref_file = 'diagnostics_ts_21600.nc'
elif test_case == 1:
    config_file = 'config_test_case_1.yaml'
    indir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
    odir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
    file = 'diagnostics_3d_ts_13200.nc'
    ref_file = None
options, update_config = sf.subfilter_options(config_file)


fname = 'test_rewrite'

odir = odir + fname +'_'+ options['FFT_type']+'/'

plot_dir = odir + 'plots/'

subfilter.global_config['test_level'] = 2
# Read data
dataset = xr.open_dataset(indir+file)

print(dataset)

if ref_file is not None:
    ref_dataset = xr.open_dataset(indir+ref_file)
else:
    ref_dataset = None

# Get model resolution values
dx, dy, options = configure_model_resolution(dataset, options)


# For plotting
#    ilev = 15
ilev = 40
#    iy = 40
iy = 95
nlevels = 40
opgrid = 'w'

filter_id = 'filter_ga03'

basename = ('.').join(file.split('.')[:-1])

derived_dataset_name = basename + '_' + fname + '.nc'
derived_dataset = xr.open_dataset(odir+derived_dataset_name)
derived_data = {'file':odir+derived_dataset_name, 'ds':derived_dataset}

print(derived_data)

filtered_dataset_name = basename + '_' + fname + '_' + filter_id + '.nc'

filtered_dataset = xr.open_dataset(odir+filtered_dataset_name)
print(filtered_dataset)

#%%
[itime, iix, iiy, iiz] = get_string_index(filtered_dataset.dims, ['time', 'x', 'y', 'z'])
[timevar, xvar, yvar, zvar] = [list(filtered_dataset.dims)[i] for i in [itime, iix, iiy, iiz]]

dbdz_r   = filtered_dataset["f(moist_dbdz_on_w)_r"]
modSsq_r = filtered_dataset["f(mod_S_sqn)_r"]

ri1 = cldm.richardson(modSsq_r, dbdz_r)

rilev = np.linspace(-1,1,41)

for it, t in enumerate(filtered_dataset.coords[timevar][:1]):

    fig, axa = plt.subplots(3,2,figsize=(10,10))

    dbdz_r.isel({timevar:it,zvar:ilev}).plot.contourf(x=xvar, y=yvar, ax=axa[0,0], levels=nlevels)
    dbdz_r.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(x=xvar, y=zvar, ax=axa[0,1], levels=nlevels)

    modSsq_r.isel({timevar:it,zvar:ilev}).plot.contourf(x=xvar, y=yvar, ax=axa[1,0], levels=nlevels)
    modSsq_r.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(x=xvar, y=zvar, ax=axa[1,1], levels=nlevels)

    #ri1 = dbdz_r / modSsq_r
    #ri1.name = 'moist_Ri'

    ri1.isel({timevar:it,zvar:ilev}).plot.contourf(x=xvar, y=yvar, ax=axa[2,0], levels=rilev)
    ri1.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(x=xvar, y=zvar, ax=axa[2,1], levels=rilev)

    fig.tight_layout()
    plt.savefig(f'{plot_dir}Richardson_lev_{ilev:03d}_x_z{iy:03d}_{filter_id}_{it:02d}{plot_type}')
    plt.close()

#%%
th_ref = get_data_on_grid(dataset, ref_dataset, derived_data, 'thref', options, grid = opgrid)
p_ref  = get_data_on_grid(dataset, ref_dataset, derived_data, 'pref',  options, grid = opgrid)
th_L = filtered_dataset["f(th_L_on_w)_r"]
qt   = filtered_dataset["f(q_total_on_w)_r"]
s_qt_qt =filtered_dataset["s(q_total,q_total)_on_w"]
s_thL_qt =filtered_dataset["s(th_L,q_total)_on_w"]
s_thL_thL =filtered_dataset["s(th_L,th_L)_on_w"]

(delta_q, qc, sig_s, cloud_fraction, qcl) = cldm.gaussian_cloud(th_L, qt, th_ref, p_ref, s_qt_qt, s_thL_qt, s_thL_thL)
#%%
for it, t in enumerate(filtered_dataset.coords[timevar]):

    fig, axa = plt.subplots(3,2,figsize=(10,10))

    sig_s.isel({timevar:it,zvar:ilev}).plot.contourf(x=xvar, y=yvar, ax=axa[0,0], levels=nlevels)
    sig_s.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(x=xvar, y=zvar, ax=axa[0,1], levels=nlevels)

    cloud_fraction.isel({timevar:it,zvar:ilev}).plot.contourf(x=xvar, y=yvar, ax=axa[1,0], levels=nlevels)
    cloud_fraction.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(x=xvar, y=zvar, ax=axa[1,1], levels=nlevels)


    qcl.isel({timevar:it,zvar:ilev}).plot.contourf(x=xvar, y=yvar, ax=axa[2,0], levels=nlevels)
    qcl.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(x=xvar, y=zvar, ax=axa[2,1], levels=nlevels)

    fig.tight_layout()
    plt.savefig(f'{plot_dir}Cloud_lev_{ilev:03d}_x_z{iy:03d}_{filter_id}_{it:02d}{plot_type}')
    plt.close()
