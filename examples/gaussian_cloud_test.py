"""
Created on Tue Feb  1 10:40:07 2022

@author: paclk
"""
#%%
import os

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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


fname = 'gaussian_cloud'

odir = odir + fname +'_'+ options['FFT_type']+'/'

plot_dir = odir + 'plots/'

# Read data
dataset = xr.open_dataset(indir+file)

print(dataset)

if ref_file is not None:
    ref_dataset = xr.open_dataset(indir+ref_file)
else:
    ref_dataset = None

# Get model resolution values
dx, dy, options = configure_model_resolution(dataset, options)
options['save_all'] = 'no'

# For plotting
#    ilev = 15
ilev = 40
#    iy = 40
iy = 95
nlevels = 40
opgrid = 'p'

filter_no = -1
if filter_no == -1:
    filter_id = 'filter_do'
else:
    filter_id = f'filter_ga{filter_no:02d}'

basename = ('.').join(file.split('.')[:-1])

derived_dataset_name = basename + '_' + fname + '.nc'
derived_dataset = xr.open_dataset(odir+derived_dataset_name)
derived_data = {'file':odir+derived_dataset_name, 'ds':derived_dataset}

print(derived_data)

filtered_dataset_name = basename + '_' + fname + '_' + filter_id + '.nc'

filtered_dataset = xr.open_dataset(odir+filtered_dataset_name)
print(filtered_dataset)

#%%

#%%
th_ref = get_data_on_grid(dataset, ref_dataset, derived_data, 'thref', options, grid = opgrid)
p_ref  = get_data_on_grid(dataset, ref_dataset, derived_data, 'pref',  options, grid = opgrid)
th_L = filtered_dataset[f"f(th_L_on_{opgrid})_r"]
qt   = filtered_dataset[f"f(q_total_on_{opgrid})_r"]
s_qt_qt =filtered_dataset[f"s(q_total,q_total)_on_{opgrid}"]
s_thL_qt =filtered_dataset[f"s(th_L,q_total)_on_{opgrid}"]
s_thL_thL =filtered_dataset[f"s(th_L,th_L)_on_{opgrid}"]

# cloud_fraction_true = np.clip(filtered_dataset[f"f(cloud_fraction_on_{opgrid})_r"],0,None)
# qcl_true = np.clip(filtered_dataset[f"f(q_cloud_liquid_mass_on_{opgrid})_r"],0,None)
cloud_fraction_true = filtered_dataset[f"f(cloud_fraction_on_{opgrid})_r"]
qcl_true = filtered_dataset[f"f(q_cloud_liquid_mass_on_{opgrid})_r"]


(delta_q, qc, sig_s, cloud_fraction, qcl) = cldm.gaussian_cloud(th_L, qt, th_ref, p_ref, s_qt_qt, s_thL_qt, s_thL_thL)

#%%

Cmax = [1, 1, 0.4, 0.2]
qmax = [1E-3, 5E-4, 1E-4, 5E-5]
Clev = np.linspace(0,Cmax[filter_no],21)
qlev = np.linspace(0,qmax[filter_no],21)

if filter_no == -1:
    [itime, iiz] = get_string_index(th_L.dims, ['time', 'z'])
    [timevar, zvar] = [list(th_L.dims)[i] for i in [itime, iiz]]

    for it, t in enumerate(filtered_dataset.coords[timevar]):
        print(f"it: {it}")
        fig1, axa = plt.subplots(3,1,figsize=(10,20))

        sig_s.isel({timevar:it,zvar:slice(1,None)}).plot(
            y=zvar, ax=axa[0], label=True)
        x = cloud_fraction.isel({timevar:it, zvar:slice(1,None)})
        y = x.coords[zvar].values
        l1 = axa[1].plot(x, y, label=x.name)

        x = cloud_fraction_true.isel({timevar:it, zvar:slice(1,None)})
        y = x.coords[zvar].values
        l2 = axa[1].plot(x, y, label=x.name)

        axa[1].set_xlabel(r"Cloud Fraction")
        axa[1].set_ylabel(zvar)
        axa[1].legend()

        x = qcl.isel({timevar:it, zvar:slice(1,None)})
        y = x.coords[zvar].values
        l1 = axa[2].plot(x, y, label=x.name)

        x = qcl_true.isel({timevar:it, zvar:slice(1,None)})
        y = x.coords[zvar].values
        l2 = axa[2].plot(x, y, label=x.name)

        axa[2].set_xlabel(r"$q_{cl}$")
        axa[2].set_ylabel(zvar)
        axa[2].legend()

        # cloud_fraction.isel({timevar:it,zvar:slice(1,None)}).plot(
        #     y=zvar, ax=axa[1])

        # cloud_fraction_true.isel({timevar:it,zvar:slice(1,None)}).plot(
        #     y=zvar, ax=axa[1])

        # qcl.isel({timevar:it,zvar:slice(1,None)}).plot(
        #     y=zvar, ax=axa[2])

        # qcl_true.isel({timevar:it,zvar:slice(1,None)}).plot(
        #     y=zvar, ax=axa[2])

        fig1.tight_layout()
        fig1.savefig(f'{plot_dir}Cloud_{filter_id}_{it:02d}{plot_type}')
        plt.close()

        fig3, axc = plt.subplots(1,2,figsize=(5,2.5))
        heatmap, xedges, yedges = np.histogram2d(cloud_fraction_true.isel({timevar:it}).values.flatten(),
                                                 cloud_fraction.isel({timevar:it}).values.flatten(),
                                                 bins=[Clev, Clev])
        heatmap[0:1,0:1] =0
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        im1 = axc[0].imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm())
        axc[0].set_xlabel('True Cloud Fraction')
        axc[0].set_ylabel('Gaussian Cloud Fraction')
        plt.colorbar(im1, ax=axc[0])

        heatmap, xedges, yedges = np.histogram2d(qcl_true.isel({timevar:it}).values.flatten(),
                                                 qcl.isel({timevar:it}).values.flatten(),
                                                 bins=[qlev, qlev])
        heatmap[0:1,0:1] =0
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        im2 = axc[1].imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm())
        axc[1].set_xlabel(r'True $q_{cl}$')
        axc[1].set_ylabel(r'Gaussian $q_{cl}$')
        plt.colorbar(im2, ax=axc[1])
        fig3.tight_layout()
        fig3.savefig(f'{plot_dir}Cloud_heatmap_{filter_id}_{it:02d}{plot_type}')
        plt.close()

    fig1, axa = plt.subplots(3,1,figsize=(6,10))

    sig_s.isel({zvar:slice(1,None)}).mean(dim=timevar).plot.line(
        y=zvar, ax=axa[0])


    x = cloud_fraction.isel({zvar:slice(1,None)}).mean(dim=timevar)
    y = x.coords[zvar].values
    l1 = axa[1].plot(x, y, label=x.name)

    x = cloud_fraction_true.isel({zvar:slice(1,None)}).mean(dim=timevar)
    y = x.coords[zvar].values
    l2 = axa[1].plot(x, y, label=x.name)

    axa[1].set_ylabel(zvar)

    axa[1].legend()

    x = qcl.isel({zvar:slice(1,None)}).mean(dim=timevar)
    y = x.coords[zvar].values
    l1 = axa[2].plot(x, y, label=x.name)

    x = qcl_true.isel({zvar:slice(1,None)}).mean(dim=timevar)
    y = x.coords[zvar].values
    l2 = axa[2].plot(x, y, label=x.name)

    axa[2].set_xlabel(r"$q_{cl}$")
    axa[2].set_ylabel(zvar)
    axa[2].legend()
    fig1.tight_layout()
    fig1.savefig(f'{plot_dir}Cloud_{filter_id}{plot_type}')
    plt.close()


else:
    [itime, iix, iiy, iiz] = get_string_index(th_L.dims, ['time', 'x', 'y', 'z'])
    [timevar, xvar, yvar, zvar] = [list(th_L.dims)[i] for i in [itime, iix, iiy, iiz]]
    # cmap = plt.get_cmap('GnBu')
    cmap = plt.get_cmap('YlGnBu')
    # cmap = plt.get_cmap('plasma')
    for it, t in enumerate(filtered_dataset.coords[timevar]):
    # for it in [5]:
        print(f"it: {it}")

        fig1, axa = plt.subplots(5,2,figsize=(10,20))

        sig_s.isel({timevar:it,zvar:ilev}).plot.contourf(
            x=xvar, y=yvar, ax=axa[0,0], levels=nlevels, cmap=cmap)
        sig_s.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(
            x=xvar, y=zvar, ax=axa[0,1], levels=nlevels, cmap=cmap)

        cloud_fraction.isel({timevar:it,zvar:ilev}).plot.contourf(
            x=xvar, y=yvar, ax=axa[1,0], levels=Clev, cmap=cmap)
        cloud_fraction.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(
            x=xvar, y=zvar, ax=axa[1,1], levels=Clev, cmap=cmap)

        cloud_fraction_true.isel({timevar:it,zvar:ilev}).plot.contourf(
            x=xvar, y=yvar, ax=axa[2,0], levels=Clev, cmap=cmap)
        cloud_fraction_true.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(
            x=xvar, y=zvar, ax=axa[2,1], levels=Clev, cmap=cmap)

        qcl.isel({timevar:it,zvar:ilev}).plot.contourf(
            x=xvar, y=yvar, ax=axa[3,0], levels=qlev, cmap=cmap)
        qcl.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(
            x=xvar, y=zvar, ax=axa[3,1], levels=qlev, cmap=cmap)

        qcl_true.isel({timevar:it,zvar:ilev}).plot.contourf(
            x=xvar, y=yvar, ax=axa[4,0], levels=qlev, cmap=cmap)
        qcl_true.isel({timevar:it,yvar:iy,zvar:slice(1,None)}).plot.contourf(
            x=xvar, y=zvar, ax=axa[4,1], levels=qlev, cmap=cmap)

        fig1.tight_layout()
        fig1.savefig(f'{plot_dir}Cloud_lev_{ilev:03d}_x_z{iy:03d}_{filter_id}_{it:02d}{plot_type}')
        plt.close()

        fig3, axc = plt.subplots(1,2,figsize=(7,3))
        heatmap, xedges, yedges = np.histogram2d(cloud_fraction_true.isel({timevar:it}).values.flatten(),
                                                 cloud_fraction.isel({timevar:it}).values.flatten(),
                                                 bins=[Clev, Clev])
        heatmap[0:1,0:1] =0
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        im1 = axc[0].imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm())
        axc[0].set_xlabel('True Cloud Fraction')
        axc[0].set_ylabel('Gaussian Cloud Fraction')
        plt.colorbar(im1, ax=axc[0])

        heatmap, xedges, yedges = np.histogram2d(qcl_true.isel({timevar:it}).values.flatten(),
                                                 qcl.isel({timevar:it}).values.flatten(),
                                                 bins=[qlev, qlev])
        heatmap[0:1,0:1] =0
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        im2 = axc[1].imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm())
        axc[1].set_xlabel(r'True $q_{cl}$')
        axc[1].set_ylabel(r'Gaussian $q_{cl}$')
        plt.colorbar(im2, ax=axc[1])
        fig3.tight_layout()
        fig3.savefig(f'{plot_dir}Cloud_heatmap_{filter_id}_{it:02d}{plot_type}')
        plt.close()

    fig1, axa = plt.subplots(3,1,figsize=(6,10))

    sig_s_mean = np.sqrt((sig_s * sig_s).mean(dim=[timevar,xvar,yvar]))
    sig_s_mean.isel({zvar:slice(1,None)}).plot.line(
        y=zvar, ax=axa[0])

    x = cloud_fraction.isel({zvar:slice(1,None)}).mean(dim=[timevar,xvar,yvar])
    y = x.coords[zvar].values
    l1 = axa[1].plot(x, y, label=x.name)

    x = cloud_fraction_true.isel({zvar:slice(1,None)}).mean(dim=[timevar,xvar,yvar])
    y = x.coords[zvar].values
    l2 = axa[1].plot(x, y, label=x.name)

    axa[1].set_ylabel(zvar)
    axa[1].set_xlabel(r"Cloud Fraction")
    axa[1].legend()

    x = qcl.isel({zvar:slice(1,None)}).mean(dim=[timevar,xvar,yvar])
    y = x.coords[zvar].values
    l1 = axa[2].plot(x, y, label=x.name)

    x = qcl_true.isel({zvar:slice(1,None)}).mean(dim=[timevar,xvar,yvar])
    y = x.coords[zvar].values
    l2 = axa[2].plot(x, y, label=x.name)

    axa[2].set_ylabel(zvar)
    axa[2].set_xlabel(r"$q_{cl}$")
    axa[2].legend()
    fig1.tight_layout()
    fig1.savefig(f'{plot_dir}Cloud_{filter_id}{plot_type}')
    plt.close()

fig2, axb = plt.subplots(1,2,figsize=(7,3))
axb[0].plot(cloud_fraction_true.values.flatten(),
            cloud_fraction.values.flatten(),
            marker='.', linestyle='none', markersize=0.01)
axb[0].set_xlim([0,Cmax[filter_no]])
axb[0].set_ylim([0,Cmax[filter_no]])
#axb[0].set_xlim([0,0.2])
#axb[0].set_ylim([0,0.2])
axb[0].set_xlabel('True Cloud Fraction')
axb[0].set_ylabel('Gaussian Cloud Fraction')


axb[1].plot(qcl_true.values.flatten(),
            qcl.values.flatten(),
            marker='.', linestyle='none', markersize=0.01)
axb[1].set_xlim([0,qmax[filter_no]])
axb[1].set_ylim([0,qmax[filter_no]])
#axb[1].set_xlim([0,0.0001])
#axb[1].set_ylim([0,0.0001])
axb[1].set_xlabel(r'True $q_{cl}$')
axb[1].set_ylabel(r'Gaussian $q_{cl}$')

fig2.tight_layout()
fig2.savefig(f'{plot_dir}Cloud_scatter_{filter_id}{plot_type}')
plt.close()


fig3, axc = plt.subplots(1,2,figsize=(7,3))
heatmap, xedges, yedges = np.histogram2d(cloud_fraction_true.values.flatten(),
                                         cloud_fraction.values.flatten(),
                                         bins=[Clev, Clev])
heatmap[0:1,0:1] =0
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

im1 = axc[0].imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm())
axc[0].set_xlabel('True Cloud Fraction')
axc[0].set_ylabel('Gaussian Cloud Fraction')
plt.colorbar(im1, ax=axc[0])

heatmap, xedges, yedges = np.histogram2d(qcl_true.values.flatten(),
                                         qcl.values.flatten(),
                                         bins=[qlev, qlev])
heatmap[0:1,0:1] =0
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

im2 = axc[1].imshow(heatmap.T, extent=extent, origin='lower', norm=LogNorm())
axc[1].set_xlabel(r'True $q_{cl}$')
axc[1].set_ylabel(r'Gaussian $q_{cl}$')
plt.colorbar(im2, ax=axc[1])
fig3.tight_layout()
fig3.savefig(f'{plot_dir}Cloud_heatmap_{filter_id}{plot_type}')
plt.close()
    #plt.show()

filtered_dataset.close()
