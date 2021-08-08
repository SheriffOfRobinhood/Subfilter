# -*- coding: utf-8 -*-
"""

  cloud_analysis.py

Created on Wed May 15 09:46:31 2019

@author: Peter Clark

"""
import os
import numpy as np
import xarray as xr

import matplotlib
import matplotlib.pyplot as plt
import scipy.special as sp
import math as m


from scipy.stats import skew
from scipy import ndimage

import dask

import subfilter.thermodynamics.thermodynamics_constants as thc
import subfilter.thermodynamics.thermodynamics as th

import subfilter.subfilter as sf
import subfilter.filters as filt

from subfilter.utils.string_utils import get_string_index
from subfilter.io.dataout import save_field
from subfilter.io.datain import get_data, get_pref, get_thref
from subfilter.io.MONC_utils import options_database



fname = 'cloud'
test_case = 0
calc_data = False

zlevel = 612.5

if test_case == 0:
    options = {
    #        'FFT_type': 'FFTconvolve',
    #        'FFT_type': 'FFT',
            'FFT_type': 'RFFT',
            'save_all': 'Yes',
            'th_ref': 300.0,
            'dx': 50.0,
            'dy': 50.0,
              }
    # indir = '/storage/silver/mft-scratch/dc915452/UG_project/BOMEX/m0050_g0128/'
    # indir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
    indir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'

#    odir = '/storage/silver/wxproc/xm904103/cloud/'
    odir = indir
    odir = odir + 'cloud/'
    file = 'diagnostics_3d_ts_21600.nc'
    ref_file = 'diagnostics_ts_21600.nc'
elif test_case == 1:
    options = {
    #        'FFT_type': 'FFTconvolve',
    #        'FFT_type': 'FFT',
            'FFT_type': 'RFFT',
            'save_all': 'Yes',
            'th_ref': 300.0,
            'dx': 5.0,
            'dy': 5.0,
              }
    indir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
    odir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
    odir = odir + 'test_dask_' + options['FFT_type']+'/'
    file = 'diagnostics_3d_ts_13200.nc'
    ref_file = None


# files_1D = glob.glob(indir+"diagnostics_ts_*.nc")

# dataset_1D = Dataset(files_1D[-1])

# thref= dataset_1D.variables['thref'][-1,...]
# pref= dataset_1D.variables['prefn'][-1,...]
# print(dataset_1D)
os.makedirs(odir, exist_ok = True)

#file = 'diagnostics_ts_18000.0.nc'
#ref_file = 'diagnostics_ts_18000.0.nc'

plot_dir = odir + 'plots/'
os.makedirs(plot_dir, exist_ok = True)

plot_type = '.png'
figshow = True
debug_label = False
def label_cyclic(mask) :
    """
    Function to label 3D objects taking account of cyclic boundary
    in x and y. Uses ndimage(label) as primary engine.

    Args:
        mask: 3D logical array with object mask (i.e. objects are
            contiguous True).

    Returns:
        Object identifiers::

            labs  : Integer array[nx,ny] of labels. -1 denotes unlabelled.
            nobjs : number of distinct objects. Labels range from 0 to nobjs-1.

    @author: Peter Clark

    """

#    print np.shape(mask)
    (nx, ny) = np.shape(mask)
    labels, nobjects = ndimage.label(mask)
    labels -=1
#    print 'labels', np.shape(labels)
    def relabel(labs, nobjs, i,j) :
#        if debug_label :
#            print('Setting label {:3d} to {:3d}'.format(j,i))
        lj = (labs == j)
        labs[lj] = i
        for k in range(j+1,nobjs) :
            lk = (labs == k)
            labs[lk] = k-1
#            if debug_label :
#                print('Setting label {:3d} to {:3d}'.format(k,k-1))
        nobjs -= 1
#        if debug_label : print('nobjects = {:d}'.format(nobjects))
        return labs, nobjs

    def find_objects_at_edge(minflag, x_or_y, n, labs, nobjs) :
#        for minflag in (True,False) :
        i = 0
        while i < (nobjs-2) :
            posi = np.where(labs == i)
#        print(np.shape(posi))
            if minflag :
                test1 = (np.min(posi[x_or_y][:]) == 0)
                border = '0'
            else:
                test1 = (np.max(posi[x_or_y][:]) == (n-1))
                border = 'n{}-1'.format(['x','y'][x_or_y])
            if test1 :
                if debug_label :
                    print('Object {:03d} on {}={} border?'.\
                          format(i,['x','y'][x_or_y],border))
                j = i+1
                while j < nobjs :
                    posj = np.where(labs == j)

                    if minflag :
                        test2 = (np.max(posj[x_or_y][:]) == (n-1))
                        border = 'n{}-1'.format(['x','y'][x_or_y])
                    else:
                        test2 = (np.min(posj[x_or_y][:]) == 0)
                        border = '0'

                    if test2 :
                        if debug_label :
                            print('Match Object {:03d} on {}={} border?'\
                                  .format(j,['x','y'][x_or_y],border))

                        if minflag :
                            ilist = np.where(posi[x_or_y][:] == 0)
                            jlist = np.where(posj[x_or_y][:] == (n-1))
                        else :
                            ilist = np.where(posi[x_or_y][:] == (n-1))
                            jlist = np.where(posj[x_or_y][:] == 0)

                        if np.size( np.intersect1d(posi[1-x_or_y][ilist], \
                                           posj[1-x_or_y][jlist]) ) :
#                            if np.size( np.intersect1d(posi[2][ilist], \
#                                           posj[2][jlist]) ) :

                           if debug_label :
                               print('Yes!',i,j)
                           labs, nobjs = relabel(labs, nobjs, i, j)
                    j += 1
            i += 1
        return labs, nobjs

    labels, nobjects = find_objects_at_edge(True,  0, nx, labels, nobjects)
    labels, nobjects = find_objects_at_edge(False, 0, nx, labels, nobjects)
    labels, nobjects = find_objects_at_edge(True,  1, ny, labels, nobjects)
    labels, nobjects = find_objects_at_edge(False, 1, ny, labels, nobjects)

    return labels, nobjects


def file_key(file):
    f1 = file.split('_')[-1]
    f2 = f1.split('.')[0]
    return float(f2)


# def read_multiple(files, var):
#     r = []
#     for f in files:
#         print(f)
#         dataset = Dataset(f)
#         d = dataset.variables[var]
#         r.append(d[...])

#     return np.vstack(r)

dask.config.set({"array.slicing.split_large_chunks": True})
dataset = xr.open_dataset(indir+file)
[itime, iix, iiy, iiz] = get_string_index(dataset.dims, ['time', 'x', 'y', 'z'])
timevar = list(dataset.dims)[itime]
xvar = list(dataset.dims)[iix]
yvar = list(dataset.dims)[iiy]
zvar = list(dataset.dims)[iiz]
max_ch = sf.subfilter_setup['chunk_size']

nch = np.min([int(dataset.dims[xvar]/(2**int(np.log(dataset.dims[xvar]
                                            *dataset.dims[yvar]
                                            *dataset.dims[zvar]
                                            /max_ch)/np.log(2)/2))),
              dataset.dims[xvar]])

print(f'nch={nch}')


npoints = dataset.dims[xvar]

dataset.close()

defn = 1
dataset = xr.open_dataset(indir+file, chunks={timevar: defn,
                                                xvar:nch, yvar:nch,
                                                'z':'auto', 'zn':'auto'})


#dataset = Dataset(indir+file, 'r') # Dataset is the class behavior to open the file
                         # and create an instance of the ncCDF4 class
ref_dataset=xr.open_dataset(indir + ref_file)

th_L = get_data(dataset, ref_dataset, 'th_L', options)
theta= get_data(dataset, ref_dataset, 'th', options)
q_cl = get_data(dataset, ref_dataset, 'q_cloud_liquid_mass', options)

q_t = get_data(dataset, ref_dataset, 'q_total', options)

w = get_data(dataset, ref_dataset, 'w', options)
w = w.interp({'z':th_L.coords['zn']}, kwargs={"fill_value": "extrapolate"})
w = w.drop_vars('z')

(pref, piref) = get_pref(dataset, ref_dataset,  options)

th_ref = get_thref(ref_dataset, options)
th_L_r = th_ref

p = 1.0E5 * (piref ** thc.rk)
T = theta * piref

#alpha_L = xr.DataArray(th.dqsatbydT(T_L_r.data, p.data), coords=th_L_r.coords)

s_N = q_t/th.qsat(T.data, p.data)-1

fig, ax = plt.subplots(1,1)

ax.hist2d(s_N.sel(zn=zlevel).data.flatten(), w.sel(zn=zlevel).data.flatten(),
           bins = [np.linspace(-15,5,101)*0.01, np.linspace(-0.5,1,101)] )

ax.set_xlabel(r'$s_N$')
ax.set_ylabel(r'$w$ m s$^{-1}$')
plt.savefig(plot_dir+f"w_v_q_t_raw"+plot_type)

plt.show()
# Q_cr = a_L*(q_tf_r - th.qsat(T_Lf_r.data, p.data))
# Q_cr = Q_cr.rename('Q_cr')
# QN = Q_cr / sigma_s
# QN = QN.rename('QN')

# q_cl_N = q_clf_r / sigma_s
# q_cl_N = q_cl_N.rename('q_cl_N')

# th_L_q_s = alpha_L_pi * th_Lf_s


# s_s = a_L * (q_tf_s - th_L_q_s)
# s_s = s_s.rename('s_s')

x = w.coords['x_p'].data
dx = x[1] - x[0]
y = w.coords['y_p'].data
dy = y[1] - y[0]

npts = len(w.coords['x_p'].data)*len(w.coords['y_p'].data)

domain_area = npts * dx * dy

#print(np.max(alpha_L),np.min(alpha_L))
#a_L = 1.0 / (1.0 + L_over_cp * alpha_L)
#print(np.shape(a_L))
#print(np.max(a_L),np.min(a_L))


q_cls = q_cl.sel(zn=zlevel)
s_Ns = s_N.sel(zn=zlevel)
ws = w.sel(zn=zlevel)
q_clt = 1E-5
wg = ws.where(q_cls >= q_clt)
wl = ws.where(q_cls < q_clt)

s_Ng = s_Ns.where(ws >  0)
s_Nl = s_Ns.where(ws <= 0)

fs = (8,12)
q_cls.plot(x='x_p', y='y_p', col='time',col_wrap=5)
plt.savefig(plot_dir+f"q_cl_raw"+plot_type)

ws.plot(x='x_p', y='y_p', col='time',col_wrap=5)
plt.savefig(plot_dir+f"w_raw"+plot_type)
plt.show()

wbins=np.linspace(-1, 2, 101)
fig0, axa = plt.subplots(3,1, sharex=True, figsize=fs)
ws.plot.hist(bins=wbins, ax=axa[0])
axa[0].set_title(r'All $w$')
wl.plot.hist(bins=wbins, ax=axa[1])
axa[1].set_title(r'$w$, $q_{cl} < 10^{-5}$ kg/kg ')
wg.plot.hist(bins=wbins, ax=axa[2])
axa[2].set_title(r'$w$, $q_{cl} \geq 10^{-5}$ kg/kg ')
plt.tight_layout()
plt.savefig(plot_dir+f"hist_w_raw"+plot_type)
plt.show()

def dist_stats(v, vname):

    print(vname)
    print(f'mean: {v.mean().values:0.2f}')
    print(f'std: {v.std().values:0.2f}')
    sk = skew(v, axis=None, nan_policy='omit', bias=False)
    print(f'sk: {sk:0.2f}')
    return

dist_stats(wg, 'w all ')

prop_wg = (np.size(wg)-np.isnan(wg).sum().values)/np.size(ws)*100.0
print(f'prop: {prop_wg:2.1f}')
dist_stats(wg, 'w cloudy ')
# print(f'mean: {wg.mean().values:0.2f}')
# print(f'std: {wg.std().values:0.2f}')
# sk = skew(wg, axis=None, nan_policy='omit', bias=False)
# print(f'sk: {sk:0.2f}')

prop_wl = (np.size(wl)-np.isnan(wl).sum().values)/np.size(ws)*100.0
print(f'prop: {prop_wl:2.1f}')
dist_stats(wl, 'w not cloudy ')
# print(f'mean: {wl.mean().values:0.2f}')
# print(f'std: {wl.std().values:0.2f}')
# sk = skew(wl, axis=None,nan_policy='omit', bias=False)
# print(f'sk: {sk:0.2f}')

sbins = np.linspace(-15,5,101)*0.01
fig0a, axa = plt.subplots(3,1, sharex=True, figsize=fs)
s_Ns.plot.hist(bins=sbins, ax=axa[0])
axa[0].set_title(r'All $s_N$')
s_Nl.plot.hist(bins=sbins, ax=axa[1])
axa[1].set_title(r'$s_N$, $w \leq 0$ m s$^{-1}$ ')
s_Ng.plot.hist(bins=sbins, ax=axa[2])
axa[2].set_title(r'$s_N$, $w > 0$ m s$^{-1}$ ')
plt.tight_layout()
plt.savefig(plot_dir+f"hist_s_N_raw"+plot_type)
plt.show()

dist_stats(s_Ns, 's_N all ')
prop_s_Ng = (np.size(s_Ng)-np.isnan(s_Ng).sum().values)/np.size(s_Ns)*100.0
print(f'prop: {prop_s_Ng:2.1f}')
dist_stats(s_Ng, 's_N w >0 ')
prop_s_Nl = (np.size(s_Nl)-np.isnan(s_Nl).sum().values)/np.size(s_Ns)*100.0
print(f'prop: {prop_s_Nl:2.1f}')
dist_stats(s_Nl, 's_N w<=0 ')

# print(f'mean: {s_Ns.mean().values:0.2f}')
# print(f'std: {s_Ns.std().values:0.2f}')
# sk = skew(s_Ns, axis=None, nan_policy='omit', bias=False)
# print(f'sk: {sk:0.2f}')

# print(f'mean: {s_Ng.mean().values:0.2f}')
# print(f'std: {s_Ng.std().values:0.2f}')
# sk = skew(s_Ng, axis=None, nan_policy='omit', bias=False)
# print(f'sk: {sk:0.2f}')

# print(f'mean: {s_Nl.mean().values:0.2f}')
# print(f'std: {s_Nl.std().values:0.2f}')
# sk = skew(s_Nl, axis=None,nan_policy='omit', bias=False)
# print(f'sk: {sk:0.2f}')

################################### Correllation #############################
def correlate( a, b, cond=None):
    if cond is None:
        af = a.data.flatten().compute()
        bf = b.data.flatten().compute()
    else:
        af = a.where(cond).data.flatten().compute()
        bf = b.where(cond).data.flatten().compute()
        af = af[~np.isnan(af)]
        bf = bf[~np.isnan(bf)]

    l = len(af)
    ma = af.mean()
    mb = bf.mean()
    ap = af - ma
    bp = bf - mb
    va = np.mean(ap*ap)
    vb = np.mean(bp*bp)
    cov = np.mean(ap*bp)

    return(l, ma, va, mb, vb, cov)

(l, mw, vw, ms, vs, cov_ws) = correlate(ws, s_Ns)
print(l, mw, vw, ms, vs, cov_ws, cov_ws/np.sqrt(vw * vs))

(lg, mwg, vwg, msg, vsg, cov_wsg) = correlate(ws, s_Ns, q_cls >= q_clt)
print(lg, mwg, vwg, msg, vsg, cov_wsg, cov_wsg/np.sqrt(vwg * vsg))
(ll, mwl, vwl, msl, vsl, cov_wsl) = correlate(ws, s_Ns, q_cls < q_clt)
print(ll, mwl, vwl, msl, vsl, cov_wsl, cov_wsl/np.sqrt(vwl * vsl))


labels, nobjects = label_cyclic(q_cls[0,:,:] >= q_clt)

plt.imshow(labels)
l_d = np.sqrt(domain_area / nobjects)
l_c = np.sqrt(np.sum(labels >= 0) * dx * dy /nobjects)

print(f'nobjects: {nobjects}')
print(f'l_c: {l_c:3.0f}')
print(f'l_d: {l_d:4.0f}')
print(r'$\sigma$:'+f'{l_c*l_c/(l_d*l_d)*100:2.2f}')

# #    sigma_list = [0.5,0.2]
sigma_list = [125.0, 250.0, 500.0, 1000.0]
#sigma_list = [125.0]
width = -1
filter_name = 'gaussian'

opgrid = 'p'
var_list = ["w", \
        "th", \
        "th_v", \
        "th_L", \
        "q_vapour", \
        "q_cloud_liquid_mass", \
        "q_total", \
        ]

var_list2 = [
     ["w","w"],
     ["w","th_L"],
     ["w","q_vapour"],
     ["w","q_cloud_liquid_mass"],
     ["w","q_total"],
     ["th_L","th_L"],
     ["th_L","q_total"],
     ["q_total","q_total"],
     ["th_L","q_vapour"],
     ["th_L","q_cloud_liquid_mass"],
   ]


derived_data, exists = \
    sf.setup_derived_data_file( indir+file, odir, fname,
                               options, override=True)
# print("Variables in derived dataset.")
# print(derived_data['ds'].variables)

# Now create list of filter definitions.

filter_list = list([])

for i,sigma in enumerate(sigma_list):
    if filter_name == 'gaussian':
        filter_id = 'filter_ga{:02d}'.format(i)
        twod_filter = filt.Filter(filter_id,
                                  filter_name,
                                  npoints=npoints,
                                  sigma=sigma, width=width,
                                  delta_x=dx)
    elif filter_name == 'wave_cutoff':
        filter_id = 'filter_wc{:02d}'.format(i)
        twod_filter = filt.Filter(filter_id,
                                  filter_name, wavenumber=np.pi/(2*sigma),
                                  npoints=npoints,
                                  width=width,
                                  delta_x=dx)
    elif filter_name == 'running_mean':
        filter_id = 'filter_rm{:02d}'.format(i)
        width = int(np.round( sigma/dx * np.pi * 2.0 / 3.0)+1)
        twod_filter = filt.Filter(filter_id,
                                  filter_name,
                                  npoints=npoints,
                                  width=width,
                                  delta_x=dx)

    print(twod_filter)
    filter_list.append(twod_filter)

for twod_filter in filter_list:

    print("Processing using filter: ")
    print(twod_filter)

    if calc_data:

        filtered_data, exists = \
            sf.setup_filtered_data_file( indir+file, odir, fname,
                                   options, twod_filter, override=True)
    else:
        filtered_data, exists = \
            sf.setup_filtered_data_file( indir+file, odir, fname,
                                   options, twod_filter, override=False)
#        print("Variables in filtered dataset.")
#        print(filtered_data['ds'].variables)

    if not exists or calc_data:
        field_list =sf.filter_variable_list(dataset, ref_dataset,
                                            derived_data, filtered_data,
                                            options,
                                            twod_filter, var_list=var_list,
                                            grid = opgrid)

        quad_field_list =sf.filter_variable_pair_list(dataset, ref_dataset,
                                     derived_data, filtered_data,
                                     options,
                                     twod_filter, var_list=var_list2,
                                     grid = opgrid)

        filtered_data['ds'] = xr.open_dataset(filtered_data['file'])
    else :
        print('Filtered data file exists' )



    wf_r = filtered_data['ds']['f(w_on_p)_r']
    wf_s = filtered_data['ds']['f(w_on_p)_s']

    q_clf_r = filtered_data['ds']['f(q_cloud_liquid_mass_on_p)_r']
    q_clf_s = filtered_data['ds']['f(q_cloud_liquid_mass_on_p)_s']

    q_tf_r = filtered_data['ds']['f(q_total_on_p)_r']
    q_tf_s = filtered_data['ds']['f(q_total_on_p)_s']

    th_Lf_r = filtered_data['ds']['f(th_L_on_p)_r']
    th_Lf_s = filtered_data['ds']['f(th_L_on_p)_s']

    wf_r_z = wf_r.sel(zn=zlevel)
    wf_s_z = wf_s.sel(zn=zlevel)

    q_clf_r_z = q_clf_r.sel(zn=zlevel)
    q_clf_s_z = q_clf_s.sel(zn=zlevel)

    wfg = wf_s_z.where(q_clf_r_z >= q_clt)
    wfl = wf_s_z.where(q_clf_r_z <  q_clt)

    q_clf_r_z.plot(x='x_p', y='y_p', col='time',col_wrap=5)
    plt.savefig(plot_dir+f"q_cl_r_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)

    q_clf_s_z.plot(x='x_p', y='y_p', col='time',col_wrap=5)
    plt.savefig(plot_dir+f"q_cl_s_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)

    wf_r_z.plot(x='x_p', y='y_p', col='time',col_wrap=5)
    plt.savefig(plot_dir+f"wr_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)
    wf_s_z.plot(x='x_p', y='y_p', col='time',col_wrap=5)
    plt.savefig(plot_dir+f"ws_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)
    plt.show()

    fig1, axa = plt.subplots(3,1, sharex=True, figsize=fs)#,figsize=(10,12)
    wf_s_z.plot.hist(bins=wbins, ax=axa[0])
    axa[0].set_title(r'All $w$')
    wfl.plot.hist(bins=wbins, ax=axa[1])
    axa[1].set_title(r'$w$, $q_{cl} < 10^{-5}$ kg/kg ')
    wfg.plot.hist(bins=wbins, ax=axa[2])
    axa[2].set_title(r'$w$, $q_{cl} \geq 10^{-5}$ kg/kg ')
    plt.tight_layout()
    plt.savefig(plot_dir+f"hist_w_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)
    plt.show()

    th_L_th_L = filtered_data['ds']['s(th_L,th_L)_on_p']
    th_L_q_t = filtered_data['ds']['s(th_L,q_total)_on_p']
    q_t_q_t = filtered_data['ds']['s(q_total,q_total)_on_p']

    T_Lf_r = th_Lf_r * piref
    alpha_L = xr.DataArray(th.dqsatbydT(T_Lf_r.data, p.data), coords=th_Lf_r.coords)
    a_L = 1.0 / (1.0 + thc.L_over_cp * alpha_L)
    alpha_L_pi = alpha_L * piref

    var_s = a_L * a_L * (q_t_q_t - 2 * alpha_L_pi * th_L_q_t + \
                          alpha_L_pi * alpha_L_pi * th_L_th_L)


    sigma_s = np.sqrt(var_s)
    print(np.max(sigma_s),np.min(sigma_s))

    Q_cr = a_L*(q_tf_r - th.qsat(T_Lf_r.data, p.data))
    Q_cr = Q_cr.rename('Q_cr')
    QN = Q_cr / sigma_s
    QN = QN.rename('QN')

    C_MY=0.5*(1+sp.erf(QN/m.sqrt(2)))
    q_cl_MY = C_MY * QN + np.exp(-(QN**2)/2)/m.sqrt(2*m.pi)
    q_cl_MY_z = q_cl_MY.sel(zn=zlevel)

    q_cl_N = q_clf_r / sigma_s
    q_cl_N = q_cl_N.rename('q_cl_N')

    th_L_q_s = alpha_L_pi * th_Lf_s


    s_s = a_L * (q_tf_s - th_L_q_s)
    s_s = s_s.rename('s_s')
    s_s_N = s_s / sigma_s
    s_s_N = s_s_N.rename('s_s_N')

    s_s_N_z = s_s_N.sel(zn=zlevel)
    QN_z = QN.sel(zn=zlevel)

    s_N_z = QN_z + s_s_N_z
    s_N_z = s_N_z.rename('s_N')

    sbins = np.linspace(-20,4,241)

    s_N_z.plot(x='x_p', y='y_p', col='time',col_wrap=5,levels = sbins)
    plt.savefig(plot_dir+f"s_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)

    QN_z.plot(x='x_p', y='y_p', col='time',col_wrap=5,levels = sbins)
    plt.savefig(plot_dir+f"sr_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)
    s_s_N_z.plot(x='x_p', y='y_p', col='time',col_wrap=5)
    plt.savefig(plot_dir+f"ss_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)

    sg = s_N_z.where(wf_r_z > 0)
    sl = s_N_z.where(wf_r_z <=0)
    fig2, axa = plt.subplots(3,1, sharex=True, figsize=fs)#,figsize=(10,12)
    s_N_z.plot.hist(bins=sbins, ax=axa[0])
    axa[0].set_title(r'All $s$ ')
    sl.plot.hist(bins=sbins, ax=axa[1])
    axa[1].set_title(r'$s$, $w \leq$ 0 m s$^{-1}$ ')
    sg.plot.hist(bins=sbins, ax=axa[2])
    axa[2].set_title(r'$s$, $w >$ 0 m s$^{-1}$ ')
    plt.tight_layout()
    plt.savefig(plot_dir+f"hist_s_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)
    plt.show()


    wfg1 = wf_s_z.where(q_cl_MY_z >= q_clt)
    wfl1 = wf_s_z.where(q_cl_MY_z <  q_clt)

    fig3, axa = plt.subplots(3,1, sharex=True, figsize=fs)#,figsize=(10,12)
    wf_s_z.plot.hist(bins=wbins, ax=axa[0])
    axa[0].set_title(r'All $w$')
    wfl1.plot.hist(bins=wbins, ax=axa[1])
    axa[1].set_title(r'$w$, $q_{cl_{MY}} < 10^{-5}$ kg/kg ')
    wfg1.plot.hist(bins=wbins, ax=axa[2])
    axa[2].set_title(r'$w$, $q_{cl_{MY}} \geq 10^{-5}$ kg/kg ')
    plt.tight_layout()
    plt.savefig(plot_dir+f"hist_w_filt_MY_{twod_filter.attributes['sigma']:04.0f}"+plot_type)
    plt.show()

    fig4, ax = plt.subplots(1,1)

    ax.hist2d(s_N_z.data.flatten(), wf_r_z.data.flatten(),
               bins = [np.linspace(-15,5,101), np.linspace(-0.5,1,151)] )

    ax.set_xlabel(r'$s_N$')
    ax.set_ylabel(r'$w$ m s$^{-1}$')
    plt.savefig(plot_dir+f"w_v_q_t_filt_{twod_filter.attributes['sigma']:04.0f}"+plot_type)

    plt.show()

    filtered_data['ds'].close()




# C_MY=0.5*(1+sp.erf(QN_c/m.sqrt(2)))
# q_cl_MY = C_MY * QN_c + np.exp(-(QN_c**2)/2)/m.sqrt(2*m.pi)


