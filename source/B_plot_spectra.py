# -*- coding: utf-8 -*-
"""

@author: Todd Jones



                    BOMEX version


"""
print("importing...")
import os 
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import glob
import re
# Ignore divide-by-zero warnings that occur when converting between wavenumbers/frequency/wavelength
np.seterr(divide='ignore')

outpath='/home/users/toddj/20210127/'


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def forward(n):
    return 2*np.pi/n # l = 2*np.pi/n

def inverse(l):
    return 2*np.pi/l # n = 2*np.pi/l

def plot_spec_z(kpo,z,vname,flt,spec,time,k_max,k_mean,k_diss,xmax=0):
    fig, ax = plt.subplots(figsize=(6,5),dpi=200)
    cs=ax.contourf(kpo,z,spec*kpo,np.arange(0,0.155,0.005),cmap='RdPu',extend='max')
    cs.cmap.set_over('red')
    cbar = fig.colorbar(cs)
    cbar.ax.set_ylabel('$k*E$_'+vname+'(k)~[$m^2~s^{-2}$]')#r'$k*E_'+vname+'(k)~[m^2~s^{-2}]$')
    csl=ax.contour(kpo,z,spec*kpo,np.arange(0,0.155,0.01),cmap='Greys')
    cbar.add_lines(csl)
    plt.xscale('log')
    if xmax != 0:
        x1,x2,y1,y2 = plt.axis()
        plt.axis((x1,xmax,y1,y2))
    plt.title(flt+', '+vname+', Spectra (2D), time = %1.2f ' %(time))
    plt.xlabel(r'$k~[rad~m^{-1}]$')
    plt.ylabel('height $[m]$')
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel(r'wavelength $[m]$')
    ax.scatter(k_max,z,s=2,zorder=2,
               label='$\lambda_{max}$') # Return position of max, then index wavenumber
    ax.scatter(k_mean,z,color='tab:red',s=2,zorder=2,
               label='$\lambda_{mean}$')
    ax.scatter(k_diss,z,color='tab:olive',s=2,zorder=2,
               label='$\lambda_{diss}$')
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig(outpath+'spec_z_'+flt+'_'+vname+'_%1.2f.png' %(time))
    plt.close()

def plot_lamMax(times,z,kxarr,vname,flt):
    fig, ax = plt.subplots(figsize=(6,2),dpi=100)
    cs=ax.contourf(times,z,2*np.pi/kxarr.T,np.arange(0,6001,500),cmap='plasma')
    cbar = fig.colorbar(cs)
    cbar.ax.set_ylabel('[m]')
    plt.title(flt+', '+vname+', $\lambda_{max}$')
    plt.xlabel('time [hr]')
    plt.ylabel('height $[m]$')
    plt.tight_layout()
    plt.savefig(outpath+'lamMax_zt_'+flt+'_'+vname+'.png')
    plt.close()

def plot_lamMean(times,z,kmarr,vname,flt):
    fig, ax = plt.subplots(figsize=(6,2),dpi=100)
    cs=ax.contourf(times,z,2*np.pi/kmarr.T,np.arange(0,6001,100),cmap='plasma')
    cbar = fig.colorbar(cs)
    cbar.ax.set_ylabel('[m]')
    plt.title(flt+', '+vname+', $\lambda_{mean}$')
    plt.xlabel('time [hr]')
    plt.ylabel('height $[m]$')
    plt.tight_layout()
    plt.savefig(outpath+'lamMean_zt_'+flt+'_'+vname+'.png')
    plt.close()

def plot_lamDiss(times,z,kdarr,vname,flt):
    fig, ax = plt.subplots(figsize=(6,2),dpi=100)
    cs=ax.contourf(times,z,2*np.pi/kdarr.T,np.arange(0,6001,100),cmap='plasma')
    cbar = fig.colorbar(cs)
    cbar.ax.set_ylabel('[m]')
    plt.title(flt+', '+vname+', $\lambda_{diss}$')
    plt.xlabel('time [hr]')
    plt.ylabel('height $[m]$')
    plt.tight_layout()
    plt.savefig(outpath+'lamDiss_zt_'+flt+'_'+vname+'.png')
    plt.close()



#################################################################################################
#################################################################################################
#################################################################################################
def main():


    
    inpath='/gws/nopw/j04/paracon_rdg/users/toddj/updates_suite/'
    text_file = open(inpath+"sBlist.txt", "r")  # sBlist.txt contains a list of BOMEX cases
    lines = [line.rstrip('\n') for line in text_file.readlines()]
    print("beginning")
    collection=[]
    for lnc in lines:

        case_label = lnc
        case_pattern = 'BOMEX*test_plot_gn_00100_spectra.nc'
        shapevar = 'spec_2d_w'  # name of a spectra variable
        
        
        # Load the full datasets
        dsf = xr.open_mfdataset(natural_sort(glob.glob(inpath+lnc+'/diagnostic_files/spectra/*w_2D.nc')), concat_dim="time")
    
    
    
        # Prep some helpers
        #rnames = dsr['variable_names'][0,...].values.astype(str)
    
    # HERE, we need to watch the varshape variable and the flt label
    # AS a time saver, you can cancel the plot_spec_z call, which makes the loop take much longer.
        fnames = dsf['variable_names'][0,...].values.astype(str)
    
    
        wn = dsf['hwaven'][0,...].values
        wl = dsf['hwavel'][0,...].values
        freq = dsf['hfreq'].values
        z = dsf['z'].values
        times = dsf['time'].values/3600.
        
        varshape = dsf[shapevar].shape
        nt=varshape[0]
        nk=varshape[1]
        nz=varshape[2]
        
        kpo=wn
        if lnc == lines[0]:
            xmax=np.max(kpo)
            xmin=np.min(kpo)
        dk=kpo[1]-kpo[0]
        flt=case_label
        kxarr=np.zeros((nt,nz))
        kmarr=np.zeros((nt,nz))
        kdarr=np.zeros((nt,nz))
        for vname in fnames:
            print(vname)
            for tnc, time in enumerate(times):
                spec=dsf['spec_2d_'+vname][tnc,:,:].values.T  #[time, hreq, z] --> [z, hfreq]
                # print(spec.shape,kpo.shape)
                # at each height, record wavenumber where power*wavenember is maximized
                k_max=kpo[np.argmax(spec*kpo,axis=1)]
                # at each height, record wavenumber associated with mean power*wavenumber
                # requested as: Power-weighted mean wavelength as a function of level. 
                #               i.e. Integral k E(k) dk / Integral E(k) dk expressed as wavelength
                k_mean=np.trapz(spec*kpo, dx=dk, axis=1)/np.trapz(spec, dx=dk, axis=1)
                # at each height, record the dissipation wavenumber
                # BOUN-D-13-00082-14_Revision2.pdf, doi:10.1007/s10546-013-9881-3:
                # Bob Beare (2014)
                # The dissipation length scale can now potentially separate the relative effects of
                # numerical dissipation in different model configurations. For example, for a given
                # grid length, a more dissipative simulation will have a larger value of ld. It is
                # important to note that this dissipation length scale is due to the TKE spectrum
                # associated with the simulation, including implicit damping from the advection
                # scheme, and not just the subgrid viscosity
                k_diss=np.sqrt(np.trapz(spec*(kpo**2), dx=dk, axis=1)/np.trapz(spec, dx=dk, axis=1))

                # store the above profiles at each time for later plotting
                kxarr[tnc,:] = k_max
                kmarr[tnc,:] = k_mean
                kdarr[tnc,:] = k_diss
                print(tnc,end=' ')

                # Plot a contour of the spectra as a funtion of height at the last time.
                if tnc == nt-1:
                    plot_spec_z(kpo,z,vname,flt,spec,time,k_max,k_mean,k_diss,xmax=xmax)

            # endfor tnc loop over times

            # Plot data contours (time/height) at all times for this case
            plot_lamMax(times,z,kxarr,vname,flt)
            plot_lamMean(times,z,kmarr,vname,flt)
            plot_lamDiss(times,z,kdarr,vname,flt)

            # Calculate and plot Plot data averaged over a time range
            #  Spectra averaged first, then lengths are calculated.
            #  updates_suite BOMEX cases record every 10 minutes, 6 entries is an hour.
            aspec=np.mean(dsf['spec_2d_'+vname][-6*4:,:,:].values, axis=0).T
            akx=kpo[np.argmax(aspec*kpo,axis=1)]
            akm=np.trapz(aspec*kpo, dx=dk, axis=1)/np.trapz(aspec, dx=dk, axis=1)
            akd=np.sqrt(np.trapz(aspec*(kpo**2), dx=dk, axis=1)/np.trapz(aspec, dx=dk, axis=1))
            plot_spec_z(kpo,z,vname,flt+'_avg',aspec,20.24,akx,akm,akd,xmax=xmax)

            # Collect profile data to plot these all together later
            # ( case_label, mean_wavenum, height [m], max_wavenum, diss_wavenum)
            collection.append( (flt, akm, z, akx, akd) )

        # endfor vname loop over fnames(variable names)


    # Plot time-averaged profiles for all cases all together
    # \lambda_{mean}
    fig, ax = plt.subplots(figsize=(6,5),dpi=200)
    for data in collection:
        ax.scatter(data[1],data[2],s=2,label=data[0])
    plt.xscale('log')
    x1,x2,y1,y2 = plt.axis()
    plt.axis((xmin,xmax,y1,y2))
    plt.title('$\lambda_{mean}$, '+vname+', Spectra (2D), time = 20-24 h')
    plt.xlabel(r'$k~[rad~m^{-1}]$')
    plt.ylabel('height $[m]$')
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel(r'wavelength $[m]$')
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig(outpath+'spec_z_BOMEXallavgLamMean_'+vname+'.png')
    plt.close()

    #\lambda_{max}
    fig, ax = plt.subplots(figsize=(6,5),dpi=200)
    for data in collection:
        ax.scatter(data[3],data[2],s=2,label=data[0])
    plt.xscale('log')
    x1,x2,y1,y2 = plt.axis()
    plt.axis((xmin,xmax,y1,y2))
    plt.title('$\lambda_{max}$, '+vname+', Spectra (2D), time = 20-24 h')
    plt.xlabel(r'$k~[rad~m^{-1}]$')
    plt.ylabel('height $[m]$')
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel(r'wavelength $[m]$')
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig(outpath+'spec_z_BOMEXallavgLamMax_'+vname+'.png')
    plt.close()

    #\lambda_{diss}
    fig, ax = plt.subplots(figsize=(6,5),dpi=200)
    for data in collection:
        ax.scatter(data[4],data[2],s=2,label=data[0])
    plt.xscale('log')
    x1,x2,y1,y2 = plt.axis()
    plt.axis((xmin,xmax,y1,y2))
    plt.title('$\lambda_{diss}$, '+vname+', Spectra (2D), time = 20-24 h')
    plt.xlabel(r'$k~[rad~m^{-1}]$')
    plt.ylabel('height $[m]$')
    secax = ax.secondary_xaxis('top', functions=(forward, inverse))
    secax.set_xlabel(r'wavelength $[m]$')
    ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig(outpath+'spec_z_BOMEXallavgLamDiss_'+vname+'.png')
    plt.close()


    # endfor lnc loop over cases
# END main


if __name__ == "__main__":
    main()
