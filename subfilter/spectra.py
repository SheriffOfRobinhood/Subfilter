
# -*- coding: utf-8 -*-
"""
    This program evaluates MONC fields with 4 dimensions (vertical, x, y, time) to produce
    horizontal power spectra at each time and vertical level written to new netcdf files.

    BY DEFAULT, each variable contained in the input file's xarray.Dataset "Data variables"
    has its horizontal power spectra evaluated.  These are all placed in the same
    output file.

    They can alternatively be placed in a list in the user settings section.

    Several options can influence the form of the final result.

    Assumes the horizontal grid dimensions are the same for each variable being analysed.
    Assumes the horizontal dimensions are named 'x' and 'y'.
    Assumes the vertical dimension is the only dimension with a 'z' in its name, but it can be
    either 'z' or 'zn'.
    The time dimension name is identified by a user-supplied string, currently: 'time'.

    "Durran" calculation based on Durran et al. (2017): https://doi.org/10.1175/MWR-D-17-0056.1

    User must supply::

       dir:    input directory (slash-agnostic)
       file:   input file
                 Suggest switching to argument input (see below)
       outtag: output file tag (appended to input file name)
                 Creates 'spectra/' directory within the given dir
       dx:     x-direction grid spacing [m]
       dy:     y-direction grid spacing [m]

    @author: Todd Jones
    @modified: Peter Clark
"""

import sys
import numpy as np
import xarray as xr
from scipy import ndimage      # Required for the radial summations
from dask.diagnostics import ProgressBar
import dask.array as da
from .utils.string_utils import get_string_index
from .io.datain import configure_model_resolution
from subfilter import executing_on_cluster

import yaml


time_dim_always_contains='time'

def spectra_options(config_file=None):
# Set computation options
    update_config = None
    options = {
           'spec_1D': True,            # Compute 2D spectra as average of
                                       # 1D spectra in both directions
           'spec_2D': True,            # Compute 2D spectra via 2D-fft
           'spec_method': 'durran',    # [Durran, ndimage] Use Durran method
                                       # (which actually also
#           'spec_method': 'ndimage',    # [Durran, ndimage] Use Durran method
                                       # (which actually also uses ndimage),
                                       # or faster, less accurate ndimage
                                       # method
           'spec_compensation': True,  # With spec_method: 'durran', use
                                       # Durran/Tanguay method to compensate
                                       # for systematic noise in the annular
                                       # summation (spectra does not preserve
                                       #energy when enabled)
           'spec_restrict': True,      # With spec_method: 'durran',
                                       # restrict the spec_2d result to values
                                       # below the Nyquist frequency
           'override': True,           # Overwrite output file if it exists.
             }

    if config_file is not None:
        with open(config_file) as c:
            update_config = yaml.load(c, Loader = yaml.FullLoader)

        options.update(update_config['options'])

    return options, update_config


def spectra_variable_list(ds, derived_dataset, options, var_list=None):
    """
    Compute and store 1D forms of 2D power spectra.

    Parameters
    ----------
    ds : xarray dataset
        Input data
    derived_dataset : dict
        Output data.    'ds':xarray dataset 'file': Path to output file.
    options : dict
        options for spectral calculations
    var_list : list[str], optional
        List of variable names in ds to derive spectra for. The default is None.

    Returns
    -------
    dso : xarray dataset
        Output data.
    outfile : str

    """

    # map for 2D spectra

    if var_list is None:
        var_list = list(ds.data_vars.keys())

    # Check if KE spectra is requested
#not yet implemented    if options['do_ke']: var_list.append('ke')

    # Get model resolution values
    dx, dy, options = configure_model_resolution(ds, options)

#    if 'dx' in ds.attrs:
#        dx = ds.attrs['dx']
#        dy = ds.attrs['dy']
#    else:
#        od = options_database(ds)
#        if od is None:
#            dx = options['dx']
#            dy = options['dy']
#            print("in spectra, NOT using options_database")
#        else:
#            dx = float(od['dxx'])
#            dy = float(od['dyy'])
#            print("in spectra, IS, IN FACT, using options_database")

    dso = derived_dataset['ds']
    outfile = derived_dataset['file']
    kmap = None

    # Loop over var_list
    for vname in var_list:

        if vname == 'options_database':
            print(f'Ignoring variable {vname}.')
            continue
        if not (vname in ds):
            print(f'[WARN] {vname} not present. skipping...')
            continue

        if options['spec_1D']:
            spectrum_ave_1D(ds, dso, vname, outfile, options, dx, dy)

        if options['spec_2D']:
            kmap = spectrum_ave_1D_radial(ds, dso, vname, outfile,
                                           options, dx, dy, kmap=kmap)

    ds.close()

    return dso

def spectrum_ave_1D(ds, dso, vname, outfile, options, dx, dy):
    """
    Compute averaged 1D spectra in x and y directions separately.
    Use real fft of anomalies, Durran et al. (2017), Eq. 13,
    and average over y results over x direction (and vice versa),
    handling Nyquist (Kr_delta).


    Parameters
    ----------
    ds : xarray Dataset
        Input data.
    dso :  xarray Dataset
        Output data.
    vname : str
        Variable name (in ds).
    outfile : str
        Path to output file.
    options : dict
        Options controlling spectrum calculations.
    dx : float
        x grid spacing.
    dy : float
        y grid spacing.

    Returns
    -------
    None.

    """
    # Ignore divide-by-zero warnings that occur when converting between wavenumbers/frequency/wavelength
    np.seterr(divide='ignore')

    print(f"Working on spectra for: {vname}.")
    # Load variable
#    var = ds[vname].load()
    var = ds[vname]

    if np.size(var.shape) != 4:
        print("Expecting 4-dimensional data for ", vname)
        print("FAIL - shape")
        sys.exit()
        # FOR OTHER KINDS OF FILES WITH MANY DIFFERENTLY STRUCTURED
        #   FIELDS, YOU MIGHT USE continue TO SKIP THOSE VARIABLES
        #   RATHER THAN sys.exit()

    # Store some parameters
    # Coordinate positions
    (tx, xx, yx, zx) = get_string_index(var.dims,
                        [time_dim_always_contains, 'x', 'y', 'z'])

    # Time and height names
    tname = var.dims[tx]
    var   = var.rename({tname: 'time'})
    tname = "time"
    times = var[tname].values
    xname = var.dims[xx]
    yname = var.dims[yx]
    zname = var.dims[zx]
    z     = ds[zname].values

    # Coordinate sizes
    ny = var.sizes[yname]
    nx = var.sizes[xname]
    nt = var.sizes[tname]
    nz = var.sizes[zname]

    # Horizontal domain lengths
    L_x = dx*nx
    L_y = dy*ny

    print(var)

    #_____________________________________________________________________________
    # Handle averaged 1D spectra option
    # Use real fft of anomalies, Durran et al. (2017), Eq. 13,
    # and average over y results over x direction (and vice versa), handling Nyquist (Kr_delta)
    #
    #  This condition could be written as a function, given:
    #    Input:
    #          yx, dy, ny (could derive)
    #          xx, dx, nx (could derive)
    #          var (as a 4D xarray DataArray, assuming data structured [time,...horiz-dims...,vertical])
    #          vname (could derive)
    #          zname
    #          times
    #    Returns:
    #          xarray Dataset containing
    #              spectral density in y
    #              spectral density in x
    #              wavenumbers, wavelengths, frequencies
    #
    #
    #    Alternatively, we could make a function to only return the spectra and deal with Dataset assembly elsewhere,
    #      in which case this would require very little input.
    #_____________________________________________________________________________
    print("Calculating 1D spectra")

    # Spectra in y-direction (mean removed)

    field = (var - var.mean(dim=(yname))).data.rechunk(chunks={yx:ny})
    temp2 = np.fft.rfft(field, axis=yx)
    # Average across all x
    ydir = np.mean((temp2*np.conj(temp2)).real, axis = xx) * (dy/(2*np.pi*ny))

    # Halve the Nyqist value
    # Now assuming that data is structured ["time","yfreq",zname]
    ydir[:,-1,:] /= 2

    # Repeat for x-direction
    field = (var - var.mean(dim=(xname))).data.rechunk(chunks={xx:nx})
    temp2 = np.fft.rfft(field, axis=xx)
    # Average across all x
    xdir = np.mean((temp2*np.conj(temp2)).real, axis = yx) * (dx/(2*np.pi*nx))
    xdir[:,-1,:] /= 2

    # Useful plotting values
    fx = np.fft.rfftfreq(nx, d = dx)     # frequencies (1/m) [n/L]
    fy = np.fft.rfftfreq(ny, d = dy)
    xwaven = fx*2*np.pi                  # wavenumbers (radians/m) [n2pi/L]
    ywaven = fy*2*np.pi
    xwavel = 2*np.pi/xwaven  # 1/fx      # wavelength (m)
    ywavel = 2*np.pi/ywaven  # 1/fy

    # Assemble the 1D Dataset
    ds1 = xr.Dataset(data_vars = {"spec_ydir_"+vname:(["time","yfreq",zname],ydir),
                                  "spec_xdir_"+vname:(["time","xfreq",zname],xdir),
                                  "ywaven":(['yfreq'],ywaven),
                                  "xwaven":(['xfreq'],xwaven),
                                  "ywavel":(['yfreq'],ywavel),
                                  "xwavel":(['xfreq'],xwavel),   },
                     coords    = {"time": (["time"],times),
                                  zname: ([zname],z),
                                  "yfreq": (["yfreq"],fy),
                                  "xfreq": (["xfreq"],fx) }  )
    ds1['spec_ydir_'+vname].attrs = {"units": "m * base_unit^2",
                              "long_name": "spectral density, y-direction"}
    ds1['spec_xdir_'+vname].attrs = {"units": "m * base_unit^2",
                              "long_name": "spectral density, x-direction"}
    ds1['yfreq'].attrs = {"units": "m-1",
                         "long_name": "spectral frequency, y-direction"}
    ds1['xfreq'].attrs = {"units": "m-1",
                         "long_name": "spectral frequency, x-direction"}
    ds1['time'].attrs = {"units": "s",
                         "long_name": "model time since start"}
    ds1[zname].attrs = {"units": "m",
                        "long_name": "height above surface"}
    ds1['ywaven'].attrs = {"units": "rad m-1",
                           "long_name": "spectral wavenumber, y-direction"}
    ds1['xwaven'].attrs = {"units": "rad m-1",
                           "long_name": "spectral wavenumber, x-direction"}
    ds1['ywavel'].attrs = {"units": "m",
                           "long_name": "spectral wavelength, y-direction"}
    ds1['xwavel'].attrs = {"units": "m",
                           "long_name": "spectral wavelength, x-direction"}

        # Append new vname spectra to the output Dataset
#    dso = xr.merge([ds1, ds2])
    d = ds1.to_netcdf(path = outfile, unlimited_dims="time", mode='a',
                      compute=False)
    if executing_on_cluster:
        results = d.compute()
    else:
        with ProgressBar():
            results = d.compute()
    ds1.close()

    return

#===================================================================
# Get PSD 1D (total radial power spectrum)
# for spec_method: ndimage
#===================================================================
def GetPSD1D(psd2D, k):
    """
    Get PSD 1D (total radial power spectrum)
    For use with option spec_method: ndimage

    Args:
        psd2D    : 2D numpy array containing 2D spectra values

    Returns:
        psd1D    : 1D numpy array containing 1D spectra ordered from
                   wavenumber 0 to highest wavenumber of shortest
                   dimension

    @author:  https://gist.github.com/TangibitStudios/47beaf24690329ac7fecddde70835ce9

    """
    h  = psd2D.shape[0]
    w  = psd2D.shape[1]
    wc = w//2
    hc = h//2

    # create an array of integer radial distances from the center
    Y, X = np.ogrid[0:h, 0:w]
    r    = np.hypot(X - wc, Y - hc).astype(np.int)

    # SUM all psd2D pixels with label 'r' for 0<=r<=wc
    # NOTE: this will miss power contributions in 'corners' r>wc
    psd2D = np.fft.fftshift(psd2D)
    psd1D = ndimage.sum(psd2D, r, index=np.arange(0, wc))

    return psd1D
#===================================================================


def prepare_map(fkx, fky, kp, dkh, Nmax):
    print("Preparing map")
    nx = len(fkx)
    ny = len(fky)
    # Prepare wavnumber grid (rmap) based on kx and ky
    # To be used in radial sum and/or applying the Tanguay/Durran correction term
    gkx=np.tile(fkx, (ny, 1))                      # x wavenumber array, repeated ny times
    gky=np.tile(np.array([fky]).T, (1, nx))        # y wavenumber array, repeated nx times
    rmap=np.sqrt(gkx**2 + gky**2,dtype='float64')  # wavenumber grid
    rlab=(rmap*0).astype(np.int)                   # grid of labels denoting index of wavenumbers (kp)
    kcount=kp*0                                    # to hold count of points in kp; sum(kcount)=(nx*ny)-1
    kpbar=kp*0                                     # to hold mean of wavenumber grid values at kp
    rindex=np.arange(0,Nmax-1)                     # list of labels to sum (all - to start with...)

    for knc,kpl in enumerate(kp):
        keep=((kpl-dkh/2 <= rmap) & (rmap < kpl+dkh/2))  # see eqn 18
        rlab[keep] = knc
        kcount[knc] = np.count_nonzero(keep)
        kpbar[knc] = np.mean(rmap[keep])

    kmap = {
            'rlab'   : rlab,
            'rindex' : rindex,
            'kcount' : kcount,
            'kpbar'  : kpbar,
           }

    return kmap


def rad_ave_with_comp(Ek, rlab, index, norm, comp):
                # Sum points
    Ekp = norm * ndimage.sum(Ek, rlab, index=index) * comp

    return Ekp

def rad_ave_without_comp(Ek, rlab, index, norm):
                # Sum points
    Ekp = norm * ndimage.sum(Ek, rlab, index=index)

    return Ekp

def spectrum_ave_1D_radial(ds, dso, vname, outfile, options, dx, dy,
                           kmap=None):
    """
    Compute averaged 2D spectra averaged to 1D.
    Use real fft of anomalies, Durran et al. (2017), Eq. 13,
    and average over y results over x direction (and vice versa),
    handling Nyquist (Kr_delta).


    Parameters
    ----------
    ds : xarray Dataset
        Input data.
    dso :  xarray Dataset
        Output data.
    vname : str
        Variable name (in ds).
    outfile : str
        Path to output file.
    options : dict
        Options controlling spectrum calculations.
    dx : float
        x grid spacing.
    dy : float
        y grid spacing.
    kmap : dict, optional
        Previously computed mapping from radial k to 2D grid. The default is None.

    Returns
    -------
    kmap : dict
        Previously computed mapping from radial k to 2D grid.
    """


    # Ignore divide-by-zero warnings that occur when converting between wavenumbers/frequency/wavelength
    np.seterr(divide='ignore')


    print("Working on spectra for: ", vname)
    # Load variable
 #   var = ds[vname].load()
    var = ds[vname]

    if np.size(var.shape) != 4:
        print("Expecting 4-dimensional data for ", vname)
        print("FAIL - shape")
        sys.exit()
        # FOR OTHER KINDS OF FILES WITH MANY DIFFERENTLY STRUCTURED
        #   FIELDS, YOU MIGHT USE continue TO SKIP THOSE VARIABLES
        #   RATHER THAN sys.exit()

    # Store some parameters
    # Coordinate positions
    (tx, xx, yx, zx) = get_string_index(var.dims,
                        [time_dim_always_contains, 'x', 'y', 'z'])

    # Time and height names
    tname = var.dims[tx]
    var   = var.rename({tname: 'time'})
    tname = "time"
    times = var[tname].values
    xname = var.dims[xx]
    yname = var.dims[yx]
    zname = var.dims[zx]
    z     = ds[zname].values

    # Coordinate sizes
    ny = var.sizes[yname]
    nx = var.sizes[xname]
    nt = var.sizes[tname]
    nz = var.sizes[zname]

    # Horizontal domain lengths
    L_x = dx*nx
    L_y = dy*ny

    print(var)
    #_____________________________________________________________________________
    # Handle averaged 2D spectra option (to 1D horizontal)
    # NOTE: Wavenumber zero is neither considered nor reported - mean removed from data.
    #
    #  This condition could be written as a function, given:
    #    Input:
    #          prepare_map logical
    #          yx, dy, ny (could derive), L_y (could derive)
    #          xx, dx, nx (could derive), L_x (could derive)
    #          var (as a 4D xarray DataArray, assuming data structured [time,...horiz-dims...,vertical])
    #          vname (could derive)
    #          zname
    #          times
    #    Returns:
    #          xarray Dataset containing
    #              1D horizontal spectral density
    #              1D horizontal wavenumbers, wavelengths, frequencies
    #
    #
    #    Alternatively, we could make a function to only return the spectra and deal with Dataset assembly elsewhere,
    #      in which case this would require very little input.
    #_____________________________________________________________________________
    print("Calculating 2D spectra")

    # Set up helpful parameters

    # NOTE: Wavenumber zero is neither considered nor reported - mean removed from data.

    fx = np.fft.fftfreq(nx, d=dx)[1:]   # frequencies (1/m)
    fy = np.fft.fftfreq(ny, d=dy)[1:]
    kx = fx*2*np.pi                     # wavenumbers (radians/m)
    ky = fy*2*np.pi
    kh = np.sqrt(kx**2 + ky**2)         # total wavenumber
    dkx = 2*np.pi/L_x                   # delta wavenumbers
    dky = 2*np.pi/L_y
    # discretize the 2D wavenumber in multiples of the maximum one-dimensional wavenumber
    dkh = np.max((dkx,dky))
    Nmax = np.int(np.ceil(np.sqrt(2)*np.max((nx/2,ny/2))))  # maximum number of points
    kp = (np.arange(Nmax-1)+1)*dkh                          # for eqn 18
    dkmin = np.min((dkx,dky))
    # Reliable Nyquist frequency (highest reliable wavenumber, shortest reliable wavelength)
    NYQwavel = 2*np.max((dy,dx))
    NYQwaven = 2*np.pi/NYQwavel
    kp_keep = kp <= NYQwaven   # wavenumbers to retain if using "restrict" option

    fkx = np.fft.fftfreq(nx, d=dx)*2*np.pi
    fky = np.fft.fftfreq(ny, d=dy)*2*np.pi
    norm = (dx*dy*dkmin)/(8*(np.pi**2)*nx*ny) # normalization factor (see eqn 24)

    # Compute the 2D fft for the full [t,y,x,z] dataset, over y and x, at once.
    # Keep as xarray to expoit dask. Means losing 'keepdims' option.
    field = (var - var.mean(dim=(xname,yname))).data
    field = field.rechunk(chunks={xx:nx, yx:ny})
    temp2 = np.fft.fft2(field, axes=(xx,yx))  # (time,...horiz...,vertical)
    Ek = (temp2*np.conj(temp2)).real

    # Populate the labels, counts, and means for each kp [for radial summation]
    if kmap is None:    # basically, this is "if this is the first time working on data with the same horizontal dimensions"
        kmap = prepare_map(fkx, fky, kp, dkh, Nmax)
        #ALWAYS ASSUMES UNIFORMITY IN FILE    # When enabled, only compute the map once.

    rlab   = kmap['rlab']
    rindex = kmap['rindex']
    kcount = kmap['kcount']
    kpbar  = kmap['kpbar']


    # Calculate the fft2 with requested computation method supplied via options
    method = options['spec_method']  # [Durran | ndimage]
    compensation = options['spec_compensation'] # [True | False]
    restrict = options['spec_restrict'] # [True | False]

    if method.upper() == 'DURRAN':
        # Prepare result array
        kpo=kp
        if restrict:   # Keep only values below Nyquist frequency
            kpo=kp[kp_keep]
            rindexo=rindex[kp_keep]
            kcounto=kcount[kp_keep]
            kpbaro=kpbar[kp_keep]
        else:
            kpo=kp
            rindexo=rindex
            kcounto=kcount
            kpbaro=kpbar


        if compensation:
            comp = (2*np.pi*kpbaro[:]/(dkmin*kcounto[:]))
            # @da.as_gufunc(signature="(i,j,k,l),(j,k),(m),(),(m),(m)->(i,m,l)",
#               output_dtypes=float, vectorize=True)
            rave = da.gufunc(rad_ave_with_comp,
                             signature="(j,k),(j,k),(m),(),(m)->(m)",
                             output_dtypes=float, vectorize=True,
                             axes=[(1,2), (0,1), (0), (), (0), (1)])
            Ekp = rave(Ek, rlab, rindexo, norm, comp)
        else:
            rave = da.gufunc(rad_ave_without_comp,
                             signature="(j,k),(j,k),(m),()->(m)",
                             output_dtypes=float, vectorize=True,
                             axes=[(1,2), (0,1), (0), (), (1)])
            Ekp = rave(Ek, rlab, rindexo, norm)

        # kplen=kpo.size
        # Ekp=np.zeros((nt,kplen,nz))
        # for tnc in np.arange(nt):
        #     for znc in np.arange(nz):
        #         # Sum points
        #         Ekp[tnc,:,znc] = norm * ndimage.sum(Ek[tnc,...,znc],
        #                                             rlab, index=rindexo)
        #         # Apply compensation  (See eqn 29)
    elif (method.upper() == 'NDIMAGE') and (nx == ny):
        # Prepare result array
        kplen=Ek[0,...,0].shape[0]//2
        Ekp=np.zeros((nt,kplen,nz))
        kpo=fkx[0:kplen]

        gpsd = da.gufunc(GetPSD1D,
                         signature="(j,k),(l)->(l)",
                         output_dtypes=float, vectorize=True,
                         axes=[(1,2), (0), (1)])

        Ekp = norm * gpsd(Ek, kpo)
        # for tnc in np.arange(nt):
        #     for znc in np.arange(nz):
        #         # Sum points
        #         Ekp[tnc,:,znc] = norm * GetPSD1D(np.fft.fftshift(Ek[tnc,...,znc]))
    else:
        print("Must supply 2D FFT spec_method: [DURRAN | NDIMAGE] and must have nx == ny for NDIMAGE")
        print("FAIL - spec_method")
        sys.exit()

    # Useful plotting values
    #kpo is radian wavenumber
    freq = kpo/(2*np.pi)
    wavel = 1/freq

    # Assemble the 2D Dataset
    ds2 = xr.Dataset(data_vars = {"spec_2d_"+vname:(["time","hfreq",zname],Ekp),
                                  "hwaven":(['hfreq'],kpo),
                                  "hwavel":(['hfreq'],wavel)  },
                     coords    = {"time": (["time"],times),
                                  zname: ([zname],z),
                                  "hfreq": (["hfreq"],freq) } )
    ds2['spec_2d_'+vname].attrs = {"units": "m * base_unit^2",
                                   "long_name": "horizontal 1D spectral density from 2D-FFT"}
    ds2['hfreq'].attrs = {"units": "m-1",
                         "long_name": "horizontal 1D spectral frequency from 2D-FFT"}
    ds2['time'].attrs = {"units": "s",
                         "long_name": "model time since start"}
    ds2[zname].attrs = {"units": "m",
                        "long_name": "height above surface"}
    ds2['hwaven'].attrs = {"units": "rad m-1",
                           "long_name": "horizontal 1D spectral wavenumber from 2D-FFT"}
    ds2['hwavel'].attrs = {"units": "m",
                           "long_name": "horizontal 1D spectral wavelength from 2D-FFT"}

# end if spec_2D
#_____________________________________________________________________________


    # Append new vname spectra to the output Dataset
    print(f"Saving {ds2} to {outfile}")

    d = ds2.to_netcdf(path = outfile, unlimited_dims="time", mode='a',
                      compute=False)
    if executing_on_cluster:
        results = d.compute()
    else:
        with ProgressBar():
            results = d.compute()
    ds2.close()

    return kmap

