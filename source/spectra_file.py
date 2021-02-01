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

    User must supply:
       dir:    input directory (slash-agnostic)
       file:   input file
                 Suggest switching to argument input (see below)
       outtag: output file tag (appended to input file name)
                 Creates 'spectra/' directory within the given dir
       dx:     x-direction grid spacing [m]
       dy:     y-direction grid spacing [m]

    @author: Todd Jones
"""

import os
import sys
import numpy as np
import xarray as xr
from scipy import ndimage      # Required for the radial summations

# Ignore divide-by-zero warnings that occur when converting between wavenumbers/frequency/wavelength
np.seterr(divide='ignore')



# ---------------------------------- USER SETTINGS ----------------------------------
# Set computation options
options = {
           'spec_1D': False,            # Compute 2D spectra as average of 1D spectra in both 
                                       #   directions
           'spec_2D': True,            # Compute 2D spectra via 2D-fft
           'spec_method': 'durran',    # [Durran, ndimage] Use Durran method (which actually also 
                                       #   uses ndimage), or faster, less accurate ndimage method
           'spec_compensation': True,  # With spec_method: 'durran', use Durran/Tanguay method to 
                                       #   compensate for systematic noise in the annular summation
                                       #   (spectra does not preserve energy when enabled)
           'spec_restrict': True       # With spec_method: 'durran', restrict the spec_2d result 
                                       #   to values below the Nyquist frequency
             }

# Single file source location and name
# dir = '/gws/nopw/j04/paracon_rdg/users/toddj/BOMEX/CA_S_SC_BOMEX_25_600/f3/test_op_RFFT/'
# dir = '/work/scratch-pw/toddj/'  # testing
# file = 'BOMEX_m0025_g0600_all_66600.0_test_plot.nc'
# -- or --
#   Modify script to receive file and/or path name
if len(sys.argv) > 2:
    dir  = sys.argv[1]
    file = sys.argv[2]
infile = os.path.join(dir, file)

# Set up outfile
outdir = os.path.join(dir, 'spectra/')
os.makedirs(outdir, exist_ok = True)  # make outdir if it doesn't already exist
outtag = "spectra_w_2D"
outfile = os.path.join(outdir,('.').join(os.path.basename(file).split('.')[:-1]) + "_"+outtag+".nc")

#dx=25.0
#dy=25.0
# Obtain dx=dy from file name '_m' encoding
mnc = infile.find('_m')
dx = float(infile[mnc+2:mnc+6])
dy = float(infile[mnc+2:mnc+6])

time_dim_always_contains='time'

var_names_spec = ['w']  # list of variable names to evaluate
        #   - leave empty to work on all present variables
        #   - populate with string variable names to work on specified list
# ---------------------------------- USER SETTINGS ----------------------------------








def get_string_index(strings, substr, exact=False):
    """
    Searches for the first index in list of strings that contains substr string
      anywhere or, if exact=True, contains an exact match

    Args:
        strings  : List of 
        substr   : Name of variable
        exact    : (optional) Whether to perform search for exact match

    Returns:
        integer  : first index of match, or None, if not found

    @author: Todd Jones
    """

    idx = None

    for idx, string in enumerate(strings):
        if exact:
            if substr == string:
                break
        else:
            if substr in string:
                break
        # end if (exact)
    return idx


#===================================================================
# Get PSD 1D (total radial power spectrum) 
# for spec_method: ndimage
#===================================================================
def GetPSD1D(psd2D):
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
    psd1D = ndimage.sum(psd2D, r, index=np.arange(0, wc))

    return psd1D
#===================================================================



def forward(n):
    return 2*np.pi/n # l = 2*np.pi/n

def inverse(l):
    return 2*np.pi/l # n = 2*np.pi/l





#################################################################################################
#################################################################################################
#################################################################################################
def main():
    '''
    Top level code
    '''

    # Open input file as xarray dataset
    ds = xr.open_dataset(infile)
    atts = ds.attrs
   
    # By default, we are going to look at all data variables from the file
    #   Save data variable names to list
    # ALTERNATIVELY, THE CODE COULD BE CONFIGURED TO PASS IN A LIST OF 
        # SPECIFIC FIELDS TO EVALUATE.
    if len(var_names_spec) == 0: 
        var_names = list(ds.data_vars.keys())
    else:
        var_names = var_names_spec.copy()

    # Initialise xarray dataset for output and append attributes with present options
    dso = xr.Dataset(data_vars = {"variable_names":(["nvars"],var_names)}, 
                     coords    = {"nvars": (["nvars"],np.arange(len(var_names)))} )
    print("Writing Dataset to ",outfile)
    # Modify the True/False options for writing.
    for inc in options:
        if options[inc] == True:
            options[inc]=1
        if options[inc] == False:
            options[inc]=0
    dso.attrs = {**atts, **options}

    #prep empty file
    dso.to_netcdf(path = outfile, mode='w')
    dso.close()

    #prep empty potential datasets
    ds1 = xr.Dataset(None)
    ds2 = xr.Dataset(None)
        # Read back names as names = dso['variable_names'].values.astype(str)

    print("Working on file: "+file)

    # Prepare map flag
    prepare_map = True
        
    # Loop over var_names
    for vname in var_names:
        
        print("Working on spectra for: ", vname)
        # Load variable
        var = ds[vname].load()

        if np.size(var.shape) != 4:
            print("Expecting 4-dimensional data for ", vname)
            print("FAIL - shape")
            sys.exit()
            # FOR OTHER KINDS OF FILES WITH MANY DIFFERENTLY STRUCTURED
            #   FIELDS, YOU MIGHT USE continue TO SKIP THOSE VARIABLES
            #   RATHER THAN sys.exit()

        # Store some parameters
        # Coordinate positions
        tx = get_string_index(var.dims,time_dim_always_contains)
        xx = get_string_index(var.dims,'x',exact = True)
        yx = get_string_index(var.dims,'y',exact = True)
        zx = get_string_index(var.dims,'z')
        
        # Time and height names
        tname = var.dims[tx]
        var   = var.rename({tname: 'time'})
        tname = "time"
        times = var[tname].values
        zname = var.dims[zx]
        z     = ds[zname].values
        
        # Coordinate sizes
        ny = var.sizes['y']
        nx = var.sizes['x']
        nt = var.sizes[tname]
        nz = var.sizes[zname]

        # Horizontal domain lengths
        L_x = dx*nx
        L_y = dy*ny
        
        
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
        if options['spec_1D']:
            print("Calculating 1D spectra")
            
            # Spectra in y-direction (mean removed)
            temp = np.fft.rfft((var.values - var.values.mean(axis = yx, keepdims = True, dtype = 'float64')), axis = yx)
            # Average across all x
            ydir = np.mean((temp*np.conj(temp)).real, axis = xx) * (dy/(2*np.pi*ny))
            # Halve the Nyqist value
            # Now assuming that data is structured ["time","yfreq",zname]
            ydir[:,-1,:] /= 2
            
            # Repeat for x-direction
            temp = np.fft.rfft((var.values - var.values.mean(axis = xx, keepdims = True, dtype = 'float64')), axis = xx)
            xdir = np.mean((temp*np.conj(temp)).real, axis = yx) * (dx/(2*np.pi*nx))
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
            ds1['yfreq'].attrs = {"units": "rad m-1",
                                 "long_name": "spectral frequency, y-direction"}
            ds1['xfreq'].attrs = {"units": "rad m-1",
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

        # end if spec_1D
        #_____________________________________________________________________________


    
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
        if options['spec_2D']:
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
            temp = np.fft.fft2((var.values - var.values.mean(axis=(yx,xx), keepdims=True, dtype='float64')),axes=(yx,xx))  # (time,...horiz...,vertical)
            Ek = (temp*np.conj(temp)).real
            
            # Populate the labels, counts, and means for each kp [for radial summation]
            if prepare_map:    # basically, this is "if this is the first time working on data with the same horizontal dimensions"
                print("Preparing map")
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
            prepare_map = False     #ALWAYS ASSUMES UNIFORMITY IN FILE    # When enabled, only compute the map once.

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

                kplen=kpo.size
                Ekp=np.zeros((nt,kplen,nz))
                for tnc in np.arange(nt):
                    for znc in np.arange(nz):
                        # Sum points
                        Ekp[tnc,:,znc] = norm * ndimage.sum(Ek[tnc,...,znc], rlab, index=rindexo)
                        # Apply compensation  (See eqn 29)
                        if compensation:
                            Ekp[tnc,:,znc] *= (2*np.pi*kpbaro[:]/(dkmin*kcounto[:]))
                
            elif (method.upper() == 'NDIMAGE') and (nx == ny):
                # Prepare result array
                kplen=Ek[tnc,...,znc].shape[0]//2
                Ekp=np.zeros((nt,kplen,nz))
                kpo=fkx[0:kplen]
                for tnc in np.arange(nt):
                    for znc in np.arange(nz):
                        # Sum points
                        Ekp[tnc,:,znc] = norm * GetPSD1D(np.fft.fftshift(Ek[tnc,...,znc]))
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
            ds2['hfreq'].attrs = {"units": "rad m-1",
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
        dso = xr.merge([ds1, ds2])
        dso.to_netcdf(path = outfile, unlimited_dims="time", mode='a')
        ds1.close()
        ds2.close()
        dso.close()

    #endfor vnc loop over var_names

    ds.close()
    print("IAMDONE")  # Tag for successful completion.  
    
# END main


if __name__ == "__main__":
    main()
