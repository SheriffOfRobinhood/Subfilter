"""

  filters.py

    - This module contains the code to generate a selection of 2-dimensional filters.

"""

import numpy as np
from sys import float_info

# Global constants
#===============================================================================
pi = np.pi # mathematical constant
n = 101 # resolution for calculation of mean of non-linear filters.
eps = float_info.min # smallest possible float

#===============================================================================
class Filter :
    '''
    Class defining a filter function.

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
        ndim (int): Number of dimensions (default=2)

    @author: Peter Clark

    '''

    def __init__(self, filter_id, filter_name,
                 delta_x=1000.0, cutoff=0.000001, npoints = None,
                 high_pass=0, wavenumber=-1, width=-1, sigma=-1,
                 ndim=2):

        if (filter_name == 'domain'):
            data = np.ones([1,1])

        elif (filter_name == 'gaussian'):
            if (sigma == -1):
                data = self.filter_error(filter_name, 'sigma')
            else:
                data = gaussian_filter(sigma, delta_x, cutoff, npoints,
                                               ndim=ndim)
        elif (filter_name == 'running_mean'):
            if (width == -1):
                data = self.filter_error(filter_name, 'width')
            else:
                data = running_mean_filter(width, npoints, ndim=ndim)
                width = np.shape(data)[0]
        elif (filter_name == 'wave_cutoff'):
            if (wavenumber == -1):
                data = self.filter_error(filter_name, 'wavenumber')
            else:
                data = wave_cutoff_filter(wavenumber, delta_x, npoints,
                                          cutoff, high_pass,
                                          ndim=ndim)
        elif (filter_name == 'circular_wave_cutoff'):
            if (wavenumber == -1):
                data = self.filter_error(filter_name, 'wavenumber')
            else:
                data = circular_wave_cutoff_filter(wavenumber, delta_x,
                                                   npoints,
                                                   cutoff, high_pass,
                                                   ndim=ndim)
        else:
            print('This filter type is not available.')
            print('Available filters are:')
            print('domain, gaussian, running_mean, wave_cutoff & '
                  'circular_wave_cutoff')
            data = -9999

        if (np.size(np.shape(data)) > 1 ) :
            self.data = data

            self.id = filter_id
            self.attributes = {'filter_type' : filter_name,
                  'ndim' : ndim,
                  'wavenumber' : wavenumber,
                  'delta_x' : delta_x,
                  'width' : width,
                  'cutoff' : cutoff,
                  'high_pass' : high_pass,
                  'sigma' : sigma}

    def __str__(self):
        rep = "Filter id: {0}\n".format(self.id)
#        rep += self.attributes.__str__()
        for attr in self.attributes:
            rep += "{0}: {1}\n".format(attr, self.attributes[attr])
        return rep

    def __repr__(self):
        rep = "filter:"
        rep += " id: {0}, data{1}, attributes{2}\n".format(self.id,\
                     np.shape(self.data), \
                     self.attributes)
        return rep

    def filter_error(filter_name, problem):
        '''
        Prints error when parameter required by filter does not exist.

        Args:
          filter_name (str): Name of filter
          problem (str): Name of parameter that has not been set

        Returns:
          filter_err (-9999): Error code for filter.
        '''
        print(f'A {filter_name:s} filter was selected, but a suitable value')
        print(f'for the {problem:s} was not chosen.')
        filter_err = -9999
        return filter_err

def running_mean_filter(width, npoints, ndim=2):
    '''
    Calculates a square 1 or 2D running mean filter with the given width

    Args:
        width (int): width of the filter
        ndim (int): Number of dimensions (default=2)

    Returns:
        ndarray: ndim dimensional array of size width in each dimension.
          Every element equals 1.0/(width**ndim)
    '''
    width = int(width)
    if npoints is None:
        npoints = width

    if ndim == 1:
        result = np.ones(width)/width
        result = np.pad(result, npoints)
    else:
        result = np.ones((width,width))/(width*width)
        result = np.pad(result,
                ((npoints-width)//2, (npoints-width-(npoints-width)//2)))
    return result


def is_npi(x, tol=0.000001):
    r = np.abs(np.pi*np.round(x/np.pi )- x) <= tol
    return r

def wave_cutoff_filter(wavenumber, delta_x=1000.0, npoints=-1, cutoff=0.000001,
                       high_pass=0, ndim=2):
    '''
    Calculates a 2D wave-cutoff filter caculated using the given wavenumber.

    Uses filter(x,y) = :math:`\sin(wavenumber * x)/x * \sin(wavenumber * y)/y`
    in 2D.
    Normalised by sum(filter(x,y)).
    Note that this returns the point sampled value of filter(x).

    Args:
      wavenumber (float):
      delta_x (float, default=1000.0): The distance between two points in the
                                    data that the filter will be applied to.
      npoints (int, default=-1): If not -1, used to explicitly set the npoints of
                               the filter.
      cutoff (float, default=0.0001): If npoints=-1, the npoints of the filter is
                                      set dynamically, and increased until the
                                      smallest value of the filter is less
                                      than the cutoff value.
      high_pass (bool, default=0): If true a high pass filter is calculated
      ndim (int): Number of dimensions (default=2)

    Returns:
      ndarray: 2D array of filter values
    '''
    if high_pass:
        print("High pass not yet coded.")
        return None
    if is_npi(wavenumber*delta_x):
        print("Use fixed npoints as wavenumber*delta_x = n * pi")
        return None
    if npoints == -1:
        if is_npi(wavenumber*delta_x):
            print("Use fixed npoints as wavenumber*delta_x = n * pi")
            return None

        half_width = 0
        if ndim == 1:
            result = np.ones((1))
            rmin = 1
            while rmin > cutoff:
                half_width += 1
                L = half_width * delta_x
                x = np.linspace(-L, L, 2*half_width+1)
                x[x == 0] = eps
                result = np.sin(wavenumber*x) / x
                result /= np.sum(result)
                rmin = np.abs(result[~is_npi(wavenumber*x)]).min()
            npoints = 2 * half_width+1
        else:
            result = np.ones((1,1))
            rmin = 1
            while rmin > cutoff:
                half_width += 1
                L = half_width * delta_x
                c = np.linspace(-L, L, 2*half_width+1)
                x, y = np.meshgrid(c, c)
                x[x == 0] = eps
                y[y == 0] = eps
                result = np.sin(wavenumber*x) / x * np.sin(wavenumber*y) / y
                result /= np.sum(result)
                rmin = np.abs(result[np.logical_and(
                                        ~is_npi(wavenumber*x),
                                        ~is_npi(wavenumber*y))]).min()
            npoints = 2 * half_width+1
    else:
        if ndim == 1:
            L = (npoints-1)/2 * delta_x
            x = np.linspace(-L, L, npoints)
            x[x == 0] = eps
            result = np.sin(wavenumber*x) / x
            result /= np.sum(result)
        else:
            L = (npoints-1)/2 * delta_x
            c = np.linspace(-L, L, npoints)
            x, y = np.meshgrid(c, c)
            x[x == 0] = eps
            y[y == 0] = eps
            result = np.sin(wavenumber*x) / x * np.sin(wavenumber*y) / y
            result /= np.sum(result)
    return result

def circular_wave_cutoff_filter(wavenumber, delta_x=1000.0, npoints=-1,
                       cutoff=0.000001, high_pass=0, ndim=2):
    '''
    Calculates a 2D wave-cutoff filter caculated using the given wavenumber.

    Uses filter(x,y) = :math:`\sin(wavenumber * x)/x * \sin(wavenumber * y)/y`
    in 2D.
    Normalised by sum(filter(x,y)).
    Note that this returns the point sampled value of filter(x).

    Args:
      wavenumber (float):
      delta_x (float, default=1000.0): The distance between two points in the data
                                    that the filter will be applied to.
      npoints (int, default=-1): If not -1, used to explicitly set the npoints of the
                               filter.
      cutoff (float, default=0.0001): If npoints=-1, the npoints of the filter is set
                                      dynamically, and increased until the
                                      smallest value of the filter is less than
                                      the cutoff value.
      high_pass (bool, default=0): If true a high pass filter is calculated
      ndim (int): Number of dimensions (default=2)

    Returns:
      ndarray: 2D array of filter values
    '''
    if high_pass:
        print("High pass not yet coded.")
        return None

    if npoints == -1:
            print("Use fixed npoints.")
            return None
    else :

        if wavenumber < (2 * np.pi) / (npoints * delta_x):
            print("Wave number too small.")
            return None

        if ndim == 1:

            k = np.fft.fftshift(np.fft.fftfreq(npoints,delta_x / (2 * np.pi)))
            filt = np.ones((npoints)) +0j
            filt[k > wavenumber] = 0.0
            filt = np.fft.ifftshift(filt)
            result = np.fft.ifftshift(np.fft.ifft(filt)).real

        else:

            frq = np.fft.fftshift(np.fft.fftfreq(npoints, delta_x /(2*np.pi)))
            kx, ky = np.meshgrid(frq, frq)
            k = np.sqrt(kx * kx + ky * ky)
            filt = np.ones((npoints, npoints)) +0j
            filt[k > wavenumber] = 0.0
            filt = np.fft.ifftshift(filt)
            result = np.fft.ifftshift(np.fft.ifft2(filt)).real

    return result


def gaussian_filter(sigma, delta_x=1000.0, cutoff=0.000001, npoints=-1,
                       ndim=2):
    '''
    Calculates a 1 or 2D Gaussian filter calculated with the given lengthscale (sigma)

    Uses filter(x,y) = :math:`\exp(-(x^2+y^2)/(2\sigma^2))` in 2D.
    Normalised by sum(filter(x)).
    Note that this returns the point sampled value of filter(x).

    Args:
      sigma (float): The lengthscale of the filter.
      delta_x (float, default=1000.0): The distance between two points in the data
                                    that the filter will be applied to.
      npoints (int, default=-1): If not -1, used to explicitly set the npoints of the
                               filter in gridpoints.
      cutoff (float, default=0.0001): If npoints=-1, the npoints of the filter is set
                                      dynamically, and increased until the
                                      smallest value of the filter is less than
                                      the cutoff value.
      ndim (int): Number of dimensions (default=2)

    Returns:
      ndarray: 2D array of filter values
    '''
    if npoints == -1:
        half_width = 0
        result = np.ones((2))
        while result.min() > cutoff:
            half_width += 1
            L = half_width * delta_x
            if ndim == 1:
                x = np.linspace(-L, L, 2*half_width+1)
                r_sq = x * x
            else:
                c = np.linspace(-L, L, 2*half_width+1)
                x, y = np.meshgrid(c, c)
                r_sq = x * x + y * y
            result = np.exp(-r_sq/(2 * (sigma**2)))
            result /= np.sum(result)
        npoints = 2 * half_width + 1
    else:
        L = (npoints-1)/2 * delta_x
        if ndim == 1:
            x = np.linspace(-L, L, npoints)
            r_sq = x * x
        else:
            c = np.linspace(-L, L, npoints)
            x, y = np.meshgrid(c, c)
            r_sq = x * x + y * y
        result = np.exp(-r_sq/(2 * (sigma**2)))
        result /= np.sum(result)

    print(f"cutoff = {cutoff}, min={np.min(result)}")

    return result
