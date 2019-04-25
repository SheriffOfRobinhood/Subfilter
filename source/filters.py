import numpy as np
from sys import float_info

# Global constants
#===============================================================================
pi = np.pi # mathematical constant
n = 101 # resolution for calculation of mean of non-linear filters.
eps = float_info.min # smallest possible float

#===============================================================================
class filter_2d :
    def __init__(self, filter_id, filter_name, wavenumber=-1, \
                         delta_x=1.0, width=-1,cutoff=0.0001, \
                         high_pass=0, sigma=-1):
        '''    Args:
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
        '''    
            
        if (filter_name == 'domain'):
            data = np.ones([1,1])
         
        elif (filter_name == 'gaussian'):
            if (sigma == -1):
                data = filter_2d_error(filter_name, 'sigma')
            else:
                data = gaussian_filter_2d(sigma, delta_x, cutoff, width)
        elif (filter_name == 'running_mean'):
            if (width == -1):
                data = filter_2d_error(filter_name, 'width')
            else:
                data = running_mean_filter_2d(width)
        elif (filter_name == 'wave_cutoff'):
            if (wavenumber == -1):
                data = filter_2d_error(filter_name, 'wavenumber')
            else:
                data = wave_cutoff_filter_2d(wavenumber, delta_x, width,
                                             cutoff, high_pass)     
        else:
            print('This filter type is not available.')
            print('Available filters are:')
            print('domain, gaussian, running_mean & wave_cutoff')
            data = -9999
            
        if (np.size(np.shape(data)) > 1 ) : 
            self.data = data
            
#            x = filter_dataset.createDimension('x'+filter_id[6:],np.shape(data)[0])
#            y = filter_dataset.createDimension('y'+filter_id[6:],np.shape(data)[1])
#            filter_data = filter_dataset.createVariable(filter_id, "f8",\
#                                    ("x"+filter_id[6:],"y"+filter_id[6:]))
#            filter_data[:,:] = data
#            filter_data.filter_name = filter_name
#            filter_data.sigma = sigma
#            filter_data.wavenumber = wavenumber
#            filter_data.delta_x = delta_x
#            filter_data.width = width
#            filter_data.cutoff = cutoff
#            filter_data.high_pass = high_pass  
            self.id = filter_id
            self.attributes = {'filter_type' : filter_name, \
                  'wavenumber' : wavenumber, \
                  'delta_x' : delta_x, \
                  'width' : width, \
                  'cutoff' : cutoff, \
                  'high_pass' : high_pass, \
                  'sigma' : sigma}
            
    def __str__(self):
        rep = "Filter (2D) id: {0}\n".format(self.id)
#        rep += self.attributes.__str__()
        for attr in self.attributes:
            rep += "{0}: {1}\n".format(attr, self.attributes[attr])
        return rep
           
    def __repr__(self):
        rep = "filter_2d:"
        rep += " id: {0}, data{1}, attributes{2}\n".format(self.id,\
                     np.shape(self.data), \
                     self.attributes)
        return rep
    
    def filter_2d_error(filter_name, problem):
        '''
        Prints error when parameter required by filter does not exist.
    
        Args:
          filter_name (str): Name of filter
          problem (str): Name of parameter that has not been set
    
        Returns:
          filter_2d (-9999): Error code for filter.
        '''
        print('A ' + filter_name + ' filter was selcted, but a suitable value')
        print('for the ' + problem + ' was not chosen')
        filter_2d = -9999
        return filter_2d
    
def running_mean_filter_2d(width):
    '''
    Calculates a square 2D running mean filter with the given width

    Args:
      width (int): width of the filter

    Returns:
      ndarray: 2D array of size width*width. Every element equals 1.0/(width**2)
    '''
    width = int(width)
    return (np.ones(width**2)/width**2).reshape(width, width)


def running_mean_filter(width):
    '''
    Calculates a 1D running mean filter with the given width

    Args:
      width (int): width of the filter

    Returns:
      ndarray: 1D array of size width. Every element equals 1.0/(width)
    '''
    width = int(width)
    return np.ones(width)/width


def wave_cutoff_filter(wavenumber, delta_x=1.0, width=-1, cutoff=0.0001,
                       high_pass=0):
    '''
    Calculates a 1D wave-cutoff filter caculated using the given wavenumber.
    
    Uses filter(x) = sin(wavenumber*x)/x. Normalised by sum(filter(x)).
    Note that this returns the average value of filter(x) over the range
    (x-delta_x/2.0, x+delta_x/2.0).
    
    Args:
      wavenumber (float):
      delta_x (float, default=1.0): The distance between two points in the data
                                    that the filter will be applied to.
      width (int, default=-1): If not -1, used to explicitly set the width of the
                               filter.
      cutoff (float, default=0.0001): If width=-1, the width of the filter is set
                                      dynamically, and increased until the
                                      smallest value of the filter is less than
                                      the cutoff value.
      high_pass (bool, default=0): If true a high pass filter is calculated

    Returns:
      ndarray: 1D array of filter values
    '''
    if high_pass == 0:
        if width == -1:
            width = 1
            result = np.ones(1)
            while abs(result).min()/pi > cutoff:
                width  += 2
                if width/2 == width/2.0:
                    x = np.concatenate((np.arange(width/2)[::-1],
                                        np.arange(width/2)))
                else:  
                    x = np.arange(width)-(width/2)
                x = x * delta_x
                result = np.zeros(width)
                for i in range(width):
                    temp = delta_x*np.arange(n)/(n-1)+x[i]-(delta_x/2.0)
#                    temp = np.zeros(1) + x[i] # If don't want to do averaging
                    temp[temp == 0] = eps
#                    result[i] = np.mean((np.sin(wavenumber*temp)/temp)/pi)
                    result[i] = np.mean(np.sin(wavenumber*temp)/temp,dtype=np.float64)
        else:
            if width/2 == width/2.0:
                x = np.concatenate((np.arange(width/2)[::-1],
                                    np.arange(width/2)))
            else:  
                x = np.arange(width)-(width/2)
            x = x * delta_x
            result = np.zeros(width)
            for i in range(width):
                temp = delta_x*np.arange(n)/(n-1)+x[i]-(delta_x/2.0)
#                temp = np.zeros(1) + x[i] # If don't want to do averaging
                temp[temp == 0] = eps
#                result[i] = np.mean((np.sin(wavenumber*temp)/temp)/pi)
                result[i] = np.mean(np.sin(wavenumber*temp)/temp,dtype=np.float64)
        result = result/np.sum(result)
        return result
    else:
        print(' high pass filter not yet coded!')
        return 1.0


def wave_cutoff_filter_2d(wavenumber, delta_x=1.0, width=-1, cutoff=0.0001,
                          high_pass=0):
    '''
    Calculates a 2D wave-cutoff filter caculated using the given wavenumber.
    
    Uses filter(x,y) = sin(wavenumber*x)/x * sin(wavenumber*y)/y. 
    Normalised by sum(filter(x,y)).
    Note that this returns the average value of filter(x) over the range
    (x-delta_x/2.0, x+delta_x/2.0).
    
    Args:
      wavenumber (float):
      delta_x (float, default=1.0): The distance between two points in the data
                                    that the filter will be applied to.
      width (int, default=-1): If not -1, used to explicitly set the width of the
                               filter.
      cutoff (float, default=0.0001): If width=-1, the width of the filter is set
                                      dynamically, and increased until the
                                      smallest value of the filter is less than
                                      the cutoff value.
      high_pass (bool, default=0): If true a high pass filter is calculated

    Returns:
      ndarray: 2D array of filter values
    '''
# Uses sin(Ax)*sin(Ay)/(x*y). Need to check this equation is correct.
    if high_pass == 0:
        if width == -1:
            width = 1
            result = np.ones(1).reshape(1,1)
            cutoff = cutoff * pi * pi # multiply cutoff by normalising factor rather than dividing filter each time
            while abs(result[:(width+1)/2,:(width+1)/2]).min() > cutoff:
                width  += 2
                if width/2 == width/2.0:
                    x = np.concatenate((np.arange(width/2)[::-1],
                                        np.arange(width/2)))
                else:  
                    x = np.arange(width)-(width/2)
                x = x * delta_x
                result = np.zeros(width**2).reshape(width,width)
                for i in range((width+1)/2):
                    for j in range((width+1)/2):
                        temp1 = np.tile((delta_x*np.arange(n)/(n-1)+x[i]
                                         -(delta_x/2.0)), (n,1))
                        temp2 = np.tile((delta_x*np.arange(n)/(n-1)+x[j]
                                         -(delta_x/2.0)), (n,1)).transpose()
#                        temp1 = np.zeros(1) + x[i] #Use if don't want to average (may be necessary to get sensible spectral filtering?)
#                        temp2 = np.zeros(1) + x[j]
                        temp1[temp1 == 0] = eps
                        temp2[temp2 == 0] = eps
                        result[i, j] = np.mean((np.sin(wavenumber*temp1)/temp1)
                                               *(np.sin(wavenumber*temp2)/temp2),dtype=np.float64)
        else:
            if width/2 == width/2.0:
                x = np.concatenate((np.arange(width/2)[::-1], np.arange(width/2)))
            else:  
                x = np.arange(width)-(width/2)
            x = x * delta_x
            result = np.zeros(width**2).reshape(width,width)
            for i in range((width+1)/2):
                for j in range((width+1)/2):
                    temp1 = np.tile((delta_x*np.arange(n)/(n-1)+x[i]
                                     -(delta_x/2.0)), (n,1))
                    temp2 = np.tile((delta_x*np.arange(n)/(n-1)+x[j]
                                     -(delta_x/2.0)), (n,1)).transpose()
#                    temp1 = np.zeros(1) + x[i]
#                    temp2 = np.zeros(1) + x[j]
                    temp1[temp1 == 0] = eps
                    temp2[temp2 == 0] = eps
                    result[i, j] = np.mean((np.sin(wavenumber*temp1)/temp1)
                                           *(np.sin(wavenumber*temp2)/temp2),dtype=np.float64)
        temp = result[0:(width+1)/2,0:(width+1)/2]
        temp = np.concatenate((temp[:,:width/2], temp[:,::-1]),axis=1)
        result = np.concatenate((temp[:width/2,:], temp[::-1,:]),axis=0)
        result = result/np.sum(result)
        return result
    else:
        print(' high pass filter not yet coded!')
        return 1.0


def gaussian_filter(sigma, delta_x=1.0, cutoff = 0.0001, width = -1):
    '''
    Calculates a 1D Gaussian filter caculated with the given lengthscale (sigma)
    
    Uses filter(x) = EXP(-x^2/2.0*sigma^2). Normalised by sum(filter(x)).
    Note that this returns the average value of filter(x) over the range
    (x-delta_x/2.0, x+delta_x/2.0).
    
    Args:
      sigma (float): The lengthscale of the filter. 
      delta_x (float, default=1.0): The distance between two points in the data
                                    that the filter will be applied to.
      width (int, default=-1): If not -1, used to explicitly set the width of the
                               filter.
      cutoff (float, default=0.0001): If width=-1, the width of the filter is set
                                      dynamically, and increased until the
                                      smallest value of the filter is less than
                                      the cutoff value.

    Returns:
      ndarray: 1D array of filter values
    '''
    if width == -1:
        width = 1
        result = np.ones(1)
        cutoff = cutoff * (np.sqrt(2.0 * pi) * sigma)
        while result.min() > cutoff:
            width += 2
            x = np.arange((width+1)/2) * delta_x
            result = np.zeros((width+1)/2)
            for i in range((width+1)/2):
                y = (delta_x*np.arange(n)/(n-1))+x[i]-(delta_x/2.0)
                result[i] = np.mean(np.exp(-y**2 / (2.0 * sigma**2)),dtype=np.float64)
            if width/2 == width/2.0:
                result = np.concatenate((result[::-1], result))
            else:
                result = np.concatenate((result[-1:0:-1], result))
    else:
        x = np.arange((width+1)/2) * delta_x
        result = np.zeros((width+1)/2)
        for i in range((width+1)/2):
            y = (delta_x*np.arange(n)/(n-1))+x[i]-(delta_x/2.0)
            result[i] = np.mean(np.exp(-y**2 / (2.0 * sigma**2)),dtype=np.float64)
        if width/2 == width/2.0:
            result = np.concatenate((result[::-1], result))
        else:
            result = np.concatenate((result[-1:0:-1], result))
    result = result / np.sum(result) # need to scale result to ensure doesn't affect mean in case cutoff is large or width is too small. Shouldn't have much effect if this is not the case.
    return result


def gaussian_filter_2d(sigma, delta_x=1.0, cutoff=0.0001, width=-1):
    '''
    Calculates a 2D Gaussian filter caculated with the given lengthscale (sigma)
    
    Uses filter(x,y) = EXP(-(x+y)^2/2.0*sigma^2). Normalised by sum(filter(x)).
    Note that this returns the average value of filter(x) over the range
    (x-delta_x/2.0, x+delta_x/2.0).
    
    Args:
      sigma (float): The lengthscale of the filter. 
      delta_x (float, default=1.0): The distance between two points in the data
                                    that the filter will be applied to.
      width (int, default=-1): If not -1, used to explicitly set the width of the
                               filter.
      cutoff (float, default=0.0001): If width=-1, the width of the filter is set
                                      dynamically, and increased until the
                                      smallest value of the filter is less than
                                      the cutoff value.

    Returns:
      ndarray: 2D array of filter values
    '''
    if width == -1:
        width = 1
        result = np.ones(1).reshape(1,1)
        cutoff = cutoff * (np.sqrt(2.0 * pi) * sigma)**2
        while result[:(width+1)//2,:(width+1)//2].min() > cutoff:
            width += 2
            if width/2 == width/2.0:
                x = np.concatenate((-np.arange(width/2)[::-1],
                                    np.arange(width/2)))
            else:  
                x = np.arange(width)-(width/2)
            x = x * delta_x
            result = np.zeros(width**2).reshape(width,width)
            for i in range((width+1)//2):
                for j in range((width+1)//2):
                    temp1 = np.tile((delta_x*np.arange(n)/(n-1)+x[i]
                                     -(delta_x/2.0))**2, (n,1))
                    temp2 = np.tile((delta_x*np.arange(n)/(n-1)+x[j]
                                     -(delta_x/2.0))**2, (n,1)).transpose()
                    y = temp1 + temp2
                    result[i,j] = np.mean(np.exp(-y / (2.0 * sigma**2)),dtype=np.float64)
    else:
        if width/2 == width/2.0:
            x = np.concatenate((-np.arange(width/2)[::-1], np.arange(width/2)))
        else:  
            x = np.arange(width)-(width/2)
        x = x * delta_x
        result = np.zeros(width**2).reshape(width,width)
        for i in range((width+1)//2):
            for j in range((width+1)//2):
                temp1 = np.tile((delta_x*np.arange(n)/(n-1)+x[i]
                                 -(delta_x/2.0))**2, (n,1))
                temp2 = np.tile((delta_x*np.arange(n)/(n-1)+x[j]
                                 -(delta_x/2.0))**2, (n,1)).transpose()
                y = temp1 + temp2
                result[i,j] = np.mean(np.exp(-y / (2.0 * sigma**2)),dtype=np.float64)
    temp = result[0:(width+1)//2,0:(width+1)//2]
    temp = np.concatenate((temp[:,:width//2], temp[:,::-1]),axis=1)
    result = np.concatenate((temp[:width//2,:], temp[::-1,:]),axis=0)
    result = result / np.sum(result)
    return result


