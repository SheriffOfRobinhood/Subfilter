# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 11:33:51 2021

@author: paclk
"""
import numpy as np
from subfilter.utils.string_utils import get_string_index

def re_chunk(f, chunks = None, xch = 'all', ych = 'all', zch = 'auto'):

    defn = 1

    if chunks is None:

        chunks={}
        sh = np.shape(f)
        for ip, dim in enumerate(f.dims):
            if 'x' in dim:
                if xch == 'all':
                    chunks[dim] = sh[ip]
                else:
                    chunks[dim] = np.min([xch, sh[ip]])
            elif 'y' in dim:
                if ych == 'all':
                    chunks[dim] = sh[ip]
                else:
                    chunks[dim] = np.min([ych, sh[ip]])
            elif 'z' in dim:
                if zch == 'all':
                    chunks[dim] = sh[ip]
                elif zch == 'auto':
                    chunks[dim] = 'auto'
                else:
                    chunks[dim] = np.min([zch, sh[ip]])
            else:
                chunks[f.dims[ip]] = defn

    f = f.chunk(chunks=chunks)

    return f

def guess_chunk(max_ch, dataset):

    [itime, iix, iiy, iiz] = get_string_index(dataset.dims,
                                              ['time', 'x', 'y', 'z'])
    timevar = list(dataset.dims)[itime]
    xvar = list(dataset.dims)[iix]
    yvar = list(dataset.dims)[iiy]
    zvar = list(dataset.dims)[iiz]

    nch = np.min([int(dataset.dims[xvar]/(2**int(np.log(dataset.dims[xvar]
                                                *dataset.dims[yvar]
                                                *dataset.dims[zvar]
                                                /max_ch)/np.log(2)/2))),
                  dataset.dims[xvar]])
    return nch