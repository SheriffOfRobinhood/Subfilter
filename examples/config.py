# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 13:20:58 2021

@author: paclk
"""
subfilter_setup = {'write_sleeptime':3,
                   'use_concat':True,
                   'chunk_size':2**22 }

dask_opts={
            'no_dask': False,
            'use_map_overlap': True
            }
