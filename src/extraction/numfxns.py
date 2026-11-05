#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 12:20:08 2016

@author: Debarany
"""
import numpy as np
from scipy.signal import butter, filtfilt
import os

def fcount(path, map = {}):
    """Count the number of folders in a directory"""
    count = 0
    for f in os.listdir(path):
        child = os.path.join(path, f)
        if os.path.isdir(child):
            child_count = fcount(child, map)
            count += child_count + 1 # unless include self
    map[path] = count
    return count

def butter_lowpass(cutoff,fs,order):
    """Specify the parameters for a lowpass butterworth filter"""
    nyq = .5 * fs
    normal_cutoff = cutoff/nyq
    b, a = butter(order, normal_cutoff, btype = 'low', analog=False) #digital filter for normally sampled data
    return b, a

def butter_lowpass_filter(data,cutoff,fs,order):
    """Apply a lowpass butterworth filter"""
    b, a = butter_lowpass(cutoff,fs,order=order)
    y = filtfilt(b,a,data)
    return y

def find_nearest_above(my_array, target):
    """Find the index of the nearest target value in an array that is above the target value"""
    diff = my_array - target
    mask = np.ma.less_equal(diff, 0)
    # We need to mask the negative differences and zero
    # since we are looking for values above
    if np.all(mask):
        return None # returns None if target is greater than any value
    masked_diff = np.ma.masked_array(diff, mask)
    return masked_diff.argmin()
    
def find_nearest(array,value):
    """Find the index of the nearest target value in an array"""
    idx = (np.abs(array-value)).argmin()
    return idx