# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 22:31:31 2021

@author: USER_2
"""
import numpy as np
from scipy import interpolate
from scipy.interpolate import UnivariateSpline

###############################################################################
# Uses a spline to estimate instantaneous velocities at a point in Cartesian co-ordinates
# KCole PEPT Flotation 050221
###############################################################################

def velocity_spline_cart(location_data, config):
    x = location_data['x']
    y = location_data['y']
    z = location_data['z']
    t = location_data['t']
    dxdt = np.zeros((len(t)))
    dydt = np.zeros((len(t)))
    dzdt = np.zeros((len(t)))
    for l in range(4,len(t)-4):
        tsample_test = t[l-4:l+5]
        xsample_test = x[l-4:l+5]
        ysample_test = y[l-4:l+5]
        zsample_test = z[l-4:l+5]
        # Test for duplicates in the t sample and interpolate sample if necessary
        tsample_set = set(tsample_test)
        if len(tsample_set) != len(tsample_test):
            tsample = np.linspace(tsample_test[1], tsample_test[-1], 9)
            f_xsample = interpolate.interp1d(tsample_test, xsample_test)
            xsample = f_xsample(tsample)
            f_ysample = interpolate.interp1d(tsample_test, ysample_test)
            ysample = f_ysample(tsample)
            f_zsample = interpolate.interp1d(tsample_test, zsample_test)
            zsample = f_zsample(tsample)
        else: 
            tsample = tsample_test
            xsample = xsample_test
            ysample = ysample_test
            zsample = zsample_test
        # Fit spline
        porder = 4
        splx = UnivariateSpline(tsample, xsample, k=porder)
        sply = UnivariateSpline(tsample, ysample, k=porder)
        splz = UnivariateSpline(tsample, zsample, k=porder)
        rmsdevx = np.sqrt(np.sum((splx(tsample)-xsample)**2))
        rmsdevy = np.sqrt(np.sum((sply(tsample)-ysample)**2))
        rmsdevz = np.sqrt(np.sum((splz(tsample)-zsample)**2))
        
        # Find velocity
        dsplx_dt = splx.derivative()
        dsply_dt = sply.derivative()
        dsplz_dt = splz.derivative()
        dxdt[l] = dsplx_dt(tsample[5])
        dydt[l] = dsply_dt(tsample[5])
        dzdt[l] = dsplz_dt(tsample[5])
        
        if np.mod(l,10000) == 0:
            print('location ' + str(l) + '/' + str(len(t)))
    sp = np.sqrt(dxdt**2 + dydt**2 + dzdt**2)
    velocity_data_cart = {'dxdt':dxdt, 'dydt':dydt, 'dzdt':dzdt, 'sp':sp, 't':t, 'porder':porder, 'rmsdevx':rmsdevx, 'rmsdevy':rmsdevy, 'rmsdevz':rmsdevz }
    return velocity_data_cart