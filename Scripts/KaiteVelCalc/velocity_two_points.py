# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 15:10:35 2021

@author: USER_2
"""
import numpy as np

###############################################################################
# Calculate velocity between consecutive points (careful of zeros)
# Forward difference scheme, first order (f(x0+ih))
# KCole PEPT Flotation 020221
###############################################################################

def velocity_two_points(location_data):    
    x = location_data['x']
    y = location_data['y']
    z = location_data['z']
    t = location_data['t']
    dxdt = np.zeros((len(t)))
    dydt = np.zeros((len(t)))
    dzdt = np.zeros((len(t)))
    
    for l in range(0,len(t)-1):    
        deltax = x[l+1]-x[l]
        deltay = y[l+1]-y[l]
        deltaz = z[l+1]-z[l]
        deltat = t[l+1]-t[l]
    
        if deltat == 0 or deltax == 0 or deltay == 0 or deltaz == 0:
            if (t[l+1]-t[l-1] != 0) and (x[l+1]-x[l-1] != 0) and (y[l+1]-y[l-1] != 0) and (z[l+1]-z[l-1] != 0):
                dxdt[l] = (x[l+1]-x[l-1])/(t[l+1]-t[l-1])
                dydt[l] = (y[l+1]-y[l-1])/(t[l+1]-t[l-1])
                dzdt[l] = (z[l+1]-z[l-1])/(t[l+1]-t[l-1])
            else:
                dxdt[l] = 0
                dydt[l] = 0
                dzdt[l] = 0          
        else:
            dxdt[l] = (x[l+1]-x[l])/(t[l+1]-t[l])
            dydt[l] = (y[l+1]-y[l])/(t[l+1]-t[l])
            dzdt[l] = (z[l+1]-z[l])/(t[l+1]-t[l])
    sp = np.sqrt(dxdt**2 + dydt**2 + dzdt**2)
    velocity_data = {'dxdt':dxdt, 'dydt':dydt, 'dzdt':dzdt, 'sp':sp, 't':t}
    return velocity_data