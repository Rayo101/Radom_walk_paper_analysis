"""
Script containing the 3 different methods for velocity comptuation.  
The 3 different methods are a 2 point, 6 point and spline velocity.  
Based on velocity_six_points.py, velocity_two_points.py and 
velocity_spline_cart.py by Katie Cole

Author: Rayhaan Perin 
Date: 26-03-2024
Packages:  numpy, pept, scipy
Optional Packages: None

Change List: 
                - 26-03-2024 - Rayhaan - Creation of Script
                - 27-03-2024 - Rayhaan - Added methods for velocity calculation and fit method
                - 27-03-2024 - Rayhaan - Handled additional argument of porder for spline method
                - 27-03-2024 - Rayhaan - Added exception handling
                - 27-03-2024 - Rayhaan - Added additional comments to class and methods

How to use:  You can place this script anywhere you'll just have to add it to path using sys or os. 
You should then be able to import the script then class and calculate velocities.

e.g. 

import sys
sys.path.append("/home/rayhaan/randomWalk_V3/Scripts/") # this is where the script is located in the file system.
from velocityCalculation import Velocity
"""
# PACKAGES
import numpy as np
from pept import PointData
from scipy import interpolate
from scipy.interpolate import UnivariateSpline

# CLASSES
class Velocity:
    """
    Class to handle all velocity calculations in cartesian coordinates.  
    Three different methods are available namely the two point method,
    six point method and spline fitting for velocity computation.
    """
    def __init__(self, locations: PointData):
        """
        Constructor for Veloicty class

        Input:
                locations: PointData - the locations of which you want to compute the velocity
        """

        if not isinstance(self.locations, PointData):
            raise TypeError("locations need to be of type PointData from the pept package!")
        self.locations = locations

    def _two_point_method(self):
        """
        Implementation of the two point method.  v = delta(X)/delta(t).  
        There is handling of whether dx, dy, dz or dt = 0.  

        Output:
                out: dictionary - Contains t, vx, vy, vz, v (speed)
        """
        
        # the positions and times
        x: np.ndarray = self.locations['x']
        y: np.ndarray = self.locations['x']
        z: np.ndarray = self.locations['z']
        t: np.ndarray = self.locations['t']

        # to store velocities (no real point in memory allocation due to run time inference)
        vx: list = []
        vy: list = []
        vz: list = []
        vt: list = []

        for i in range(len(x) - 1):
            dx = x[i+1] - x[i]
            dy = y[i+1] - y[i]
            dz = z[i+1] - z[i]
            dt = t[i+1] - t[i]
            averageT = (t[i+1] + t[i])/2

            if (dx == 0) or (dy == 0) or (dz == 0) or (dt == 0):
                if (t[i+1]-t[i-1] != 0) and (x[i+1]-x[i-1] != 0) and (y[i+1]-y[i-1] != 0) and (z[i+1]-z[i-1] != 0):
                    dxdt = (x[i+1]-x[i-1])/(t[i+1]-t[i-1])
                    dydt = (y[i+1]-y[i-1])/(t[i+1]-t[i-1])
                    dzdt = (z[i+1]-z[i-1])/(t[i+1]-t[i-1])

                    vt.append(averageT)
                    vx.append(dxdt)
                    vy.append(dydt)
                    vz.append(dzdt)
                else:
                    vt.append(averageT)
                    vx.append(0)
                    vy.append(0)
                    vz.append(0)
            else:
                vx.append(dx/dt)
                vy.append(dy/dt)
                vz.append(dz/dt)
                vt.append(averageT)
        
        vx = np.array(vx, dtype = np.float32)
        vy = np.array(vy, dtype = np.float32)
        vz = np.array(vz, dtype = np.float32)
        vt = np.array(vt, dtype = np.float32)

        speed = np.sqrt(vx**2 + vy**2 + vz**2)

        out = {"t": vt, "vx": vx, "vy": vy, "vz": vz, "v": speed}
        return out
    
    def _six_point_method(self):
        """
        Implementation of the six point method. See Stewart et al. (2001)

        Output:
                out: dictionary - Contains t, vx, vy, vz, v (speed)
        """

        x = self.locations['x']
        y = self.locations['y']
        z = self.locations['z']
        t = self.locations['t']

        t_reduced = t[5:-6]

        dxdt = 0.1*((x[10:-1]-x[5:-6])/(t[10:-1]-t[5:-6])) + 0.15*((x[9:-2]-x[4:-7])/(t[9:-2]-t[4:-7])) + 0.25*((x[8:-3]-x[3:-8])/(t[8:-3]-t[3:-8])) + 0.25*((x[7:-4]-x[2:-9])/(t[7:-4]-t[2:-9])) + 0.15*((x[6:-5]-x[1:-10])/(t[6:-5]-t[1:-10])) + 0.1*((x[5:-6]-x[0:-11])/(t[5:-6]-t[0:-11]))
        dydt = 0.1*((y[10:-1]-y[5:-6])/(t[10:-1]-t[5:-6])) + 0.15*((y[9:-2]-y[4:-7])/(t[9:-2]-t[4:-7])) + 0.25*((y[8:-3]-y[3:-8])/(t[8:-3]-t[3:-8])) + 0.25*((y[7:-4]-y[2:-9])/(t[7:-4]-t[2:-9])) + 0.15*((y[6:-5]-y[1:-10])/(t[6:-5]-t[1:-10])) + 0.1*((y[5:-6]-y[0:-11])/(t[5:-6]-t[0:-11]))
        dzdt = 0.1*((z[10:-1]-z[5:-6])/(t[10:-1]-t[5:-6])) + 0.15*((z[9:-2]-z[4:-7])/(t[9:-2]-t[4:-7])) + 0.25*((z[8:-3]-z[3:-8])/(t[8:-3]-t[3:-8])) + 0.25*((z[7:-4]-z[2:-9])/(t[7:-4]-t[2:-9])) + 0.15*((z[6:-5]-z[1:-10])/(t[6:-5]-t[1:-10])) + 0.1*((z[5:-6]-z[0:-11])/(t[5:-6]-t[0:-11]))

        sp = np.sqrt(dxdt**2 + dydt**2 + dzdt**2)

        velocity_data = {'t':t_reduced, 'vx':dxdt, 'vy':dydt, 'vz':dzdt, 'v':sp}
        return velocity_data
    
    def _spline_method(self, porder):
        """
        Implementation of spline interpolation in order to calculate
        instantaneous veloicty given input times of locations.  

        Output:
                out: dictionary - Contains t, vx, vy, vz, v (speed), porder, rmsdevx, rmsdevy, rmsdevz
        """

        x = self.locations['x']
        y = self.locations['y']
        z = self.locations['z']
        t = self.locations['t']

        # doubt this will lead to speedup but for time I'm leaving it becasue it works - R Perin 27-03-2024
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
            
            # if np.mod(l,10000) == 0:
            #     print('location ' + str(l) + '/' + str(len(t)))
        sp = np.sqrt(dxdt**2 + dydt**2 + dzdt**2)
        velocity_data_cart = {'t':t, 'vx':dxdt, 'vy':dydt, 'vz':dzdt, 'v':sp, 'porder':porder, 'rmsdevx':rmsdevx, 'rmsdevy':rmsdevy, 'rmsdevz':rmsdevz }
        return velocity_data_cart
    
    def fit(self, method: str, **kwargs):
        """
        Implementation of fit algorithm.  This simply allows you to 
        select what algorithm you wish to apply to the data in order
        to compute the velocities.  

        Input:
                method: str - The velocity method you wish to use (two point, six point or spline)
                **kwargs - Additional keyword arguments for porder (int)
        """

        if method == "two point":
            velocity_info = self._two_point_method()
            return velocity_info

        elif method == "six point":
            velocity_info = self._six_point_method()
            return velocity_info
        
        elif method == "spline":

            if "porder" in kwargs.keys():
                porder = kwargs["porder"]

                if not isinstance(porder, int):
                    raise TypeError("porder must be an integer!")
                velocity_info = self._spline_method(porder = porder)

                return velocity_info
            
            else: 
                raise TypeError("Please add polynomial order argument for spline velocity calculation! e.g. porder = 3")
        else:
            raise TypeError("invalid selection of method.  The avavailable methods are: \"two point\", \"six point\" or \"spline\"")

# FUNCTIONS