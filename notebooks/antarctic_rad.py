#!/usr/bin/env python
# coding: utf-8

# In[1]:


from __future__ import print_function, division


# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import qgrid

import numpy as np
from numpy import diff
import math

import pandas as pd
import xarray as xr

import climlab
from climlab.solar.insolation import daily_insolation
from climlab.radiation import DailyInsolation
from climlab.radiation import FixedInsolation
from climlab.process import TimeDependentProcess


# In[4]:


###### create a Turbulence class #######
class Turbulence(TimeDependentProcess):
    #class turbulence that decays with height computed from:
    # theta = Tatm(P0/Patm)^R/Cp
    # k0 = Net surface flux/(dtheta/dz) at surface 
    # ka = k0*exp(-z/d) in atmosphere
    # -ka * dtheta/dz is turbulent flux in the atmosphere 
    
    def __init__(self, **kwargs):
        #define initial state
        super(Turbulence, self).__init__(**kwargs)
        self.experiment = experiment
        self.time_type = 'explicit'
        self._hr = {}
        for var in self.state:
            self._hr[var] = 0. *self.state[var]
        ### define new levels
        #n+1 pressure level (lev_boundsin climlab) reworked so that we are no longer using 0 and 1000 as our TOA and surface
        self.lev_bounds[0] = self.lev_bounds[1]-2*(self.lev_bounds[1] - self.lev[0])
        self.lev_bounds[-1] = self.lev_bounds[-2]+2*(self.lev[-1] - self.lev_bounds[-2])
        #altitude (n)
        self._z = (self.experiment['Altitude'].sel(month = m))- self.experiment['Altitude'].sel(month = m, level=0)
        self._z_ground = self.experiment['Altitude'].sel(month = m, level = 0)
        #altitude bounds (n+1)
        self._z_bounds = np.zeros(len(state['Tatm']) +1)
        self._z_bounds[0] =  (ds['Altitude'].sel(month = m, level = len(self.experiment.level)+1)) -self._z_ground
        for n in range(len(state['Tatm'])-1):
            self._z_bounds[n+1] = ((self.experiment['Altitude'].sel(month = m)[n]+self.experiment['Altitude'].sel(month = m)[n+1])/2) - self._z_ground
        self._z_bounds[-1] = -self._z_bounds[-2]
        ### calculate turbulent flux
        #surface flux (LW + SW)
        self._total_sfc_flux_init = (rad.diagnostics['LW_flux_net_clr'] + rad.diagnostics['SW_flux_net_clr'])[-1]
        #gas constant / specific heat of air
        R = 287 #gas constant of air J/kg*K (at 250 K)
        cp = 1003 #specific heat of air (at 250 K) J/kg*K
        self.R_cp = R/cp #theta (length of n)
        self._theta_init = self.state['Tatm']*((self.lev[-1]/self.lev)**(self.R_cp)) 
        #dtheta_dz_init (just need the first value of dtheta/dz)
        self._dtheta_dz_surf_init = (np.diff(self._theta_init)/np.diff(self._z))[-1]
        #surface diffk
        self._surface_diffk = -self._total_sfc_flux_init/self._dtheta_dz_surf_init
        #atmospheric diffk
        scale_factor = 1000
        self._atm_diffk = self._surface_diffk * (np.exp(-self._z_bounds/scale_factor)) 
        

    def _compute(self):
        #constants
        density = 1.05 #density of air kg/m3
        cp = 1003 #specific heat of air (at 250 K) J/kg*K
        #surface flux (LW + SW)
        self._total_sfc_flux = (rad.diagnostics['LW_flux_net_clr'] + rad.diagnostics['SW_flux_net_clr'])[-1]
        #theta (length of n)
        self._theta = self.state['Tatm']*((self.lev[-1]/self.lev)**(self.R_cp)) 
        #dtheta_dz_init (length of n+1)
        self._dtheta_dz = np.zeros(len(state['Tatm']) +1)
        self._dtheta_dz[1:-1] = np.diff(self._theta)/np.diff(self._z)
        self._dtheta_dz[0] =  0 #TOA flux is zero
        self._dtheta_dz[-1] =  -self._total_sfc_flux/self._surface_diffk #total ground flux is zero, so turbulent flux + LW/SW flux = 0
        #calculate the turbulent flux
        self._turbulent_flux = -self._atm_diffk * self._dtheta_dz
        # calculate heating rate (flux convergence) from flux and convert into K/day (which is the heating rate output in climlab)
        self._hr = (np.diff(self._turbulent_flux)/np.diff(self._z_bounds)) * (1/(density*cp)* climlab.constants.seconds_per_day)
        tendencies = {'Tatm' : self._hr}
        return tendencies


# In[ ]:





# In[ ]:




