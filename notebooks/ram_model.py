import xarray as xr
import climlab
from attrdict import AttrDict
from climlab.process import TimeDependentProcess
import numpy as np
import warnings

class Turbulence(TimeDependentProcess):
    """
    class turbulence that decays with height computed from:
    theta = Tatm(P0/Patm)^R/Cp
    k0 = Net surface flux/(dtheta/dz) at surface 
    ka = k0*exp(-z/d) in atmosphere
    -ka * dtheta/dz is turbulent flux in the atmosphere
    Surface level is at lev_bound and z_bound [-1]; lowest atm level at z and lev [-1]
    """    
    def __init__(self, surface_diffk = None, **kwargs):
        #define initial state
        super(Turbulence, self).__init__(**kwargs) #initialize our model
        self.dataset = kwargs.get('ds') #name the dataset as part of the model
        self.time_type = 'explicit' #explicit time calculations
        ### add new diagnostics
        self.add_diagnostic('total_sfc_flux_init', 0. * self.Ts)
        self.add_diagnostic('total_sfc_flux', 0. * self.Ts)
        self.add_diagnostic('theta_init', 0. * (self.Tatm))
        self.add_diagnostic('dtheta_dz_surf_init', 0. * self.Ts)
        self.add_diagnostic('surface_diffk', 0. * self.Ts)
        self.add_diagnostic('theta', 0. * (self.Tatm))
        self.add_diagnostic('dtheta_dz', 0. * (self.Tatm))
        self.add_diagnostic('atm_diffk', 0. * (self.Tatm))
        self.add_diagnostic('atm_turbulent_flux', 0. * (self.Tatm+1))
        self.add_diagnostic('turb_atm_hr', 0. * self.Tatm)
        self.add_diagnostic('turb_ground_hr', 0. * self.Ts)
        self.add_diagnostic('turb_hr', 0. * (self.Tatm+1))
        self.add_diagnostic('z', 0. * self.Tatm)
        self.add_diagnostic('z_bounds', 0. * (self.Tatm +1)) #include z_bound lowest as surface level
        
        ### read in our kwargs
        m = kwargs.get('m')
        state = kwargs.get('state')
        rad = kwargs.get('rad')
        self.rad = rad
        self.state = state
        self.m = m
        z_ground = -.5
        ### define new levels
        #n+1 pressure level (lev_bounds in climlab) reworked so that we are no longer using 0 and 1000 as our TOA and surface
        self.lev_bounds[0] = self.lev_bounds[1]-2*(self.lev_bounds[1] - self.lev[0])
        self.lev_bounds[-1] =  self.lev_bounds[-2]+2*(self.lev[-1] - self.lev_bounds[-2])
        #altitude (n)
        self.z[:-1] = (self.dataset['Altitude'].sel(month = m) - self.dataset['Altitude'].sel(month = m).isel(level = -1))[:-1]
        self.z[-1] = self.z[-2]/2
        self.z[:-1] += self.z[-1]
        #altitude bounds (n)
        self.z_bounds = np.zeros(len(self.state['Tatm']) + 1)
        for n in range(len(self.state['Tatm'])-1):
            self.z_bounds[n+1] = ((self.z[n]+self.z[n+1])/2)
        self.z_bounds[0] = (
            (self.dataset.attrs['extra_Altitude'].sel(month = m) + 
             self.dataset['Altitude'].sel(month = m).sel(level = len(self.dataset['Altitude'].sel(month = m))-1))/2 - 
            self.dataset['Altitude'].sel(month = m).isel(level = -1)
        ) + self.z[-1]
        self.z_bounds[-1] = 0
        
        ### calculate turbulent flux
        #surface flux (LW + SW)
        self.total_sfc_flux_init = (rad.diagnostics['SW_sfc_clr']- rad.diagnostics['LW_sfc_clr']) #W/m^2
        #gas constant / specific heat of air
        R = 287 #gas constant of air J/kg*K (at 250 K)
        cp_air = 1004 #specific heat of air (at 250 K) J/(kg*K)
        dp_sfc = (self.lev_bounds[-1] - self.lev[-1]) #pressure difference between lowest bounds (surface) and the lowest atm layer
        dz_sfc = (self.z_bounds[-1] - self.z[-1]) #altitude difference between lowest bounds (surface) and the lowest atm layer
        self.R_cp = R/cp_air 
        self.theta_init = self.state['Tatm']*(np.asarray(self.dataset['Pressure'].sel(month = m).isel(level = -1))/self.lev)**(self.R_cp)#K
        #dtheta_dz_init (just need the first value of dtheta/dz)
        self.dtheta_dz_surf_init = (self.theta_init[-1]-self.Ts)/(dz_sfc) #K/m
        if surface_diffk == None:
            self.surface_diffk = (-self.total_sfc_flux_init/self.dtheta_dz_surf_init)/(climlab.utils.heat_capacity.atmosphere(dp_sfc)/(dz_sfc)) #m^2/s
        else:
            self.surface_diffk = surface_diffk
        #atmospheric diffk (n)
        scale_factor = 500 #100 #m
        self.atm_diffk = np.ones_like(self.rad.diagnostics['LW_flux_down'])
        self.atm_diffk = self.surface_diffk * (np.exp(-(self.z_bounds)[:-1]/scale_factor)) #m^2/s

    def _compute(self):
        #constants
        dz_ground = 1.
        z_ground = -.5
        dp_sfc = (self.lev[-1] - self.lev_bounds[-1]) #pressure difference between lowest bounds (surface) and the lowest atm layer
        dz_sfc = (self.z[-1] - self.z_bounds[-1]) #altitude difference between lowest bounds (surface) and the lowest atm layer
        cp_air = 1004.  #specific heat of air (at 250 K) J/(kg*K)
        g = 9.81 #m/s2
        #surface flux (LW + SW)
        self.total_sfc_flux = (self.rad.diagnostics['SW_sfc_clr'] - self.rad.diagnostics['LW_sfc_clr']) #W/m^2
        #theta (length of n)
        self.theta = self.state['Tatm'] * (np.asarray(self.dataset['Pressure'].sel(month = self.m).isel(level = -1))/self.lev)**(self.R_cp)#K
        #dtheta_dz_init (length of n+1)
        self.dtheta_dz = np.zeros_like(self.state['Tatm'])
        self.dtheta_dz[1:] = np.diff(self.theta)/np.diff(self.z)
        self.dtheta_dz[0] =  0 #TOA flux is zero so we set this to zero so that when they are multiplied TOA flux = 0
        #calculate the atmospheric turbulent flux
        dz_toa = (self.z[0] + (self.z[0] - self.z[1])) - self.z[0]
        dp_toa = (self.lev[0] + (self.lev[0] - self.lev[1])) - self.z[0]
        # calculate heating rate (flux convergence) from flux and convert into K/sec (which is the heating rate output in climlab)
        self.atm_turbulent_flux = np.zeros_like(self.lev_bounds)
        self.atm_turbulent_flux[1:-1] = -self.atm_diffk[1:] * self.dtheta_dz[1:] *(-climlab.utils.heat_capacity.atmosphere(np.diff(self.lev))/(np.diff(self.z))) #W/m^2
        self.atm_turbulent_flux[0] = -self.atm_diffk[0] * self.dtheta_dz[0] *(-climlab.utils.heat_capacity.atmosphere(dp_toa)/(dz_toa))
        self.atm_turbulent_flux[-1] = np.squeeze(-self.surface_diffk * ((self.theta[-1]-self.Ts)/(dz_sfc)) *(-climlab.utils.heat_capacity.atmosphere(dp_sfc)/(dz_sfc))) #W/m^2
        
        #calculate heating rate (flux convergence) from flux and convert into K/sec (which is the heating rate output in climlab)
        self.turb_atm_hr = np.zeros_like(self.state['Tatm'])
        #ignore bottom z_bounds since this is the surface level
        self.turb_atm_hr = -np.diff(self.atm_turbulent_flux)/(-climlab.utils.heat_capacity.atmosphere(np.diff((self.lev_bounds)))) #K/sec 
        self.turb_ground_hr= [-(np.asarray(self.atm_turbulent_flux[-1] - 0)/dz_ground)*(dz_ground/(climlab.utils.heat_capacity.ocean(dz_ground)))] #K/sec
        self.turb_hr = np.concatenate([self.turb_atm_hr,self.turb_ground_hr])
        tendencies = {'Tatm' : self.turb_atm_hr, 'Ts' : self.turb_ground_hr}
        
        # check that CFL condition is met
        self._CFL = self.surface_diffk*(self.time['timestep']/(np.diff(self.z_bounds)**2))[-1] 
        if self._CFL > 1.:
            warnings.warn(f"CFL Condition not met, {self._CFL}, timestep too large or lower level z difference too small for {self.m}, CO2 kg/kg")
        return tendencies

def init_ram(
        ds, m, CO2, timestep, turbulence_on, advection_on, advection = None,
        surface_diffk = None, albedo = .8
    ):
    #create two domains: atm and surface
    sfc, atm=climlab.domain.single_column(lev=ds['Pressure'].sel(month=m).values, water_depth=1.);
    #change the level bounds
    atm.axes['lev'].bounds[0] = atm.axes['lev'].bounds[1]-2*(atm.axes['lev'].bounds[1] - atm.axes['lev'].points[0])
    atm.axes['lev'].bounds[-1] = atm.axes['lev'].bounds[-2]+2*(atm.axes['lev'].points[-1] - atm.axes['lev'].bounds[-2])
    #change the level delta
    atm.axes['lev'].delta = atm.axes['lev'].bounds[1:] - atm.axes['lev'].bounds[0:-1]
    #update heat capacity
    atm.set_heat_capacity()
    #create an atmospheric state
    state = AttrDict()
    #set up a surface temperature profile 
    Ts_dict = {}
    for month in ds['month'].values:
        Ts_dict[month] = (ds['Temperature'].sel(month = month)[-1] - 52*(ds['Temperature'].sel(month = month)[-2] - ds['Temperature'].sel(month = month)[-1])/
                     (ds['Altitude'].sel(month = month)[-2] - ds['Altitude'].sel(month = month)[-1])).values    
    T_s=climlab.Field(Ts_dict[m], domain=sfc);
    state['Ts']=T_s #K
    #set up an atmospheric temperature profile
    T_atm=climlab.Field(ds['Temperature'].sel(month=m).values, domain=atm);
    state['Tatm']=T_atm #K

    #radiation model setup
    rad = climlab.radiation.RRTMG(name='Radiation(all gases)', state = state,
                                      specific_humidity = ds['spec_humidity'].sel(month = m).values,
                                      albedo = albedo,
                                      timestep = timestep,
                                      ozone_file = None,
                                      S0 = 1365.2,
                                      insolation = ds['monthly_insolation'].sel(month = m).values,
                                      isolvar = 1 #see https://climlab.readthedocs.io/en/latest/_modules/climlab/radiation/rrtm/rrtmg_sw.html#RRTMG_SW._compute_heating_rates 
                                                  #1 = Solar variability (using NRLSSI2  solar
                                                    # model) with solar cycle contribution
                                                    # determined by fraction of solar cycle
                                                    # with facular and sunspot variations
                                                    # fixed to their mean variations over the
                                                    # average of Solar Cycles 13-24;
                                                    # two amplitude scale factors allow
                                                    # facular and sunspot adjustments from
                                                    # mean solar cycle as defined by indsolvar
                                    )
    rad.absorber_vmr['O3'] = ds['O3'].sel(month = m).values #kg/kg
    rad.absorber_vmr['CO2'] = CO2 #kg/kg
    #create ram
    ram = climlab.TimeDependentProcess(state = state, timestep = timestep)
    #add latitude axis
    lat = climlab.domain.axis.Axis(axis_type='lat', points=-90.)
    ram.domains['Ts'].axes['lat'] = lat
    #add radiation
    ram.add_subprocess('Radiation', rad)
    #compute ram
    ram.compute()
    if turbulence_on == True:
        #turbulence model setup (coupled to rad model, ds, and month)
        turb = Turbulence(surface_diffk = surface_diffk, name = 'Turbulence', state=state, rad = rad, m = m, ds = ds, timestep = timestep)
        #add turbulence
        ram.add_subprocess('Turbulence', turb) #add insolation subprocess
        #compute ram
        ram.compute()
    if turbulence_on == False:
        turb = None
    if advection_on:
        #advective model setup (coupled to rad model)
        adv = climlab.process.external_forcing.ExternalForcing(state = state, ram = ram, turb = turb)
        if turbulence_on:
            normal_advection = -((ram.TdotSW_clr + ram.TdotLW_clr)/climlab.constants.seconds_per_day + turb.turb_atm_hr) #(K/day + K/day)/(sec/day) + K/sec
        if advection == None:
            adv.forcing_tendencies['Tatm'] = normal_advection 
        else:
            adv.forcing_tendencies['Tatm'] = np.copy(advection[m])
        #add advection
        ram.add_subprocess('Advection', adv)
        #compute ram
        ram.compute()
    
    return ram


























def init_ram_no_advection(
        ds, m, CO2, timestep,
        surface_diffk = None, albedo = .8
    ):
    #create two domains: atm and surface
    sfc, atm=climlab.domain.single_column(lev=ds['Pressure'].sel(month=m).values, water_depth=1.);
    #create an atmospheric state
    state = AttrDict()
    #set up a surface temperature profile 
    Ts_dict = {}
    for month in ds['month'].values:
        Ts_dict[month] = (ds['Temperature'].sel(month = month)[-1] - 52*(ds['Temperature'].sel(month = month)[-2] - ds['Temperature'].sel(month = month)[-1])/
                     (ds['Altitude'].sel(month = month)[-2] - ds['Altitude'].sel(month = month)[-1])).values    
    T_s=climlab.Field(Ts_dict[m], domain=sfc);
    state['Ts']=T_s #K
    #set up an atmospheric temperature profile
    T_atm=climlab.Field(ds['Temperature'].sel(month=m).values, domain=atm);
    state['Tatm']=T_atm #K

    #radiation model setup
    rad = climlab.radiation.RRTMG(name='Radiation(all gases)', state = state,
                                      specific_humidity = ds['spec_humidity'].sel(month = m).values,
                                      albedo = albedo,
                                      timestep = timestep,
                                      ozone_file = None,
                                      S0 = 1365.2,
                                      insolation = ds['monthly_insolation'].sel(month = m).values,
                                      isolvar = 1 #see https://climlab.readthedocs.io/en/latest/_modules/climlab/radiation/rrtm/rrtmg_sw.html#RRTMG_SW._compute_heating_rates 
                                                  #1 = Solar variability (using NRLSSI2  solar
                                                    # model) with solar cycle contribution
                                                    # determined by fraction of solar cycle
                                                    # with facular and sunspot variations
                                                    # fixed to their mean variations over the
                                                    # average of Solar Cycles 13-24;
                                                    # two amplitude scale factors allow
                                                    # facular and sunspot adjustments from
                                                    # mean solar cycle as defined by indsolvar
                                    )
    rad.absorber_vmr['O3'] = ds['O3'].sel(month = m).values #kg/kg
    rad.absorber_vmr['CO2'] = CO2 #kg/kg
    #create ram
    ram = climlab.TimeDependentProcess(state = state, timestep = timestep)
    #add latitude axis
    lat = climlab.domain.axis.Axis(axis_type='lat', points=-90.)
    ram.domains['Ts'].axes['lat'] = lat
    #add radiation
    ram.add_subprocess('Radiation', rad)
    #compute ram
    ram.compute()
    #turbulence model setup (coupled to rad model, ds, and month)
    turb = Turbulence(surface_diffk = surface_diffk, name = 'Turbulence', state=state, rad = rad, m = m, ds = ds, timestep = timestep)
    #add turbulence
    ram.add_subprocess('Turbulence', turb) #add insolation subprocess
    #compute ram
    ram.compute()
    
    return ram
def annual_mean_sfc_diffk(ds, ram_dict):
    sum_surface_diffk = 0
    months_for_mean = 0
    for m in ds['month'].values:
        if ram_dict[0.00038][m].surface_diffk >0:
            months_for_mean +=1
            sum_surface_diffk += ram_dict[0.00038][m].surface_diffk
    average_surface_diffk = sum_surface_diffk/months_for_mean
    return average_surface_diffk

def fill_ensemble(ds, ram_dict, timestep, turbulence_on, advection_on, advection, surface_diffk):
    for CO2 in ds['CO2_list'].values:
        ram_dict[CO2] = {}
        for m in np.asarray(ds['month']):
            ram = init_ram(ds = ds, m = m, CO2 = CO2, timestep = timestep, turbulence_on = turbulence_on, advection_on = advection_on, advection = advection, surface_diffk = surface_diffk)
            ram_dict[CO2][m] = {}
            ram_dict[CO2][m] = ram
    return ram_dict


def fill_ensemble_no_advection(ds, ram_dict, timestep, surface_diffk):
    for CO2 in ds['CO2_list'].values:
        ram_dict[CO2] = {}
        for m in np.asarray(ds['month']):
            ram = init_ram_no_advection(ds = ds, m = m, CO2 = CO2, timestep = timestep, surface_diffk = surface_diffk)
            ram_dict[CO2][m] = {}
            ram_dict[CO2][m] = ram
    return ram_dict