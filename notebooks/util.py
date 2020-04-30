import pandas as pd
import xarray as xr
import numpy as np
from climlab.solar.insolation import daily_insolation

single_level_process = {'LW_sfc_clr':r'Surface Longwave Flux ($\frac{W}{m^2}$)', 
                        'Ts': 'Surface Temperature (K)', 
                        'ASRclr':r'ASR ($\frac{W}{m^2}$)', 
                        'SW_sfc_clr':r'Surface Shortwave Flux ($\frac{W}{m^2}$)', 
                        'total_sfc_flux':r'Total Radiative Surface Flux ($\frac{W}{m^2}$)', 
                        'dtheta_dz_near_surf_init':r'$\frac{d\theta}{dz}$ ($\frac{K}{m}$)', 
                        'surface_diffk':r'Surface $\kappa$ ($\frac{m^2}{s}$)', 
                        'sfc_turbulent_flux':r'Surface Turbulent Flux ($\frac{W}{m^2}$)', 
                        'turb_sfc_hr':r'Surface Turbulent Heating Rate ($\frac{K}{s}$)', 
                        'advection_Ts':r'Surface Advection ($\frac{K}{s}$)', 
                        'OLR':r'OLR ($\frac{W}{m^2}$)'}

def load_soundings():
    atm_df = pd.concat(pd.read_excel('../schmit_data/schmit_compiled_temp_data_final.xlsx', #read in atm data
                                     sheet_name=None,   # import all sheets in the excel
                                     na_values=np.nan), # convert nan to numpy nan values
                       axis=0) #concatenate along 0th axis
    gas_df = pd.concat(pd.read_excel('../schmit_data/schmit_compiled_gas_data.xlsx', #read in gas data
                                     sheet_name=None,   # import all sheets in the excel
                                     na_values=np.nan), # convert nan to numpy nan values
                       axis=0) #concatenate along 0th axis and flip
    df = pd.concat([atm_df,gas_df],axis=1, sort=True) #concatenate gas and atm data into one dataframe

    ds = xr.Dataset(df).unstack() #create dataset and unstack our coordinates
    ds = ds.rename({'dim_0_level_0': 'month'}) #rename coordinate
    ds = ds.rename({'dim_0_level_1': 'level'}) #rename coordinate
    for var_idx, var_name in enumerate(ds.variables): #transpose our data
        ds[var_name] = ds[var_name].transpose()
    ds['spec_humidity']=(ds['H2O']/(1+ds['H2O'])) #add in specific humidity
    ds = ds.reindex(level=ds.level[::-1]) #reindex our data so it reads from top to bottom
    ds.__delitem__('CO2') #delete unused variables
    ds.__delitem__('CO') #delete unused variables
    ds.__delitem__('CH4') #delete unused variables
    ds.__delitem__('Height') #delete unused variables
    ds['CO2_list'] = np.array([0.,100.,200.,380.,760.,1000.,1500.])/(1e6) #list of CO2 values we will use

    ds.attrs['name'] = 'all'
    
    return ds

def add_monthly_insolation(ds):
    month_days = [31,28,31,30,31,30,31,31,30,31,30,31] #list of the number of days in a month
    day_of_month = [
        np.cumsum(month_days) - np.array(month_days),
        np.cumsum(month_days)
    ] #list of first and last day of each month
    day_num = np.arange(1,366) #number of days in a year
    month_num  = np.arange(0,12) #number of months in a year

    insolation = np.zeros(len(day_num)) #insolation for each day
    monthly_insol = np.zeros(len(month_num)) #insolation for each month
    for day in day_num: #calculate daily insolation
        insolation[day-1] = daily_insolation(lat = -90, day = day)
    for month in month_num: #calculate monthly mean insolation
        monthly_insol[month] = insolation[day_of_month[0][month]:day_of_month[1][month]].mean()
    ds['monthly_insolation'] = xr.DataArray(data = monthly_insol, dims = 'month') #add into our dataset
    
def adjust_lev(ds, nlevels = 96, name = 'strat'):
    #modify ds to includes up to nlevels
    ds.attrs['extra_Altitude'] = ds['Altitude'].sel(level=nlevels+1)
    ds = ds.sel(level=slice(nlevels,0))
    ds.attrs['name'] = name
    dropped_lev = [1,2,3,4,5,6,7,9,10,11,12,14,15,17,18,20,22,24,26]
    #modify our dataset to skip levels near the surface so that it is more spaced out
    ds = ds.drop(level = [1,2,3,4,5,6,7,9,10,11,12,14,15,17,18,20,22,24,26]) #drops to one every 100m
    #ds = ds.drop(level = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,24]) #drops to one every 300m
    #ds = ds.drop(level = dropped_lev) #drops to one every 500m for the first 1000m
    #ds_strat = ds_strat.drop(level = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]) #drops to one every 1000m
    print(f'Dropped levels {dropped_lev}')
    ds['level'] = np.flip(np.arange(ds['level'].size))
    ds['Altitude'].sel(month = 'January')
    return(ds)