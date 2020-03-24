import matplotlib.pyplot as plt
import numpy as np




#### plot diagnostics to see turbulence, theta, dtheta/dz ####
def plot_diagnostics(ram_dict, experiment, CO2, month):
    fig, ax = plt.subplots(4,2,figsize = (14,14), constrained_layout=True)
    fig.suptitle(month, fontsize =16)
    ax[0,0].plot(ram_dict[experiment][CO2][month].subprocess.Turbulence._atm_diffk, ram_dict[experiment][CO2][month].subprocess.Turbulence._z_bounds) 
    ax[0,0].set_title('Atmospheric Kappa (m^2/s)')
    ax[0,0].set_xlabel('Kappa (m^2/s)')
    ax[0,0].set_ylabel('Altitude (m)')
    
    ax[0,1].plot(ram_dict[experiment][CO2][month].subprocess.Turbulence._theta_init, ram_dict[experiment][CO2][month].subprocess.Turbulence._z)
    ax[0,1].set_title('Theta (K)')
    ax[0,1].set_xlabel('Theta (K)')
    ax[0,1].set_ylabel('Altitude (m)')

    ax[1,0].plot(ram_dict[experiment][CO2][month].subprocess.Turbulence._dtheta_dz, ram_dict[experiment][CO2][month].subprocess.Turbulence._z_bounds)
    ax[1,0].set_title('dtheta/dz (K/m)')
    ax[1,0].set_xlabel('dtheta/dz (K/m)')
    ax[1,0].set_ylabel('Altitude (m)')
    
    ax[1,1].plot(ram_dict[experiment][CO2][month].subprocess.Advection.forcing_tendencies['Tatm'], ram_dict[experiment][CO2][month].subprocess.Turbulence._z[:-1])
    ax[1,1].set_title('Advective Heating Rate (K/sec)')
    ax[1,1].set_xlabel('Advection (K/sec)')
    ax[1,1].set_ylabel('Altitude (m)')
    
    ax[2,1].plot(ram_dict[experiment][CO2][month].subprocess.Turbulence._hr, ram_dict[experiment][CO2][month].subprocess.Turbulence._z)
    ax[2,1].set_title('Turbulent Heating Rate (K/sec)')
    ax[2,1].set_xlabel('Turbulent Heating Rate (K/sec)')
    ax[2,1].set_ylabel('Altitude (m)')
    ax[2,1].set_ylim(-100,1000)

    ax[2,0].plot((np.pad(np.asarray(ram_dict[experiment][CO2][month].TdotSW_clr),pad_width = (0, 1))-np.pad(np.asarray(ram_dict[experiment][CO2][month].TdotLW_clr),pad_width = (0, 1)) - ram_dict[experiment][CO2][month].subprocess.Turbulence._hr), ram_dict[experiment][CO2][month].subprocess.Turbulence._z)
    ax[2,0].set_title('Radiative and Turbulent Heating Rate (K/day)')
    ax[2,0].set_xlabel('Radiative Heating Rate (K/day)')
    ax[2,0].set_ylabel('Altitude (m)')
    
    ax[3,0].plot((ram_dict[experiment][CO2][month].TdotLW_clr), ram_dict[experiment][CO2][month].subprocess.Turbulence._z[:-1])
    ax[3,0].set_title('LW Heating Rate (K/day)')
    ax[3,0].set_xlabel('Heating Rate (K/day)')
    ax[3,0].set_ylabel('Altitude (m)')
    
    ax[3,1].plot((ram_dict[experiment][CO2][month].TdotSW_clr), ram_dict[experiment][CO2][month].subprocess.Turbulence._z[:-1])
    ax[3,1].set_title('SW Heating Rate (K/day)')
    ax[3,1].set_xlabel('Heating Rate (K/day)')
    ax[3,1].set_ylabel('Altitude (m)')