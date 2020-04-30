import matplotlib.pyplot as plt
import numpy as np
import climlab 
from matplotlib.pyplot import cm
import util
    
def plot_co2_dif(results_dict, outputs, output_dict, CO2_conc1, CO2_conc2, time1, month, figsize, rows, columns, ylim, diff_only = False):
    fig = plt.figure(figsize = figsize)
    atm_process_dict = {'Tatm':'Atmospheric Temperature (K)',
                    'atm_diffk':r'Atmospheric $\kappa$ ($\frac{m^2}{s}$)',
                    'dtheta_dz':r'$\frac{d\theta}{dz}$ ($\frac{K}{m}$)',
                    'theta':r'$\theta$',
                    'SW_flux_net_clr':r'Shortwave Flux ($\frac{W}{m^2}$)',
                    'LW_flux_net_clr':r'Longwave Flux ($\frac{W}{m^2}$)', 
                    'TdotLW_clr':r'Longwave Heating Rate ($\frac{K}{s}$)',
                    'TdotSW_clr':r'Shortwave Heating Rate ($\frac{K}{s}$)',
                    'turb_atm_hr':r'Turbulent Heating Rate ($\frac{K}{s}$)',
                    'advection_Tatm':r'Advective Heating Rate ($\frac{K}{s}$)',
                   'atm_turbulent_flux':'Turbulent Flux'}
    colors = cm.twilight(np.linspace(0,1,7))
    for idx, output in enumerate(outputs):  
        ax = fig.add_subplot(rows, columns, idx +1)
        if (output == 'atm_diffk' or output == 'dtheta_dz' or output == 'atm_turbulent_flux'):
            y1 = np.asarray(results_dict[time1][output_dict['bounds'][output]][CO2_conc1][month][:-1])
        else:
            y1 = np.asarray(results_dict[time1][output_dict['bounds'][output]][CO2_conc1][month])
        yidx = y1 < ylim
        y1 = y1[yidx]
        
        if (output == 'TdotLW_clr' or output == 'TdotSW_clr'):
            x1 = np.asarray(results_dict[time1][output][CO2_conc1][month])/86400
        else: 
            x1 = np.asarray(results_dict[time1][output][CO2_conc1][month])
        x1 = x1[yidx]

        if (output == 'atm_diffk' or output == 'dtheta_dz' or output == 'atm_turbulent_flux'):
            y2 = np.asarray(results_dict[time1][output_dict['bounds'][output]][CO2_conc2][month][:-1])
        else:
            y2 = np.asarray(results_dict[time1][output_dict['bounds'][output]][CO2_conc2][month])  
        y2 = y2[yidx]

        if (output == 'TdotLW_clr' or output == 'TdotSW_clr'):
            x2 = np.asarray(results_dict[time1][output][CO2_conc2][month])/86400
        else: 
            x2 = np.asarray(results_dict[time1][output][CO2_conc2][month])
        x2 = x2[yidx]

        if diff_only:
            plt.plot(x2-x1, y1, color = colors[4], linewidth = 2, linestyle = '--');
            plt.title(f'Difference between {CO2_conc1*(1e6)} ppm and {CO2_conc2*(1e6)} ppm \n {atm_process_dict[output]} in {month} at {time1/ climlab.constants.seconds_per_day} days')
        else:
            plt.plot(x1, y1, color = colors[2], linewidth = 2, label = f'{output} at {CO2_conc1*(1e6)} ppm');
            plt.plot(x2, y2, color = colors[4], linewidth = 2,label = f'{output} at {CO2_conc2*(1e6)} ppm');
            plt.title(f'{atm_process_dict[output]} in {month} at CO2 {CO2_conc1*(1e6)} ppm and {CO2_conc2*(1e6)} ppm \n at {np.round(time1 / climlab.constants.seconds_per_day, 3)} days')
            plt.legend()
        
        plt.ylabel('Altitude')
        plt.ylim([0,ylim])
        plt.xlabel(f'{output}')
        plt.xticks(rotation = 45)
        plt.tight_layout()
        
def single_level_plot(results_dict, output_list, CO2_conc1, CO2_conc2, month_list, figsize):
    fig, axes = plt.subplots(len(output_list), len(month_list), figsize = figsize)
    colors = cm.twilight(np.linspace(0,1,7))
    for idx_m, month in enumerate(month_list):
        for idx_o, output in enumerate(output_list):
            ax = axes[idx_o, idx_m]
            for t in list(results_dict.keys())[0::10]:
                if output == 'total_sfc_flux':
                    ax.plot(t / climlab.constants.seconds_per_day,results_dict[t][output][CO2_conc1][month]*-1, c = colors[2], marker = '.', linestyle = '')
                    ax.plot(t / climlab.constants.seconds_per_day,results_dict[t][output][CO2_conc2][month]*-1, c = colors[5], marker = '.', linestyle = '')
                else:
                    ax.plot(t / climlab.constants.seconds_per_day,results_dict[t][output][CO2_conc1][month], c = colors[2], marker = '.', linestyle = '')
                    ax.plot(t / climlab.constants.seconds_per_day,results_dict[t][output][CO2_conc2][month], c = colors[5], marker = '.', linestyle = '')
                axes[-1,idx_m].set_xlabel('Time (days)', fontsize = 14)
                #ax.set_xticks(rotation = 45)
                axes[idx_o, 0].set_ylabel(f'{util.single_level_process[output]}', fontsize = 14, wrap=True)
                axes[0,idx_m].set_title(month, fontsize = 20)
        fig.legend([f'{CO2_conc1*(1e6)}',f'{CO2_conc2*(1e6)}'], title = 'ppm', bbox_to_anchor=(1.1, .95))
    plt.tight_layout()


def plot_time_dif(results_dict, output, output_dict, CO2_conc1, time1, timesteps, month, figsize, ylim, diff_only=False):
    fig = plt.figure(figsize = figsize)
    atm_process_dict = {'Tatm':'Atmospheric Temperature (K)',
                    'atm_diffk':r'Atmospheric $\kappa$ ($\frac{m^2}{s}$)',
                    'dtheta_dz':r'$\frac{d\theta}{dz}$ ($\frac{K}{m}$)',
                    'theta':r'$\theta$',
                    'SW_flux_net_clr':r'Shortwave Flux ($\frac{W}{m^2}$)',
                    'LW_flux_net_clr':r'Longwave Flux ($\frac{W}{m^2}$)', 
                    'TdotLW_clr':r'Longwave Heating Rate ($\frac{K}{s}$)',
                    'TdotSW_clr':r'Shortwave Heating Rate ($\frac{K}{s}$)',
                    'turb_atm_hr':r'Turbulent Heating Rate ($\frac{K}{s}$)',
                    'advection_Tatm':r'Advective Heating Rate ($\frac{K}{s}$)',
                       'atm_turbulent_flux': r'Atmospheric Turbulent Flux ($\frac{W}{m^2}$)'}
    ax = fig.add_subplot()
    color = iter(cm.twilight(np.linspace(0,1,len(timesteps)+1)))
    if (output == 'atm_diffk' or output == 'dtheta_dz' or output == 'atm_turbulent_flux'):
        y1 = np.asarray(results_dict[time1][output_dict['bounds'][output]][CO2_conc1][month][:-1])
    else:
        y1 = np.asarray(results_dict[time1][output_dict['bounds'][output]][CO2_conc1][month])          
    yidx = y1 < ylim
    y1 = y1[yidx]
        
    if (output == 'TdotLW_clr' or output == 'TdotSW_clr'):
        x1 = np.asarray(results_dict[time1][output][CO2_conc1][month])/86400
    else: 
        x1 = np.asarray(results_dict[time1][output][CO2_conc1][month])
    x1 = x1[yidx]
    if diff_only == False:
        plt.plot(x1, y1, label = f'{np.round(time1 / climlab.constants.seconds_per_day, 0)}');
    for time2 in timesteps:
        c=next(color)
        if (output == 'atm_diffk' or output == 'dtheta_dz' or output == 'atm_turbulent_flux'):
            y2 = np.asarray(results_dict[time2][output_dict['bounds'][output]][CO2_conc1][month][:-1])
        else:
            y2 = np.asarray(results_dict[time2][output_dict['bounds'][output]][CO2_conc1][month])
        y2 = y2[yidx]

        if (output == 'TdotLW_clr' or output == 'TdotSW_clr'):
            x2 = np.asarray(results_dict[time2][output][CO2_conc1][month])/86400
        else: 
            x2 = np.asarray(results_dict[time2][output][CO2_conc1][month])
        x2 = x2[yidx]

        if diff_only:
            plt.plot(x2-x1, y1, c=c, label  = f'Difference at {np.round(time2 / climlab.constants.seconds_per_day, 0)} days' );
            plt.title(f'Difference between {np.round(time1 / climlab.constants.seconds_per_day, 0)} and {np.round(time2 / climlab.constants.seconds_per_day, 0)} \n {atm_process_dict[output]} in {month} at CO2 {CO2_conc1*(1e6)} ppm')
            plt.legend()

        else:
            plt.plot(x2, y2, c=c, label = f'{np.round(time2 / climlab.constants.seconds_per_day, 0)}');
            plt.title(f'{atm_process_dict[output]} in {month} at CO2 {CO2_conc1*(1e6)} ppm')
            plt.legend(loc='best', title = 'Days')
        plt.ylabel('Altitude')
        plt.ylim([0,ylim])
        plt.xlabel(f'{output}')
        plt.xticks(rotation = 45)
        plt.tight_layout()

            

def plot_temp(ds, results_dict):
    fig, ax = plt.subplots(figsize = [6,8])
    color=iter(cm.twilight(np.linspace(0,1,13)))
    for idx, month in enumerate(ds['month'].values):
        c=next(color)
        x = np.asarray(results_dict[0]['Tatm'][0.00038][month])
        x = np.append(x, results_dict[0]['Ts'][0.00038][month])
        y = results_dict[0]['z'][0.00038][month]/1000
        y = np.append(y, 0.)
        plt.plot(x, y, c = c, label = month, lw = 2)
        plt.xlabel('Temperature (K)', fontsize = 14)
        plt.ylabel('Altitude relative to surface level (km)', fontsize = 14)
        plt.yscale('log')
        ax.yaxis.set_ticklabels([0, 0,.1,1, 10])
        plt.title('Monthly Temperature Profiles', fontsize = 20)
        plt.legend(title = 'Months')

def plot_temp_timestepped(ds, results_dict, month, CO2, time1, timesteps, ylim, diff_only):
    fig, ax = plt.subplots(figsize = [6,8])
    color=iter(cm.twilight(np.linspace(0,1,len(timesteps)+2)))
    y1 = results_dict[time1]['z'][CO2][month]
    y1 = np.append(y1, 0.)
    yidx = y1 < ylim
    y1 = y1[yidx]
    x1 = np.asarray(results_dict[time1]['Tatm'][CO2][month])
    x1 = np.append(x1, results_dict[time1]['Ts'][CO2][month])
    x1 = x1[yidx]
    c=next(color)
    if diff_only == False:
        plt.plot(x1, y1, c = c, label = np.round(time1 / climlab.constants.seconds_per_day, 0), lw = 2)
    for idx, time in enumerate(timesteps):
        c=next(color)
        y2 = results_dict[time]['z'][CO2][month]
        y2 = np.append(y2, 0.)
        yidx = y2 < ylim
        y2 = y2[yidx]
        x2 = np.asarray(results_dict[time]['Tatm'][CO2][month])
        x2 = np.append(x2, results_dict[time]['Ts'][CO2][month])
        x2 = x2[yidx]
        if diff_only:
                plt.plot(x2-x1, y2, c = c, label = np.round(time / climlab.constants.seconds_per_day, 0), lw = 2)
        else:
            plt.plot(x2, y2, c = c, label = np.round(time / climlab.constants.seconds_per_day, 0), lw = 2)

    plt.xlabel('Temperature (K)', fontsize = 14)
    plt.ylabel('Altitude relative to surface level (km)', fontsize = 14)
    plt.title(f'Temperature Profiles {month}', fontsize = 20)
    plt.legend(title="Days")
        
def plot_turbulent_hr(ds, results_dict, month, CO2, timesteps, ylim):
    fig, ax = plt.subplots(figsize = [6,8])
    color=iter(cm.twilight(np.linspace(0,1,len(timesteps)+1)))
    for idx, time in enumerate(timesteps):
        c=next(color)
        y = results_dict[time]['z_bounds'][CO2][month]
        yidx = y < ylim
        y = y[yidx]
        x = np.asarray(results_dict[time]['turb_atm_hr'][CO2][month])
        x = np.append(x, results_dict[time]['turb_sfc_hr'][CO2][month])
        x = x[yidx]
        plt.plot(x, y, c = c, label = np.round(time / climlab.constants.seconds_per_day, 0), lw = 2)
        plt.xlabel('Turbulent Heating Rate', fontsize = 14)
        plt.ylabel('Altitude relative to surface level (km)', fontsize = 14)
        plt.title(f'Turbulent Heating Rate Profiles {month}', fontsize = 20)
        plt.legend(title="Days")
        plt.xticks(rotation = '45')

def plot_turbulent_flux(ds, results_dict, month, CO2, timesteps, ylim):
    fig, ax = plt.subplots(figsize = [6,8])
    color=iter(cm.twilight(np.linspace(0,1,len(timesteps)+1)))
    for idx, time in enumerate(timesteps):
        c=next(color)
        y = results_dict[time]['z_bounds'][CO2][month]
        yidx = y < ylim
        y = y[yidx]
        x = np.asarray(results_dict[time]['atm_turbulent_flux'][CO2][month])
        x = np.append(x, results_dict[time]['sfc_turbulent_flux'][CO2][month])
        x = x[yidx]
        plt.plot(x, y, c = c, label = np.round(time / climlab.constants.seconds_per_day, 0), lw = 2)
        plt.xlabel('Turbulent Flux', fontsize = 14)
        plt.ylabel('Altitude relative to surface level (km)', fontsize = 14)
        plt.title(f'Turbulent Flux Profiles {month}', fontsize = 20)
        plt.legend(title="Days")
        plt.xticks(rotation = '45')

def plot_sfc_TOA_process(ds, results_dict, process, months, single_level_process):
    fig, axes = plt.subplots(1,2, figsize = [12,4])
    color=iter(cm.twilight(np.linspace(0,1,8)))
    for CO2_conc in results_dict[0][process].keys():
        c=next(color)
        for idx, month in enumerate(months):
            ax = axes[idx]
            x = CO2_conc*10e5
            y = results_dict[0][process][CO2_conc][month]
            ax.plot(x, y,'o', c = c,mec ='k',ms = 9)
            ax.set_xlabel('$CO_2$ Concentration (ppm)', fontsize = 14)
            ax.set_ylabel(single_level_process[process], fontsize = 14)
            plt.tight_layout()
            ax.set_title(f'{month}', fontsize = 14)
            
            
def plot_sum_hr(results_dict, month, CO2_conc, timesteps):
    fig = plt.figure(figsize = [8,8])
    ax = fig.add_subplot(1, 1, 1)
    color=iter(cm.twilight(np.linspace(0,1,len(timesteps)+1)))
    for time in timesteps:
        c = next(color)
        x1 = np.asarray(results_dict[time]['TdotLW_clr'][CO2_conc][month]/86400) + np.asarray(results_dict[time]['TdotSW_clr'][CO2_conc][month]/86400) + np.asarray(results_dict[time]['turb_atm_hr'][CO2_conc][month]) + np.asarray((results_dict[time]['advection_Tatm'][CO2_conc][month]))

        y1 = results_dict[time]['z_bounds'][CO2_conc][month][:-1]

        plt.plot(x1, y1, c= c, label = np.round(time/climlab.constants.seconds_per_day, 0));

        plt.legend(title = 'Days')
    plt.ylabel('Altitude', fontsize = 14)
    plt.xlabel(f'Heating Rate (K/s)', fontsize = 14)
    plt.title(f'Sum of heating rates in {month} at CO2 {CO2_conc * 1e6} ppm', fontsize = 16);

    
def plot_adv_LW_hr(results_dict, month, CO2_conc1, CO2_conc2, time, diff_only = False):
    fig = plt.figure(figsize = [8,8])
    ax = fig.add_subplot(1, 1, 1)
    color=cm.twilight(np.linspace(0,1,5))
    x1 = np.asarray(results_dict[time]['TdotLW_clr'][CO2_conc1][month]/86400) + np.asarray((results_dict[time]['advection_Tatm'][CO2_conc1][month]))
    y1 = results_dict[time]['z_bounds'][CO2_conc1][month][:-1]
        
    x2 = np.asarray(results_dict[time]['TdotLW_clr'][CO2_conc2][month]/86400) + np.asarray((results_dict[time]['advection_Tatm'][CO2_conc2][month]))

    y2 = results_dict[time]['z_bounds'][CO2_conc2][month][:-1]
    if diff_only:
        plt.plot(x2-x1, y2, c= color[2], label = f"Difference between {CO2_conc1}, {CO2_conc2}")
    else:
        plt.plot(x1, y1, c= color[2], label = CO2_conc1);
        plt.plot(x2, y2, c= color[4], label = CO2_conc2);
    plt.legend(title = 'Days')
    plt.ylabel('Altitude', fontsize = 14)
    plt.xlabel(f'Heating Rate (K/s)', fontsize = 14)
    plt.xticks(rotation = '45')
    plt.title(f'Sum of LW and Adv HR in {month} at {time} ppm', fontsize = 16);
    
    
def rad_sfc_HR_plot(rad_sfc_HR, results_dict, CO2_conc1, CO2_conc2, month_list, figsize):

    fig = plt.figure(figsize = figsize)
    colors = cm.twilight(np.linspace(0,1,7))
    for idx_m, month in enumerate(month_list):
        plt.subplot(1, 2, idx_m+1)
        for t in list(results_dict.keys())[0::10]:
            plt.plot(t / climlab.constants.seconds_per_day,-rad_sfc_HR[t][CO2_conc1][month], c = colors[2], marker = '.', linestyle = '')
            plt.plot(t / climlab.constants.seconds_per_day,-rad_sfc_HR[t][CO2_conc2][month], c = colors[5], marker = '.', linestyle = '')
            plt.title(f'{month} Surface Radiative Heating Rate')
    fig.legend([f'{CO2_conc1*(1e6)}',f'{CO2_conc2*(1e6)}'], title = 'ppm', bbox_to_anchor=(1.1, .95))
    plt.tight_layout()