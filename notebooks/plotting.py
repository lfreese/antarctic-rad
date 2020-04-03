import matplotlib.pyplot as plt
import numpy as np
import climlab 
       
    
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
                    'atm_hr':r'Turbulent Heating Rate ($\frac{K}{s}$)',
                    'advection_Tatm':r'Advective Heating Rate ($\frac{K}{s}$)'}
    for idx, output in enumerate(outputs):  
        ax = fig.add_subplot(rows, columns, idx +1)
        if (output == 'TdotLW_clr' or output == 'TdotSW_clr' or output == 'Tatm' or output == 'atm_hr' or output == 'advection_Tatm'):
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

        if (output == 'TdotLW_clr' or output == 'TdotSW_clr' or output == 'Tatm' or output == 'atm_hr' or output == 'advection_Tatm'):
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
            plt.plot(x2-x1, y1, color = 'C0');
            plt.title(f'Difference between {CO2_conc1*(1e6)} ppm and {CO2_conc1*(1e6)} ppm \n {atm_process_dict[output]} in {month} at {time1/ climlab.constants.seconds_per_day} days')
        else:
            plt.plot(x1, y1, color = 'C0', label = f'{output} at {CO2_conc1*(1e6)} ppm');
            plt.plot(x2, y2, color = 'C1',label = f'{output} at {CO2_conc2*(1e6)} ppm');
            plt.title(f'{atm_process_dict[output]} in {month} at CO2 {CO2_conc1*(1e6)} ppm and {CO2_conc2*(1e6)} ppm \n at {np.round(time1 / climlab.constants.seconds_per_day, 3)} days')
            plt.legend()
        plt.ylabel('Altitude')
        plt.ylim([0,ylim])
        plt.xlabel(f'{output}')
        plt.xticks(rotation = 45)
        plt.tight_layout()
        
def single_level_plot(results_dict, output, CO2_conc1, CO2_conc2, month):
    fig, ax = plt.subplots(figsize = [10,4])
    for t in list(results_dict.keys())[0::10]:
        plt.plot(t / climlab.constants.seconds_per_day,results_dict[t][output][CO2_conc1][month],'C0.')
        plt.plot(t / climlab.constants.seconds_per_day,results_dict[t][output][CO2_conc2][month],'C1.')
        plt.xlabel('Time (days)')
        plt.xticks(rotation = 45)
        plt.ylabel(f'{output}')
        plt.legend([f'{CO2_conc1*(1e6)} ppm',f'{CO2_conc2*(1e6)} ppm'])
        

def plot_time_dif(results_dict, outputs, output_dict, CO2_conc1, time1, timesteps, month, figsize, rows, columns, ylim, diff_only=False):
    fig = plt.figure(figsize = figsize)
    atm_process_dict = {'Tatm':'Atmospheric Temperature (K)',
                    'atm_diffk':r'Atmospheric $\kappa$ ($\frac{m^2}{s}$)',
                    'dtheta_dz':r'$\frac{d\theta}{dz}$ ($\frac{K}{m}$)',
                    'theta':r'$\theta$',
                    'SW_flux_net_clr':r'Shortwave Flux ($\frac{W}{m^2}$)',
                    'LW_flux_net_clr':r'Longwave Flux ($\frac{W}{m^2}$)', 
                    'TdotLW_clr':r'Longwave Heating Rate ($\frac{K}{s}$)',
                    'TdotSW_clr':r'Shortwave Heating Rate ($\frac{K}{s}$)',
                    'atm_hr':r'Turbulent Heating Rate ($\frac{K}{s}$)',
                    'advection_Tatm':r'Advective Heating Rate ($\frac{K}{s}$)'}
    for idx, output in enumerate(outputs):
        ax = fig.add_subplot(rows, columns, idx +1)
        
        if (output == 'TdotLW_clr' or output == 'TdotSW_clr' or output == 'Tatm' or output == 'atm_hr' or output == 'advection_Tatm'):
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
            plt.plot(x1, y1, label = f'{output} at {time1 / climlab.constants.seconds_per_day} days');
        for time2 in timesteps:
            if (output == 'TdotLW_clr' or output == 'TdotSW_clr' or output == 'Tatm' or output == 'atm_hr' or output == 'advection_Tatm'):
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
                plt.plot(x2-x1, y1, label  = f'Difference in {output} at {time2 / climlab.constants.seconds_per_day} days' );
                plt.title(f'Difference between {time1 / climlab.constants.seconds_per_day} and {time2 / climlab.constants.seconds_per_day} \n {atm_process_dict[output]} in {month} at CO2 {CO2_conc1*(1e6)} ppm')
                plt.legend()

            else:
                plt.plot(x2, y2,label = f'{output} at {time2 / climlab.constants.seconds_per_day} days');
                plt.title(f'{atm_process_dict[output]} in {month} at CO2 {CO2_conc1*(1e6)} ppm')
                plt.legend()
            plt.ylabel('Altitude')
            plt.ylim([0,ylim])
            plt.xlabel(f'{output}')
            plt.xticks(rotation = 45)
            plt.tight_layout()
