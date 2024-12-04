# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 16:26:14 2024

@author: ljh3218
"""

#%% Import things
from Lattice import lattice
from Lattice import Parameters
from datetime import datetime
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 

now = datetime.now()
pwd = os.getcwd()

save_path = pwd + '/plots/' + now.strftime('%d%b%Y-1/') # where you want to save them
if not os.path.exists(save_path): os.makedirs(save_path)

#%% Calculate states vs parameter of interest for 0 field
parameter_to_vary = 't0'
parameter_array = [1e-3, 2e-3, 3e-3, 4e-3, 5e-3]
labels = ['1 meV', '2 meV', '3 meV', '4 meV', '5 meV']

L0_dict = {}

spacing = 10   # Angstrom
num_sites_coupled =  1.45 #1.45, not one here as I want to count diagonal neighbours as CT states
size = 4

j0 = 1.5
r0j = 0.1
min_dist = min([size, r0j*(j0/0.0257 - 1)])

save = 1
for i in range(len(parameter_array)):
#kout is the extraction rate of the separted charges
#distance over which charges are considered 'separated' is determined by dist_CS_min parameter in L0.buildHam
    params = Parameters(Temp = 300, 
                    Epeak = 0.16,
                    kout = 1e11,
                    reorE_inner = 0.202, reorE_outer = 0.048,
                    j0 = j0, r0j = r0j*spacing, 
                    E_singlet = 1.4,
                    recom = 1)

#Build the lattice
    L0 = lattice()
    L0.generate_uniform(size = size,
                        HOMO = 0, LUMO = 1.8,
                        dist_sites_A = spacing, 
                        min_dist_near_neighbour = num_sites_coupled*spacing + 0.01,
                        t0_homo = parameter_array[i], t0_lumo = parameter_array[i],
                        d0 = 5e-3, r0d = 0.1*spacing,
                        V_ex = 0.02, V_CT = 0.001,
                        CT_abs_prop = 0, recom = 1)
    
    time_start = datetime.now()
    L0.buildHam(params,
                F = [0,0,0],
                min_dist_near_neighbour = (num_sites_coupled*spacing) + 0.01,
                dist_CS_min = min_dist*spacing,
                disorder_site_ene = 0.05,
                random_seed = 42)
    time_ham = datetime.now()
    print('\n' + 'build ham: ' + str(time_ham - time_start))   
    L0.states_from_Ham_v2(params, max_energy_diff = 1.5)
    time_states = datetime.now()
    print('get states: ' + str(time_states - time_ham))
    #calculate the rates of transitions between states using Redfield 
    L0.rates_mat_2D_v7(params) 
    time_rates = datetime.now()
    print('get rates: ' + str(time_rates - time_states))
    #solve for steady state population of the different states 
    L0.solve_steady(params)
    time_steady = datetime.now()
    print('solve steady: ' + str(time_steady - time_rates))
    print('total: ' + str(time_steady - time_ham)) 
    print( "charge generation efficiency: {val:.2f}".format(val = L0.charge_gen) + "\nratio of recombination Ex: {val:.2f}".format(val = L0.recombination_Ex) + "\nratio of recombination CT: {val:.2f}".format(val = L0.recombination_CT))
  
    L0_dict[i] = L0
        
#%%
#Plot recombination rates
save = 0

fig = plt.figure(facecolor = 'white', figsize = (8,5), dpi = 600)
ax = fig.add_axes([0.15, 0.25, 0.5, 0.7])
cmap = mpl.colormaps['viridis']  
values = np.linspace(0, 1, len(parameter_array))
colours = [cmap(i) for i in values]

for i in range(len(parameter_array)):
    ax.plot(L0_dict[i].states.dis_eh, L0_dict[i].states.En_eV, label = labels[i], color=colours[i], marker = 'x', ls = ' ', zorder = len(parameter_array)-i)
    
ax.set_xlabel('r$_{e-h}$ (Lattice Sites)', fontsize = 16)
ax.set_ylabel('Energy (eV)', fontsize = 16)
ax.tick_params(axis = 'both', direction = 'in', top = True, right = True, labelsize = 14)
ax.legend(fontsize = 14, loc = 'center left', bbox_to_anchor = [1.03, 0.5])

if save == 1:
        plt.savefig(save_path + 'reh-Energy-Vary-t0', bbox_inches = 'tight', dpi=600)

