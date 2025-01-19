"""Various plots to visualise the eigenstates of the lattice."""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def plot_energies(lattice_dict:dict, parameter_array:list[np.float32], 
                  labels:list[str], input_parameter:str, save:bool, 
                  save_path:str = " "):
    """Plots the energies of the eigenstates versus the electron-hole separation.

    Args:
        lattice_dict: A dictionary of solutions as produced by the sweep_parameter 
            function in RunLattice.py
        parameter_array: A list containing the sample values for the parameter being 
            varied.
        labels: A list containing strings to be used in the legend of the plot.
        input_parameter: A string indicating which parameter was being varied. This will
        be used as the title of the legend. 
        save: A boolean indicating if you want to save the plot.
        save_path: A string indicating the file in which you want to save the plot. Only
            required if save = True. 
    """
    fig = plt.figure(facecolor="white", figsize=(6,4), dpi=600)
    ax = fig.add_axes([0.15, 0.25, 0.5, 0.7])
    cmap = mpl.colormaps["viridis"]
    values = np.linspace(0, 1, len(parameter_array))
    colours = [cmap(i) for i in values]
    
    for i in range(len(parameter_array)):
        ax.plot(
            lattice_dict[i].states.dis_eh,
            lattice_dict[i].states.energies,
            label=labels[i],
            color=colours[i],
            marker="x",
            ls=" ",
            zorder=len(parameter_array) - i,
        )
    
    ax.set_xlabel("r$_{e-h}$ ($\AA$)", fontsize=16)
    ax.set_ylabel("Energy (eV)", fontsize=16)
    ax.tick_params(axis="both", direction="in", top=True, right=True, labelsize=14)
    ax.legend(fontsize=14, loc="center left", bbox_to_anchor=[1.03, 0.5],
              title = input_parameter, title_fontsize = 14)
    
    if save:
        plt.savefig(save_path + f"Energy-vs-reh-Vary-{input_parameter}.png", 
                    bbox_inches="tight", dpi=300)
        
def plot_state_distribution(L0,state_num:int,save:bool,sort_by_energy:bool=True, 
                            save_path:str = " "):
    """Plots the distribution of the eigenstates over the lattice.
    
    This function creates a plot showing how the electron and hole probability density 
    of a given eigenstate are distributed over the lattice. 

    Args:
        L0: An instance of the lattice class for which the eigenstates of the 
            Hamiltonian have been found (i.e., you have run the states_from_ham 
            function).
        state_num: The number of the state for which you want to plot the distribution.
            By default, the states are ranked by their energy so state_num = 0 would be 
            the lowest energy eigenstate. 
        save: A boolean indicating if you want to save the plot.
        sort_by_energy: A boolean indicating if you want the states to be sorted by 
            their energy (True) or by their electron-hole separation (False). The 
            default is True.  
        save_path: A string indicating the file in which you want to save the plot. Only
            required if save = True. 
    """
    States = L0.states
    if not sort_by_energy:
        #Sort states in order of reh from low to high
        States = States.sort_values(by=['dis_eh']).reset_index(drop=True)
    elif sort_by_energy:
        #Sort states in order of energy from low to high
        States = States.sort_values(by=['energies']).reset_index(drop=True)
    num_states = len(L0.basis)
    
    electron_density = np.zeros(int(num_states**0.5))
    hole_density = np.zeros(int(num_states**0.5))
    
    for i in range(num_states):
        el_loc = L0.basis[i][1]
        hole_loc = L0.basis[i][0]
        electron_density[el_loc] += States['occupation_probability'].iloc[state_num][i]
        hole_density[hole_loc] += States['occupation_probability'].iloc[state_num][i]
    electron_density = np.resize(electron_density, 
                                 (int(num_states**0.25), int(num_states**0.25)))
    hole_density = np.resize(hole_density, 
                             (int(num_states**0.25), int(num_states**0.25)))

    fig, axes = plt.subplots(1,2, facecolor = 'white', figsize = (10,10), dpi = 600)
    sns.heatmap(electron_density.transpose(), cmap = 'Blues', robust = True, 
                annot = False, cbar = False, ax = axes[0], square = True,
                linecolor = 'black', linewidths = 1)
    sns.heatmap(hole_density.transpose(), cmap = 'Reds', robust = True, 
                annot = False, cbar = False, ax = axes[1], square = True, 
                linecolor = 'black', linewidths = 1)
    
    for ax in axes:
        ax.set_xticks([], [])
        ax.set_yticks([], [])  
    
    axes[0].set_xlabel('Electron', fontsize = 20)
    axes[1].set_xlabel('Hole', fontsize = 20)
    
    if save:
        name = "state_distribution_{}_eV.png".format(States['energies'].iloc[state_num])
        plt.savefig(save_path + name, 
                    bbox_inches="tight", dpi=300)