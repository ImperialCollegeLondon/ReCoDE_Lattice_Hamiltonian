# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 18:12:40 2025

@author: ljh3218
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# %%
# Plot recombination rates
save = 0

def plot_energies(lattice_dict:dict, parameter_array:list[np.float32], 
                  labels:list[str], save:bool, save_path:str = " "):
    fig = plt.figure(facecolor="white", figsize=(8, 5), dpi=600)
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
    
    ax.set_xlabel("r$_{e-h}$ (Lattice Sites)", fontsize=16)
    ax.set_ylabel("Energy (eV)", fontsize=16)
    ax.tick_params(axis="both", direction="in", top=True, right=True, labelsize=14)
    ax.legend(fontsize=14, loc="center left", bbox_to_anchor=[1.03, 0.5])
    
    if save == 1:
        plt.savefig(save_path + "reh-Energy-Vary-t0", bbox_inches="tight", dpi=300)