"""Calculate eigenstates of the lattice for several values of a given input parameter.

This file allows the user to see how the properties of the system's eigenstates depend
on the value of the input parameter t0. This paraemter controls the electronic coupling
between lattice sites. Once the eigenstates have been found for 5 values of t0, the code
then plots the eigenstates' energies as a function of the electron-hole separation of
the eigenstate.'
"""

# %% Import things
import os
from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from Lattice import Lattice, Parameters

now = datetime.now()
pwd = os.getcwd()

save_path = pwd + "/plots/" + now.strftime("%d%b%Y-1/")
if not os.path.exists(save_path):
    os.makedirs(save_path)

# %% Calculate states vs parameter of interest for 0 field
parameter_to_vary = "t0"
parameter_array = [1e-3, 2e-3, 3e-3, 4e-3, 5e-3]
labels = ["1 meV", "2 meV", "3 meV", "4 meV", "5 meV"]

lattice_dict = {}

spacing = 10
num_sites_coupled = 1.45
size = 4

j0 = 1.5
r0j = 0.1
min_dist = min([size, r0j * (j0 / 0.0257 - 1)])

save = 1
for i in range(len(parameter_array)):
    params = Parameters(
        temp=300,
        e_peak=0.16,
        lambda_inner=0.202,
        lambda_outer=0.048,
        j0=j0,
        r0j=r0j * spacing,
        e_singlet=1.4,
        const_recombination=False,
    )

    lattice = Lattice()
    lattice.generate_uniform(
        size=size,
        HOMO=0,
        LUMO=1.8,
        dist_sites=spacing,
        min_dist_near_neighbour=num_sites_coupled * spacing + 0.01,
        t0_homo=parameter_array[i],
        t0_lumo=parameter_array[i],
        d0=5e-3,
        r0d=0.1 * spacing,
        v_ex=0.02,
        const_recombination=False,
    )

    time_start = datetime.now()
    lattice.build_ham(
        params,
        F=[0, 0, 0],
        min_dist_near_neighbour=(num_sites_coupled * spacing) + 0.01,
        dist_cs_min=min_dist * spacing,
        disorder_site_ene=0.05,
        random_seed=42,
    )
    time_ham = datetime.now()
    print("\n" + "build ham: " + str(time_ham - time_start))
    lattice.states_from_ham(params, max_energy_diff=1.5)
    time_states = datetime.now()
    print("get states: " + str(time_states - time_ham))
    # calculate the rates of transitions between states using Redfield
    lattice.get_rate_mat(params)
    time_rates = datetime.now()
    print("get rates: " + str(time_rates - time_states))
    # solve for steady state population of the different states
    lattice.solve_steady(params)
    time_steady = datetime.now()
    print("solve steady: " + str(time_steady - time_rates))
    print("total: " + str(time_steady - time_ham))

    lattice_dict[i] = lattice

# %%
# Plot recombination rates
save = 0

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
    plt.savefig(save_path + "reh-Energy-Vary-t0", bbox_inches="tight", dpi=600)
