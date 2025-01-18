"""Calculate eigenstates of the lattice for several values of a given input parameter.

This file finds the eigenstates and occupancies of the lattice as a function of a user
specified input parameter."
"""

import numpy as np
from Lattice import Lattice, Parameters

def SweepParameter(parameter_to_vary:str, parameter_array:list[np.float32],
                   parameter_dict:dict):
       
    lattice_dict = {}
    
    spacing = 10
    num_sites_coupled = 1.45
    size = parameter_dict["size"]
    
    j0 = parameter_dict["j0"]
    r0j = parameter_dict["r0j"]
    min_dist = min([size, r0j * (j0 / 0.0257 - 1)])
    
    for i in range(len(parameter_array)):
        parameter_dict[parameter_to_vary] = parameter_array[i]
        params = Parameters(
            temp=300,
            e_peak=0.16,
            lambda_inner=0.202,
            lambda_outer=0.048,
            j0=j0,
            r0j=r0j * spacing,
            e_singlet=parameter_dict["e_singlet"],
            const_recombination=parameter_dict["const_recombination"],
        )
    
        lattice = Lattice()
        if not parameter_dict["const_recombination"]:
            lattice.generate_uniform(
                size=size,
                HOMO=0,
                LUMO=1.8,
                dist_sites=spacing,
                min_dist_near_neighbour=num_sites_coupled * spacing + 0.01,
                t0_homo=parameter_dict["t0"],
                t0_lumo=parameter_dict["t0"],
                d0=parameter_dict["d0"],
                r0d=parameter_dict["r0d"] * spacing,
                v_ex=parameter_dict["v_ex"],
                const_recombination=False,
            )
        else:
            lattice.generate_uniform(
                size=size,
                HOMO=0,
                LUMO=1.8,
                dist_sites=spacing,
                min_dist_near_neighbour=num_sites_coupled * spacing + 0.01,
                t0_homo=parameter_dict["t0"],
                t0_lumo=parameter_dict["t0"],
                d0=parameter_dict["d0"],
                r0d=parameter_dict["r0d"] * spacing,
                krec_ex=parameter_dict["krec_ex"],
                const_recombination=False,
            )
    
        lattice.build_ham(
            params,
            F=parameter_dict["F"],
            min_dist_near_neighbour=(num_sites_coupled * spacing) + 0.01,
            dist_cs_min=min_dist * spacing,
            disorder_site_ene=parameter_dict["disorder"],
            random_seed=42,
        )

        lattice.states_from_ham(params, max_energy_diff=1.5)
        # calculate the rates of transitions between states using Redfield
        lattice.get_rate_mat(params)
        # solve for steady state population of the different states
        lattice.solve_steady(params)    
        lattice_dict[i] = lattice
        
        return lattice_dict


