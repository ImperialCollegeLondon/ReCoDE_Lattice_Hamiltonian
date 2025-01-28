"""Calculate eigenstates of the lattice for several values of a given input parameter.

This file finds the eigenstates of the lattice and their occupancies as a function of 
a user specified input parameter. It returns a dictionary containing the solutions."
"""

import numpy as np

from Lattice import Lattice, Parameters


def sweep_parameter(parameter_to_vary:str, parameter_array:list[np.float32],
                   parameter_dict:dict):
    """A function to solve the Hamiltonian for several values of one parameter.

    Args:
        parameter_to_vary: A string defining the name of the parameter being varied. The
            string should match the name of the argument which is relevant to the 
            parameter being varied e.g., if varying the exciton energy, one would 
            have parameter_to_vary = "e_singlet". The parameters which can be varied 
            are: "size", "j0", "r0j", "e_singlet", "const_recombination", "krec_ex",
            "v_ex", "d0", "r0d", "t0", "disorder_site_ene", "F".
        parameter_array: A list containing the sample values for the parameter being 
            varied.
        parameter_dict: A dictionary containing the values of the parameters being held
            constant. This should contain values for each of the parameters: size, j0, 
            r0j, e_singlet, const_recombination, krec_ex, v_ex, d0, r0d, t0, 
            disorder_site_ene and F except for whichever of these is being varied. 

    Returns:
        L0_dict: A dictionary containing the solutions for each value of the parameter 
            being varied. The keys of the dictionary are numerical and start from 0 
            i.e., the solution for the first value of the parameter being varied can be 
            accessed by calling L0_dict[0]. 
    """
    lattice_dict = {}
    
    spacing = 10
    num_sites_coupled = 1.45
    size = parameter_dict["size"]
    
    j0 = parameter_dict["j0"]
    r0j = parameter_dict["r0j"]
    
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
            v_ex = parameter_dict["v_ex"]
            krec_ex = 0
        else:
            v_ex = 0
            krec_ex = parameter_dict["krec_ex"]
            
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
            v_ex=v_ex,
            krec_ex=krec_ex,
            const_recombination=parameter_dict["const_recombination"],
        )

    
        lattice.build_ham(
            params,
            F=parameter_dict["F"],
            min_dist_near_neighbour=(num_sites_coupled * spacing) + 0.01,
            disorder_site_ene=parameter_dict["disorder_site_ene"],
            random_seed=42,
        )

        # calculate the eigenstates of the lattice
        lattice.states_from_ham(params, max_energy_diff=1.5)
        # calculate the rates of transitions between states 
        lattice.get_rate_mat(params)
        # solve for steady state population of the different states
        lattice.solve_steady(params)    
        lattice_dict[i] = lattice
        
    return lattice_dict


