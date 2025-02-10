"""Calculate eigenstates of the lattice for several values of a given input parameter.

This file finds the eigenstates of the lattice and their occupancies as a function of
a user specified input parameter. It returns a dictionary containing the solutions."
"""

import numpy as np

from lattice_hamiltonian.lattice import Lattice, Parameters


def sweep_parameter(
    parameter_to_vary: str, parameter_array: list[np.float32], parameter_dict: dict
):
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

    valid_parameters = {
        "size": [2, 10],
        "j0": [1, 3],
        "r0j": [0.01, 1],
        "v_ex": [1e-3, 0.1],
        "krec_ex": [1e8, 1e11],
        "t0": [1e-3, 50e-3],
        "d0": [1e-3, 50e-3],
        "r0d": [0.01, 1],
        "disorder_site_ene": [0, 100e-3],
        "e_singlet": [0.5, 1.75],
        "F": [0, 20e-3],
    }

    bounds = valid_parameters[parameter_to_vary]
    for par in parameter_array:
        if par > bounds[1] or par < bounds[0]:
            raise ValueError(
                f"{par} is not a physical value of the parameter {parameter_to_vary}.\
                \nPlease choose a value in the range {bounds[0]} to {bounds[1]}."
            )
    
    for key in parameter_dict:
        if key != 'const_recombination':
            par = parameter_dict[key]
            bounds = valid_parameters[key]
            if par > bounds[1] or par < bounds[0]:
                raise ValueError(
                    f"{par} is not a physical value of the parameter {key}.\
                    \nPlease choose a value in the range {bounds[0]} to {bounds[1]}."
                )

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
            F=[0, parameter_dict["F"], 0],
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

def single_parameter(parameter_dict: dict):
    """A function to solve the Hamiltonian for the inputs given in parameter_dict.

    Args:
        parameter_dict: A dictionary containing the values of all the parameters. This 
            should contain values for each of the parameters: size, j0,
            r0j, e_singlet, const_recombination, krec_ex, v_ex, d0, r0d, t0,
            disorder_site_ene and F.

    Returns:
        L0: A instance of the lattice class for which the eigenstates and their 
            occupations have been calculated based upon the parameter values given in 
            parameter_dict. 
    """

    valid_parameters = {
        "size": [2, 10],
        "j0": [1, 3],
        "r0j": [0.01, 1],
        "v_ex": [1e-3, 0.1],
        "krec_ex": [1e8, 1e11],
        "t0": [1e-3, 50e-3],
        "d0": [1e-3, 50e-3],
        "r0d": [0.01, 1],
        "disorder_site_ene": [0, 100e-3],
        "e_singlet": [0.5, 1.75],
        "F": [0, 20e-3],
    }

    for key in parameter_dict:
        if key != 'const_recombination':
            par = parameter_dict[key]
            bounds = valid_parameters[key]
            if par > bounds[1] or par < bounds[0]:
                raise ValueError(
                    f"{par} is not a physical value of the parameter {key}.\
                    \nPlease choose a value in the range {bounds[0]} to {bounds[1]}."
                )

    spacing = 10
    num_sites_coupled = 1.45
    size = parameter_dict["size"]

    j0 = parameter_dict["j0"]
    r0j = parameter_dict["r0j"]

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
        F=[0, parameter_dict["F"], 0],
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

    return lattice
