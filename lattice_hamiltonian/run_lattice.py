"""Calculate eigenstates of the lattice for several values of a given input parameter.

This file finds the eigenstates of the lattice and their occupancies as a function of
a user specified input parameter. It returns a dictionary containing the solutions."
"""

import numpy as np

from lattice_hamiltonian.lattice import Lattice, Parameters


def check_valid_parameters(parameter_dict: dict) -> None:
    """A function which checks that the simulation parameters are physically reasonable.

    Args:
        parameter_dict: A dictionary containing the values of the parameters which are
            to be used in the simulation.
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
        "e_singlet": [0.5, 3],
        "F": [0, 20e-3],
        "temp": [10, 400],
        "e_peak": [0.1, 0.2],
        "lambda_inner": [0.01, 0.5],
        "lambda_outer": [0.01, 0.5],
        "HOMO": [0, 1],
        "LUMO": [0.5, 3],
        "spacing": [1, 15],
        "num_sites_coupled": [1, 3],
    }

    for key in parameter_dict:
        if key == "const_recombination":
            continue
        par = parameter_dict[key]
        bounds = valid_parameters[key]
        if not bounds[0] <= par <= bounds[1]:
            raise ValueError(
                f"{par} is not a physical value of the parameter {key}."
                f"Please choose a value in the range {bounds[0]} to {bounds[1]}."
            )

    if (parameter_dict["HOMO"] >= parameter_dict["LUMO"]) or *parameter_dict["HOMO"] >=  parameter_dict["e_singlet"]):
        raise ValueError(
            "The HOMO energy must be less than e_singlet and the LUMO energy."
        )
    if parameter_dict["e_singlet"] >= parameter_dict["LUMO"]:
        raise ValueError("e_singlet must be less than the LUMO energy.")


def check_parameter_dict(parameter_dict: dict) -> dict:
    """A function which checks that all the necessary parameters are specified.

    This function checks the parameter_dict provided by the user to see which input
    parameters have not been determined. Any undetermined parameters will be set to the
    default values contained in the default_parameters dictionary.

    Args:
        parameter_dict: A dictionary containing the values of the parameters which are
            to be used in the simulation.
    """
    default_parameters = {
        "size": 6,
        "j0": 1.5,
        "r0j": 0.1,
        "v_ex": 0.02,
        "krec_ex": 1e9,
        "t0": 5e-3,
        "d0": 5e-3,
        "r0d": 0.1,
        "disorder_site_ene": 0.05,
        "e_singlet": 1.4,
        "F": 0,
        "temp": 300,
        "e_peak": 0.16,
        "lambda_inner": 0.202,
        "lambda_outer": 0.048,
        "HOMO": 0,
        "LUMO": 1.8,
        "spacing": 10,
        "num_sites_coupled": 1.45,
        "const_recombination": 0,
    }
    for key in default_parameters.keys():
        if key not in parameter_dict.keys():
            parameter_dict[key] = default_parameters[key]
    return parameter_dict


def single_parameter(parameter_dict: dict, max_energy_diff: float = 1.5):
    """A function to solve the Hamiltonian for the inputs given in parameter_dict.

    Args:
        parameter_dict: A dictionary containing the values of the parameters to be set
            by the user. Any input parameters not specified by the user will be set
            to the default values provided in the check_parameter_dict function.
        max_energy_diff: Eigenstates with an energy greater than max_energy_diff above
                the lowest energy state will be discarded. Set to be 1.5 eV by default.

    Returns:
        L0: A instance of the lattice class for which the eigenstates and their
            occupations have been calculated based upon the parameter values given in
            parameter_dict.
    """
    parameter_dict = check_parameter_dict(parameter_dict)
    check_valid_parameters(parameter_dict)
    spacing = parameter_dict["spacing"]
    num_sites_coupled = parameter_dict["num_sites_coupled"]
    size = parameter_dict["size"]

    j0 = parameter_dict["j0"]
    r0j = parameter_dict["r0j"]

    params = Parameters(
        temp=parameter_dict["temp"],
        e_peak=parameter_dict["e_peak"],
        lambda_inner=parameter_dict["lambda_inner"],
        lambda_outer=parameter_dict["lambda_outer"],
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
        HOMO=parameter_dict["HOMO"],
        LUMO=parameter_dict["LUMO"],
        dist_sites=spacing,
        min_dist_near_neighbour=(num_sites_coupled * spacing) + 0.01,
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
    lattice.states_from_ham(params, max_energy_diff=max_energy_diff)
    # calculate the rates of transitions between states
    lattice.get_rate_mat(params)
    # solve for steady state population of the different states
    lattice.solve_steady(params)

    return lattice


def sweep_parameter(
    parameter_to_vary: str,
    parameter_array: list[np.float32],
    parameter_dict: dict,
    max_energy_diff: float = 1.5,
) -> dict:
    """A function to solve the Hamiltonian for several values of one parameter.

    Args:
        parameter_to_vary: A string defining the name of the parameter being varied. The
            string should match the name of the argument which is relevant to the
            parameter being varied e.g., if varying the exciton energy, one would
            have parameter_to_vary = "e_singlet".
        parameter_array: A list containing the sample values for the parameter being
            varied.
        parameter_dict: A dictionary containing the values of the parameters being held
            constant. Any input parameters not specified by the user will be set to the
            default values provided in the check_parameter_dict function.
        max_energy_diff: Eigenstates with an energy greater than max_energy_diff above
            the lowest energy state will be discarded. Set to be 1.5 eV by default.

    Returns:
        L0_dict: A dictionary containing the solutions for each value of the parameter
            being varied. The keys of the dictionary are numerical and start from 0
            i.e., the solution for the first value of the parameter being varied can be
            accessed by calling L0_dict[0].
    """
    lattice_dict = {}

    parameter_dict = check_parameter_dict(parameter_dict)

    for i in range(len(parameter_array)):
        parameter_dict[parameter_to_vary] = parameter_array[i]
        lattice_dict[i] = single_parameter(parameter_dict, max_energy_diff)
    return lattice_dict
