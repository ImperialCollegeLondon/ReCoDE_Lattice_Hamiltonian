"""Calculate eigenstates and populations for a 2D lattice of molecules.

The purpose of this module is to build the Hamiltonian which describes the (optically)
excited states of a 2D lattice of organic molecules. After building the Hamiltonian, the
code finds its eigenstates and the corresponding eigenvalues. Additionally, the
populations of the eigenstates are calculated under conditions of constant illumination.
"""

import math
from copy import deepcopy
from dataclasses import dataclass, field
from itertools import combinations, product

import numpy as np
import pandas as pd
import scipy.constants as sp
from numpy.typing import NDArray
from scipy import linalg
from scipy.sparse import csr_array

import lattice_hamiltonian.recombination as recombination
import lattice_hamiltonian.redfield as redfield


@dataclass
class Parameters:
    """Stores global parameters of the system."""

    temp: float
    """The temperature of the lattice in Kelvin."""
    e_peak: float
    """The energy of the peak of the spectral density function in eV. Typically 0.16 eV 
    for organic molecules."""
    lambda_outer: float
    """The outer reorganisation energy of each molecule in the lattice in eV."""
    lambda_inner: float
    """The inner reorganisation energy of each molecule in the lattice in eV."""
    j0: float
    """A parameter used in determining the shape of the Mataga potenial. Units are 
    eV."""
    r0j: float
    """A parameter used in determining the shape of the Mataga potenial. Units are 
    angstrom."""
    e_singlet: float
    """The energy of a singlet exciton in eV. This should be lower than the bandgap of 
    the material."""
    const_recombination: bool
    """A switch which determines whether recombination rates are treated as a constant 
    (True) or are calculated using generalised Marcus-Levich-Jortner (False)."""

    @property
    def kT(self) -> float:
        """Define the typical thermal energy of the system."""
        return sp.k / sp.e * self.temp


@dataclass
class Site:
    """Stores the properties of each lattice site."""

    coordinate: list[float]
    """A list of the (x,y,z) coordinates of the lattice site. Units are angstrom. """
    HOMO: float
    """Energy of the lattice site's highest occupied molecular orbital in eV."""
    LUMO: float
    """Energy of the lattice site's lowest unoccupied molecular orbital in eV."""
    id: int
    """The ID of the lattice site."""
    const_recombination: bool
    """A switch which determines whether recombination rates are treated as a constant 
    (True) or are calculated using generalised Marcus-Levich-Jortner (False)."""
    transition_dipole_ex: float
    """Number between zero and one which controls the rate at which excitons are 
    generated in the lattice."""
    nearest_neighbour: list = field(default_factory=list)
    """List containing the IDs of all the lattice sites to include when calculating the 
    size of the dipole-dipole coupling."""
    LUMO_coupling: list = field(default_factory=list)
    """List containing the electronic coupling between the LUMO of the lattice site and 
    the LUMOs of all its adjacent neighbours."""
    HOMO_coupling: list = field(default_factory=list)
    """List containing the electronic coupling between the HOMO of the lattice site and 
    the HOMOs of all its adjacent neighbours."""
    dipole_coupling: list = field(default_factory=list)
    """List containing the dipole coupling between the current lattice site and all 
    the lattice sites included in the list Nearest_Neighbour."""
    krec_ex: float = field(default=0, kw_only=True)
    """The decay rate of excitonic states in s-1."""
    v_ex: float = field(default=0, kw_only=True)
    """The strength of the coupling between the ground state and the exciton state in 
    eV."""

    # Check they are non-zero
    def __post_init__(self) -> None:
        """Check that a value has been set for the recombination parameter.

        If clause so that the code checks the value of the recombination parameter
        which is relevant for the chosen value of const_recombination.
        """
        if self.const_recombination:
            if self.krec_ex == 0:
                raise ValueError(
                    "Recombination rate for exciton basis states cannot \
                                 be zero. Set krec_ex to a finite value."
                )
        else:
            if self.v_ex == 0:
                raise ValueError(
                    "Recombination rate for exciton basis states cannot \
                                 be zero. Set v_ex to a finite value."
                )


class Lattice:
    """Class used to constuct the system's Hamiltonian and find its properties."""

    def __init__(self) -> None:
        """Declare the attributes of the Lattice class.

        Assign attributes to the lattice class which either contain information about
        the system that needs to be passed from one function to another (e.g., to build 
        the Hamiltonian, you need to know what the properties of the sites are as these 
        make up your basis set) or which contain information which is a desired output 
        of the simulation (e.g., the states DataFrame contains information about the 
        eigenstates of the system).
        """
        self.sites = []
        """A list containing instances of the Site class which will be populated by the 
        generate_uniform function and contain information about the lattice sites"""
        self.const_recombination: bool
        """A switch which determines whether recombination rates are treated as a 
        constant (True) or are calculated using generalised Marcus-Levich-Jortner 
        (False)."""
        self.ham: csr_array
        """The Hamiltionian of the system, which is stored as a sparse matrix."""
        self.basis: list
        """The basis elements of the system. Each element of basis is a two element 
        list where the first number indicates the lattice site on which the electron is 
        located and the second number indicates the lattice site on which the hole is 
        located."""
        self.transdip_vec_ex: list[float]
        """A list containing the transition dipole of each basis state. This is used to 
        calculate the transition dipole of the eigenstates and thus the rate at which 
        they are generated."""
        self.krec_vec_ex: list[float]
        """A list containing either the recombination rate of each basis state 
        (if const_recombination = True) ot the value of the coupling stregth between 
        the basis state and the ground state (if const_recombination = False). This is 
        used to calculate the recombination rate of the eigenstates and thus the
        rate at which they decay."""
        self.dist_he: list[float]
        """A list containing the electron-hole separation of each eigenstate. This is 
        used to calculate the expectation value of the electron-hole separation of the 
        eigenstates"""
        self.is_ex: list[int]
        """A list of 1s and 0s indicating if a given basis state is excitonic or not. 
        This information is used in the calc_IPR function."""
        self.states: pd.DataFrame()
        """A dataframe containing useful information about the system's eigenstates."""
        self.rates_mat: NDArray[np.float32]
        """A 2D matrix containing the rates of transitions between the eigenstates as 
        calculated in the get_rate_mat function."""

    def generate_uniform(
        self,
        size: int,
        HOMO: float,
        LUMO: float,
        dist_sites: float,
        min_dist_near_neighbour: float,
        t0_homo: float,
        t0_lumo: float,
        d0: float,
        r0d: float,
        const_recombination: bool = True,
        krec_ex: float = 0,
        v_ex: float = 0,
    ) -> None:
        """Build the lattice and assign properties to lattice sites.

        Args:
            size: Determines the total number of lattice sites. The lattice will be a
                size x size square.
            HOMO: Energy of the material's highest occupied molecular orbital in eV.
            LUMO: Energy of the material's lowest unoccupied molecular orbital in eV.
            dist_sites: The separation between adjacent lattice sites. Units are
                angstom.
            min_dist_near_neighbour: The maximum separtion between two lattice sites
                for which excitonic coupling is considered to be non-zero. Units are
                angstrom.
            t0_homo: The electronic coupling between the HOMOs of two adjacent lattice
                sites. Units are eV.
            t0_lumo: The electronic coupling between the LUMOs of two adjacent lattice
                sites. Units are eV.
            d0: Parameter determining the strength of the dipole-dipole coupling
                between excitons. Units are eV.
            r0d: Parameter determining the strength of the dipole-dipole coupling
                between excitons. Units are eV.
            const_recombination: A switch which determines whether recombination rates
                are treated as a constant (True) or are calculated using generalised
                Marcus-Levich-Jortner (False).
            krec_ex: Rate at which excitons decay to the ground state. Defaults to 0,
                but must be non-zero if const_recombination = True. Units are inverse
                seconds.
            v_ex: Coupling of excitons to the ground state. Defaults to zero, but must
                be non-zero if const_recombination = False. Units are eV.
        """
        counter_id = 0
        self.const_recombination = const_recombination

        for x, y in product(range(size), repeat=2):
            # The 1 here is the transition dipole of the exciton.
            if const_recombination:
                self.sites.append(
                    Site(
                        np.array([x * dist_sites, y * dist_sites, 0]),
                        HOMO,
                        LUMO,
                        counter_id,
                        const_recombination,
                        1,
                        krec_ex=krec_ex,
                    )
                )
            else:
                self.sites.append(
                    Site(
                        np.array([x * dist_sites, y * dist_sites, 0]),
                        HOMO,
                        LUMO,
                        counter_id,
                        const_recombination,
                        1,
                        v_ex=v_ex,
                    )
                )
            counter_id = counter_id + 1
        for ii, jj in combinations(range(counter_id), 2):
            # Find distance from site ii to site jj.
            distance_bet_sites = math.dist(
                self.sites[ii].coordinate, self.sites[jj].coordinate
            )
            # Implement cut-off distance below which sites are considered coupled.
            if distance_bet_sites < min_dist_near_neighbour:
                self.sites[ii].nearest_neighbour.append(self.sites[jj].id)
                self.sites[jj].nearest_neighbour.append(self.sites[ii].id)
                if math.isclose(distance_bet_sites - dist_sites, 0, abs_tol=1e-3):
                    self.sites[ii].LUMO_coupling.append(t0_lumo)
                    self.sites[jj].LUMO_coupling.append(t0_lumo)
                    self.sites[jj].HOMO_coupling.append(t0_homo)
                    self.sites[ii].HOMO_coupling.append(t0_homo)
                else:
                    self.sites[ii].LUMO_coupling.append(0)
                    self.sites[jj].LUMO_coupling.append(0)
                    self.sites[jj].HOMO_coupling.append(0)
                    self.sites[ii].HOMO_coupling.append(0)
                self.sites[jj].dipole_coupling.append(
                    d0 * ((distance_bet_sites - dist_sites) / r0d + 1) ** -3
                )
                self.sites[ii].dipole_coupling.append(
                    d0 * ((distance_bet_sites - dist_sites) / r0d + 1) ** -3
                )

    def build_ham(
        self,
        params,
        F: list[float],
        min_dist_near_neighbour: float,
        disorder_site_ene: float,
        random_seed: int = 0,
    ) -> None:
        """Build the Hamiltonian for the lattice defined by self.sites.

        Args:
            params: An instance of the parameters class.
            F: The electric field applied to the lattice in the form of a list contaning
                the x, y and z components of the field i.e., F = [Fx, Fy, Fz].
            min_dist_near_neighbour: The maximum separtion between two lattice sites
                for which excitonic coupling is considered to be non-zero. Units are
                angstrom.
            disorder_site_ene: The static disorder associated with the energies of the
                basis states. Diagonal elements of the Hamiltonian will have this have
                an energy added to them which is randomly drawn from a Gaussian centred
                about zero and with a standard deviation of disorder_site_ene. Units
                are eV.
            random_seed: Optional parameter which can be used to set the value of the
                random seed. Setting this value allows one to rerun a given instance of
                the lattice. Default value is 0, meaning no user control of the random
                seed.
        """
        row = []
        col = []
        data = []
        basis = []
        const_recombination = self.const_recombination
        if random_seed != 0:
            np.random.seed(seed=random_seed)
        # Create empty lists to store the values of various quantities of interest.
        transdip_vec_ex = []
        krec_vec_ex = []
        dist_he = []
        e_singlet = params.e_singlet
        is_ex = []
        counter = 0
        for elec_pos in self.sites:
            for hol_pos in self.sites:
                distance_ele_hole = math.dist(elec_pos.coordinate, hol_pos.coordinate)
                dist_vec = elec_pos.coordinate - hol_pos.coordinate
                dist_he.append(distance_ele_hole)
                coul_j = params.j0 / (1 + distance_ele_hole / params.r0j)
                # Setting generation and decay rates for excitons.
                if distance_ele_hole == 0:
                    transdip_vec_ex.append(elec_pos.transition_dipole_ex)
                    is_ex.append(1)
                    if const_recombination:
                        krec_vec_ex.append(elec_pos.krec_ex)
                    else:
                        krec_vec_ex.append(elec_pos.v_ex)
                else:
                    transdip_vec_ex.append(0)
                    is_ex.append(0)
                    krec_vec_ex.append(0)
                # Calculate the values of the diagonal elements of the Hamiltonian.
                row.append(counter)
                col.append(counter)
                if distance_ele_hole == 0:
                    data.append(
                        e_singlet + disorder_site_ene * np.random.standard_normal(1)[0]
                    )
                else:
                    data.append(
                        elec_pos.LUMO
                        - hol_pos.HOMO
                        + disorder_site_ene * np.random.standard_normal(1)[0]
                        - coul_j
                        + np.dot(dist_vec, F)
                    )
                counter += 1
                basis.append([elec_pos.id, hol_pos.id])

        # Calculate the off diagonal elements of the Hamiltonian.
        # using combinations so will never have case where ii = jj
        for ii, jj in combinations(range(len(basis)), 2):
            if (
                basis[jj][0] in self.sites[basis[ii][0]].nearest_neighbour
                and basis[ii][1] == basis[jj][1]
            ):
                # index returns the index of the first occurance of the specified value.
                idx_el = self.sites[basis[ii][0]].nearest_neighbour.index(basis[jj][0])
                row.append(ii)
                col.append(jj)
                # hopping for electron
                data.append(self.sites[basis[ii][0]].LUMO_coupling[idx_el])
                # same again as H is a symmetric matrix
                row.append(jj)
                col.append(ii)
                data.append(self.sites[basis[ii][0]].LUMO_coupling[idx_el])
            if (
                basis[ii][0] == basis[ii][1]
                and basis[jj][0] == basis[jj][1]
                and basis[jj][0] in self.sites[basis[ii][0]].nearest_neighbour
            ):
                idx_ex = self.sites[basis[ii][1]].nearest_neighbour.index(basis[jj][1])
                row.append(ii)
                col.append(jj)
                data.append(self.sites[basis[ii][0]].dipole_coupling[idx_ex])
                row.append(jj)
                col.append(ii)
                data.append(self.sites[basis[ii][0]].dipole_coupling[idx_ex])
            if (
                basis[jj][1] in self.sites[basis[ii][1]].nearest_neighbour
                and basis[ii][0] == basis[jj][0]
            ):
                idx_h = self.sites[basis[ii][1]].nearest_neighbour.index(basis[jj][1])
                row.append(ii)
                col.append(jj)
                data.append(self.sites[basis[ii][1]].HOMO_coupling[idx_h])
                row.append(jj)
                col.append(ii)
                data.append(self.sites[basis[ii][1]].HOMO_coupling[idx_h])

        self.ham = csr_array(
            (data, (row, col)), shape=(len(self.sites) ** 2, len(self.sites) ** 2)
        )
        self.basis = basis
        self.transdip_vec_ex = transdip_vec_ex
        self.krec_vec_ex = krec_vec_ex
        self.dist_he = dist_he
        self.is_ex = is_ex

    def states_from_ham(self, params, max_energy_diff: float) -> None:
        """Find the eigenstates of the Hamiltonian and calculate their properties.

        If const_recombination = False, decay rates are found using generalised Marcus-
        Levich-Jortner following Taylor and Kassal 2018 and D'Avino et al. J. Phys.
        Chem. Lett. 2016, 7, 536-540.

        Args:
            params: An instance of the parameters class.
            max_energy_diff: An energy cut-off used to limit the number of eigenstates
                which are analysed. Eigenstates with an energy greater than that of the
                lowest energy eigenstate plus max_energy_diff are discarded.
        """
        lambda_inner, lambda_outer = (
            params.lambda_inner,
            params.lambda_outer,
        )
        evals, evecs = linalg.eigh(self.ham.toarray())
        evals = np.real(evals)
        evecs = np.real(evecs)
        # remove eigenvectors corresponding to states above the energy cut off.
        evecs = evecs[:, evals < evals.min() + max_energy_diff]
        # remove eigenvalues corresponding to states above the energy cut off.
        evals = evals[evals < evals.min() + max_energy_diff]
        # Declare a empty lists which will be populated by the function.
        list_evecs = []
        dis_st = []
        IPR = []
        ex_char = []
        transdip_ex = []
        occupation_prob = []
        krec_ex = []
        is_ex = np.array(self.is_ex)
        basis = np.array(self.basis)
        if not params.const_recombination:
            v_eff_ex = []
        for i in range(len(evals)):
            list_evecs.append(evecs[:, i])
            occupation_probability = evecs[:, i] ** 2
            occupation_prob.append(occupation_probability)
            dis_st.append(occupation_probability @ self.dist_he)
            transdip_ex.append(occupation_probability @ self.transdip_vec_ex)
            ex_char.append(occupation_probability @ is_ex)
            if params.const_recombination:
                IPR.append(calc_IPR(evecs[:, i], basis, is_ex))
                krec_ex.append(occupation_probability @ self.krec_vec_ex)
            else:
                IPR.append(calc_IPR(evecs[:, i], basis, is_ex))
                effective_coupling_ex = evecs[:, i] @ self.krec_vec_ex
                v_eff_ex.append(effective_coupling_ex)
                effective_inner_lambda_ex = lambda_inner * (1 / IPR[i][0])
                effective_outer_lambda_ex = lambda_outer * (1 / IPR[i][0])
                krec_ex.append(
                    recombination.decay_rate(
                        effective_inner_lambda_ex,
                        params.e_peak,
                        effective_outer_lambda_ex,
                        evals[i],
                        effective_coupling_ex,
                    )
                )

        states = pd.DataFrame(
            {
                "energies": evals,
                "dis_eh": dis_st,
                "IPR": IPR,
                "ex_char": ex_char,
                "transdip_ex": transdip_ex,
                "krec_ex": krec_ex,
                "occupation_probability": occupation_prob,
            }
        )
        # Make it so the states are numbered 1 to n, not 0 to n-1.
        states = states.sort_values(by=["energies"])
        states.insert(0, "state", np.arange(1, len(evals) + 1))
        if not params.const_recombination:
            states = states.assign(v_effective_ex=v_eff_ex)
        if len(self.states) == 0:
            self.states = states
        else:
            # adding newly calculated values onto the pre-existing states dataframe
            self.states = pd.concat([self.states, states], ignore_index=True)

    def get_rate_mat(self, params) -> None:
        """Calculate rates of population transfer between eigenstates.

        Rates are calculates using Redfield theory within the secular approximation as
        in Marcus and Renger 2002 & Quantum biology revisited 2020.

        Args:
            params: An instance of the parameters class.
        """
        e_peak, kT = params.e_peak, params.kT
        (
            lambda_outer,
            lambda_inner,
        ) = params.lambda_outer, params.lambda_inner
        lambda_inner = redfield.norm_spectral_density(
            lambda_outer + lambda_inner, e_peak
        )
        energies = np.array(self.states.energies)
        num_states = len(energies)
        occupation_prob = np.array(
            [self.states.occupation_probability[i] for i in range(num_states)]
        )
        rates_mat = np.zeros((num_states, num_states), dtype=np.float32)
        inds = np.tril_indices_from(rates_mat, k=-1)
        len_inds = int(((num_states) * (num_states - 1)) / 2)
        w = np.array(
            [energies[inds[0][k]] - energies[inds[1][k]] for k in range(len_inds)],
            dtype=np.float32,
        )
        C = redfield.correlation_function_real_part(w, lambda_inner, e_peak, kT)
        # Here I pick the kth element from inds[0], which gives the row index, and the
        # kth element from inds[1], which gives the column index. This corresponds to a
        # transition between the inds[0][k] and inds[1][k] eigenstates.
        gamma_Cw_real_tot = np.array(
            [
                occupation_prob[inds[0][k]] @ occupation_prob[inds[1][k]] * C[k]
                for k in range(len_inds)
            ]
        )
        rates_mat[inds] = gamma_Cw_real_tot
        rates_mat[(inds[1], inds[0])] = rates_mat[inds] * np.exp(-w / params.kT)
        self.rates_mat = rates_mat

    def solve_steady(self, params) -> None:
        """Solve Pauli's Master Equation to find the populations of the eingenstates.

        Args:
            params: An instance of the parameters class.
        """
        isteady_n = np.zeros(len(self.states))
        a = deepcopy(self.rates_mat)
        a = np.transpose(a)
        self.states["gen"] = self.states.transdip_ex
        for i in range(len(self.states)):
            a[i, i] = -sum(a[:, i])
            a[i, i] = a[i, i] - self.states.iloc[i].krec_ex
            isteady_n[i] = self.states.iloc[i].gen
        b = -isteady_n
        z = linalg.solve(a, b)
        self.states["z"] = z


def calc_IPR(
    eigvec: list[float], basis: list[float], is_ex: list[float]
) -> list[float]:
    """Calculate the inverse participation ratio of the eigenstates.

    These are calculated following method in D'Avino et al. J. Phys. Chem. Lett.
    2016, 7, 536-540.

    Args:
        eigvec: List of the eigenvectors of the system.
        basis: List of the system's basis states.
        is_ex: List identifying which elements of basis correspond to excitonic
            states.

    Returns:
        list[float]: List containing the inverse participation ratio (IPR) of the
            eigenstate. The elements of the list correspond to the IPR of the
            excitonic, electron and hole contributions to the eigenstate.
    """
    num_states = len(eigvec)
    num_sites = int(basis[-1][0] + 1)
    pop = eigvec**2

    frac_ex = sum(pop * is_ex)
    frac_ct = sum(pop * (1 - is_ex))
    IPR_ex = frac_ex**2 / sum((pop * is_ex) ** 2)

    # Loop over lattice sites.
    e_pop = np.zeros(num_sites)
    h_pop = np.zeros(num_sites)
    for k in range(num_sites):
        h_temp = 0
        e_temp = 0
        c_h_temp = 0
        c_e_temp = 0
        # Loop over basis elements.
        for j in range(num_states):
            # If electron is on site k, add to sum as need to find probability
            # electron on k for any hole position
            if basis[j][0] == k:
                # Put this line here to avoid including basis elements which
                # correspond to excitons in the sum
                if basis[j][0] != basis[j][1]:
                    e_temp += pop[j]
                c_e_temp += pop[j]
            if basis[j][1] == k:
                if basis[j][0] != basis[j][1]:
                    h_temp += pop[j]
                c_h_temp += pop[j]
        e_pop[k] = e_temp
        h_pop[k] = h_temp
        # loc_e += c_e
    IPR_e = frac_ct**2 / sum(e_pop**2)
    IPR_h = frac_ct**2 / sum(h_pop**2)
    return [IPR_ex, IPR_e, IPR_h]
