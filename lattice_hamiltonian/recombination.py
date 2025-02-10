"""Calculate the decay rates of eigenstates using generalised Marcus-Levich-Jornter.

If const_recombination is set to False, the decay_rate function of this module is
called by the states_from_ham function of the Lattice class in Lattice.py to calculate
the decay rate of the eigenstates.
"""

from itertools import product

import numpy as np
import scipy.constants as const


def laguerre(alpha: int, n: int, x: float) -> float:
    """Calculate values of generalised Laguerre polynomials using a recursion relation.

    Generalised Laguerre polynomials are solution of the differential equation:
        xy'' + (alpha + 1 - x)y' + ny = 0

    See https://en.wikipedia.org/wiki/Laguerre_polynomials for more information about
    Laguerre polynomials.

    Args:
        alpha: Parameter determining which Laguerre polynomial is needed.
        n: Parameter determining which Laguerre polynomial is needed.
        x: Huang-Rhys factor of the high frequency phonon mode assocaited with each
            lattice site.

    Returns:
        L: The value of the generalised Laguerre polynomial corresponding to alpha and n
            and sampled at x.
    """
    L_0 = 1
    L_1 = 1 + alpha - x
    if n == 0:
        return L_0
    elif n == 1:
        return L_1
    else:
        L = (1 / (n)) * (
            (2 * n - 1 + alpha - x) * laguerre(alpha, n - 1, x)
            - (n + alpha - 1) * laguerre(alpha, n - 2, x)
        )
        return L


def calc_decay_rate_single_vibronic_mode(
    n: int,
    m: int,
    lambda_inner: float,
    e_peak: float,
    lambda_outer: float,
    w: float,
    kT: float,
) -> float:
    """Calculates decay rate from excited to ground state for given values of m and n.

    More specifically, this function calculates the decay rate from an excited state in
    vibrational level m to the ground state in vibrational level n. It uses the
    formula found in https://journals.aps.org/prx/pdf/10.1103/PhysRevX.8.031055

    Args:
        n: Vibrational energy level of the excited state.
        m: Vibrational energy level of the ground state.
        lambda_inner: The inner reorganisation energy of each molecule in the lattice.
            Units are eV.
        e_peak: The energy of the peak of the spectral density function in eV.
            Typically 0.16 eV for organic molecules.
        lambda_outer: The outer reorganisation energy of each molecule in the lattice.
            Units are eV.
        w: The energy of the excited state, negelecting vibronic contributions.
        kT: The avaialable thermal energy in eV.

    Returns:
        Float: The decay rate from an excited state in vibrational level m to the
            ground state in vibrational level n.
    """
    # Huang-Rhys factor
    S = lambda_inner / e_peak
    # Vibronic integral for any n,m
    prefactor = np.exp(-S) * S ** (n - m) * np.math.factorial(m) / np.math.factorial(n)
    lag = laguerre(n - m, m, S)
    activation = np.exp(
        -((w - (n - m) * e_peak - lambda_outer) ** 2) / (4 * lambda_outer * kT)
    )
    # prefactor is for normalisation
    thermal_pop = (1 - np.exp(-e_peak / kT)) * np.exp(-(m * e_peak) / kT)
    return prefactor * lag**2 * activation * thermal_pop


def calc_decay_rate_all_vibronic_modes(
    lambda_inner: float,
    e_peak: float,
    lambda_outer: float,
    w: float,
    kT: float,
    N: int = 20,
    M: int = 6,
) -> float:
    """Calculates the total decay rate of an excited state to the ground state.

    N and M have been set to default values which are reasonable for kT = 0.0257 eV. If
    interested in a lower temperature, N and M could be lowered.

    Args:
        lambda_inner: The inner reorganisation energy of each molecule in the lattice.
            Units are eV.
        e_peak: The energy of the peak of the spectral density function in eV.
            Typically 0.16 eV for organic molecules.
        lambda_outer: The outer reorganisation energy of each molecule in the lattice.
            Units are eV.
        w: The energy of the excited state, negelecting vibronic contributions.
        kT: The avaialable thermal energy in eV.
        N: The total number of vibronic modes to consider for the excited state. Default
            value is 20.
        M: The total number of vibronic modes to consider for the ground state. Defualt
            value is 6.

    Returns:
        Float: The total decay rate of an excited state to the ground state, summing
            over multiple vibronic modes.
    """
    vib_overlap_total = 0
    for n, m in product(range(N), range(M), repeat=1):
        vib_overlap_total += calc_decay_rate_single_vibronic_mode(
            n, m, lambda_inner, e_peak, lambda_outer, w, kT
        )
    # Factor of 1/e to convert eV to joules
    return (
        (1 / const.e) * (1 / np.sqrt(4 * np.pi * lambda_outer * kT)) * vib_overlap_total
    )


def decay_rate(
    lambda_inner: float,
    e_peak: float,
    lambda_outer: float,
    w: float,
    v: float,
    kT: float,
) -> float:
    """Calculates the decay rate of an excited state to the ground state.

    Args:
        lambda_inner: The inner reorganisation energy of each molecule in the lattice.
            Units are eV.
        e_peak: The energy of the peak of the spectral density function in eV.
            Typically 0.16 eV for organic molecules.
        lambda_outer: The outer reorganisation energy of each molecule in the lattice.
            Units are eV.
        w: The energy of the excited state, negelecting vibronic contributions.
        v: The The strength of the coupling between the ground state and the exciton
            state in eV.
        kT: The avaialable thermal energy in eV.

    Returns:
        Float:The rate at which the eigenstate with energy w and coupling v decays to
            the ground state.
    """
    # Convert v into joules
    v *= const.e
    vib_overlap_0 = calc_decay_rate_all_vibronic_modes(
        lambda_inner, e_peak, lambda_outer, w, kT
    )
    return (2 * np.pi * v**2 * vib_overlap_0) / const.hbar
