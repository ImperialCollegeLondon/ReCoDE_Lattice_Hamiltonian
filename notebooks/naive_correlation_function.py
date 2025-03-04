"""Calculate rates of population transfer between eigenstates.

Rates are calculates using Redfield theory within the secular approximation as in
Marcus and Renger 2002 & Quantum biology revisited 2020.
"""

import numpy as np
import scipy.constants as const


def correlation_function_real_part_single(
    w: float, lambda_total: float, e_peak: float, kT: float
) -> float:
    """Calculates the real part of the correlation function of the lattice.

    Args:
        w: The energy difference between the two eigenstates of interest.
        lambda_total: The total reorganisation energy associated with each molecule.
            Defined to be the sum of the outer and inner reorganisation energies. Units
            are eV.
        e_peak: The energy of the peak of the spectral density function in eV.
            Typically 0.16 eV for organic molecules.
        kT: The typical thermal energy of the lattice in eV.

    Returns:
        float: The value of the correlation function evaluted at w.
    """
    C_w = (
        (1 + bose_einstein_occupancy_single(w, kT))
        * (
            (spectral_density_single(w, lambda_total, e_peak))
            - (spectral_density_single(-w, lambda_total, e_peak))
        )
    ) / (const.hbar / const.e)
    return 2 * np.pi * C_w


def bose_einstein_occupancy_single(w: float, kT: float) -> float:
    """Bose-Einstein occupancy function.

    Args:
        w: The energy difference between the two eigenstates of interest.
        kT: The typical thermal energy of the lattice in eV.

    Returns:
        n: The value of the Bose-Einstein occupation function for the input w and kT.
    """
    if w == 0:
        n = np.inf
    else:
        n = 1 / (np.exp(w / kT) - 1)
    return n


def spectral_density_single(
    w: float,
    lambda_total: float,
    e_peak: float,
    n: int = 15,
) -> float:
    """The spectral density function of the molecules' phonon modes.

    Args:
        w: The energy difference between the two eigenstates of interest.
        lambda_total: The total reorganisation energy associated with each molecule.
            Defined to be the sum of the outer and inner reorganisation energies.
            Units are eV.
        e_peak: The energy of the peak of the spectral density function in eV.
            Typically 0.16 eV for organic molecules.
        n: Parameter determining how rapidly the spectral denisty function decays for
            w > e_peak. Default values is 15.

    Returns:
        J_w: The value of the spectral density evaluted at w.
    """
    if w >= 0:
        a = e_peak / (3 / n) ** (1 / n)
        J_w = lambda_total * w**3 * np.exp((-w / a) ** n)
    else:
        J_w = 0
    return J_w
