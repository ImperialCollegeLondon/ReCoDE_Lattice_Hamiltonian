"""Calculate rates of population transfer between eigenstates.

Rates are calculates using Redfield theory within the secular approximation as in
Marcus and Renger 2002 & Quantum biology revisited 2020.
"""

import numpy as np
import scipy.constants as const
import scipy.integrate as integrate
from numpy.typing import NDArray


def C_re_2D_array(
    w: NDArray[np.float32], lambda_total: float, e_peak: float, kT: float
) -> NDArray[np.float32]:
    """The correlation function of the lattice.

    Args:
        w: A numpy array containing the energy differences between all the eigenstates
            of the lattice.
        lambda_total: The total reorganisation energy associated with each molecule.
            Defined to be the sum of the outer and inner reorganisation energies. Units
            are eV.
        e_peak: The energy of the peak of the spectral density function in eV.
            Typically 0.16 eV for organic molecules.
        kT: The typical thermal energy of the lattice in eV.

    Returns:
        NDArray[np.float32]: A numpy array containing the value of the correlation
            function evaluted at every energy contained in w.
    """
    C_w = (
        (1 + nw(w, kT))
        * (
            (J_Renger(w, lambda_total, e_peak, kT=kT))
            - (J_Renger(-w, lambda_total, e_peak, kT=kT))
        )
    ) / (const.hbar / const.e)
    return 2 * np.pi * C_w


def nw(w: NDArray[np.float32], kT: float) -> NDArray[np.float32]:
    """Bose-Einstein occupancy function.

    Args:
        w: A numpy array containing the energy differences between all the eigenstates
            of the lattice.
        kT: The typical thermal energy of the lattice in eV.

    Returns:
        n: The value of the Bose-Einstein occupation function for the input w and kT.
    """
    if not hasattr(w, "__len__"):
        if w == 0:
            n = np.inf
        else:
            n = 1 / (np.exp(w / kT) - 1)
    elif hasattr(w, "__len__"):
        w[w == 0] = 1e-7
        n = 1 / (np.exp(w / kT) - 1)
    return n


def J_Renger(
    w: NDArray[np.float32],
    lambda_total: float,
    e_peak: float,
    kT: float = 0.0257,
    n: int = 15,
) -> NDArray[np.float32]:
    """The spectral density function of the molecules' phonon modes.

    Args:
        w: A numpy array containing the energy differences between all the eigenstates
            of the lattice.
        lambda_total: The total reorganisation energy associated with each molecule.
            Defined to be the sum of the outer and inner reorganisation energies.
            Units are eV.
        e_peak: The energy of the peak of the spectral density function in eV.
            Typically 0.16 eV for organic molecules.
        kT: The typical thermal energy of the lattice in eV.
        n: Parameter determining how rapidly the spectral denisty function decays for
            w > e_peak. Default values is 15.

    Returns:
        J_w: A numpy array containing the value of the spectral density evaluted at
            every energy contained in w.
    """
    a = e_peak / (3 / n) ** (1 / n)
    J_w = lambda_total * w**3 * np.exp((-w / a) ** n)
    J_w[w < 0] = 0
    return J_w


def norm_J_Renger(lambda_total: float, e_peak: float, n: int = 15) -> float:
    """Ensures that the spectral density function is correctly normalised.

    Args:
        lambda_total: The total reorganisation energy associated with each molecule.
            Defined to be the sum of the outer and inner reorganisation energies.
            Units are eV.
        e_peak: The energy of the peak of the spectral density function in eV.
            Typically 0.16 eV for organic molecules.
        n: Parameter determining how rapidly the spectral denisty function decays for
            w > e_peak.

    Returns:
        lambda_total_input: The reorganisation energy which must be input to the
            function J_Renger so that the spectral density function is corretly
            normalised.
    """

    def func(w, a, n):
        return w**2 * np.exp(-((w / a) ** n))

    w_end = 0.3
    a = e_peak / (3 / n) ** (1 / n)
    i0 = integrate.quad(lambda x: func(x, a, n), 0, w_end)
    lambda_total_input = lambda_total / i0[0]
    return lambda_total_input
