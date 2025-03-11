# Description of the Spectral Density Function

We define the spectral density function of the basis states using the following functional, previously employed by Renger et al. $^{1}$

$J(\omega) = A\omega \mathrm{exp}\left(-\left(\frac{\omega}{\Omega}\right)^n\right)$

This function is shown in the [Introduction](01_ModelDescription.md). The parameter $\Omega$ is chosen such that the maximum of $J(\omega)$ occurs at $\hbar\omega$ = 0.16 eV, a typical energy for common intra-molecular vibrations in conjugated organic molecules, such as C-C stretching bonds. We use a high value of $n$ ($n$ = 15) to ensure that $J(\omega)$ rapidly decays to zero for $\omega > \Omega$ to reflect the fact that organic molecules do not couple strongly to phonon modes with energies greater than ~ 0.2 eV (1600 cm $^{-1}$). Lastly, $A$ is calculated so that the following normalisation condition is satisfied

$\frac{\lambda}{\hbar} = \int_{0}^{\infty} d\omega \ \frac{J(\omega)}{\omega}$

where $\lambda$ is the total input reorganisation energy of the basis states. 

## References 

1. Renger, T. & Marcus, R. A. On the relation of protein dynamics and exciton relaxation in pigment–protein complexes: An estimation of the spectral density and a theory for the calculation of optical spectra. *The Journal of Chemical Physics* **116**, 9997–10019 (2002).
