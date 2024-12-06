In this ReCoDe exemplar, we construct a Hamiltonian which describes the photoexcited state of an organic molecular crystal, such as rubrene. We then find the eigenstates of this Hamiltonian and their occupation under steady-state conditions. From this, we can calculate macroscopic properties of the crystal, such as its absorption spectra. 

To model the system, we consider a 2D lattice in which each site represents a molecule. When a molecule is photoexcited, one electron will move from the highest occupied molecular orbital (HOMO) to the lowest unoccupied molecular orbital (LUMO). Due to the low dielectric constant of organic materials, the resulting electron-hole pair will form a tightly bound state called a Frenkel exciton. Thus, our model contains two types of basis states: excitons, when the electron and hole occupy the same lattice site and electron-hole pairs, when the electron and hole are located on different lattice sites. The energy of excitonic states, $E_{\mathrm{ex}}$, is given by 

$E_{\mathrm{ex}}=E_{\mathrm{LUMO}}-E_{\mathrm{HOMO}}-E_{\mathrm{B}}$

In which $E_{\mathrm{LUMO}}$ and $E_{\mathrm{HOMO}}$ and the energy levels of the LUMO and HOMO, respectively and $E_{\mathrm{B}}$ is the exciton binding energy. The energy of electron-hole pairs is given by

$E(|\vec{r}|)=E_{\mathrm{LUMO}}-E_{\mathrm{HOMO}}-\frac{J_{0}}{1-|\vec{r}|/r_{\mathrm{0,j}}}$

where $\vec{r}$ is the vector connecting the electron-hole pair, and the second term represents the electrostatic attraction between them, which we have modelled using the Mataga potential. The model can also include the effects of uniform electric field, $\vec{F}$, by altering the energies of electron-hole basis elements by an amount $\Delta E= \vec{r}∙\vec{F}$.

These energies form the diagonal elements of the Hamiltonian. To simulate the disorder commonly found in the site energies of organic crystals, the model includes a Gaussian disorder term, which is added to the diagonal elements of the Hamiltonian. 

In the off diagonal elements, we consider interactions between the basis states. Each exciton is weakly coupled to other excitons in the lattice, and we model this using a dipole-dipole interaction which takes the functional form

$d(|\vec{r}|)=\frac{d_{0}}{(1+(|\vec{r}|-a)/r_{\mathrm{0,d}})^{3}}$ 

where $d_{0}$ and $r_{\mathrm{0,d}}$ are the two parameters which characterise the strength of the dipole-dipole interaction and $a$ is the spacing between lattice sites. We note that this form of the dipole-dipole interaction is equivalent to assuming that the transition dipole moments of the basis states are all oriented parallel to one another. 

In addition to this coupling between excitonic states, it is also possible for the electron (hole) to hop to a neighbouring lattice site, leaving the initial molecule with a positive (negative) charge and giving one of its neighbours a negative (positive) charge. The size of the coupling describing this interaction depends on the wavefunction overlap between the molecules involved, which typically decays exponentially with distance. Thus, we assume that this type of interaction can only take place between directly adjacent lattice sites and it is parameterised by the electronic coupling $t_{\mathrm{0,LUMO}}$ ($t_{\mathrm{0,HOMO}}$). 

Putting these parts together, we can write the Hamiltonian describing the electronic states as

$\hat{H}_{\mathrm{el}} = \sum \limits _{k} E_{\mathrm{ex}}|k,k⟩⟨k,k| + \sum \limits _{i \neq j}(E(|\vec{r}|)-q\vec{r}∙\vec{F})|i,j⟩⟨i,j|
+ \sum \limits _{k \neq k'}d(|\vec{r}|)|k,k⟩⟨k',k'| + \sum \limits _{i \neq j, j'} t_{\mathrm{0,HOMO}}|i,j⟩⟨i,j'| + \sum \limits _{i, i' \neq j}t_{\mathrm{0,LUMO}}|i,j⟩⟨i',j|$

where basis states indexed with k are excitonic in character, while those indexed using i and j describe electron-hole pairs. The i index refers to the lattice site on which the electron is located and the j index to that on which the hole is localised. We then diagonalise this Hamiltonian to find the eigenstates and energy levels of the system.