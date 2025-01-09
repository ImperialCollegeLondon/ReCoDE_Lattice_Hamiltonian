# Calculating Properties of the Eigenstates

The eigenstates of the system are found by solving the eigenvalue problem

$\hat{H}\psi = E\psi$

Each excited state, $\alpha$, is characterised by an energy, $E_{\alpha}$, and eigenvector

$\psi^{(\alpha)}= \sum_{\mathrm{k}}c_{\mathrm{kk}}^{(\alpha)}|k,k\rangle + \sum_{\mathrm{i\neq j}}c_{\mathrm{ij}}^{(\alpha)}|i,j\rangle$

The electron-hole separation of the state is calculated as

```math
r^{(\alpha)}= \langle\psi^{(\alpha)}||\vec{r}_{i} - \vec{r}_{j}||\psi^{(\alpha)}\rangle=\sum_{\mathrm{i\neq j}}|c_{\mathrm{ij}}^{(\alpha)}|^{2}|\vec{r}_{i} - \vec{r}_{j}|
```
where $\vec{r}_{i}$ is the vector describing the position of the basis state $i$ in the lattice. Other eigenstate properties, such as the transition dipole moment, can be calculated as expectation values in a similar manner. 

The excitonic and charge transfer character of each eigenstate state can be evaluated using the contributions from excitonic $|k,k\rangle$ and charge transfer states $|i,j\rangle$ basis elements, respectively, using the expressions

$\rho_{\mathrm{ex}}=\sum_{\mathrm{k}}|c_{\mathrm{kk}}^{(\alpha)}|^{2}$ 

$\rho_{\mathrm{CT}}=\sum_{\mathrm{i\neq j}}|c_{\mathrm{ij}}^{(\alpha)}|^{2}$

The calculated excited states will be delocalised over the whole basis, where the degree of delocalization is given by the inverse participation ratio, which is generally defined as

$IPR^{(\alpha)} = \sum_{\mathrm{i}}|c_{\mathrm{i}}^{(\alpha)}|^{-4}$

For a given eigenstate, the inverse participation ratio of the exciton, electron and hole is calculated by first defining a reduced wavefunction for each of these species. For example, in the case of the electron

$\psi_{\mathrm{el}}^{(\alpha)} = \frac{1}{\sqrt{\rho_{\mathrm{CT}}}}\sum_{i}\sum_{j,i\neq j}c_{\mathrm{ij}}^{(\alpha)}|i,j\rangle =  \frac{1}{\sqrt{\rho_{\mathrm{CT}}}}\sum_{i}\tilde{c}_{\mathrm{i}}^{(\alpha)}|i\rangle $

where the prefactor ensures that the new wavefunction is properly normalised. This wavefunction is then substituted into the definition of the inverse participation ratio to get 

$IPR_{\mathrm{el}}^{(\alpha)} = \frac{\rho_{\mathrm{CT}}^{2}} {\sum_{i}|\tilde{c}_{\mathrm{i}}^{(\alpha)}|^{4}}$

The same reasoning can be used to find the IPR of the exciton and hole contributions to the eigenstate's wavefunction i.e., 

$IPR_{\mathrm{ex}}^{(\alpha)} = \frac{\rho_{\mathrm{ex}}^{2}} {\sum_{kk}|{c}_{\mathrm{kk}}^{(\alpha)}|^{4}}$

$IPR_{\mathrm{hole}}^{(\alpha)} = \frac{\rho_{\mathrm{CT}}^{2}} {\sum_{j}|\tilde{c}_{\mathrm{j}}^{(\alpha)}|^{4}}$
