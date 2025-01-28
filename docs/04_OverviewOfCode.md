# Summary

<img src="assets/Codeoverview.png" alt="Alt Text" style="width:90%; height:auto;">

The structure of the code is illustrated in the above schematic. The core of the code is the `Lattice` class in the `lattice.py` file. This class contains all the subfunctions needed to construct the Hamiltonian (`generate_uniform` and `build_ham`), find its eigenstates ('states_from_ham') and calculate their steady-state populations ('get_rates_mat' and 'solve_steady'). The functions `states_from_ham` and `get_rates_mat` call on helper functions from the files `Recombination.py` and `Redfield.py`, respectively, though `Recombination.py` is only used if `const_recombination` is set to `False` (see [Section Three](CalculatingEigenstateProperties.md)). 
In addition to the `Lattice` class, the `lattice.py` file also contains the `Parameter` and `Site` dataclasses. The former contains the values of global variables, such as temperature, which are used in all the other functions and the latter contains the properties of each individual site in the lattice and is called by the `generate_uniform` function. 
