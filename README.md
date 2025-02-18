# Modelling Organic Crystals Using a Lattice Hamiltonian

## Description

In this project, we use a simplified model of an organic crystal to caluclate the systems' excited states under illumination and their populations. We then investigate how changing the input parameters of the model change the nature of the excited states.

## Learning Outcomes

- Techniques to speed up ```for``` loops
- How to run ```for``` loops in parallel uing the ```multiprocessing``` package
- How to plot heatmaps using the ```seaborn``` package

<!-- How long should they spend reading and practising using your Code.
Provide your best estimate -->

| Task       | Time    |
| ---------- | ------- |
| Reading    | 3 hours |
| Practising | 3 hours |

## Requirements

<!--
If your exemplar requires students to have a background knowledge of something
especially this is the place to mention that.

List any resources you would recommend to get the students started.

If there is an existing exemplar in the ReCoDE repositories link to that.
-->

### Academic
- Strong, undergraduate-level understanding of quantum mechanics
- Basic familiarity with solid state physics
- Intermediate-level python ability

### System
- See the [pyproject.toml](pyproject.toml) file

## Getting Started

Start by reading through sections 1-4 which describe the physics underlying this exemplar and the structure of the code. 

Once you have been through this, you can work through the next four sections. In the first of these, we walk you through how to use the code and investigate how the eigenstates of the system change when we change the strength of the coupling between lattice sites. The next three sections focus in detail on short extracts from the code which are relevant to the learning outcomes of this exemplar. 

## Project Structure

```log
.
├── docs
│   ├── 01_ModelDescription.md
│   ├── 02_FindingSteadyStatePopulations.md
│   ├── 03_CalculatingEigenstateProperties.md
│   ├── 04_OverviewOdCode.md
│   ├── 06_SpeedingUp-for-Loops.md
│   ├── 08_Using-seaborn.md
│   └── ...
├── lattice_hamiltonian
|   ├── lattice.py
|   ├── lattice_plots.py
|   ├── recombination.py
|   ├── redfield.py
│   └── run_lattice.py
├── notebooks
|   ├── plots
|   ├── 05_WorkedExample.ipynb
│   └── 07_Using-multiprocessing.ipynb
└── ...
```

<!-- Change this to your License. Make sure you have added the file on GitHub -->

## License

This project is licensed under the [BSD-3-Clause license](LICENSE.md)
