# Using ```seaborn``` to make heatmaps

Here, we will look at the ```plot_state_distribution``` in the ```lattice_plots.py``` file. The purpose of this function is to plot how the electron and hole probability densities are distributed across the lattice. This is useful as it allows you to visualise how delocalised the eigenstates are across the lattice. To see how this function can be used, work through the ```WorkedExample``` jupyter notebook and we show some examples of its output in Figure One, below. 

A heatmap provides a convenient way to represent the electron and hole probability distributions for a 2D lattice. Although it is possible to make a heatmap using ```matplotlib```, the ```seaborn``` package provides an alternative method and the plots generated tend to look nicer for less effort than those obtained using ```matplotlib```. Consequently, we will discuss the ```seaborn``` package here. 

We start by going over the relevant section of the ```plot_state_distribution`` function and annotating what each part of the code is doing:

<img src="assets/seaborn_code.png" alt="Alt Text" style="width:100%; height:auto;">



