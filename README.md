# CR7T_Connectivity
This folder contains all the code used in our manuscript:
"Network topology of symbolic and nonsymbolic number comparison" - Benjamin Conrad, Eric Wilkey, Darren Yeo, and Gavin Price.

This was a 7 Tesla imaging project looking at functional connectivity and modular topology during a symbolic/nonsymbolic comparison task. Participants saw either a set of dots or a digit (intermixed) for each trial, and decided whether the numerosity was greater than or less than 5. Stimulus numerosities were either 2, 4, 6, 8.

Data preprocessing and beta series estimation is conducted primarily with AFNI tools, called from MATLAB, thus requires AFNI to be installed on the machine.

Connectivity matrix creation and statistical analyses are conducted in MATLAB, with algorithmic implementations of Louvain modularity/consensus-clustering/etc from the Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/).

External functions and toolboxes are included.

The manuscript is open access and has been published in Network Neuroscience at https://doi.org/10.1162/netn_a_00144.

The input data for much of this code is available at https://osf.io/sb5v2/.

<a href="https://zenodo.org/badge/latestdoi/210896655"><img src="https://zenodo.org/badge/210896655.svg" alt="DOI"></a>

