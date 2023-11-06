# Constraint-based modeling of photorespiratory mutants using metabolomics in TMFA
This repository contains all code that has been generated for the constraint-based modeling of photorespiratory mutants in the following publication:

von Bismarck, T., Wendering, P., Perez de Souza, L. et al. Growth in fluctuating light buffers plants against photorespiratory perturbations. Nat Commun 14, 7052 (2023). [https://doi.org/10.1038/s41467-023-42648-x](https://doi.org/10.1038/s41467-023-42648-x)

## Prerequisites
* MATLAB (tested with releases 2020a and 2020b)
* [IBM CPLEX solver](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer) (v12.9)
* matTFA toolbox (please use the [enh/parallelization branch](https://github.com/pwendering/matTFA/tree/enh/parallelization) of my GitHub fork)
* [COBRA toolbox](https://opencobra.github.io/cobratoolbox/stable/installation.html) (v3.0)
* the code was tested and run under both Ubuntu 20.04.3 LTS and Windows 10

## Setup
1) Clone this repository

```git clone https://github.com/pwendering/model-prm```

2) Clone the enh/parallelization branch of the matTFA fork

```git clone -b enh/parallelization https://github.com/pwendering/matTFA```

## Steps to replicate the results
The script at [Code/model_prm.m](https://github.com/pwendering/model-prm/blob/master/Code/model_prm.m) performs the complete analysis.

To run any of the scripts, please change your working directory to `model-prm/Code`.

Please note that you may have to adjust the following variables before you run the model_prm.m script:

#### COBRA_PATH
Path to the top directory of the COBRA toolbox.

#### CPLEX_PATH
Path to the MATLAB functions that come with the CPLEX solver.

#### MATTFA_PATH
Path to the top directory of the matTFA toolbox.

#### N_CPU
Set to one if you do not want to run the code parallelized or set to the desired number of cpus.
This affects both variability analysis and sampling.

#### N_SAMPLES
Set the number of flux samples.


