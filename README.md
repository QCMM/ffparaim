ffparaim
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/QCMM/ffparaim/workflows/CI/badge.svg)](https://github.com/QCMM/ffparaim/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/QCMM/ffparaim/branch/master/graph/badge.svg)](https://codecov.io/gh/QCMM/ffparaim/branch/master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


Derivation of non-bonded non-polarizable force field parameters from Atom-in-Molecules density partition.

### Minimal setup

Required dependencies:

- openff-toolkit https://github.com/openforcefield/openff-toolkit
- openmmtools https://github.com/choderalab/openmmtools
- qc-denspart https://github.com/theochem/denspart
- qc-iodata: https://github.com/theochem/iodata
- qc-gbasis: https://github.com/theochem/gbasis

### Install (with dependencies):


    conda env create -n ffp -f environment.yml
    conda activate ffp
    pip install git+https://github.com/qcmm/ffparaim.git

(There are no releases yet.)


### Copyright

Copyright (c) 2022, QCMM


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
