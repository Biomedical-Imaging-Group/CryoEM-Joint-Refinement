# Cryo-EM-Refinement
Matlab implementation of TIP 2020 paper: Joint Angular Refinement and Reconstruction for Single-Particle Cryo-EM

by:Mona Zehni, Laurene Donati, Emmanuel Soubies, Zhizhen Zhao, and Michael Unser

### Pre-requisites
- GlobalBioIm: 
https://biomedical-imaging-group.github.io/GlobalBioIm/index.html
- ASPIRE: http://spr.math.princeton.edu/

### Code structure
- `experiments`: Contains the main file to run.

- `functions`: A set of utility functions.

- `LinOpCryo`: Operators used in cryo-EM SPR.

- `volumes`: A set of 3D volumes in `mrc` format.

### User manual
- Run `RunMe.m` to add the paths of all folders and required libraries.
- Run `./experiments/main.m`. You can change the dataset and the parameters of the experiments.

### Quick demo
For a quick demo, we try two small volumes of size 44 x 44 x 44. In this demo experiment, we consider high SNR and use 200 projection images. 

