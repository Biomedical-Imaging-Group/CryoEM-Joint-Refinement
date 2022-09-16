# Cryo-EM-Refinement
Matlab implementation of TIP 2020 paper: Joint Angular Refinement and Reconstruction for Single-Particle Cryo-EM

by: Mona Zehni, Laurene Donati, Emmanuel Soubies, Zhizhen Zhao, and Michael Unser

### Dependencies
- `GlobalBioIm`: 
https://biomedical-imaging-group.github.io/GlobalBioIm/index.html
- `ASPIRE`: http://spr.math.princeton.edu/

### Code structure
- `experiments`: Contains the main file to run.

- `functions`: A set of utility functions.

- `LinOpCryo`: Operators used in cryo-EM SPR.

- `volumes`: A set of 3D volumes in `mrc` format.

### User manual
- Run `RunMe.m` to add the paths of all folders and required libraries.
- Navigate to `./experiments/` folder and run `main.m`. You can change the dataset and the parameters of the experiments.

### Compilation of mex files
- Pre-compiled mex files of projection and adjoint of the projection operators are provided in `./LinOpCryo/functions_opt`.
- To re-compile these functions, navigate to `./LinOpCryo/mex_files/` and run `mex filename.c -o ../functions_opt/`.
- To avoid re-writing the provided pre-compiled mex files, the name of the `*.c` files slightly differ from the name of the pre-compiled mex files. Once you run any compilations, make sure the name of the functions being used in `LinOpCryo/LinOpPBTSShift.m` matches your compiled mex files.

### Citation
If you find this repositry helpful in your publications, please consider citing the following papers:
```r
@ARTICLE{zehni2020,
    author={Zehni, Mona and Donati, Laurène and Soubies, Emmanuel and Zhao, Zhizhen and Unser, Michael},
    journal={IEEE Transactions on Image Processing}, 
    title={Joint Angular Refinement and Reconstruction for Single-Particle Cryo-EM}, 
    year={2020},
    volume={29},
    number={},
    pages={6151-6163},
    doi={10.1109/TIP.2020.2984313}}
  
@article{donati2018,
    author = {Laurène Donati and Masih Nilchian and Carlos Oscar S. Sorzano and Michael Unser},
    title = {Fast multiscale reconstruction for Cryo-EM},
    journal = {Journal of Structural Biology},
    volume = {204},
    number = {3},
    pages = {543-554},
    year = {2018},
    issn = {1047-8477},
    doi = {https://doi.org/10.1016/j.jsb.2018.09.008},
    url = {https://www.sciencedirect.com/science/article/pii/S1047847718302697},
}
```
