# Cryo-EM-Refinement
Matlab implementation of TIP 2020 paper: Joint Angular Refinement and Reconstruction for Single-Particle Cryo-EM

by:Mona Zehni, Laurene Donati, Emmanuel Soubies, Zhizhen Zhao, and Michael Unser

### Dependencies
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
    title = {Fast multiscale reconstruction for Cryo-EM},
    journal = {Journal of Structural Biology},
    volume = {204},
    number = {3},
    pages = {543-554},
    year = {2018},
    issn = {1047-8477},
    doi = {https://doi.org/10.1016/j.jsb.2018.09.008},
    url = {https://www.sciencedirect.com/science/article/pii/S1047847718302697},
    author = {Laurène Donati and Masih Nilchian and Carlos Oscar S. Sorzano and Michael Unser},
}
```
