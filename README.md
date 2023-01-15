# Q2DmhdFoam

`Q2DmhdFoam` is a custom OpenFOAM solver based on the built-in boussinesqBuoyantFoam solver that solves magnetohydrodynamic flows with buoyancy effects. To this end, the simplification of Sommeria and Moreau [1] is used. This allows us to run magnetohydrodynamic cases much faster than with the rigorous three-dimensional approach, thus being useful to get preliminary results before spending large amounts of time computing complex cases. The C++ solver can be found in the [`Q2DmhdFoam`](https://github.com/iraola/Q2DmhdFoam/tree/main/Q2DmhdFoam) directory.

The repository also includes some [`tutorials`](https://github.com/iraola/Q2DmhdFoam/tree/main/tutorials), [`validation`](https://github.com/iraola/Q2DmhdFoam/tree/main/validation) cases, and [`Python`](https://github.com/iraola/Q2DmhdFoam/tree/main/python) scripts to automate runs and OpenFOAM data postprocessing.

`Q2DmhdFoam` was developed as a project for Hakan Nilsson's PhD course for CFD with OpenSource Software [2]. For the full report and further information, access Chalmers University of Technology webpage [here](http://www.tfd.chalmers.se/~hani/kurser/OS_CFD/#YEAR_2020).

## References

[1] Sommeria, J., & Moreau, R. (1982). Why, how, and when, MHD turbulence becomes two-dimensional. Journal of Fluid Mechanics, 118(May), 507–518. https://doi.org/10.1017/S0022112082001177

[2] Proceedings of CFD with OpenSource Software, 2020, Edited by Nilsson H. http://dx.doi.org/10.17196/OS_CFD#YEAR_2020
