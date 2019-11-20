# Fluxo

FLUXO is a numerical tool to solve linear and nonlinear advection diffusion equations, for example the compressible Navier-Stokes equations, 
the resistive magneto-hydrodynamic equations and others. 

FLUXO is based on a novel Discontinuous Galerkin Spectral Element
Method (DGSEM), which allows for high-order of accuracy 
and fully unstructured hexahedral meshes. The volume terms are formulated in a specific subcell flux based form, which allows the construction of 
entropy stable as well as kinetic energy preserving discretisations. 
The solver is parallelized very efficiently and scales up
to hundreds of thousand cores.

The main contributers are
* from Mathematical Institute, University of Cologne: Gregor Gassner
* from Max Planck Institute for Plasma Physics, Garching: Florian Hindenlang
* from Department of Mathematics, Linköping University: Andrew Winters

FLUXO is a spin-off of the FLEXI code (https://github.com/flexi-framework/flexi) 
that is developed at the Institute for Aero- and Gasdynamics at University Stuttgart 
in the group of Claus-Dieter Munz.

If you use FLUXO, please cite one of the articles mentioned in [REFERENCES.md](REFERENCES.md)

## Installation

For installation instructions see [INSTALL.md](INSTALL.md)

## License 

FLUXO is released under the terms of the GNU General Public License v3.0. 
For the full license terms see the included license file [LICENSE.md](LICENSE.md).

## Used libraries

FLUXO uses several external libraries as well as auxilliary functions from open source projects, including:
* [HDF5](https://www.hdfgroup.org/)
* [MPI](http://www.mcs.anl.gov/research/projects/mpi/)
* [LAPACK](http://www.netlib.org/lapack/)

## List of Contributors

Numerous people have worked and are working on the code and here is the full list of contributors: [CONTRIBUTORS.md](CONTRIBUTORS.md)

We would like to thank all for their effort.

## GIT cheatsheet
You find some help with git in [doc/GITcheatsheet.md](doc/GITcheatsheet.md)
