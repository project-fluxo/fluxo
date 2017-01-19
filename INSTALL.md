## FLUXO installation procedure

## Prerequisites

The required packages for the Ubuntu Linux distributions are listed in table.
Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git
    

| Package          | Ubuntu 14.04    | Ubuntu 16.04    |
|:----------------:|:---------------:|:---------------:|
| git              | x               |      x          |
| cmake            | x               |      x          |
| cmake-curses-gui | x               |      x          |
| liblapack3       | x               |      x          |
| liplapack-dev    | x               |      x          |
| gfortran         | x               |      x          |
| g++              | x               |      x          |
|  mpi-default-dev | x               |      x          |
| zlib1g-dev       | -               |     x           |

Table: Required debian packages under Ubuntu.


## Compiling the code

* Open a terminal
* Change into the FLUXO directory
* Create a new subdirectory and use CMake to configure and compile the code
```
mkdir build; cd build
cmake ../.
```
Custom configuration of compiler options may be done using
```
ccmake ../
```
and use c (configure) t( toggle to all options) and g ( generate makefiles after configuration). 
After sucessfull generation of makefile, type
```
make -j
```
to compile the code, where -j option is for parallel compilation.
The executable **fluxo**  is contained in your FLUXO directory in `build/bin/`.


