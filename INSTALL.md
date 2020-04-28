# FLUXO installation procedure

## Prerequisites for Ubuntu Linux

The required packages are listed in the table below.
Under Ubuntu, they can be obtained using the apt environment:

    sudo apt-get install git
    

| Package          | Ubuntu 14.04    | Ubuntu 16.04    |
|:----------------:|:---------------:|:---------------:|
| git              | x               |      x          |
| cmake            | x               |      x          |
| cmake-curses-gui | x               |      x          |
| liblapack3       | x               |      x          |
| liblapack-dev    | x               |      x          |
| gfortran         | x               |      x          |
| g++              | x               |      x          |
|  mpi-default-dev | x               |      x          |
| zlib1g-dev       | -               |     x           |

Table: Required debian packages under Ubuntu.


## Prerequisites for Mac OSX

For OSX, there are a few steps to install the necessary packages. The first is to get the GNU compiler suite:

1. Install *Xcode* from the App Store (this is a fairly large download). After it is downloaded and installed, launch it and agree to the Xcode license to finalize the basic installation. Note, Xcode contains many packages already, like git and LAPACK.
2. Install the *command lines tools*. To do so, open a terminal and enter:

        sudo xcode-select --install

3. Xcode comes with gcc/g++ functionality but not Fortran. So, next, install *gfortran*. There are many ways to do this. An easy way is provided by the maintainers of gfortran who offer [Apple-style installers for macOS](https://github.com/fxcoudert/gfortran-for-macOS/releases). To verify the installation of gfortran type:

        gfortran --version

     This should return the expected compiler version. 
 
The installation of other necessary packages is eased greatly with the *Macports* tool (alternatively, *homebrew* could also be used albeit with slightly different syntax). Macports provides macOS with a Synaptic Package Manager type environment which facilitates installation and updates of software libraries.

1. Install [Macports](https://www.macports.org/install.php).
2. After installation, it is recommended to run a *self-update* to ensure that the ports available are all current:


        sudo port -v selfupdate

      This should be done periodically to keep the Macports system up-to-date.

3. Use Macports to install the remaining packages listed in the table below, which are obtained through the port environment:

        sudo port install cmake

     The port environment will also install any necessary supporting packages that are required.

  | Package |
  |:-------:|
  | cmake   |
  | ctags   |
  | hdf5    |
  | mpich   |
  | openmpi |

Table: Remaining packages to be installed with Macports.

A list of ports that are installed as well as their version is provided with the command

    port installed


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


