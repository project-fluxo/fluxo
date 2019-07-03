# =========================================================================
# FLUXO CMake configuration
# =========================================================================
# NOTE: use minimum CMake version required by tools/libs (Paraview, HDF5, FFTW)
CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
# Policies (see https://cmake.org/cmake/help/v3.8/manual/cmake-policies.7.html)
cmake_policy (SET CMP0045 OLD)
cmake_policy (SET CMP0046 OLD)
# Folder for custom CMake source files
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)

# =========================================================================
# Check whether we are using Git ssh or https
# =========================================================================
EXECUTE_PROCESS(COMMAND git ls-remote --get-url OUTPUT_VARIABLE GIT_ORIGIN)
STRING(FIND "${GIT_ORIGIN}" "git@" USESSH)
IF(USESSH EQUAL 0)
  SET(FLUXO_DLPATH "git@github.com:project-fluxo")
ELSE()
  SET(FLUXO_DLPATH "https://github.com/project-fluxo")
ENDIF()
SET (HDF5DOWNLOAD https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git )
SET (HDF5TAG "hdf5-1_10_1" CACHE STRING "HDF5 version tag")
SET (P4ESTDOWNLOAD https://github.com/cburstedde/p4est.git )
SET (P4ESTSHA "7ce0525ffcd3f9cd0f5834604434b1d9b4e1ba57" CACHE STRING "P4EST commit SHA")

# =========================================================================
# Detect machine environements
# =========================================================================
SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "eslogin" ONHORNET)
STRING(FIND "${CMAKE_HOSTNAME}" "marconi" ON_MARCONI)
# Set options that need to be enabled BEFORE ANY PROJECT OR LANGUAGE COMMAND
IF (ONHORNET GREATER -1)
  SET(CMAKE_C_COMPILER       cc)
  SET(CMAKE_CXX_COMPILER     CC)
  SET(CMAKE_Fortran_COMPILER ftn)
  SET(HDF5_USE_STATIC_LIBRARIES TRUE)
  MARK_AS_ADVANCED(C_PATH CXX_PATH Fortran_PATH)
  MESSAGE(STATUS "Compiling on Hazelhen")
  ADD_DEFINITIONS(-DLUSTRE)
ENDIF()
IF (ON_MARCONI GREATER -1)
  MESSAGE(STATUS "Compiling on Marconi")
ENDIF()

# =========================================================================
# Project definitions and output paths
# =========================================================================
PROJECT(Fluxo)
ENABLE_LANGUAGE(Fortran C CXX)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Set custom install dir (needs to be done after project has been defined!)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)

# =========================================================================
# Options
# =========================================================================

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

#generate lists of subfolders to be used as choices for the variables of the configuration
SUBDIRLIST(TESTCASES  ${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/)
SUBDIRLIST(EQNSYSTEMS ${CMAKE_CURRENT_SOURCE_DIR}/src/equation/)
SUBDIRLIST(LIFTINGS   ${CMAKE_CURRENT_SOURCE_DIR}/src/dg/)


OPTION(FLUXO_BUILD_MPI "Compile SINGLE or MPI version" ON)
#OPTION(FLUXO_BUILD_MKL "Use MKL libraries instead of BLAS/LAPACK" OFF)
OPTION(FLUXO_BUILD_PAPI "Use PAPI library to perform performance measurements (e.g. flop counts)." OFF)
OPTION(FLUXO_BUILD_HDF5 "Compile and build HDF5 library" OFF)
OPTION(FLUXO_AMR        "Compile and build with p4est / AMR support" ON)

OPTION(FLUXO_PARABOLIC "Enable parabolic terms of equation-system" OFF)

SET(   FLUXO_EQNSYSNAME "linearscalaradvection" CACHE STRING
                        "Choice of equation system (linearscalaradvection/navierstokes/mhd/maxwell)" )
SET_PROPERTY(CACHE FLUXO_EQNSYSNAME PROPERTY STRINGS ${EQNSYSTEMS})

SET(   FLUXO_TESTCASE   "default"      CACHE STRING
                        "Name of testcase which is used. Must be identical to directory name in testcase folder.")
SET_PROPERTY(CACHE FLUXO_TESTCASE PROPERTY STRINGS ${TESTCASES})

SET(FLUXO_POLYNOMIAL_DEGREE "N" CACHE STRING 
                        "N: compile with variable polynomial degree. 1,2,...: fix polynomial degree (untested).")
ADD_DEFINITIONS(-DPP_N=${FLUXO_POLYNOMIAL_DEGREE})
MARK_AS_ADVANCED(FLUXO_POLYNOMIAL_DEGREE)

SET(_BUILD_FLUXO_POST  OFF  CACHE BOOL "switch to also build the fluxo_post executable (reads a collection of states and calls analyze).")

# =========================================================================
# PREPROCESSOR FLAGS
# =========================================================================
SET(FLUXO_DISCTYPE "2" CACHE STRING "disc type for volint. 1: standard form, 2: flux-average form")
SET_PROPERTY(CACHE FLUXO_DISCTYPE  PROPERTY STRINGS "1" "2")
ADD_DEFINITIONS(-DPP_DiscType=${FLUXO_DISCTYPE})

SET(FLUXO_DISC_CARTESIANFLUX  OFF  CACHE BOOL "switch on only if cartesian meshes are used.")
IF(FLUXO_DISC_CARTESIANFLUX)
  ADD_DEFINITIONS(-DCARTESIANFLUX=${FLUXO_DISC_CARTESIANFLUX} )
ENDIF()

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  UNSET(FLUXO_OPTIMIZED CACHE) # unset if debug
ELSE()
  SET(FLUXO_OPTIMIZED  ON  CACHE BOOL "switch on/off some optimized implementations.")
ENDIF()
IF(FLUXO_OPTIMIZED)
  ADD_DEFINITIONS(-DOPTIMIZED)
ENDIF()

##### Node type
SET(FLUXO_DISC_NODETYPE "GAUSS-LOBATTO" CACHE STRING "Set of nodes used for representing the solution (GAUSS/GAUSS-LOBATTO).")
SET_PROPERTY(CACHE FLUXO_DISC_NODETYPE PROPERTY STRINGS "GAUSS" "GAUSS-LOBATTO")
IF(FLUXO_DISC_NODETYPE STREQUAL "GAUSS")
  ADD_DEFINITIONS(-DPP_NodeType=1)
ELSEIF(FLUXO_DISC_NODETYPE STREQUAL "GAUSS-LOBATTO")
  ADD_DEFINITIONS(-DPP_NodeType=2)
ELSE()
  MESSAGE(SEND_ERROR "No such node type known.")
ENDIF()

IF(FLUXO_DISC_NODETYPE STREQUAL "GAUSS")
  IF(FLUXO_DISCTYPE EQUAL 2)
    MESSAGE(SEND_ERROR " Split form discretization (DiscType=2) + Gauss nodes is not possible! Choose either DiscType=2 + Gauss-Lob., or DiscType=1 + Gauss / Gauss-Lob. ")
  ENDIF()
ENDIF()

IF(FLUXO_AMR)
  ADD_DEFINITIONS(-DUSE_AMR=1)
  CMAKE_DEPENDENT_OPTION(FLUXO_BUILD_P4EST "Build p4est library" ON "FLUXO_AMR" OFF)
ELSE()
  ADD_DEFINITIONS(-DUSE_AMR=0)
ENDIF()

# =========================================================================
# Equation system
# =========================================================================
# variables set in an if clause needs to be unset by default

SET(EQNSYSNAME_FOUND false ) #to be set in each equation system
                             # so that its clear that the equation system has been accounted for here!

IF(FLUXO_EQNSYSNAME STREQUAL "navierstokes")
  SET(EQNSYSNAME_FOUND true )
  ADD_DEFINITIONS(-DPP_nVar=5 -Dnavierstokes)
ENDIF()

IF(FLUXO_EQNSYSNAME STREQUAL "linearscalaradvection")
  SET(EQNSYSNAME_FOUND true )
  ADD_DEFINITIONS(-DPP_nVar=1 -Dlinearscalaradvection)
ENDIF()

IF(FLUXO_EQNSYSNAME STREQUAL "mhd")
  SET(EQNSYSNAME_FOUND true )
  ADD_DEFINITIONS( -Dmhd)
  SET(FLUXO_EQN_GLM  ON  CACHE BOOL "switch on/off divergence cleaing for mhd.")
  IF(FLUXO_EQN_GLM)
    ADD_DEFINITIONS(-DPP_GLM=T -DPP_nVar=9)
  ELSE()
    ADD_DEFINITIONS( -DPP_nVar=8)
  ENDIF()
  SET(FLUXO_EQN_ANISO_HEAT  OFF  CACHE BOOL "switch on anisotropic heat diffusion in mhd equation.")
  IF(FLUXO_EQN_ANISO_HEAT)
    ADD_DEFINITIONS(-DPP_ANISO_HEAT=1)
  ENDIF()
ELSE()
  UNSET(FLUXO_EQN_GLM CACHE ) # if not mhd
  UNSET(FLUXO_EQN_ANISO_HEAT CACHE ) #if not mhd
ENDIF()

IF(FLUXO_EQNSYSNAME STREQUAL "maxwell")
  SET(EQNSYSNAME_FOUND true )
  ADD_DEFINITIONS(-DPP_nVar=8 -Dmaxwell)
  IF (FLUXO_PARABOLIC)
    MESSAGE(SEND_ERROR " Do not use parabolic flag with maxwell equation system!")
  ENDIF()
ENDIF()

IF(NOT EQNSYSNAME_FOUND ) 
  # no given equation system
  MESSAGE(SEND_ERROR " This equation system does not appear in CMakeLists.txt!")
ENDIF()


# variables set in an if clause needs to be unset by default

IF (FLUXO_PARABOLIC)
  ADD_DEFINITIONS(-DPARABOLIC=1)

  #### Set the lifting
  SET(FLUXO_PARABOLIC_LIFTING "br2" CACHE STRING "Lifting method used for PARABOLIC")
  SET_PROPERTY(CACHE FLUXO_PARABOLIC_LIFTING PROPERTY STRINGS ${LIFTINGS})

  IF(FLUXO_PARABOLIC_LIFTING STREQUAL "br1")
    ADD_DEFINITIONS(-DPP_Lifting=1)
  ELSEIF(FLUXO_PARABOLIC_LIFTING STREQUAL "br2")
    ADD_DEFINITIONS(-DPP_Lifting=2)
  ELSE()
    MESSAGE(SEND_ERROR " No such LIFTING known.")
  ENDIF()

  #### Set the law for viscosity
  IF(FLUXO_EQNSYSNAME STREQUAL "navierstokes")
    SET(FLUXO_EQN_VISCOSITY "constant" CACHE STRING 
                      "Method to calculate viscosity (constant/sutherland/powerlaw)")
    SET_PROPERTY(CACHE FLUXO_EQN_VISCOSITY PROPERTY STRINGS constant sutherland powerlaw)
    IF(FLUXO_EQN_VISCOSITY STREQUAL "constant")
      ADD_DEFINITIONS(-DPP_VISC=0)
    ELSEIF(FLUXO_EQN_VISCOSITY STREQUAL "sutherland")
      ADD_DEFINITIONS(-DPP_VISC=1)
    ELSEIF(FLUXO_EQN_VISCOSITY STREQUAL "powerlaw")
      ADD_DEFINITIONS(-DPP_VISC=2)
    ELSE()
      MESSAGE(SEND_ERROR " No such VISCOSITY known.")
    ENDIF()
  ELSE()
    UNSET(FLUXO_EQN_VISCOSITY CACHE ) #if parabolic & not navierstokes
  ENDIF()
ELSE()
  UNSET(FLUXO_EQN_VISCOSITY CACHE ) #if not parabolic & not navierstokes
  UNSET(FLUXO_PARABOLIC_LIFTING CACHE ) #if not parabolic
  ADD_DEFINITIONS(-DPARABOLIC=0)
ENDIF()

# =========================================================================
# CMake option output
# =========================================================================

# A function to get all user defined variables with a specified prefix
function (getListOfVarsStartingWith _prefix _varResult)
    GET_CMAKE_PROPERTY(_vars CACHE_VARIABLES)
    STRING(REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    SET(${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

# A function to get all user defined variables with a specified prefix
SET(configuration ${CMAKE_Fortran_MODULE_DIRECTORY}/configuration-cmake.f90)
FILE(WRITE ${configuration} "!This file is included during compilation in fluxo.f90 to show the configuration on the screen output. \n")
FILE(APPEND ${configuration} "  WRITE(UNIT_stdOut,'(A)') \" THIS IS AN OUTPUT OF THE FLUXO OPTIONS USED IN CMAKE:\" \n")
GETLISTOFVARSSTARTINGWITH("FLUXO_" fluxoVars)
FOREACH (_var IN LISTS fluxoVars)
  GET_PROPERTY(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
  IF(${_var})
    STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
  ENDIF()
  FILE(APPEND ${configuration} "  WRITE(UNIT_stdOut,'(A18,A30,A3,A40)')  \"   CMAKE-OPTION : \",\" ${_var} \",\" : \",\"${${_var}}\"\n")
ENDFOREACH()


# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug / Release / Profile."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_DEFINITIONS("-DDEBUG")
ENDIF()

# =========================================================================
# Location of binary and filenames
# =========================================================================

# append relative filename-macro for __FILENAME__ in Stamp of abort function (see defines.h)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# add basedir for location of corresponding userblock-file
ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

# =========================================================================
# COMPILER FLAGS
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  #environment specific variables
  IF(ONHORNET GREATER -1)
    SET(FLUXO_INSTRUCTION "-march=core-avx2")
  ELSE()
    SET(FLUXO_INSTRUCTION "-march=native")
  ENDIF()
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 ${FLUXO_INSTRUCTION} -finline-functions ") #gives segfault: -fstack-arrays
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${FLUXO_INSTRUCTION} -finline-functions ") #gives segfault: -fstack-arrays
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -ffpe-trap=invalid -fbounds-check -finit-real=snan -fbacktrace -Wall -Wextra")
  # add flags only for compiling not linking!
  SET (FLUXO_COMPILE_FLAGS "-xf95-cpp-input -fPIC")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   #environment specific variables
   IF(ONHORNET GREATER -1)
      SET(FLUXO_INSTRUCTION "-xCORE-AVX2")
   ELSEIF(ON_MARCONI GREATER -1)
      SET(FLUXO_INSTRUCTION "-axMIC-AVX512 ") #for BW and KNL
   ELSE()
      SET(FLUXO_INSTRUCTION "-xHost")
   ENDIF()
   # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 ${FLUXO_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${FLUXO_INSTRUCTION} -qopt-report0 -qopt-report-phase=vec -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
   # add flags only for compiling not linking!
   SET (FLUXO_COMPILE_FLAGS "-fpp -allow nofpp_comments -assume bscc")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -f free -s real64 -em -DCRAY")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -p . -rm")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -h profile_generate -p . -rm")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -eD -rm")
  # add flags only for compiling not linking!
  SET (FLUXO_COMPILE_FLAGS "-F")

ELSE()
  MESSAGE(SEND_ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# Add the libraries
# =========================================================================
IF(FLUXO_BUILD_MPI)
  SET(FLUXO_FC  ${MPI_Fortran_COMPILER})
  SET(FLUXO_CC  ${MPI_C_COMPILER})
  SET(FLUXO_CXX ${MPI_CXX_COMPILER})
ELSE()
  SET(FLUXO_FC  ${CMAKE_Fortran_COMPILER})
  SET(FLUXO_CC  ${CMAKE_C_COMPILER} )
  SET(FLUXO_CXX ${CMAKE_CXX_COMPILER} )
ENDIF()
SET(FLUXO_COMPILERS FC=${FLUXO_FC} CC=${FLUXO_CC} CXX=${FLUXO_CXX})


IF(FLUXO_BUILD_MPI)
  SET(FLUXO_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID}-MPI)
ELSE()
  SET(FLUXO_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID})
ENDIF()
MARK_AS_ADVANCED(FLUXO_EXTERNAL_LIB_DIR)

# =========================================================================
# p4est library
# =========================================================================

IF(FLUXO_AMR)
  IF(FLUXO_BUILD_P4EST)
    SET(FLUXO_P4EST_DLDIR ${FLUXO_EXTERNAL_LIB_DIR}/P4EST)
    SET(FLUXO_P4EST_DIR  ${FLUXO_P4EST_DLDIR}/build)
    IF (NOT EXISTS "${FLUXO_P4EST_DIR}/lib/libp4est.a")
      IF(FLUXO_BUILD_MPI)
        SET(FLUXO_P4ESTPARALLEL --enable-mpi)
      ELSE()
        UNSET(FLUXO_P4ESTPARALLEL)
      ENDIF()
      ExternalProject_Add(P4EST
        GIT_REPOSITORY ${P4ESTDOWNLOAD}
        GIT_TAG ${P4ESTTAG}
        PREFIX ${FLUXO_P4EST_DIR}
        UPDATE_COMMAND "${FLUXO_P4EST_DIR}/src/P4EST/bootstrap"
        CONFIGURE_COMMAND ${FLUXO_COMPILERS} ${FLUXO_P4EST_DIR}/src/P4EST/configure --prefix=${FLUXO_P4EST_DIR} --with-pic ${FLUXO_P4ESTPARALLEL}
        BUILD_COMMAND ${MAKE}
      )
    ENDIF()
    UNSET(P4EST_LIBRARIES)
    UNSET(P4EST_INCLUDE_DIRS)
    LIST(APPEND P4EST_LIBRARIES    ${FLUXO_P4EST_DIR}/lib/libp4est.a; ${FLUXO_P4EST_DIR}/lib/libsc.a)
    LIST(APPEND P4EST_INCLUDE_DIRS ${FLUXO_P4EST_DIR}/include )
  ELSE()
    FIND_PACKAGE(p4est)
  ENDIF()
  
  INCLUDE_DIRECTORIES(${P4EST_INCLUDE_DIRS})
  LIST(APPEND linkedlibs ${P4EST_LIBRARIES} -ldl)
ENDIF()

# =========================================================================
# HDF5 library
# =========================================================================

IF(NOT FLUXO_BUILD_HDF5)
  FIND_PACKAGE(HDF5 COMPONENTS C Fortran QUIET) #if built with configure script
  IF (HDF5_FOUND)
    SET(FLUXO_BUILD_HDF5_CMAKE FALSE)
  ELSE()
    FIND_PACKAGE(HDF5 COMPONENTS C Fortran static NO_MODULE QUIET PATH_SUFFIXES share/cmake )
    IF (HDF5_FOUND)
      SET(FLUXO_BUILD_HDF5_CMAKE TRUE)
      MESSAGE (STATUS "HDF5 C libs:${HDF5_FOUND} static:${HDF5_static_C_FOUND} and shared:${HDF5_shared_C_FOUND}")
      MESSAGE (STATUS "HDF5 Fortran libs: static:${HDF5_static_Fortran_FOUND} and shared:${HDF5_shared_Fortran_FOUND}")
      MESSAGE (STATUS "h5diff found:${HDF5_DIFF_EXECUTABLE}")
    ENDIF()
  ENDIF()
  UNSET(FLUXO_HDF5_DLDIR CACHE ) #if not build hdf5
  UNSET(FLUXO_HDF5_DIR CACHE ) #if not build hdf5
ELSE()

  SET(FLUXO_HDF5_DLDIR ${FLUXO_EXTERNAL_LIB_DIR}/HDF5)
  SET(FLUXO_HDF5_DIR  ${FLUXO_HDF5_DLDIR}/build)
  IF (NOT EXISTS "${FLUXO_HDF5_DIR}/lib/libhdf5.a")
    IF(FLUXO_BUILD_MPI)
      SET(FLUXO_HDF5PARALLEL --enable-parallel)
    ELSE()
      UNSET(FLUXO_HDF5PARALLEL)
    ENDIF()
    #NOTE: old 1.8.X versions of HDF5 require additional flag --enable-fortran2003
    ExternalProject_Add(HDF5
      GIT_REPOSITORY ${HDF5DOWNLOAD}
      GIT_TAG ${HDF5TAG}
      PREFIX ${FLUXO_HDF5_DIR}
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND ${FLUXO_COMPILERS} ${FLUXO_HDF5_DIR}/src/HDF5/configure --prefix=${FLUXO_HDF5_DIR} --with-pic --enable-fortran ${FLUXO_HDF5PARALLEL} --disable-shared --libdir=${FLUXO_HDF5_DIR}/lib
      BUILD_COMMAND ${MAKE}
    )
  ENDIF()

  SET(HDF5_hdf5_LIBRARY_RELEASE         ${FLUXO_HDF5_DIR}/lib/libhdf5.a)
  SET(HDF5_hdf5_fortran_LIBRARY_RELEASE ${FLUXO_HDF5_DIR}/lib/libhdf5_fortran.a)
  SET(HDF5_C_INCLUDE_DIR                ${FLUXO_HDF5_DIR}/include)
  SET(HDF5_Fortran_INCLUDE_DIR          ${FLUXO_HDF5_DIR}/include)
  UNSET(HDF5_LIBRARIES)
  UNSET(HDF5_INCLUDE_DIR_FORTRAN)
  UNSET(HDF5_INCLUDE_DIR)
  UNSET(HDF5_DIFF_EXECUTABLE)
  SET(HDF5_DIFF_EXECUTABLE              ${FLUXO_HDF5_DIR}/bin/h5diff)
  LIST(APPEND HDF5_INCLUDE_DIR_FORTRAN ${HDF5_Fortran_INCLUDE_DIR} ${HDF5_C_INCLUDE_DIR})
  LIST(APPEND HDF5_INCLUDE_DIR  ${HDF5_C_INCLUDE_DIR})
  FIND_LIBRARY(HDF5_z_LIBRARY_RELEASE z)
  LIST(APPEND HDF5_LIBRARIES ${HDF5_hdf5_fortran_LIBRARY_RELEASE} ${HDF5_hdf5_LIBRARY_RELEASE} ${HDF5_z_LIBRARY_RELEASE} -ldl)
ENDIF()

IF(FLUXO_BUILD_HDF5_CMAKE)
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR_FORTRAN} ${HDF5_INCLUDE_DIR})
  LIST(APPEND linkedlibs ${HDF5_FORTRAN_STATIC_LIBRARY} ${HDF5_LIBRARIES})
ELSE()
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR_FORTRAN} ${HDF5_INCLUDE_DIR})
  MESSAGE(STATUS "Found HDF5 Libs: " ${HDF5_LIBRARIES})
  LIST(APPEND linkedlibs ${HDF5_LIBRARIES} )
ENDIF()

# =========================================================================
# LAPACK/MKL library
# =========================================================================
# use Lapack/Blas for GNU
#MESSAGE(STATUS "Using BLAS/Lapack library")
#IF (FLUXO_BUILD_MKL)
#  SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
#ENDIF()
#FIND_PACKAGE(LAPACK QUIET)
#IF (LAPACK_FOUND)
#  LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
#ENDIF()

# =========================================================================
# MPI
# =========================================================================
IF(FLUXO_BUILD_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF (NOT MPI_Fortran_NO_INTERROGATE)
    FOREACH(DIR ${MPI_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    FOREACH(DIR ${MPI_Fortran_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    LIST(APPEND linkedlibs ${MPI_Fortran_LIBRARIES})
  ENDIF()
  MARK_AS_ADVANCED(MPI_LIBRARY MPI_EXTRA_LIBRARY) # don't show these ones
  MESSAGE(STATUS "Building Fluxo with MPI")
  MESSAGE(STATUS "MPI Compiler: " ${MPI_Fortran_COMPILER})
  ADD_DEFINITIONS(-DMPI=1)
ELSE()
  ADD_DEFINITIONS(-DMPI=0)
ENDIF()

# ========================================================================= 
# PAPI library 
# ========================================================================= 
IF(FLUXO_BUILD_PAPI)
  FIND_PACKAGE(PAPI REQUIRED)
  ADD_DEFINITIONS(-DPAPI) 
  LIST(APPEND linkedlibs ${PAPI_LIBRARIES})
  INCLUDE_DIRECTORIES(${PAPI_INCLUDE_DIRS})
  MESSAGE(STATUS "Building Fluxo with PAPI benchmark support.") 
ENDIF()

# =========================================================================
# Fluxo sources & targets
# =========================================================================

function(add_exec target)
  ADD_EXECUTABLE(${target} ${ARGN})
  TARGET_LINK_LIBRARIES(${target} fluxolib ${linkedlibs})
  #TARGET_LINK_LIBRARIES(${target} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.o fluxolib ${linkedlibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${FLUXO_COMPILE_FLAGS}) 
  SET_TARGET_PROPERTIES(${target} PROPERTIES LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
endfunction()

# include defines.h 
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_Fortran_MODULE_DIRECTORY}) 

FILE(GLOB_RECURSE fluxoF90 ./src/globals/*.f90
                           ./src/equation/${FLUXO_EQNSYSNAME}/*.f90
                           ./src/testcase/${FLUXO_TESTCASE}/*.f90
                           ./src/readintools/*.f90
                           ./src/interpolation/*.f90
                           ./src/mpi/*.f90
                           ./src/restart/*.f90
                           ./src/io_hdf5/*.f90 
                           ./src/analyze/*.f90
                           ./src/output/*.f90
                           ./src/mortar/*.f90
                           ./src/mesh/*.f90
                           ./src/timedisc/*.f90            )


# collect F90 files of 'dg' directory (non-recursive) in a new list 'dgF90'
FILE(GLOB dgF90            ./src/dg/*.f90)
# append F90 files of lifting to 'dgF90'-list
IF (FLUXO_PARABOLIC)
  FILE(GLOB liftingF90     ./src/dg/${FLUXO_PARABOLIC_LIFTING}/*.f90)
  LIST(APPEND dgF90 ${liftingF90})
ENDIF()

# append 'dgF90'-list to 'fluxoF90'-list
LIST(APPEND fluxoF90 ${dgF90})

IF (FLUXO_AMR)
  INCLUDE_DIRECTORIES(./src/amr/)
  FILE(GLOB_RECURSE amrF90 ./src/amr/*.f90)
  FILE(GLOB_RECURSE amrC   ./src/amr/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/globals/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/globals/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/analyze/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/analyze/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/interpolation/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/interpolation/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/mpi/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/mpi/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/p4est/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/p4est/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/refine/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/refine/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/io_hdf5/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/io_hdf5/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/mesh/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/mesh/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/output/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/output/*.c)
#  FILE(GLOB amrF90 ./src/amr/hopest/readintools/*.f90)
#  FILE(GLOB amrC   ./src/amr/hopest/readintools/*.c)
  LIST(APPEND fluxoF90 ${amrF90})
  LIST(APPEND fluxoC   ${amrC})
ENDIF()

# define libs
ADD_LIBRARY(fluxolibF90 OBJECT ${fluxoF90} )
ADD_DEPENDENCIES(fluxolibF90 HDF5)
SET_TARGET_PROPERTIES(fluxolibF90 PROPERTIES COMPILE_FLAGS ${FLUXO_COMPILE_FLAGS}) 
LIST(APPEND fluxoobj $<TARGET_OBJECTS:fluxolibF90>)
IF (fluxoC)
  ADD_LIBRARY(fluxolibC OBJECT ${fluxoC} )
  ADD_DEPENDENCIES(fluxolibC HDF5 P4EST)
  SET_TARGET_PROPERTIES(fluxolibC PROPERTIES COMPILE_FLAGS "-fPIC") # TODO: add C default flags
  LIST(APPEND fluxoobj $<TARGET_OBJECTS:fluxolibC>)
ENDIF()

ADD_LIBRARY(fluxolib STATIC ${fluxoobj})
LIST(APPEND fluxodeps fluxolib HDF5)
IF (FLUXO_AMR)
  LIST(APPEND fluxodeps P4EST)
ENDIF()

add_exec(fluxo ./src/fluxo.f90)
ADD_DEPENDENCIES(fluxo ${fluxodeps})

IF (_BUILD_FLUXO_POST)
  add_exec(fluxo_post ./src/fluxo_post.f90)
  ADD_DEPENDENCIES(fluxo_post ${fluxodeps})
ENDIF()

# special compiler flags to avoid warnings for empty functions of testcase.f90 and overwrite of intrinsic abort in globals.f90
#...# IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#...#   SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/${FLUXO_TESTCASE}/testcase.f90 PROPERTIES COMPILE_FLAGS "-Wno-unused-dummy-argument")
#...#   SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/globals/globals.f90 PROPERTIES COMPILE_FLAGS "-Wno-intrinsic-shadow")
#...# ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#...#   SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/${FLUXO_TESTCASE}/testcase.f90 PROPERTIES COMPILE_FLAGS "-warn nounused")
#...# ENDIF()

# -------------------------------------------------------------------------
# echo COMPILE_DEFINITIONS
# -------------------------------------------------------------------------
GET_DIRECTORY_PROPERTY(comp_defs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
SET(comp_defs_str "COMPILE_DEFINITIONS = ")
FOREACH (def ${comp_defs})
   SET(comp_defs_str "${comp_defs_str} -D${def}")
ENDFOREACH()
ADD_CUSTOM_TARGET(preproc_defines COMMAND echo ${comp_defs_str})

# =========================================================================
# Install
# =========================================================================

INSTALL(FILES    ${CMAKE_Fortran_MODULE_DIRECTORY}/configuration-cmake.f90 DESTINATION include )
