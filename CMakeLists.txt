# =========================================================================
# FLUXO CMake configuration
# =========================================================================
# NOTE: use minimum CMake version required by tools/libs (Paraview, HDF5, FFTW)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1)
# Policies (see https://cmake.org/cmake/help/v3.13/manual/cmake-policies.7.html)
# VERSION_GREATER_EQUAL is incompatible with CMAKE versions < 3.7
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.6.0"))
  SET(GITSHALLOW "GIT_SHALLOW ON")
ENDIF()
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.12.0"))
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()
# Folder for custom CMake source files
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)



# =========================================================================
# Detect machine environments
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
ENABLE_LANGUAGE(Fortran)

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

SUBDIRLIST(EQNSYSTEMS ${CMAKE_CURRENT_SOURCE_DIR}/src/equation/)
SUBDIRLIST(TESTCASES  ${CMAKE_CURRENT_SOURCE_DIR}/src/testcase/)
SUBDIRLIST(LIFTINGS   ${CMAKE_CURRENT_SOURCE_DIR}/src/dg/)


OPTION(FLUXO_BUILD_MPI "Compile SINGLE or MPI version" ON)
#OPTION(FLUXO_BUILD_MKL "Use MKL libraries instead of BLAS/LAPACK" OFF)
OPTION(FLUXO_BUILD_PAPI "Use PAPI library to perform performance measurements (e.g. flop counts)." OFF)
OPTION(FLUXO_BUILD_HDF5 "Compile and build hdf5 library" OFF)

OPTION(FLUXO_PARABOLIC "Enable parabolic terms of equation-system" OFF)

SET(   FLUXO_EQNSYSNAME "linearscalaradvection" CACHE STRING
                        "Choice of equation system (lienarscalaradvection/navierstokes/mhd/maxwell)" )
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

SET(FLUXO_DISC_NODETYPE "GAUSS-LOBATTO" CACHE STRING "Set of nodes used for representing the solution (GAUSS/GAUSS-LOBATTO).")
SET_PROPERTY(CACHE FLUXO_DISC_NODETYPE PROPERTY STRINGS GAUSS GAUSS-LOBATTO)
IF(FLUXO_DISC_NODETYPE STREQUAL "GAUSS")
  ADD_DEFINITIONS(-DPP_NodeType=1)
ELSEIF(FLUXO_DISC_NODETYPE STREQUAL "GAUSS-LOBATTO")
  ADD_DEFINITIONS(-DPP_NodeType=2)
ELSE()
  MESSAGE(SEND_ERROR " No such node type known.")
ENDIF()

IF(FLUXO_DISC_NODETYPE STREQUAL "GAUSS")
  IF(FLUXO_DISCTYPE EQUAL 2)
    MESSAGE(SEND_ERROR " Split form discretization (DiscType=2) + Gauss nodes is not possible! Choose either DiscType=2 + Gauss-Lob., or DiscType=1 + Gauss / Gauss-Lob. ")
  ENDIF()
ENDIF()

# Enable or disable unit tests
#OPTION(FLUXO_UNITTESTS "Enable unit tests after build" ON)
#MARK_AS_ADVANCED(FLUXO_UNITTESTS)

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

  #### Set the viscosity
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
#ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

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
  SET (FLUXO_COMPILE_FLAGS "-xf95-cpp-input")

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
   SET (FLUXO_COMPILE_FLAGS "-fpp -assume bscc")

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

IF(FLUXO_BUILD_HDF5)
  IF(FLUXO_BUILD_MPI)
    SET(FLUXO_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID}-MPI)
  ELSE()
    SET(FLUXO_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID})
  ENDIF()
ELSE()
  UNSET(FLUXO_EXTERNAL_LIB_DIR CACHE ) 
ENDIF()
MARK_AS_ADVANCED(FLUXO_EXTERNAL_LIB_DIR)

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
  UNSET(HDF5TAG CACHE ) #if not build hdf5
  UNSET(HDF5DOWNLOAD CACHE ) #if not build hdf5
  UNSET(FLUXO_HDF5_DLDIR CACHE ) #if not build hdf5
  UNSET(FLUXO_HDF5_DIR CACHE ) #if not build hdf5
  UNSET(FLUXO_HDF5FC CACHE ) 
  UNSET(FLUXO_HDF5CC CACHE ) 
ELSE()
  # define here the path to download from
  SET (HDF5DOWNLOAD https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git )
  SET (HDF5TAG "hdf5-1_10_0-patch1" CACHE STRING "HDF5 version tag")

  SET(FLUXO_HDF5_DLDIR ${FLUXO_EXTERNAL_LIB_DIR}/HDF5)
  SET(FLUXO_HDF5_DIR  ${FLUXO_HDF5_DLDIR}/build)
  IF (NOT EXISTS "${FLUXO_HDF5_DIR}/lib/libhdf5.a")
    IF(FLUXO_BUILD_MPI)
      SET(FLUXO_HDF5PARALLEL --enable-parallel)
      SET(FLUXO_HDF5FC ${MPI_Fortran_COMPILER})
      SET(FLUXO_HDF5CC ${MPI_C_COMPILER})
    ELSE()
      UNSET(FLUXO_HDF5PARALLEL)
      SET(FLUXO_HDF5FC ${CMAKE_Fortran_COMPILER})
      SET(FLUXO_HDF5CC ${CMAKE_C_COMPILER} )
    ENDIF()
    #NOTE: old 1.8.X versions of HDF5 require additional flag --enable-fortran2003
    ExternalProject_Add(HDF5
      GIT_REPOSITORY ${HDF5DOWNLOAD}
      GIT_TAG ${HDF5TAG}
      PREFIX ${FLUXO_HDF5_DIR}
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND FC=${FLUXO_HDF5FC} CC=${FLUXO_HDF5CC} ${FLUXO_HDF5_DIR}/src/HDF5/configure --prefix=${FLUXO_HDF5_DIR} --enable-fortran ${FLUXO_HDF5PARALLEL} --disable-shared --libdir=${FLUXO_HDF5_DIR}/lib
      BUILD_COMMAND ${MAKE}
    )
    LIST(APPEND SELFBUILTEXTERNALS HDF5)
    # CMake HDF5 is fast but not that reliable
    #ExternalProject_Add(HDF5
    #  GIT_REPOSITORY ${HDF5DOWNLOAD}
    #  GIT_TAG ${HDF5TAG}
    #  PREFIX ${FLUXO_HDF5_DIR}
    #  UPDATE_COMMAND ""
    #  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${FLUXO_HDF5_DIR} -DCMAKE_C_COMPILER=${FLUXO_HDF5CC} -DCMAKE_Fortran_COMPILER=${FLUXO_HDF5FC} -DBUILD-TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_PARALLEL=ON
    #  BUILD_COMMAND ${MAKE}
    #)
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
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_FORTRAN})
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



FILE(GLOB dgF90            ./src/dg/*.f90)
IF (FLUXO_PARABOLIC)
  FILE(GLOB liftingF90     ./src/dg/${FLUXO_PARABOLIC_LIFTING}/*.f90)
  LIST(APPEND dgF90 ${liftingF90})
ENDIF()


# define libs
ADD_LIBRARY(fluxolibF90 OBJECT ${fluxoF90} ${dgF90})


SET_TARGET_PROPERTIES(fluxolibF90 PROPERTIES COMPILE_FLAGS ${FLUXO_COMPILE_FLAGS}) 
#ADD_LIBRARY(userblocklib OBJECT ./src/output/read_userblock.c)
#SET_TARGET_PROPERTIES(userblocklib PROPERTIES LINKER_LANGUAGE C)
#ADD_LIBRARY(fluxolib STATIC $<TARGET_OBJECTS:fluxolibF90> $<TARGET_OBJECTS:userblocklib>)
ADD_LIBRARY(fluxolib STATIC $<TARGET_OBJECTS:fluxolibF90> )
IF(SELFBUILTEXTERNALS)
  ADD_DEPENDENCIES(fluxolibF90 ${SELFBUILTEXTERNALS})
ENDIF()

add_exec(fluxo ./src/fluxo.f90)
ADD_DEPENDENCIES(fluxo fluxolib ${INTERNALLIBS})

IF (_BUILD_FLUXO_POST)
  add_exec(fluxo_post ./src/fluxo_post.f90)
  ADD_DEPENDENCIES(fluxo_post fluxolib ${INTERNALLIBS})
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

#INSTALL(FILES    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt DESTINATION bin)


# =========================================================================
# USERBLOCK + Preproc_flags
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
#  SET(boolian (${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF")) 
#  IF(${boolian})
#    #FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
#    FILE(APPEND ${configuration} "  WRITE(UNIT_stdOut,*)  \"   CMAKE,OPTION: ${_var} ${${_var}} \" \n")
#  ELSE()
    IF(${_var})
      STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
    ENDIF()
    #FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
    FILE(APPEND ${configuration} "  WRITE(UNIT_stdOut,'(A18,A30,A3,A40)')  \"   CMAKE-OPTION : \",\" ${_var} \",\" : \",\"${${_var}}\"\n")
#  ENDIF()
ENDFOREACH()

# =========================================================================
# Add FLUXO2VTK tool
# =========================================================================
#add_exec(fluxo2vtk ./src/fluxo2vtk.f90)
#ADD_DEPENDENCIES(fluxo2vtk fluxolib HDF5)

# =========================================================================
# RegressionCheck tool && Examples
# =========================================================================
#add_exec(regressioncheck ./regressioncheck/regressioncheck.f90
#                         ./regressioncheck/regressioncheck_vars.f90
#                         ./regressioncheck/regressioncheck_build.f90
#                         ./regressioncheck/regressioncheck_compare.f90
#                         ./regressioncheck/regressioncheck_run.f90
#                         ./regressioncheck/regressioncheck_tools.f90 )
#ADD_DEFINITIONS(-DH5DIFF="${HDF5_DIFF_EXECUTABLE}")
#SET(FLUXO_REGRESSIONCHECK_EXAMPLES "${CMAKE_CURRENT_SOURCE_DIR}/regressioncheck/examples")


# =========================================================================
# Ctags
# =========================================================================
#SET(EQNS_WO_ACTUAL ${EQNSYSTEMS})
#LIST(REMOVE_ITEM EQNS_WO_ACTUAL ${FLUXO_EQNSYSNAME})
#SET(EXLUDES "")
#FOREACH(item ${EQNS_WO_ACTUAL})
#   SET(EXLUDES "--exclude=src/equations/${item}/*" ${EXLUDES})
#ENDFOREACH()
#SET(LIFTINGS_WO_ACTUAL ${LIFTINGS})
#LIST(REMOVE_ITEM LIFTINGS_WO_ACTUAL ${FLUXO_PARABOLIC_LIFTING})
#FOREACH(item ${LIFTINGS_WO_ACTUAL})
#   SET(EXLUDES "--exclude=src/dg/${item}/*" ${EXLUDES})
#ENDFOREACH()
#SET(TESTCASES_WO_ACTUAL ${TESTCASES})
#LIST(REMOVE_ITEM TESTCASES_WO_ACTUAL ${FLUXO_TESTCASE})
#FOREACH(item ${TESTCASES_WO_ACTUAL})
#   SET(EXLUDES "--exclude=src/testcase/${item}/*" ${EXLUDES})
#ENDFOREACH()
#ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ctags -R --fortran-kinds=+i+L -f ctags.txt ${EXLUDES} src)
#ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
#FIND_PROGRAM(CTAGS_PATH ctags)
#IF(CTAGS_PATH)
#    MESSAGE(STATUS "Found ctags: ${CTAGS_PATH}")
#    ADD_DEPENDENCIES(fluxo tags)
#ENDIF(CTAGS_PATH)

# =========================================================================
# Tests
# =========================================================================

# Deactivate tests on hornet since no programs are allowed to run on frontend
#IF(ONHORNET GREATER -1)
#  SET(FLUXO_UNITTESTS OFF CACHE BOOL "Enable unit tests after build" FORCE)
#ENDIF()
#IF(FLUXO_UNITTESTS)
#  INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/unitTests/CMakeLists.txt)
#ENDIF()

