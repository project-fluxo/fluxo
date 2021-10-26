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
STRING(FIND "${CMAKE_HOSTNAME}" "cobra" ON_COBRA)
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
ENABLE_LANGUAGE(Fortran  C CXX)

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


OPTION(FLUXO_ENABLE_MPI "Compile SINGLE or MPI version" ON)
#SET(FLUXO_ENABLE_LAPACK "off" CACHE STRING "link with LAPACK, (off/on/MKL)")
#SET_PROPERTY(CACHE FLUXO_ENABLE_LAPACK PROPERTY STRINGS off on MKL )
OPTION(FLUXO_BUILD_PAPI "Use PAPI library to perform performance measurements (e.g. flop counts)." OFF)
OPTION(FLUXO_BUILD_HDF5 "Compile and build hdf5 library" OFF)
OPTION(FLUXO_AMR        "Compile and build with p4est / AMR support" OFF)
OPTION(FLUXO_PARABOLIC "Enable parabolic terms of equation-system" OFF)

OPTION(FLUXO_SHOCKCAPTURE "Enable shock capturing for non-linear equations" OFF)
OPTION(FLUXO_POSITIVITYPRES "Enable positivity preservation for non-linear equations" OFF)

SET(   FLUXO_EQNSYSNAME "navierstokes" CACHE STRING
                        "Choice of equation system (lienarscalaradvection/navierstokes/mhd/maxwell)" )
SET_PROPERTY(CACHE FLUXO_EQNSYSNAME PROPERTY STRINGS ${EQNSYSTEMS})

IF(FLUXO_EQNSYSNAME STREQUAL "maxwell")
  UNSET(FLUXO_PARABOLIC CACHE ) 
ELSE()
  OPTION(FLUXO_PARABOLIC "Enable parabolic terms of equation-system" OFF)
ENDIF()

SET(   FLUXO_TESTCASE   "default"      CACHE STRING
                        "Name of testcase which is used. Must be identical to directory name in testcase folder.")
SET_PROPERTY(CACHE FLUXO_TESTCASE PROPERTY STRINGS ${TESTCASES})

SET(FLUXO_POLYNOMIAL_DEGREE "N" CACHE STRING 
                        "N: compile with variable polynomial degree. 1,2,...: fix polynomial degree (untested).")
ADD_DEFINITIONS(-DPP_N=${FLUXO_POLYNOMIAL_DEGREE})
#MARK_AS_ADVANCED(FLUXO_POLYNOMIAL_DEGREE)

SET(_BUILD_FLUXO_POST  OFF  CACHE BOOL "switch to also build the fluxo_post executable (reads a collection of states and calls analyze).")

# =========================================================================
# PREPROCESSOR FLAGS
# =========================================================================
SET(FLUXO_DISCTYPE "2" CACHE STRING "disc type for volint. 1: standard form, 2: flux-average form")
SET_PROPERTY(CACHE FLUXO_DISCTYPE  PROPERTY STRINGS "1" "2")
ADD_DEFINITIONS(-DPP_DiscType=${FLUXO_DISCTYPE})


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

IF(FLUXO_AMR)
  ADD_DEFINITIONS(-DUSE_AMR=1)
  SET(FLUXO_BUILD_P4EST  ON  CACHE BOOL "Build p4est library.")
ELSE()
  UNSET(FLUXO_BUILD_P4EST CACHE )
  ADD_DEFINITIONS(-DUSE_AMR=0)
ENDIF()

SET(FLUXO_JESSE_MORTAR  OFF  CACHE BOOL "switch on/off the mortar implementation of Jesse's paper.")
IF(FLUXO_JESSE_MORTAR)
  ADD_DEFINITIONS(-DJESSE_MORTAR=1)
  IF(NOT FLUXO_DISCTYPE EQUAL 2)
    MESSAGE(SEND_ERROR " JESSE_MORTAR ONLY with Split form discretization (DiscType=2)!")
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
  IF(FLUXO_SHOCKCAPTURE)
	  MESSAGE(SEND_ERROR " Do not use shock capturing with linear advection equation!")
  ENDIF()
  IF(FLUXO_POSITIVITYPRES)
    MESSAGE(SEND_ERROR " Do not use positivity preservation with linear advection equation!")
  ENDIF()
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
  IF(FLUXO_DISCTYPE EQUAL 2)
    SET(FLUXO_EQN_NONCONS  ON  CACHE BOOL "switch on/off non-conservative source term for mhd.")
    IF(FLUXO_EQN_NONCONS)
      SET(FLUXO_EQN_NONCONS_TYPE "Powell" CACHE STRING "type of non-conservative source term for mhd: Powell/Brackbill/Janhunen")
      SET_PROPERTY(CACHE FLUXO_EQN_NONCONS_TYPE PROPERTY STRINGS Powell Brackbill Janhunen )
      IF(FLUXO_EQN_NONCONS_TYPE STREQUAL "Powell")
        ADD_DEFINITIONS(-DNONCONS=1)
      ELSEIF(FLUXO_EQN_NONCONS_TYPE STREQUAL "Brackbill")
        ADD_DEFINITIONS(-DNONCONS=2)
      ELSEIF(FLUXO_EQN_NONCONS_TYPE STREQUAL "Janhunen")
        ADD_DEFINITIONS(-DNONCONS=3)
      ENDIF()
      IF(FLUXO_EQN_GLM)
        SET(FLUXO_EQN_NONCONS_GLM  OFF  CACHE BOOL "switch on/off non-conservative GLM variable dependent source terms for mhd.")
        IF(FLUXO_EQN_NONCONS_GLM)
          ADD_DEFINITIONS(-DPP_NC_GLM=1)
        ENDIF()
      ELSE()
        UNSET(FLUXO_EQN_NONCONS_GLM  CACHE ) #if GLM=OFF
      ENDIF()
    ELSE()
      UNSET(FLUXO_EQN_NONCONS_TYPE  CACHE ) #if  NONCONS=OFF
      UNSET(FLUXO_EQN_NONCONS_GLM  CACHE ) #if  NONCONS=OFF
    ENDIF()
  ELSE()
    UNSET(FLUXO_EQN_NONCONS CACHE ) #if not DISCTYPE=2
    UNSET(FLUXO_EQN_NONCONS_TYPE  CACHE ) #if not DISCTYPE=2
    UNSET(FLUXO_EQN_NONCONS_GLM  CACHE ) #if not DISCTYPE=2
  ENDIF()
  SET(FLUXO_EQN_ANISO_HEAT  OFF  CACHE BOOL "switch on anisotropic heat diffusion in mhd equation.")
  IF(FLUXO_EQN_ANISO_HEAT)
    ADD_DEFINITIONS(-DPP_ANISO_HEAT=1)
  ENDIF()
ELSE()
  UNSET(FLUXO_EQN_GLM CACHE ) # if not mhd
  UNSET(FLUXO_EQN_ANISO_HEAT CACHE ) #if not mhd
  UNSET(FLUXO_EQN_NONCONS CACHE ) #if not mhd
  UNSET(FLUXO_EQN_NONCONS_TYPE  CACHE ) #if not mhd
  UNSET(FLUXO_EQN_NONCONS_GLM  CACHE ) #if not mhd
ENDIF()

IF(FLUXO_EQNSYSNAME STREQUAL "maxwell")
  SET(EQNSYSNAME_FOUND true )
  ADD_DEFINITIONS(-DPP_nVar=8 -Dmaxwell)
  IF (FLUXO_PARABOLIC)
    MESSAGE(SEND_ERROR " Do not use parabolic flag with maxwell equation system!")
  ELSEIF(FLUXO_SHOCKCAPTURE)
    MESSAGE(SEND_ERROR " Do not use shock capturing with maxwell equation system!")
  ENDIF()
  IF(FLUXO_POSITIVITYPRES)
    MESSAGE(SEND_ERROR " Do not use positivity preservation with maxwell equation system!")
  ENDIF()
ENDIF()

IF(NOT EQNSYSNAME_FOUND ) 
  # no given equation system
  MESSAGE(SEND_ERROR " This equation system does not appear in CMakeLists.txt!")
ENDIF()

IF(FLUXO_DISCTYPE EQUAL 2)
  SET(FLUXO_EQN_VOLFLUX "-1" CACHE STRING 
          "Default: -1 (procedure pointer set in init). Specify two-point flux for compilation (eqn/flux_average module), for optimization!")
  ADD_DEFINITIONS(-DPP_VolFlux=${FLUXO_EQN_VOLFLUX}  )
ELSE() #if not disctype=2
  UNSET(FLUXO_EQN_VOLFLUX CACHE )
ENDIF()


# variables set in an if clause needs to be unset by default

IF (FLUXO_PARABOLIC)
  ADD_DEFINITIONS(-DPARABOLIC=1)

  #### Set the lifting
  SET(FLUXO_PARABOLIC_LIFTING "br1" CACHE STRING "Lifting method used for PARABOLIC")
  SET_PROPERTY(CACHE FLUXO_PARABOLIC_LIFTING PROPERTY STRINGS ${LIFTINGS})

  IF(FLUXO_EQNSYSNAME STREQUAL "navierstokes")
    LIST(APPEND LIFTING_VARS "entropy_var" "cons_var" "prim_var" )
  ELSEIF(FLUXO_EQNSYSNAME STREQUAL "linearscalaradvection")
    LIST(APPEND LIFTING_VARS "cons_var")
  ELSEIF(FLUXO_EQNSYSNAME STREQUAL "mhd")
    LIST(APPEND LIFTING_VARS "entropy_var" "cons_var" "prim_var" )
  ENDIF()
  
  LIST(GET LIFTING_VARS 0 LIFTVARDEF )
  SET(FLUXO_PARABOLIC_LIFTING_VAR ${LIFTVARDEF} CACHE STRING "variable used for BR1 lifting, cons/prim/entropy_var: conservative/primitive/entropy variables.")
  SET_PROPERTY(CACHE FLUXO_PARABOLIC_LIFTING_VAR PROPERTY STRINGS ${LIFTING_VARS} )
  IF(FLUXO_PARABOLIC_LIFTING_VAR STREQUAL "cons_var")
    ADD_DEFINITIONS(-DPP_Lifting_Var=1)
  ELSEIF(FLUXO_PARABOLIC_LIFTING_VAR STREQUAL "prim_var")
    ADD_DEFINITIONS(-DPP_Lifting_Var=2)
  ELSEIF(FLUXO_PARABOLIC_LIFTING_VAR STREQUAL "entropy_var")
    ADD_DEFINITIONS(-DPP_Lifting_Var=3)
  ENDIF()
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
  #if not parabolic
  UNSET(FLUXO_PARABOLIC_LIFTING_VAR CACHE )
  UNSET(FLUXO_EQN_VISCOSITY CACHE ) 
  UNSET(FLUXO_PARABOLIC_LIFTING CACHE )
  ADD_DEFINITIONS(-DPARABOLIC=0)
ENDIF()

# =========================================================================
# Shock capturing settings
# =========================================================================

IF (FLUXO_SHOCKCAPTURE)
  ADD_DEFINITIONS(-DSHOCKCAPTURE=1)
  
  #### NFVSE shock capturing
  OPTION(FLUXO_SHOCKCAP_NFVSE "Use NFVSE for shock capturing" ON)
  
  IF(FLUXO_SHOCKCAP_NFVSE)
    # NFVSE only available for Gauss-Lobatto!
    IF(FLUXO_DISC_NODETYPE STREQUAL "GAUSS")
      MESSAGE(SEND_ERROR " The NFVSE shock capturing is only available for Gauss-Lobatto nodes!")
    ENDIF()
    
    ADD_DEFINITIONS(-DSHOCK_NFVSE=1)
    OPTION(FLUXO_SHOCK_NFVSE_CORR "Blending Correction" OFF)
    IF(FLUXO_SHOCK_NFVSE_CORR)
      ADD_DEFINITIONS(-DNFVSE_CORR=1)
      OPTION(NFVSE_LOCAL_ALPHA "Local alpha??" OFF)
      IF(NFVSE_LOCAL_ALPHA)
        ADD_DEFINITIONS(-DLOCAL_ALPHA=1)
      ENDIF()
      OPTION(IDP_CHECKBOUNDS "Check boudns after IDP correction??" OFF)
      IF(IDP_CHECKBOUNDS)
        ADD_DEFINITIONS(-DIDP_CHECKBOUNDS=1)
      ENDIF()
    ELSE()
      UNSET(NFVSE_LOCAL_ALPHA CACHE )
      UNSET(IDP_CHECKBOUNDS CACHE )
    ENDIF()
  ELSE()
    UNSET(FLUXO_SHOCK_NFVSE_CORR CACHE )
    UNSET(NFVSE_LOCAL_ALPHA CACHE )
    UNSET(IDP_CHECKBOUNDS CACHE )
  ENDIF()
  
  #### None selected?
  IF((NOT FLUXO_SHOCKCAP_NFVSE))
    MESSAGE(SEND_ERROR " No shock capturing method selected.")
  ENDIF()
  
ELSE()
  #if no shock capturing
  UNSET(FLUXO_SHOCKCAP_NFVSE CACHE )
  UNSET(FLUXO_SHOCK_NFVSE_CORR CACHE )
  UNSET(NFVSE_LOCAL_ALPHA CACHE )
  UNSET(IDP_CHECKBOUNDS CACHE )
  ADD_DEFINITIONS(-DSHOCKCAPTURE=0)
ENDIF()

# =========================================================================
# Positivity preservation
# =========================================================================

IF (FLUXO_POSITIVITYPRES)
  ADD_DEFINITIONS(-DPOSITIVITYPRES=1)
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
  SET(COMPILER_ADD_FLAGS " -finline-limit=1500 " CACHE STRING "additional user defined compiler flags") 
  #environment specific variables
  IF(ONHORNET GREATER -1)
    SET(FLUXO_INSTRUCTION "-march=core-avx2")
  ELSE()
    SET(FLUXO_INSTRUCTION "-march=native")
  ENDIF()
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 ${FLUXO_INSTRUCTION} ${COMPILER_ADD_FLAGS} -finline-functions ") #gives segfault: -fstack-arrays
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 ${FLUXO_INSTRUCTION} ${COMPILER_ADD_FLAGS} -finline-functions ") #gives segfault: -fstack-arrays
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -ffpe-trap=invalid -fbounds-check -finit-real=snan -fbacktrace -Wall -Wextra")
  # add flags only for compiling not linking!
  SET (FLUXO_COMPILE_FLAGS "-xf95-cpp-input")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   SET(COMPILER_ADD_FLAGS "-inline-max-size=1500 -no-inline-max-total-size" CACHE STRING "additional user defined compiler flags") 
   #environment specific variables
   IF(ONHORNET GREATER -1)
      SET(FLUXO_INSTRUCTION "-xCORE-AVX2")
   ELSEIF(ON_MARCONI GREATER -1)
      #SET(FLUXO_INSTRUCTION "-axMIC-AVX512 ") #for BW and KNL
      SET(FLUXO_INSTRUCTION "-xCore-AVX512 ") #for SKL
      #SET(FLUXO_INSTRUCTION "-xCommon-AVX512 ") #for SKL
   ELSEIF(ON_COBRA GREATER -1)
      SET(FLUXO_INSTRUCTION "-xCore-AVX512")
   ELSE()
      SET(FLUXO_INSTRUCTION "-xHost")
   ENDIF()
   # set Flags
  SET(FLUXO_OPT_REPORT "-qopt-report=5 -qopt-report-phase=all")
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3  ${FLUXO_INSTRUCTION} ${COMPILER_ADD_FLAGS} ${FLUXO_OPT_REPORT} -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3  ${FLUXO_INSTRUCTION} ${COMPILER_ADD_FLAGS} ${FLUXO_OPT_REPORT} -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g  -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
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
IF(FLUXO_ENABLE_MPI)
  SET(FLUXO_FC  ${MPI_Fortran_COMPILER})
  SET(FLUXO_CC  ${MPI_C_COMPILER})
  SET(FLUXO_CXX ${MPI_CXX_COMPILER})
ELSE()
  SET(FLUXO_FC  ${CMAKE_Fortran_COMPILER})
  SET(FLUXO_CC  ${CMAKE_C_COMPILER} )
  SET(FLUXO_CXX ${CMAKE_CXX_COMPILER} )
ENDIF()
SET(FLUXO_COMPILERS FC=${FLUXO_FC} CC=${FLUXO_CC} CXX=${FLUXO_CXX})


# IF(FLUXO_ENABLE_MPI)
#   SET(FLUXO_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID}-MPI)
# ELSE()
#   SET(FLUXO_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID})
# ENDIF()
# MARK_AS_ADVANCED(FLUXO_EXTERNAL_LIB_DIR)

IF(FLUXO_ENABLE_MPI)
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
    SET (P4ESTDOWNLOAD https://github.com/cburstedde/p4est.git )
    # SET (P4ESTSHA "7ce0525ffcd3f9cd0f5834604434b1d9b4e1ba57" CACHE STRING "P4EST commit SHA")
    SET (P4ESTSHA "7752b1e12b596663667d5949aa922e805a42b9b6" CACHE STRING "P4EST commit SHA")
    SET(FLUXO_P4EST_DLDIR ${FLUXO_EXTERNAL_LIB_DIR}/P4EST)
    SET(FLUXO_P4EST_DIR  ${FLUXO_P4EST_DLDIR}/build)
    IF (NOT EXISTS "${FLUXO_P4EST_DIR}/lib/libp4est.so")
      IF(FLUXO_ENABLE_MPI)
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
      LIST(APPEND SELFBUILTEXTERNALS P4EST)
    ENDIF()
    UNSET(P4EST_LIBRARIES)
    UNSET(P4EST_INCLUDE_DIRS)
    LIST(APPEND P4EST_LIBRARIES    ${FLUXO_P4EST_DIR}/lib/libp4est.so; ${FLUXO_P4EST_DIR}/lib/libsc.so)
    LIST(APPEND P4EST_INCLUDE_DIRS ${FLUXO_P4EST_DIR}/include )
  ELSE()
    UNSET(P4ESTDOWNLOAD CACHE)
    UNSET(P4ESTSHA CACHE)
    
    IF(DEFINED ENV{p4est_DIR})
      SET(p4est_DIR $ENV{p4est_DIR} CACHE STRING "Installation directory of p4est")
    ELSEIF(EXISTS "${FLUXO_EXTERNAL_LIB_DIR}/P4EST/build/lib/libp4est.so")
      SET(p4est_DIR ${FLUXO_EXTERNAL_LIB_DIR}/P4EST/build CACHE STRING "Installation directory of p4est")
    ENDIF()
    IF(p4est_DIR)
      UNSET(P4EST_LIBRARIES)
      UNSET(P4EST_INCLUDE_DIRS)
      LIST(APPEND P4EST_LIBRARIES    ${p4est_DIR}/lib/libp4est.so; ${p4est_DIR}/lib/libsc.so)
      LIST(APPEND P4EST_INCLUDE_DIRS ${p4est_DIR}/include )
    ELSE()
      FIND_PACKAGE(p4est)
    ENDIF()
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
    IF(FLUXO_ENABLE_MPI)
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
#IF(NOT (FLUXO_ENABLE_LAPACK STREQUAL "off") )
#  # use Lapack/Blas for GNU
#  MESSAGE(STATUS "Using BLAS/Lapack library")
#  IF (FLUXO_ENABLE_LAPACK STREQUAL "MKL" )
#    SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
#  ENDIF()
#  FIND_PACKAGE(LAPACK QUIET)
#  IF (LAPACK_FOUND)
#    LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
#  ENDIF()
#ENDIF() 

# =========================================================================
# MPI
# =========================================================================
IF(FLUXO_ENABLE_MPI)
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
                           ./src/timedisc/*.f90   
                           ./src/indicators/*.f90
                           ./src/shockcapturing/*.f90   
       			   ./src/positivitypreservation/*.f90  )



FILE(GLOB dgF90            ./src/dg/*.f90)
IF (FLUXO_PARABOLIC)
  FILE(GLOB liftingF90     ./src/dg/${FLUXO_PARABOLIC_LIFTING}/*.f90)
  LIST(APPEND dgF90 ${liftingF90})
ENDIF()

# define libs
LIST(APPEND fluxoF90 ${dgF90})

IF (FLUXO_AMR)
  INCLUDE_DIRECTORIES(./src/amr/)
  FILE(GLOB_RECURSE amrF90 ./src/amr/*.f90)
  FILE(GLOB_RECURSE amrC   ./src/amr/*.c)
  LIST(APPEND fluxoF90 ${amrF90})
  LIST(APPEND fluxoC   ${amrC})
ENDIF()

# define libs
ADD_LIBRARY(fluxolibF90 OBJECT ${fluxoF90} )

SET_TARGET_PROPERTIES(fluxolibF90 PROPERTIES COMPILE_FLAGS ${FLUXO_COMPILE_FLAGS}) 
#ADD_LIBRARY(userblocklib OBJECT ./src/output/read_userblock.c)
#SET_TARGET_PROPERTIES(userblocklib PROPERTIES LINKER_LANGUAGE C)
#ADD_LIBRARY(fluxolib STATIC $<TARGET_OBJECTS:fluxolibF90> $<TARGET_OBJECTS:userblocklib>)

# ADD_LIBRARY(fluxolib STATIC $<TARGET_OBJECTS:fluxolibF90> )
# ADD_DEPENDENCIES(fluxolibF90 HDF5)
LIST(APPEND fluxoobj $<TARGET_OBJECTS:fluxolibF90>)
IF (fluxoC)
  ADD_LIBRARY(fluxolibC OBJECT ${fluxoC} )
  ##ADD_DEPENDENCIES(fluxolibC HDF5 P4EST)
  SET_TARGET_PROPERTIES(fluxolibC PROPERTIES COMPILE_FLAGS "-fPIC") # TODO: add C default flags
  LIST(APPEND fluxoobj $<TARGET_OBJECTS:fluxolibC>)
ENDIF()
ADD_LIBRARY(fluxolib STATIC ${fluxoobj})
LIST(APPEND fluxodeps fluxolib HDF5)
IF (FLUXO_AMR)
  LIST(APPEND fluxodeps P4EST)
ENDIF()

IF(SELFBUILTEXTERNALS)
  ADD_DEPENDENCIES(fluxolibF90 ${SELFBUILTEXTERNALS})
  ADD_DEPENDENCIES(fluxolib ${SELFBUILTEXTERNALS})
  ADD_DEPENDENCIES(fluxolibC ${SELFBUILTEXTERNALS})
ENDIF()

add_exec(fluxo ./src/fluxo.f90)

ADD_DEPENDENCIES(fluxo fluxolib ${INTERNALLIBS})
##ADD_DEPENDENCIES(fluxo ${fluxodeps} ${INTERNALLIBS})

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

