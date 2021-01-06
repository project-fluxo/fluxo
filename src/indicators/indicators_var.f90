#include "defines.h"

!==================================================================================================================================
!> This module only initializes the equation specific parameters and computes analytical functions and the evaluation of the source
!==================================================================================================================================
MODULE MOD_Indicators_vars
    IMPLICIT NONE
    ! PRIVATE
    PUBLIC
    SAVE
    
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
    REAL,ALLOCATABLE :: sVdm_Leg(:,:)                        !< 1D inverse Vandermondematrix to Legendre polynomials
    integer          :: NumberOfPerssonInd =0

 CONTAINS

END MODULE MOD_Indicators_vars
