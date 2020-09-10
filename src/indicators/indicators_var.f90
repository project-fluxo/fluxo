#include "defines.h"

!==================================================================================================================================
!> This module only initializes the equation specific parameters and computes analytical functions and the evaluation of the source
!==================================================================================================================================
MODULE MOD_Indicators_vars
    
    ! MODULES
    IMPLICIT NONE
    ! PRIVATE
    PUBLIC
    SAVE

    ! > Dummy interface for Shock-Indicator function pointer
    ! ABSTRACT INTERFACE
    !   FUNCTION ShockCaptureIndicator(U)
    !     REAL, INTENT(IN) :: U(:,:,:,:)
       
    !   END FUNCTION
    ! END INTERFACE

    
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
    LOGICAL          :: IndicatorsInitIsDone=.FALSE.     !< marks whether the shock capturing init routine is complete
    REAL,ALLOCATABLE :: sVdm_Leg(:,:)                        !< 1D inverse Vandermondematrix to Legendre polynomials
    
    ! PROCEDURE(ShockCaptureIndicator),POINTER :: Indicator !< pointer to shock capture routine, depends on U

 CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
END MODULE MOD_Indicators_vars