#include "defines.h"

!==================================================================================================================================
!> This module only initializes the equation specific parameters and computes analytical functions and the evaluation of the source
!==================================================================================================================================
MODULE MOD_Indicators
    
    ! MODULES
    IMPLICIT NONE
    PUBLIC
    SAVE 

abstract interface
  pure subroutine i_sub_GetIndicator(U,ind)
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: ind
  end subroutine i_sub_GetIndicator
end interface

procedure(i_sub_GetIndicator), pointer :: CustomIndicator

PUBLIC :: InitIndicator

CONTAINS

 SUBROUTINE InitIndicator()
    !============================================================================================================================
    ! MODULES
    USE MOD_Globals
    USE MOD_PreProc
    USE MOD_Indicators_vars
    USE MOD_ReadInTools
    USE MOD_Mesh_Vars         ,ONLY: nElems,nSides,firstSlaveSide,LastSlaveSide, isMortarMesh
    USE MOD_Interpolation_Vars,ONLY:xGP,InterpolationInitIsDone
   
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    !----------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    !----------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    !----------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    integer :: whichIndicator
    !============================================================================================================================
    IF (IndicatorsInitIsDone .OR. (.NOT.InterpolationInitIsDone)) THEN
      SWRITE(*,*) "InitIndicator not ready to be called or already called."
      RETURN
    END IF
    IF (PP_N.LT.2) THEN
      CALL abort(__STAMP__,'Polynomial Degree too small for Indicator!',999,999.)
      RETURN
    END IF
    
    ! shock caturing parameters
      
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT SHOCKCAPTURING INDICATOR...'
    CALL InitBasisTrans(PP_N,xGP)
  
#if navierstokes || mhd
    whichIndicator = GETINT('ShockIndicatorAMR','1')
    select case (whichIndicator)
      case(1)
        CustomIndicator => GetDensity
        SWRITE(UNIT_StdOut,'(A)') '    USING DENSITY AS SHOCK INDICATOR!'
      case(2)
        CustomIndicator => GetPressure
        SWRITE(UNIT_StdOut,'(A)') '    USING PRESSURE AS SHOCK INDICATOR!'
      case(3)
        CustomIndicator => GetDensityTimesPressure
        SWRITE(UNIT_StdOut,'(A)') '    USING DENSITY TIMES PRESSURE AS SHOCK INDICATOR!'
      case(4)
        CustomIndicator => GetKinEnergy
        SWRITE(UNIT_StdOut,'(A)') '    USING KINTETIC ENERGY AS SHOCK INDICATOR!'
    end select
#endif /*navierstokes || mhd*/
    if (isMortarMesh) then
      SWRITE(UNIT_stdOut,'(A)')' WARNING: Shock capturing coefficients are not transferred correctly across mortars!'
    end if
    
    
    IndicatorsInitIsDone = .TRUE.
    SWRITE(UNIT_stdOut,'(A)')' INIT SHOCKCAPTURING INDICATOR DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitIndicator
    
    


 SUBROUTINE InitBasisTrans(N_in,xGP)
    !===================================================================================================================================
    !> Initialize Vandermodematrix for basis transformation
    !===================================================================================================================================
    ! MODULES
    USE MOD_Indicators_vars         ,ONLY:sVdm_Leg
    USE MOD_Basis                   ,ONLY :BuildLegendreVdm
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! INPUT VARIABLES
    INTEGER,INTENT(IN)                         :: N_in
    REAL,INTENT(IN),DIMENSION(0:N_in)          :: xGP
    REAL,DIMENSION(0:N_in,0:N_in)              :: Vdm_Leg
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! OUTPUT VARIABLES
    !-----------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    !===================================================================================================================================
    !  NODAL <--> MODAL
    ! Compute the 1D Vandermondematrix, needed to tranform the nodal basis into a modal (Legendre) basis
    ALLOCATE(sVdm_Leg(0:N_in,0:N_in))
    CALL BuildLegendreVdm(N_in,xGP,Vdm_Leg,sVdm_Leg)
    END SUBROUTINE InitBasisTrans
    
!============================================================================================================================
!> Modal shock sensor of Persson and Peraire
!> Persson, P. O.; Peraire, J. (2006). "Sub-cell shock capturing for discontinuous Galerkin methods". In 44th AIAA Aerospace Sciences Meeting and Exhibit (p. 112).
!============================================================================================================================
FUNCTION ShockSensor_PerssonPeraire(U) RESULT (eta)
    USE MOD_PreProc
    use MOD_ChangeBasis                     ,only: ChangeBasis3D
    use MOD_Mesh_Vars                       ,only: nElems
    use MOD_Indicators_vars                 ,only: sVdm_Leg
    implicit none
    ! Arguments
    !---------------------------------------------------------------------------------------------------------------------------------
    real,dimension(PP_nVar,0:PP_N,0:PP_N,0:PP_N), intent(in)  :: U
    real                                                :: eta
    ! Local variables
    !---------------------------------------------------------------------------------------------------------------------------------
    real, dimension(1:1,0:PP_N,0:PP_N,0:PP_N) :: Uind,Umod
    real                                      :: LU,LUM1,LUM2,LU_N,LU_NM1
    integer                                   :: l
    !---------------------------------------------------------------------------------------------------------------------------------
    
 
      ! Get the indicator variable (Uind)
      call GetIndicator_3D(Uind,U(:,:,:,:))
      
      ! Transform Uind into modal Legendre interpolant Umod
      CALL ChangeBasis3D(1,PP_N,PP_N,sVdm_Leg,Uind,Umod)
      
      ! Compute (truncated) error norms
      LU     = SUM(Umod(1,0:PP_N  ,0:PP_N  ,0:PP_N  )**2)
      LUM1   = SUM(Umod(1,0:PP_N-1,0:PP_N-1,0:PP_N-1)**2)
      LUM2   = SUM(Umod(1,0:PP_N-2,0:PP_N-2,0:PP_N-2)**2)
      LU_N   = LU-LUM1
      LU_NM1 = LUM1-LUM2
  
      ! DOF energy indicator
      eta = MAX(LU_N/LU,LU_NM1/LUM1)
      
      if (eta < epsilon(eta)) eta = epsilon(eta)
   
END FUNCTION ShockSensor_PerssonPeraire

  
  !============================================================================================================================
!> Get the shock indicator quantity for all degrees of freedom of an element
!============================================================================================================================
  pure subroutine GetIndicator_3D(Uind,U)
    USE MOD_PreProc
    implicit none
    !-arguments-------------------------------------------
    real, intent(in)  :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N)
    real, intent(out) :: Uind (1:1,0:PP_N,0:PP_N,0:PP_N)
    !-local-variables-------------------------------------
    integer :: i,j,k
    !-----------------------------------------------------
    
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      call CustomIndicator(U(:,i,j,k),Uind(1,i,j,k))
    end do       ; end do       ; end do
  end subroutine GetIndicator_3D
  

  

#if navierstokes || mhd

!============================================================================================================================
!============================================================================================================================
! FOLLOWING ROUTINES COMPUTE PHYSICAL QUANTITIES... THEY COULD BE MOVED TO equation

!============================================================================================================================
!> Get Density
!============================================================================================================================
pure subroutine GetDensity(U,rho)
implicit none
real, intent(in)  :: U(PP_nVar)
real, intent(out) :: rho

rho = U(1)

end subroutine GetDensity
!============================================================================================================================
!> Get Pressure
!============================================================================================================================
pure subroutine GetPressure(U,p)
#ifdef mhd
use MOD_Equation_Vars      , only: s2mu_0
#endif /*mhd*/
use MOD_Equation_Vars      , only: KappaM1
implicit none
real, intent(in)  :: U(PP_nVar)
real, intent(out) :: p

p = KappaM1*(U(5)-0.5*(SUM(U(2:4)*U(2:4))/U(1)))
#ifdef mhd
p = p - KappaM1*s2mu_0*SUM(U(6:8)*U(6:8))
#ifdef PP_GLM
p = p - 0.5*KappaM1*U(9)*U(9)
#endif /*PP_GLM*/
#endif /*mhd*/
end subroutine GetPressure
!============================================================================================================================
!> Get Density Times Pressure
!============================================================================================================================
pure subroutine GetDensityTimesPressure(U,rhop)
implicit none
real, intent(in)  :: U(PP_nVar)
real, intent(out) :: rhop
!-----------------------------------------
real :: p
!-----------------------------------------

call GetPressure(U,p)
rhop = U(1) * p

end subroutine GetDensityTimesPressure
!============================================================================================================================
!> Get Density Times Pressure
!============================================================================================================================
pure subroutine GetKinEnergy(U,kinen)
#ifdef mhd
use MOD_Equation_Vars      , only: s2mu_0
#endif /*mhd*/
implicit none
real, intent(in)  :: U(PP_nVar)
real, intent(out) :: kinen

kinen = SUM(U(2:4)*U(2:4))/U(1)    
#ifdef mhd
kinen = kinen + s2mu_0*SUM(U(6:8)*U(6:8))
#endif /*mhd*/
end subroutine GetKinEnergy
#endif /*navierstokes || mhd*/

END MODULE MOD_Indicators
