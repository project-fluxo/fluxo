!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2020 Andrés Rueda
! Copyright (c) 2010 - 2016 Claus-Dieter Munz (github.com/flexi-framework/flexi)
!
! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 
! of the License, or (at your option) any later version.
!
! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
!==================================================================================================================================
#include "defines.h"

!==================================================================================================================================
!> Module for the shock capturing routines
!==================================================================================================================================
MODULE MOD_ShockCapturing
#if SHOCKCAPTURE
! MODULES
IMPLICIT NONE
PRIVATE
! ----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersShockCapturing
   MODULE PROCEDURE DefineParametersShockCapturing
END INTERFACE

INTERFACE InitShockCapturing
   MODULE PROCEDURE InitShockCapturing
END INTERFACE

INTERFACE FinalizeShockCapturing
   MODULE PROCEDURE FinalizeShockCapturing
END INTERFACE

abstract interface
  ! Interface for the indicator
  pure subroutine i_sub_GetIndicator(U,ind)
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: ind
  end subroutine i_sub_GetIndicator
  ! Interface for the Blending coefficient subroutine
  subroutine i_sub_CalcBlendingCoefficient(U)
    use MOD_PreProc
    use MOD_Mesh_Vars          , only: nElems
    real, intent(in)  :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
  end subroutine i_sub_CalcBlendingCoefficient
end interface

! Routines that can be pointed
procedure(i_sub_GetIndicator)           , pointer :: CustomIndicator => null()
procedure(i_sub_CalcBlendingCoefficient), pointer :: CalcBlendingCoefficient => null()

PUBLIC :: DefineParametersShockCapturing
PUBLIC :: InitShockCapturing
#if SHOCK_NFVSE
public :: CalcBlendingCoefficient
#endif /*SHOCK_NFVSE*/
PUBLIC :: FinalizeShockCapturing
public :: GetPressure ! TODO: Move this to equation
!==================================================================================================================================
! local definitions for inlining / optimizing routines
#if PP_Indicator_Var==0
#  define PP_GetIndicator CustomIndicator
#elif PP_Indicator_Var==1
#  define PP_GetIndicator GetDensity
#elif PP_Indicator_Var==2
#  define PP_GetIndicator GetPressure
#elif PP_Indicator_Var==3
#  define PP_GetIndicator GetDensityTimesPressure
#elif PP_Indicator_Var==4
#  define PP_GetIndicator GetKinEnergy
#endif
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersShockCapturing()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("ShockCapturing")

CALL prms%CreateIntOption(     "ShockIndicator",  " Specifies the quantity to be used as shock-indicator "//&
                                              "  1: Density"//&
                                              "  2: Pressure"//&
                                              "  3: Density times Pressure"//&
                                              "  4: Kinetic Energy"//&
                                              "  5: Randomly assign the blending coefficients"&
                                             ,"3")

END SUBROUTINE DefineParametersShockCapturing

SUBROUTINE InitShockCapturing()
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ShockCapturing_Vars
USE MOD_ReadInTools
USE MOD_Mesh_Vars         ,ONLY: nElems,firstSlaveSide,LastSlaveSide, isMortarMesh
USE MOD_Interpolation_Vars,ONLY: xGP,InterpolationInitIsDone
#if SHOCK_NFVSE
use MOD_NFVSE             , only: InitNFVSE
#endif /*SHOCK_NFVSE*/
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
IF (ShockCapturingInitIsDone.OR.(.NOT.InterpolationInitIsDone)) THEN
  SWRITE(*,*) "InitShockCapturing not ready to be called or already called."
  RETURN
END IF
IF (PP_N.LT.2) THEN
  CALL abort(__STAMP__,'Polynomial Degree too small for Shock Capturing!',999,999.)
  RETURN
END IF

nu_max = 0.

#if SHOCK_NFVSE
allocate ( alpha(nElems) )
allocate ( alpha_Master(firstSlaveSide:LastSlaveSide) ) ! Only allocating on slave sides (no BCs needed, and mortars not considered yet -TODO!)
allocate ( alpha_Slave (firstSlaveSide:LastSlaveSide) )
alpha        = 0.d0
alpha_Master = 0.d0
alpha_Slave  = 0.d0


!~ threshold = 0.5d0 * 10.d0 ** (-1.8d0 * (PP_N + 1)**4)      ! Old Sebastian: This is almost cero ALWAYS
!~ threshold = 0.5d0 * (PP_N + 1)**(-1.8d0 * 4)               ! Mine
threshold = 0.5d0 * 10.d0 ** (-1.8d0 * (PP_N + 1.d0)**0.25d0) ! New Sebastian

call InitNFVSE()
#endif /*SHOCK_NFVSE*/

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SHOCKCAPTURING...'
CALL InitBasisTrans(PP_N,xGP)

#if PP_Indicator_Var==0
whichIndicator = GETINT('ShockIndicator','3')
#else
whichIndicator = PP_Indicator_Var
#endif

CalcBlendingCoefficient => CalcBlendingCoefficient_indicator ! Default
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
  case(5)
    CalcBlendingCoefficient => CalcBlendingCoefficient_random ! Override CalcBlendingCoefficient
    SWRITE(UNIT_StdOut,'(A)') '    USING RANDOM BLENDING COEFFICIENTS!'
end select

if (isMortarMesh) then
  SWRITE(UNIT_stdOut,'(A)')' WARNING: Shock capturing coefficients are not transferred correctly across mortars!'
end if

ShockCapturingInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SHOCKCAPTURING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitShockCapturing


SUBROUTINE InitBasisTrans(N_in,xGP)
!===================================================================================================================================
!> Initialize Vandermodematrix for basis transformation
!===================================================================================================================================
! MODULES
USE MOD_ShockCapturing_Vars,ONLY:sVdm_Leg
#if SHOCK_LOC_ARTVISC
USE MOD_ShockCapturing_Vars,ONLY:FilterMat
#endif /*SHOCK_LOC_ARTVISC*/
USE MOD_Basis, ONLY :BuildLegendreVdm
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

#if SHOCK_NFVSE
!===================================================================================================================================
!> Routines to compute the blending coefficient for NFVSE
!> -> See Hennemann and Gassner (2020). "Entropy stable shock capturing for the discontinuous galerkin spectral element
!>                                          method with native finite volume sub elements"
!> -> This routine computes the sensor, makes the correction (with alpha_min and alpha_max), and sends the information with MPI
!> -> No propagation is done yet (MPI informationmust be received).
!===================================================================================================================================
subroutine CalcBlendingCoefficient_indicator(U)
  use MOD_PreProc
  use MOD_ShockCapturing_Vars
  use MOD_Mesh_Vars          , only: nElems
  use MOD_NFVSE_MPI          , only: ProlongBlendingCoeffToFaces
  implicit none
  ! Arguments
  !---------------------------------------------------------------------------------------------------------------------------------
  real,dimension(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems), intent(in)  :: U
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Local variables
  real ::  eta(nElems)
  !---------------------------------------------------------------------------------------------------------------------------------
  
  ! Shock indicator
  call ShockSensor_PerssonPeraire(U,eta)
  
  ! Compute and correct alpha
  alpha = 1.0 / (1.0 + exp(-sharpness * (eta - threshold)/threshold ))
  
  where (alpha < alpha_min)
    alpha = 0.0
  elsewhere (alpha >= alpha_max)
    alpha = alpha_max
  end where
  
  call ProlongBlendingCoeffToFaces()
  
end subroutine CalcBlendingCoefficient_indicator
!===================================================================================================================================
!> This routine selects the blending coefficient randomly (and sends it with MPI)
!===================================================================================================================================
subroutine CalcBlendingCoefficient_random(U)
  use MOD_PreProc
  use MOD_ShockCapturing_Vars
  use MOD_Mesh_Vars          , only: nElems
  use MOD_NFVSE_MPI          , only: ProlongBlendingCoeffToFaces
  implicit none
  ! Arguments
  !---------------------------------------------------------------------------------------------------------------------------------
  real,dimension(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems), intent(in)  :: U
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Local variables
  real ::  eta(nElems)
  integer :: eID
  !---------------------------------------------------------------------------------------------------------------------------------
  
  do eID=1, nElems
    call RANDOM_NUMBER(alpha(eID))
  end do
  
  call ProlongBlendingCoeffToFaces()
end subroutine CalcBlendingCoefficient_random
#endif /*SHOCK_NFVSE*/

!============================================================================================================================
!> Modal shock sensor of Persson and Peraire
!> Persson, P. O.; Peraire, J. (2006). "Sub-cell shock capturing for discontinuous Galerkin methods". In 44th AIAA Aerospace Sciences Meeting and Exhibit (p. 112).
!============================================================================================================================
subroutine ShockSensor_PerssonPeraire(U,eta)
  USE MOD_PreProc
  use MOD_ChangeBasis        , only: ChangeBasis3D
  use MOD_Mesh_Vars          , only: nElems
  use MOD_ShockCapturing_Vars, only: sVdm_Leg
  implicit none
  ! Arguments
  !---------------------------------------------------------------------------------------------------------------------------------
  real,dimension(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems), intent(in)  :: U
  real                                               , intent(out) :: eta(nElems)
  ! Local variables
  !---------------------------------------------------------------------------------------------------------------------------------
  real, dimension(1:1,0:PP_N,0:PP_N,0:PP_N) :: Uind,Umod
  real                                      :: LU,LUM1,LUM2,LU_N,LU_NM1
  integer                                   :: l
  !---------------------------------------------------------------------------------------------------------------------------------
  
  do l=1, nElems
    ! Get the indicator variable (Uind)
    call GetIndicator_3D(Uind,U(:,:,:,:,l))
    
    ! Transform Uind into modal Legendre interpolant Umod
    CALL ChangeBasis3D(1,PP_N,PP_N,sVdm_Leg,Uind,Umod)
    
    ! Compute (truncated) error norms
    LU     = SUM(Umod(1,0:PP_N  ,0:PP_N  ,0:PP_N  )**2)
    LUM1   = SUM(Umod(1,0:PP_N-1,0:PP_N-1,0:PP_N-1)**2)
    LUM2   = SUM(Umod(1,0:PP_N-2,0:PP_N-2,0:PP_N-2)**2)
    LU_N   = LU-LUM1
    LU_NM1 = LUM1-LUM2

    ! DOF energy indicator
    eta(l) = MAX(LU_N/LU,LU_NM1/LUM1)
    
  end do !l
end subroutine ShockSensor_PerssonPeraire
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
      call PP_GetIndicator(U(:,i,j,k),Uind(1,i,j,k))
    end do       ; end do       ; end do
  end subroutine GetIndicator_3D
  
SUBROUTINE FinalizeShockCapturing()
!============================================================================================================================
!> Deallocate all global shock capturing variables.
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_ShockCapturing_Vars
#if SHOCK_NFVSE
use MOD_NFVSE             , only: FinalizeNFVSE
#endif /*SHOCK_NFVSE*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!============================================================================================================================
IF (.NOT.ShockCapturingInitIsDone) THEN
  WRITE(UNIT_stdOut,*) "InitShockCapturing was not called before."
  RETURN
END IF
ShockCapturingInitIsDone = .FALSE.
SDEALLOCATE(sVdm_Leg)
SDEALLOCATE(alpha)

#if SHOCK_NFVSE
call FinalizeNFVSE()
#endif /*SHOCK_NFVSE*/

END SUBROUTINE FinalizeShockCapturing
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

#endif /*SHOCKCAPTURE*/
END MODULE MOD_ShockCapturing
