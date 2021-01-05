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
public :: InitShockCapturingAfterAdapt
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
use MOD_NFVSE       ,only: DefineParametersNFVSE
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

CALL prms%CreateIntOption(     "ModalThreshold",  " Threshold to be used for the indicator "//&
                                              "  1: 0.5 * 10.0 ** (-1.8 * (PP_N+1)**0.25)"//&
                                              "  2: 0.5 * 10.0 ** (-1.8 * PP_N**0.25)"&
                                             ,"1")

#if NFVSE_CORR
CALL prms%CreateRealOption(   "PositCorrFactor",  " The correction factor for NFVSE", "0.1")
CALL prms%CreateIntOption(       "PositMaxIter",  " Maximum number of iterations for positivity limiter", "10")
#endif /*NFVSE_CORR*/

call DefineParametersNFVSE()

END SUBROUTINE DefineParametersShockCapturing

SUBROUTINE InitShockCapturing()
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ShockCapturing_Vars
USE MOD_ReadInTools
USE MOD_Mesh_Vars         ,ONLY: nElems,nSides,firstSlaveSide,LastSlaveSide, isMortarMesh, firstMortarInnerSide
USE MOD_Interpolation_Vars,ONLY:xGP,InterpolationInitIsDone
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
integer :: i,j,k    ! DOF counters
integer :: eID      ! Element counter
real    :: conMetrics(3,3)  ! Contravariant metric tensor in each element
!============================================================================================================================
SDEALLOCATE(alpha)
SDEALLOCATE(alpha_Master)
SDEALLOCATE(alpha_Slave)

IF (ShockCapturingInitIsDone.OR.(.NOT.InterpolationInitIsDone)) THEN
  SWRITE(*,*) "InitShockCapturing not ready to be called or already called."
  RETURN
END IF
IF (PP_N.LT.2) THEN
  CALL abort(__STAMP__,'Polynomial Degree too small for Shock Capturing!',999,999.)
  RETURN
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT SHOCKCAPTURING...'

ModalThreshold = GETINT('ModalThreshold',"1")

#if SHOCK_NFVSE
allocate ( alpha(nElems) )
allocate ( alpha_Master(firstMortarInnerSide:nSides ) )
allocate ( alpha_Slave (firstSlaveSide:LastSlaveSide) )
alpha        = 0.d0
alpha_Master = 0.d0
alpha_Slave  = 0.d0

select case (ModalThreshold)
  case(1) ; threshold = 0.5d0 * 10.d0 ** (-1.8d0 * (PP_N+1.)**0.25d0) ! Sebastian's thresold (Euler paper)
  case(2) ; threshold = 0.5 * 10.0 ** (-1.8 * PP_N**0.25) ! New threshold (MHD paper)
end select

call InitNFVSE()

#endif /*SHOCK_NFVSE*/

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

#if NFVSE_CORR
PositCorrFactor = GETREAL('PositCorrFactor','0.1')
PositMaxIter = GETINT('PositMaxIter','10')
SWRITE(UNIT_stdOut,'(A,ES16.7)') '    *NFVSE correction activated with PositCorrFactor=', PositCorrFactor
#endif /*NFVSE_CORR*/

ShockCapturingInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT SHOCKCAPTURING DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitShockCapturing
!===================================================================================================================================
!> Reinitializes all variables that need reinitialization after the h-adaptation
!> ATTENTION: At the moment only for NFVSE
!===================================================================================================================================
SUBROUTINE InitShockCapturingAfterAdapt(ChangeElem,nElemsOld,nSidesOld,firstSlaveSideOld,LastSlaveSideOld,firstMortarInnerSideOld)
  use MOD_NFVSE               , only: InitNFVSEAfterAdaptation
  use MOD_NFVSE_Vars          , only: TimeRelFactor
  USE MOD_ShockCapturing_Vars
  USE MOD_ReadInTools
  USE MOD_Mesh_Vars         ,ONLY: nElems,nSides,firstSlaveSide,LastSlaveSide,firstMortarInnerSide
  implicit none
  !-arguments-----------------------------------
  integer, intent(in) :: ChangeElem(8,nElems)
  integer, intent(in) :: nElemsOld,nSidesOld,firstSlaveSideOld,LastSlaveSideOld,firstMortarInnerSideOld
  !-local-variables-----------------------------
  integer          :: eID
  real,allocatable,target :: alphaNew(:)
  !---------------------------------------------
  
#if SHOCK_NFVSE
  
! Reallocate storage
! ------------------
  
  if ( (firstMortarInnerSideOld .ne. firstMortarInnerSide) .or. (nSides .ne. nSidesOld) ) then
    SDEALLOCATE(alpha_Master)
    allocate ( alpha_Master(firstMortarInnerSide:nSides ) )
  end if
  if ( (firstSlaveSide .ne. firstSlaveSideOld) .or. (LastSlaveSide .ne. LastSlaveSideOld) ) then
    SDEALLOCATE(alpha_Slave)
    allocate ( alpha_Slave (firstSlaveSide:LastSlaveSide) )
  end if
  
! Initialize values
! -----------------
  alpha_Master = 0.0
  alpha_Slave  = 0.0
  
  if (TimeRelFactor < alpha_min/alpha_max) then
    ! The time relaxation has no effect, alpha can be set to 0
    if (nElems /= nElemsOld) then
      SDEALLOCATE(alpha)
      allocate(alpha(nElems))
    end if
    alpha = 0.0
  else
    allocate ( alphaNew(nElems) )
    ! Set with old values
    do eID=1, nElems
      if (ChangeElem(1,eID) < 0) then
        ! refinement
        alphaNew(eID) = alpha(-ChangeElem(1,eID))
      elseif (ChangeElem(2,eID) > 0) then
        ! coarsening
        alphaNew(eID) = maxval(alpha(ChangeElem(1:8,eID)))
      else
        ! simple reasignment
        alphaNew(eID) = alpha(ChangeElem(1,eID))
      endif
    end do
    call move_alloc(alphaNew,alpha)
  end if
  
  call InitNFVSEAfterAdaptation(ChangeElem,nElemsOld)
#endif /*SHOCK_NFVSE*/
END SUBROUTINE InitShockCapturingAfterAdapt

SUBROUTINE InitBasisTrans(N_in,xGP)
!===================================================================================================================================
!> Initialize Vandermodematrix for basis transformation
!===================================================================================================================================
! MODULES
USE MOD_ShockCapturing_Vars,ONLY:sVdm_Leg
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
real,DIMENSION(0:N_in,0:N_in)              :: Filter
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
  use MOD_NFVSE_MPI          , only: ProlongBlendingCoeffToFaces, PropagateBlendingCoeff
  use MOD_NFVSE_Vars         , only: SpacePropSweeps, TimeRelFactor, RECONS_NEIGHBOR
  ! For reconstruction on boundaries
#if MPI
  use MOD_Mesh_Vars          , only: firstSlaveSide, lastSlaveSide
  use MOD_NFVSE_Vars         , only: ReconsBoundaries, MPIRequest_Umaster
  use MOD_MPI                , only: StartReceiveMPIData,StartSendMPIData
  USE MOD_MPI_Vars
  use MOD_DG_Vars            , only: U_master
#endif /*MPI*/
  implicit none
  ! Arguments
  !---------------------------------------------------------------------------------------------------------------------------------
  real,dimension(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems), intent(in)  :: U
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Local variables
  real ::  eta(nElems)
  integer :: eID, sweep
  !---------------------------------------------------------------------------------------------------------------------------------
  
  ! If we do reconstruction on boundaries, we need to send the U_master
#if MPI
  if (ReconsBoundaries >= RECONS_NEIGHBOR) then
    ! receive the master
    call StartReceiveMPIData(U_master(:,:,:,firstSlaveSide:lastSlaveSide), DataSizeSide, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_Umaster(:,1), SendID=1) ! Receive YOUR  (sendID=1) 
    
    ! Send the master
    call StartSendMPIData   (U_master(:,:,:,firstSlaveSide:lastSlaveSide), DataSizeSide, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_Umaster(:,2),SendID=1) 
  end if
#endif /*MPI*/
  
  ! Shock indicator
  call ShockSensor_PerssonPeraire(U,eta)
  
  ! Compute and correct alpha
  do eID=1, nElems
    alpha(eID) = max(TimeRelFactor*alpha(eID), 1.0 / (1.0 + exp(-sharpness * (eta(eID) - threshold)/threshold )) )
  end do

  where (alpha < alpha_min)
    alpha = 0.0
  elsewhere (alpha >= alpha_max)
    alpha = alpha_max
  end where
  
  ! Start first space propagation sweep (MPI-optimized)
  if (SpacePropSweeps > 0) call ProlongBlendingCoeffToFaces()
  
end subroutine CalcBlendingCoefficient_indicator
!===================================================================================================================================
!> This routine selects the blending coefficient randomly (and sends it with MPI)
!===================================================================================================================================
subroutine CalcBlendingCoefficient_random(U)
  use MOD_PreProc
  use MOD_ShockCapturing_Vars
  use MOD_Mesh_Vars          , only: nElems
  use MOD_NFVSE_MPI          , only: ProlongBlendingCoeffToFaces, PropagateBlendingCoeff
  use MOD_NFVSE_Vars         , only: SpacePropSweeps, TimeRelFactor
#if MPI
  use MOD_Mesh_Vars          , only: firstSlaveSide, lastSlaveSide
  use MOD_NFVSE_Vars         , only: ReconsBoundaries, MPIRequest_Umaster, RECONS_NEIGHBOR
  use MOD_MPI                , only: StartReceiveMPIData,StartSendMPIData
  USE MOD_MPI_Vars
  use MOD_DG_Vars            , only: U_master
#endif /*MPI*/
  implicit none
  ! Arguments
  !---------------------------------------------------------------------------------------------------------------------------------
  real,dimension(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems), intent(in)  :: U
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Local variables
  real ::  eta(nElems)
  integer :: eID, sweep
  !---------------------------------------------------------------------------------------------------------------------------------
  
  ! If we do reconstruction on boundaries, we need to send the U_master
#if MPI
  if (ReconsBoundaries >= RECONS_NEIGHBOR) then
    ! receive the master
    call StartReceiveMPIData(U_master(:,:,:,firstSlaveSide:lastSlaveSide), DataSizeSide, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_Umaster(:,1), SendID=1) ! Receive YOUR  (sendID=1) 
    
    ! Send the master
    call StartSendMPIData   (U_master(:,:,:,firstSlaveSide:lastSlaveSide), DataSizeSide, firstSlaveSide, lastSlaveSide, &
                             MPIRequest_Umaster(:,2),SendID=1) 
  end if
#endif /*MPI*/
  
  do eID=1, nElems
    call RANDOM_NUMBER(alpha(eID))
  end do
  
  if (SpacePropSweeps > 0) then
    ! Do first sweep (MPI-optimized)
    call ProlongBlendingCoeffToFaces()
    
    ! Do remaining sweeps (overhead in MPI communication)
    do sweep=2, SpacePropSweeps
      call PropagateBlendingCoeff()
      call ProlongBlendingCoeffToFaces()
    end do
  end if
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
