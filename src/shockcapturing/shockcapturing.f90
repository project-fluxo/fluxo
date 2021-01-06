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

PUBLIC :: DefineParametersShockCapturing
PUBLIC :: InitShockCapturing
#if SHOCK_NFVSE
public :: CalcBlendingCoefficient
#endif /*SHOCK_NFVSE*/
PUBLIC :: FinalizeShockCapturing
public :: InitShockCapturingAfterAdapt
!===================================================================================================================================
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersShockCapturing()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
use MOD_NFVSE       ,only: DefineParametersNFVSE
use MOD_Equation_Vars, only: IndicatorQuantityNames, nIndVar
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
integer            :: i
character(len=255) :: IndicatorQuantities,fmt
!==================================================================================================================================
CALL prms%SetSection("ShockCapturing")

CALL prms%CreateIntOption(     "ShockIndicator",  " Specifies the quantity to be used as shock-indicator:\n"//&
                                              "   1: Persson-Peraire indicator\n"//&
                                              "   2: Randomly assign the blending coef.\n"//&
                                              "   3: Fixed blending coef. (alpha=ShockBlendCoef)"&
                                             ,"3")

write(fmt,'(A,I0,A)') '(',nIndVar,'A)'
write(IndicatorQuantities,fmt) ('  * '//trim(IndicatorQuantityNames(i))//'\n', i=1, nIndVar)

CALL prms%CreateStringOption("ShockIndicatorQuantity",  " Specifies the quantity to be used for the shock indicator. One of the following:\n"//&
                                              trim(IndicatorQuantities)&
                                             ,"DensityTimesPressure")

CALL prms%CreateRealOption(   "ShockBlendCoef",  " Fixed blending coefficient to be used with ShockIndicator=3", "0.0")

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
alpha        = 0.0
alpha_Master = 0.0
alpha_Slave  = 0.0

select case (ModalThreshold)
  case(1) ; threshold = 0.5 * 10.0 ** (-1.8 * (PP_N+1.)**0.25) ! Sebastian's thresold (Euler and MHD paper)
  case(2) ; threshold = 0.5 * 10.0 ** (-1.8 * PP_N**0.25)      ! New threshold
end select

call InitNFVSE()

#endif /*SHOCK_NFVSE*/

CALL InitBasisTrans(PP_N,xGP)

ShockIndicator         = GETINT('ShockIndicator','1')
ShockIndicatorQuantity = GETSTR('ShockIndicatorQuantity','DensityTimesPressure')
ShockBlendCoef         = GETREAL('ShockBlendCoef','0.0')

if (ShockIndicator==1) then
  call Shock_Indicator % construct(ShockIndicatorQuantity)
elseif (ShockIndicator>3) then
  CALL abort(__STAMP__,'Shock indicator not defined!',999,999.)
  RETURN
end if

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
subroutine CalcBlendingCoefficient(U)
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
! -------------------------------------------------------------------
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
  
! Compute the blending coefficients
! ---------------------------------
  select case (ShockIndicator)
    case(1)   ! Persson-Peraire indicator
      ! Shock indicator
      eta = Shock_Indicator % compute(U)
      
      ! Compute and correct alpha
      do eID=1, nElems
        alpha(eID) = max(TimeRelFactor*alpha(eID), 1.0 / (1.0 + exp(-sharpness * (eta(eID) - threshold)/threshold )) )
      end do
      
    case(2)   ! Random indicator
      do eID=1, nElems
        call RANDOM_NUMBER(alpha(eID))
      end do
      
    case(3)   ! Fixed blending coefficients
      alpha = ShockBlendCoef
  end select
  
! Impose alpha_max and alpha_min
! ------------------------------
  where (alpha < alpha_min)
    alpha = 0.0
  elsewhere (alpha >= alpha_max)
    alpha = alpha_max
  end where
  
  
! Start first space propagation sweep (MPI-optimized)
! ---------------------------------------------------
  if (SpacePropSweeps > 0) call ProlongBlendingCoeffToFaces()
  
end subroutine CalcBlendingCoefficient
#endif /*SHOCK_NFVSE*/
  
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

call Shock_Indicator % destruct
#if SHOCK_NFVSE
call FinalizeNFVSE()
#endif /*SHOCK_NFVSE*/

END SUBROUTINE FinalizeShockCapturing

#endif /*SHOCKCAPTURE*/
END MODULE MOD_ShockCapturing
