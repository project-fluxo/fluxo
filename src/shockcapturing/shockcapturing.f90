!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
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

#if SHOCK_ARTVISC
INTERFACE CalcArtificialViscosity
   MODULE PROCEDURE CalcArtificialViscosity
END INTERFACE
#endif /*SHOCK_ARTVISC*/

INTERFACE FinalizeShockCapturing
   MODULE PROCEDURE FinalizeShockCapturing
END INTERFACE

abstract interface
  pure subroutine i_sub_GetIndicator(U,ind)
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: ind
  end subroutine i_sub_GetIndicator
end interface

procedure(i_sub_GetIndicator), pointer :: CustomIndicator

PUBLIC :: DefineParametersShockCapturing
PUBLIC :: InitShockCapturing
#if SHOCK_ARTVISC
PUBLIC :: CalcArtificialViscosity
#endif /*SHOCK_ARTVISC*/
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
                                              "  4: Kinetic Energy"&
                                             ,"3")
CALL prms%CreateRealOption(     "ShockCorrFactor",  " The correction factor for NFVSE")
END SUBROUTINE DefineParametersShockCapturing

SUBROUTINE InitShockCapturing()
!============================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_ShockCapturing_Vars
USE MOD_ReadInTools
USE MOD_Mesh_Vars         ,ONLY: nElems,nSides,firstSlaveSide,LastSlaveSide, isMortarMesh
USE MOD_Interpolation_Vars,ONLY:xGP,InterpolationInitIsDone
#if SHOCK_NFVSE
use MOD_NFVSE             , only: InitNFVSE
#endif /*SHOCK_NFVSE*/
#if SHOCK_LOC_ARTVISC
use MOD_Mesh_Vars         , only: Metrics_fTilde, Metrics_gTilde, Metrics_hTilde, sJ
use MOD_Basis             , only: INV33
use MOD_Sensors           , only: SENS_NUM
#endif /*SHOCK_LOC_ARTVISC*/
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
real    :: Mh(3,3)  ! Metric tensor in each element
!============================================================================================================================
IF (ShockCapturingInitIsDone.OR.(.NOT.InterpolationInitIsDone)) THEN
  SWRITE(*,*) "InitShockCapturing not ready to be called or already called."
  RETURN
END IF
IF (PP_N.LT.2) THEN
  CALL abort(__STAMP__,'Polynomial Degree too small for Shock Capturing!',999,999.)
  RETURN
END IF

! shock caturing parameters
#if SHOCK_ARTVISC
ALLOCATE(nu(nElems),nu_Master(nSides),nu_Slave(firstSlaveSide:LastSlaveSide))
nu     = 0.
nu_max = 0.
nu_Master = 0.
nu_Slave = 0.
#endif /*SHOCK_ARTVISC*/

#if SHOCK_LOC_ARTVISC
! Allocate vectors
allocate ( artVisc(SENS_NUM,0:PP_N,0:PP_N,0:PP_N,nElems) )
allocate ( Mh_inv      (3,3,0:PP_N,0:PP_N,0:PP_N,nElems) )

! Inverse of metric tensor
do eID=1, nElems
  do k=0, PP_N  ; do j=0, PP_N  ; do i=0, PP_N
    Mh(:,1) = Metrics_fTilde(:,i,j,k,eID) * sJ(i,j,k,eID)
    Mh(:,2) = Metrics_gTilde(:,i,j,k,eID) * sJ(i,j,k,eID)
    Mh(:,3) = Metrics_hTilde(:,i,j,k,eID) * sJ(i,j,k,eID)
    call INV33(Mh,Mh_inv(:,:,i,j,k,eID))
  end do        ; end do        ; end do ! ijk
end do !eID
#endif /*SHOCK_LOC_ARTVISC*/

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

if (isMortarMesh) then
  SWRITE(UNIT_stdOut,'(A)')' WARNING: Shock capturing coefficients are not transferred correctly across mortars!'
end if

#if NFVSE_CORR
beta = GETREAL('ShockCorrFactor','0.1')
SWRITE(UNIT_stdOut,'(A,ES16.7)') '    *NFVSE correction activated with beta=', beta
#endif /*NFVSE_CORR*/

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

#if SHOCK_ARTVISC
SUBROUTINE CalcArtificialViscosity(U)
!===================================================================================================================================
!> Use framework of Persson and Peraire to measure shocks with DOF energy indicator and calculate artificial viscosity, if necessary
!===================================================================================================================================
! MODULES
USE MOD_PreProc
USE MOD_ShockCapturing_Vars, ONLY: nu,nu_max
USE MOD_Mesh_Vars          , ONLY: nElems,sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars      , ONLY: ConsToPrim
USE MOD_Equation_Vars      , ONLY: FastestWave3D
#if SHOCK_LOC_ARTVISC
use MOD_Sensors            , only: SensorsByFernandezEtAl
USE MOD_ShockCapturing_Vars, ONLY: artVisc, Mh_inv
use MOD_Lifting_Vars       , only: gradPx, gradPy, gradPz
#endif /*SHOCK_LOC_ARTVISC*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems),INTENT(IN) :: U
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                     :: eta_dof(nElems),eta_min,eta_max,eps0
REAL                                     :: v(3),Prim(1:PP_nVar),cf,Max_Lambda(6),h,lambda_max2
INTEGER                                  :: l,i,j,k
!===================================================================================================================================

nu    =0.
nu_max=0.

call ShockSensor_PerssonPeraire(U,eta_dof)

!<temp Compute the physics based sensors to debug
#if SHOCK_LOC_ARTVISC
do l=1, nElems
  do k=0, PP_N  ; do j=0, PP_N  ; do i=0, PP_N
    call SensorsByFernandezEtAl ( artVisc(:,i,j,k,l), &
                                        U(:,i,j,k,l), &
                                   gradPx(:,i,j,k,l), &
                                   gradPy(:,i,j,k,l), &
                                   gradPz(:,i,j,k,l), &
                                 Mh_inv(:,:,i,j,k,l) )
  end do        ; end do        ; end do ! ijk
end do !l
#endif /*SHOCK_LOC_ARTVISC*/
!temp>

eta_dof = log10(eta_dof)
DO l=1,nElems
  ! Artificial Viscosity
  eta_min = -9.0
  eta_max = -3.0
  eps0 = 0.01

  IF (eta_dof(l).GE.eta_max) THEN
    nu(l) = eps0
  ELSE IF (eta_dof(l).LE.eta_min) THEN
    nu(l) = 0.
  ELSE
    nu(l) = 0.5*eps0*(1.0+SIN(PP_Pi*(eta_dof(l)-0.5*(eta_max+eta_min))/(eta_max-eta_min)))
  END IF
  
  ! Save max artificial viscosity for DFL timestepping
  nu_max = MAX(nu_max,nu(l))

  ! Get (transformed) max eigenvalue:
  Max_Lambda = 0.0
  DO i=0,PP_N
    DO j=0,PP_N
      DO k=0,PP_N
        CALL ConsToPrim(Prim,U(:,i,j,k,l))
        CALL FastestWave3D(Prim,cf)
        v(:)=Prim(2:4) 
        Max_Lambda(1)=MAX(Max_Lambda(1),sJ(i,j,k,l)*(ABS(SUM(Metrics_fTilde(:,i,j,k,l)*v)) + &
                        cf*SQRT(SUM(Metrics_fTilde(:,i,j,k,l)*Metrics_fTilde(:,i,j,k,l)))))
        Max_Lambda(2)=MAX(Max_Lambda(2),sJ(i,j,k,l)*(ABS(SUM(Metrics_gTilde(:,i,j,k,l)*v)) + &
                        cf*SQRT(SUM(Metrics_gTilde(:,i,j,k,l)*Metrics_gTilde(:,i,j,k,l)))))
        Max_Lambda(3)=MAX(Max_Lambda(3),sJ(i,j,k,l)*(ABS(SUM(Metrics_hTilde(:,i,j,k,l)*v)) + &
                        cf*SQRT(SUM(Metrics_hTilde(:,i,j,k,l)*Metrics_hTilde(:,i,j,k,l)))))
        Max_Lambda(4)=MAX(Max_Lambda(4),ABS(v(1)) + cf)
        Max_Lambda(5)=MAX(Max_Lambda(5),ABS(v(2)) + cf)
        Max_Lambda(6)=MAX(Max_Lambda(6),ABS(v(3)) + cf)
      END DO
    END DO
  END DO

  v(1) = MAXVAL(Metrics_fTilde(:,:,:,:,l))
  v(2) = MAXVAL(Metrics_gTilde(:,:,:,:,l))
  v(3) = MAXVAL(Metrics_hTilde(:,:,:,:,l))
  eps0 = 1.0/SQRT(MAXVAL(v))

!  h=2.0*eps0/MINVAL(sJ(:,:,:,l))
!  lambda_max=MAXVAL(Max_Lambda(1:3))*h
!  h=0.5/MINVAL(sJ(:,:,:,l))
  h=(8.0/MINVAL(sJ(:,:,:,l)))**(1.0/3.0)
  lambda_max2=MAXVAL(Max_Lambda(4:6))*h

  ! Scaling of artificial viscosity
  nu(l) = nu(l)*lambda_max2/(REAL(PP_N))
!  nu(l) = 0.007

END DO ! l

END SUBROUTINE CalcArtificialViscosity
#endif /*SHOCK_ARTVISC*/
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
  alpha = 1.d0 / (1 + exp(-sharpness * (eta - threshold)/threshold ))
  
  where (alpha < alpha_min)
    alpha = 0.d0
  elsewhere (alpha >= alpha_max)
    alpha = alpha_max
  end where
  
  call ProlongBlendingCoeffToFaces()
  
end subroutine CalcBlendingCoefficient
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
nu_max = 0.
SDEALLOCATE(sVdm_Leg)
SDEALLOCATE(nu)
SDEALLOCATE(nu_Master)
SDEALLOCATE(nu_Slave)
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
