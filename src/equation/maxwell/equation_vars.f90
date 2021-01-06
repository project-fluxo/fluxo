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
!> Contains the parameters for the Maxwell system (including the hyperbolic divergence cleaning variable psi)
!==================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: doCalcSource            !< logical to define if a source term (e.g. exactfunc) is added
REAL                :: c_corr
REAL                :: c_corr2    !< c_corr^2
REAL                :: c_corr_c   !< c_corr*c
REAL                :: c_corr_c2  !< c_corr*c^2
REAL                :: centralflux !=1. purely central flux, =0. purely upwind flux
REAL                :: fDamping
REAL                :: eta_c      !< (c_corr -1 )*c
REAL                :: scr        !< constant for damping in divcorr
INTEGER             :: IniExactFunc
INTEGER             :: BCType(6)=-999
INTEGER             :: BoundaryCondition(6,2)
! Boundary condition arrays
REAL,ALLOCATABLE    :: BCData(:,:,:,:)
INTEGER,ALLOCATABLE :: nBCByType(:)
INTEGER,ALLOCATABLE :: BCSideID(:,:)

CHARACTER(LEN=255),DIMENSION(4),PARAMETER ::  StrVarNames(8)=&
   (/ CHARACTER(LEN=255) :: 'ElectricFieldX', &
                            'ElectricFieldY', &
                            'ElectricFieldZ', &
                            'MagneticFieldX', &
                            'MagneticFieldY', &
                            'MagneticFieldZ', &
                            'Phi           ', &
                            'Psi           ' /)

integer, parameter :: nIndVar = 1
character(len=255) :: IndicatorQuantityNames(nIndVar) = (/character(len=132) :: 'MagneticFieldMag'/)

!REAL,PARAMETER      :: c=299792458.
!REAL,PARAMETER      :: eps0=8.8541878176E-12
REAL                :: c  
REAL                :: c_inv
REAL                :: c2  ! c^2
REAL                :: c2_inv
REAL                :: eps0
REAL                :: mu0 
REAL                :: smu0
INTEGER             :: alpha_shape
REAL                :: shapeFuncPrefix
REAL                :: rCutoff

LOGICAL             :: EquationInitIsDone=.FALSE.
#if (PP_DiscType==2)
INTEGER             :: WhichVolumeFlux
PROCEDURE(i_sub_VolumeFluxAverageVec),POINTER :: VolumeFluxAverageVec
ABSTRACT INTERFACE
  PURE SUBROUTINE i_sub_VolumeFluxAverageVec (UL,UR,metric_L,metric_R,Fstar)
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
    REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
    REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
    REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar          !< transformed central flux
  END SUBROUTINE i_sub_VolumeFluxAverageVec
END INTERFACE
#endif /*PP_DiscType==2*/
ABSTRACT INTERFACE
  pure subroutine i_indicatorFunction(U,ind)
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: ind
  end subroutine i_indicatorFunction
END INTERFACE

contains

!===================================================================================================================================
!> Subroutine to initialize a procedure pointer to a specific indicator function
!===================================================================================================================================
subroutine SetIndicatorFunction(ind_id,IndicatorFunc)
  use MOD_Globals
  implicit none
  !-arguments---------------------------------------
  character(len=255)        , intent(in) :: ind_id
  procedure(i_indicatorFunction),pointer :: IndicatorFunc
  !-------------------------------------------------
  
  select case (trim(ind_id))
    case('MagneticFieldMag') ; IndicatorFunc => Get_MagFieldMag
    case default
      CALL abort(__STAMP__,'Indicator quantity "'//trim(ind_id)//'" is not defined for this equation!',999,999.)
      RETURN
  end select
  
end subroutine SetIndicatorFunction

!===================================================================================================================================
!> Returns the solution
!===================================================================================================================================
pure subroutine Get_MagFieldMag(U,sol)
  implicit none
  !-arguments---------------------------------------
  real,intent(in)  :: U(PP_nVar)
  real,intent(out) :: sol
  !-------------------------------------------------
  
  sol = sqrt(sum(U(4:6)*U(4:6)))
end subroutine Get_MagFieldMag

END MODULE MOD_Equation_Vars
