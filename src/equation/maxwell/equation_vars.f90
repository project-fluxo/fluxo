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

!===================================================================================================================================
!> Contains the parameters for the Maxwell system (including the hyperbolic divergence cleaning variable psi)
!===================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: doCalcSource            !< logical to define if a source term (e.g. exactfunc) is added
REAL                :: c_corr
REAL                :: c_corr2    !< c_corr^2
REAL                :: c_corr_c   !< c_corr*c
REAL                :: c_corr_c2  !< c_corr*c^2
REAL                :: fDamping
REAL                :: eta_c      !< (c_corr -1 )*c
REAL                :: scr        !< constant for damping in divcorr
REAL                :: Pi
INTEGER             :: IniExactFunc
INTEGER             :: BCType(6)=-999
INTEGER             :: BoundaryCondition(6,2)
INTEGER,PARAMETER   :: nAuxVar=1
! Boundary condition arrays
REAL,ALLOCATABLE    :: BCData(:,:,:,:)
INTEGER,ALLOCATABLE :: nBCByType(:)
INTEGER,ALLOCATABLE :: BCSideID(:,:)

CHARACTER(LEN=255),DIMENSION(4),PARAMETER :: StrVarNames(8)=(/ CHARACTER(LEN=255) :: 'ElectricFieldX', &
                                                                                     'ElectricFieldY', &
                                                                                     'ElectricFieldZ', &
                                                                                     'MagneticFieldX', &
                                                                                     'MagneticFieldY', &
                                                                                     'MagneticFieldZ', &
                                                                                     'Phi           ', &
                                                                                     'Psi           ' /)

LOGICAL             :: EquationInitIsDone=.FALSE.
!REAL,PARAMETER      :: c=299792458.
!REAL,PARAMETER      :: eps0=8.8541878176E-12
REAL,PARAMETER      :: c=1.  !normalized
REAL,PARAMETER      :: c2=1.  !normalized c^2
REAL,PARAMETER      :: eps0=1.  !normalized
INTEGER             :: alpha_shape
REAL                :: shapeFuncPrefix
REAL                :: rCutoff

INTEGER             :: WhichVolumeFlux
PROCEDURE(),POINTER :: VolumeFluxAverageVec

END MODULE MOD_Equation_Vars
