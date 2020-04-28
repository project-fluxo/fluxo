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

!==================================================================================================================================
!> Contains global variables used by the Analyze modules.
!==================================================================================================================================
MODULE MOD_Analyze_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER              :: NWriteData=1                      !< data output (writing/visualizeing the solution, timeaverages etc.) 
                                                          !< is performed every multiple of Analyze_dt
REAL                 :: Analyze_dt                        !< time intervall at which analysis routines are called
REAL                 :: WriteData_dt                      !< time intervall at which solution data is written
REAL                 :: tWriteData                        !< actual time at which next solution IO will be performed
! precomputed variables
REAL,ALLOCATABLE     :: wGPSurf(:,:)                      !< wGPSurf(i,j)=wGP(i)*wGP(j)
REAL,ALLOCATABLE     :: wGPVol(:,:,:)                     !< wGPVol(i,j,k)=wGP(i)*wGP(j)*wGP(k)
REAL,ALLOCATABLE     :: Surf(:)                           !< surface of each analyze set (e.g. of each boundary condition)
REAL,ALLOCATABLE     :: ElemVol(:)                        !< volume of each element
REAL                 :: Vol                               !< volume of the domain
! Analyze features
LOGICAL              :: doCalcErrorNorms  =.FALSE.        !< input, flag to compute L2 and Linf error norms
LOGICAL              :: doCalcBulk        =.FALSE.        !<handle for analyze: compute Bulk integral of all variables
LOGICAL              :: doCalcMeanFlux    =.FALSE.        !< input, flag to compute the Mean Flux at boundaries
INTEGER              :: AnalyzeExactFunc                  !< input, default = IniExactFunc 

! Analyze to file
LOGICAL              :: doAnalyzeToFile   =.FALSE.        !< marks whether error norms should be written to a file
REAL                 :: iterRestart=0                     !< contains iteration count of previous computation in case a restart is
                                                          !< performed. No restart: 0
REAL                 :: calcTimeRestart=0.                !< contains simulation time at which a restart has been performed
                                                          !< No restart: 0
INTEGER                        :: A2F_ioUnit              !< unit for analyze out file
CHARACTER(LEN=255),ALLOCATABLE :: A2F_VarNames(:)         !< Varnames for analyze to file, setduring InitAnalyze and IniAnalyzeEquation
INTEGER                        :: A2F_nVars               !< Total Number of Variables to be written
INTEGER                        :: A2F_iVar                !< Global Counter, always reset after output 
REAL,ALLOCATABLE               :: A2F_Data(:)             !< one line of data for the Analyze File

! Variables for the specific analyze routines
! ErrorNorms
INTEGER              :: NAnalyze                          !< polynomial degree analysis is performed at (e.g. computation of L2
                                                          !< norms). Number of points: NAnalyze+1
REAL,ALLOCATABLE     :: wGPVolAnalyze(:,:,:)              !< product of GL integration weights used for analyze routines
REAL,ALLOCATABLE     :: Vdm_GaussN_NAnalyze(:,:)          !< Vandermonde for interpolating the solution to analyze points


LOGICAL              :: AnalyzeInitIsDone = .FALSE.       !< marks whether analyze routines have been inittialized
!==================================================================================================================================
END MODULE MOD_Analyze_Vars
