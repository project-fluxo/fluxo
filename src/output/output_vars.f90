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

!==================================================================================================================================
!> Contains global variables provided by the output routines
!==================================================================================================================================
MODULE MOD_Output_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                         :: NVisu                 !< Number of points at which solution is sampled for visualization
REAL,ALLOCATABLE                :: Vdm_GaussN_NVisu(:,:) !< Vandermonde for direct interpolation from computation grid to visu grid
REAL,PARAMETER                  :: FileVersion=0.1       !< version written into output file
CHARACTER(LEN=255),PARAMETER    :: ProgramName='Fluxo'   !< name of program written into output file
CHARACTER(LEN=255)              :: ProjectName           !< Name of the current simulation (mandatory).
CHARACTER(LEN=255)              :: vtuPath               !< Path for vtu files
INTEGER                         :: outputFormat=0        !< File format for visualization. <=0: no visualization, 1: Tecplot binary,
                                                      !< 2: Tecplot ASCII, 3: Paraview binary. Note: Tecplot output is currently
                                                      !< unavailable due to licensing issues.
INTEGER                         :: ASCIIOutputFormat=0   !< File format for ASCII output. 0: CSV, 1: Tecplot
LOGICAL                         :: OutputInitIsDone=.FALSE.  !< marks whether output routines have been initialized
LOGICAL                         :: doPrintStatusLine     !< flag indicating if status line should be printed
LOGICAL                         :: PrimVisuDefault       !< flag indicating if the visualization routines output the primitive variables
INTEGER                         :: nBoundingBoxes        !< number of visualization bounding boxes (default=0)
REAL,ALLOCATABLE                :: VisuBoundingBox(:,:)  !< bounding boxes from input file, size (6,nBoundingBoxes)
                                                      !< with (xmin,xmax,ymin,ymax,zmin,zmax)
integer                         :: nOutVars              !< Default number of output variables for visualize routine
character(LEN=255), allocatable :: strvarnames_tmp(:)   !< Default variable names for visualize routine
!==================================================================================================================================
END MODULE MOD_Output_Vars
