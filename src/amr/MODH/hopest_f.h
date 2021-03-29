!===================================================================================================================================
! Copyright (c) 2018 - 2020 Alexander Astanin
!
! This file is part of FLUXO (github.com/project-fluxo/fluxo). FLUXO is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! FLUXO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLUXO. If not, see <http://www.gnu.org/licenses/>.
!===================================================================================================================================
!===================================================================================================================================
! Here, preprocessor variables for different equation systems and abbreviations for specific expressions are defined
!===================================================================================================================================
#include"defines.h"

! Change if p4est changes its types
! Ideally, this should be in a different file
#define P4EST_F90_TOPIDX INTEGER(KIND=C_INT32_T)
#define P4EST_F90_QCOORD INTEGER(KIND=C_INT32_T)
#define P4EST_F90_LOCIDX INTEGER(KIND=C_INT32_T)
#define P4EST_F90_GLOIDX INTEGER(KIND=C_INT64_T)
#define P4EST_F90_QLEVEL INTEGER(KIND=C_INT8_T)

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#define ADEALLOCATE(A) IF(ASSOCIATED(A)) DEALLOCATE(A)

! Predefined "PARAMETER-like" variables

! Entry position in BC: ATTENTION - TODO: These quantities differ from the ones in defines.h
#define BC_SIZE   4
#define BC_TYPE   1
#define BC_CURVED 2
#define BC_STATE  3
#define BC_ALPHA  4

!entry positions in ElemInfo 
#define ElemInfoSize      6        /*number of entry in each line of ElemInfo*/
#define ELEM_Type         1        /*entry position in ElemInfo */
#define ELEM_Zone         2           
#define ELEM_FirstSideInd 3
#define ELEM_LastSideInd  4
#define ELEM_FirstNodeInd 5
#define ELEM_LastNodeInd  6

!entry positions in SideInfo 
#define SideInfoSize      5        /*number of entry in each line of SideInfo*/
#define SIDE_ID           2
#define SIDE_nbElemID     3
#define SIDE_Flip         4
#define SIDE_BCID         5
