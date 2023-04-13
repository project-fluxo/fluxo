!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2022 Andrés Rueda
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
!> This module contains the routines to evaluate the Loehner indicator (under construction)
!==================================================================================================================================
module MOD_Indicator_Loehner
  use MOD_Indicators_vars, only: Indicator_Generic
  implicit none
  
  private
  public :: Indicator_Loehner
  
!============================================================================================================================
!> Class for the modal shock sensor of Persson and Peraire
!> Persson, P. O.; Peraire, J. (2006). "Sub-cell shock capturing for discontinuous Galerkin methods". In 44th AIAA Aerospace Sciences Meeting and Exhibit (p. 112).
!============================================================================================================================
  type, extends(Indicator_Generic) :: Indicator_Loehner
    contains
      procedure :: construct  => Loehner_Construct
      procedure :: compute    => Loehner_Compute
      procedure :: destruct   => Loehner_Destruct
  end type Indicator_Loehner

contains
!===================================================================================================================================
!> Constructor for Indicator_Loehner 
!===================================================================================================================================
  subroutine Loehner_Construct(this,IndicatorQuantity,ReturnLog)
    USE MOD_Globals
    USE MOD_PreProc
    ! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
    !-arguments------------------------------------------------------------------------------------------------------------------
    class(Indicator_Loehner), intent(inout) :: this
    character(len=255)             , intent(in)    :: IndicatorQuantity
    logical, optional              , intent(in)    :: ReturnLog
    !============================================================================================================================
    
    ! Do generic initialization
    call this % Indicator_Generic % construct (IndicatorQuantity,ReturnLog)
    
    IF (PP_N.LT.2) THEN
      CALL abort(__STAMP__,'Polynomial Degree too small for Indicator!',999,999.)
      RETURN
    END IF
    
    ! All is done
    this % InitIsDone = .TRUE.
    
  end subroutine Loehner_Construct
!===================================================================================================================================
!> Destructor for Indicator_Loehner 
!===================================================================================================================================
  subroutine Loehner_Destruct(this)
    implicit none
    class(Indicator_Loehner), intent(inout) :: this
    
    if (.not. this % InitIsDone) return
    
    ! Call generic destructor
    call this % Indicator_Generic % destruct()
    
    this % InitIsDone =.FALSE.
    
  end subroutine Loehner_Destruct
    
!============================================================================================================================
!> Compute the indicator of Löhner
!============================================================================================================================
  function Loehner_Compute(this,U) RESULT (eta)
    use MOD_PreProc
    use MOD_Mesh_Vars,only: nElems !number of local elements
    implicit none
    ! Arguments
    !---------------------------------------------------------------------------------------------------------------------------------
    class(Indicator_Loehner), intent(in) :: this
    real                           , intent(in) :: U(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
    real                                        :: eta(nElems)
    ! Local variables
    !---------------------------------------------------------------------------------------------------------------------------------
    real, dimension(1:1,0:PP_N,0:PP_N,0:PP_N) :: Uind
    integer                                   :: eID,i,j,k
    real                                      :: estimate,u0,up,um,num,den
    !---------------------------------------------------------------------------------------------------------------------------------

    do eID=1, nElems
      
      ! Get the indicator variable (Uind)
      call this % GetIndicator_3D(Uind,U(:,:,:,:,eID))
      
      estimate = 0.0

      do k=0, PP_N ; do j=0, PP_N ; do i=1, PP_N-1
        ! x direction
        u0 = Uind(1,i,j,k)
        up = Uind(1,i+1,j,k)
        um = Uind(1,i-1,j,k)
        num = ABS(up - 2 * u0 + um)
        den = ABS(up - u0) + ABS(u0 - um) + 0.2 * (ABS(up) + 2*ABS(u0) + ABS(um))
        estimate = MAX(estimate, num/den)
      end do       ; end do       ; end do

      do k=0, PP_N ; do j=1, PP_N-1 ; do i=0, PP_N
        ! y direction
        u0 = Uind(1,i,j,k)
        up = Uind(1,i,j+1,k)
        um = Uind(1,i,j-1,k)
        num = ABS(up - 2 * u0 + um)
        den = ABS(up - u0) + ABS(u0 - um) + 0.2 * (ABS(up) + 2*ABS(u0) + ABS(um))
        estimate = MAX(estimate, num/den)
      end do       ; end do    ; end do

      do k=1, PP_N-1 ; do j=0, PP_N ; do i=0, PP_N
        ! z direction
        u0 = Uind(1,i,j,k)
        up = Uind(1,i,j,k+1)
        um = Uind(1,i,j,k-1)
        num = ABS(up - 2 * u0 + um)
        den = ABS(up - u0) + ABS(u0 - um) + 0.2 * (ABS(up) + 2*ABS(u0) + ABS(um))
        estimate = MAX(estimate, num/den)
      end do       ; end do    ; end do
      
      eta(eID) = estimate
      
    end do
    
    where (eta < epsilon(1.0)) eta = epsilon(eta(eID))
    
    if (this % ReturnLog) eta = log10(eta)
    
  end function Loehner_Compute

end module MOD_Indicator_Loehner
