!==================================================================================================================================
! Copyright (c) 2012 - 2019 Florian Hindenlang
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
!> Contains all linear algebra routines that are required for the solution of the sparse linear systems of the type A x = b 
!> arising when using implicit temporal discretisations
!>
!===================================================================================================================================
MODULE MOD_LinAlg
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
!no interfaces, because of different vector dimensions
!INTERFACE MatVec
!  MODULE PROCEDURE MatVec
!END INTERFACE

PUBLIC::MatVec
PUBLIC::Norm_2
PUBLIC::DotProd
!===================================================================================================================================


CONTAINS

!===================================================================================================================================
!> computes the matrix-vector product y = A * x, where A is the system matrix of the spatial operator
!>
!===================================================================================================================================
SUBROUTINE MATVEC(x, y, t)
! MODULES
USE MOD_PreProc   
#if MPI
USE MOD_Globals
#endif /*MPI*/
USE MOD_EigenValues_Vars, ONLY:U0,FD_degree,Nx,Ny,Nz   ! Variables for matrix-vector product
USE MOD_DG,               ONLY:DGTimeDerivative   ! DG Operator
USE MOD_TimeDisc_Vars,    ONLY: TimeStep
USE MOD_DG_Vars,          ONLY:U,Ut               ! DG Operator Variables
USE MOD_Mesh_Vars,        ONLY:nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)  :: t                                               !< current time
  REAL, INTENT(IN)  :: x(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)      !< input vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  REAL, INTENT(OUT) :: y(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)      !< result of matvec
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL              :: eps                                             !< eps for matrix-vector product
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  REAL, DIMENSION(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems):: yp1,yp2,ym1
!===================================================================================================================================
! Compute eps
  eps = SUM(ABS(x*x))
#if MPI
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,eps,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  ELSE
    CALL MPI_REDUCE(eps         ,eps,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  END IF
#endif /*MPI*/
  IF(MPIroot) eps = SQRT(EPSILON(0.)/eps)
#if MPI
  CALL MPI_BCAST(eps,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
#endif /*MPI*/
  SELECT CASE(FD_degree)
  CASE(1) !y~x exact,i forward difference y= d/dx ut ~ (ut(u0+eps*x)-ut(u0))/(eps)
    !U   = U0+eps*x
    CALL UpEpsX(( eps),x)
    CALL DGTimeDerivative(t)
    yp1 = Ut(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)
    !U   = U0
    CALL UpEpsX((0.),x)
    CALL DGTimeDerivative(t)
    ! Compute dU=A*x
    y   = (yp1-Ut(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems))/(eps)
  CASE(2) !y~x^2 exact, central difference y= d/dx ut ~ (ut(u0+eps*x)-ut(u0-eps*x))/(2*eps)
    !U   = U0+eps*x
    CALL UpEpsX(( eps),x)
    CALL DGTimeDerivative(t)
    yp1 = Ut(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)
    !U   = U0-eps*x
    CALL UpEpsX((-eps),x)
    CALL DGTimeDerivative(t)
    ! Compute dU=A*x
    y   = (yp1-Ut(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems))/(2.*eps)
  CASE(20) !fully discrete, one timestep!!
           !y~x^2 exact, central difference y= d/dx (U^(n+1)) ~ (u^(n+1)(u0+eps*x)-u^(n+1)(u0-eps*x))/(2*eps)
    !U   = U0+eps*x
    CALL UpEpsX(( eps),x)
    CALL TimeStep(t)
    yp1 = U(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)
    !U   = U0-eps*x
    CALL UpEpsX((-eps),x)
    CALL TimeStep(t)
    ! Compute dU=A*x
    y   = (yp1-U(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems))/(2.*eps)
  CASE(4) !y~x^4 exact, central difference y= d/dx ut ~ (-ut(+2eps)+8*ut(+eps)-8*ut(-eps)+ut(-2eps))/(12*eps)
    !U   = U0+eps*x
    CALL UpEpsX(( eps),x)
    CALL DGTimeDerivative(t)
    yp1 = Ut(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)
    !U   = U0-eps*x
    CALL UpEpsX((-eps),x)
    CALL DGTimeDerivative(t) 
    ym1 = Ut(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)
    !U   = U0+(2.*eps)*x
    CALL UpEpsX(( 2.*eps),x)
    CALL DGTimeDerivative(t) 
    yp2 = Ut(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)
    !U   = U0-(2.*eps)*x
    CALL UpEpsX((-2.*eps),x)
    CALL DGTimeDerivative(t) 
    ! Compute dU=A*x
    y   = (-yp2+8.*(yp1-ym1) + Ut(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems))/(12.*eps)
  CASE DEFAULT
    STOP 'FD degree not implemented'
  END SELECT !FDorder
END SUBROUTINE MatVec

!===================================================================================================================================
!> 
!===================================================================================================================================
SUBROUTINE UpEpsX(eps_in,x_in)
! MODULES
USE MOD_Preproc
USE MOD_EigenValues_Vars, ONLY:U0,Nx,Ny,Nz   ! Variables for matrix-vector product
USE MOD_DG_Vars,          ONLY:U
USE MOD_Mesh_Vars,        ONLY:nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
  REAL, INTENT(IN)  :: eps_in
  REAL, INTENT(IN)  :: x_in(1:PP_nVar,0:Nx,0:Ny,0:Nz,1:nElems)      !< input vector
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: j,k
!===================================================================================================================================
  IF(Nz.EQ.0)THEN
    IF(Ny.EQ.0)THEN
      DO k=0,PP_N; DO j=0,PP_N
        U(:,:,j,k,:)   = U0(:,:,0,0,:)+eps_in*x_in(:,:,0,0,:)
      END DO;END DO
    ELSE
      DO k=0,PP_N
        U(:,:,:,k,:)   = U0(:,:,:,0,:)+eps_in*x_in(:,:,:,0,:)
      END DO
    END IF
  ELSE
    U   = U0+eps_in*x_in
  END IF
END SUBROUTINE UpEpsX

!===================================================================================================================================
!> computes the euclidean norm of the vector x 
!>
!===================================================================================================================================
FUNCTION Norm_2(x)
! MODULES
#if MPI
USE MOD_Globals
#endif /*MPI*/
USE MOD_EigenValues_Vars, ONLY: localDOF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL    :: x(localDOF)  !< vector to multiply with
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL    :: Norm_2 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  Norm_2 = SUM(x*x)
#if MPI
  IF(MPIroot)THEN
    CALL MPI_REDUCE(MPI_IN_PLACE,norm_2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  ELSE
    CALL MPI_REDUCE(norm_2      ,norm_2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  END IF
#endif /*MPI*/
  IF(MPIroot) Norm_2 = SQRT(Norm_2)
#if MPI
  CALL MPI_BCAST(norm_2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
#endif /*MPI*/
END FUNCTION Norm_2


!===================================================================================================================================
!> computes the dot product of the vectors x and y
!>
!===================================================================================================================================
FUNCTION DotProd(x, y)
! MODULES
#if MPI
USE MOD_Globals
#endif /*MPI*/
USE MOD_EigenValues_Vars, ONLY: localDOF
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL    :: x(localDOF)  !< vector to multiply with
REAL    :: y(localDOF)  !< vector to multiply with
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL    :: DotProd
!===================================================================================================================================
DotProd = SUM(x*y)
#if MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,DotProd,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
#endif /*MPI*/
END FUNCTION DotProd

END MODULE MOD_LinAlg
