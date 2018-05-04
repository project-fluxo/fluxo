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
!> Routines for the constraint transport 

!> Contains the routines to 
!> - initialize and finalize the CT global variables
!> - compute the curl of a projected At
!> 
!==================================================================================================================================
MODULE MOD_CT
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitCT
  MODULE PROCEDURE InitCT
END INTERFACE

INTERFACE InitB_CT
  MODULE PROCEDURE InitB_CT
END INTERFACE

INTERFACE CT_TimeDerivative
  MODULE PROCEDURE CT_TimeDerivative
END INTERFACE

!INTERFACE swapBt
!  MODULE PROCEDURE swapBt
!END INTERFACE

INTERFACE FinalizeCT
  MODULE PROCEDURE FinalizeCT
END INTERFACE

PUBLIC:: InitCT
PUBLIC:: InitB_CT
PUBLIC:: CT_TimeDerivative
PUBLIC:: swapBt
PUBLIC:: FinalizeCT
PUBLIC:: DefineParametersCT
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersCT()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("CT")
CALL prms%CreateIntOption('use_CT', &
     "switch on the Constraint transport, -1:off, 0: only init,1: init and Ut", '-1')
CALL prms%CreateIntOption(     "CT_Mover"  , " M for gauss quadrature of projection step, includes +/-1 point (default=(N-1)+2)")
END SUBROUTINE DefineParametersCT

!==================================================================================================================================
!> Allocate all global CT variables
!> 
!==================================================================================================================================
SUBROUTINE InitCT()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_CT_Vars
USE MOD_Interpolation_Vars, ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: MeshInitIsDone,nElems
USE MOD_Mesh_Vars,          ONLY: Metrics_ftilde,Metrics_gtilde,Metrics_htilde,dXGL_N,Elem_xGP
USE MOD_Mesh_vars,          ONLY: nBCSides
USE MOD_Basis,              ONLY: Inv33
USE MOD_ChangeBasis,        ONLY: ChangeBasis3D
USE MOD_ReadInTools,        ONLY: GETINT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: i,j,k,iElem
REAL    :: InvJacMat(3,3)
CHARACTER(LEN=3) :: tmpStr
!===================================================================================================================================
  IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone))THEN
     CALL abort(__STAMP__, &
     'InitCT not ready to be called or already called.',999,999.)
  END IF
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' INIT CT...'
  use_CT=GETINT('use_CT','-1') !switch: 1: use CT also in time update, 0: use CT only for initialization, -1: off
  IF(use_CT.EQ.-1)THEN
    SWRITE(UNIT_stdOut,'(A)')' INIT CT DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')
    RETURN
  END IF
  IF(nBCSides.GT.0)THEN
     CALL abort(__STAMP__, &
      'CT not implemented with BC, only for periodic,cartmesh,1proc')
  END IF
#if MPI
  IF(nProcessors.GT.1)THEN
     CALL abort(__STAMP__, &
         'CT not implemented for MPI, only for periodic,cartmesh,1proc')
  END IF
#endif
  
  ALLOCATE(curlA( 3,0:PP_N,0:PP_N,0:PP_N,nElems))
  ALLOCATE(curlAt(3,0:PP_N,0:PP_N,0:PP_N,nElems))
  curlA=0.
  curlAt=0.
  
  nTotalcurlA=3*(PP_N+1)**3*nElems
  
  ALLOCATE(JacMat(3,3,0:PP_N,0:PP_N,0:PP_N,nElems))
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      InvJacMat(1,:)=Metrics_ftilde(:,i,j,k,iElem)
      InvJacMat(2,:)=Metrics_gtilde(:,i,j,k,iElem)
      InvJacMat(3,:)=Metrics_htilde(:,i,j,k,iElem)
      CALL INV33(InvJacMat,JacMat(:,:,i,j,k,iElem)) !JacMat=(InvJacMat)^-1
    END DO; END DO; END DO !i,j,k=0,PP_N
  END DO !iElem=1,nElems
  
#if PP_NodeType != 2
  STOP 'constraint transport not implemented for Gauss points!'
#endif
  WRITE(tmpStr,'(I3)')PP_N+1
  Mover = GETINT('CT_Mover',tmpStr)
  IF(Mover.LT.PP_N+1) STOP 'Mover must be >=N+1'
  N_ned = PP_N-1
  
  CALL InitProjectionBasis()

  ALLOCATE(CovVec1_int(3,0:Mover,0:Mover,0:Mover,nElems))
  ALLOCATE(CovVec2_int(3,0:Mover,0:Mover,0:Mover,nElems))
  ALLOCATE(CovVec3_int(3,0:Mover,0:Mover,0:Mover,nElems))
  ALLOCATE(Elem_xint(3,0:Mover,0:Mover,0:Mover,nElems))


  DO iElem=1,nElems
    CALL ChangeBasis3D(3,PP_N,Mover,Vdm_GLN_int,dXGL_N(1,:,:,:,:,iElem),CovVec1_int(:,:,:,:,iElem))
    CALL ChangeBasis3D(3,PP_N,Mover,Vdm_GLN_int,dXGL_N(2,:,:,:,:,iElem),CovVec2_int(:,:,:,:,iElem))
    CALL ChangeBasis3D(3,PP_N,Mover,Vdm_GLN_int,dXGL_N(3,:,:,:,:,iElem),CovVec3_int(:,:,:,:,iElem))
    CALL ChangeBasis3D(3,PP_N,Mover,Vdm_GLN_int,Elem_xGP(:,:,:,:,iElem),  Elem_xint(:,:,:,:,iElem))
  END DO !iElem

  
  
  doEvalSource_A=.TRUE. !default
  
  CTInitIsDone=.TRUE.
  SWRITE(UNIT_stdOut,'(A)')' INIT CT DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitCT


SUBROUTINE InitProjectionBasis()
!===================================================================================================================================
! basis functions for projections, are of degree N+1, so they are evaluated first  on (N+2) GL points
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_CT_Vars
USE MOD_Interpolation_Vars,ONLY:NodeTypeG,xGP,wBary
USE MOD_Interpolation     ,ONLY:GetNodesAndWeights
USE MOD_Basis             ,ONLY:InitializeVandermonde
USE MOD_Basis             ,ONLY:JacobiP
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j
REAL    :: x_e(0:N_ned)
REAL    :: w_e(0:N_ned)
REAL    :: wBary_e(0:N_ned)
REAL    :: check
REAL    :: JacPol11_GL(  0:PP_N ,0:N_ned-1)
REAL    :: JacPol01_P_GL(0:PP_N)
REAL    :: JacPol10_P_GL(0:PP_N)
REAL    :: x_GMover(1:Mover-1)
REAL    :: w_GMover(1:Mover-1)
REAL    :: JacPol11_GMover(  1:Mover-1,0:N_ned-1)
REAL    :: JacPol01_P_GMover(1:Mover-1)
REAL    :: JacPol10_P_GMover(1:Mover-1)
REAL    :: e_GMover(1:Mover-1,0:N_ned)
REAL    :: phi_GMover(1:Mover-1,0:N_ned+1)
REAL    :: JP10_pm(2),JP01_pm(2)
!===================================================================================================================================
#if PP_NodeType != 2
  STOP 'InitProjectionbasis only works for Gauss-Lobatto points!'
#endif
  !higher integration rule with gauss Mover-2
  CALL GetNodesAndWeights(Mover-2,NodeTypeG,x_GMover,wIP=w_GMover)
  
  !add boundaries for evaluation
  ALLOCATE(x_int(0:Mover),w_int(0:Mover))
  x_int(1:Mover-1)=x_GMover(:)
  x_int(0)=-1.; x_int(Mover)= 1.

  ALLOCATE(Vdm_GLN_int(0:Mover ,0:PP_N))

  !interpolate GL N to integration points 
  CALL InitializeVandermonde(PP_N  ,Mover,wBary ,xGP ,x_int(:),Vdm_GLN_int(:,:))

  !N_ned=PP_N-1
  IF(N_ned.NE.PP_N-1) STOP 'N_ned must be set to N-1'
  !edge base
  CALL GetNodesAndWeights(N_ned,NodeTypeG,x_e,wIP=w_e,wIPBary=wBary_e) !gauss edge base
 
  !evaluate edge gauss base on GL ponts
  ALLOCATE(Vdm_e_GLN(0:PP_N,0:N_ned))
  CALL InitializeVandermonde(N_ned  ,PP_N,wBary_e    ,x_e  ,xGP,Vdm_e_GLN(:,:))
  DO j=0,N_ned
    Vdm_e_GLN(:,j)=Vdm_e_GLN(:,j)/w_e(j)  !M^-1 of edge gauss base
  END DO

  !integration onto edge: from 1..Mover-1 ONLY!
  ALLOCATE(proj_int_e(0:N_ned,1:Mover-1))
  ! edge base evaluated at Mover 
  CALL InitializeVanderMonde(N_ned,Mover-2,wBary_e,x_e,x_GMover(:),e_GMover(:,:))

  DO j=1,Mover-1
    proj_int_e(:,j)=e_Gmover(j,:)*w_GMover(j)
  END DO

  DO i=0,N_ned
    e_GMover(:,i)=e_GMover(:,i)/w_e(i)
  END DO

  !WRITE(*,*)'check edge'
  !DO i=0,N_ned
  !  WRITE(*,*) MATMUL(proj_int_e(i,:),e_Gmover)
  !END DO
  check=MAXVAL(ABS(MATMUL(proj_int_e(:,:),e_Gmover(:,:))-EYE(N_ned+1)))
  IF(check.GT.1.0e-12)THEN
    CALL abort(__STAMP__, &
           'check edge in initprojectionbasis failed!',-1,check)   
  END IF

  !Jacobi Polynomial (1,1) of degree N_ned-1
  CALL JacobiP(PP_N+1 ,xGP,1,1,0,N_ned-1,JacPol11_GL(:,:))!number of nodes, nodes [-1,1], alpha,beta
                                                          !min./max.degree, JP(i,j) at node i and degree j
  CALL JacobiP(Mover-1,x_GMover,1,1, 0,N_ned-1,JacPol11_GMover(:,:))
  
  !Jacobi Polynomial (1,0)/(0,1) of degree N_ned, evaluated at GL_N
  CALL JacobiP(PP_N+1,xGP,   1, 0, N_ned,N_ned,JacPol10_P_GL(:)) 
  CALL JacobiP(PP_N+1,xGP,   0, 1, N_ned,N_ned,JacPol01_P_GL(:)) 
  !Jacobi Polynomial (1,0)/(0,1) of degree N_ned,evaluated at +/- 1
  CALL JacobiP(  2,(/-1.,1./), 1, 0, N_ned,N_ned,JP10_pm         ) 
  CALL JacobiP(  2,(/-1.,1./), 0, 1, N_ned,N_ned,JP01_pm         )
  
  !Jacobi Polynomial (1,0)/(0,1) of degree N_ned, evaluated at GMover
  CALL JacobiP(Mover-1,x_GMover,   1, 0, N_ned,N_ned,JacPol10_P_GMover(:)) 
  CALL JacobiP(Mover-1,x_GMover,   0, 1, N_ned,N_ned,JacPol01_P_GMover(:)) 
  
  !Jacobi Polynomial (0,1) of degree N_ned
  !phiminus = (1-x)*P^(1,0)_p(x) /(2*P^(1,0)_p(-1) )
  !phiplus  = (1+x)*P^(0,1)_p(x) /(2*P^(0,1)_p(+1) )
  
  
  ALLOCATE(Vdm_phi_GLN(0:PP_N,0:N_ned+1))
  Vdm_phi_GLN(:,0)=(1.-xGP(:))*JacPol10_P_GL(:)/(2.*JP10_pm(1))     ! -1 edge basis in volume
  DO i=1,N_ned
    Vdm_phi_GLN(:,i)=(1.-xGP(:))*(1.+xGP(:))*JacPol11_GL(:,i-1)
  END DO
  Vdm_phi_GLN(:,N_ned+1)=(1.+xGP(:))*JacPol01_P_GL(:)/(2.*JP01_pm(2)) ! +1 edge basis in volume


  phi_GMover( :,0)=(1.-x_GMover(:))*JacPol10_P_GMover(:)/(2.*JP10_pm(1)) 
  DO i=1,N_ned
    phi_Gmover( :,i)=(1.-x_Gmover(:))*(1.+x_Gmover(:))*JacPol11_Gmover(:,i-1) !face / volume base
  END DO
  phi_GMover( :,N_ned+1)=(1.+x_GMover(:))*JacPol01_P_GMover(:)/(2.*JP01_pm(2)) 


  !include integration weight here, and a copy of the boudnary values at (-1,1)
  ALLOCATE(proj_int_phi(0:N_ned+1,0:Mover))
  proj_int_phi=0.
  proj_int_phi(0,0)=1.  !boundary copy
  proj_int_phi(N_ned+1,Mover)=1.
  DO i=1,N_ned
    DO j=1,Mover-1
      proj_int_phi(i,j)=JacPol11_GMover(j,i-1)*w_GMover(j)
    END DO
  END DO

  !WRITE(*,*)'check face'
  !DO i=1,N_ned
  !  WRITE(*,*) MATMUL(proj_int_phi(i,1:Mover-1),phi_Gmover(1:Mover-1,:))
  !END DO
  check=MAXVAL(ABS(MATMUL(proj_int_phi(1:N_ned,1:Mover-1),phi_Gmover(1:Mover-1,1:N_ned))-EYE(N_ned)))
  IF(check.GT.1.0e-12)THEN
    CALL abort(__STAMP__, &
           'check face in initprojectionbasis failed!',-1,check)   
  END IF
  check=MAXVAL(ABS(MATMUL(proj_int_phi(1:N_ned,1:Mover-1),phi_Gmover(1:Mover-1,0))))
  IF(check.GT.1.0e-12)THEN
    CALL abort(__STAMP__, &
           'check face Bc -1 in initprojectionbasis failed!',-1,check)   
  END IF
  check=MAXVAL(ABS(MATMUL(proj_int_phi(1:N_ned,1:Mover-1),phi_Gmover(1:Mover-1,N_ned+1))))
  IF(check.GT.1.0e-12)THEN
    CALL abort(__STAMP__, &
           'check face Bc +1 in initprojectionbasis failed!',-1,check)   
  END IF

  ! GERRITSMA edge basis function is defined as 
  ! e_j(xi) = - sum_k=0^(j) d/dxi l_k(xi)  , j=0...N-1 
  !
  ! evalute at GL N points xi_i:
  ! e_j(xi_i) = - sum_k=0^(j)  [d/dxi l_k(xi) ]_(xi_i) 
  !
  ! where we can use the definition of the D matrix: D_ik = [d/dxi l_k(xi) ]_(xi_i) 

  !CALL PolynomialDerivativeMatrix(PP_N,xGP,D_GL)
  !DO i=0,PP_N; DO j=0,PP_N-1
  !  Vdm_e_GLN(i,j)=-SUM(D_GL(i,0:j))
  !END DO; END DO ! i,j=0,PP_N

  CONTAINS
  
  FUNCTION EYE(N)
  INTEGER,INTENT(IN) :: N
  REAL               :: EYE(N,N)
  EYE(:,:)=0.
  DO i=1,N
    EYE(i,i)=1.
  END DO
  END  FUNCTION EYE

 
END SUBROUTINE InitProjectionBasis

!==================================================================================================================================
!> Initialize divergence free magnetic field
!> 
!==================================================================================================================================
SUBROUTINE InitB_CT(ExactFunc,U_inout)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_CT_Vars
USE MOD_Mesh_Vars,          ONLY: nElems,Elem_xGP
USE MOD_ChangeBasis, ONLY: ChangeBasis3D_XYZ
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: exactFunc
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)    :: U_inout(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER :: i,j,k,iElem
REAL    :: Acov(3,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL    :: Acov_e(3,0:PP_N-1,0:PP_N  ,0:PP_N  ,nElems)
REAL    :: Acart_loc(3),curlA_loc(3),constB(3)
REAL    :: Acov_int(3,0:Mover,0:Mover,0:Mover)
LOGICAL :: discreteDivB
!===================================================================================================================================
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      CALL EvalExactFunc_A(ExactFunc,Elem_xGP(:,i,j,k,iElem),Acart_loc(:), &
                                                  curlA(:,i,j,k,iElem),constB,discreteDivB)
    END DO; END DO; END DO !i,j,k=0,PP_N
  END DO !iElem=1,nElems
  IF(.NOT.discreteDivB)THEN
    DO iElem=1,nElems
      DO k=0,Mover; DO j=0,Mover; DO i=0,Mover
        CALL EvalExactFunc_A(ExactFunc,Elem_xint(:,i,j,k,iElem),Acart_loc(:), &
                                                    curlA_loc(:),constB,discreteDivB)
        Acov_int(1,i, j,k)=SUM(CovVec1_int(:,i,j,k,iElem)*Acart_loc(:))
        Acov_int(2,j, i,k)=SUM(CovVec2_int(:,i,j,k,iElem)*Acart_loc(:))
        Acov_int(3,k, i,j)=SUM(CovVec3_int(:,i,j,k,iElem)*Acart_loc(:))
      END DO; END DO; END DO !i,j,k=0,PP_N
      CALL ChangeBasis3D_XYZ(3,Mover-2,PP_N-1,proj_int_e, &
                               Mover,PP_N  ,proj_int_phi, &
                               Mover,PP_N  ,proj_int_phi, &
                               Acov_int(:,1:Mover-1, :,:),Acov_e(:,:,:,:,iElem))

    END DO !iElem=1,nElems
    ! stitch edges/faces together 
    CALL ProjectToNedelec(ACov_e,ACov)
     
    CALL ComputeCartCurl(ACov,curlA)
    curlA(1,:,:,:,:)=curlA(1,:,:,:,:)+constB(1)
    curlA(2,:,:,:,:)=curlA(2,:,:,:,:)+constB(2)
    curlA(3,:,:,:,:)=curlA(3,:,:,:,:)+constB(3)
  END IF
  WRITE(*,*)'CHECK B=curlA_init...'
  CALL CheckB(curlA)
  !replace B in U_in
  CALL swapB(U_inout,curlA)
   
END SUBROUTINE InitB_CT

!==================================================================================================================================
!> compute ODE time derivative of A, A_t=-(B x v [+eta J] )
!> project to Necelec base (in covariant components!), compute curl back to cartesian components. 
!==================================================================================================================================
SUBROUTINE CT_TimeDerivative()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars    ,ONLY: U       !! input  !!
USE MOD_CT_Vars 
USE MOD_Equation_Vars    ,ONLY: IniExactFunc
USE MOD_Mesh_Vars  ,ONLY: nElems
USE MOD_ChangeBasis, ONLY: ChangeBasis3D,ChangeBasis3D_XYZ
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: AtCov(3,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL    :: AtCov_e(3,0:PP_N-1,0:PP_N,0:PP_N,nElems)
REAL    :: At(3)
REAL    :: AtCov_int(3,0:Mover,0:Mover,0:Mover)
REAL    :: At_source(3,0:Mover,0:Mover,0:Mover)
REAL    :: U_int(PP_nVar,0:Mover,0:Mover,0:Mover)
!==================================================================================================================================
  At_source=0.
  DO iElem=1,nElems
    IF(doEvalSource_A)THEN
      DO k=0,Mover;DO j=0,Mover; DO i=0,Mover
        CALL EvalSource_A(IniExactFunc,Elem_xint(:,i,j,k,iElem),At_source(:,i,j,k))
      END DO; END DO; END DO 
    END IF !doSource
    CALL ChangeBasis3D(PP_nVar,PP_N,Mover,Vdm_GLN_int,U(:,:,:,:,iElem),U_int(:,:,:,:))
    DO k=0,Mover;DO j=0,Mover; DO i=0,Mover
      At = CROSS(U_int(2:4,i,j,k)/U_int(1,i,j,k),U_int(6:8,i,j,k))+At_source(:,i,j,k)  !=v x B = - B x v
      AtCov_int(1,i, j,k)= SUM(covVec1_int(:,i,j,k,iElem)*At(:))
      AtCov_int(2,j, i,k)= SUM(covVec2_int(:,i,j,k,iElem)*At(:))
      AtCov_int(3,k, i,j)= SUM(covVec3_int(:,i,j,k,iElem)*At(:))
    END DO; END DO; END DO 
    CALL ChangeBasis3D_XYZ(3,Mover-2,PP_N-1,proj_int_e, &
                             Mover,PP_N  ,proj_int_phi, &
                             Mover,PP_N  ,proj_int_phi, &
                             Atcov_int(:,1:Mover-1, :,:),Atcov_e(:,:,:,:,iElem))
  END DO !iElem=1,nElems
  
  CALL ProjectToNedelec(AtCov_e,AtCov)
   
  CALL ComputeCartCurl(AtCov,curlAt)
  
  !WRITE(*,*)'check B=curlA_t  !!!!!!!!!'
  !CALL checkB(curlAt)
END SUBROUTINE CT_TimeDerivative

!==================================================================================================================================
!> compute the cartesian components of the curl, only locally per element,
!> such that the transformed divergence in the volume is discretely zero (1/J sum_i d/dxi^i Ja^i . v), using the inverse of Ja^i .
!> input can be cartesian or covariant 
!==================================================================================================================================
SUBROUTINE ComputeCartCurl(vec_cov,curl_cart)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_DG_Vars    ,ONLY: D
USE MOD_CT_Vars    ,ONLY: JacMat
USE MOD_Mesh_Vars  ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)      :: vec_cov(3,0:PP_N,0:PP_N,0:PP_N,nElems) !< vector in covariant components 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)     :: curl_cart(3,0:PP_N,0:PP_N,0:PP_N,nElems) !< curl of vector in cartesian components 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,l,iElem
REAL    :: dvdxi(3),dvdeta(3),dvdzeta(3)
!==================================================================================================================================
  DO iElem=1,nElems
    
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      dvdxi(:)   = D(i,0)*vec_cov(:,0,j,k,iElem)
      dvdeta(:)  = D(j,0)*vec_cov(:,i,0,k,iElem)
      dvdzeta(:) = D(k,0)*vec_cov(:,i,j,0,iElem)
      DO l=1,PP_N
        dvdxi(:)   = dvdxi(:)   +D(i,l)*vec_cov(:,l,j,k,iElem)
        dvdeta(:)  = dvdeta(:)  +D(j,l)*vec_cov(:,i,l,k,iElem)
        dvdzeta(:) = dvdzeta(:) +D(k,l)*vec_cov(:,i,j,l,iElem)
      END DO !l
      !cartesian components from (J B)^i = (Ja^i.Bcart) scaled contravariant components
      curl_cart(:,i,j,k,iElem)= MATMUL(JacMat(:,:,i,j,k,iElem),(/(dvdeta( 3)-dvdzeta(2)), &
                                                                 (dvdzeta(1)-dvdxi(  3)), &
                                                                 (dvdxi(  2)-dvdeta( 1)) /))
!           D(l,j)*vec_cov(3,i,l,k)-D(l,k)*vec_cov(2,i,j,l)
!           D(l,k)*vec_cov(1,i,j,l)-D(l,i)*vec_cov(3,l,j,k)
!           D(l,i)*vec_cov(2,l,j,k)-D(l,j)*vec_cov(1,i,l,k)
    END DO; END DO; END DO !i,j,k=0,PP_N                      
  END DO !iElem=1,nElems

END SUBROUTINE ComputeCartCurl


!==================================================================================================================================
!> project discontinuous covariant vector represented by LGL polynomials to Nedelec Finite Element LGL polynomial basis
!> !!! WORKS ONLY FOR 1 PROC &  PERIODIC STRUCTURED BLOCK MESH !!!
!==================================================================================================================================
SUBROUTINE projectToNedelec(vec_in,vec)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars,ONLY: nElems,SideToElem,ElemToSide,Elem_xGP
USE MOD_CT_Vars,ONLY: Vdm_e_GLN,Vdm_phi_GLN
USE MOD_Changebasis,ONLY: ChangeBasis2D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: vec_in(3,0:PP_N-1,0:PP_N,0:PP_N,nElems) !<  covariant components as edge dof, (1,m,j,k), (2,m,i,k), (3,m,i,j)
REAL,INTENT(OUT):: vec(3,0:PP_N,0:PP_N,0:PP_N,nElems)  !< covariant components stiched together and interpolated to GL N 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL       :: vec1(0:PP_N-1,0:PP_N,0:PP_N)
REAL       :: vec2(0:PP_N-1,0:PP_N,0:PP_N)
REAL       :: vec3(0:PP_N-1,0:PP_N,0:PP_N)
INTEGER    :: i,j,k,iElem,nb(-1:1,-1:1),xx,yy,zz
!==================================================================================================================================
  DO iElem=1,nElems
    !A1: C0 in y,z
    nb( 0, 0) =iElem  ! 
    nb( 0,-1) =GetnbElemID(ZETA_MINUS,iElem)      !    z- 
    nb( 0, 1) =GetnbElemID(ZETA_PLUS ,iElem)      !    z+ 
    nb(-1, 0) =GetnbElemID( ETA_MINUS,iElem)      ! y-
    nb( 1, 0) =GetnbElemID( ETA_PLUS ,iElem)      ! y+
    nb(-1,-1) =GetnbElemID( ETA_MINUS,nb(0,-1))   ! y-,z-
    nb( 1,-1) =GetnbElemID( ETA_PLUS ,nb(0,-1))   ! y+,z-
    nb(-1, 1) =GetnbElemID( ETA_MINUS,nb(0, 1))   ! y-,z+
    nb( 1, 1) =GetnbElemID( ETA_PLUS ,nb(0, 1))   ! y+,z+
    !corners
    DO yy=0,1; DO zz=0,1
      vec1(:,yy*PP_N,zz*PP_N) = 0.25*( vec_in(1,:,   0,   0,nb(yy  ,zz  )) &
                                      +vec_in(1,:,   0,PP_N,nb(yy  ,zz-1)) &
                                      +vec_in(1,:,PP_N,   0,nb(yy-1,zz  )) &
                                      +vec_in(1,:,PP_N,PP_N,nb(yy-1,zz-1)) )
    END DO; END DO !yy,zz
    !inner surface points
    DO k=1,PP_N-1 
      vec1( :,   0,k)=0.5*( vec_in(1,:,   0,k,nb( 0, 0)) &
                           +vec_in(1,:,PP_N,k,nb(-1, 0)) ) !y-
      vec1( :,PP_N,k)=0.5*( vec_in(1,:,PP_N,k,nb( 0, 0)) &
                           +vec_in(1,:,   0,k,nb( 1, 0)) ) !y+
    END DO
    DO j=1,PP_N-1 
      vec1( :,j,   0)=0.5*( vec_in(1,:,j,   0,nb( 0, 0)) &
                           +vec_in(1,:,j,PP_N,nb( 0,-1)) ) !z-
      vec1( :,j,PP_N)=0.5*( vec_in(1,:,j,PP_N,nb( 0, 0)) &
                           +vec_in(1,:,j,   0,nb( 0, 1)) ) !z+
    END DO
    !inner volume points
    DO k=1,PP_N-1 
      DO j=1,PP_N-1 
        vec1(:,j,k)=vec_in(1,:,j,k,iElem)
      END DO
    END DO
    !A2: C0 in x,z
    nb( 0, 0) =iElem  ! 
    nb( 0,-1) =GetnbElemID(ZETA_MINUS,iElem)    !    z- 
    nb( 0, 1) =GetnbElemID(ZETA_PLUS ,iElem)    !    z+ 
    nb(-1, 0) =GetnbElemID(  XI_MINUS,iElem)    ! x-
    nb( 1, 0) =GetnbElemID(  XI_PLUS ,iElem)    ! x+
    nb(-1,-1) =GetnbElemID(  XI_MINUS,nb(0,-1)) ! x-,z-
    nb( 1,-1) =GetnbElemID(  XI_PLUS ,nb(0,-1)) ! x+,z-
    nb(-1, 1) =GetnbElemID(  XI_MINUS,nb(0, 1)) ! x-,z+
    nb( 1, 1) =GetnbElemID(  XI_PLUS ,nb(0, 1)) ! x+,z+
    !corners
    DO xx=0,1; DO zz=0,1
      vec2(:,xx*PP_N,zz*PP_N)= 0.25*( vec_in(2,:,   0,   0,nb(xx  ,zz  )) &
                                     +vec_in(2,:,   0,PP_N,nb(xx  ,zz-1)) &
                                     +vec_in(2,:,PP_N,   0,nb(xx-1,zz  )) &
                                     +vec_in(2,:,PP_N,PP_N,nb(xx-1,zz-1)) )
    END DO; END DO !xx,zz
    !inner surface points
    DO k=1,PP_N-1 
      vec2(:,   0, k)=0.5*( vec_in(2,:,   0,k,nb( 0, 0)) &
                           +vec_in(2,:,PP_N,k,nb(-1, 0)) ) !x-
      vec2(:,PP_N, k)=0.5*( vec_in(2,:,PP_N,k,nb( 0, 0)) &
                           +vec_in(2,:,   0,k,nb( 1, 0)) ) !x+
    END DO
    DO i=1,PP_N-1 
      vec2(:, i,   0)=0.5*( vec_in(2,:,   i,   0,nb( 0, 0)) &
                           +vec_in(2,:,   i,PP_N,nb( 0,-1)) ) !z-
      vec2(:, i,PP_N)=0.5*( vec_in(2,:,   i,PP_N,nb( 0, 0)) &
                           +vec_in(2,:,   i,   0,nb( 0, 1)) ) !z+
    END DO
    !inner volume points
    DO k=1,PP_N-1 
      DO i=1,PP_N-1 
        vec2(:,i,k)=vec_in(2,:,i,k,iElem)
      END DO
    END DO
    !A3: C0 in x,y
    nb( 0, 0) =iElem  ! 
    nb( 0,-1) =GetnbElemID( ETA_MINUS,iElem)       !    y- 
    nb( 0, 1) =GetnbElemID( ETA_PLUS ,iElem)       !    y+ 
    nb(-1, 0) =GetnbElemID(  XI_MINUS,iElem)       ! x-
    nb( 1, 0) =GetnbElemID(  XI_PLUS ,iElem)       ! x+
    nb(-1,-1) =GetnbElemID(  XI_MINUS,nb(0,-1))    ! x-,y-
    nb( 1,-1) =GetnbElemID(  XI_PLUS ,nb(0,-1))    ! x+,y-
    nb(-1, 1) =GetnbElemID(  XI_MINUS,nb(0, 1))    ! x-,y+
    nb( 1, 1) =GetnbElemID(  XI_PLUS ,nb(0, 1))    ! x+,y+
    !corners
    DO xx=0,1; DO yy=0,1
      vec3(:,xx*PP_N,yy*PP_N)= 0.25*( vec_in(3,:,   0,   0,nb(xx  ,yy  )) &
                                     +vec_in(3,:,   0,PP_N,nb(xx  ,yy-1)) &
                                     +vec_in(3,:,PP_N,   0,nb(xx-1,yy  )) &
                                     +vec_in(3,:,PP_N,PP_N,nb(xx-1,yy-1)) )
    END DO; END DO !xx,yy
    !inner surface points
    DO j=1,PP_N-1 
      vec3(:,   0, j)=0.5*( vec_in(3,:,   0,   j,nb( 0, 0)) &
                           +vec_in(3,:,PP_N,   j,nb(-1, 0)) ) !x-
      vec3(:,PP_N, j)=0.5*( vec_in(3,:,PP_N,   j,nb( 0, 0)) &
                           +vec_in(3,:,   0,   j,nb( 1, 0)) ) !x+
    END DO
    DO i=1,PP_N-1 
      vec3(:, i,   0)=0.5*( vec_in(3,:,   i,   0,nb( 0, 0)) &
                           +vec_in(3,:,   i,PP_N,nb( 0,-1)) ) !y-
      vec3(:, i,PP_N)=0.5*( vec_in(3,:,   i,PP_N,nb( 0, 0)) &
                           +vec_in(3,:,   i,   0,nb( 0, 1)) ) !y+
    END DO
    !inner volume points
    DO j=1,PP_N-1 
      DO i=1,PP_N-1 
        vec3(:,i,j)=vec_in(3,:,i,j,iElem)
      END DO
    END DO
    
    !change basis from Nedelec to GLN vec1,vec2,vec3 --> vec(3,...,iElem)
    CALL ChangeBasis2D(PP_N,PP_N,PP_N,Vdm_phi_GLN,vec1,vec1)
    DO k=0,PP_N; DO j=0,PP_N
        vec(1,:,j,k,iElem)=MATMUL(Vdm_e_GLN,vec1(:,j,k))
    END DO; END DO
    CALL ChangeBasis2D(PP_N,PP_N,PP_N,Vdm_phi_GLN,vec2,vec2)
    DO k=0,PP_N; DO i=0,PP_N
        vec(2,i,:,k,iElem)=MATMUL(Vdm_e_GLN,vec2(:,i,k))
    END DO; END DO
    CALL ChangeBasis2D(PP_N,PP_N,PP_N,Vdm_phi_GLN,vec3,vec3)
    DO j=0,PP_N; DO i=0,PP_N
        vec(3,i,j,:,iElem)=MATMUL(Vdm_e_GLN,vec3(:,i,j))
    END DO; END DO

  END DO !iElem=1,nElems


CONTAINS

  FUNCTION GETnbElemID(locSide,elemID)
    !----------------------------------
    INTEGER,INTENT(IN) :: locSide
    INTEGER,INTENT(IN) :: elemID
    INTEGER            :: GETnbElemID
    !----------------------------------
    INTEGER            :: SideID,Flip
    !----------------------------------
    SideID=ElemToSide(E2S_SIDE_ID,locSide,ElemID)
    Flip  =ElemToSide(E2S_FLIP   ,locSide,ElemID)
    IF(flip.EQ.0)THEN
      GETnbElemID=SideToElem(S2E_NB_ELEM_ID,SideID)
    ELSE
      GETnbElemID=SideToElem(S2E_ELEM_ID,SideID)
    END IF

  END FUNCTION GETnbElemID

END SUBROUTINE projectToNedelec


!==================================================================================================================================
!> evaluate the magnetic vector potential A, belonging to the magnetic field
!> of a given IniExactfunction 
!> B_1= d(A_3)/dy - d(A_2)/dz
!> B_2= d(A_1)/dz - d(A_3)/dx
!> B_3= d(A_2)/dx - d(A_1)/dy
!==================================================================================================================================
SUBROUTINE EvalExactFunc_A(ExactFunc,x,Acart,Bcart,constB,discreteDivB)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars  ,ONLY: s2mu_0,IniHalfwidth,IniFrequency,IniAmplitude,IniWaveNumber
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: ExactFunc
REAL   ,INTENT(IN)  :: x(3) 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Acart(3)
REAL,INTENT(OUT)    :: Bcart(3)
REAL,INTENT(OUT)    :: constB(3)
LOGICAL,INTENT(OUT) :: discreteDivB   !< set to true if divergence of Bcart is already DISCRETELY zero!
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: Omega,e,nx,ny,phi_alv
REAL  :: q0,a,r0,va
!==================================================================================================================================
  constB=0.
  SELECT CASE(ExactFunc)
  CASE(3) !alfven wave domain [-1,1]^3, A is periodic
    Omega=2.*PP_Pi*IniFrequency
    !r=2,sqr=1,teval=0.
    e   = 0.2
    nx  = 1./SQRT(2**2+1.)
    ny  = 2./SQRT(2**2+1.)
    phi_alv = omega*(nx/ny*(x(1)-1.) + (x(2)-1.))
  
  
    Acart(1) = 0.
    Acart(2) = e*ny/(nx*Omega)*COS(phi_alv)
    Acart(3) = e*ny/Omega*SIN(phi_alv) 
  
    constB(1)=nx
    constB(2)=ny
    constB(3)=0.
  
    discreteDivB=.FALSE.
  
    Bcart(1) =nx+ e*ny*COS(phi_alv)
    Bcart(2) =ny -e*nx*COS(phi_alv)
    Bcart(3) =   -e*SIN(phi_alv) 

  CASE(31,32,33) ! linear shear alfven wave , linearized MHD,|B|>=1 , p,rho from inirefstate 
           !IniWavenumber=(k_x,k_yk_z): k_parallel=k_x*e_x+k_y*e_y, k_perp=k_z*e_z
           !IniAmplitude should be small compare to density and pressure (1e-08)
           !31: backward moving, 32: forward moving, 33: standing (31+32!)
    constB(3)=0.
    IF(IniWavenumber(2).LT.0.01) THEN
      constB(1:2)=(/1.,0./)
    ELSEIF(IniWavenumber(1).LT.0.01) THEN
      constB(1:2)=(/0.,1./)
    ELSE
      constB(1:2)=(/MIN(1.,IniWavenumber(1)/IniWavenumber(2)),MIN(1.,IniWavenumber(2)/IniWavenumber(1)) /)
    END IF
    q0=SQRT(SUM(constB(:)**2))        !=|B_0|
    a=1.  !SQRT(mu_0*rho_0)  !=|B_0|/va = sqrt(mu_0*rho_0) ASSUMING rho_0=11,,mu_0=1
    IF(Exactfunc.EQ.32) a=-a ! case(32) -va!
    va=q0/a      !(+va)=|B_0|/sqrt(mu_0*rho_0)
    IF(Exactfunc.EQ.33)THEN
      r0=0. !switch for standing wave
    ELSE
      r0=1.
    END IF
    !xc(1:3)=x(1:3)-r0*constB(1:3)/a*tEval ! b_0/a = B_0/|B_0|*va
    IF(IniWaveNumber(1).LT.IniWaveNumber(2))THEN 
      Acart(3)=IniAmplitude/(2.0*PP_Pi*IniWaveNumber(2))*COS(2.0*PP_Pi*SUM(x(:)*IniWavenumber(:)))
    ELSE
      Acart(3)=IniAmplitude/(2.0*PP_Pi*IniWaveNumber(1))*COS(2.0*PP_Pi*SUM(x(:)*IniWavenumber(:)))
    END IF
    Acart(3)=Acart(3)*r0*a/q0
    Acart(1) =0. 
    Acart(2) =0.
  
    discreteDivB=.FALSE.
    e=IniAmplitude*COS(2.0*PP_Pi*SUM(x(:)*IniWavenumber(:)))*(-va*2.0*PP_Pi*SUM(constB(:)*IniWavenumber(:)))
    Bcart(1) =constB(1) + e*r0*a/q0*constB(2)
    Bcart(2) =constB(2) - e*r0*a/q0*constB(1)
    Bcart(3) = 0.

  CASE(3001) !alfven wave domain [-1,1]^3, A is periodic, rotated
    Omega=2.*PP_Pi*IniFrequency
    !r=2,sqr=1,teval=0.
    e   = 0.2
    nx  = 1./SQRT(2**2+1.)
    ny  = 2./SQRT(2**2+1.)
    phi_alv = omega*(nx/ny*(x(1)-1.) + (x(3)-1.))
  
  
    Acart(1) = 0.
    Acart(2) =-e*ny/Omega*SIN(phi_alv) 
    Acart(3) =-e*ny/(nx*Omega)*COS(phi_alv)
  
    constB(1)=nx
    constB(2)=0
    constB(3)=ny
  
    discreteDivB=.FALSE.
  
    Bcart(1) =nx+ e*ny*COS(phi_alv)
    Bcart(2) =   -e*SIN(phi_alv) 
    Bcart(3) =ny -e*nx*COS(phi_alv)
  
  CASE(7)
    Omega=PP_Pi*IniFrequency
    Acart(1)= 2.*IniAmplitude*SIN(Omega*SUM(x(:)))
    Acart(2)=    IniAmplitude*SIN(Omega*SUM(x(:)))
    Acart(3)=  - IniAmplitude*SIN(Omega*SUM(x(:)))
    discreteDivB=.FALSE.
  
    !B1=dA3/dy-dA2/dz = -2*IniAmplitude*Omega*COS(Omega*SUM(x))
    !B2=dA1/dz-dA3/dx =  3*IniAmplitude*Omega*COS(Omega*SUM(x))
    !B3=dA2/dx-dA1/dy = -1*IniAmplitude*Omega*COS(Omega*SUM(x))
  
    Bcart(1) = -2*IniAmplitude*Omega*COS(Omega*SUM(x))
    Bcart(2) =  3*IniAmplitude*Omega*COS(Omega*SUM(x))
    Bcart(3) = -  IniAmplitude*Omega*COS(Omega*SUM(x))
  CASE(73) !IS NOT PERIODIC!!!
    Bcart(1)=TANH((ABS(x(2))-0.5)/IniHalfwidth)
    Bcart(2)=0.
    Bcart(3)=1.
    discreteDivB=.TRUE. !A,constB are not used...
    Acart(1) = 0.
    Acart(2) = - x(3)*Bcart(1) !NOT PERIODIC!!!!
    Acart(3) = 0. 
    !constant Bfield,to be added
    constB(3)=1.0
  CASE(74) !IS NOT PERIODIC !!!
    Bcart(1)=TANH((ABS(x(2))-0.5)/IniHalfwidth)
    Bcart(2)=0.
    Bcart(3)=SQRT(1.0-Bcart(1)**2)
    discreteDivB=.TRUE.  ! A,constB are not used...
    Acart(1) = 0.
    Acart(2) = x(1)*Bcart(3)-x(3)*Bcart(1) !NOT PERIODIC
    Acart(3) = 0. 
  
  CASE DEFAULT
    STOP 'this exactfunc is not implemented in EvalExactFunc_A'
  END SELECT !ExactFunc

END SUBROUTINE EvalExactFunc_A

!==================================================================================================================================
!> evaluate the magnetic vector potential A, belonging to the magnetic field
!> of a given IniExactfunction 
!> B_1= d(A_3)/dy - d(A_2)/dz
!> B_2= d(A_1)/dz - d(A_3)/dx
!> B_3= d(A_2)/dx - d(A_1)/dy
!==================================================================================================================================
SUBROUTINE EvalSource_A(ExactFunc,x,Atcart)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_CT_Vars,        ONLY: doEvalSource_A
!USE MOD_Equation_Vars  ,ONLY: s2mu_0,IniHalfwidth,IniFrequency,IniAmplitude
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: ExactFunc
REAL   ,INTENT(IN)  :: x(3) 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Atcart(3)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL  :: Omega
!==================================================================================================================================
  Atcart=0.
  SELECT CASE (ExactFunc)
  CASE(4,5,6) 
    STOP 'you need to provide a source for the magnetic vector potential At here!!'
  CASE DEFAULT
    ! No source -> do nothing
    doEvalSource_A=.FALSE. 
  END SELECT ! ExactFunction

END SUBROUTINE EvalSource_A

!==================================================================================================================================
!> Replaces B of current U by curlA, and changing total energy and keeping the same pressure
!==================================================================================================================================
SUBROUTINE swapB(U_inout,curlA_in)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars  ,ONLY: s2mu_0
USE MOD_Mesh_Vars      ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: curlA_in(3,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: U_inout(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: EmagU,EmagCurlA
!REAL    :: diffB(3),mindiffB(3),maxdiffB(3) 
!==================================================================================================================================
  !maxdiffB=-1.0e20
  !mindiffB=1.0e20
  !DO iElem=1,nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
  !  diffB=U_inout(6:8,i,j,k,iElem)-curlA_in(:,i,j,k,iElem)
  !  mindiffB=MIN(mindiffB,diffB)
  !  maxdiffB=MAX(maxdiffB,diffB)
  !END DO; END DO; END DO; END DO ! i,j,k,iElem

  !WRITE(*,*)'SWAPB: min(B-curlA)',mindiffB
  !WRITE(*,*)'       max(B-curlA)',maxdiffB
  
  !WRITE(*,*)'CHECK B(U)...'
  !CALL CheckB(U_inout(6:8,:,:,:,:))
  !WRITE(*,*)'CHECK B=curlA...'
  !CALL CheckB(curlA_in)
  
  DO iElem=1,nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
    EmagU=SUM(U_inout(6:8,i,j,k,iElem)**2)
    EmagcurlA=SUM(curlA_in(:,i,j,k,iElem)**2)
    U_inout(5,i,j,k,iElem)=U_inout(5,i,j,k,iElem)+s2mu_0*(EmagcurlA-EmagU) !add new magn. energy, substract old energy
    U_inout(6:8,i,j,k,iElem)=curlA_in(:,i,j,k,iElem)
  END DO; END DO; END DO; END DO ! i,j,k,iElem

END SUBROUTINE swapB


!==================================================================================================================================
!> Replaces Bt of current Ut by curlAt, and change total energy Et to remain consistent with the entropy
!> inequality:
!> S_t=w^T Ut = ((gamma-s)/(gamma-1)-(rho/p)1/2|v|^2,(rho/p)vvec,-(rho/p),(rho/p)Bvec) ^T 
!>               (rho_t, (rho*vvec)_t,E_t,Bvec_t)
!>            =  ... (rho/p) (-E_t + Bvec*Bvec_t  + [ Bvec*(curlA_t-curlA_t) ] ) 
!>
!> in new Ut, the B_t is replaced with curlA_t, yielding
!>        =>  ... (rho/p) (-E*_t + Bvec*curlA_t)
!> if we choose:
!> E*_t =E_t + Bvec*(curlA_t-B_t)
!> we recover the original form, such that (w^T *U_t) remains the same!

!==================================================================================================================================
SUBROUTINE swapBt(U_in,Ut_inout)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Equation_Vars  ,ONLY: smu_0
USE MOD_Mesh_Vars      ,ONLY: nElems
USE MOD_CT_Vars      ,ONLY: curlAt
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: U_in(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Ut_inout(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k,iElem
REAL    :: EmagUt,EmagCurlAt
!==================================================================================================================================
  DO iElem=1,nElems; DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N 
    EmagUt=SUM(U_in(6:8,i,j,k,iElem)*Ut_inout(6:8,i,j,k,iElem))
    EmagcurlAt=SUM(U_in(6:8,i,j,k,iElem)*curlAt(:,i,j,k,iElem))
    Ut_inout(5,i,j,k,iElem)=Ut_inout(5,i,j,k,iElem)+smu_0*(EmagcurlAt-EmagUt) !to keep entropy inequality
    Ut_inout(6:8,i,j,k,iElem)=curlAt(:,i,j,k,iElem)
  END DO; END DO; END DO; END DO ! i,j,k,iElem

END SUBROUTINE swapBt

!==================================================================================================================================
!> checks divergence and surface jumps of cartesian vector field 
!> ONLY 1 PROC & STRUCTURED MESH 
!==================================================================================================================================
SUBROUTINE checkB(B_in)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Mesh_Vars      ,ONLY: nElems
USE MOD_Mesh_Vars      ,ONLY: sJ,metrics_ftilde,metrics_gtilde,metrics_htilde
USE MOD_Mesh_Vars      ,ONLY: ElemToSide,nSides
USE MOD_DG_Vars        ,ONLY: D 
!USE MOD_DG_Vars        ,ONLY: nTotal_IP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: B_in(3,0:PP_N,0:PP_N,0:PP_N,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!INTEGER :: i
REAL    :: Bn_plus(0:PP_N,0:PP_N,nSides),Bn_minus(0:PP_N,0:PP_N,nSides)
INTEGER :: i,j,k,l,iElem
REAL    :: mindivB,maxDivB,divB_loc,Btilde(3,0:PP_N,0:PP_N,0:PP_N)
!==================================================================================================================================

  mindivB=1.0e20
  maxdivB=-1.0e20
  Bn_plus=9999.
  Bn_minus=9999.
  DO iElem=1,nElems
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Btilde(1,i,j,k)=SUM(Metrics_ftilde(:,i,j,k,iElem)*B_in(:,i,j,k,iElem))
      Btilde(2,i,j,k)=SUM(Metrics_gtilde(:,i,j,k,iElem)*B_in(:,i,j,k,iElem))
      Btilde(3,i,j,k)=SUM(Metrics_htilde(:,i,j,k,iElem)*B_in(:,i,j,k,iElem))
    END DO; END DO; END DO ! i,j,k
    Bn_plus( :,:,ElemToSide(E2S_SIDE_ID,   XI_PLUS,iElem))=Btilde(1,PP_N,:,:)
    Bn_plus( :,:,ElemToSide(E2S_SIDE_ID,  ETA_PLUS,iElem))=Btilde(2,:,PP_N,:)
    Bn_plus( :,:,ElemToSide(E2S_SIDE_ID, ZETA_PLUS,iElem))=Btilde(3,:,:,PP_N)
    Bn_minus(:,:,ElemToSide(E2S_SIDE_ID,  XI_MINUS,iElem))=Btilde(1,0,:,:)
    Bn_minus(:,:,ElemToSide(E2S_SIDE_ID, ETA_MINUS,iElem))=Btilde(2,:,0,:)
    Bn_minus(:,:,ElemToSide(E2S_SIDE_ID,ZETA_MINUS,iElem))=Btilde(3,:,:,0)
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      divB_loc=0.
      DO l=0,PP_N
        divB_loc=divB_loc + D(i,l)*Btilde(1,l,j,k) + &
                            D(j,l)*Btilde(2,i,l,k) + &
                            D(k,l)*Btilde(3,i,j,l)
      END DO ! l
      divB_loc = sJ(i,j,k,iElem) * divB_loc 
      mindivB=MIN(mindivB,divB_loc)
      maxdivB=MAX(maxdivB,divB_loc)
  
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem

WRITE(*,*)'..check ..min/max divB    : ',mindivB,maxdivB
WRITE(*,*)'        ..min/max [[B*n]] : ',MINVAL((Bn_plus-Bn_minus)),MAXVAL((Bn_plus-Bn_minus))
!IF(MAXVAL(ABS(Bn_plus-Bn_minus)).GT.1.0e-10)THEN
!  DO i=1,nSides
!    WRITE(*,*)'SideID',i,'locSide master',SideToElem(S2E_LOC_SIDE_ID,i) &
!                        ,'locSide slave ',SideToElem(S2E_NB_LOC_SIDE_ID,i)               
!                                                 
!    WRITE(*,*)i,'     ..min/max [[B*n]] : ',MINVAL((Bn_plus(:,:,i)-Bn_minus(:,:,i))),MAXVAL((Bn_plus(:,:,i)-Bn_minus(:,:,i))),MAXVAL(Bn_plus(:,:,i))
!  END DO
!  READ(*,*)
!END IF

END SUBROUTINE checkB


!==================================================================================================================================
!> Finalizes global variables of the module.
!> Deallocate allocatable arrays, nullify pointers, set *InitIsDone = .FALSE.
!==================================================================================================================================
SUBROUTINE FinalizeCT()
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_CT_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
  SDEALLOCATE(CurlAt)
  SDEALLOCATE(CurlA)
  SDEALLOCATE(JacMat)
  SDEALLOCATE(Vdm_GLN_int)
  SDEALLOCATE(x_int)
  SDEALLOCATE(w_int)
  SDEALLOCATE(Vdm_e_GLN)
  SDEALLOCATE(Vdm_phi_GLN)
  SDEALLOCATE(proj_int_e)
  SDEALLOCATE(proj_int_phi)
  SDEALLOCATE(covVec1_int)
  SDEALLOCATE(covVec2_int)
  SDEALLOCATE(covVec3_int)
  SDEALLOCATE(Elem_xint)
  CTInitIsDone = .FALSE.
END SUBROUTINE FinalizeCT

END MODULE MOD_CT
