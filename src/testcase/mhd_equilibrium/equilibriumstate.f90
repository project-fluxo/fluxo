!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
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


MODULE MOD_EquilibriumState
!==================================================================================================================================
!
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================

INTERFACE InitEquilibriumState
  MODULE PROCEDURE InitEquilibriumState 
END INTERFACE

PUBLIC::InitEquilibriumState
!==================================================================================================================================

CONTAINS


SUBROUTINE InitEquilibriumState()
!==================================================================================================================================
! store initial gradients of the magnetic field 
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_MPI
USE MOD_MPI_Vars
USE MOD_Testcase_Vars
USE MOD_Testcase_ExactFunc ,ONLY: TestcaseExactFunc
USE MOD_CalcTimeStep       ,ONLY: CalcTimeStep
USE MOD_DG                 ,ONLY: DGTimeDerivative
USE MOD_DG_Vars            ,ONLY: DGInitIsDone
USE MOD_DG_Vars            ,ONLY: U,Ut,U_Minus,U_Plus
USE MOD_ProlongToFace      ,ONLY: ProlongToFace
USE MOD_Equation           ,ONLY: ExactFunc
USE MOD_Equation_Vars      ,ONLY: IniExactFunc
USE MOD_Equation_Vars      ,ONLY: nBCByType,BCSideID,BCdata
USE MOD_Equation_Vars      ,ONLY: ConsToPrimVec,PrimToConsVec
USE MOD_Equation_Vars      ,ONLY: StrVarNames
USE MOD_Restart_Vars       ,ONLY: DoRestart
USE MOD_Mesh_Vars          ,ONLY: MeshFile
USE MOD_Mesh_Vars          ,ONLY: nBCSides,nBCs,BoundaryType
USE MOD_Mesh_Vars          ,ONLY: Elem_xGP
USE MOD_Mesh_Vars          ,ONLY: NormVec,TangVec1,TangVec2
USE MOD_Mesh_Vars          ,ONLY: nBCSides,nInnerSides,nMPISides_MINE
USE MOD_Mesh_Vars          ,ONLY: SideID_plus_lower,SideID_plus_upper
USE MOD_Output             ,ONLY: VisualizeAny
USE MOD_HDF5_Output        ,ONLY: WriteAnyStateToHDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i,j,k,iElem
INTEGER             :: iBC,iSide,SideID,BCType,nBCLoc
REAL                :: dt_min
REAL                :: resu_t(PP_nVar),resu_tt(PP_nVar)
REAL                :: U_tmp(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
REAL,ALLOCATABLE    :: Apot(:,:,:,:,:),Bdivfree(:,:,:,:,:)
INTEGER             :: errType
INTEGER             :: nTotal
INTEGER             :: EquilibriumDisturbFunc
CHARACTER(LEN=255)  :: errMsg
CHARACTER(LEN=255)  :: FileTypeStr
CHARACTER(LEN=255)  :: EquilibriumStateFile
LOGICAL             :: EquilibriumDivBcorr
!CHECK
REAL                :: deltaB(3),maxjmp_B(3)
!==================================================================================================================================
IF(.NOT.DGInitIsDone)THEN
   CALL abort(__STAMP__, &
   'InitEquilibriumState not ready to be called or already called.',999,999.)
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT EQUILIBRIUM STATE...'


ALLOCATE(Ueq(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
ALLOCATE(Ueq_BC(PP_nVar,0:PP_N,0:PP_N,nBCSides))
ALLOCATE(Uteq(PP_nVar,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
Ueq  =0.
Ueq_BC=0.
Uteq =0.
IF(EquilibriumStateIni.EQ.0) THEN
  EquilibriumStateIni=IniExactFunc
  SWRITE(UNIT_StdOut,'(A,A33,A3,I22)') ' | EquilibriumStateIni changed!   | ', &
                                             'set to IniExactFunc'    ,' | ',EquilibriumStateIni
END IF
IF(EquilibriumDivBcorr) THEN
  ALLOCATE(Apot(3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  Apot=-999.
END IF

IF(EquilibriumStateIni.GT.0)THEN
  evalEquilibrium=.TRUE. !for exactfunctestcase
  IF(EquilibriumDivBcorr)THEN
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N;  DO i=0,PP_N
        CALL TestcaseExactFunc(EquilibriumStateIni,0., &
                               Elem_xGP(1:3,i,j,k,iElem),Ueq(1:PP_nVar,i,j,k,iElem), &
                               resu_t,resu_tt,Apot(:,i,j,k,iElem))
      END DO; END DO; END DO!i,j,k
    END DO ! iElem=1,PP_nElems
  ELSE
    DO iElem=1,PP_nElems
      DO k=0,PP_N; DO j=0,PP_N;  DO i=0,PP_N
        CALL ExactFunc(EquilibriumStateIni,0.,Elem_xGP(1:3,i,j,k,iElem),Ueq(1:PP_nVar,i,j,k,iElem))
      END DO; END DO; END DO!i,j,k
    END DO ! iElem=1,PP_nElems
  END IF !EquilibriumDivBcorr
ELSEIF(EquilibriumStateIni.EQ.-2) THEN
  !Read Ueq from Mesh
  IF(EquilibriumDivBcorr)THEN
    CALL ReadEquilibriumFromMesh(Ueq,Apot)
  ELSE
    CALL ReadEquilibriumFromMesh(Ueq)
  END IF
ELSEIF(EquilibriumStateIni.EQ.-3) THEN
  !EquilibriumStateFile=GETSTR('EquilibriumStateFile')
  EquilibriumDivBcorr=.FALSE. !Apot not there
  CALL ReadEquilibriumFromMesh(Ueq)
  CALL ReadEquilibriumFromState(EquilibriumStateFile,Ueq) !overwrite Ueq with U from Statefile
END IF !EquilibriumState

IF(EquilibriumDivBcorr)THEN
  ALLOCATE(Bdivfree(3,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  !==== overwrite B FROM POTENTIAL
  CALL ComputeBfromPotential(Apot,Ueq(6:8,:,:,:,:),Bdivfree)
  nTotal=(PP_N+1)**3*PP_nElems
  !PrimToCons before overwriting magnetic field, to keep same pressure!!
  CALL ConsToPrimVec(nTotal,U_tmp,Ueq)
  U_tmp(6:8,:,:,:,:)=BdivFree(1:3,:,:,:,:)
  CALL PrimToConsVec(nTotal,U_tmp,Ueq)
  DEALLOCATE(Apot,Bdivfree)

END IF !EquilibriumDivBcorr

!CHECK continuity of tangential B-field on FACES
CALL ProlongToFace(Ueq,U_Minus,U_Plus,doMPISides=.FALSE.)
#ifdef MPI
! Prolong to face for MPI sides - send direction
CALL StartReceiveMPIData(U_Plus,DataSizeSide,SideID_plus_lower,SideID_plus_upper,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE
CALL ProlongToFace(Ueq,U_Minus,U_Plus,doMPiSides=.TRUE.)
CALL StartSendMPIData(   U_Plus,DataSizeSide,SideID_plus_lower,SideID_plus_upper,MPIRequest_U(:,RECV),SendID=2) ! Send YOUR
! Complete send / receive
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U) !Send YOUR - receive MINE
#endif /*MPI*/
maxjmp_B(:)=0.
DO SideID=nBCsides+1,nBCsides+nInnerSides+nMPISides_MINE
  DO j=0,PP_N ; DO i=0,PP_N
    deltaB(:) = U_Plus(6:8,i,j,SideID)-U_Minus(6:8,i,j,SideID)
    maxjmp_B(1)  = MAX(maxjmp_B(1),ABS(SUM(NormVec( :,i,j,SideID)*deltaB(:))))
    maxjmp_B(2)  = MAX(maxjmp_B(2),ABS(SUM(TangVec1(:,i,j,SideID)*deltaB(:))))
    maxjmp_B(3)  = MAX(maxjmp_B(3),ABS(SUM(TangVec2(:,i,j,SideID)*deltaB(:))))
  END DO; END DO !i,j
END DO !SideID
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,maxjmp_B  ,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  WRITE(*,'(A,3E21.11)') 'CHECK MAX jmp_Bn,jmp_Bt1,jmp_Bt2',maxjmp_B(:)
ELSE
  CALL MPI_REDUCE(maxjmp_B  ,0,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
END IF

!save U to U_tmp for restart
IF(doRestart) U_tmp=U
!overwrite initial U for dg_timederivative 
U=Ueq

!set state state Dirichlet BC 21: equilibrium state evaluated at the Boundary!
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  IF(BCType.NE.21) CYCLE
  nBCLoc =nBCByType(iBC)
  IF(nBCLoc.EQ.0) CYCLE
  ! FOR BCType 21, use equilibrium state as BC
  CALL ProlongToFace(U,U_Minus,U_Plus,doMPISides=.FALSE.)
  DO iSide=1,nBCLoc
    SideID=BCSideID(iBC,iSide)
    BCdata(:,:,:,SideID)=U_Minus(:,:,:,SideID)
  END DO !iSide=1,nBCloc
  U_Minus=0.
  U_Plus=0.
END DO !iBC=1,nBCs

! Call DG operator to fill face data, fluxes, gradients for analyze
dt_Min=CALCTIMESTEP(errType,errMsg) 
CALL DGTimeDerivative(0.)

!store time derivative of Ueq
Uteq= Ut 
Ut=0.
!check time derivative norms
CALL EvalNorms(Uteq)


IF(doRestart)THEN
  U=U_tmp !copy U back
ELSE
  EquilibriumDisturbFunc=GETINT('EquilibriumDisturbFunc','0')
  IF(EquilibriumDisturbFunc.EQ.0) THEN
    EquilibriumDisturbFunc=IniExactFunc
    SWRITE(UNIT_StdOut,'(A,A33,A3,I22)') ' | EquilibriumDisturbFunc changed!| ', &
                                               'set to IniExactFunc'    ,' | ',EquilibriumDisturbFunc
  END IF
  !then disturb the Initial state U
  CALL DisturbU(EquilibriumDisturbFunc)
  FileTypeStr='EquilibriumSource'
  CALL VisualizeAny(Uteq,FileTypeStr,0.,PrimVisu=.FALSE.)
  CALL WriteAnyStateToHDF5(PP_nVar,Uteq,MeshFile,FileTypeStr,StrVarNames,0.)
  FileTypeStr='EquilibriumState'
  CALL VisualizeAny(Ueq,FileTypeStr,0.,PrimVisu=.TRUE.)
  CALL WriteAnyStateToHDF5(PP_nVar,Ueq,MeshFile,FileTypeStr,StrVarNames,0.)
END IF !doRestart

SWRITE(UNIT_stdOut,'(A)')' INIT EQUILIBRIUM STATE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquilibriumState


SUBROUTINE ComputeBfromPotential(Apot,Bexact,Bdivfree)
!==================================================================================================================================
! Calculates B=curl(A) locally for one element
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars, ONLY: sJ,Metrics_ftilde,Metrics_gtilde,Metrics_htilde
USE MOD_Mesh_Vars, ONLY: Elem_xGP
USE MOD_DG_Vars,   ONLY: D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN   )              :: Apot(3,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
REAL,INTENT(IN   )              :: Bexact(3,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: Bdivfree(3,0:PP_N,0:PP_N,0:PP_N,PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: i,j,k,l,iElem
REAL                            :: dAdxi(3),dAdeta(3),dAdzeta(3)
REAL                            :: B_f,B_g,B_h
REAL                            :: aCov(3,0:PP_N,0:PP_N,0:PP_N,3)  ! last index a_1, a_2, a_3, first index x,y,z
REAL                            :: APotCov(3,0:PP_N,0:PP_N,0:PP_N)
REAL                            :: divB(0:PP_N,0:PP_N,0:PP_N)
REAL                            :: maxDivB,maxBdiff(3),maxBex(3)
REAL                            :: x(3) 
INTEGER                         :: method
!==================================================================================================================================
SWRITE(UNIT_stdOut,'(A)') ' COMPUTE B FROM POTENTIAL...'

!!TEST !!!!!!!!!!!!!!
!! OVERWRITE Apot and Bexact
!DO iElem=1,PP_nElems
!  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
!      x(:)=Elem_xGP(:,i,j,k,iElem)-1.
!!      Apot(1,i,j,k,iElem)  =0.
!!      Apot(2,i,j,k,iElem)  =0.
!!      Apot(3,i,j,k,iElem)  =-0.5*(x(1)*x(1)+x(2)*x(2)) 
!!      Bexact(1,i,j,k,iElem)=-x(2)                   !=d/dy Az - d/dz Ay
!!      Bexact(2,i,j,k,iElem)= x(1)                   !=d/dz Ax - d/dx Az
!!      Bexact(3,i,j,k,iElem)= 0.                     !=d/dx Ay - d/dy Ax
!      Apot(1,i,j,k,iElem)  =x(3)
!      Apot(2,i,j,k,iElem)  =x(1)
!      Apot(3,i,j,k,iElem)  =-x(2)
!      Bexact(1,i,j,k,iElem)=-1.                  !=d/dy Az - d/dz Ay
!      Bexact(2,i,j,k,iElem)= 1.                  !=d/dz Ax - d/dx Az
!      Bexact(3,i,j,k,iElem)= 1.                  !=d/dx Ay - d/dy Ax
!!      Apot(1,i,j,k,iElem)  =1.
!!      Apot(2,i,j,k,iElem)  =2.
!!      Apot(3,i,j,k,iElem)  =3.
!!      Bexact(1,i,j,k,iElem)=0.                  !=d/dy Az - d/dz Ay
!!      Bexact(2,i,j,k,iElem)=0.                  !=d/dz Ax - d/dx Az
!!      Bexact(3,i,j,k,iElem)=0.                  !=d/dx Ay - d/dy Ax
!  END DO; END DO; END DO!i,j,k
!END DO !iElem
!!TEST !!!!!!!!!!!!!!


Bdivfree=0.
DO iElem=1,PP_nElems
  !compute covariant vectors
  aCov=0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      aCov(:,i,j,k,1) = aCov(:,i,j,k,1) +D(i,l)*Elem_xGP(:,l,j,k,iElem)
      aCov(:,i,j,k,2) = aCov(:,i,j,k,2) +D(j,l)*Elem_xGP(:,i,l,k,iElem)
      aCov(:,i,j,k,3) = aCov(:,i,j,k,3) +D(k,l)*Elem_xGP(:,i,j,l,iElem)
    END DO !l
  END DO; END DO; END DO!i,j,k
!  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
!    aCov(:,i,j,k,1) = sJ(i,j,k,iElem)*CROSS(Metrics_gtilde(:,i,j,k,iElem),Metrics_htilde(:,i,j,k,iElem))
!    aCov(:,i,j,k,2) = sJ(i,j,k,iElem)*CROSS(Metrics_htilde(:,i,j,k,iElem),Metrics_ftilde(:,i,j,k,iElem))
!    aCov(:,i,j,k,3) = sJ(i,j,k,iElem)*CROSS(Metrics_ftilde(:,i,j,k,iElem),Metrics_gtilde(:,i,j,k,iElem))
!  END DO; END DO; END DO!i,j,k
  
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    ApotCov(1,i,j,k)=SUM(Apot(:,i,j,k,iElem)*aCov(:,i,j,k,1))
    ApotCov(2,i,j,k)=SUM(Apot(:,i,j,k,iElem)*aCov(:,i,j,k,2))
    ApotCov(3,i,j,k)=SUM(Apot(:,i,j,k,iElem)*aCov(:,i,j,k,3))
  END DO; END DO; END DO!i,j,k
  
  ! Compute the covariant components JB^i 
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    dAdxi=0.
    dAdeta=0.
    dAdzeta=0.
    DO l=0,PP_N
      dAdxi(:)   = dAdxi(:)   +D(i,l)*ApotCov(:,l,j,k)
      dAdeta(:)  = dAdeta(:)  +D(j,l)*ApotCov(:,i,l,k)
      dAdzeta(:) = dAdzeta(:) +D(k,l)*ApotCov(:,i,j,l)
    END DO !l
    Bdivfree(1,i,j,k,iElem)=(dAdeta( 3)-dAdzeta(2))
    Bdivfree(2,i,j,k,iElem)=(dAdzeta(1)-dAdxi(  3))
    Bdivfree(3,i,j,k,iElem)=(dAdxi(  2)-dAdeta( 1))
    !cartesian components 1/J(sum_i JB^i*a_i))
    Bdivfree(:,i,j,k,iElem)= (  Bdivfree(1,i,j,k,iElem)*aCov(:,i,j,k,1) &
                              + Bdivfree(2,i,j,k,iElem)*aCov(:,i,j,k,2) &
                              + Bdivfree(3,i,j,k,iElem)*aCov(:,i,j,k,3) )*sJ(i,j,k,iElem)
  END DO; END DO; END DO!i,j,k
  
END DO !iElem


maxDivB=0.
!CHECK STRONG DIVERGENCE div B = 1/J sum_i d/dxi^i (JB^i)
DO iElem=1,PP_nElems
  divB=0.
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    ! Compute the covariant components 
    B_f = SUM(Metrics_fTilde(:,i,j,k,iElem)*Bdivfree(:,i,j,k,iElem))
    B_g = SUM(Metrics_gTilde(:,i,j,k,iElem)*Bdivfree(:,i,j,k,iElem))
    B_h = SUM(Metrics_hTilde(:,i,j,k,iElem)*Bdivfree(:,i,j,k,iElem))
    DO l=0,PP_N
      divB(l,j,k) = divB(l,j,k)+D(l,i)*B_f
      divB(i,l,k) = divB(i,l,k)+D(l,j)*B_g
      divB(i,j,l) = divB(i,j,l)+D(l,k)*B_h
    END DO !l
  END DO; END DO; END DO!i,j,k
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    divB(i,j,k) = divB(i,j,k)*sJ(i,j,k,iElem) 
  END DO; END DO; END DO!i,j,k
  maxDivB    =MAX(maxDivB,MAXVAL(ABS(divB)))
END DO !iElem
maxBex(1)=MAXVAL(ABS(Bexact(1,:,:,:,:)))
maxBex(2)=MAXVAL(ABS(Bexact(2,:,:,:,:)))
maxBex(3)=MAXVAL(ABS(Bexact(3,:,:,:,:)))
maxBdiff(1)=MAXVAL(ABS(Bdivfree(1,:,:,:,:)-Bexact(1,:,:,:,:)))
maxBdiff(2)=MAXVAL(ABS(Bdivfree(2,:,:,:,:)-Bexact(2,:,:,:,:)))
maxBdiff(3)=MAXVAL(ABS(Bdivfree(3,:,:,:,:)-Bexact(3,:,:,:,:)))
#ifdef MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,maxBex  ,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,maxBdiff,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,maxDivB ,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  WRITE(*,'(A,3E21.11)')'   -> max|B-Bini| x,y,z',maxBdiff
  WRITE(*,'(A,3E21.11)')'   -> max|B-Bini|/max(Bini) x,y,z',maxBdiff/maxBex
  WRITE(*,'(A,E21.11)') '   -> Max. element-local DIV B',maxDivB
  WRITE(UNIT_stdOut,'(A)')' ... DONE!'
ELSE
  CALL MPI_REDUCE(maxBex  ,0,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(maxBdiff,0,3,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(maxDivB ,0,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/

END SUBROUTINE ComputeBfromPotential


SUBROUTINE DisturbU(EquilibriumDisturbFunc) 
!==================================================================================================================================
! add a disturbance to the solution U, depending on the EquilibriumDisturbFunc 
!==================================================================================================================================
! MODULES
USE MOD_Preproc
USE MOD_Globals               ,ONLY: Abort
USE MOD_Testcase_Vars         ,ONLY: InputEq
USE MOD_DG_Vars               ,ONLY: U
USE MOD_Mesh_Vars             ,ONLY: Elem_xGP
USE MOD_Equation_Vars         ,ONLY: Pi,sSqrt4Pi,Kappa,sKappaM1,AdvVel,RefStateCons,RefStatePrim,IniRefState
USE MOD_Equation_Vars         ,ONLY: smu_0,mu_0
USE MOD_Equation_Vars         ,ONLY: IniCenter,IniFrequency,IniAmplitude,IniHalfwidth,IniWaveNumber
USE MOD_Equation_Vars         ,ONLY: IniDisturbance
USE MOD_Equation_Vars         ,ONLY: PrimToCons,ConsToPrim
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)              :: EquilibriumDisturbFunc   
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: Prim(1:PP_nVar),x(3)
INTEGER                         :: i,j,k,l,m,iElem
REAL                            :: a,vtor,psinorm
!==================================================================================================================================
! Determine the value, the first and the second time derivative
SELECT CASE (EquilibriumDisturbFunc)
CASE(1) !no disturbance
  !resu=resu
CASE(10071) !Tearing mode instability, of paper Landi et al. , domain [0,6*pi]x[-pi/2,pi/2]x[0:2Pi]
        ! "Three-dimensional simulations of compressible tearing instability"
        ! rho_0=1, p0=0.5*beta (choose with refstate)
        ! Re_eta=5000, mu=0.,kappa=5/3 1/delta=0.1(=IniHalfwidth)  IniAmplitude=1.0E-04
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N;  DO i=0,PP_N
      x=Elem_xGP(:,i,j,k,iElem)
      CALL ConsToPrim(Prim,U(:,i,j,k,iElem))
        DO m=0,NINT(IniWaveNumber(3))
          DO l=0,NINT(IniWaveNumber(1))
            a=(1.+0.8*REAL(l)+0.9*REAL(m))/(1+0.8*IniWaveNumber(1)+0.9*IniWaveNumber(3))
            Prim(3)=Prim(3)+SIN(x(1)/3.*REAL(l)+ x(3)*REAL(m)+2*PP_Pi*a)
          END DO
        END DO
        Prim(3)=IniDisturbance*Prim(6)*Prim(8)*Prim(3)
      !Prim(6:8)=sSqrt4pi*Prim(6:8) ! scaling with sqrt(4pi)!?!
      CALL PrimToCons(Prim,U(:,i,j,k,iElem))
    END DO; END DO; END DO !i,j,k
  END DO ! iElem=1,PP_nElems

CASE(10090,10091) !cylindrical equilibrium for ideal MHD for current hole (Czarny, JCP, 2008), current Jz in z direction is given:
         ! cylindrical domain r[0,1], z[0,20] (from q(r=1)=4.4 =2*pi*B0/(Lz*Bphi(r=1)) => B0/Lz=0.364 B0~7.44, L0~20.)
         ! Jz=j1*(1-r^4)-j2*(1-r^2)^8, j1=0.2, j2=0.266
         ! from J=rotB (Br=0) follows
         ! Bphi(r) =mu_0 1/r  \int_0^r r*Jc dr
         ! pressure difference from gradp=J x B:
         ! dp(r) = -smu_0 ( 0.5*Bphi^2 + \int_0^r Bphi^2/r dr)
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N;  DO i=0,PP_N
      x=Elem_xGP(:,i,j,k,iElem)
      CALL ConsToPrim(Prim,U(:,i,j,k,iElem))
      Prim(2)=IniDisturbance*(1.-SUM(x(1:2)**2))**2
      Prim(3)=x(2)*Prim(2)
                  
      CALL PrimToCons(Prim,U(:,i,j,k,iElem))
    END DO; END DO; END DO !i,j,k
  END DO ! iElem=1,PP_nElems
CASE(101) ! R=10, a=1, 1st toroidal n=1&2 mode of toroidal velocity, 
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N;  DO i=0,PP_N
      x=Elem_xGP(:,i,j,k,iElem)
      CALL ConsToPrim(Prim,U(:,i,j,k,iElem))
      a=ATAN2(x(2),x(1)) !toroidal angle
      !psinorm=SQRT(InputEq(6,i,j,k,iElem))
      !vtor=IniDisturbance*(4.*(1-psinorm)*psinorm)**3*COS(a)
      psinorm=InputEq(6,i,j,k,iElem)
      vtor=IniDisturbance*(1.-psinorm)*psinorm*(COS(a)+0.5*SIN(2*a))
      prim(2)= -vtor*SIN(a)
      prim(3)=  vtor*COS(a)
      prim(4)= 0.
      CALL PrimToCons(Prim,U(:,i,j,k,iElem))
    END DO; END DO; END DO !i,j,k
  END DO ! iElem=1,PP_nElems
CASE(999) ! TESTING
  DO iElem=1,PP_nElems
    DO k=0,PP_N; DO j=0,PP_N;  DO i=0,PP_N
      x=Elem_xGP(:,i,j,k,iElem)
      CALL ConsToPrim(Prim,U(:,i,j,k,iElem))
      a=ATAN2(x(2),x(1)) !toroidal angle
      !vtor=(10.+2*COS(a)+0.03*SIN(2*a))/SQRT(Prim(1))
      vtor=(10+3*cos(a-0.3)+4.*sin(a+0.5)+0.3*cos(2*a-0.1) + 0.4*sin(2*a+0.3) + 0.01*cos(3*a+0.25) -0.04*sin(3*a-0.7))/SQRT(prim(1)
      prim(2)= -vtor*SIN(a)
      prim(3)=  vtor*COS(a)
      prim(4)= 0.
      CALL PrimToCons(Prim,U(:,i,j,k,iElem))
    END DO; END DO; END DO !i,j,k
  END DO ! iElem=1,PP_nElems
CASE DEFAULT
  CALL abort(__STAMP__,' DisturbU function not specified!')
END SELECT ! ExactFunction
END SUBROUTINE DisturbU


SUBROUTINE EvalNorms(Uin)
!==================================================================================================================================
! Calculates L_infinfity and L_2 norms of Ut 
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: sJ
USE MOD_Interpolation_Vars, ONLY: wGP
USE MOD_Equation_Vars,      ONLY: StrVarNames
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                 :: Uin(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                            :: L_2(PP_nVar),L_Inf(PP_nVar)
INTEGER                         :: iElem,i,j,k
REAL                            :: IntegrationWeight,vol
!==================================================================================================================================
! Calculate error norms
L_Inf(:)=-1.E10
L_2(:)=0.
vol=0.
! Interpolate values of Error-Grid from GP's
DO iElem=1,PP_nElems
   DO k=0,PP_N
     DO j=0,PP_N
       DO i=0,PP_N
         L_Inf = MAX(L_Inf,abs(Uin(:,i,j,k,iElem)))
         IntegrationWeight = wGP(i)*wGP(j)*wGP(k)/sJ(i,j,k,iElem)
         ! To sum over the elements, We compute here the square of the L_2 error
         L_2 = L_2+Uin(:,i,j,k,iElem)**2*IntegrationWeight
         vol = vol+IntegrationWeight
       END DO ! k
     END DO ! l
   END DO ! m
END DO ! iElem=1,PP_nElems

#ifdef MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,L_2  ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,L_Inf,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(MPI_IN_PLACE,vol  ,      1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
ELSE
  CALL MPI_REDUCE(L_2  ,0           ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(L_Inf,0           ,PP_nVar,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,iError)
  CALL MPI_REDUCE(vol  ,0           ,      1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,iError)
END IF
#endif /*MPI*/

! We normalize the L_2 norm with the Volume of the domain and take into account that we have to use the square root
IF(MPIroot)THEN
    L_2 = SQRT(L_2/Vol)
    WRITE(UNIT_StdOut,'(66("-"))')
    WRITE(UNIT_StdOut,'(A20,2(2X,A21))')'Variablename: ','L_2  (Ut_eq): ',' L_inf(Ut_eq): '
    DO i=1,PP_nVar
      WRITE(UNIT_StdOut,'(A20,2(2X,ES21.15))')TRIM(StrVarNames(i))//' : ',L_2(i),L_inf(i)
    END DO
    WRITE(UNIT_StdOut,'(66("-"))')
END IF
END SUBROUTINE EvalNorms


SUBROUTINE ReadEquilibriumFromMesh(Ueq,Apot)
!==================================================================================================================================
! READ equilibrium solution from mesh file 
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars             ,ONLY: MeshFile,Ngeo
USE MOD_Mesh_Vars             ,ONLY: offsetElem
USE MOD_Interpolation_Vars    ,ONLY: NodeTypeGL
USE MOD_Testcase_Vars         ,ONLY: InputEq
USE MOD_Interpolation         ,ONLY: GetVandermonde 
USE MOD_Equation_Vars         ,ONLY: PrimToConsVec
USE MOD_HDF5_input            ,ONLY: OpenDataFile,CloseDataFile,ReadArray,DataSetExists,GetDataSize
USE MOD_HDF5_input            ,ONLY: File_ID
USE MOD_ChangeBasis           ,ONLY: ChangeBasis3D
USE MOD_IO_HDF5               ,ONLY: Hsize,nDims
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: Ueq(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
REAL,INTENT(OUT),OPTIONAL :: Apot(3,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,ALLOCATABLE   :: MHDEQdata_GL(:,:,:,:,:)
REAL               :: Vdm_GLNgeo_N(0:PP_N,0:Ngeo)
REAL               :: maxPsi(2),minPsi(2)
INTEGER            :: nTotal,iElem,nVarMHDEQ
LOGICAL            :: exists
!==================================================================================================================================
  CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.)
  CALL DataSetExists(File_ID,'MHDEQdata_GL',exists)
  IF(exists)THEN 
    CALL GetDataSize(File_ID,'MHDEQdata_GL',nDims,HSize)
    nVarMHDEQ=HSize(1)
    DEALLOCATE(HSize)
    IF((nVarMHDEQ.NE.7).AND.(nVarMHDEQ.NE.10)) THEN 
      CALL abort(__STAMP__, &
      'ERROR: MHDEQdata_GL has not the right size (7/10) in Meshfile: '//TRIM(Meshfile))
    END IF
    IF((nVarMHDEQ.EQ.7).AND.(PRESENT(Apot)) ) THEN 
      CALL abort(__STAMP__, &
      'ERROR: MHDEQdata_GL has not the magnetic potential included, in Meshfile: '//TRIM(Meshfile))
    END IF
    ALLOCATE(MHDEQdata_GL(nVarMHDEQ,0:Ngeo,0:Ngeo,0:Ngeo,PP_nElems))
    CALL ReadArray('MHDEQdata_GL',5,(/nVarMHDEQ,Ngeo+1,Ngeo+1,Ngeo+1,PP_nElems/),OffsetElem,5,RealArray=MHDEQdata_GL)
  END IF
  CALL CloseDataFile() 
  IF(.NOT.exists) CALL Abort(__STAMP__, &
          'ERROR: MHDEQdata_GL does not exist in Meshfile: '//TRIM(Meshfile))
  !INPUT variables, on XGL_Ngeo points: 
  ! 1 Density
  ! 2 Pressure
  ! 3-5 (Bx,By,Bz)
  ! 6 normalized poloidal flux
  ! 7 normalized toroidal flux
  ! 8-10 (Ax,Ay,Az) 

  !interpolate to computing nodes
  ALLOCATE(InputEq(1:nVarMHDEQ,0:PP_N,0:PP_N,0:PP_N,PP_nElems))
  CALL GetVandermonde(Ngeo,NodeTypeGL,PP_N,NodeType,Vdm_GLNgeo_N)
  DO iElem=1,PP_nElems
    CALL ChangeBasis3D(nVarMHDEQ,Ngeo,PP_N,Vdm_GLNgeo_N,MHDEQdata_GL(:,:,:,:,iElem),InputEq(:,:,:,:,iElem))
  END DO !iElem
  Ueq=0.
  Ueq(  1,:,:,:,:) =InputEq(  1,:,:,:,:) !density
  Ueq(  5,:,:,:,:) =InputEq(  2,:,:,:,:) !pressure
  Ueq(6:8,:,:,:,:) =InputEq(3:5,:,:,:,:) !magnetic field (x,y,z components)
  nTotal=(PP_N+1)**3*PP_nElems
  ! input data is in primitive variables, change to cons 
  CALL PrimToConsVec(nTotal,Ueq,Ueq)
  !if potential is given
  IF(PRESENT(Apot))THEN
    Apot(:,:,:,:,:) =InputEq(8:10,:,:,:,:) !magnetic potential (x,y,z components)
  END IF !Apot present
  DEALLOCATE(MHDEQdata_GL)
  !Renormalize poloidal & toroidal flux
  maxPsi(1)=MAXVAL(InputEq(6,:,:,:,:))
  maxPsi(2)=MAXVAL(InputEq(7,:,:,:,:))
  minPsi(1)=MINVAL(InputEq(6,:,:,:,:))
  minPsi(2)=MINVAL(InputEq(7,:,:,:,:))
#ifdef MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxPsi,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,iError)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,minPsi,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#endif /*MPI*/
  InputEq(6,:,:,:,:)= (InputEq(6,:,:,:,:)-minPsi(1))/(maxPsi(1)-minPsi(1))
  InputEq(7,:,:,:,:)= (InputEq(7,:,:,:,:)-minPsi(2))/(maxPsi(2)-minPsi(2))
  
END SUBROUTINE ReadEquilibriumFromMesh 


SUBROUTINE ReadEquilibriumFromState(StateFile,Ueq)
!==================================================================================================================================
! READ equilibrium solution from mesh file 
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars             ,ONLY: offsetElem
USE MOD_HDF5_input            ,ONLY: OpenDataFile,CloseDataFile,ReadArray,DataSetExists,GetDataSize
USE MOD_HDF5_input            ,ONLY: File_ID
USE MOD_IO_HDF5               ,ONLY: Hsize,nDims
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=255)  :: StateFile
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Ueq(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL               :: U_HDF5(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:PP_nElems)
INTEGER            :: nVarHDF5,N_HDF5
LOGICAL            :: exists
!==================================================================================================================================
  CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.)
  CALL DataSetExists(File_ID,'DG_Solution',exists)
  IF(exists)THEN 
    CALL GetDataSize(File_ID,'DG_Solution',nDims,HSize)
    nVarHDF5=HSize(1)
    N_HDF5=HSize(2)-1
    DEALLOCATE(HSize)
    IF(nVarHDF5.NE.PP_nVar) CALL abort(__STAMP__, &
          'ERROR: nVar of DG_Solution /=PP_nVar in Statefile: '//TRIM(StateFile))
    IF(N_HDF5.NE.PP_N) CALL abort(__STAMP__, &
          'ERROR: N of DG_Solution /=PP_N in Statefile: '//TRIM(StateFile))
    CALL ReadArray('DG_Solution',5,(/PP_nVar,PP_N+1,PP_N+1,PP_N+1,PP_nElems/),OffsetElem,5,RealArray=U_HDF5)
  END IF
  CALL CloseDataFile() 
  IF(.NOT.exists) CALL Abort(__STAMP__, &
          'ERROR: DG_Solution does not exist in Statefile: '//TRIM(Statefile))
  Ueq(1:8,:,:,:,:)=U_HDF5(1:8,:,:,:,:)
  Ueq(9,:,:,:,:)=0. !reset psi 
END SUBROUTINE ReadEquilibriumFromState 

END MODULE MOD_EquilibriumState
