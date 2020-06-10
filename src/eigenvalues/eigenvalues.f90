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
!> Contains the eigenvalue intialization and output routines 
!>
!===================================================================================================================================
MODULE MOD_EigenValues
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
INTERFACE DefineParametersEigenValues
  MODULE PROCEDURE DefineParametersEigenValues
END INTERFACE

INTERFACE InitEigenValues
  MODULE PROCEDURE InitEigenValues
END INTERFACE

INTERFACE EigenValues
  MODULE PROCEDURE EigenValues
END INTERFACE

INTERFACE FinalizeEigenValues
  MODULE PROCEDURE FinalizeEigenValues
END INTERFACE

PUBLIC  :: DefineParametersEigenValues
PUBLIC  :: InitEigenValues
PUBLIC  :: EigenValues
PUBLIC  :: FinalizeEigenValues
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!>
!==================================================================================================================================
SUBROUTINE DefineParametersEigenValues()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("EigenValues")
CALL prms%CreateLogicalOption('calcEigVals', &
     "switch on computation of eigenvalues", 'F')
CALL prms%CreateIntOption('stableDtfromEV', &
   "=-1: off, =0: find stable timestep from Eigenvalues & charact. polynomial of time integator, =1: apply new timestep", '0')
CALL prms%CreateIntOption('EV_reducedDim', &
     " leave out 0:PP_N in higher dimensions: 1: only xi (1D grid in x), 2: xi,eta (2D grid in x,y),3: 3D", '3')
CALL prms%CreateIntOption('EV_maxIter', &
     "Max. Iteration count of Arnoldi -1: set to Problem size, >0 compute only maxIter number of eigenvalues", '-1')
CALL prms%CreateIntOption(     "FD_degree"  , " polynomial degree for approximation of Finite Diff. of Matvec (d/dU R(U))*v")
END SUBROUTINE DefineParametersEigenValues

!===================================================================================================================================
!> Initializes all module variables
!>
!===================================================================================================================================
SUBROUTINE InitEigenValues()
! MODULES
USE MOD_Globals                                           ! UNIT_stdOut
USE MOD_PreProc                                           ! all PP_*** variables
USE MOD_ReadInTools,ONLY:GETLOGICAL,GETINT
USE MOD_EigenValues_Vars
USE MOD_Mesh_Vars, ONLY: nGlobalElems,nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER ::EV_reducedDim
!===================================================================================================================================
  SWRITE(UNIT_StdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,'(A)') ' INIT EIGENVALUES...'

  calcEigVals=GETLOGICAL('calcEigVals','.TRUE.')

  IF(.NOT.calcEigVals) THEN
    SWRITE(UNIT_stdOut,'(A)') ' INIT EIGENVALUES DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')
    RETURN
  END IF

  stableDtfromEV=GETINT('stableDtfromEV','-1')
  FD_degree = GETINT('FD_degree','4')
  ! Compute total number of degrees of freedom
  !localDOF  = PP_nVar*(PP_N+1)**3*nElems
  !GlobalDOF = PP_nVar*(PP_N+1)**3*nGlobalElems
  EV_reducedDim = GETINT('EV_reducedDim','3')
  SELECT CASE(EV_reducedDim)
  CASE(1)
    Nx=PP_N; Ny=0;  Nz=0
  CASE(2)
    Nx=PP_N; Ny=PP_N;  Nz=0
  CASE DEFAULT
    Nx=PP_N; Ny=PP_N;  Nz=PP_N
  END SELECT
  localDOF  = PP_nVar*(Nx+1)*(Ny+1)*(Nz+1)*nElems
  GlobalDOF = PP_nVar*(Nx+1)*(Ny+1)*(Nz+1)*nGlobalElems
  
  EV_maxIter=GETINT('EV_maxIter','-1')
  IF(EV_maxIter.LE.0) THEN
    EV_maxIter=GlobalDOF
  ELSE
    EV_maxIter=MIN(EV_maxIter,GlobalDOF)
  END IF
  !only allocate hesseberg matrix on root!
  IF(MPIroot) ALLOCATE(Hmat(EV_maxIter,EV_maxIter)) 

  ! Allocate additional global fields required for the matrix-free matrix-vector product routine
  ALLOCATE(U0( 1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems))
  
  InitEigenValuesDone=.TRUE.
  SWRITE(UNIT_stdOut,'(A)') ' INIT EIGENVALUES DONE!'
  SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEigenValues


!===================================================================================================================================
!> Computes the Eigenvalues of an upper Hessenberg Matrix
!>
!===================================================================================================================================
SUBROUTINE EigenValues(t,dt,dt_scale)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_EigenValues_Vars
USE MOD_DG_vars           , ONLY: U
USE MOD_Output_Vars       , ONLY: ProjectName
USE MOD_Mesh_Vars         , ONLY: nElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)     :: t
REAL,INTENT(IN)      :: dt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)    :: dt_scale
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                 :: maxReal,AbsImag_maxReal,minReal,AbsImag_minReal,minRealonaxis,AbsImag_minRealonaxis
REAL                 :: maxAbsImag,Real_maxAbsImag
REAL                 :: tS, tE
REAL                 :: z
REAL                 :: G0,G,dtmax,dtnew
INTEGER              :: subDim
INTEGER              :: lWork
INTEGER              :: info
REAL,ALLOCATABLE     :: RealPart(:)   !< all  eigenvalues (real part)
REAL,ALLOCATABLE     :: ImagPart(:)   !< all  eigenvalues (imaginary part)
REAL, ALLOCATABLE    :: Work(:)
REAL                 ::U_tmp( 1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)

INTEGER              :: i 
LOGICAL              :: fullMat
!===================================================================================================================================
!default:
dtnew=dt
dt_scale=1.
IF(.NOT.calcEigVals) RETURN
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,*) 'Compute Eigenvalues:' 
U0=U
tS = FLUXOTIME()
fullMat=(EV_maxIter.EQ.GlobalDOF)
IF(fullMat)THEN
SWRITE(UNIT_stdOut,'(A,I8)') 'Computing full Linearized Matrix ...  globalDOF: ' , GlobalDOF
  subDim=GlobalDOF
  CALL LinearizedMatrix(t) !--> computes Hmat (for MPIroot)
ELSE
SWRITE(UNIT_stdOut,'(A,I8,A,I8)') 'Computing Upper Hessenberg Matrix ...  ',EV_maxIter," , globalDOF: " , GlobalDOF
  CALL arnoldi(t,SubDim) !--> computes Hmat (for MPIroot)
END IF


!save back
U=U0

tE = FLUXOTIME()
SWRITE(UNIT_stdOut,'(A,F10.5,A)') 'Matrix computed in ', tE - tS, ' s'
tS = FLUXOTIME()

IF(MPIroot)THEN
  WRITE(UNIT_stdOut,'(A)') 'Computing Eigenvalues ...'
  lWork = 11*SubDim
  ALLOCATE(RealPart(subDim))
  ALLOCATE(ImagPart(subDim))

  !call lapack routine
  ALLOCATE(work(lWork))
  IF(fullMat)THEN
    call dgeev('N','N', GlobalDOF, Hmat, GlobalDOF, RealPart, ImagPart, z,1,z,1,work,lwork,info)
  ELSE
    CALL dhseqr('E','N',SubDim,1,SubDim,Hmat(1:SubDim,1:SubDim),SubDim,RealPart(1:SubDim),ImagPart(1:SubDim),z,1,work,lwork,info)
  END IF
  IF(info.NE.0) &
    CALL abort(__STAMP__,&
       "Eigenvalue computation failed in Lapack routine dhseqr. Info=",info,999.)
  DEALLOCATE(Hmat,work)

END IF
tE = FLUXOTIME()
IF(MPIroot)THEN
  WRITE(UNIT_stdOut,'(A,F10.5,A)') 'Eigenvalues computed in ', tE - tS, ' s'
  WRITE(*,'(A,F10.5)') 'Eigenvalues of the spatial operator at t = ', t
  
  IF((MAXVAL(RealPart)*dt).GT. 1.0e-06)THEN
    G0=evalCharPoly(dt)
    WRITE(UNIT_stdOut,'(A,E25.15)') 'FOUND EIGENVALUES WITH POSITIVE REAL PART (Re*dt>1.0e-06)! max(|G|-1)=',G0-1.
    stableDTfromEV=-1 !switch off dt computation!
  END IF
  IF(stableDtfromEV.NE.-1)THEN
    !Find max stable timestep
    G0=evalCharPoly(dt)
    IF((G0-1.).GT.1.0E-8) THEN
      WRITE(UNIT_stdOut,'(A,E25.15)') 'timestep unstable!'
      dtmax=dt
      DO 
        G=evalCharPoly(dtmax)
        IF((G-1).LT.1E-08) EXIT
        dtmax=0.995*dtmax
        IF(dtmax.LT.1E-10) THEN
           WRITE(*,*) 'WARNING, dtmin<1E-10  !!!'
           EXIT
        END IF
      END DO
    ELSE
      !Find max stable timestep
      dtmax=dt
      DO
        G=evalCharPoly(dtmax)
        IF((G-1).GT.1E-08) EXIT
        dtmax=1.005*dtmax
      END DO
    END IF
 
    WRITE(UNIT_stdOut,'(A,E21.15)') '|G(dt)-1|       : ',G
    WRITE(UNIT_stdOut,'(A,F21.15)') 'dt_scale needed : ',dtmax/dt
    IF(stableDtfromEV.EQ.1)THEN !apply new timestep
      dt_scale=dtnew/dt 
      dtnew=0.99*dtmax  
      WRITE(UNIT_stdOut,'(A,E21.15)') 'new dt is set to : ',0.99*dtmax
    ELSE
      WRITE(UNIT_stdOut,'(A,E21.15)') 'new dt would be  : ',0.99*dtmax
    END IF
  END IF !stableDtfromEV.NE.-1
  
  OPEN(44,FILE=TRIM(ProjectName)//'_eigenvalues.dat',Status='UNKNOWN')
    WRITE(44,'(A)')'TITLE=""'
    WRITE(44,'(A)')'VARIABLES='
    WRITE(44,'(A)')'"Re","Im","Re_dt","Im_dt"'
    WRITE(44,'(A,I8)')'ZONE T="'//TRIM(ProjectName)//'",I=',SubDim
    maxReal=-1.E12
    minReal= 1.E12
    minRealonaxis= 1.E12
    maxAbsImag=-1.E12
    DO i = 1, SubDim
      WRITE(44,'(E25.15,3((","),E25.15))') RealPart(i),ImagPart(i),dtnew*RealPart(i),dtnew*ImagPart(i)
      !analyze
      IF(maxReal.LT.RealPart(i))THEN
        maxReal=RealPart(i)
        AbsImag_maxReal=ABS(ImagPart(i))
      END IF
      IF(minReal.GT.RealPart(i))THEN
        minReal=RealPart(i)
        AbsImag_minReal=ABS(ImagPart(i))
      END IF
      IF(maxAbsImag.LT.ABS(ImagPart(i)))THEN
        maxAbsImag=ABS(ImagPart(i))
        Real_maxAbsImag=RealPart(i)
      END IF
      IF(ABS(ImagPart(i)).LT.1.E-5)THEN
        IF(minRealonaxis.GT.RealPart(i))THEN
           minRealonaxis=RealPart(i)
        END IF
      END IF
    END DO
  CLOSE(44)
  DEALLOCATE(RealPart)
  DEALLOCATE(ImagPart)
  
  WRITE(UNIT_stdOut,'(A,E25.15,A,E25.15)') 'min Real*dt  : ',dtnew*minReal   ,' at |Imag|*dt ',dtnew*AbsImag_minReal
  WRITE(UNIT_stdOut,'(A,E25.15,A,E25.15)') 'max Real*dt  : ',dtnew*maxReal   ,' at |Imag|*dt ',dtnew*AbsImag_maxReal
  WRITE(UNIT_stdOut,'(A,E25.15,A,E25.15)') 'max |Imag|*dt: ',dtnew*maxAbsImag,' at  Real*dt  ',dtnew*Real_maxAbsImag
  WRITE(UNIT_stdOut,'(A,E25.15,A)') 'min Real*dt  : ',dtnew*minRealonaxis ,' on axis   '
  
END IF !MPIRoot

  CONTAINS

  !=================================================================================================================================
  !> Evaluate characteristic polynomial of chosen time integrator 
  !>
  !=================================================================================================================================
  FUNCTION evalCharPoly(dt_in)
  ! MODULES
  USE MOD_TimeDisc_Vars, ONLY: TimeDiscMethod
  ! uses subDim, RealPart, ImagPart from calling subroutine
  ! IMPLICIT VARIABLE HANDLING
  IMPLICIT NONE
  !---------------------------------------------------------------------------------------------------------------------------------
  ! INPUT VARIABLES
  REAL,INTENT(IN) :: dt_in
  !---------------------------------------------------------------------------------------------------------------------------------
  ! OUTPUT VARIABLES
  REAL            :: evalCharPoly
  !---------------------------------------------------------------------------------------------------------------------------------
  ! LOCAL VARIABLES 
  COMPLEX,DIMENSION(subDim) :: dtlambda
  !=================================================================================================================================
    dtlambda=dt_in*CMPLX(RealPart,ImagPart)
    SELECT CASE (TRIM(TimeDiscMethod))
    CASE('standardrk3-3') 
      dtlambda=1. + dtlambda+(0.5)*dtlambda**2+(1/6)*dtlambda**3
    CASE('carpenterrk4-5')
      dtlambda=1. + 1.       * dtlambda       &
                  + 0.5      * dtlambda ** 2  &
                  + (1./6.)  * dtlambda ** 3  &
                  + (1./24.) * dtlambda ** 4  &
                  + 0.005    * dtlambda ** 5
    CASE DEFAULT
      SWRITE(UNIT_StdOut,'(A,A)') 'characteristic polynomial not implemented for TimeDiscMethod: ' ,TRIM(TimeDiscMethod)
    END SELECT
    evalCharPoly=MAXVAL(ABS(dtlambda))
  END FUNCTION evalCharPoly

END SUBROUTINE EigenValues


!===================================================================================================================================
!> Arnoldi Algoriothm for building an upper hessenberg matrix
!>
!===================================================================================================================================
SUBROUTINE Arnoldi(t,subDim)
! MODULES
USE MOD_Globals
USE MOD_EigenValues_Vars,ONLY: localDOF,globalDOF,Hmat,EV_maxIter  !only allocated for MPIroot
USE MOD_LinAlg, ONLY: Norm_2,DotProd,MatVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)     :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT) :: subDim
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                 :: W(localDOF)
REAL                 :: V(localDOF,EV_maxIter)

INTEGER              :: i,j 
INTEGER              :: norm_w_eps

REAL                 :: EPS = 1.E-08
REAL                 :: wdotv,norm_w 
!===================================================================================================================================
v(:,1) = 0.
v(1,1) = 1.
DO j = 1, EV_maxIter
! Construct new Krylov Basis Vector
  CALL Matvec(v(:,j),w,t)
  !modified Gram-Schmidt
  DO i = 1, j
    wdotv = DotProd(w,v(:,i))
    w  = w - wdotv*v(:,i)
    IF(MPIroot) Hmat(i,j)=wdotv
  END DO ! i

  norm_w = Norm_2(w)

  IF(MPIroot) norm_w_eps =MERGE(1,0,norm_w>EPS)
#if MPI
  CALL MPI_BCAST(norm_w_eps,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
#endif /*MPI*/

  subDim=j
  IF ((norm_w_eps .EQ. 1).AND.(j<EV_maxIter)) THEN
    IF(MPIroot) Hmat(j+1,j) = norm_w
    v(:,j+1) = w / norm_w
  ELSE
    IF(MPIroot)THEN
      WRITE(UNIT_stdOut,'(A,I10)') 'Arnoldi has finished. Iteration number: ', j
      IF(j<globalDOF) WRITE(UNIT_stdOut,'(A,I10)') 'WARNING: STOPPED BEFORE ALL DOFS WERE CONSIDERED!!'   
      WRITE(UNIT_stdOut,'(A,E21.15)') '      Last entry H(j,j-1)= ', Hmat(j,j-1)
      WRITE(UNIT_stdOut,'(A,E21.15)') '      Last Norm   ||w||_2= ', norm_w
    END IF
    RETURN
  END IF
  IF (MOD(j,200) == 0) THEN
    SWRITE(UNIT_stdOut,'(A,I10,A,I10)') 'Iteration', j, ' out of ', EV_maxIter
  END IF
END DO ! j
!===================================================================================================================================
END SUBROUTINE Arnoldi

!===================================================================================================================================
!> Compute linearized operator matrix by using matrix-vector product with unit vectors
!>
!===================================================================================================================================
SUBROUTINE linearizedMatrix(t)
! MODULES
USE MOD_Globals
USE MOD_EigenValues_Vars,ONLY: localDOF,globalDOF,Hmat,EV_maxIter  !only allocated for MPIroot
USE MOD_LinAlg, ONLY: Norm_2,DotProd,MatVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)     :: t
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                 :: W(localDOF)
REAL                 :: V(localDOF)

INTEGER              :: i,j 
INTEGER              :: norm_w_eps

REAL                 :: EPS = 1.E-08
REAL                 :: wdotv,norm_w 
INTEGER              :: procDOF_MPI(1:nProcessors)
INTEGER              :: offsetDOF(0:nProcessors)
!===================================================================================================================================
#if MPI
CALL MPI_GATHER(localDOF,1,MPI_INTEGER,procDOF_MPI,1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
IF(MPIroot)THEN
  offsetDOF(0)=0
  DO i=1,nProcessors
    offsetDOF(i)=offsetDOF(i-1)+procDOF_MPI(i)
  END DO
END IF !MPIroot
CALL MPI_BCAST(offsetDOF,nProcessors+1,MPI_INTEGER,0,MPI_COMM_WORLD,iError)
!IF(offsetDOF(nProcessors).NE.GlobalDOF) &
!  CALL abort(__STAMP__,&
!             "problem with offsetDOF in linearizedMatrix")
#else
offsetDOF(0)=0
offsetDOF(1)=GlobalDOF
#endif /*MPI*/
DO j = 1, EV_maxIter
  v=0.0d0
  IF((j.GT.offsetDOF(myRank)).AND.(j.LE.offsetDOF(myRank+1))) THEN
    v(j-offsetDOF(myRank))=1.0d0
  END IF
  CALL Matvec(v(:),w(:),t)
#if MPI
  CALL MPI_GATHERV(w,localDOF,MPI_DOUBLE_PRECISION,Hmat(:,j),procDOF_MPI,offsetDOF(0:nProcessors-1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iError)
!  IF(MPIroot)THEN
!    IF(SUM(ABS(Hmat(1:localDOF,j)-w)).GT.1.0e-12*REAL(localDOF)) &
!      CALL abort(__STAMP__, &
!                 "problem with offset in Gatherv in linearizedMatrix")
!  END IF
#else
  Hmat(:,j)=w
#endif /*MPI*/
  IF (MOD(j,200) == 0) THEN
    SWRITE(UNIT_stdOut,'(A,I10,A,I10)') 'Iteration', j, ' out of ', EV_maxIter
  END IF
END DO ! j
!===================================================================================================================================
END SUBROUTINE linearizedMatrix

!===================================================================================================================================
!> Deallocates the additional global fields required by the matrix-free method
!>
!===================================================================================================================================
SUBROUTINE FinalizeEigenValues()
! MODULES
USE MOD_EigenValues_Vars 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
SDEALLOCATE(U0)
SDEALLOCATE(Hmat)
END SUBROUTINE FinalizeEigenValues

END MODULE MOD_EigenValues

