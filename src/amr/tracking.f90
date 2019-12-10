#include "defines.h"
#include "amr_f.h"

#if USE_AMR
!==================================================================================================================================
!> Contains routines for AMR as trecking function
!==================================================================================================================================
MODULE MOD_AMR_tracking

INTERFACE ShockCapturingAMR
  MODULE PROCEDURE ShockCapturingAMR
END INTERFACE


! IMPLICIT NONE
  INTEGER :: doLBalance = 0
CONTAINS

SUBROUTINE ShockCapturingAMR()
!   USE MOD_AMR_vars,            ONLY: P4EST_PTR, CONNECTIVITY_PTR
    USE MOD_PreProc
  USE MOD_Globals,              ONLY: MPIroot, myrank
  USE MOD_DG_Vars,            ONLY: U
  USE MOD_AMR,                 ONLY: RunAMR,LoadBalancingAMR, SaveMesh,LoadBalancingAMRold;
  USE MOD_Mesh_Vars,           ONLY:nElems, Elem_xGP, nGlobalElems
  USE MOD_Interpolation_Vars, ONLY: xGP
  USE MOD_Basis,         ONLY: BuildLegendreVdm
  USE MOD_ChangeBasis,        ONLY: ChangeBasis3D
  ! USE MOD_Equation_Vars,      ONLY: kappaM1, RefStatePrim, IniRefState
  IMPLICIT NONE
  ! SAVE
    !Local variables
  INTEGER, ALLOCATABLE, TARGET  :: ElemToRefineAndCoarse(:) ! positive Number - refine, negative - coarse, 0 - do nothing
  INTEGER :: l, MaxLevel, MinLevel
  ! REAL :: 
  REAL,DIMENSION(1:1,0:PP_N,0:PP_N,0:PP_N) :: Uind,Umod
  REAL,DIMENSION(0:PP_N,0:PP_N)              :: Vdm_Leg, sVdm_Leg
  REAL                                     :: LU,LUM1,LUM2,LU_N,LU_NM1,eta_dof,eta_min,eta_max,eps0, RhoInf, Pinf, RhoMax, RhoMin, Xmin(3), Xmax(3), Abst
  INTEGER               :: iXMax(3), iXMin(3)
  LOGICAL :: doBalance = .TRUE.

MaxLevel = 4;
MinLevel = 2;

!MaxLevel = 3;

!MinLevel = 0;
ALLOCATE(ElemToRefineAndCoarse(1:nElems))! 
    ElemToRefineAndCoarse=0;
! ALLOCATE(sVdm_Leg(0:PP_N,0:PP_N))
! ALLOCATE(Vdm_Leg(0:PP_N,0:PP_N))
CALL BuildLegendreVdm(PP_N,xGP,Vdm_Leg,sVdm_Leg)

!print *, RefStatePrim(IniRefState,:)
! RhoInf = RefStatePrim(IniRefState,1)
!CALL EXIT()
DO l=1,nElems
!     ! if (l .EQ. 1) U()
     Uind(1,:,:,:) = U(1,:,:,:,l) 
! !    Uind(1,:,:,:) = KappaM1*(U(5,:,:,:,l)-0.5*(SUM(U(2:4,:,:,:,l)*U(2:4,:,:,:,l))/U(1,:,:,:,l)))
! !    Uind(1,:,:,:) = Uind(1,:,:,:)*U(1,:,:,:,l)
     CALL ChangeBasis3D(1,PP_N,PP_N,sVdm_Leg,Uind,Umod)
!       ! Compute (truncated) error norms
   LU     = SUM(Umod(1,:,:,:)**2)
   LUM1   = SUM(Umod(1,0:PP_N-1,0:PP_N-1,0:PP_N-1)**2)
   LUM2   = SUM(Umod(1,0:PP_N-2,0:PP_N-2,0:PP_N-2)**2)
   LU_N   = LU-LUM1
   LU_NM1 = LUM1-LUM2


  ! DOF energy indicator
  eta_dof = LOG10(MAX(LU_N/LU,LU_NM1/LUM1, TINY(1.0)))
!   print *, " Umod(1,:,1,1) = ", Umod(1,:,1,1)
!   print *, " Umod(1,1,:,1) = ", Umod(1,1,:,1)
!   print *, " Umod(1,1,1,:) = ", Umod(1,1,1,:)
! print *, " LU = ", LU
! print *, " LUM1 = ", LUM1
! print *, " LUM2 = ", LUM2
! print *, " LU_N = ", LU_N
! print *, " LU_NM1 = ", LU_NM1
! print *, " LU = ", LU
!eta_dof = MAXVAL(ABS(U(1,:,:,:,l)- RhoInf)/RhoInf)
!eta_dof= LOG10(eta_dof)
!iXmax = 0
!iXmin = 0
!RhoMax = lOG10(MAXVAL(U(1,:,:,:,l)))
!RhoMin = LOG10(MINVAL(U(1,:,:,:,l)))
! iXMAX = MAXLOC(U(1,:,:,:,l))
! iXMIN = MINLOC(U(1,:,:,:,l))

 !xMin = Elem_xGP(:,iXmin(1), iXmin(2),iXmin(3),l)
! xMax = Elem_xGP(:,iXmax(1), iXmax(2),iXmax(3),l)

 !Abst = SQRT(SUM((xMax(:) - xMin(:))**2))
 !print *, eta_dof
 !IF ((ixMax(1) .EQ. iXMin(1)) .AND. (ixMax(2) .EQ. iXMin(2)) .AND. (ixMax(3) .EQ. iXMin(3)) )  THEN
 !   eta_dof = 0.
 !ELSE 
  !  eta_dof = ABS( (RhoMax - RhoMin))/Abst
    !print *, "eta_dof = ", eta_dof, ABST
 !ENDIF

  ! Artificial Viscosity
  ! eta_min = -9.0
  !eta_max = -3.0
  ! for N-4
 ! eta_min = -13.5
  !eta_max = -12.7
!  eta_min = -14.5
!  eta_max = -11.7
  !for N=5
  eta_min = -15.5
  eta_max = -10.0
  !eta_min = -8.
  !eta_max = -6.99
  
  !eta_min = 0.0001/250.
  !eta_max = 0.1/100.
  ! eps0 = 0.01
  IF (eta_dof.GE.eta_max) THEN
    ! nu(l) = eps0
    ElemToRefineAndCoarse(l) = MaxLevel
    !Print *, "Refine!!!", eta_dof
    ELSE IF (eta_dof.LE.eta_min) THEN
    ! nu(l) = 0.
    ElemToRefineAndCoarse(l) = -MinLevel -1 
   ! print *, "Coarse!", eta_dof
  ELSE
    ! nu(l) = 0.5*eps0*(1.0+SIN(PP_Pi*(eta_dof-0.5*(eta_max+eta_min))/(eta_max-eta_min)))
    ElemToRefineAndCoarse(l) = 0
  END IF
 
 
ENDDO
!CALL EXIT()
!DO l=1,nElems
 !   IF (ElemToRefineAndCoarse(l)>1) THEN
!
 !   Print *, ElemToRefineAndCoarse
 !   CALL EXIT()
 !   ENDIF
!ENDDO
! PRINT *, "RunAMRs"
!ElemToRefineAndCoarse = 0
!IF (MPIRoot) THEN
!  ElemToRefineAndCoarse(1) = 2
!ENDIF
CALL RunAMR(ElemToRefineAndCoarse);
! CALL RunAMR(ElemToRefineAndCoarse);
! IF (MPIRoot) THEN
!   Print *, "AMR RUN!"
! ENDIF
Deallocate(ElemToRefineAndCoarse)
doLBalance = doLBalance+1
IF (doLBalance .EQ. 1) THEN
    doLBalance = 0;
    ! IF (MPIRoot)  print *, "Balance Loading START: Number of Global Elements = ", nGlobalElems
    ! PRINT *, "BEFORE: MPIRANK = ", myrank, "nELems = ", nElems
    IF (doBalance) CALL LoadBalancingAMR()
    IF (MPIRoot) THEN
        print *, "LoadBalance: Done! nGlobalElems =", nGlobalElems
    ENDIF
ENDIF

END SUBROUTINE ShockCapturingAMR

SUBROUTINE InitVortex()
    USE MOD_DG_Vars,            ONLY: U
    USE MOD_Mesh_Vars,          ONLY: Elem_xGP, nElems
    USE MOD_Equation_Vars,      ONLY: kappam1, skappam1, PrimToCons, kappa, IniExactFunc
    USE MOD_Equation,           ONLY: ExactFunc
    USE MOD_AMR,                 ONLY: RunAMR
    IMPLICIT NONE
    INTEGER :: iElem, i,j,k,PP, nVar, Iter
    REAL :: X0, Y0, X, Y, Rho, Ux,Vy, F, OMEGA,R, S, MachInf, Beta,alfa, PInf, Tinf, sigma, P, Prim(5)

    PP = size(U(1,:,0,0,1))-1;
    nVar = size(U(:,0,0,0,1));
    X0=5.!0.5
    Y0=5.!0.5
    R=3./2. !sqrt(0.5) !0.005
    sigma = 1.
    MachInf = 0.4 !sqrt(skappa) !0.5
    beta = Machinf*27./3.14159/4.*exp(2./9.)!Machinf*5./3.14159/2.*exp(1.)!1./5.
    alfa=3.14159/2.
    DO  Iter= 1, 5
    DO iElem = 1, nElems
	DO i=0,PP; DO j=0,PP; DO k =0,PP;
          
    
    ! X = Elem_xGP(1,i,j,k,iElem)
    ! 	    Y = Elem_xGP(2,i,j,k,iElem)
    ! 	    X = (x-x0)
    ! 	    Y = (y-y0)
    ! 	    F = -1./2./sigma/Sigma*(x*x/R/R + y*y/R/R)
    ! 	    Omega = beta * exp(f)
    ! 	    Rho = (1. - ((kappam1)/2.)*Omega*Omega)**skappam1
    ! 	    Ux = Machinf*Cos(alfa) - y/r*Omega
    ! 	    Vy = MachInf*sin(alfa) + x/r*Omega
    ! 	    P = skappa*(1 - kappam1/2. * Omega*Omega)**(kappa*skappam1)
    !         Prim = (/Rho, Ux, Vy, 0., P/)
            ! Prim = (/1.29, 0.,0., 0., 100000./)
        ! CALL PrimToCons(Prim, U(:,i,j,k,iElem))    	    
    !SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu) 
		CALL ExactFunc(IniExactFunc, 0., Elem_xGP(:,i,j,k,iElem), U(:,i,j,k,iElem));
        ENDDO; ENDDO; ENDDO;! i,j,k
	ENDDO ! nElems
    ! Call RunAMR();
    ENDDO !Iter

END SUBROUTINE

END MODULE MOD_AMR_tracking
#endif