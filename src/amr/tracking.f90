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
    INTEGER :: Count = 0
CONTAINS

! FUNCTION GetShockCapturing(Uin) result(eta_dof)
!     USE MOD_PreProc
!     USE MOD_ChangeBasis,            ONLY : ChangeBasis3D
!     IMPLICIT NONE
!     REAL, DIMENSION(PP_nVar, 0:PP_N, 0:PP_N, 0:PP_N), INTENT(IN) :: Uin
!     REAL, DIMENSION(1:1, 0:PP_N, 0:PP_N, 0:PP_N) :: Umod
!     REAL, DIMENSION(0:PP_N, 0:PP_N) :: Vdm_Leg, sVdm_Leg
!     REAL :: LU, LUM1, LUM2, LU_N, LU_NM1!, eta_dof, eta_min, eta_max, eps0, RhoInf, Pinf, RhoMax, RhoMin, Xmin(3), Xmax(3), Abst
!     REAL eta_dof
!     CALL ChangeBasis3D(1, PP_N, PP_N, sVdm_Leg, Uin(1:1,:,:,:), Umod)
!     LU = SUM(Umod(1, :, :, :)**2)
!     LUM1 = SUM(Umod(1, 0:PP_N - 1, 0:PP_N - 1, 0:PP_N - 1)**2)
!     LUM2 = SUM(Umod(1, 0:PP_N - 2, 0:PP_N - 2, 0:PP_N - 2)**2)
!     LU_N = LU - LUM1
!     LU_NM1 = LUM1 - LUM2
!     ! DOF energy indicator
!     eta_dof = LOG10(MAX(LU_N / LU, LU_NM1 / LUM1, TINY(1.0)))
    
! END FUNCTION 

!==================================================================================================================================
!> Specifies the initial AMR refinement
!==================================================================================================================================
subroutine InitialAMRRefinement()
  USE MOD_PreProc     , only: PP_N
  use MOD_AMR_Vars    , only: InitialRefinement, UseAMR, MaxLevel, MinLevel
  use MOD_AMR         , only: RunAMR
  use MOD_Mesh_Vars   , only: nElems, Elem_xGP
  implicit none
  !-local-variables-----------------------------------------
  real    :: r
  logical :: RefineElem
  integer :: iter
  integer :: iElem, i,j,k
  integer, allocatable :: ElemToRefineAndCoarse(:)
  !---------------------------------------------------------
  
  if (.not. UseAMR) return
  
  select case (InitialRefinement)  
    case default ! Use the default indicator up to max-level
      do iter = 1,MaxLevel
        call ShockCapturingAMR()
        call InitData()
      end do
    
    case(1) ! Refine any element containing a node in the spherewith radius r=0.2 to the MaxLevel
      do iter = 1,MaxLevel
        allocate (ElemToRefineAndCoarse(1:nElems))!
        
        ! Fill ElemToRefineAndCoarse
        ! --------------------------
        do iElem=1, nElems  
          ! Check if the element has a node on the desired region
          RefineElem = .FALSE.
          do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
            r = sqrt(sum(Elem_xGP(:,i,j,k,iElem)**2))
            if (r <= 0.2) then
              RefineElem = .TRUE.
              exit ; exit ; exit
            end if
          end do       ; end do       ; end do
          
          ! Set that element for refinement
          if (RefineElem) then
            ElemToRefineAndCoarse(iElem) = MaxLevel
          else
            ElemToRefineAndCoarse(iElem) = MinLevel
          end if
          
        end do
        
        ! Refine
        ! ------
        CALL RunAMR(ElemToRefineAndCoarse)
        
        deallocate (ElemToRefineAndCoarse)
      end do
      call InitData()
  end select
  

end subroutine


    SUBROUTINE ShockCapturingAMR()
        !   USE MOD_AMR_vars,            ONLY: P4EST_PTR, CONNECTIVITY_PTR
        USE MOD_PreProc
        USE MOD_Globals,                ONLY : MPIroot
        USE MOD_DG_Vars,                ONLY : U
        USE MOD_AMR,                    ONLY : RunAMR, LoadBalancingAMR, SaveMesh;
        USE MOD_Mesh_Vars,              ONLY : nElems
        USE MOD_Basis,                  ONLY : BuildLegendreVdm
        USE MOD_Indicators,             ONLY : ShockSensor_PerssonPeraire
        USE MOD_AMR_Vars,               ONLY : MinLevel, MaxLevel, RefineVal, CoarseVal
        USE MOD_P4EST,                  ONLY: SaveP4est
        ! USE MOD_Equation_Vars,      ONLY: kappaM1, RefStatePrim, IniRefState
        IMPLICIT NONE
        ! SAVE
        !Local variables
        INTEGER, ALLOCATABLE, TARGET :: ElemToRefineAndCoarse(:) ! positive Number - refine, negative - coarse, 0 - do nothing
        INTEGER :: iElem
         
        ! REAL, DIMENSION(0:PP_N, 0:PP_N) :: Vdm_Leg, sVdm_Leg
        ! REAL :: LU, LUM1, LUM2, LU_N, LU_NM1, eta_dof, eta_min, eta_max, eps0, RhoInf, Pinf, RhoMax, RhoMin, Xmin(3), Xmax(3), Abst
        ! INTEGER :: iXMax(3), iXMin(3),i,j,k
        REAL    :: eta_dof
            
        ALLOCATE(ElemToRefineAndCoarse(1:nElems))!
        ElemToRefineAndCoarse = MinLevel;
        
    !! < ----- Commented for the production ----- >
 
         DO iElem = 1, nElems
             eta_dof = LOG10( ShockSensor_PerssonPeraire(U(:,:,:,:,iElem)))
             IF (eta_dof .GE. RefineVal) THEN
                 ElemToRefineAndCoarse(iElem) = MaxLevel
           
             ELSE IF (eta_dof .LE. CoarseVal) THEN
                 ElemToRefineAndCoarse(iElem) = -MinLevel - 1
             ELSE
                 ElemToRefineAndCoarse(iElem) = MinLevel
             END IF
        ENDDO
        !! < ----- Commented for the production ----- >

        ! IF (Count .EQ. 0 ) THEN
            ! COUNT = 1; 
        ! DO iElem = 1,nElems

        !     ! ==== > StressTest
        !     ! CALL RANDOM_NUMBER(R)
        !     ! ! PRINT *, "R = ", R
        !     ! If (R .LE. 0.2) THEN 
        !     !     ElemToRefineAndCoarse(iElem) = MaxLevel  
        !     !     ! PRINT *, "Refine = ", R       
        !     ! ELSE
        !     !     ElemToRefineAndCoarse(iElem) = -MinLevel - 1          
        !     !     ! PRINT *, "Coarse = ", R       
        !     ! ENDIF
        !     ! ==== < StressTest
        !     ! IF ((Elem_xGP(1,0,0,0,iElem)) .LE. 0.1499 .OR. &
        !     !     (Elem_xGP(2,0,0,0,iElem)) .LE. 0.1499 .OR. &
        !     !     ! (Elem_xGP(3,0,0,0,l)) .LE. 0.1499 .OR. &
        !     !     (Elem_xGP(1,PP_N,PP_N,PP_N,iElem)) .GE. 1.-0.1499 .OR. &
        !     !     (Elem_xGP(2,PP_N,PP_N,PP_N,iElem)) .GE. 1.-0.1499) THEN ! .OR. &
        !     !     ! (Elem_xGP(3,PP_N,PP_N,PP_N,l)) .GE. 1.-0.1499 &
                
        !     !         ElemToRefineAndCoarse(iElem) = MaxLevel          
        !     !     ELSE 
        !     !         ElemToRefineAndCoarse(iElem) = MinLevel          
        !     !     ! PRINT *, "REFINE!!!!!!!!!!!!!!!!"
        !     !     ! PRINT *, "=>>>>", Minval(Elem_xGP(1,:,:,:,iElem))
        !     !     ! CALL EXIT()
        !     ! ENDIF
        ! ENDDO
          ! IF (MPIRoot) THEN
            ! ElemToRefineAndCoarse(1) = MaxLevel
            ! ENDIF
        ! ENDIF
        ! CALL EXIT()
        ! ElemToRefineAndCoarse(1) = 1 

        CALL RunAMR(ElemToRefineAndCoarse);
        ! IF ((Count .EQ. 1) .OR. (Count .EQ. 0)) THEN
        !     COUNT = 1; 
        !    CALL InitData();
        !     IF (MPIRoot) PRINT *, "InitData"
        ! ENDIF
        Deallocate(ElemToRefineAndCoarse)
    !     PRINT *, "REFINED"
    !     CALL SaveMesh("new_mesh.h5")
    !     Call SaveP4est("new.p4est")

    !    CALL EXIT()
       

    END SUBROUTINE ShockCapturingAMR

    SUBROUTINE InitData()
        USE MOD_PreProc
        USE MOD_DG_Vars, ONLY : U

        USE MOD_Equation_Vars, ONLY :IniExactFunc
        USE MOD_Equation, ONLY : FillIni
        USE MOD_Restart_Vars, only: DoRestart
        USE MOD_AMR, ONLY : RunAMR
        IMPLICIT NONE


        
        if (.not. DoRestart) call FillIni(IniExactFunc,U)
        
!#        PP = size(U(1, :, 0, 0, 1)) - 1;
!#        nVar = size(U(:, 0, 0, 0, 1));
!#        X0 = 0.5
!#        Y0 = 0.5
!#        Mejecta = 0.5 !0.5
!#        SigmaEjecta = 15.e-2!3.e-2
!#        Eblast = 0*0.1 !1.
!#        Sigmablast = 5.e-2!2.e-2
!#        ! R = 3. / 2. !sqrt(0.5) !0.005
!#        ! sigma = 1.
!#        ! MachInf = 0.4 !sqrt(skappa) !0.5
!#        ! beta = Machinf * 27. / 3.14159 / 4. * exp(2. / 9.)!Machinf*5./3.14159/2.*exp(1.)!1./5.
!#        ! alfa = 3.14159 / 2.
!#!        PRINT *, "IniExactFunc", IniExactFunc
!#        ! DO l = 1, nElems
!#        !     IF ((Elem_xGP(1,0,0,0,l)) .LE. 0.1499 .OR. &
!#        !         (Elem_xGP(2,0,0,0,l)) .LE. 0.1499 .OR. &
!#        !         (Elem_xGP(3,0,0,0,l)) .LE. 0.1499 .OR. &
!#        !         (Elem_xGP(1,PP_N,PP_N,PP_N,l)) .GE. 1.-0.1499 .OR. &
!#        !         (Elem_xGP(2,PP_N,PP_N,PP_N,l)) .GE. 1.-0.1499 .OR. &
!#        !         (Elem_xGP(3,PP_N,PP_N,PP_N,l)) .GE. 1.-0.1499 ) THEN !
!#        !             U(1,:,:,:,l ) = 5.
!#        !         ELSE
!#        !             U(1,:,:,:,l ) = 1.
!#        !             ! PRINT *," ====>", Maxval(Elem_xGP(1,:,:,:,l))
!#        !     ENDIF
!#        ! ENDDO
!#        ! DO  Iter = 1, 5
!#            DO iElem = 1, nElems
!#                DO i = 0, PP; DO j = 0, PP; DO k = 0, PP;

!#                    X = Elem_xGP(1,i,j,k,iElem)
!#                    Y = Elem_xGP(2,i,j,k,iElem)
!#                    X = (x-x0)
!#                    Y = (y-y0)

!#                    Ux = 1.
!#                    Vy = 1.
!#                    Rho = 2. + Mejecta / 2. / PP_Pi / SigmaEjecta / SigmaEjecta * &
!#                        exp(-0.5 * (x * x + y * y) / SigmaEjecta / SigmaEjecta)
!#                    P = 1.e+2 * skappam1 + Eblast/(2 * PP_Pi * SigmaBlast * SigmaBlast) * &
!#                        exp(-0.5 * (x * x + y * y) / Sigmablast / Sigmablast)
!#                    	    ! F = -1./2./sigma/Sigma*(x*x/R/R + y*y/R/R)
!#                    	    ! Omega = beta * exp(f)
!#                    	    ! Rho = (1. - ((kappam1)/2.)*Omega*Omega)**skappam1
!#                    	    ! Ux = Machinf*Cos(alfa) - y/r*Omega
!#                    	    ! Vy = MachInf*sin(alfa) + x/r*Omega
!#                    	    ! P = skappa*(1 - kappam1/2. * Omega*Omega)**(kappa*skappam1)
!#                      Prim = (/Rho, Ux, Vy, 0., P/)
!#                    ! Prim = (/1.29, 0.,0., 0., 100000./)
!#                    ! CALL PrimToCons(Prim, U(:,i,j,k,iElem))
!#                    ! SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu)
!#                    CALL ExactFunc(IniExactFunc, 0., Elem_xGP(:, i, j, k, iElem), U(:, i, j, k, iElem));
!#                ENDDO;
!#                ENDDO;
!#                ENDDO;
!#                ! i,j,k
!#            ENDDO ! nElems
!#        !     ! Call RunAMR();
!#        ! ENDDO !Iter

    END SUBROUTINE InitData

END MODULE MOD_AMR_tracking
#endif
