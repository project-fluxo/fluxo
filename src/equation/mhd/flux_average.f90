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

!==================================================================================================================================
!> Routines to compute two-point average fluxes for the volint when using the split-form (DiscType=2)
!==================================================================================================================================
MODULE MOD_Flux_Average
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE EvalEulerFluxAverage3D_eqn
  MODULE PROCEDURE EvalEulerFluxAverage3D_eqn
END INTERFACE

!INTERFACE EvalEulerFluxTilde3D_eqn
!  MODULE PROCEDURE EvalEulerFluxTilde3D_eqn
!END INTERFACE

#if NONCONS
INTERFACE AddNonConsFluxTilde3D
  MODULE PROCEDURE AddNonConsFluxTilde3D
END INTERFACE

INTERFACE AddNonConsFluxVec
  MODULE PROCEDURE AddNonConsFluxVec
END INTERFACE
#endif /*NONCONS*/

INTERFACE EvalUaux
  MODULE PROCEDURE EvalUaux
END INTERFACE

INTERFACE StandardDGFlux
  MODULE PROCEDURE StandardDGFlux
END INTERFACE

INTERFACE StandardDGFluxVec
  MODULE PROCEDURE StandardDGFluxVec
END INTERFACE

INTERFACE StandardDGFluxDealiasedMetricVec
  MODULE PROCEDURE StandardDGFluxDealiasedMetricVec
END INTERFACE

INTERFACE EntropyandKinEnergyConservingFlux
  MODULE PROCEDURE EntropyandKinEnergyConservingFlux
END INTERFACE

INTERFACE EntropyandKinEnergyConservingFluxVec
  MODULE PROCEDURE EntropyandKinEnergyConservingFluxVec
END INTERFACE

INTERFACE LN_MEAN 
  MODULE PROCEDURE LN_MEAN
END INTERFACE


#if (PP_DiscType==2)
PUBLIC:: VolInt_splitForm_eqn
#endif /*PP_DiscType==2*/
PUBLIC::EvalEulerFluxAverage3D_eqn
PUBLIC::EvalEulerFluxTilde3D_eqn
#if NONCONS
PUBLIC::AddNonConsFluxTilde3D
PUBLIC::AddNonConsFluxVec
#endif /*NONCONS*/
PUBLIC::EvalUaux
PUBLIC::StandardDGFlux
PUBLIC::StandardDGFluxVec
PUBLIC::StandardDGFluxDealiasedMetricVec
PUBLIC::EntropyandKinEnergyConservingFlux
PUBLIC::EntropyandKinEnergyConservingFluxVec
PUBLIC::LN_MEAN

!==================================================================================================================================
! local definitions for inlining / optimizing routines,
! depending on PP_VolFlux 
#if PP_VolFlux==0
#  define PP_VolumeFluxAverageMat StandardDGFluxMat 
#elif PP_VolFlux==10
#  define PP_VolumeFluxAverageMat EntropyAndKinEnergyConservingFluxMat
#else
!default, using PP_VolumeFluxAverageVec inside
#  define PP_VolumeFluxAverageMat GeneralFluxMat 
#endif

#if PP_VolFlux==0
#  define PP_VolumeFluxAverageVec StandardDGFluxVec 
#elif PP_VolFlux==1
#  define PP_VolumeFluxAverageVec StandardDGFluxDealiasedMetricVec
#elif PP_VolFlux==10
#  define PP_VolumeFluxAverageVec EntropyandKinEnergyConservingFluxVec
#else
#  define PP_VolumeFluxAverageVec VolumeFluxAverageVec
#endif

!==================================================================================================================================

CONTAINS


#if (PP_DiscType==2)
!==================================================================================================================================
!> Computes the volume integral using flux differencing 
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut=0. and is updated with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolInt_SplitForm_eqn(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars   ,ONLY:DvolSurf_T,U
USE MOD_Mesh_Vars ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Mesh_Vars ,ONLY:nElems
!#if PARABOLIC
!USE MOD_Flux      ,ONLY:EvalDiffFluxTilde3D
!USE MOD_DG_Vars   ,ONLY:D_Hat_T
!#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)                                :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N):: rtilde
INTEGER                                           :: i,j,k,l,iElem
!==================================================================================================================================

DO iElem=1,nElems


!  !opt_v2, with larger calls 
  DO k=0,PP_N; DO j=0,PP_N
    !diagonal (consistent) part
    CALL PP_VolumeFluxAverageMat(              U(:,:,j,k,iElem), &
                                  Metrics_fTilde(:,:,j,k,iElem), &
                                          rtilde(:,:,:)         ) 
    
    DO i=0,PP_N; DO l=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf_T(l,i)*rtilde(:,l,i)
    END DO; END DO !i,l
  END DO; END DO ! j,k
  DO k=0,PP_N; DO i=0,PP_N
    CALL PP_VolumeFluxAverageMat(             U(:,i,:,k,iElem), &
                                 Metrics_gTilde(:,i,:,k,iElem), &
                                         rtilde(:,:,:)          )
    DO j=0,PP_N; DO l=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf_T(l,j)*rtilde(:,l,j)
    END DO; END DO !j,l
  END DO; END DO ! i,k
  DO j=0,PP_N; DO i=0,PP_N
    CALL PP_VolumeFluxAverageMat(             U(:,i,j,:,iElem), &
                                 Metrics_hTilde(:,i,j,:,iElem), &
                                         rtilde(:,:,:)          )
    DO k=0,PP_N; DO l=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf_T(l,k)*rtilde(:,l,k)
    END DO; END DO !j,l
  END DO; END DO ! i,j

!SAME SPEED
!  DO k=0,PP_N; DO j=0,PP_N
!    !diagonal (consistent) part
!    CALL TwoPointEntropyConservingFluxMat(              U(:,:,j,k,iElem), &
!!                                                     Uaux(:,:,j,k)      , &
!                                           Metrics_fTilde(:,:,j,k,iElem), &
!                                                   ftilde(:,:,:,j,k)         ) 
!    
!  END DO; END DO ! j,k
!  DO k=0,PP_N; DO i=0,PP_N
!    CALL TwoPointEntropyConservingFluxMat(             U(:,i,:,k,iElem), &
!!                                                    Uaux(:,i,:,k)      , &
!                                          Metrics_gTilde(:,i,:,k,iElem), &
!                                                  gtilde(:,:,:,i,k)          )
!  END DO; END DO ! i,k
!  DO j=0,PP_N; DO i=0,PP_N
!    CALL TwoPointEntropyConservingFluxMat(             U(:,i,j,:,iElem), &
!!                                                    Uaux(:,i,j,:)      , &
!                                          Metrics_hTilde(:,i,j,:,iElem), &
!                                                  htilde(:,:,:,i,j)          )
!  END DO; END DO ! i,j
!  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO l=0,PP_N
!     Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf_T(l,i)*ftilde(:,l,i,j,k) &
!                                           + Dvolsurf_T(l,j)*gtilde(:,l,j,i,k) &
!                                           + Dvolsurf_T(l,k)*htilde(:,l,k,i,j)
!  END DO; END DO; END DO; END DO !j,l

!#if PARABOLIC
!  CALL EvalDiffFluxTilde3D(iElem,ftildeDiff,gtildeDiff,htildeDiff)
!  !compute Diffusion flux contribution of 
!  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO l=0,PP_N
!      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) +    D_Hat_T(l,i)*ftildeDiff(:,l,j,k)  &
!                                            +    D_Hat_T(l,j)*gtildeDiff(:,i,l,k)  &
!                                            +    D_Hat_T(l,k)*htildeDiff(:,i,j,l)
!  END DO; END DO; END DO; END DO !j,l
!#endif /*PARABOLIC*/

END DO ! iElem
END SUBROUTINE VolInt_SplitForm_eqn
#endif /*PP_DiscType==2*/

!==================================================================================================================================
!> Compute flux differences in 3D, making use of the symmetry and appling also directly the metrics  
!==================================================================================================================================
SUBROUTINE EvalEulerFluxAverage3D_eqn(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars        ,ONLY:U
USE MOD_Mesh_Vars      ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
#if PP_VolFlux==-1
USE MOD_Equation_Vars  ,ONLY:VolumeFluxAverageVec !pointer to flux averaging routine
#endif
USE MOD_Equation_Vars  ,ONLY:nAuxVar
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem  !< current element ID from volint
REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde,gtilde,htilde !< 4D transformed fluxes (iVar,i,,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N):: ftilde_c,gtilde_c,htilde_c !central euler flux at ijk 
REAL,DIMENSION(nAuxVar,0:PP_N,0:PP_N,0:PP_N)  :: Uaux                       !auxiliary variables
INTEGER             :: i,j,k,l
!==================================================================================================================================

!opt_v1
CALL EvalEulerFluxTilde3D_eqn(              U(:,:,:,:,iElem) &
                              ,Metrics_fTilde(:,:,:,:,iElem) &
                              ,Metrics_gTilde(:,:,:,:,iElem) &
                              ,Metrics_hTilde(:,:,:,:,iElem) &
                              ,ftilde_c,gtilde_c,htilde_c, Uaux)
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part
  ftilde(:,i,i,j,k)=ftilde_c(:,i,j,k) 
  DO l=i+1,PP_N
    CALL PP_VolumeFluxAverageVec(             U(:,i,j,k,iElem),              U(:,l,j,k,iElem), &
                                           Uaux(:,i,j,k)      ,           Uaux(:,l,j,k)      , &
                                 Metrics_fTilde(:,i,j,k,iElem), Metrics_fTilde(:,l,j,k,iElem), &
                                       ftilde(:,l,i,j,k)                                       )
    ftilde(:,i,l,j,k)=ftilde(:,l,i,j,k) !symmetric
  END DO!l=i+1,N
END DO; END DO; END DO ! i,j,k
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part
  gtilde(:,j,i,j,k)=gtilde_c(:,i,j,k) 
  DO l=j+1,PP_N
    CALL PP_VolumeFluxAverageVec(             U(:,i,j,k,iElem),              U(:,i,l,k,iElem), &
                                           Uaux(:,i,j,k)      ,           Uaux(:,i,l,k)      , &
                                 Metrics_gTilde(:,i,j,k,iElem), Metrics_gTilde(:,i,l,k,iElem), &
                                       gtilde(:,l,i,j,k)                                       )
    gtilde(:,j,i,l,k)=gtilde(:,l,i,j,k) !symmetric
  END DO!l=j+1,N
END DO; END DO; END DO ! i,j,k
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part
  htilde(:,k,i,j,k)=htilde_c(:,i,j,k) 
  DO l=k+1,PP_N
    CALL PP_VolumeFluxAverageVec(             U(:,i,j,k,iElem),              U(:,i,j,l,iElem), &
                                           Uaux(:,i,j,k)      ,           Uaux(:,i,j,l)      , &
                                 Metrics_hTilde(:,i,j,k,iElem), Metrics_hTilde(:,i,j,l,iElem), &
                                       htilde(:,l,i,j,k)                                       )
    htilde(:,k,i,j,l)=htilde(:,l,i,j,k) !symmetric
  END DO!l=k+1,N
END DO; END DO; END DO ! i,j,k

#if NONCONS
CALL AddNonConsFluxTilde3D(iElem,Uaux,ftilde,gtilde,htilde)
#endif /*NONCONS*/

END SUBROUTINE EvalEulerFluxAverage3D_eqn

!==================================================================================================================================
!> Compute MHD transformed fluxes using conservative variables and derivatives for every volume Gauss point.
!> directly apply metrics and output the tranformed flux 
!==================================================================================================================================
SUBROUTINE EvalEulerFluxTilde3D_eqn(U_in,M_f,M_g,M_h,ftilde,gtilde,htilde,Uaux)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:nAuxVar,kappaM1,kappaM2,smu_0,s2mu_0
#ifdef PP_GLM
USE MOD_Equation_vars ,ONLY:GLM_ch
#endif /*PP_GLM*/
USE MOD_DG_Vars       ,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN )   :: U_in(PP_nVar,1:nTotal_vol)
REAL,INTENT(IN )   :: M_f(       3,1:nTotal_vol) 
REAL,INTENT(IN )   :: M_g(       3,1:nTotal_vol) 
REAL,INTENT(IN )   :: M_h(       3,1:nTotal_vol) 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(1:PP_nVar,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(1:PP_nVar,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(1:PP_nVar,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
REAL,INTENT(OUT)   :: Uaux(    nAuxVar,1:nTotal_vol) !< auxiliary variables:(srho,v1,v2,v3,p,|v|^2)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:PP_nVar) :: f,g,h                             ! Cartesian fluxes (iVar)
REAL                :: srho                                    ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,v_2,pt                         ! velocity and pressure(including magnetic pressure
REAL                :: bb2,vb                                  ! magnetic field, bb2=|bvec|^2, v dot b
REAL                :: Ep                                      ! E + p
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  ASSOCIATE(rho   =>U_in(1,i), &
            rhov1 =>U_in(2,i), &
            rhov2 =>U_in(3,i), &
            rhov3 =>U_in(4,i), &
#ifdef PP_GLM
            Etotal=>U_in(5,i)-0.5*smu_0*U_in(9,i)**2, &
#else
            Etotal=>U_in(5,i), &
#endif /*def PP_GLM*/
            b1    =>U_in(6,i), &
            b2    =>U_in(7,i), &
            b3    =>U_in(8,i)  ) 
  ! auxiliary variables
  srho = 1. / rho ! 1/rho
  v1   = rhov1*srho 
  v2   = rhov2*srho 
  v3   = rhov3*srho 
  v_2  = v1*v1+v2*v2+v3*v3 
  bb2  = (b1*b1+b2*b2+b3*b3)
  vb   = (b1*v1+b2*v2+b3*v3)
  !p = ptilde (includes magnetic pressure)
  pt   = kappaM1*(Etotal-0.5*rho*(v_2))-KappaM2*s2mu_0*bb2
  Ep   = (Etotal + pt)
  
  Uaux(:,i)=(/srho,v1,v2,v3,pt,v_2,bb2,vb/)
  ! Advection part
  ! Advection fluxes x-direction
  f(1)=rhov1                     ! rho*u
  f(2)=rhov1*v1+pt -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
  f(3)=rhov1*v2    -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
  f(4)=rhov1*v3    -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
  f(5)=Ep*v1       -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
  f(6)=0.
  f(7)=v1*b2-b1*v2
  f(8)=v1*b3-b1*v3
  ! Advection fluxes y-direction
  g(1)=rhov2                     ! rho*v      
  g(2)=f(3)                      ! rho*u*v      -1/mu_0*b2*b1
  g(3)=rhov2*v2+pt -smu_0*b2*b2  ! rho*v²+p     -1/mu_0*b2*b2
  g(4)=rhov2*v3    -smu_0*b2*b3  ! rho*v*w      -1/mu_0*b2*b3
  g(5)=Ep*v2       -smu_0*b2*vb  ! (rho*e+p)*v  -1/mu_0*b2*(v dot B)
  g(6)=-f(7)                     ! (v2*b1-b2*v1)
  g(7)=0.
  g(8)=v2*b3-b2*v3
  ! Advection fluxes z-direction
  h(1)=rhov3                     ! rho*v
  h(2)=f(4)                      ! rho*u*w      -1/mu_0*b3*b1
  h(3)=g(4)                      ! rho*v*w      -1/mu_0*b3*b2
  h(4)=rhov3*v3+pt -smu_0*b3*b3  ! rho*v²+p     -1/mu_0*b3*b3
  h(5)=Ep*v3       -smu_0*b3*vb  ! (rho*e+p)*w  -1/mu_0*b3*(v dot B)
  h(6)=-f(8)                     ! v3*b1-b3*v1 
  h(7)=-g(8)                     ! v3*b2-b3*v2
  h(8)=0.

#ifdef PP_GLM
  f(5) = f(5)+smu_0*GLM_ch*b1*U_in(9,i)
  f(6) = f(6)+GLM_ch*U_in(9,i)
  f(9) = GLM_ch*b1

  g(5) = g(5)+smu_0*GLM_ch*b2*U_in(9,i)
  g(7) = g(7)+GLM_ch*U_in(9,i)
  g(9) = GLM_ch*b2

  h(5) = h(5)+smu_0*GLM_ch*b3*U_in(9,i)
  h(8) = h(8)+GLM_ch*U_in(9,i)
  h(9) = GLM_ch*b3
#endif /* PP_GLM */

  END ASSOCIATE ! rho,rhov1,rhov2,rhov3,Etotal,b1,b2,b3

  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(:,i) =   f(:)*M_f(1,i) + g(:)*M_f(2,i) + h(:)*M_f(3,i)
  gtilde(:,i) =   f(:)*M_g(1,i) + g(:)*M_g(2,i) + h(:)*M_g(3,i)
  htilde(:,i) =   f(:)*M_h(1,i) + g(:)*M_h(2,i) + h(:)*M_h(3,i)
END DO ! i
END SUBROUTINE EvalEulerFluxTilde3D_eqn


#if NONCONS
!==================================================================================================================================
!> Compute transformed nonconservative MHD fluxes 
!==================================================================================================================================
SUBROUTINE AddNonConsFluxTilde3D(iElem,Uaux,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars ,ONLY:nAuxVar
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: iElem  !< Determines the actual element
REAL   ,INTENT(IN)           :: Uaux(nAuxVar,0:PP_N,0:PP_N,0:PP_N)   !auxiliary variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: ftilde !< add to transformed flux f(iVar,l,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: gtilde !< add to transformed flux g(iVar,l,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: htilde !< add to transformed flux h(iVar,l,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,l
!REAL :: phi(PP_nVar) 
REAL :: phi_s4(2:8) 
#ifdef PP_GLM
REAL :: phi_GLM_f_s2(2),phi_GLM_g_s2(2),phi_GLM_h_s2(2) 
#endif
!==================================================================================================================================
!phi=0.
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  ! Powell
  Phi_s4(2:4)=0.25* U(   6:8,i,j,k,iElem) ! B
  Phi_s4(  5)=0.25* Uaux(  8,i,j,k)      ! vB
  Phi_s4(6:8)=0.25* Uaux(2:4,i,j,k)    ! v
#ifdef PP_GLM
  phi_GLM_f_s2(1:2) = (0.5*SUM(Metrics_ftilde(:,i,j,k,iElem)*Uaux(2:4,i,j,k)))*(/U(9,i,j,k,iElem),1./)
  phi_GLM_g_s2(1:2) = (0.5*SUM(Metrics_gtilde(:,i,j,k,iElem)*Uaux(2:4,i,j,k)))*(/U(9,i,j,k,iElem),1./)
  phi_GLM_h_s2(1:2) = (0.5*SUM(Metrics_htilde(:,i,j,k,iElem)*Uaux(2:4,i,j,k)))*(/U(9,i,j,k,iElem),1./)
#endif /*PP_GLM*/
  DO l=0,PP_N
    ftilde(    2:8,l,i,j,k) = ftilde(2:8,l,i,j,k)+(SUM(( Metrics_ftilde(:,i,j,k,iElem) &
                                                        +Metrics_ftilde(:,l,j,k,iElem))*U(6:8,l,j,k,iElem)))*Phi_s4(2:8) 
#ifdef PP_GLM
    !nonconservative term to restore galilein invariance for GLM term: (grad\psi) (0,0,0,0,vec{v}\psi, 0,0,0, \vec{v})
    ! => 5/9. component: 1/2 vec{Ja^d}_{(l,i),jk} . vec{v}_ijk \psi_ljk (\psi_ijk,1)
    !
    ftilde((/5,9/),l,i,j,k) = ftilde((/5,9/),l,i,j,k)                +U(9,l,j,k,iElem) *phi_GLM_f_s2(1:2)
#endif /*PP_GLM*/
  END DO !l=0,PP_N
  DO l=0,PP_N
    gtilde(    2:8,l,i,j,k) = gtilde(2:8,l,i,j,k)+(SUM(( Metrics_gtilde(:,i,j,k,iElem) & 
                                                        +Metrics_gtilde(:,i,l,k,iElem))*U(6:8,i,l,k,iElem)))*Phi_s4(2:8) 
#ifdef PP_GLM
    gtilde((/5,9/),l,i,j,k) = gtilde((/5,9/),l,i,j,k)                +U(9,i,l,k,iElem) *phi_GLM_g_s2(1:2)
#endif /*PP_GLM*/
  END DO !l=0,PP_N
  DO l=0,PP_N
    htilde(2:8,    l,i,j,k) = htilde(2:8,l,i,j,k)+(SUM(( Metrics_htilde(:,i,j,k,iElem) & 
                                                        +Metrics_htilde(:,i,j,l,iElem))*U(6:8,i,j,l,iElem)))*Phi_s4(2:8) 
#ifdef PP_GLM
    htilde((/5,9/),l,i,j,k) = htilde((/5,9/),l,i,j,k)                +U(9,i,j,l,iElem) *phi_GLM_h_s2(1:2)
#endif /*PP_GLM*/
  END DO !l=0,PP_N
END DO; END DO; END DO ! i,j,k

END SUBROUTINE AddNonConsFluxTilde3D

!==================================================================================================================================
!> Compute transformed nonconservative MHD flux given left and right states and metrics 
!==================================================================================================================================
SUBROUTINE AddNonConsFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
REAL,DIMENSION(8),INTENT(IN)        :: UauxL          !< left auxiliary variables
REAL,DIMENSION(8),INTENT(IN)        :: UauxR          !< right auxiliary variables
REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(INOUT) :: Fstar   !< added to flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
  ! Powell
  !Phi(2:4)=UL(6:8)
  !Phi(5)  =UauxL(8)      ! vB
  !Phi(6:8)=UauxL(2:4)    ! v
  Fstar(2:8) = Fstar(2:8) &
               +(0.25*SUM((metric_L(:)+metric_R(:))*UR(6:8))) &
                *(/UL(6:8),UauxL(8),UauxL(2:4)/)
#ifdef PP_GLM
  !nonconservative term to restore galilein invariance for GLM term, 1/2 cancels with 2*Dmat 
  ! grad\psi (0,0,0,0,vec{v}\phi, 0,0,0, \vec{v}) => vec{Ja^d}_{i,j,k} . vec{v}_ijk \psi_l,j,k
  Fstar((/5,9/)) = Fstar((/5,9/))  &
                   +(UR(9)*0.5*SUM(metric_L(:)*UauxL(2:4))) &
                   *(/UL(9),1./)
#endif /*PP_GLM*/

END SUBROUTINE AddNonConsFluxVec
#endif /*NONCONS*/

!==================================================================================================================================
!> computes auxiliary nodal variables (1/rho,v_1,v_2,v_3,p_t,|v|^2) 
!==================================================================================================================================
SUBROUTINE EvalUaux(iElem,Uaux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Equation_Vars ,ONLY:nAuxVar
USE MOD_Equation_Vars ,ONLY:kappaM1,KappaM2,s2mu_0
#ifdef OPTIMIZED
USE MOD_DG_Vars,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                        :: iElem !< current element index in volint
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(nAuxVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: Uaux   !< auxiliary variables:(srho,v1,v2,v3,p_t,|v|^2,|B|^2,v*b
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i 
#ifndef OPTIMIZED
INTEGER             :: j,k
#endif
!==================================================================================================================================
#ifdef OPTIMIZED
DO i=0,nTotal_vol-1
#else /*OPTIMIZED*/
DO k=0,PP_N;  DO j=0,PP_N; DO i=0,PP_N
#endif /*OPTIMIZED*/
  ! auxiliary variables
  Uaux(1  ,PP_IJK) = 1./U(1,PP_IJK,iElem)                      ! 1/rho
  Uaux(2:4,PP_IJK) = Uaux(1,PP_IJK)*U(2:4,PP_IJK,iElem)        ! vec{rho*v}/rho
  Uaux(6  ,PP_IJK) = SUM(Uaux(2:4,PP_IJK)*Uaux(2:4,PP_IJK))                  ! |v|^2
  Uaux(7  ,PP_IJK)  =SUM(U(6:8,PP_IJK,iElem)**2)               ! |B|^2
  Uaux(8  ,PP_IJK)  =SUM(Uaux(2:4,PP_IJK)*U(6:8,PP_IJK,iElem)) ! v*B
  !total pressure=gas pressure + magnetic pressure
  Uaux(5  ,PP_IJK)=kappaM1*(U(5,PP_IJK,iElem) -0.5*U(1,PP_IJK,iElem)*Uaux(6,PP_IJK) &
#ifdef PP_GLM
                                              -s2mu_0*U(9,PP_IJK,iElem)**2 &
#endif /*PP_GLM*/
                                                    )-kappaM2*s2mu_0*Uaux(7,PP_IJK) !p_t 
#ifdef OPTIMIZED
END DO ! i
#else /*OPTIMIZED*/
END DO; END DO; END DO ! i,j,k
#endif /*OPTIMIZED*/
END SUBROUTINE EvalUaux


!==================================================================================================================================
!> General FluxMat routine for  VolInt_SplitForm_eqn, which is basically a wrapper for the FluxVec point-wise operations, 
!> since not all flux averages have been implemented in the optimized FluxMat way. 
!> uses PP_VolumeFluxAverageVec, being either the pointer the flux averaging routine set at runtime, or a precompiled routine 
!> defined by PP_VolFlux parameter.
!==================================================================================================================================
PURE SUBROUTINE GeneralFluxMat(U_in,metric_in,Fstar) 
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,kappaM2,smu_0,s2mu_0
#ifdef PP_GLM
USE MOD_Equation_vars ,ONLY:GLM_ch
#endif /*PP_GLM*/
#if PP_VolFlux==-1
USE MOD_Equation_Vars  ,ONLY:VolumeFluxAverageVec !pointer to flux averaging routine
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: U_in(PP_nVar,0:PP_N)   !< right state
REAL,INTENT(IN) :: metric_in(3,0:PP_N)   !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: Fstar(PP_nVar,0:PP_N,0:PP_N)          !< transformed flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: l,i
REAL,DIMENSION(0:PP_N)            :: qvHat,qBHat,srho,v1,v2,v3,b1,b2,b3,pt,vv,vb,bb
#ifdef PP_GLM
REAL,DIMENSION(0:PP_N)              :: Etotal
#endif 
#if NONCONS
REAL                                  :: btilde
#endif /*NONCONS*/
!==================================================================================================================================
srho(:)  = 1./U_in(1,:)
v1(:)    = U_in(2,:)*srho(:)
v2(:)    = U_in(3,:)*srho(:)
v3(:)    = U_in(4,:)*srho(:)
b1(:)    = U_in(6,:)
b2(:)    = U_in(7,:)
b3(:)    = U_in(8,:)
vv(:)    = v1(:)*v1(:)+v2(:)*v2(:)+v3(:)*v3(:)
vb(:)    = v1(:)*b1(:)+v2(:)*b2(:)+v3(:)*b3(:)
bb(:)    = b1(:)*b1(:)+b2(:)*b2(:)+b3(:)*b3(:)
#ifdef PP_GLM
Etotal(:)=U_in(5,:)-s2mu_0*U_in(9,:)**2
#else
Etotal(:)=U_in(5,:)
#endif /*def PP_GLM*/
pt(:)    = kappaM1*(Etotal(:)-0.5*U_in(1,:)*vv(:)) &
           -kappaM2*s2mu_0*bb(:) ! = -kappaM1*s2mu_0*bb + s2mu_0*bb, include magn. press
qvHat(:) =   v1(:)*metric_in(1,:) +   v2(:)*metric_in(2,:) +   v3(:)*metric_in(3,:) 
qbHat(:) =   b1(:)*metric_in(1,:) +   b2(:)*metric_in(2,:) +   b3(:)*metric_in(3,:) 

DO i=0,PP_N
  !consistency euler flux f(Ui,Ui)=f(Ui)
  Fstar(1,i,i) = (U_in(1,i)*qvHat(i)  )
  Fstar(2,i,i) = (U_in(2,i)*qvHat(i)  + metric_in(1,i)*pt(i) -smu_0*(qbHat(i)*b1(i)) )
  Fstar(3,i,i) = (U_in(3,i)*qvHat(i)  + metric_in(2,i)*pt(i) -smu_0*(qbHat(i)*b2(i)) )
  Fstar(4,i,i) = (U_in(4,i)*qvHat(i)  + metric_in(3,i)*pt(i) -smu_0*(qbHat(i)*b3(i)) )

#ifdef PP_GLM
  Fstar(5,i,i) = (Etotal(i) + pt(i))*qvHat(i) -smu_0*(qbHat(i)*vb(i)) + GLM_ch*U_in(9,i)*qbHat(i)             
  Fstar(6,i,i) = (qvHat(i)*b1(i)-qbHat(i)*v1(i))                      + GLM_ch*U_in(9,i)*metric_in(1,i) 
  Fstar(7,i,i) = (qvHat(i)*b2(i)-qbHat(i)*v2(i))                      + GLM_ch*U_in(9,i)*metric_in(2,i) 
  Fstar(8,i,i) = (qvHat(i)*b3(i)-qbHat(i)*v3(i))                      + GLM_ch*U_in(9,i)*metric_in(3,i) 
  Fstar(9,i,i) =                                                        GLM_ch* qbHat(i)                   
#else
  Fstar(5,i,i) = (Etotal(i) + pt(i))*qvHat(i) -smu_0*(qbHat(i)*vb(i)) 
  Fstar(6,i,i) = (qvHat(i)*b1(i)-qbHat(i)*v1(i))
  Fstar(7,i,i) = (qvHat(i)*b2(i)-qbHat(i)*v2(i))
  Fstar(8,i,i) = (qvHat(i)*b3(i)-qbHat(i)*v3(i))
#endif /* PP_GLM */

  DO l=i+1,PP_N
    ! call fluxvec (general but not the fastest) 
    CALL PP_VolumeFluxAverageVec(  U_in(:,i)   ,      U_in(:,l), &
                                   (/srho(i),v1(i),v2(i),v3(i),pt(i),vv(i),bb(i),vb(i)/) , & !Uaux(i)
                                   (/srho(l),v1(l),v2(l),v3(l),pt(l),vv(l),bb(l),vb(l)/) , & !Uaux(l)
                                 Metric_in(:,i), Metric_in(:,l), &
                                 Fstar(:,l,i)                                       )
!    !symmetry
    Fstar(:,i,l)=Fstar(:,l,i)
  END DO !l=i+1,PP_N
END DO !i=0,PP_N
!add nonconservative (non-symmetric!) term
#if NONCONS
#ifdef PP_GLM
qvHat(:)    =   0.5*qvHat(:)
#endif /*PP_GLM*/
DO i=0,PP_N
  DO l=0,PP_N
    btilde=0.25*(  (metric_in(1,i)+metric_in(1,l))*b1(i)  &
                  +(metric_in(2,i)+metric_in(2,l))*b2(i)  &
                  +(metric_in(3,i)+metric_in(3,l))*b3(i)  )

    !powell term
    Fstar(2,l,i) = Fstar(2,l,i)+ btilde*b1(i)
    Fstar(3,l,i) = Fstar(3,l,i)+ btilde*b2(i)
    Fstar(4,l,i) = Fstar(4,l,i)+ btilde*b3(i)
    Fstar(6,l,i) = Fstar(6,l,i)+ btilde*v1(i)
    Fstar(7,l,i) = Fstar(7,l,i)+ btilde*v2(i)
    Fstar(8,l,i) = Fstar(8,l,i)+ btilde*v3(i)
#ifdef PP_GLM
    !nonconservative term to restore galilein invariance for GLM term, 1/2 cancels with 2*Dmat 
    ! grad\psi (0,0,0,0,vec{v}\phi, 0,0,0, \vec{v}) => vec{Ja^d}_{i,j,k} . vec{v}_ijk \psi_l,j,k
    Fstar(5,l,i) = Fstar(5,l,i)+ btilde*vb(i) + U_in(9,l)*qvHat(i)*U_in(9,i)
    Fstar(9,l,i) = Fstar(9,l,i)               + U_in(9,l)*qvHat(i)
#else
    Fstar(5,l,i) = Fstar(5,l,i)+ btilde*vb(i)
#endif /*PP_GLM*/
  END DO !l=0,PP_N
END DO !i=0,PP_N
#endif /*NONCONS*/

END SUBROUTINE GeneralFluxMat


!==================================================================================================================================
!> Computes the standard flux in x-direction for the hyperbolic part ( normally used with a rotated state)
!==================================================================================================================================
PURE SUBROUTINE StandardDGFlux(UL,UR,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,kappaM2,smu_0,s2mu_0
#ifdef PP_GLM
USE MOD_Equation_Vars,ONLY:GLM_ch
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar !< 1D flux in x direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: rhoqL,rhoqR
REAL                                :: sRho_L,sRho_R,v1_L,v2_L,v3_L,v1_R,v2_R,v3_R,bb2_L,bb2_R,vb_L,vb_R,pt_L,pt_R
!==================================================================================================================================
! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>UL(1),  rho_R =>UR(1), &
           rhoU_L =>UL(2), rhoU_R =>UR(2), &
           rhoV_L =>UL(3), rhoV_R =>UR(3), &
           rhoW_L =>UL(4), rhoW_R =>UR(4), &
#ifdef PP_GLM
              E_L =>UL(5)-0.5*smu_0*UL(9)**2, E_R =>UR(5)-0.5*smu_0*UR(9)**2, &
#else
              E_L =>UL(5),    E_R =>UR(5), &
#endif
             b1_L =>UL(6),   b1_R =>UR(6), &
             b2_L =>UL(7),   b2_R =>UR(7), &
             b3_L =>UL(8),   b3_R =>UR(8)  )
sRho_L = 1./rho_L
sRho_R = 1./rho_R
v1_L = sRho_L*rhoU_L;  v2_L = sRho_L*rhoV_L ; v3_L = sRho_L*rhoW_L
v1_R = sRho_R*rhoU_R;  v2_R = sRho_R*rhoV_R ; v3_R = sRho_R*rhoW_R
bb2_L  = (b1_L*b1_L+b2_L*b2_L+b3_L*b3_L)
bb2_R  = (b1_R*b1_R+b2_R*b2_R+b3_R*b3_R)

pt_L    = (kappaM1)*(E_L - 0.5*(rhoU_L*v1_L + rhoV_L*v2_L + rhoW_L*v3_L))-kappaM2*s2mu_0*bb2_L
pt_R    = (kappaM1)*(E_R - 0.5*(rhoU_R*v1_R + rhoV_R*v2_R + rhoW_R*v3_R))-kappaM2*s2mu_0*bb2_R

vb_L  = v1_L*b1_L+v2_L*b2_L+v3_L*b3_L
vb_R  = v1_R*b1_R+v2_R*b2_R+v3_R*b3_R
  
! Standard DG flux
rhoqL    = rho_L*v1_L
rhoqR    = rho_R*v1_R
Fstar(1) = 0.5*(rhoqL      + rhoqR)
Fstar(2) = 0.5*(rhoqL*v1_L + rhoqR*v1_R +(pt_L + pt_R)-smu_0*(b1_L*b1_L+b1_R*b1_R))
Fstar(3) = 0.5*(rhoqL*v2_L + rhoqR*v2_R               -smu_0*(b2_L*b1_L+b2_R*b1_R))
Fstar(4) = 0.5*(rhoqL*v3_L + rhoqR*v3_R               -smu_0*(b3_L*b1_L+b3_R*b1_R))
Fstar(5) = 0.5*((E_L + pt_L)*v1_L + (E_R + pt_R)*v1_R- smu_0*(b1_L*vb_L+b1_R*vb_R))
Fstar(6) = 0.
Fstar(7) = 0.5*(v1_L*b2_L-b1_L*v2_L + v1_R*b2_R-b1_R*v2_R)
Fstar(8) = 0.5*(v1_L*b3_L-b1_L*v3_L + v1_R*b3_R-b1_R*v3_R)
#ifdef PP_GLM
Fstar(5) = Fstar(5)+0.5*smu_0*GLM_ch*(b1_L*UL(9)+b1_R*UR(9))
Fstar(6) = Fstar(6)+0.5      *GLM_ch*(     UL(9)+     UR(9))
Fstar(9) =          0.5      *GLM_ch*(b1_L      +b1_R      )
#endif /* PP_GLM */
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE StandardDGFlux


!==================================================================================================================================
!> Computes the standard DG flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the advection
!> part of the MHD equations
!> for curved metrics, no dealiasing is done (exactly = standard DG )!
!==================================================================================================================================
PURE SUBROUTINE StandardDGFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:smu_0
#ifdef PP_GLM
USE MOD_Equation_vars ,ONLY:GLM_ch
#endif /*PP_GLM*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
REAL,DIMENSION(8),INTENT(IN)        :: UauxL          !< left auxiliary variables
REAL,DIMENSION(8),INTENT(IN)        :: UauxR          !< right auxiliary variables
REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< transformed central flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: qv_L,qv_R,qb_L,qb_R
!==================================================================================================================================

! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
#ifdef PP_GLM
             E_L =>UL(5)-0.5*smu_0*UL(9)**2, E_R =>UR(5)-0.5*smu_0*UR(9)**2, &
#else
             E_L =>UL(5), E_R =>UR(5), &
#endif
             b1_L =>   UL(6),   b1_R =>   UR(6), &
             b2_L =>   UL(7),   b2_R =>   UR(7), &
             b3_L =>   UL(8),   b3_R =>   UR(8), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
             v1_L =>UauxL(2),   v1_R =>UauxR(2), &
             v2_L =>UauxL(3),   v2_R =>UauxR(3), &
             v3_L =>UauxL(4),   v3_R =>UauxR(4), &
             pt_L =>UauxL(5),   pt_R =>UauxR(5), & !total pressure = gas pressure+magnetic pressure
           ! v2_L =>UauxL(6),   v2_R =>UauxR(6), &
           ! b2_L =>UauxL(7),   b2_R =>UauxR(7), &
             vb_L =>UauxL(8),   vb_R =>UauxR(8)  )

!without metric dealiasing (=standard DG weak form on curved meshes)
qv_L = v1_L*metric_L(1) + v2_L*metric_L(2) + v3_L*metric_L(3)
qb_L = b1_L*metric_L(1) + b2_L*metric_L(2) + b3_L*metric_L(3)

qv_R = v1_R*metric_R(1) + v2_R*metric_R(2) + v3_R*metric_R(3)
qb_R = b1_R*metric_R(1) + b2_R*metric_R(2) + b3_R*metric_R(3)

! Standard DG flux
!without metric dealiasing (=standard DG weak form on curved meshes)
Fstar(1) = 0.5*( rho_L*qv_L +  rho_R*qv_R )
Fstar(2) = 0.5*(rhoU_L*qv_L + rhoU_R*qv_R + metric_L(1)*pt_L+metric_R(1)*pt_R -smu_0*(qb_L*b1_L+qb_R*b1_R) )
Fstar(3) = 0.5*(rhoV_L*qv_L + rhoV_R*qv_R + metric_L(2)*pt_L+metric_R(2)*pt_R -smu_0*(qb_L*b2_L+qb_R*b2_R) )
Fstar(4) = 0.5*(rhoW_L*qv_L + rhoW_R*qv_R + metric_L(3)*pt_L+metric_R(3)*pt_R -smu_0*(qb_L*b3_L+qb_R*b3_R) )
Fstar(5) = 0.5*((E_L + pt_L)*qv_L  + (E_R + pt_R)*qv_R      -smu_0*(qb_L*vb_L+qb_R*vb_R) )
Fstar(6) = 0.5*(qv_L*b1_L-qb_L*v1_L + qv_R*b1_R-qb_R*v1_R)
Fstar(7) = 0.5*(qv_L*b2_L-qb_L*v2_L + qv_R*b2_R-qb_R*v2_R)
Fstar(8) = 0.5*(qv_L*b3_L-qb_L*v3_L + qv_R*b3_R-qb_R*v3_R)

#ifdef PP_GLM
!without metric dealiasing (=standard DG weak form on curved meshes)
Fstar(5) = Fstar(5) + 0.5*GLM_ch*(qb_L*UL(9)             + qb_R*UR(9))
Fstar(6) = Fstar(6) + 0.5*GLM_ch*(     UL(9)*metric_L(1) +      UR(9)*metric_R(1))
Fstar(7) = Fstar(7) + 0.5*GLM_ch*(     UL(9)*metric_L(2) +      UR(9)*metric_R(2))
Fstar(8) = Fstar(8) + 0.5*GLM_ch*(     UL(9)*metric_L(3) +      UR(9)*metric_R(3))
Fstar(9) =            0.5*GLM_ch*(qb_L                   +qb_R                   )

#endif /* PP_GLM */

END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE StandardDGFluxVec

!==================================================================================================================================
!> Computes the standard DG euler flux transformed with the metrics 
!> fstar=1/2((fL*metric1L+gL*metric2L+h*metric3L)+(fR*metric1R+gR*metric2R+h*metric3R)  )
!==================================================================================================================================
PURE SUBROUTINE standardDGFluxMat(U_in,metric_in,Fstar) 
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,kappaM2,smu_0,s2mu_0
#ifdef PP_GLM
USE MOD_Equation_vars ,ONLY:GLM_ch
#endif /*PP_GLM*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: U_in(PP_nVar,1:PP_N+1)   !< right state
REAL,INTENT(IN) :: metric_in(3,1:PP_N+1)   !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: Fstar(PP_nVar,1:PP_N+1,1:PP_N+1)          !< transformed flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                             :: l,i
REAL,DIMENSION(1:PP_N+1)              :: qvHat,qBHat,srho,v1,v2,v3,pt,vb,bb,flux1,flux2,flux3,flux4,flux5,flux6,flux7,flux8 
#ifdef PP_GLM
REAL,DIMENSION(1:PP_N+1)              :: Etotal,flux9
#endif 
#if NONCONS
REAL                                  :: btilde
#endif /*NONCONS*/
!==================================================================================================================================
ASSOCIATE(  rhov1 =>U_in(2,:), &
            rhov2 =>U_in(3,:), &
            rhov3 =>U_in(4,:), &
#ifndef PP_GLM
            Etotal=>U_in(5,:), &
#endif
            b1    =>U_in(6,:), &
            b2    =>U_in(7,:), &
            b3    =>U_in(8,:)  ) 
srho(:)  = 1./U_in(1,:)
v1(:)    = rhov1*srho(:)
v2(:)    = rhov2*srho(:)
v3(:)    = rhov3*srho(:)
vb(:)    = v1(:)*b1(:)+v2(:)*b2(:)+v3(:)*b3(:)
bb(:)    = b1(:)*b1(:)+b2(:)*b2(:)+b3(:)*b3(:)
#ifdef PP_GLM
Etotal(:)=U_in(5,:)-s2mu_0*U_in(9,:)**2
#endif /*def PP_GLM*/
pt(:)    = kappaM1*(Etotal(:)-0.5*(rhov1(:)*v1(:)+rhov2(:)*v2(:)+rhov3(:)*v3(:))) &
           -kappaM2*s2mu_0*bb(:) ! = -kappaM1*s2mu_0*bb + s2mu_0*bb, include magn. press
qvHat(:) =   v1(:)*metric_in(1,:) +   v2(:)*metric_in(2,:) +   v3(:)*metric_in(3,:) 
qbHat(:) =   b1(:)*metric_in(1,:) +   b2(:)*metric_in(2,:) +   b3(:)*metric_in(3,:) 

Flux1(:) = (U_in(1,:)*qvHat(:)  )
Flux2(:) = ( rhov1(:)*qvHat(:)  + metric_in(1,:)*pt(:) -smu_0*(qbHat(:)*b1(:)) )
Flux3(:) = ( rhov2(:)*qvHat(:)  + metric_in(2,:)*pt(:) -smu_0*(qbHat(:)*b2(:)) )
Flux4(:) = ( rhov3(:)*qvHat(:)  + metric_in(3,:)*pt(:) -smu_0*(qbHat(:)*b3(:)) )

#ifdef PP_GLM
Flux5(:) = (Etotal(:) + pt(:))*qvHat(:) -smu_0*(qbHat(:)*vb(:)) + GLM_ch*U_in(9,:)*qbHat(:)             
Flux6(:) = (qvHat(:)*b1(:)-qbHat(:)*v1(:))                      + GLM_ch*U_in(9,:)*metric_in(1,:) 
Flux7(:) = (qvHat(:)*b2(:)-qbHat(:)*v2(:))                      + GLM_ch*U_in(9,:)*metric_in(2,:) 
Flux8(:) = (qvHat(:)*b3(:)-qbHat(:)*v3(:))                      + GLM_ch*U_in(9,:)*metric_in(3,:) 
Flux9(:) =                                                        GLM_ch* qbHat(:)                   
#else
Flux5(:) = (Etotal(:) + pt(:))*qvHat(:) -smu_0*(qbHat(:)*vb(:)) 
Flux6(:) = (qvHat(:)*b1(:)-qbHat(:)*v1(:))
Flux7(:) = (qvHat(:)*b2(:)-qbHat(:)*v2(:))
Flux8(:) = (qvHat(:)*b3(:)-qbHat(:)*v3(:))
#endif /* PP_GLM */


DO i=1,PP_N+1
  !consistency euler flux f(Ui,Ui)=f(Ui)
  Fstar(1,i,i) = flux1(i)
  Fstar(2,i,i) = flux2(i)
  Fstar(3,i,i) = flux3(i)
  Fstar(4,i,i) = flux4(i)
  Fstar(5,i,i) = flux5(i)
  Fstar(6,i,i) = flux6(i)
  Fstar(7,i,i) = flux7(i)
  Fstar(8,i,i) = flux8(i)
#ifdef PP_GLM
  Fstar(9,i,i) = flux9(i)
#endif
  DO l=i+1,PP_N+1
    ! 
    Fstar(1,l,i) = 0.5*(flux1(l)+flux1(i))
    Fstar(2,l,i) = 0.5*(flux2(l)+flux2(i))
    Fstar(3,l,i) = 0.5*(flux3(l)+flux3(i))
    Fstar(4,l,i) = 0.5*(flux4(l)+flux4(i))
    Fstar(5,l,i) = 0.5*(flux5(l)+flux5(i))
    Fstar(6,l,i) = 0.5*(flux6(l)+flux6(i))
    Fstar(7,l,i) = 0.5*(flux7(l)+flux7(i))
    Fstar(8,l,i) = 0.5*(flux8(l)+flux8(i))
#ifdef PP_GLM
    Fstar(9,l,i) = 0.5*(flux9(l)+flux9(i))
#endif
    !symmetry
    Fstar(:,i,l)=Fstar(:,l,i)
  END DO !l=i+1,PP_N+1
END DO !i=1,PP_N+1

!add nonconservative (non-symmetric!) term
#if NONCONS
#ifdef PP_GLM
qvHat(:)    =   0.5*qvHat(:)
#endif /*PP_GLM*/
DO i=1,PP_N+1
  DO l=1,PP_N+1
    btilde=0.25*(  (metric_in(1,i)+metric_in(1,l))*b1(i)  &
                  +(metric_in(2,i)+metric_in(2,l))*b2(i)  &
                  +(metric_in(3,i)+metric_in(3,l))*b3(i)  )

    !powell term
    Fstar(2,l,i) = Fstar(2,l,i)+ btilde*b1(i)
    Fstar(3,l,i) = Fstar(3,l,i)+ btilde*b2(i)
    Fstar(4,l,i) = Fstar(4,l,i)+ btilde*b3(i)
    Fstar(6,l,i) = Fstar(6,l,i)+ btilde*v1(i)
    Fstar(7,l,i) = Fstar(7,l,i)+ btilde*v2(i)
    Fstar(8,l,i) = Fstar(8,l,i)+ btilde*v3(i)
#ifdef PP_GLM
    !nonconservative term to restore galilein invariance for GLM term, 1/2 cancels with 2*Dmat 
    ! grad\psi (0,0,0,0,vec{v}\phi, 0,0,0, \vec{v}) => vec{Ja^d}_{i,j,k} . vec{v}_ijk \psi_l,j,k
    Fstar(5,l,i) = Fstar(5,l,i)+ btilde*vb(i) + U_in(9,l)*qvHat(i)*U_in(9,i)
    Fstar(9,l,i) = Fstar(9,l,i)               + U_in(9,l)*qvHat(i)
#else
    Fstar(5,l,i) = Fstar(5,l,i)+ btilde*vb(i)
#endif /*PP_GLM*/
  END DO !l=1,PP_N+1
END DO !i=1,PP_N+1
#endif /*NONCONS*/
END ASSOCIATE

END SUBROUTINE standardDGFluxMat

!==================================================================================================================================
!> Computes the standard DG flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the advection
!> part of the MHD equations
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
PURE SUBROUTINE StandardDGFluxDealiasedMetricVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:smu_0
#ifdef PP_GLM
USE MOD_Equation_vars ,ONLY:GLM_ch
#endif /*PP_GLM*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
REAL,DIMENSION(8),INTENT(IN)        :: UauxL          !< left auxiliary variables
REAL,DIMENSION(8),INTENT(IN)        :: UauxR          !< right auxiliary variables
REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< transformed central flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: qv_L,qv_R,qb_L,qb_R
#ifdef PP_GLM
REAL                                :: phiHat
#endif /*PP_GLM*/
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)

! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
#ifdef PP_GLM
             E_L =>UL(5)-0.5*smu_0*UL(9)**2, E_R =>UR(5)-0.5*smu_0*UR(9)**2, &
#else
             E_L =>UL(5), E_R =>UR(5), &
#endif
             b1_L =>   UL(6),   b1_R =>   UR(6), &
             b2_L =>   UL(7),   b2_R =>   UR(7), &
             b3_L =>   UL(8),   b3_R =>   UR(8), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
             v1_L =>UauxL(2),   v1_R =>UauxR(2), &
             v2_L =>UauxL(3),   v2_R =>UauxR(3), &
             v3_L =>UauxL(4),   v3_R =>UauxR(4), &
             pt_L =>UauxL(5),   pt_R =>UauxR(5), & !total pressure = gas pressure+magnetic pressure
           !vv2_L =>UauxL(6),  vv2_R =>UauxR(6), &
           !bb2_L =>UauxL(7),  bb2_R =>UauxR(7), &
             vb_L =>UauxL(8),   vb_R =>UauxR(8)  )

qv_L = v1_L*metric(1) + v2_L*metric(2) + v3_L*metric(3)
qb_L = b1_L*metric(1) + b2_L*metric(2) + b3_L*metric(3)

qv_R = v1_R*metric(1) + v2_R*metric(2) + v3_R*metric(3)
qb_R = b1_R*metric(1) + b2_R*metric(2) + b3_R*metric(3)

! Standard DG flux
Fstar(1) = 0.5*( rho_L*qv_L +  rho_R*qv_R )
Fstar(2) = 0.5*(rhoU_L*qv_L + rhoU_R*qv_R + metric(1)*(pt_L+pt_R) -smu_0*(qb_L*b1_L+qb_R*b1_R) )
Fstar(3) = 0.5*(rhoV_L*qv_L + rhoV_R*qv_R + metric(2)*(pt_L+pt_R) -smu_0*(qb_L*b2_L+qb_R*b2_R) )
Fstar(4) = 0.5*(rhoW_L*qv_L + rhoW_R*qv_R + metric(3)*(pt_L+pt_R) -smu_0*(qb_L*b3_L+qb_R*b3_R) )
Fstar(5) = 0.5*((E_L + pt_L)*qv_L  + (E_R + pt_R)*qv_R      -smu_0*(qb_L*vb_L+qb_R*vb_R) )
Fstar(6) = 0.5*(qv_L*b1_L-qb_L*v1_L + qv_R*b1_R-qb_R*v1_R)
Fstar(7) = 0.5*(qv_L*b2_L-qb_L*v2_L + qv_R*b2_R-qb_R*v2_R)
Fstar(8) = 0.5*(qv_L*b3_L-qb_L*v3_L + qv_R*b3_R-qb_R*v3_R)

#ifdef PP_GLM
Fstar(5) = Fstar(5) + 0.5*GLM_ch*(qb_L*UL(9)+qb_R*UR(9))
phiHat   = 0.5*GLM_ch*(UL(9)+UR(9))
Fstar(6) = Fstar(6) + phiHat*metric(1)
Fstar(7) = Fstar(7) + phiHat*metric(2)
Fstar(8) = Fstar(8) + phiHat*metric(3)
Fstar(9) =            0.5*GLM_ch*(qb_L+qb_R)
#endif /* PP_GLM */

END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE StandardDGFluxDealiasedMetricVec



!==================================================================================================================================
!> entropy conservation for MHD, kinetric Energy conservation only in the Euler case
!> following D.Dergs et al."a novel Entropy consistent nine-wave field divergence diminishing ideal MHD system" 
!> mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t) 
!==================================================================================================================================
PURE SUBROUTINE EntropyAndKinEnergyConservingFlux(UL,UR,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,skappaM1,smu_0
#ifdef PP_GLM
USE MOD_Equation_Vars,ONLY:GLM_ch
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: betaLN,beta_R,beta_L
REAL            :: rhoLN,B2_L,B2_R,v2_L,v2_R
REAL            :: pTilde,p_L,p_R
REAL            :: v_L(3),v_R(3)
REAL            :: BAvg(3),vAvg(3)
REAL            :: v1_B2Avg
REAL            :: vB_Avg
#ifdef PP_GLM
REAL            :: psiAvg
#endif
!==================================================================================================================================
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoV_L => UL(2:4), rhoV_R => UR(2:4), &
#ifdef PP_GLM
              E_L =>UL(5)-0.5*smu_0*UL(9)**2, E_R =>UR(5)-0.5*smu_0*UR(9)**2, &
            psi_L =>UL(9)   ,  psi_R =>UR(9), &
#else
              E_L =>UL(5)   ,    E_R =>UR(5), &
#endif
              B_L => UL(6:8),    B_R => UR(6:8)  )
! Get the inverse density, velocity, and pressure on left and right
v_L = rhoV_L(:)/rho_L
v_R = rhoV_R(:)/rho_R

v2_L = SUM(v_L(:)*v_L(:))
v2_R = SUM(v_R(:)*v_R(:))
B2_L = SUM(B_L(:)*B_L(:))
B2_R = SUM(B_R(:)*B_R(:))

!beta=rho/(2*p)
p_L    = kappaM1*(E_L - 0.5*(rho_L*v2_L+smu_0*B2_L))
p_R    = kappaM1*(E_R - 0.5*(rho_R*v2_R+smu_0*B2_R))
beta_L = 0.5*rho_L/p_L
beta_R = 0.5*rho_R/P_R

! Get the averages for the numerical flux

rhoLN      = LN_MEAN( rho_L, rho_R)
betaLN     = LN_MEAN(beta_L,beta_R)
vAvg       = 0.5 * ( v_L +  v_R)
BAvg       = 0.5 * ( B_L +  B_R)
!B2Avg      = 0.5 * (B2_L + B2_R)
v1_B2Avg   = 0.5 * (v_L(1)*B2_L              + v_R(1)*B2_R)
vB_Avg     = 0.5 * (SUM(V_L(:)*B_L(:))+ SUM(V_R(:)*B_R(:)))
                                                                   
pTilde     = 0.5*((rho_L+rho_R)/(beta_L+beta_R)+smu_0*0.5*(B2_L+B2_R)) !rhoLN/(2*betaLN)+1/(2mu_0)({{|B|^2}}...)
#ifdef PP_GLM
psiAvg     = 0.5*(psi_L+psi_R)
#endif

! Entropy conserving and kinetic energy conserving flux
Fstar(1) = rhoLN*vAvg(1)
Fstar(2) = Fstar(1)*vAvg(1) - smu_0*BAvg(1)*BAvg(1) + pTilde
Fstar(3) = Fstar(1)*vAvg(2) - smu_0*BAvg(1)*BAvg(2)
Fstar(4) = Fstar(1)*vAvg(3) - smu_0*BAvg(1)*BAvg(3)
Fstar(7) = vAvg(1)*Bavg(2) - BAvg(1)*vAvg(2)
Fstar(8) = vAvg(1)*Bavg(3) - BAvg(1)*vAvg(3)
#ifdef PP_GLM
Fstar(6) = GLM_ch*psiAvg
Fstar(9) = GLM_ch*BAvg(1)
#else
Fstar(6) =0.
#endif

Fstar(5) = Fstar(1)*0.5*(skappaM1/betaLN - 0.5*(v2_L+v2_R))  &
           + SUM(vAvg(:)*Fstar(2:4)) &
           +smu_0*( SUM(BAvg(:)*Fstar(6:8)) &
                   -0.5*v1_B2Avg +BAvg(1)*vB_Avg &
#ifdef PP_GLM
                   +Fstar(9)*psiAvg-GLM_ch*0.5*(psi_L*B_L(1)+psi_R*B_R(1))    &
#endif
                   )

END ASSOCIATE 
END SUBROUTINE EntropyAndKinEnergyConservingFlux


!==================================================================================================================================
!> entropy conservation for MHD, kinetric Energy conservation only in the Euler case
!> following D.Dergs et al."a novel Entropy consistent nine-wave field divergence diminishing ideal MHD system" 
!> mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t) 
!> firectly compute tranformed flux: fstar=f*metric1+g*metric2+h*metric3
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
PURE SUBROUTINE EntropyAndKinEnergyConservingFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
USE MOD_Equation_Vars,ONLY:sKappaM1,s2mu_0,smu_0
#ifdef PP_GLM
USE MOD_Equation_Vars,ONLY:GLM_ch
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL             !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR             !< right state
REAL,DIMENSION(nAuxVar),INTENT(IN)  :: UauxL          !< left auxiliary variables
REAL,DIMENSION(nAuxVar),INTENT(IN)  :: UauxR          !< right auxiliary variables
REAL,INTENT(IN)                     :: metric_L(3)    !< left metric
REAL,INTENT(IN)                     :: metric_R(3)    !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< transformed flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                   :: vAvg(3),BAvg(3)
REAL                   :: rhoLN,betaLN
REAL                   :: beta_R,beta_L
REAL                   :: vm,Bm,pTilde
#ifdef PP_GLM
REAL                   :: PsiAvg
#endif /*PP_GLM*/
REAL                   :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)

ASSOCIATE(  rho_L => UL(1)  ,  rho_R => UR(1)    , &
              B_L => UL(6:8),    B_R => UR(6:8)  , &
#ifdef PP_GLM
            psi_L => UL(9)  ,  psi_R => UR(9)    , &
#endif /*PP_GLM*/
              v_L =>UauxL(2:4),  v_R =>UauxR(2:4), &
             pt_L =>UauxL(5),   pt_R =>UauxR(5)  , & !pt=p+1/(2mu_0)|B|^2
             v2_L =>UauxL(6),   v2_R =>UauxR(6)  , & !|v|^2 left/right
             B2_L =>UauxL(7),   B2_R =>UauxR(7)  , & !|B|^2 left/right
             vB_L =>UauxL(8),   vB_R =>UauxR(8)  )

!p_L=pt_L -s2mu_0*B2_L
!p_R=pt_R -s2mu_0*B2_R
beta_L = 0.5*rho_L/(pt_L-s2mu_0*B2_L) !0.5*rho_L/p_L
beta_R = 0.5*rho_R/(pt_R-s2mu_0*B2_R) !0.5*rho_R/p_R 

! Get the averages for the numerical flux

!rho_MEAN  = 0.5*(   rho_L+rho_R)
!beta_MEAN = 0.5*(    beta_L+beta_R)
rhoLN     = LN_MEAN( rho_L, rho_R)
betaLN    = LN_MEAN(beta_L,beta_R)
vAvg      = 0.5*( v_L(:)+ v_R(:))
BAvg      = 0.5*( B_L(:)+ B_R(:))
!B2Avg     = 0.5*(B2_L+B2_R)
pTilde    = 0.5*((rho_L+rho_R)/(beta_L+beta_R)+smu_0*0.5*(B2_L+B2_R))  !rho_MEAN/(2*beta_MEAN) + 1/(2mu_0){{|B|^2}}

vm=SUM(vAvg(:)*metric(:))
Bm=SUM(BAvg(:)*metric(:))
#ifdef PP_GLM
PsiAvg = 0.5*(Psi_L+Psi_R)
#endif /*PP_GLM*/

! Entropy conserving and kinetic energy conserving flux
Fstar(  1) = rhoLN*vm
Fstar(2:4) = Fstar(1)*vAvg(1:3)-(smu_0*Bm)*BAvg(:) + pTilde*metric(1:3)
#ifdef PP_GLM
Fstar(6:8) = vm*BAvg(1:3) - Bm*vAvg(1:3) + (GLM_ch*PsiAvg)*metric(1:3)
Fstar(  9) = GLM_ch*Bm
#else
Fstar(6:8) = vm*BAvg(1:3) - Bm*vAvg(1:3)
#endif /*PP_GLM*/


Fstar(5) = Fstar(1)*0.5*(skappaM1/betaLN - 0.5*(v2_L+v2_R)) &
           + SUM(vAvg(:)*Fstar(2:4))  &
           +smu_0*( SUM(BAvg(:)*Fstar(6:8))                        &
                   - 0.25*SUM((B2_L*v_L(:)+B2_R*v_R(:))*metric(:)) & ! -0.5* {{|B|^2v(:)}}.{{m(:)}}
                   + 0.5*(vb_L+vb_R)*bm                            & !{{(v.B)}}{{ B(:) }} . {{m(:)}}
#ifdef PP_GLM
                   +Fstar(9)*PsiAvg - (GLM_ch*0.5)*SUM((Psi_L*B_L(:)+Psi_R*B_R(:))*metric(:))  & !c_h{{psi B}}.{{m}}
#endif
                  )
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE EntropyAndKinEnergyConservingFluxVec

!==================================================================================================================================
!> Computes the standard DG euler flux transformed with the metrics 
!> fstar=1/2((fL*metric1L+gL*metric2L+h*metric3L)+(fR*metric1R+gR*metric2R+h*metric3R)  )
!==================================================================================================================================
PURE SUBROUTINE EntropyAndKinEnergyConservingFluxMat(U_in,metric_in,Fstar) 
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:sKappaM1,kappaM1,kappaM2,smu_0,s2mu_0
#ifdef PP_GLM
USE MOD_Equation_vars ,ONLY:GLM_ch
#endif /*PP_GLM*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: U_in(PP_nVar,1:PP_N+1)   !< right state
REAL,INTENT(IN) :: metric_in(3,1:PP_N+1)   !< right metric
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: Fstar(PP_nVar,1:PP_N+1,1:PP_N+1)          !< transformed flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: l,i
REAL,DIMENSION(1:PP_N+1) :: srho,v1,v2,v3,pt,vv,vb,bb,beta,qvHat
#ifdef PP_GLM
REAL,DIMENSION(1:PP_N+1) :: Etotal
#endif 
REAL                   :: qBHat,v1Avg,v2Avg,v3Avg,B1Avg,B2Avg,B3Avg
REAL                   :: rhoLN,betaLN
REAL                   :: vm,Bm,pTilde
#ifdef PP_GLM
REAL                   :: PsiAvg
#endif /*PP_GLM*/
REAL                   :: metric1,metric2,metric3
#if NONCONS
REAL                   :: btilde
#endif /*NONCONS*/
!==================================================================================================================================
ASSOCIATE(  rhov1 =>U_in(2,:), &
            rhov2 =>U_in(3,:), &
            rhov3 =>U_in(4,:), &
#ifndef PP_GLM
            Etotal=>U_in(5,:), &
#endif
            b1    =>U_in(6,:), &
            b2    =>U_in(7,:), &
            b3    =>U_in(8,:)  ) 
srho(:)  = 1./U_in(1,:)
v1(:)    = rhov1*srho(:)
v2(:)    = rhov2*srho(:)
v3(:)    = rhov3*srho(:)
qvHat(:) = v1(:)*metric_in(1,:) +   v2(:)*metric_in(2,:) +   v3(:)*metric_in(3,:) 
vv(:)    = v1(:)*v1(:)+v2(:)*v2(:)+v3(:)*v3(:)
vb(:)    = v1(:)*b1(:)+v2(:)*b2(:)+v3(:)*b3(:)
bb(:)    = b1(:)*b1(:)+b2(:)*b2(:)+b3(:)*b3(:)
#ifdef PP_GLM
Etotal(:)=U_in(5,:)-s2mu_0*U_in(9,:)**2
#endif /*def PP_GLM*/
pt(:)    = kappaM1*(Etotal(:)-0.5*U_in(1,:)*vv(:)) &
           -kappaM2*s2mu_0*bb(:) ! = -kappaM1*s2mu_0*bb + s2mu_0*bb, include magn. press
beta(:) = 0.5*U_in(1,:)/(pt-s2mu_0*bb(:)) !0.5*rho/p

DO i=1,PP_N+1
  qbHat    =   b1(i)*metric_in(1,i) +   b2(i)*metric_in(2,i) +   b3(i)*metric_in(3,i) 
  !consistency euler flux f(Ui,Ui)=f(Ui)
  Fstar(1,i,i) = (U_in(1,i)*qvHat(i)  )
  Fstar(2,i,i) = ( rhov1(i)*qvHat(i)  + metric_in(1,i)*pt(i) -smu_0*(qbHat*b1(i)) )
  Fstar(3,i,i) = ( rhov2(i)*qvHat(i)  + metric_in(2,i)*pt(i) -smu_0*(qbHat*b2(i)) )
  Fstar(4,i,i) = ( rhov3(i)*qvHat(i)  + metric_in(3,i)*pt(i) -smu_0*(qbHat*b3(i)) )

#ifdef PP_GLM
  Fstar(5,i,i) = (Etotal(i) + pt(i))*qvHat(i) -smu_0*(qbHat*vb(i)) + GLM_ch*U_in(9,i)*qbHat             
  Fstar(6,i,i) = (qvHat(i)*b1(i)-qbHat*v1(i))                      + GLM_ch*U_in(9,i)*metric_in(1,i) 
  Fstar(7,i,i) = (qvHat(i)*b2(i)-qbHat*v2(i))                      + GLM_ch*U_in(9,i)*metric_in(2,i) 
  Fstar(8,i,i) = (qvHat(i)*b3(i)-qbHat*v3(i))                      + GLM_ch*U_in(9,i)*metric_in(3,i) 
  Fstar(9,i,i) =                                                  GLM_ch*qbHat                   
#else
  Fstar(5,i,i) = (Etotal(i) + pt(i))*qvHat(i) -smu_0*(qbHat*vb(i)) 
  Fstar(6,i,i) = (qvHat(i)*b1(i)-qbHat*v1(i))
  Fstar(7,i,i) = (qvHat(i)*b2(i)-qbHat*v2(i))
  Fstar(8,i,i) = (qvHat(i)*b3(i)-qbHat*v3(i))

#endif /* PP_GLM */
  DO l=i+1,PP_N+1
    rhoLN     = LN_MEAN( U_in(1,i), U_in(1,l))
    betaLN    = LN_MEAN(beta(i),beta(l))
    v1Avg     = 0.5*( v1(i)+ v1(l))
    v2Avg     = 0.5*( v2(i)+ v2(l))
    v3Avg     = 0.5*( v3(i)+ v3(l))
    B1Avg     = 0.5*( b1(i)+ b1(l))
    B2Avg     = 0.5*( b2(i)+ b2(l))
    B3Avg     = 0.5*( b3(i)+ b3(l))
    pTilde    = 0.5*((U_in(1,i)+U_in(1,l))/(beta(i)+beta(l))+smu_0*0.5*(bb(i)+bb(l)))  !rho_MEAN/(2*beta_MEAN) + 1/(2mu_0){{|B|^2}}

    metric1 = 0.5*(metric_in(1,i)+metric_in(1,l))
    metric2 = 0.5*(metric_in(2,i)+metric_in(2,l))
    metric3 = 0.5*(metric_in(3,i)+metric_in(3,l))

    vm=v1Avg*metric1+v2Avg*metric2+v3Avg*metric3
    Bm=B1Avg*metric1+B2Avg*metric2+B3Avg*metric3

    Fstar(1,l,i) = rhoLN*vm
    Fstar(2,l,i) = Fstar(1,l,i)*v1Avg-(smu_0*Bm)*B1Avg + pTilde*metric1
    Fstar(3,l,i) = Fstar(1,l,i)*v2Avg-(smu_0*Bm)*B2Avg + pTilde*metric2
    Fstar(4,l,i) = Fstar(1,l,i)*v3Avg-(smu_0*Bm)*B3Avg + pTilde*metric3
#ifdef PP_GLM
    PsiAvg = 0.5*(U_in(9,i)+U_in(9,l))
    Fstar(6,l,i) = vm*B1Avg - Bm*v1Avg + (GLM_ch*PsiAvg)*metric1
    Fstar(7,l,i) = vm*B2Avg - Bm*v2Avg + (GLM_ch*PsiAvg)*metric2
    Fstar(8,l,i) = vm*B3Avg - Bm*v3Avg + (GLM_ch*PsiAvg)*metric3
    Fstar(9,l,i) = GLM_ch*Bm
#else
    Fstar(6,l,i) = vm*B1Avg - Bm*v1Avg
    Fstar(7,l,i) = vm*B2Avg - Bm*v2Avg
    Fstar(8,l,i) = vm*B3Avg - Bm*v3Avg
#endif /*PP_GLM*/
    Fstar(5,l,i) = Fstar(1,l,i)*0.5*(skappaM1/betaLN - 0.5*(vv(i)+vv(l))) &
                   +v1Avg*Fstar(2,l,i)+v2Avg*Fstar(3,l,i)+v3Avg*Fstar(4,l,i)  &
                   +smu_0*( B1Avg*Fstar(6,l,i)+B2Avg*Fstar(7,l,i)+B3Avg*Fstar(8,l,i)    &
                           - 0.25*((bb(i)*v1(i)+bb(l)*v1(l))*metric1+   &
                                   (bb(i)*v2(i)+bb(l)*v2(l))*metric2+   &
                                   (bb(i)*v3(i)+bb(l)*v3(l))*metric3 )  &! -0.5* {{|B|^2v(:)}}.{{m(:)}}
                           + 0.5*(vb(i)+vb(l))*bm                       & !{{(v.B)}}{{ B(:) }} . {{m(:)}}
#ifdef PP_GLM
                           +Fstar(9,l,i)*PsiAvg  &
                           - (GLM_ch*0.5)*((U_in(9,i)*b1(i)+U_in(9,l)*b1(l))*metric1+ &
                                           (U_in(9,i)*b2(i)+U_in(9,l)*b2(l))*metric2+ &
                                           (U_in(9,i)*b3(i)+U_in(9,l)*b3(l))*metric3) &!c_h{{psi B}}.{{m}}
#endif
                          )
    !symmetry
    Fstar(:,i,l)=Fstar(:,l,i)
  END DO !l=i+1,PP_N+1
END DO !i=1,PP_N+1
!add nonconservative (non-symmetric!) term
#if NONCONS
#ifdef PP_GLM
qvHat(:) =   0.5*qvHat(:)
#endif /*PP_GLM*/
DO i=1,PP_N+1
  DO l=1,PP_N+1
    btilde=0.25*(  (metric_in(1,i)+metric_in(1,l))*b1(i)  &
                  +(metric_in(2,i)+metric_in(2,l))*b2(i)  &
                  +(metric_in(3,i)+metric_in(3,l))*b3(i)  )

    !powell term
    Fstar(2,l,i) = Fstar(2,l,i)+ btilde*b1(i)
    Fstar(3,l,i) = Fstar(3,l,i)+ btilde*b2(i)
    Fstar(4,l,i) = Fstar(4,l,i)+ btilde*b3(i)
    Fstar(6,l,i) = Fstar(6,l,i)+ btilde*v1(i)
    Fstar(7,l,i) = Fstar(7,l,i)+ btilde*v2(i)
    Fstar(8,l,i) = Fstar(8,l,i)+ btilde*v3(i)
#ifdef PP_GLM
    !nonconservative term to restore galilein invariance for GLM term, 1/2 cancels with 2*Dmat 
    ! grad\psi (0,0,0,0,vec{v}\phi, 0,0,0, \vec{v}) => vec{Ja^d}_{i,j,k} . vec{v}_ijk \psi_l,j,k
    Fstar(5,l,i) = Fstar(5,l,i)+ btilde*vb(i) + U_in(9,l)*qvHat(i)*U_in(9,i)
    Fstar(9,l,i) = Fstar(9,l,i)               + U_in(9,l)*qvHat(i)
#else
    Fstar(5,l,i) = Fstar(5,l,i)+ btilde*vb(i)
#endif /*PP_GLM*/
  END DO !l=1,PP_N+1
END DO !i=1,PP_N+1
#endif /*NONCONS*/

END ASSOCIATE

END SUBROUTINE EntropyAndKinEnergyConservingFluxMat

!==================================================================================================================================
!> Computes the logarithmic mean: (aR-aL)/(LOG(aR)-LOG(aL)) = (aR-aL)/LOG(aR/aL)
!> Problem: if aL~= aR, then 0/0, but should tend to --> 0.5*(aR+aL)
!>
!> introduce xi=aR/aL and f=(aR-aL)/(aR+aL) = (xi-1)/(xi+1) 
!> => xi=(1+f)/(1-f) 
!> => Log(xi) = log(1+f)-log(1-f), and for smaRl f (f^2<1.0E-02) :
!>
!>    Log(xi) ~=     (f - 1/2 f^2 + 1/3 f^3 - 1/4 f^4 + 1/5 f^5 - 1/6 f^6 + 1/7 f^7)
!>                  +(f + 1/2 f^2 + 1/3 f^3 + 1/4 f^4 + 1/5 f^5 + 1/6 f^6 + 1/7 f^7)
!>             = 2*f*(1           + 1/3 f^2           + 1/5 f^4           + 1/7 f^6)
!>  (aR-aL)/Log(xi) = (aR+aL)*f/(2*f*(1 + 1/3 f^2 + 1/5 f^4 + 1/7 f^6)) = (aR+aL)/(2 + 2/3 f^2 + 2/5 f^4 + 2/7 f^6)
!>  (aR-aL)/Log(xi) = 0.5*(aR+aL)*(105/ (105+35 f^2+ 21 f^4 + 15 f^6)
!==================================================================================================================================
PURE FUNCTION LN_MEAN(aL,aR)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: aL  !< left value
REAL,INTENT(IN) :: aR  !< right value
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL            :: LN_MEAN  !< result
!----------------------------------------------------------------------------------------------------------------------------------
! LOCaR VaLIABLES
REAL           :: Xi,u
REAL,PARAMETER :: eps=1.0E-4  ! tolerance for f^2, such that switch is smooth in double precision 
!==================================================================================================================================
Xi = aR/aL
u=(Xi*(Xi-2.)+1.)/(Xi*(Xi+2.)+1.) !u=f^2, f=(aR-aL)/(aR+aL)=(xi-1)/(xi+1)
LN_MEAN=MERGE((aL+aR)*52.5d0/(105.d0 + u*(35.d0 + u*(21.d0 +u*15.d0))), & !u <eps (test true)
              (aR-aL)/LOG(Xi)                                         , & !u>=eps (test false)
              (u.LT.eps)                                              )   !test
END FUNCTION LN_MEAN


END MODULE MOD_Flux_Average
