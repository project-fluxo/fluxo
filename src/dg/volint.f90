!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2020 - 2021 Andrés Rueda
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

MODULE MOD_VolInt
!==================================================================================================================================
! Containes the different DG volume integrals
! Computes the volume integral contribution based on U and updates Ut
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
! Public Part ---------------------------------------------------------------------------------------------------------------------
#if (PP_DiscType==1)
INTERFACE VolInt
  MODULE PROCEDURE VolInt_weakForm
END INTERFACE

PUBLIC::VolInt
PUBLIC::VolInt_adv
#elif (PP_DiscType==2)
INTERFACE VolInt_Adv_SplitForm
  !optimized VolInt with general 4D flux array from flux_average.f90, only advection part!
  MODULE PROCEDURE VolInt_SplitForm
END INTERFACE

PUBLIC::VolInt_Adv_SplitForm
#endif /*PP_DiscType*/

#if PARABOLIC
INTERFACE VolInt_visc
  MODULE PROCEDURE VolInt_visc
END INTERFACE

PUBLIC::VolInt_visc
#endif /*PARABOLIC*/

!==================================================================================================================================


CONTAINS

#if (PP_DiscType==1)
!==================================================================================================================================
!>  Computes the volume integral of the weak DG form a la Kopriva
!>  Attention 1: 1/J(i,j,k) is not yet accounted for
!>  Attention 2: input Ut=0. and is updated with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolInt_weakForm(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars   ,ONLY:D_hat_T,U
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY:gradPx,gradPy,gradPz
#endif /*PARABOLIC*/
USE MOD_Flux      ,ONLY:EvalFluxTilde3D  
USE MOD_Mesh_Vars ,ONLY:nElems,metrics_ftilde,metrics_gtilde,metrics_htilde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)                                :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N)      :: ftilde,gtilde,htilde ! transformed volume fluxes at all Gauss points
INTEGER                                           :: i,j,k,l,iElem
!==================================================================================================================================
DO iElem=1,nElems
  ! Compute for all Gauss point values already the transformed fluxes
  CALL EvalFluxTilde3D(             U(:,:,:,:,iElem), &
                       metrics_fTilde(:,:,:,:,iElem), &
                       metrics_gTilde(:,:,:,:,iElem), &
                       metrics_hTilde(:,:,:,:,iElem), &
#if PARABOLIC
                               gradPx(:,:,:,:,iElem), &
                               gradPy(:,:,:,:,iElem), &
                               gradPz(:,:,:,:,iElem), &
#endif /*PARABOLIC*/
                               ftilde,gtilde,htilde)
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
  DO l=0,PP_N
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat_T(l,i)*ftilde(:,l,j,k) + &
                                              D_hat_T(l,j)*gtilde(:,i,l,k) + &
                                              D_hat_T(l,k)*htilde(:,i,j,l)
    END DO; END DO; END DO ! i,j,k
  END DO ! l
END DO ! iElem
END SUBROUTINE VolInt_weakForm

!==================================================================================================================================
!>  Computes the volume integral of ONLY the advective terms a la Kopriva
!>  Attention 1: 1/J(i,j,k) is not yet accounted for
!>  Attention 2: input Ut=0. and is updated with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolInt_adv(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars   ,ONLY:D_hat_T,U
#if LOCAL_ALPHA
use MOD_NFVSE_Vars,only: ftilde_DG, gtilde_DG, htilde_DG
use MOD_Interpolation_Vars , only: wGP
USE MOD_DG_Vars   ,ONLY:D_T
#endif /*LOCAL_ALPHA*/
USE MOD_Flux      ,ONLY:EvalAdvFluxTilde3D  
USE MOD_Mesh_Vars ,ONLY:nElems,metrics_ftilde,metrics_gtilde,metrics_htilde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)                                :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N)      :: ftilde,gtilde,htilde ! transformed volume fluxes at all Gauss points
INTEGER                                           :: i,j,k,l,iElem
!==================================================================================================================================
DO iElem=1,nElems
  ! Compute for all Gauss point values already the transformed fluxes
  CALL EvalAdvFluxTilde3D(          U(:,:,:,:,iElem), &
                       metrics_fTilde(:,:,:,:,iElem), &
                       metrics_gTilde(:,:,:,:,iElem), &
                       metrics_hTilde(:,:,:,:,iElem), &
                               ftilde,gtilde,htilde)
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
#if LOCAL_ALPHA
    
    do i=0, PP_N-1
      ftilde_DG(:,i,:,:,iElem) = 0.
      gtilde_DG(:,:,i,:,iElem) = 0.
      htilde_DG(:,:,:,i,iElem) = 0.
      do l=0, i
        do j=0, PP_N
          ftilde_DG(:,i,:,:,iElem) = ftilde_DG(:,i,:,:,iElem) - wGP(j)*D_T(l,j)*ftilde(:,j,:,:)
          gtilde_DG(:,:,i,:,iElem) = gtilde_DG(:,:,i,:,iElem) - wGP(j)*D_T(l,j)*gtilde(:,:,j,:)
          htilde_DG(:,:,:,i,iElem) = htilde_DG(:,:,:,i,iElem) - wGP(j)*D_T(l,j)*htilde(:,:,:,j)
        end do
      end do
    end do
#endif /*LOCAL_ALPHA*/
  
  ! xi derivatives
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N; DO l=0,PP_N
    Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat_T(l,i)*ftilde(:,l,j,k)
  END DO; END DO; END DO; END DO
  
  ! eta derivatives
  DO k=0,PP_N; DO j=0,PP_N; DO l=0,PP_N; DO i=0,PP_N
    Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat_T(l,j)*gtilde(:,i,l,k)
  END DO; END DO; END DO; END DO
  
  ! zeta derivatives
  DO k=0,PP_N; DO l=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat_T(l,k)*htilde(:,i,j,l)
  END DO; END DO; END DO; END DO
  
END DO ! iElem
END SUBROUTINE VolInt_adv
#endif /*PP_DiscType==1*/



#if (PP_DiscType==2)
!==================================================================================================================================
!> Computes the volume integral using flux differencing 
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut=0. and is updated with the volume flux derivatives
!> Attention 3: If we are using ES Gauss collocation methods, we store Uaux here, since we'll need it for the surface integral
!==================================================================================================================================
SUBROUTINE VolInt_SplitForm(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,            ONLY: DvolSurf_T,U
#if (PP_NodeType==1 & defined(PP_u_aux_exist))
USE MOD_DG_Vars,            ONLY: Uaux
#endif /*(PP_NodeType==1 & defined(PP_u_aux_exist))*/
USE MOD_Mesh_Vars,          ONLY: nElems,metrics_ftilde,metrics_gtilde,metrics_htilde
USE MOD_Flux_Average,       ONLY: EvalAdvFluxAverage3D
#if LOCAL_ALPHA
#if NONCONS
USE MOD_Flux_Average,       ONLY: EvalAdvFluxAverage3D_separate
use MOD_NFVSE_Vars,         only: ftildeR_DG, gtildeR_DG, htildeR_DG
#endif NONCONS
use MOD_NFVSE_Vars,         only: ftilde_DG, gtilde_DG, htilde_DG, sWGP
use MOD_Interpolation_Vars, only: wGP
#endif /*LOCAL_ALPHA*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)                                :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N):: ftilde,gtilde,htilde !transformed flux differences, one more dimension!
                                                                           ! ftilde(:,l,i,j,k) ={{metrics1}}.vecF(U_ljk,U_ijk)
                                                                           ! gtilde(:,l,i,j,k) ={{metrics2}}.vecF(U_ilk,U_ijk)
                                                                           ! htilde(:,l,i,j,k) ={{metrics3}}.vecF(U_ijl,U_ijk)
#if NONCONS
#if defined(PP_GLM) && defined (PP_NC_GLM)
INTEGER,PARAMETER:: nnonc = 2
#else
INTEGER,PARAMETER:: nnonc = 1
#endif /*PP_GLM and PP_NC_GLM*/
REAL,DIMENSION(            nnonc,0:PP_N,0:PP_N,0:PP_N,0:PP_N) :: f_noncons,g_noncons,h_noncons !< 4D transformed symmetric part of non-conservative terms  (iVar,i,,k)
REAL,DIMENSION(1:PP_nVar,3,nnonc,       0:PP_N,0:PP_N,0:PP_N) :: phi                           !< 4D transformed symmetric part of non-conservative terms (iVar,i,,k)
REAL :: phi_tildex(nnonc,-1:PP_N, 0:PP_N, 0:PP_N) 
REAL :: phi_tildey(nnonc, 0:PP_N,-1:PP_N, 0:PP_N) 
REAL :: phi_tildez(nnonc, 0:PP_N, 0:PP_N,-1:PP_N) 
#endif /*NONCONS*/
real :: Ut_dbg(PP_nVar) ! Debug!!
INTEGER                                           :: i,j,k,l,m,iElem,term
!==================================================================================================================================
#if LOCAL_ALPHA
ftilde_DG=0.0
htilde_DG=0.0
gtilde_DG=0.0
#if NONCONS
phi_tildex=0.0
phi_tildey=0.0
phi_tildez=0.0
#endif /*NONCONS*/
#endif /*LOCAL_ALPHA*/

DO iElem=1,nElems
  
#if LOCAL_ALPHA && NONCONS
  ! Get the conservative and nonconservative fluxes separately
  CALL EvalAdvFluxAverage3D_separate(    U(:,:,:,:,iElem), &
#if (PP_NodeType==1 & defined(PP_u_aux_exist))
                                      Uaux(:,:,:,:,iElem), &
#endif /*(PP_NodeType==1 & defined(PP_u_aux_exist))*/
                            metrics_fTilde(:,:,:,:,iElem), &
                            metrics_gTilde(:,:,:,:,iElem), &
                            metrics_hTilde(:,:,:,:,iElem), &
                                   ftilde,gtilde,htilde,f_noncons,g_noncons,h_noncons,phi)
#else /*LOCAL_ALPHA && NONCONS*/
  !compute advective+nonconservative contribution of the fluxes
  CALL EvalAdvFluxAverage3D(             U(:,:,:,:,iElem), &
#if (PP_NodeType==1 & defined(PP_u_aux_exist))
                                      Uaux(:,:,:,:,iElem), &
#endif /*(PP_NodeType==1 & defined(PP_u_aux_exist))*/
                            metrics_fTilde(:,:,:,:,iElem), &
                            metrics_gTilde(:,:,:,:,iElem), &
                            metrics_hTilde(:,:,:,:,iElem), &
                                   ftilde,gtilde,htilde)
#endif /*LOCAL_ALPHA && NONCONS*/

#if LOCAL_ALPHA
    ! Conservative part of high-order "staggered" fluxes:
    ! ***************************************************
    do j=0, PP_N-1
      ftilde_DG(:,j,:,:,iElem) = ftilde_DG(:,j-1,:,:,iElem)
      gtilde_DG(:,:,j,:,iElem) = gtilde_DG(:,:,j-1,:,iElem)
      htilde_DG(:,:,:,j,iElem) = htilde_DG(:,:,:,j-1,iElem)
#if NONCONS
      phi_tildex(:,j,:,:) = phi_tildex(:,j-1,:,:)
      phi_tildey(:,:,j,:) = phi_tildey(:,:,j-1,:)
      phi_tildez(:,:,:,j) = phi_tildez(:,:,:,j-1)
#endif /*NONCONS*/
      
      do i=0, PP_N
        ftilde_DG(:,j,:,:,iElem) = ftilde_DG(:,j,:,:,iElem) + wGP(j)*Dvolsurf_T(i,j)*ftilde(:,i,j,:,:)
        gtilde_DG(:,:,j,:,iElem) = gtilde_DG(:,:,j,:,iElem) + wGP(j)*Dvolsurf_T(i,j)*gtilde(:,i,:,j,:)
        htilde_DG(:,:,:,j,iElem) = htilde_DG(:,:,:,j,iElem) + wGP(j)*Dvolsurf_T(i,j)*htilde(:,i,:,:,j)
#if NONCONS
        phi_tildex(:,j,:,:) = phi_tildex(:,j,:,:) + wGP(j)*Dvolsurf_T(i,j)*f_noncons(:,i,j,:,:)
        phi_tildey(:,:,j,:) = phi_tildey(:,:,j,:) + wGP(j)*Dvolsurf_T(i,j)*g_noncons(:,i,:,j,:)
        phi_tildez(:,:,:,j) = phi_tildez(:,:,:,j) + wGP(j)*Dvolsurf_T(i,j)*h_noncons(:,i,:,:,j)
#endif /*NONCONS*/
      end do
    end do
#if NONCONS
    ! Copy conservative part into the right non-conservative fluxes
    ftildeR_DG(:,:,:,:,iElem) = ftilde_DG(:,:,:,:,iElem)
    gtildeR_DG(:,:,:,:,iElem) = gtilde_DG(:,:,:,:,iElem)
    htildeR_DG(:,:,:,:,iElem) = htilde_DG(:,:,:,:,iElem)
#endif /*NONCONS*/
#endif /*LOCAL_ALPHA*/
  
  !only euler
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    
#if LOCAL_ALPHA
#if NONCONS
    ! Add nonconservative part to fluxes
    do term=1,nnonc
      ! Fluxes on the left
      ftilde_DG(:,i-1,j,k,iElem) = ftilde_DG(:,i-1,j,k,iElem) + phi(:,1,term,i,j,k) * phi_tildex(term,i-1,j,k)
      gtilde_DG(:,i,j-1,k,iElem) = gtilde_DG(:,i,j-1,k,iElem) + phi(:,2,term,i,j,k) * phi_tildey(term,i,j-1,k)
      htilde_DG(:,i,j,k-1,iElem) = htilde_DG(:,i,j,k-1,iElem) + phi(:,3,term,i,j,k) * phi_tildez(term,i,j,k-1)
      ! Fluxes on the right
      ftildeR_DG(:,i,j,k,iElem) = ftildeR_DG(:,i,j,k,iElem) + phi(:,1,term,i,j,k) * phi_tildex(term,i,j,k)
      gtildeR_DG(:,i,j,k,iElem) = gtildeR_DG(:,i,j,k,iElem) + phi(:,2,term,i,j,k) * phi_tildey(term,i,j,k)
      htildeR_DG(:,i,j,k,iElem) = htildeR_DG(:,i,j,k,iElem) + phi(:,3,term,i,j,k) * phi_tildez(term,i,j,k)
    end do
    ! Use flux differencing formula to get Ut
    Ut(:,i,j,k,iElem) =  sWGP(i) * ( ftildeR_DG(:,i,j,k,iElem) - ftilde_DG(:,i-1,j  ,k  ,iElem) ) &
                       + sWGP(j) * ( gtildeR_DG(:,i,j,k,iElem) - gtilde_DG(:,i  ,j-1,k  ,iElem) ) &
                       + sWGP(k) * ( htildeR_DG(:,i,j,k,iElem) - htilde_DG(:,i  ,j  ,k-1,iElem) )
#else /*NONCONS*/
    ! Use flux differencing formula to get Ut
    Ut(:,i,j,k,iElem) =  sWGP(i) * ( ftilde_DG(:,i,j,k,iElem) - ftilde_DG(:,i-1,j  ,k  ,iElem) ) &
                       + sWGP(j) * ( gtilde_DG(:,i,j,k,iElem) - gtilde_DG(:,i  ,j-1,k  ,iElem) ) &
                       + sWGP(k) * ( htilde_DG(:,i,j,k,iElem) - htilde_DG(:,i  ,j  ,k-1,iElem) )
#endif /*NONCONS*/
#else /*LOCAL_ALPHA*/
    DO l=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf_T(l,i)*ftilde(:,l,i,j,k)  &
                                            + Dvolsurf_T(l,j)*gtilde(:,l,i,j,k)  &
                                            + Dvolsurf_T(l,k)*htilde(:,l,i,j,k)
    END DO ! l
#endif /*LOCAL_ALPHA*/
    
  END DO; END DO; END DO ! i,j,k
  
END DO ! iElem
END SUBROUTINE VolInt_SplitForm
#endif /*PP_DiscType==2*/


#if PARABOLIC
!==================================================================================================================================
!> Computes the volume integral using flux differencing 
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut=0. and is updated with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolInt_visc(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY:D_Hat_T,U
USE MOD_Lifting_Vars ,ONLY:gradPx,gradPy,gradPz
USE MOD_Flux         ,ONLY:EvalDiffFluxTilde3D
USE MOD_Mesh_Vars    ,ONLY:nElems,metrics_ftilde,metrics_gtilde,metrics_htilde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)                                :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N)      :: ftildeDiff,gtildeDiff,htildeDiff !transformed diffusion fluxes
#endif /*PARABOLIC*/
INTEGER                                           :: i,j,k,l,iElem
!==================================================================================================================================

DO iElem=1,nElems
  !compute Diffusion flux contribution of 
  CALL EvalDiffFluxTilde3D(             U(:,:,:,:,iElem), &
                           metrics_fTilde(:,:,:,:,iElem), &
                           metrics_gTilde(:,:,:,:,iElem), &
                           metrics_hTilde(:,:,:,:,iElem), &
                                   gradPx(:,:,:,:,iElem), &
                                   gradPy(:,:,:,:,iElem), &
                                   gradPz(:,:,:,:,iElem), &
                                   ftildeDiff,gtildeDiff,htildeDiff)
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
  ! euler fluxes in flux differencing form: strong, but surface parts included 
  ! diffusion fluxes are accouted in the standard weak form
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) +    D_Hat_T(l,i)*ftildeDiff(:,l,j,k)  &
                                            +    D_Hat_T(l,j)*gtildeDiff(:,i,l,k)  &
                                            +    D_Hat_T(l,k)*htildeDiff(:,i,j,l)
    END DO ! l
  END DO; END DO; END DO ! i,j,k
END DO ! iElem
END SUBROUTINE VolInt_visc
#endif /*PARABOLIC*/

END MODULE MOD_VolInt
