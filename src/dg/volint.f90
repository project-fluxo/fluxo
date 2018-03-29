!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
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
INTERFACE VolInt
#if (PP_DiscType==1)
  MODULE PROCEDURE VolInt_weakForm
#elif (PP_DiscType==2)
  MODULE PROCEDURE VolInt_SplitForm
  !new optimized version, other versions are still here, see below
  !MODULE PROCEDURE VolInt_SplitForm3 
#endif /*PP_DiscType*/
END INTERFACE

PUBLIC::VolInt
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
USE MOD_DG_Vars   ,ONLY:D_hat
USE MOD_Flux      ,ONLY:EvalFluxTilde3D  
USE MOD_Mesh_Vars ,ONLY:nElems
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
  ! Compute for all Gauss point values already the tranformed fluxes
  CALL EvalFluxTilde3D(iElem,ftilde,gtilde,htilde)
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
  DO l=0,PP_N
    DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + D_hat(i,l)*ftilde(:,l,j,k) + &
                                              D_hat(j,l)*gtilde(:,i,l,k) + &
                                              D_hat(k,l)*htilde(:,i,j,l)
    END DO; END DO; END DO ! i,j,k
  END DO ! l
END DO ! iElem
!PRINT*,'inside volint: ftilde', ftilde
!PRINT*,'inside volint: gtilde', gtilde
!PRINT*,'inside volint: htilde', htilde
END SUBROUTINE VolInt_weakForm
#endif /*PP_DiscType==1*/



#if (PP_DiscType==2)
!==================================================================================================================================
!> Computes the volume integral using flux differencing 
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut=0. and is updated with the volume flux derivatives
!==================================================================================================================================
SUBROUTINE VolInt_SplitForm(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars   ,ONLY:DvolSurf
USE MOD_Mesh_Vars ,ONLY:nElems
#if PARABOLIC
USE MOD_Flux      ,ONLY:EvalDiffFluxTilde3D
USE MOD_DG_Vars   ,ONLY:D_Hat
#endif /*PARABOLIC*/
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
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N)      :: ftildeDiff,gtildeDiff,htildeDiff !transformed diffusion fluxes
#endif /*PARABOLIC*/
INTEGER                                           :: i,j,k,l,iElem
!==================================================================================================================================

DO iElem=1,nElems
  !compute Euler contribution of the fluxes, 
  CALL EvalEulerFluxAverage3D(iElem,ftilde,gtilde,htilde)
#if PARABOLIC
  !compute Diffusion flux contribution of 
  CALL EvalDiffFluxTilde3D(iElem,ftildeDiff,gtildeDiff,htildeDiff)
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
  ! euler fluxes in flux differencing form: strong, but surface parts included 
  ! diffusion fluxes are accouted in the standard weak form
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(i,l)*  ftilde(:,l,i,j,k)  &
                                            + Dvolsurf(j,l)*  gtilde(:,l,i,j,k)  &
                                            + Dvolsurf(k,l)*  htilde(:,l,i,j,k)  &
                                            +    D_Hat(i,l)*ftildeDiff(:,l,j,k)  &
                                            +    D_Hat(j,l)*gtildeDiff(:,i,l,k)  &
                                            +    D_Hat(k,l)*htildeDiff(:,i,j,l)
    END DO ! l
  END DO; END DO; END DO ! i,j,k
#else  
  !only euler
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(i,l)*ftilde(:,l,i,j,k)  &
                                            + Dvolsurf(j,l)*gtilde(:,l,i,j,k)  &
                                            + Dvolsurf(k,l)*htilde(:,l,i,j,k)
    END DO ! l
  END DO; END DO; END DO ! i,j,k
#endif /*PARABOLIC*/
END DO ! iElem
END SUBROUTINE VolInt_SplitForm


!==================================================================================================================================
!> Compute flux differences in 3D, making use of the symmetry and appling also directly the metrics  
!==================================================================================================================================
SUBROUTINE EvalEulerFluxAverage3D(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars        ,ONLY:U
USE MOD_Mesh_Vars      ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars  ,ONLY:VolumeFluxAverageVec !pointer to flux averaging routine
USE MOD_Equation_Vars  ,ONLY:nAuxVar
USE MOD_Flux_Average   ,ONLY:EvalEulerFluxTilde3D
#if NONCONS
USE MOD_Flux_Average   ,ONLY:AddNonConsFluxTilde3D
#endif /*NONCONS*/
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
! due to consisteny, if left and right are the same, its just the transformed eulerflux
CALL EvalEulerFluxTilde3D(iElem,ftilde_c,gtilde_c,htilde_c,Uaux)

DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !symmetric part
  ftilde(:,i,i,j,k)=ftilde_c(:,i,j,k) 
  gtilde(:,j,i,j,k)=gtilde_c(:,i,j,k) 
  htilde(:,k,i,j,k)=htilde_c(:,i,j,k) 
  DO l=i+1,PP_N
    CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,l,j,k,iElem), &
                                           Uaux(:,i,j,k)      ,           Uaux(:,l,j,k)      , &
                                 Metrics_fTilde(:,i,j,k,iElem), Metrics_fTilde(:,l,j,k,iElem), &
                                       ftilde(:,l,i,j,k)                                       )
    ftilde(:,i,l,j,k)=ftilde(:,l,i,j,k) !symmetric
  END DO!l=i+1,N
  DO l=j+1,PP_N
    CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,i,l,k,iElem), &
                                           Uaux(:,i,j,k)      ,           Uaux(:,i,l,k)      , &
                                 Metrics_gTilde(:,i,j,k,iElem), Metrics_gTilde(:,i,l,k,iElem), &
                                       gtilde(:,l,i,j,k)                                       )
    gtilde(:,j,i,l,k)=gtilde(:,l,i,j,k) !symmetric
  END DO!l=j+1,N
  DO l=k+1,PP_N
    CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,i,j,l,iElem), &
                                           Uaux(:,i,j,k)      ,           Uaux(:,i,j,l)      , &
                                 Metrics_hTilde(:,i,j,k,iElem), Metrics_hTilde(:,i,j,l,iElem), &
                                       htilde(:,l,i,j,k)                                       )
    htilde(:,k,i,j,l)=htilde(:,l,i,j,k) !symmetric
  END DO!l=k+1,N
END DO; END DO; END DO ! i,j,k

#if NONCONS
CALL AddNonConsFluxTilde3D(iElem,Uaux,ftilde,gtilde,htilde)
#endif /*NONCONS*/

END SUBROUTINE EvalEulerFluxAverage3D


!==================================================================================================================================
!> Computes the volume integral using flux differencing 
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut=0. and is updated with the volume flux derivatives
!>
!> VERSION 2: low memory consumption, do not store 4D arrays ftilde,gtilde,htilde and update Ut directly,
!>            stupid variant without taking symmetry and diagonal into account
!>            => 2x slower than version 1 or version 3
!==================================================================================================================================
SUBROUTINE VolInt_SplitForm2(Ut)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars            ,ONLY:DvolSurf
USE MOD_DG_Vars            ,ONLY:U
USE MOD_Mesh_Vars          ,ONLY:nElems
USE MOD_Mesh_Vars          ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars      ,ONLY:VolumeFluxAverageVec !pointer to flux Averaging routine
USE MOD_Equation_Vars      ,ONLY:nAuxVar
USE MOD_Flux_Average       ,ONLY:EvalUaux
#if PARABOLIC
USE MOD_Flux               ,ONLY:EvalDiffFluxTilde3D
USE MOD_DG_Vars            ,ONLY:D_Hat
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)                                :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N):: ftilde,gtilde,htilde !transformed flux differences, one more dimension!
                                                                           ! ftilde(:,l,i,j,k) ={{metrics1}}.vecF(U_ljk,U_ijk)
                                                                           ! gtilde(:,l,i,j,k) ={{metrics2}}.vecF(U_ilk,U_ijk)
                                                                           ! htilde(:,l,i,j,k) ={{metrics3}}.vecF(U_ijl,U_ijk)
REAL,DIMENSION(PP_nVar)                           :: ftilde,gtilde,htilde
REAL,DIMENSION(nAuxVar,0:PP_N,0:PP_N,0:PP_N)      :: Uaux
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N)      :: ftildeDiff,gtildeDiff,htildeDiff !transformed diffusion fluxes
#endif /*PARABOLIC*/
INTEGER                                           :: i,j,k,l,iElem
!==================================================================================================================================

DO iElem=1,nElems

  !compute Euler contribution of the fluxes, 
  !compute Diffusion flux contribution of 
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
  ! euler fluxes in flux differencing form: strong, but surface parts included 
  ! diffusion fluxes are accouted in the standard weak form
  CALL EvalUaux(iElem,Uaux)
#if PARABOLIC
  CALL EvalDiffFluxTilde3D(iElem,ftildeDiff,gtildeDiff,htildeDiff)
#endif /*PARABOLIC*/

  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,l,j,k,iElem), &
                                             Uaux(:,i,j,k)      ,           Uaux(:,l,j,k)      , &
                                   Metrics_fTilde(:,i,j,k,iElem), Metrics_fTilde(:,l,j,k,iElem), &
                                           ftilde(:)                                             )

      CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,i,l,k,iElem), &
                                             Uaux(:,i,j,k)      ,           Uaux(:,i,l,k)      , &
                                   Metrics_gTilde(:,i,j,k,iElem), Metrics_gTilde(:,i,l,k,iElem), &
                                           gtilde(:)                                             )

      CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,i,j,l,iElem), &
                                             Uaux(:,i,j,k)      ,           Uaux(:,i,j,l)      , &
                                   Metrics_hTilde(:,i,j,k,iElem), Metrics_hTilde(:,i,j,l,iElem), &
                                           htilde(:)                                             )

#if PARABOLIC
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(i,l)*  ftilde(:)          &
                                            + Dvolsurf(j,l)*  gtilde(:)          &
                                            + Dvolsurf(k,l)*  htilde(:)          &
                                            +    D_Hat(i,l)*ftildeDiff(:,l,j,k)  &
                                            +    D_Hat(j,l)*gtildeDiff(:,i,l,k)  &
                                            +    D_Hat(k,l)*htildeDiff(:,i,j,l)
#else
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(i,l)*  ftilde(:)          &
                                            + Dvolsurf(j,l)*  gtilde(:)          &
                                            + Dvolsurf(k,l)*  htilde(:) 
#endif /*PARABOLIC*/
    END DO ! l=0,N
  END DO; END DO; END DO ! i,j,k
END DO ! iElem
END SUBROUTINE VolInt_SplitForm2


!==================================================================================================================================
!> Computes the volume integral using flux differencing 
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut=0. and is updated with the volume flux derivatives
!>
!> VERSION 3: Lower memory consumption, do not store 4D arrays ftilde,gtilde,htilde and update Ut directly,
!>            clever variant taking symmetry and diagonal into account (like version 1, also performance-wise) 
!==================================================================================================================================
SUBROUTINE VolInt_SplitForm3(Ut)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars            ,ONLY:DvolSurf
USE MOD_DG_Vars            ,ONLY:U
USE MOD_Mesh_Vars          ,ONLY:nElems
USE MOD_Mesh_Vars          ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Equation_Vars      ,ONLY:VolumeFluxAverageVec !pointer to flux Averaging routine
USE MOD_Equation_Vars      ,ONLY:nAuxVar
USE MOD_Flux_Average       ,ONLY:EvalEulerFluxTilde3D
#if PARABOLIC
USE MOD_Flux               ,ONLY:EvalDiffFluxTilde3D
USE MOD_DG_Vars            ,ONLY:D_Hat
#endif /*PARABOLIC*/
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
!< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N):: ftilde,gtilde,htilde !transformed flux differences, one more dimension!
                                                                           ! ftilde(:,l,i,j,k) ={{metrics1}}.vecF(U_ljk,U_ijk)
                                                                           ! gtilde(:,l,i,j,k) ={{metrics2}}.vecF(U_ilk,U_ijk)
                                                                           ! htilde(:,l,i,j,k) ={{metrics3}}.vecF(U_ijl,U_ijk)
REAL,DIMENSION(PP_nVar)                           :: ftilde,gtilde,htilde
REAL,DIMENSION(nAuxVar,0:PP_N,0:PP_N,0:PP_N)      :: Uaux
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N)      :: ftilde_c,gtilde_c,htilde_c !transformed diffusion fluxes
#if PARABOLIC
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N)      :: ftildeDiff,gtildeDiff,htildeDiff !transformed diffusion fluxes
#endif /*PARABOLIC*/
INTEGER                                           :: i,j,k,l,iElem
!==================================================================================================================================

DO iElem=1,nElems

  !compute Euler contribution of the fluxes, 
  !compute Diffusion flux contribution of 
  ! Update the time derivative with the spatial derivatives of the transformed fluxes
  ! euler fluxes in flux differencing form: strong, but surface parts included 
  ! diffusion fluxes are accouted in the standard weak form
  CALL EvalEulerFluxTilde3D(iElem,ftilde_c,gtilde_c,htilde_c,Uaux)
#if PARABOLIC
  CALL EvalDiffFluxTilde3D(iElem,ftildeDiff,gtildeDiff,htildeDiff)
#endif /*PARABOLIC*/

  DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
    DO l=i+1,PP_N
      CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,l,j,k,iElem), &
                                             Uaux(:,i,j,k)      ,           Uaux(:,l,j,k)      , &
                                   Metrics_fTilde(:,i,j,k,iElem), Metrics_fTilde(:,l,j,k,iElem), &
                                         ftilde(:)                                               )
#if PARABOLIC
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(i,l)*  ftilde(:) +  D_Hat(i,l)*ftildeDiff(:,l,j,k)
      Ut(:,l,j,k,iElem) = Ut(:,l,j,k,iElem) + Dvolsurf(l,i)*  ftilde(:) +  D_Hat(l,i)*ftildeDiff(:,i,j,k)
#else
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(i,l)*  ftilde(:)
      Ut(:,l,j,k,iElem) = Ut(:,l,j,k,iElem) + Dvolsurf(l,i)*  ftilde(:)
#endif /*PARABOLIC*/
    END DO !l=i+1,N
    DO l=j+1,PP_N
      CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,i,l,k,iElem), &
                                             Uaux(:,i,j,k)      ,           Uaux(:,i,l,k)      , &
                                   Metrics_gTilde(:,i,j,k,iElem), Metrics_gTilde(:,i,l,k,iElem), &
                                         gtilde(:)                                               )
#if PARABOLIC
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(j,l)*  gtilde(:) +  D_Hat(j,l)*gtildeDiff(:,i,l,k) 
      Ut(:,i,l,k,iElem) = Ut(:,i,l,k,iElem) + Dvolsurf(l,j)*  gtilde(:) +  D_Hat(l,j)*gtildeDiff(:,i,j,k) 
#else
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(j,l)*  gtilde(:) 
      Ut(:,i,l,k,iElem) = Ut(:,i,l,k,iElem) + Dvolsurf(l,j)*  gtilde(:) 
#endif /*PARABOLIC*/
    END DO !l=j+1,N
    DO l=k+1,PP_N
      CALL VolumeFluxAverageVec(                U(:,i,j,k,iElem),              U(:,i,j,l,iElem), &
                                             Uaux(:,i,j,k)      ,           Uaux(:,i,j,l)      , &
                                   Metrics_hTilde(:,i,j,k,iElem), Metrics_hTilde(:,i,j,l,iElem), &
                                         htilde(:)                                               )
#if PARABOLIC
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(k,l)*  htilde(:) +  D_Hat(k,l)*htildeDiff(:,i,j,l) 
      Ut(:,i,j,l,iElem) = Ut(:,i,j,l,iElem) + Dvolsurf(l,k)*  htilde(:) +  D_Hat(l,k)*htildeDiff(:,i,j,k) 
#else
      Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(k,l)*  htilde(:)
      Ut(:,i,j,l,iElem) = Ut(:,i,j,l,iElem) + Dvolsurf(l,k)*  htilde(:) 
#endif /*PARABOLIC*/
    END DO ! l=k+1,N
#if PARABOLIC
    Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(i,i)*  ftilde_c(:,i,j,k)  &
                                          + Dvolsurf(j,j)*  gtilde_c(:,i,j,k)  &
                                          + Dvolsurf(k,k)*  htilde_c(:,i,j,k)  &
                                          +    D_Hat(i,i)*ftildeDiff(:,i,j,k)  &
                                          +    D_Hat(j,j)*gtildeDiff(:,i,j,k)  &
                                          +    D_Hat(k,k)*htildeDiff(:,i,j,k)
#else
    Ut(:,i,j,k,iElem) = Ut(:,i,j,k,iElem) + Dvolsurf(i,i)*  ftilde_c(:,i,j,k)  &
                                          + Dvolsurf(j,j)*  gtilde_c(:,i,j,k)  &
                                          + Dvolsurf(k,k)*  htilde_c(:,i,j,k) 
#endif /*PARABOLIC*/
  END DO; END DO; END DO ! i,j,k
END DO ! iElem
END SUBROUTINE VolInt_SplitForm3

#endif /*PP_DiscType==2*/

END MODULE MOD_VolInt
