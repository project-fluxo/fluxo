!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2021 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
! Copyright (c) 2020 - 2021 Andrés Rueda
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

#if (PP_DiscType==2)
INTERFACE EvalAdvFluxAverage3D
  MODULE PROCEDURE EvalAdvFluxAverage3D
END INTERFACE
#endif /*PP_DiscType==2*/

#if NONCONS
INTERFACE AddNonConsFluxVec
  MODULE PROCEDURE AddNonConsFluxVec
END INTERFACE
#endif /*NONCONS*/

INTERFACE StandardDGFlux
  MODULE PROCEDURE StandardDGFlux
END INTERFACE

INTERFACE StandardDGFluxVec
  MODULE PROCEDURE StandardDGFluxVec
END INTERFACE

INTERFACE StandardDGFluxDealiasedMetricVec
  MODULE PROCEDURE StandardDGFluxDealiasedMetricVec
END INTERFACE

INTERFACE EntropyandKinEnergyConservingFlux_Derigs
  MODULE PROCEDURE EntropyandKinEnergyConservingFlux_Derigs
END INTERFACE

INTERFACE EntropyandKinEnergyConservingFluxVec_Derigs
  MODULE PROCEDURE EntropyandKinEnergyConservingFluxVec_Derigs
END INTERFACE

INTERFACE EntropyandKinEnergyConservingFlux_FloGor
  MODULE PROCEDURE EntropyandKinEnergyConservingFlux_FloGor
END INTERFACE

INTERFACE EntropyandKinEnergyConservingFluxVec_FloGor
  MODULE PROCEDURE EntropyandKinEnergyConservingFluxVec_FloGor
END INTERFACE

INTERFACE LN_MEAN 
  MODULE PROCEDURE LN_MEAN
END INTERFACE


#if (PP_DiscType==2)
PUBLIC::EvalAdvFluxAverage3D
PUBLIC::EvalAdvFluxAverage
PUBLIC::EvalUaux
#endif /*PP_DiscType==2*/
#if NONCONS
PUBLIC::AddNonConsFluxVec
#endif /*NONCONS*/
PUBLIC::StandardDGFlux
PUBLIC::StandardDGFluxVec
PUBLIC::StandardDGFluxDealiasedMetricVec
PUBLIC::EntropyandKinEnergyConservingFlux_Derigs
PUBLIC::EntropyandKinEnergyConservingFluxVec_Derigs
PUBLIC::EntropyandKinEnergyConservingFlux_FloGor
PUBLIC::EntropyandKinEnergyConservingFluxVec_FloGor
PUBLIC::LN_MEAN

!==================================================================================================================================
! local definitions for inlining / optimizing routines,depending on PP_VolFlux, DEFAULT=-1: USE POINTER defined at runtime!
#if PP_VolFlux==-1
#  define PP_VolumeFluxAverageVec VolumeFluxAverageVec
#elif PP_VolFlux==0
#  define PP_VolumeFluxAverageVec StandardDGFluxVec 
#elif PP_VolFlux==1
#  define PP_VolumeFluxAverageVec StandardDGFluxDealiasedMetricVec
#elif PP_VolFlux==10
#  define PP_VolumeFluxAverageVec EntropyandKinEnergyConservingFluxVec_Derigs
#elif PP_VolFlux==12
#  define PP_VolumeFluxAverageVec EntropyandKinEnergyConservingFluxVec_FloGor
#endif

!==================================================================================================================================

CONTAINS


#if (PP_DiscType==2)
!==================================================================================================================================
!> Compute flux differences in 3D, making use of the symmetry and appling also directly the metrics  
!==================================================================================================================================
SUBROUTINE EvalAdvFluxAverage3D(U_in,&
#if (PP_NodeType==1)
                                Uaux, &
#endif /*(PP_NodeType==1)*/
                                     M_f,M_g,M_h,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
#if PP_VolFlux==-1
USE MOD_Equation_Vars  ,ONLY:VolumeFluxAverageVec !pointer to flux averaging routine
#endif
USE MOD_Equation_Vars  ,ONLY:nAuxVar
USE MOD_DG_Vars       ,ONLY:nTotal_vol
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN ) :: U_in        !< solution
REAL,DIMENSION(1:3      ,0:PP_N,0:PP_N,0:PP_N),INTENT(IN ) :: M_f,M_g,M_h !< metrics
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde,gtilde,htilde !< 4D transformed fluxes (iVar,i,,k)
#if (PP_NodeType==1)
REAL,DIMENSION(1:nAuxVar,0:PP_N,0:PP_N,0:PP_N)       ,INTENT(OUT) :: Uaux                 !auxiliary variables
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#else
REAL           :: Uaux(nAuxVar,  0:PP_N,0:PP_N,0:PP_N)  !auxiliary variables
#endif /*(PP_NodeType==1)*/
INTEGER   :: i,j,k,l
!==================================================================================================================================

!opt_v1
CALL EvalUaux(nTotal_vol,U_in,Uaux)
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part not needed since diagonal of DvolSurfMat is zero!
  !ftilde(:,i,i,j,k)=ftilde_c(:,i,j,k) 
  ftilde(:,i,i,j,k)=0.
  DO l=i+1,PP_N
    CALL PP_VolumeFluxAverageVec(U_in(:,i,j,k),U_in(:,l,j,k), &
                                 Uaux(:,i,j,k),Uaux(:,l,j,k), &
                                  M_f(:,i,j,k), M_f(:,l,j,k), &
                             ftilde(:,l,i,j,k)                )
    ftilde(:,i,l,j,k)=ftilde(:,l,i,j,k) !symmetric
  END DO!l=i+1,N
END DO; END DO; END DO ! i,j,k
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part not needed since diagonal of DvolSurfMat is zero!
  !gtilde(:,j,i,j,k)=gtilde_c(:,i,j,k) 
  gtilde(:,j,i,j,k)=0.
  DO l=j+1,PP_N
    CALL PP_VolumeFluxAverageVec(U_in(:,i,j,k),U_in(:,i,l,k), &
                                 Uaux(:,i,j,k),Uaux(:,i,l,k), &
                                  M_g(:,i,j,k), M_g(:,i,l,k), &
                             gtilde(:,l,i,j,k)                )
    gtilde(:,j,i,l,k)=gtilde(:,l,i,j,k) !symmetric
  END DO!l=j+1,N
END DO; END DO; END DO ! i,j,k
DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
  !diagonal (consistent) part not needed since diagonal of DvolSurfMat is zero!
  !htilde(:,k,i,j,k)=htilde_c(:,i,j,k) 
  htilde(:,k,i,j,k)=0.
  DO l=k+1,PP_N
    CALL PP_VolumeFluxAverageVec(U_in(:,i,j,k),U_in(:,i,j,l), &
                                 Uaux(:,i,j,k),Uaux(:,i,j,l), &
                                  M_h(:,i,j,k), M_h(:,i,j,l), &
                             htilde(:,l,i,j,k)                )
    htilde(:,k,i,j,l)=htilde(:,l,i,j,k) !symmetric
  END DO!l=k+1,N
END DO; END DO; END DO ! i,j,k

#if NONCONS
CALL AddNonConsFluxTilde3D(U_in,Uaux,M_f,M_g,M_h,ftilde,gtilde,htilde)
#endif /*NONCONS*/

END SUBROUTINE EvalAdvFluxAverage3D

!==================================================================================================================================
!> Compute volumetric flux differences (advective and non-conservative contributions) between two points appling also directly the metrics  
!==================================================================================================================================
SUBROUTINE EvalAdvFluxAverage(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
! IMPLICIT VARIABLE HANDLING
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
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar          !< transformed central flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

! Compute advective flux
CALL PP_VolumeFluxAverageVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)

END SUBROUTINE EvalAdvFluxAverage

!==================================================================================================================================
!> computes auxiliary nodal variables (1/rho,v_1,v_2,v_3,p_t,|v|^2) 
!==================================================================================================================================
PURE SUBROUTINE EvalUaux(np,Uin,Uaux)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER                     ,INTENT(IN)  :: np !size of input/output arrays
REAL,DIMENSION(PP_nVar,1:np),INTENT(IN)  :: Uin
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(nAuxVar,1:np),INTENT(OUT) :: Uaux   !< auxiliary variables:(srho,v1,v2,v3,p_t,|v|^2,|B|^2,v*b
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i 
REAL                :: srho,vel(1:3),v2,B2
!==================================================================================================================================
DO i=1,np
  ! auxiliary variables
  srho = 1./Uin(IRHO1,i) 
  vel  = Uin(IRHOU:IRHOW,i)*srho
  v2   = SUM(vel*vel)
  B2   = SUM(Uin(IB1:IB3,i)*Uin(IB1:IB3,i))
  Uaux(ISRHO,i) = srho
  Uaux(IU:IW,i) = vel
  Uaux(IVV  ,i) = v2
  Uaux(IBB  ,i) = B2
  Uaux(IVB  ,i)  =SUM(vel(:)*Uin(IB1:IB3,i)) ! v*B
  !total pressure=gas pressure + magnetic pressure
  Uaux(IP   ,i)=kappaM1*(Uin(IRHOE,i) -0.5*Uin(IRHO1,i)*v2 &
#ifdef PP_GLM
                                   -s2mu_0*Uin(IPSI ,i)**2 &
#endif /*PP_GLM*/
                                   )-kappaM2*s2mu_0*B2 !p_t 
END DO ! i
END SUBROUTINE EvalUaux
#endif /*PP_DiscType==2*/



#if NONCONS
!==================================================================================================================================
!> Compute transformed nonconservative MHD fluxes 
!==================================================================================================================================
PURE SUBROUTINE AddNonConsFluxTilde3D(U_in,Uaux,M_f,M_g,M_h,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN) :: U_in   !< solution
REAL,DIMENSION(nAuxVar,0:PP_N,0:PP_N,0:PP_N),INTENT(IN) :: Uaux   !< auxiliary variables
REAL,DIMENSION(1:3    ,0:PP_N,0:PP_N,0:PP_N),INTENT(IN) :: M_f,M_g,M_h !< metrics
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: f !< add to transformed flux f(iVar,l,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: g !< add to transformed flux g(iVar,l,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,0:PP_N),INTENT(INOUT) :: h !< add to transformed flux h(iVar,l,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i,j,k,l
#if NONCONS==1 /*Powell*/
INTEGER,PARAMETER:: vs=IRHOU
INTEGER,PARAMETER:: ve=IB3
#elif NONCONS==2 /*Brackbill*/
INTEGER,PARAMETER:: vs=IRHOU
INTEGER,PARAMETER:: ve=IRHOW
#elif NONCONS==3 /*Janhunen*/
INTEGER,PARAMETER:: vs=IB1
INTEGER,PARAMETER:: ve=IB3
#endif /*NONCONSTYPE*/
REAL :: Phi_MHD_s2(vs:ve)
REAL :: Phi_GLM_s2(2)
REAL :: Bhat_L
!==================================================================================================================================

DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N

#if NONCONS==1 /*Powell*/
  ! Powell: Phi(IRHOU:IB3) =B,vB,v
  Phi_MHD_s2(vs:ve) = 0.5 * (/U_in(IB1:IB3,i,j,k),Uaux(IVB,i,j,k),Uaux(IU:IW,i,j,k)/)
#elif NONCONS==2 /*Brackbill*/
  ! Brackbill: Phi(IRHOU:IRHOW) =B
  Phi_MHD_s2(vs:ve) = 0.5 * U_in(IB1:IB3,i,j,k)
#elif NONCONS==3 /*Janhunen*/
  ! Janhunen: Phi(IB1:IB3) =v
  Phi_MHD_s2(vs:ve) = 0.5 * Uaux(IU:IW,i,j,k)
#endif /*NONCONSTYPE*/
  
  ! Non-conservative terms in xi
  !-----------------------------
#if defined(PP_GLM) && defined (PP_NC_GLM)
  Phi_GLM_s2 = 0.5*SUM(M_f(:,i,j,k)*Uaux(IU:IW,i,j,k)) *(/U_in(IPSI,i,j,k),1./)
#endif /*PP_GLM and PP_NC_GLM*/
  
  ! First metrics dealiasing term
  Bhat_L = dot_product(M_f(:,i,j,k),U_in(IB1:IB3,i,j,k))
  
  DO l=0,PP_N
    f(vs:ve,l,i,j,k) = f(vs:ve,l,i,j,k) + Phi_MHD_s2 * (Bhat_L + dot_product(0.5*(M_f(:,i,j,k)+M_f(:,l,j,k)),U_in(IB1:IB3,l,j,k)))
#if defined(PP_GLM) && defined (PP_NC_GLM)
    !nonconservative term to restore Galilean invariance for GLM term
    f((/IRHOE,IPSI/),l,i,j,k) = f((/IRHOE,IPSI/),l,i,j,k) + (U_in(IPSI,i,j,k)+U_in(IPSI,l,j,k)) * Phi_GLM_s2
#endif /*PP_GLM and PP_NC_GLM*/
  END DO !l=0,PP_N
  
  ! Non-conservative terms in eta
  !------------------------------
#if defined(PP_GLM) && defined (PP_NC_GLM)
  Phi_GLM_s2 = 0.5*SUM(M_g(:,i,j,k)*Uaux(IU:IW,i,j,k)) *(/U_in(IPSI,i,j,k),1./)
#endif /*PP_GLM and PP_NC_GLM*/
  
  ! First metrics dealiasing term
  Bhat_L = dot_product(M_g(:,i,j,k),U_in(IB1:IB3,i,j,k))
  
  DO l=0,PP_N
    g(vs:ve,l,i,j,k) = g(vs:ve,l,i,j,k) + Phi_MHD_s2 * (Bhat_L + dot_product(0.5*(M_g(:,i,j,k)+M_g(:,i,l,k)),U_in(IB1:IB3,i,l,k)))
#if defined(PP_GLM) && defined (PP_NC_GLM)
    !nonconservative term to restore Galilean invariance for GLM term
    g((/IRHOE,IPSI/),l,i,j,k) = g((/IRHOE,IPSI/),l,i,j,k) + (U_in(IPSI,i,j,k)+U_in(IPSI,i,l,k)) * Phi_GLM_s2
#endif /*PP_GLM and PP_NC_GLM*/
  END DO !l=0,PP_N
  
  ! Non-conservative terms in zeta
  !-------------------------------
#if defined(PP_GLM) && defined (PP_NC_GLM)
  Phi_GLM_s2 = 0.5*SUM(M_h(:,i,j,k)*Uaux(IU:IW,i,j,k)) *(/U_in(IPSI,i,j,k),1./)
#endif /*PP_GLM and PP_NC_GLM*/
  
  ! First metrics dealiasing term
  Bhat_L = dot_product(M_h(:,i,j,k),U_in(IB1:IB3,i,j,k))
  
  DO l=0,PP_N
    h(vs:ve,l,i,j,k) = h(vs:ve,l,i,j,k) + Phi_MHD_s2 * (Bhat_L + dot_product(0.5*(M_h(:,i,j,k)+M_h(:,i,j,l)),U_in(IB1:IB3,i,j,l)))
#if defined(PP_GLM) && defined (PP_NC_GLM)
    !nonconservative term to restore Galilean invariance for GLM term
    h((/IRHOE,IPSI/),l,i,j,k) = h((/IRHOE,IPSI/),l,i,j,k) + (U_in(IPSI,i,j,k)+U_in(IPSI,i,j,l)) * Phi_GLM_s2
#endif /*PP_GLM and PP_NC_GLM*/
  END DO !l=0,PP_N
END DO; END DO; END DO ! i,j,k

END SUBROUTINE AddNonConsFluxTilde3D
!==================================================================================================================================
!> Compute transformed nonconservative MHD flux given left and right states and metrics:
!> phi^◇ = 0.5*(B_L·metrics_L+B_R·{metrics})*phi_L^MHD + {psi}*metrics_L·phi_L^GLM
!>
!> phi^MHD is the Powell, Brackbill or Janhunen nonconservative term:
!> * Powell:    phi^MHD = (0,B_1,B_2,B_3,v·B,v_1,v_2,v_3,0)
!> * Brackbill: phi^MHD = (0,B_1,B_2,B_3,0,0,0,0)
!> * Janhunen:  phi^MHD = (0,0,0,0,0,v_1,v_2,v_3,0)
!> phi^GLM is the GLM nonconservative term
!> * phi^GLM = (0,0,0,0,psi*v,0,0,0,v)
!==================================================================================================================================
PURE SUBROUTINE AddNonConsFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
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
REAL,DIMENSION(PP_nVar),INTENT(INOUT) :: Fstar   !< added to flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
real :: Bhat
!==================================================================================================================================

Bhat = 0.5*(dot_product(metric_L,UL(IB1:IB3))+dot_product(0.5*(metric_L+metric_R),UR(IB1:IB3)))

#if NONCONS==1 /*Powell*/
  ! Powell: Phi(IRHOU:IB3) =B,vB,v
  Fstar(IRHOU:IB3) = Fstar(IRHOU:IB3) +Bhat*(/UL(IB1:IB3),UauxL(IVB),UauxL(IRHOU:IRHOW)/)
#elif NONCONS==2 /*Brackbill*/
  ! Brackbill: Phi(IRHOU:IRHOW) =B
  Fstar(IRHOU:IRHOW) = Fstar(IRHOU:IRHOW) +Bhat*UL(IB1:IB3)
#elif NONCONS==3 /*Janhunen*/
  ! Janhunen: Phi(IB1:IB3) =v
  Fstar(IB1:IB3) = Fstar(IB1:IB3) +Bhat*UauxL(IU:IW)
#endif /*NONCONSTYPE*/
#if defined(PP_GLM) && defined (PP_NC_GLM)
  !nonconservative term to restore Galilean invariance for GLM term
  Fstar((/IRHOE,IPSI/)) = Fstar((/IRHOE,IPSI/)) +(0.5*(UL(IPSI)+UR(IPSI))*SUM(metric_L*UauxL(IU:IW))) *(/UL(IPSI),1./)
#endif /*PP_GLM and PP_NC_GLM*/

END SUBROUTINE AddNonConsFluxVec
#endif /*NONCONS*/

!==================================================================================================================================
!> Computes the standard flux in x-direction for the hyperbolic part ( normally used with a rotated state)
!==================================================================================================================================
PURE SUBROUTINE StandardDGFlux(UL,UR,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
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
ASSOCIATE(  rho_L =>UL(IRHO1),  rho_R =>UR(IRHO1), &
           rhoU_L =>UL(IRHOU), rhoU_R =>UR(IRHOU), &
           rhoV_L =>UL(IRHOV), rhoV_R =>UR(IRHOV), &
           rhoW_L =>UL(IRHOW), rhoW_R =>UR(IRHOW), &
#ifdef PP_GLM
              E_L =>UL(IRHOE)-0.5*smu_0*UL(IPSI)**2, E_R =>UR(IRHOE)-0.5*smu_0*UR(IPSI)**2, &
#else
              E_L =>UL(IRHOE),    E_R =>UR(IRHOE), &
#endif
             b1_L =>UL(IB1),   b1_R =>UR(IB1), &
             b2_L =>UL(IB2),   b2_R =>UR(IB2), &
             b3_L =>UL(IB3),   b3_R =>UR(IB3)  )
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
Fstar(IRHO1) = 0.5*(rhoqL      + rhoqR)
Fstar(IRHOU) = 0.5*(rhoqL*v1_L + rhoqR*v1_R +(pt_L + pt_R)-smu_0*(b1_L*b1_L+b1_R*b1_R))
Fstar(IRHOV) = 0.5*(rhoqL*v2_L + rhoqR*v2_R               -smu_0*(b2_L*b1_L+b2_R*b1_R))
Fstar(IRHOW) = 0.5*(rhoqL*v3_L + rhoqR*v3_R               -smu_0*(b3_L*b1_L+b3_R*b1_R))
Fstar(IRHOE) = 0.5*((E_L + pt_L)*v1_L + (E_R + pt_R)*v1_R- smu_0*(b1_L*vb_L+b1_R*vb_R))
Fstar(IB1) = 0.
Fstar(IB2) = 0.5*(v1_L*b2_L-b1_L*v2_L + v1_R*b2_R-b1_R*v2_R)
Fstar(IB3) = 0.5*(v1_L*b3_L-b1_L*v3_L + v1_R*b3_R-b1_R*v3_R)
#ifdef PP_GLM
Fstar(IRHOE) = Fstar(IRHOE)+0.5*smu_0*GLM_ch*(b1_L*UL(IPSI)+b1_R*UR(IPSI))
Fstar(IB1) = Fstar(IB1)+0.5      *GLM_ch*(     UL(IPSI)+     UR(IPSI))
Fstar(IPSI) =          0.5      *GLM_ch*(b1_L      +b1_R      )
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
USE MOD_Equation_Vars
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
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< transformed central flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: qv_L,qv_R,qb_L,qb_R
!==================================================================================================================================

! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>   UL(IRHO1),  rho_R =>   UR(IRHO1), &
           rhoU_L =>   UL(IRHOU), rhoU_R =>   UR(IRHOU), &
           rhoV_L =>   UL(IRHOV), rhoV_R =>   UR(IRHOV), &
           rhoW_L =>   UL(IRHOW), rhoW_R =>   UR(IRHOW), &
#ifdef PP_GLM
             E_L =>UL(IRHOE)-0.5*smu_0*UL(IPSI)**2, E_R =>UR(IRHOE)-0.5*smu_0*UR(IPSI)**2, &
#else
             E_L =>UL(IRHOE), E_R =>UR(IRHOE), &
#endif
             b1_L =>   UL(IB1),   b1_R =>   UR(IB1), &
             b2_L =>   UL(IB2),   b2_R =>   UR(IB2), &
             b3_L =>   UL(IB3),   b3_R =>   UR(IB3), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
             v1_L =>UauxL(IU),   v1_R =>UauxR(IU), &
             v2_L =>UauxL(IV),   v2_R =>UauxR(IV), &
             v3_L =>UauxL(IW),   v3_R =>UauxR(IW), &
             pt_L =>UauxL(IP),   pt_R =>UauxR(IP), & !total pressure = gas pressure+magnetic pressure
           ! v2_L =>UauxL(IVV),   v2_R =>UauxR(IVV), &
           ! b2_L =>UauxL(IBB),   b2_R =>UauxR(IBB), &
             vb_L =>UauxL(IVB),   vb_R =>UauxR(IVB)  )

!without metric dealiasing (=standard DG weak form on curved meshes)
qv_L = v1_L*metric_L(1) + v2_L*metric_L(2) + v3_L*metric_L(3)
qb_L = b1_L*metric_L(1) + b2_L*metric_L(2) + b3_L*metric_L(3)

qv_R = v1_R*metric_R(1) + v2_R*metric_R(2) + v3_R*metric_R(3)
qb_R = b1_R*metric_R(1) + b2_R*metric_R(2) + b3_R*metric_R(3)

! Standard DG flux
!without metric dealiasing (=standard DG weak form on curved meshes)
Fstar(IRHO1) = 0.5*( rho_L*qv_L +  rho_R*qv_R )
Fstar(IRHOU) = 0.5*(rhoU_L*qv_L + rhoU_R*qv_R + metric_L(1)*pt_L+metric_R(1)*pt_R -smu_0*(qb_L*b1_L+qb_R*b1_R) )
Fstar(IRHOV) = 0.5*(rhoV_L*qv_L + rhoV_R*qv_R + metric_L(2)*pt_L+metric_R(2)*pt_R -smu_0*(qb_L*b2_L+qb_R*b2_R) )
Fstar(IRHOW) = 0.5*(rhoW_L*qv_L + rhoW_R*qv_R + metric_L(3)*pt_L+metric_R(3)*pt_R -smu_0*(qb_L*b3_L+qb_R*b3_R) )
Fstar(IRHOE) = 0.5*((E_L + pt_L)*qv_L  + (E_R + pt_R)*qv_R      -smu_0*(qb_L*vb_L+qb_R*vb_R) )
Fstar(IB1) = 0.5*(qv_L*b1_L-qb_L*v1_L + qv_R*b1_R-qb_R*v1_R)
Fstar(IB2) = 0.5*(qv_L*b2_L-qb_L*v2_L + qv_R*b2_R-qb_R*v2_R)
Fstar(IB3) = 0.5*(qv_L*b3_L-qb_L*v3_L + qv_R*b3_R-qb_R*v3_R)

#ifdef PP_GLM
!without metric dealiasing (=standard DG weak form on curved meshes)
Fstar(IRHOE) = Fstar(IRHOE) + 0.5*smu_0*GLM_ch*(qb_L*UL(IPSI)             + qb_R*UR(IPSI))
Fstar(IB1) = Fstar(IB1) + 0.5      *GLM_ch*(     UL(IPSI)*metric_L(1) +      UR(IPSI)*metric_R(1))
Fstar(IB2) = Fstar(IB2) + 0.5      *GLM_ch*(     UL(IPSI)*metric_L(2) +      UR(IPSI)*metric_R(2))
Fstar(IB3) = Fstar(IB3) + 0.5      *GLM_ch*(     UL(IPSI)*metric_L(3) +      UR(IPSI)*metric_R(3))
Fstar(IPSI) =            0.5      *GLM_ch*(qb_L                   +qb_R                   )

#endif /* PP_GLM */

END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE StandardDGFluxVec

!==================================================================================================================================
!> Computes the standard DG flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the advection
!> part of the MHD equations
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
PURE SUBROUTINE StandardDGFluxDealiasedMetricVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
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
ASSOCIATE(  rho_L =>   UL(IRHO1),  rho_R =>   UR(IRHO1), &
           rhoU_L =>   UL(IRHOU), rhoU_R =>   UR(IRHOU), &
           rhoV_L =>   UL(IRHOV), rhoV_R =>   UR(IRHOV), &
           rhoW_L =>   UL(IRHOW), rhoW_R =>   UR(IRHOW), &
#ifdef PP_GLM
             E_L =>UL(IRHOE)-0.5*smu_0*UL(IPSI)**2, E_R =>UR(IRHOE)-0.5*smu_0*UR(IPSI)**2, &
#else
             E_L =>UL(IRHOE), E_R =>UR(IRHOE), &
#endif
             b1_L =>   UL(IB1),   b1_R =>   UR(IB1), &
             b2_L =>   UL(IB2),   b2_R =>   UR(IB2), &
             b3_L =>   UL(IB3),   b3_R =>   UR(IB3), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
             v1_L =>UauxL(IU),   v1_R =>UauxR(IU), &
             v2_L =>UauxL(IV),   v2_R =>UauxR(IV), &
             v3_L =>UauxL(IW),   v3_R =>UauxR(IW), &
             pt_L =>UauxL(IP),   pt_R =>UauxR(IP), & !total pressure = gas pressure+magnetic pressure
           !vv2_L =>UauxL(IVV),  vv2_R =>UauxR(IVV), &
           !bb2_L =>UauxL(IBB),  bb2_R =>UauxR(IBB), &
             vb_L =>UauxL(IVB),   vb_R =>UauxR(IVB)  )

qv_L = v1_L*metric(1) + v2_L*metric(2) + v3_L*metric(3)
qb_L = b1_L*metric(1) + b2_L*metric(2) + b3_L*metric(3)

qv_R = v1_R*metric(1) + v2_R*metric(2) + v3_R*metric(3)
qb_R = b1_R*metric(1) + b2_R*metric(2) + b3_R*metric(3)

! Standard DG flux
Fstar(IRHO1) = 0.5*( rho_L*qv_L +  rho_R*qv_R )
Fstar(IRHOU) = 0.5*(rhoU_L*qv_L + rhoU_R*qv_R + metric(1)*(pt_L+pt_R) -smu_0*(qb_L*b1_L+qb_R*b1_R) )
Fstar(IRHOV) = 0.5*(rhoV_L*qv_L + rhoV_R*qv_R + metric(2)*(pt_L+pt_R) -smu_0*(qb_L*b2_L+qb_R*b2_R) )
Fstar(IRHOW) = 0.5*(rhoW_L*qv_L + rhoW_R*qv_R + metric(3)*(pt_L+pt_R) -smu_0*(qb_L*b3_L+qb_R*b3_R) )
Fstar(IRHOE) = 0.5*((E_L + pt_L)*qv_L  + (E_R + pt_R)*qv_R      -smu_0*(qb_L*vb_L+qb_R*vb_R) )
Fstar(IB1) = 0.5*(qv_L*b1_L-qb_L*v1_L + qv_R*b1_R-qb_R*v1_R)
Fstar(IB2) = 0.5*(qv_L*b2_L-qb_L*v2_L + qv_R*b2_R-qb_R*v2_R)
Fstar(IB3) = 0.5*(qv_L*b3_L-qb_L*v3_L + qv_R*b3_R-qb_R*v3_R)

#ifdef PP_GLM
Fstar(IRHOE) = Fstar(IRHOE) + 0.5*smu_0*GLM_ch*(qb_L*UL(IPSI)+qb_R*UR(IPSI))
phiHat   = 0.5*GLM_ch*(UL(IPSI)+UR(IPSI))
Fstar(IB1) = Fstar(IB1) + phiHat*metric(1)
Fstar(IB2) = Fstar(IB2) + phiHat*metric(2)
Fstar(IB3) = Fstar(IB3) + phiHat*metric(3)
Fstar(IPSI) =            0.5*GLM_ch*(qb_L+qb_R)
#endif /* PP_GLM */

END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE StandardDGFluxDealiasedMetricVec



!==================================================================================================================================
!> entropy conservation for MHD, kinetric Energy conservation only in the Euler case
!> following D.Dergs et al."a novel Entropy consistent nine-wave field divergence diminishing ideal MHD system" 
!> mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t) 
!==================================================================================================================================
PURE SUBROUTINE EntropyAndKinEnergyConservingFlux_Derigs(UL,UR,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
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
ASSOCIATE(  rho_L =>   UL(IRHO1),  rho_R =>   UR(IRHO1), &
           rhoV_L => UL(IRHOU:IRHOW), rhoV_R => UR(IRHOU:IRHOW), &
#ifdef PP_GLM
              E_L =>UL(IRHOE)-0.5*smu_0*UL(IPSI)**2, E_R =>UR(IRHOE)-0.5*smu_0*UR(IPSI)**2, &
            psi_L =>UL(IPSI)   ,  psi_R =>UR(IPSI), &
#else
              E_L =>UL(IRHOE)   ,    E_R =>UR(IRHOE), &
#endif
              B_L => UL(IB1:IB3),    B_R => UR(IB1:IB3)  )
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
Fstar(IRHO1) = rhoLN*vAvg(1)
Fstar(IRHOU) = Fstar(IRHO1)*vAvg(1) - smu_0*BAvg(1)*BAvg(1) + pTilde
Fstar(IRHOV) = Fstar(IRHO1)*vAvg(2) - smu_0*BAvg(1)*BAvg(2)
Fstar(IRHOW) = Fstar(IRHO1)*vAvg(3) - smu_0*BAvg(1)*BAvg(3)
Fstar(IB2) = vAvg(1)*Bavg(2) - BAvg(1)*vAvg(2)
Fstar(IB3) = vAvg(1)*Bavg(3) - BAvg(1)*vAvg(3)
#ifdef PP_GLM
Fstar(IB1) = GLM_ch*psiAvg
Fstar(IPSI) = GLM_ch*BAvg(1)
#else
Fstar(IB1) =0.
#endif

Fstar(IRHOE) = Fstar(IRHO1)*0.5*(skappaM1/betaLN - 0.5*(v2_L+v2_R))  &
           + SUM(vAvg(:)*Fstar(IRHOU:IRHOW)) &
           +smu_0*( SUM(BAvg(:)*Fstar(IB1:IB3)) &
                   -0.5*v1_B2Avg +BAvg(1)*vB_Avg &
#ifdef PP_GLM
                   +Fstar(IPSI)*psiAvg-GLM_ch*0.5*(psi_L*B_L(1)+psi_R*B_R(1))    &
#endif
                   )

END ASSOCIATE 
END SUBROUTINE EntropyAndKinEnergyConservingFlux_Derigs


!==================================================================================================================================
!> entropy conservation for MHD, kinetric Energy conservation only in the Euler case
!> following D.Dergs et al."a novel Entropy consistent nine-wave field divergence diminishing ideal MHD system" 
!> mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t) 
!> firectly compute tranformed flux: fstar=f*metric1+g*metric2+h*metric3
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
PURE SUBROUTINE EntropyAndKinEnergyConservingFluxVec_Derigs(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
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

ASSOCIATE(  rho_L => UL(IRHO1)  ,  rho_R => UR(IRHO1)    , &
              B_L => UL(IB1:IB3),    B_R => UR(IB1:IB3)  , &
#ifdef PP_GLM
            psi_L => UL(IPSI)  ,  psi_R => UR(IPSI)    , &
#endif /*PP_GLM*/
              v_L =>UauxL(IU:IW),  v_R =>UauxR(IU:IW), &
             pt_L =>UauxL(IP),   pt_R =>UauxR(IP)  , & !pt=p+1/(2mu_0)|B|^2
             v2_L =>UauxL(IVV),   v2_R =>UauxR(IVV)  , & !|v|^2 left/right
             B2_L =>UauxL(IBB),   B2_R =>UauxR(IBB)  , & !|B|^2 left/right
             vB_L =>UauxL(IVB),   vB_R =>UauxR(IVB)  )

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
Fstar(IRHO1) = rhoLN*vm
Fstar(IRHOU:IRHOW) = Fstar(IRHO1)*vAvg(1:3)-(smu_0*Bm)*BAvg(:) + pTilde*metric(1:3)
#ifdef PP_GLM
Fstar(IB1:IB3) = vm*BAvg(1:3) - Bm*vAvg(1:3) + (GLM_ch*PsiAvg)*metric(1:3)
Fstar(IPSI) = GLM_ch*Bm
#else
Fstar(IB1:IB3) = vm*BAvg(1:3) - Bm*vAvg(1:3)
#endif /*PP_GLM*/


Fstar(IRHOE) = Fstar(IRHO1)*0.5*(skappaM1/betaLN - 0.5*(v2_L+v2_R)) &
           + SUM(vAvg(:)*Fstar(IRHOU:IRHOW))  &
           +smu_0*( SUM(BAvg(:)*Fstar(IB1:IB3))                        &
                   - 0.25*SUM((B2_L*v_L(:)+B2_R*v_R(:))*metric(:)) & ! -0.5* {{|B|^2v(:)}}.{{m(:)}}
                   + 0.5*(vb_L+vb_R)*bm                            & !{{(v.B)}}{{ B(:) }} . {{m(:)}}
#ifdef PP_GLM
                   +Fstar(IPSI)*PsiAvg - (GLM_ch*0.5)*SUM((Psi_L*B_L(:)+Psi_R*B_R(:))*metric(:))  & !c_h{{psi B}}.{{m}}
#endif
                  )
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE EntropyAndKinEnergyConservingFluxVec_Derigs

!==================================================================================================================================
!> entropy conservation for MHD, kinetric Energy conservation only in the Euler case
!> following D.Dergs et al."a novel Entropy consistent nine-wave field divergence diminishing ideal MHD system" 
!> mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t) 
!==================================================================================================================================
PURE SUBROUTINE EntropyAndKinEnergyConservingFlux_FloGor(UL,UR,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
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
REAL            :: srho_L,srho_R,rhoLN,B2_L,B2_R,v2_L,v2_R
REAL            :: p_L,p_R,p_avg
REAL            :: v_L(3),v_R(3)
REAL            :: BAvg(3),vAvg(3)
REAL            :: in_e_L,in_e_R
REAL            :: B2_ZIP,v2_ZIP
!==================================================================================================================================
ASSOCIATE(  rho_L =>   UL(IRHO1),  rho_R =>   UR(IRHO1), &
           rhoV_L => UL(IRHOU:IRHOW), rhoV_R => UR(IRHOU:IRHOW), &
#ifdef PP_GLM
              E_L =>UL(IRHOE)-0.5*smu_0*UL(IPSI)**2, E_R =>UR(IRHOE)-0.5*smu_0*UR(IPSI)**2, &
            psi_L =>UL(IPSI)   ,  psi_R =>UR(IPSI), &
#else
              E_L =>UL(IRHOE)   ,    E_R =>UR(IRHOE), &
#endif
              B_L => UL(IB1:IB3),    B_R => UR(IB1:IB3)  )
! Get the inverse density, velocity, and pressure on left and right
srho_L = 1./rho_L
srho_R = 1./rho_R
v_L = rhoV_L(:)*srho_L
v_R = rhoV_R(:)*srho_R


v2_L = SUM(v_L(:)*v_L(:))
v2_R = SUM(v_R(:)*v_R(:))
B2_L = SUM(B_L(:)*B_L(:))
B2_R = SUM(B_R(:)*B_R(:))

p_L    = kappaM1*(E_L - 0.5*(rho_L*v2_L+smu_0*B2_L))
p_R    = kappaM1*(E_R - 0.5*(rho_R*v2_R+smu_0*B2_R))

!specific inner energy
in_e_L = p_L*srho_L*skappaM1
in_e_R = p_R*srho_R*skappaM1

! Get the averages for the numerical flux

rhoLN      = LN_MEAN( rho_L, rho_R)
vAvg       = 0.5 * ( v_L +  v_R)
v2_ZIP   = (v_L(1)*v_R(1)+v_L(2)*v_R(2)+v_L(3)*v_R(3))
BAvg       = 0.5 * ( B_L +  B_R)
B2_ZIP   = (B_L(1)*B_R(1)+B_L(2)*B_R(2)+B_L(3)*B_R(3))
                                                                   
p_avg     = 0.5*(p_L+p_R)

#define ZIP(a,b,c,d) 0.5*(a*d+b*c)
! Entropy conserving and kinetic energy conserving flux
Fstar(IRHO1) = rhoLN*vAvg(1)

Fstar(IRHOU) = Fstar(IRHO1)*vAvg(1)+p_avg+s2mu_0*B2_ZIP- smu_0*ZIP(B_L(1),B_R(1),B_L(1),B_R(1))
Fstar(IRHOV) = Fstar(IRHO1)*vAvg(2)                    - smu_0*ZIP(B_L(1),B_R(1),B_L(2),B_R(2))
Fstar(IRHOW) = Fstar(IRHO1)*vAvg(3)                    - smu_0*ZIP(B_L(1),B_R(1),B_L(3),B_R(3))

!opt without ZIP
!Fstar(IRHOU) = Fstar(IRHO1)*vAvg(1)+p_avg+s2mu_0*B2_ZIP- s2mu_0*(B_L(1)*B_R(1)+B_R(1)*B_L(1)) !1/2 of ZIP in s2mu_0=smu_0*0.5
!Fstar(IRHOV) = Fstar(IRHO1)*vAvg(2)                    - s2mu_0*(B_L(1)*B_R(2)+B_R(1)*B_L(2))
!Fstar(IRHOW) = Fstar(IRHO1)*vAvg(3)                    - s2mu_0*(B_L(1)*B_R(3)+B_R(1)*B_L(3))

Fstar(IRHOE) = Fstar(IRHO1)*(0.5*v2_ZIP+in_e_L*in_e_R/LN_MEAN(in_e_L,in_e_R))+ZIP(p_L,p_R,v_L(1),v_R(1))  &
         + smu_0*(ZIP(v_L(1)*B_L(2),v_R(1)*B_R(2),B_L(2),B_R(2))-ZIP(v_L(2)*B_L(1),v_R(2)*B_R(1),B_L(2),B_R(2)) &
         +        ZIP(v_L(1)*B_L(3),v_R(1)*B_R(3),B_L(3),B_R(3))-ZIP(v_L(3)*B_L(1),v_R(3)*B_R(1),B_L(3),B_R(3)) &
#ifdef PP_GLM
                  +GLM_ch*ZIP(B_L(1),B_R(1),psi_L,psi_R)  &
#endif /*PP_GLM*/
                  )

!opt without ZIP
!Fstar(IRHOE) = Fstar(IRHO1)*(0.5*v2_ZIP+in_e_L*in_e_R/LN_MEAN(in_e_L,in_e_R))+0.5*(p_L*v_R(1)+p_R*v_L(1))  &
!         + s2mu_0*((v_L(1)+v_R(1))*B2_ZIP                              &  !1/2 of ZIP in s2mu_0=smu_0*0.5
!                   -B_L(1)*(v_L(1)*B_R(1)+v_L(2)*B_R(2)+v_L(3)*B_R(3)) &
!                   -B_R(1)*(B_L(1)*v_R(1)+B_L(2)*v_R(2)+B_L(3)*v_R(3)) &
!#ifdef PP_GLM
!                   +GLM_ch*(B_L(1)*psi_R+B_R(1)*psi_L)  &
!#endif /*PP_GLM*/
!                  )

#ifdef PP_GLM
Fstar(IB1) = GLM_ch*0.5*(psi_L+psi_R)
Fstar(IPSI) = GLM_ch*BAvg(1)
#else
Fstar(IB1) = 0.
#endif /*PP_GLM*/
Fstar(IB2) = 0.5* ((v_L(1)*B_L(2)-v_L(2)*B_L(1)) + (v_R(1)*B_R(2)-v_R(2)*B_R(1)))
Fstar(IB3) = 0.5* ((v_L(1)*B_L(3)-v_L(3)*B_L(1)) + (v_R(1)*B_R(3)-v_R(3)*B_R(1)))

#undef ZIP

END ASSOCIATE 
END SUBROUTINE EntropyAndKinEnergyConservingFlux_FloGor


!==================================================================================================================================
!> entropy conservation for MHD, kinetric Energy conservation only in the Euler case
!> following D.Dergs et al."a novel Entropy consistent nine-wave field divergence diminishing ideal MHD system" 
!> mu_0 added, total energy contribution is 1/(2mu_0)(|B|^2+psi^2), in energy flux: 1/mu_0*(B.B_t + psi*psi_t) 
!> firectly compute tranformed flux: fstar=f*metric1+g*metric2+h*metric3
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
PURE SUBROUTINE EntropyAndKinEnergyConservingFluxVec_FloGor(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars
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
REAL                   :: rhoLN
REAL                   :: vm_L,vm_R,Bm_L,Bm_R
REAL                   :: p_L,p_R,in_e_L,in_e_R
REAL                   :: p_avg,v2_ZIP,B2_ZIP
REAL                   :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)

ASSOCIATE(  rho_L => UL(IRHO1)  ,  rho_R => UR(IRHO1)    , &
              B_L => UL(IB1:IB3),    B_R => UR(IB1:IB3)  , &
#ifdef PP_GLM
            psi_L => UL(IPSI)  ,  psi_R => UR(IPSI)    , &
#endif /*PP_GLM*/
              v_L =>UauxL(IU:IW),  v_R =>UauxR(IU:IW), &
             pt_L =>UauxL(IP),   pt_R =>UauxR(IP)  , & !pt=p+1/(2mu_0)|B|^2
             v2_L =>UauxL(IVV),   v2_R =>UauxR(IVV)  , & !|v|^2 left/right
             B2_L =>UauxL(IBB),   B2_R =>UauxR(IBB)  , & !|B|^2 left/right
             vB_L =>UauxL(IVB),   vB_R =>UauxR(IVB)  )

p_L=pt_L -s2mu_0*B2_L
p_R=pt_R -s2mu_0*B2_R

! specific inner energy
in_e_L = p_L/rho_L*skappaM1
in_e_R = p_R/rho_R*skappaM1

! Get the averages for the numerical flux
rhoLN     = LN_MEAN( rho_L, rho_R)
vAvg      = 0.5*( v_L(:)+ v_R(:))
BAvg      = 0.5*( B_L(:)+ B_R(:))
v2_ZIP   = (v_L(1)*v_R(1)+v_L(2)*v_R(2)+v_L(3)*v_R(3))
B2_ZIP   = (B_L(1)*B_R(1)+B_L(2)*B_R(2)+B_L(3)*B_R(3))
p_avg      = 0.5*(p_L+p_R)


vm_L=SUM(v_L(:)*metric(:))
vm_R=SUM(v_R(:)*metric(:))
Bm_L=SUM(B_L(:)*metric(:))
Bm_R=SUM(B_R(:)*metric(:))

#define ZIP(a,b,c,d) 0.5*(a*d+b*c)
! Entropy conserving and kinetic energy conserving flux
Fstar(IRHO1) = rhoLN*0.5*(vm_L+vm_R)

Fstar(IRHOU) = Fstar(IRHO1)*vAvg(1)+metric(1)*(p_avg+s2mu_0*B2_ZIP)-smu_0*ZIP(B_L(1),B_R(1),Bm_L,Bm_R)
Fstar(IRHOV) = Fstar(IRHO1)*vAvg(2)+metric(2)*(p_avg+s2mu_0*B2_ZIP)-smu_0*ZIP(B_L(2),B_R(2),Bm_L,Bm_R)
Fstar(IRHOW) = Fstar(IRHO1)*vAvg(3)+metric(3)*(p_avg+s2mu_0*B2_ZIP)-smu_0*ZIP(B_L(3),B_R(3),Bm_L,Bm_R)

!opt without ZIP
!Fstar(IRHOU) = Fstar(IRHO1)*vAvg(1)+metric(1)*(p_avg+s2mu_0*B2_ZIP)-s2mu_0*(B_L(1)*Bm_R+B_R(1)*Bm_L) !1/2 of ZIP in s2mu_0=smu_0*0.5
!Fstar(IRHOV) = Fstar(IRHO1)*vAvg(2)+metric(2)*(p_avg+s2mu_0*B2_ZIP)-s2mu_0*(B_L(2)*Bm_R+B_R(2)*Bm_L)
!Fstar(IRHOW) = Fstar(IRHO1)*vAvg(3)+metric(3)*(p_avg+s2mu_0*B2_ZIP)-s2mu_0*(B_L(3)*Bm_R+B_R(3)*Bm_L)

Fstar(IRHOE) = Fstar(IRHO1)*(0.5*v2_ZIP+in_e_L*in_e_R/LN_MEAN(in_e_L,in_e_R))&
         + ZIP(p_L,p_R,vm_L,vm_R) &
         + smu_0*(  ZIP(vm_L*B_L(1),vm_R*B_R(1),B_L(1),B_R(1))-ZIP(v_L(1)*Bm_L,v_R(1)*Bm_R,B_L(1),B_R(1)) &
                  + ZIP(vm_L*B_L(2),vm_R*B_R(2),B_L(2),B_R(2))-ZIP(v_L(2)*Bm_L,v_R(2)*Bm_R,B_L(2),B_R(2)) &
                  + ZIP(vm_L*B_L(3),vm_R*B_R(3),B_L(3),B_R(3))-ZIP(v_L(3)*Bm_L,v_R(3)*Bm_R,B_L(3),B_R(3)) &
#ifdef PP_GLM
                  + GLM_ch*ZIP(Bm_L,Bm_R,Psi_L,Psi_R)                  & 
#endif /*PP_GLM*/
                 )

!opt without ZIP
!Fstar(IRHOE) = Fstar(IRHO1)*(0.5*v2_ZIP+in_e_L*in_e_R/LN_MEAN(in_e_L,in_e_R))   + 0.5*(p_L*vm_R+p_R*vm_L) &
!         + s2mu_0*((vm_L+vm_R)*B2_ZIP                                &  !1/2 of ZIP in s2mu_0=smu_0*0.5
!                   -Bm_L*(v_L(1)*B_R(1)+v_L(2)*B_R(2)+v_L(3)*B_R(3)) &
!                   -Bm_R*(B_L(1)*v_R(1)+B_L(2)*v_R(2)+B_L(3)*v_R(3)) &
!#ifdef PP_GLM
!                   + GLM_ch*(Bm_L*Psi_R+Bm_R*Psi_L)                  & 
!#endif /*PP_GLM*/
!                 )
#undef ZIP

#ifdef PP_GLM
Fstar(IB1:IB3) = 0.5* ((vm_L*B_L(1:3)-v_L(1:3)*Bm_L) + (vm_R*B_R(1:3)-v_R(1:3)*Bm_R)  + GLM_ch*(Psi_L+Psi_R)*metric(1:3))
Fstar(IPSI) = GLM_ch*0.5*(Bm_L+Bm_R)
#else                                                            
Fstar(IB1:IB3) = 0.5* ((vm_L*B_L(1:3)-v_L(1:3)*Bm_L) + (vm_R*B_R(1:3)-v_R(1:3)*Bm_R)) 
#endif /*PP_GLM*/
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE EntropyAndKinEnergyConservingFluxVec_FloGor


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

#undef PP_VolumeFluxAverageVec

END MODULE MOD_Flux_Average
