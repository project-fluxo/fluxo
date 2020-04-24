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
!> Computes two-point average fluxes for the volint when using the split-form (DiscType=2)
!==================================================================================================================================
MODULE MOD_Flux_Average
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE EvalEulerFluxTilde3D
  MODULE PROCEDURE EvalEulerFluxTilde3D
END INTERFACE

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

INTERFACE TwoPointEntropyConservingFlux
  MODULE PROCEDURE TwoPointEntropyConservingFlux
END INTERFACE

INTERFACE TwoPointEntropyConservingFluxVec
  MODULE PROCEDURE TwoPointEntropyConservingFluxVec
END INTERFACE

INTERFACE KennedyAndGruberFlux1
  MODULE PROCEDURE KennedyAndGruberFlux1
END INTERFACE

INTERFACE KennedyAndGruberFluxVec1
  MODULE PROCEDURE KennedyAndGruberFluxVec1
END INTERFACE

INTERFACE KennedyAndGruberFlux2
  MODULE PROCEDURE KennedyAndGruberFlux2
END INTERFACE

INTERFACE KennedyAndGruberFluxVec2
  MODULE PROCEDURE KennedyAndGruberFluxVec2
END INTERFACE

INTERFACE DucrosFlux
  MODULE PROCEDURE DucrosFlux
END INTERFACE

INTERFACE DucrosFluxVec
  MODULE PROCEDURE DucrosFluxVec
END INTERFACE

INTERFACE MorinishiFlux
  MODULE PROCEDURE MorinishiFlux
END INTERFACE

INTERFACE MorinishiFluxVec
  MODULE PROCEDURE MorinishiFluxVec
END INTERFACE

INTERFACE ggFlux
  MODULE PROCEDURE ggFlux
END INTERFACE

INTERFACE ggFluxVec
  MODULE PROCEDURE ggFluxVec
END INTERFACE

INTERFACE EntropyAndEnergyConservingFlux
  MODULE PROCEDURE EntropyAndEnergyConservingFlux
END INTERFACE

INTERFACE EntropyAndEnergyConservingFluxVec
  MODULE PROCEDURE EntropyAndEnergyConservingFluxVec
END INTERFACE

INTERFACE EntropyAndEnergyConservingFlux2
  MODULE PROCEDURE EntropyAndEnergyConservingFlux2
END INTERFACE

INTERFACE EntropyAndEnergyConservingFluxVec2
  MODULE PROCEDURE EntropyAndEnergyConservingFluxVec2
END INTERFACE

INTERFACE GassnerWintersWalchFlux
   MODULE PROCEDURE GassnerWintersWalchFlux
END INTERFACE

INTERFACE GassnerWintersWalchFluxVec
   MODULE PROCEDURE GassnerWintersWalchFluxVec
END INTERFACE

INTERFACE LN_MEAN 
   MODULE PROCEDURE LN_MEAN
END INTERFACE


PUBLIC:: EvalEulerFluxTilde3D
PUBLIC:: EvalUaux
PUBLIC:: StandardDGFlux
PUBLIC:: StandardDGFluxVec
PUBLIC:: StandardDGFluxDealiasedMetricVec
PUBLIC:: TwoPointEntropyConservingFlux
PUBLIC:: TwoPointEntropyConservingFluxVec
PUBLIC:: KennedyAndGruberFlux1
PUBLIC:: KennedyAndGruberFluxVec1
PUBLIC:: KennedyAndGruberFlux2
PUBLIC:: KennedyAndGruberFluxVec2
PUBLIC:: DucrosFlux
PUBLIC:: DucrosFluxVec
PUBLIC:: MorinishiFlux
PUBLIC:: MorinishiFluxVec
PUBLIC:: EntropyAndEnergyConservingFlux
PUBLIC:: EntropyAndEnergyConservingFluxVec
PUBLIC:: EntropyAndEnergyConservingFlux2
PUBLIC:: EntropyAndEnergyConservingFluxVec2
PUBLIC:: ggflux
PUBLIC:: ggfluxVec
PUBLIC:: GassnerWintersWalchFlux
PUBLIC:: GassnerWintersWalchFluxVec
PUBLIC:: LN_MEAN

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute Euler fluxes using the conservative variables and derivatives for every volume Gauss point.
!> directly apply metrics and output the tranformed flux 
!==================================================================================================================================
SUBROUTINE EvalEulerFluxTilde3D(iElem,ftilde,gtilde,htilde,Uaux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Equation_Vars ,ONLY:nAuxVar,kappaM1
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
#ifdef OPTIMIZED
USE MOD_DG_Vars       ,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: iElem                              !< current element index in volint
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(    5,0:PP_N,0:PP_N,0:PP_N) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(    5,0:PP_N,0:PP_N,0:PP_N) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(    5,0:PP_N,0:PP_N,0:PP_N) !< transformed flux h(iVar,i,j,k)
REAL,INTENT(OUT)   :: Uaux(nAuxVar,0:PP_N,0:PP_N,0:PP_N) !< auxiliary variables:(srho,v1,v2,v3,p,|v|^2)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: f(5),g(5),h(5)                    !cartesi an fluxes 
REAL                :: srho                              ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,v_2,p                    ! auxiliary variables
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
  ASSOCIATE(rho   =>U(1,PP_IJK,iElem), &
            rhov1 =>U(2,PP_IJK,iElem), &
            rhov2 =>U(3,PP_IJK,iElem), &
            rhov3 =>U(4,PP_IJK,iElem), &
            rhoE  =>U(5,PP_IJK,iElem)  )
  ! auxiliary variables
  srho = 1./rho
  v1   = rhov1*srho 
  v2   = rhov2*srho 
  v3   = rhov3*srho 
  v_2  = v1*v1+v2*v2+v3*v3 
  p    = kappaM1*(rhoE-0.5*rho*v_2)
  Uaux(:,PP_IJK)=(/srho,v1,v2,v3,p,v_2/)
  ! Euler part
  ! Euler fluxes x-direction
  f(1)=rhov1         
  f(2)=rhov1*v1+p    
  f(3)=rhov1*v2      
  f(4)=rhov1*v3      
  f(5)=(rhoE+p)*v1         
  ! Euler fluxes y-direction
  g(1)=rhov2
  g(2)=f(3)                      ! rho*u*v
  g(3)=rhov2*v2+p  
  g(4)=rhov2*v3
  g(5)=(rhoE+p)*v2 
  ! Euler fluxes z-direction
  h(1)=rhov3
  h(2)=f(4)               ! rho*v1*v3
  h(3)=g(4)               ! rho*v2*v3  
  h(4)=rhov3*v3+p    
  h(5)=(rhoE+p)*v3  
  END ASSOCIATE !rho,rhov1,rhov2,rhov3,rhoE
  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(:,PP_IJK) =   f(:)*Metrics_fTilde(1,PP_IJK,iElem)  &
                     + g(:)*Metrics_fTilde(2,PP_IJK,iElem)  &
                     + h(:)*Metrics_fTilde(3,PP_IJK,iElem)
  gtilde(:,PP_IJK) =   f(:)*Metrics_gTilde(1,PP_IJK,iElem)  &
                     + g(:)*Metrics_gTilde(2,PP_IJK,iElem)  &
                     + h(:)*Metrics_gTilde(3,PP_IJK,iElem)
  htilde(:,PP_IJK) =   f(:)*Metrics_hTilde(1,PP_IJK,iElem)  &
                     + g(:)*Metrics_hTilde(2,PP_IJK,iElem)  &
                     + h(:)*Metrics_hTilde(3,PP_IJK,iElem)

#ifdef OPTIMIZED
END DO ! i
#else /*OPTIMIZED*/
END DO; END DO; END DO ! i,j,k
#endif /*OPTIMIZED*/
END SUBROUTINE EvalEulerFluxTilde3D


!==================================================================================================================================
!> computes auxiliary nodal variables (1/rho,v_1,v_2,v_3,p,|v|^2) from state U
!==================================================================================================================================
SUBROUTINE EvalUaux(iElem,Uaux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Equation_Vars ,ONLY:nAuxVar,kappaM1
#ifdef OPTIMIZED
USE MOD_DG_Vars,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: iElem                               !<current element in volint
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: Uaux(nAuxVar,0:PP_N,0:PP_N,0:PP_N)  !<auxiliary variables:(srho,v1,v2,v3,p,|v|^2)
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
  Uaux(1  ,PP_IJK) = 1./U(1,PP_IJK,iElem)
  Uaux(2:4,PP_IJK) = Uaux(1,PP_IJK)*U(2:4,PP_IJK,iElem) 
  Uaux(6  ,PP_IJK) = Uaux(2,PP_IJK)*Uaux(2,PP_IJK)+Uaux(3,PP_IJK)*Uaux(3,PP_IJK)+Uaux(4,PP_IJK)*Uaux(4,PP_IJK)
  Uaux(5  ,PP_IJK)=kappaM1*(U(5,PP_IJK,iElem)-0.5*U(1,PP_IJK,iElem)*Uaux(6,PP_IJK))
#ifdef OPTIMIZED
END DO ! i
#else /*OPTIMIZED*/
END DO; END DO; END DO ! i,j,k
#endif /*OPTIMIZED*/
END SUBROUTINE EvalUaux


!==================================================================================================================================
!> Computes the standard flux in x-direction for the Euler equations ( normally used with a rotated state)
!==================================================================================================================================
SUBROUTINE StandardDGFlux(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !< central flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: rhoqL,rhoqR
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,p_L,p_R,VelV_L,VelV_R,VelW_L,VelW_R
!==================================================================================================================================
! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>UL(1),  rho_R =>UR(1), &
           rhoU_L =>UL(2), rhoU_R =>UR(2), &
           rhoV_L =>UL(3), rhoV_R =>UR(3), &
           rhoW_L =>UL(4), rhoW_R =>UR(4), &
           rhoE_L =>UL(5), rhoE_R =>UR(5)  )
sRho_L = 1./rho_L
sRho_R = 1./rho_R
VelU_L = sRho_L*rhoU_L;  VelV_L = sRho_L*rhoV_L ; VelW_L = sRho_L*rhoW_L
VelU_R = sRho_R*rhoU_R;  VelV_R = sRho_R*rhoV_R ; VelW_R = sRho_R*rhoW_R

p_L    = (kappaM1)*(rhoE_L - 0.5*(rhoU_L*VelU_L + rhoV_L*VelV_L + rhoW_L*VelW_L))
p_R    = (kappaM1)*(rhoE_R - 0.5*(rhoU_R*VelU_R + rhoV_R*VelV_R + rhoW_R*VelW_R))

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = 0.5*(rho_L + rho_R)
uHat   = 0.5*(VelU_L + VelU_R)
vHat   = 0.5*(VelV_L + VelV_R)
wHat   = 0.5*(VelW_L + VelW_R)
p1Hat  = 0.5*(p_L + p_R)
aHat   = SQRT(kappa*p1Hat/rhoHat)
HHat   = kappa*p1Hat/(rhoHat*kappaM1) + 0.5*(uHat*uHat + vHat*vHat + wHat*wHat)
! Standard DG flux
rhoqL    = rho_L*VelU_L
rhoqR    = rho_R*VelU_R
Fstar(1) = 0.5*(rhoqL        + rhoqR)
Fstar(2) = 0.5*(rhoqL*VelU_L + rhoqR*VelU_R +(p_L + p_R))
Fstar(3) = 0.5*(rhoqL*VelV_L + rhoqR*VelV_R             )
Fstar(4) = 0.5*(rhoqL*VelW_L + rhoqR*VelW_R             )
Fstar(5) = 0.5*((rhoE_L + p_L)*VelU_L + (rhoE_R + p_R)*VelU_R)
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE StandardDGFlux


!==================================================================================================================================
!> Computes the standard DG euler flux transformed with the metrics 
!> fstar=1/2((fL*metric1L+gL*metric2L+h*metric3L)+(fR*metric1R+gR*metric2R+h*metric3R)  )
!==================================================================================================================================
SUBROUTINE StandardDGFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
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
REAL                                :: q_L,q_R
!==================================================================================================================================

! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5)  )


!curved, without metric dealiasing (=standard DG weak form on curved meshes)
q_L   = VelU_L*metric_L(1) + VelV_L*metric_L(2) + VelW_L*metric_L(3)
q_R   = VelU_R*metric_R(1) + VelV_R*metric_R(2) + VelW_R*metric_R(3)

Fstar(1) = 0.5*( rho_L*q_L  +  rho_R*q_R                         )
Fstar(2) = 0.5*(rhoU_L*q_L  + rhoU_R*q_R + metric_L(1)*p_L+metric_R(1)*p_R   )
Fstar(3) = 0.5*(rhoV_L*q_L  + rhoV_R*q_R + metric_L(2)*p_L+metric_R(2)*p_R   )
Fstar(4) = 0.5*(rhoW_L*q_L  + rhoW_R*q_R + metric_L(3)*p_L+metric_R(3)*p_R   )
Fstar(5) = 0.5*((rhoE_L + p_L)*q_L    + (rhoE_R + p_R)*q_R )

END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE StandardDGFluxVec


!==================================================================================================================================
!> Computes the standard DG euler flux with dealiased metrics (fstar=f*metric1+g*metric2+h*metric3 ) 
!> where for curved metrics, metric=1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE StandardDGFluxDealiasedMetricVec(UL,UR,UauxL,UauxR,metric_L,metric_R, Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
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
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar          !< transformed flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: q_L,q_R
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)

! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5)  )


q_L   = VelU_L*metric(1) + VelV_L*metric(2) + VelW_L*metric(3)
q_R   = VelU_R*metric(1) + VelV_R*metric(2) + VelW_R*metric(3)

! Standard DG flux
Fstar(1) = 0.5*( rho_L*q_L  +  rho_R*q_R                         )
Fstar(2) = 0.5*(rhoU_L*q_L  + rhoU_R*q_R + metric(1)*(p_L+p_R)   )
Fstar(3) = 0.5*(rhoV_L*q_L  + rhoV_R*q_R + metric(2)*(p_L+p_R)   )
Fstar(4) = 0.5*(rhoW_L*q_L  + rhoW_R*q_R + metric(3)*(p_L+p_R)   )
Fstar(5) = 0.5*((rhoE_L + p_L)*q_L    + (rhoE_R + p_R)*q_R )


END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE StandardDGFluxDealiasedMetricVec


!==================================================================================================================================
!> Computes the entropy conserving numerical 3D flux (in the normal direction) for the Euler equations derived by Ismail and Roe
!==================================================================================================================================
SUBROUTINE TwoPointEntropyConservingFlux(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,kappaP1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: z1L,z1R,z2L,z2R,z3L,z3R,z4L,z4R,z5L,z5R
REAL                                :: sz1Mean,z1Mean,z4Mean,z5LN,z3Mean,z1LN,z2Mean,z5Mean
REAL                                :: kappa_p2Hat
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,p_L,p_R,VelV_L,VelV_R,VelW_L,VelW_R
!==================================================================================================================================
! Get the inverse density, velocity, and pressure on left and right
sRho_L = 1./UL(1); VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1); VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

p_L    = (kappaM1)*(UL(5) - 0.5*(UL(2)*VelU_L + UL(3)*VelV_L + UL(4)*VelW_L))
p_R    = (kappaM1)*(UR(5) - 0.5*(UR(2)*VelU_R + UR(3)*VelV_R + UR(4)*VelW_R))

! z1 = √(rho/pressure) values on left and right, arithmatic and logorithmic means
z1L    = SQRT(UL(1)/p_L)
z1R    = SQRT(UR(1)/p_R)
z1Mean = 0.5*(z1L+z1R)
z1LN  = LN_Mean(z1L,z1R)
! z2 = u*√(rho/pressure) values on left and right and arithmatic mean
z2L    = z1L*VelU_L
z2R    = z1R*VelU_R
z2Mean = 0.5*(z2L+z2R)
! z3 = v*√(rho/pressure) values on left and right and arithmatic mean
z3L    = z1L*VelV_L
z3R    = z1R*VelV_R
z3Mean = 0.5*(z3L+z3R)
! z4 = w*√(rho/pressure) values on left and right and arithmatic mean
z4L    = z1L*VelW_L
z4R    = z1R*VelW_R
z4Mean = 0.5*(z4L+z4R)
! z5 = √(rho*pressure) values on left and right, arithmatic and logorithmic means
z5L    = p_L*z1L !SQRT(UL(1)*p_L)
z5R    = p_R*z1R !SQRT(UR(1)*p_R)
z5Mean = 0.5*(z5L+z5R)
z5LN   = LN_Mean(z5L,z5R)

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = z1Mean*z5LN
sz1Mean = 1./z1Mean
uHat   = z2Mean*sz1Mean
vHat   = z3Mean*sz1Mean
wHat   = z4Mean*sz1Mean
p1Hat  = z5Mean*sz1Mean
kappa_p2Hat  = 0.5*(kappaP1*(z5LN/z1LN) + kappaM1*p1Hat)
aHat   = SQRT(kappa_p2Hat/rhoHat)
HHat   = kappa_p2Hat/(rhoHat*kappaM1) + 0.5*(uHat*uHat + vHat*vHat + wHat*wHat)
! Entropy conserving flux
Fstar(1) = rhoHat*uHat
Fstar(2) = Fstar(1)*uHat + p1Hat
Fstar(3) = Fstar(1)*vHat 
Fstar(4) = Fstar(1)*wHat 
Fstar(5) = Fstar(1)*HHat
END SUBROUTINE TwoPointEntropyConservingFlux


!==================================================================================================================================
!> Computes the entropy conserving flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the Euler equations
!> derived by Ismail and Roe
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE TwoPointEntropyConservingFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
USE MOD_Equation_Vars,ONLY:skappaM1,kappaP1
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
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar          !< transformed flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: z1L,z1R,z2L,z2R,z3L,z3R,z4L,z4R,z5L,z5R
REAL                                :: sz1Mean,z1Mean,z4Mean,z5LN,z3Mean,z1LN,z2Mean,z5Mean
REAL                                :: rhoHat,uHat,vHat,wHat,p1Hat,rhoHat_HHat
REAL                                :: qHat
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5)  )

! z1 = √(rho/pressure) values on left and right, arithmatic and logorithmic means
z1L    = SQRT(rho_L/p_L)
z1R    = SQRT(rho_R/p_R)
z1Mean = 0.5*(z1L+z1R)
z1LN   = LN_Mean(z1L,z1R)
! z2 = u*√(rho/pressure) values on left and right and arithmatic mean
z2L    = z1L*VelU_L
z2R    = z1R*VelU_R
z2Mean = 0.5*(z2L+z2R)
! z3 = v*√(rho/pressure) values on left and right and arithmatic mean
z3L    = z1L*VelV_L
z3R    = z1R*VelV_R
z3Mean = 0.5*(z3L+z3R)
! z4 = w*√(rho/pressure) values on left and right and arithmatic mean
z4L    = z1L*VelW_L
z4R    = z1R*VelW_R
z4Mean = 0.5*(z4L+z4R)
! z5 = √(rho*pressure) values on left and right, arithmatic and logorithmic means
z5L    = z1L*p_L !SQRT(rho_L*p_L)
z5R    = z1R*p_R !SQRT(rho_R*p_R)
z5Mean = 0.5*(z5L+z5R)
z5LN   = LN_Mean(z5L,z5R)

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = z1Mean*z5LN
sz1Mean=1./z1Mean
uHat   = z2Mean*sz1Mean
vHat   = z3Mean*sz1Mean
wHat   = z4Mean*sz1Mean
p1Hat  = z5Mean*sz1Mean
rhoHat_HHat   = 0.5*(kappaP1*sKappaM1*(z5LN/z1LN)+p1Hat + rhoHat*(uHat*uHat + vHat*vHat + wHat*wHat)) !=rhoHat*HHat

qHat=uHat*metric(1)+vHat*metric(2)+wHat*metric(3)

! Entropy conserving flux
Fstar(1) = rhoHat*qHat
Fstar(2) = Fstar(1)*uHat + metric(1)*p1Hat
Fstar(3) = Fstar(1)*vHat + metric(2)*p1Hat
Fstar(4) = Fstar(1)*wHat + metric(3)*p1Hat
Fstar(5) = rhoHat_HHat*qHat

!! Entropy conserving flux
!Fstar(1) = rhoHat*uHat
!Fstar(2) = Fstar(1)*uHat + p1Hat
!Fstar(3) = Fstar(1)*vHat
!Fstar(4) = Fstar(1)*wHat
!Fstar(5) = uHat*rhoHat_HHat
!
!Gstar(1) = rhoHat*vHat
!Gstar(2) = Fstar(3)              !rhoHat*vHat*uHat 
!Gstar(3) = Gstar(1)*vHat + p1Hat
!Gstar(4) = Gstar(1)*wHat
!Gstar(5) = vHat*rhoHat_HHat
!
!Hstar(1) = rhoHat*wHat
!Hstar(2) = Fstar(4)              !rhoHat*wHat*uHat 
!Hstar(3) = Gstar(4)              !rhoHat*wHat*vHat
!Hstar(4) = Hstar(1)*wHat + p1Hat
!Hstar(5) = wHat*rhoHat_HHat

END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE TwoPointEntropyConservingFluxVec


!==================================================================================================================================
!> Computes the skew-symmetric 3D flux (in the normal direction) of Kennedy and Gruber for the Euler equations
!> Attention 1: Note that normal in this instance is always xHat, yHat, or zHat
!==================================================================================================================================
SUBROUTINE KennedyAndGruberFlux1(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: eHat
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,VelV_L,VelV_R,VelW_L,VelW_R
!==================================================================================================================================
! Get the inverse density, velocity, and pressure on left and right
sRho_L = 1./UL(1); VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1); VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

!p_L    = kappaM1*(UL(5) - 0.5*(UL(2)*VelU_L + UL(3)*VelV_L + UL(4)*VelW_L))
!p_R    = kappaM1*(UR(5) - 0.5*(UR(2)*VelU_R + UR(3)*VelV_R + UR(4)*VelW_R))
!e_L    = sRho_L*UL(5)
!e_R    = sRho_R*UR(5)

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = 0.5*(UL(1) + UR(1))
uHat   = 0.5*(VelU_L + VelU_R)
vHat   = 0.5*(VelV_L + VelV_R)
wHat   = 0.5*(VelW_L + VelW_R)
p1Hat  = 0.5*kappaM1*(UL(5)+UR(5)-0.5*( UL(2)*VelU_L + UL(3)*VelV_L + UL(4)*VelW_L  & 
                                       +UR(2)*VelU_R + UR(3)*VelV_R + UR(4)*VelW_R  )) !   0.5*(p_L + p_R)
eHat   = 0.5*(sRho_L*UL(5)+sRho_R*UR(5)) !0.5*(e_L + e_R)
aHat   = SQRT(kappa*p1Hat/rhoHat)
HHat   = kappa*p1Hat/(rhoHat*kappaM1) + 0.5*(uHat*uHat + vHat*vHat + wHat*wHat)
! Kennedy and Gruber skew-symmetric flux
Fstar(1) = rhoHat*uHat
Fstar(2) = Fstar(1)*uHat + p1Hat
Fstar(3) = Fstar(1)*vHat 
Fstar(4) = Fstar(1)*wHat
Fstar(5) = Fstar(1)*eHat + p1Hat*uHat
END SUBROUTINE KennedyAndGruberFlux1


!==================================================================================================================================
!> Computes the skew-symmetric flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the Euler equations
!> derived by Kennedy and Gruber
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE KennedyAndGruberFluxVec1(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
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
REAL                                :: rhoHat,uHat,vHat,wHat,p1Hat
REAL                                :: eHat,qHat
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
! Get the inverse density, velocity, and pressure on left and right
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
           srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5)  )

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = 0.5*( rho_L +  rho_R)
uHat   = 0.5*(VelU_L + VelU_R)
vHat   = 0.5*(VelV_L + VelV_R)
wHat   = 0.5*(VelW_L + VelW_R)
p1Hat  = 0.5*(p_L + p_R)
eHat   = 0.5*(sRho_L*rhoE_L+sRho_R*rhoE_R) !0.5*(e_L + e_R)

qHat=uHat*metric(1)+vHat*metric(2)+wHat*metric(3)

! Kennedy and Gruber skew-symmetric flux
Fstar(1) = rhoHat*qHat
Fstar(2) = Fstar(1)*uHat + metric(1)*p1Hat
Fstar(3) = Fstar(1)*vHat + metric(2)*p1Hat
Fstar(4) = Fstar(1)*wHat + metric(3)*p1Hat
Fstar(5) = Fstar(1)*eHat + p1Hat*qHat ! =\sum_i=1^3 (0.5*rhoHat*(eL+eR)*0.5*(u_iL+u_iR)+p1Hat*0.5*(u_iL+u_iR))*n_i

!! Kennedy and Gruber skew-symmetric flux
!Fstar(1) = rhoHat*uHat
!Fstar(2) = Fstar(1)*uHat + p1Hat
!Fstar(3) = Fstar(1)*vHat 
!Fstar(4) = Fstar(1)*wHat 
!Fstar(5) = Fstar(1)*eHat + p1Hat*uHat !=0.5*rhoHat*(eL+eR)*0.5*(uL+uR)+p1Hat*0.5*(uL+uR)
!
!Gstar(1) = rhoHat*vHat
!Gstar(2) = Fstar(3)               ! rhoHat*vHat*uHat 
!Gstar(3) = Gstar(1)*vHat + p1Hat
!Gstar(4) = Gstar(1)*wHat 
!Gstar(5) = Gstar(1)*eHat + p1Hat*vHat !=0.5*rhoHat*(eL+eR)*0.5*(vL+vR)+p1Hat*0.5*(vL+vR)
!
!Hstar(1) = rhoHat*wHat
!Hstar(2) = Fstar(4)               !rhoHat*wHat*uHat 
!Hstar(3) = Gstar(4)               !rhoHat*wHat*vHat
!Hstar(4) = Hstar(1)*wHat + p1Hat 
!Hstar(5) = Hstar(1)*eHat + p1Hat*wHat !=0.5*rhoHat*(eL+eR)*0.5*(wL+wR)+p1Hat*0.5*(wL+wR)

END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE KennedyAndGruberFluxVec1


!==================================================================================================================================
!> Computes Pirozzoli's misinterpretation of Kennedy and Gruber flux for the Euler equations
!> Attention 1: Note that normal in this instance is always xHat, yHat, or zHat
!==================================================================================================================================
SUBROUTINE KennedyAndGruberFlux2(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,p_L,p_R,VelV_L,VelV_R,VelW_L,VelW_R
!==================================================================================================================================
! Get the inverse density, velocity, and pressure on left and right
sRho_L = 1./UL(1); VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1); VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

p_L    = kappaM1*(UL(5) - 0.5*(UL(2)*VelU_L + UL(3)*VelV_L + UL(4)*VelW_L))
p_R    = kappaM1*(UR(5) - 0.5*(UR(2)*VelU_R + UR(3)*VelV_R + UR(4)*VelW_R))
!H_L    = (UL(5)+p_L)*sRho_L  
!H_R    = (UR(5)+p_R)*sRho_R

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = 0.5*(UL(1) + UR(1))
uHat   = 0.5*(VelU_L + VelU_R)
vHat   = 0.5*(VelV_L + VelV_R)
wHat   = 0.5*(VelW_L + VelW_R)
p1Hat  = 0.5*(   p_L +    p_R)
aHat   = SQRT(kappa*p1Hat/rhoHat)
HHat   = 0.5*( (UL(5)+p_L)*sRho_L + (UR(5)+p_R)*sRho_R ) !0.5*(H_L + H_R)
! Kennedy and Gruber skew-symmetric flux
Fstar(1) = rhoHat*uHat
Fstar(2) = Fstar(1)*uHat + p1Hat
Fstar(3) = Fstar(1)*vHat
Fstar(4) = Fstar(1)*wHat
Fstar(5) = Fstar(1)*HHat
END SUBROUTINE KennedyAndGruberFlux2


!==================================================================================================================================
!> Computes Pirozzoli's misinterpretation of Kennedy and Gruber flux for the Euler equations
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE KennedyAndGruberFluxVec2(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
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
REAL                                :: rhoHat,uHat,vHat,wHat,HHat,p1Hat,qHat
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
           srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5)  )

!H_L    = (rhoE_L+p_L)*sRho_L  
!H_R    = (rhoE_R+p_R)*sRho_R

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = 0.5*( rho_L +  rho_R)
uHat   = 0.5*(VelU_L + VelU_R)
vHat   = 0.5*(VelV_L + VelV_R)
wHat   = 0.5*(VelW_L + VelW_R)
p1Hat  = 0.5*(   p_L +    p_R)
HHat   = 0.5*((rhoE_L+p_L)*sRho_L + (rhoE_R+p_R)*sRho_R ) !0.5*(H_L + H_R)

qHat=uHat*metric(1)+vHat*metric(2)+wHat*metric(3)

! Kennedy and Gruber skew-symmetric flux
Fstar(1) = rhoHat*qHat
Fstar(2) = Fstar(1)*uHat + metric(1)*p1Hat
Fstar(3) = Fstar(1)*vHat + metric(2)*p1Hat
Fstar(4) = Fstar(1)*wHat + metric(3)*p1Hat
Fstar(5) = Fstar(1)*HHat

!! Kennedy and Gruber skew-symmetric flux
!Fstar(1) = rhoHat*uHat
!Fstar(2) = Fstar(1)*uHat + p1Hat
!Fstar(3) = Fstar(1)*vHat 
!Fstar(4) = Fstar(1)*wHat 
!Fstar(5) = Fstar(1)*HHat
!
!Gstar(1) = rhoHat*vHat
!Gstar(2) = Fstar(3)               !rhoHat*vHat*uHat 
!Gstar(3) = Gstar(1)*vHat + p1Hat
!Gstar(4) = Gstar(1)*wHat 
!Gstar(5) = Gstar(1)*HHat
!
!Hstar(1) = rhoHat*wHat
!Hstar(2) = Fstar(4)               !rhoHat*wHat*uHat 
!Hstar(3) = Gstar(4)               !rhoHat*wHat*vHat 
!Hstar(4) = Hstar(1)*wHat + p1Hat
!Hstar(5) = Hstar(1)*HHat
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE KennedyAndGruberFluxVec2


!==================================================================================================================================
!> Computes the skew-symmetric flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the Euler equations
!> derived by Ducros
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE DucrosFlux(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: qLR
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,VelV_L,VelV_R,VelW_L,VelW_R
!==================================================================================================================================
! Get the inverse density, velocity, and pressure on left and right

sRho_L = 1./UL(1); VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1); VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

!p_L    = kappaM1*(UL(5) - 0.5*(UL(2)*VelU_L + UL(3)*VelV_L + UL(4)*VelW_L))
!p_R    = kappaM1*(UR(5) - 0.5*(UR(2)*VelU_R + UR(3)*VelV_R + UR(4)*VelW_R))

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = 0.5*(UL(1) + UR(1))
uHat   = 0.5*(VelU_L + VelU_R)
vHat   = 0.5*(VelV_L + VelV_R)
wHat   = 0.5*(VelW_L + VelW_R)
p1Hat  = 0.5*kappaM1*(UL(5)+UR(5)-0.5*( UL(2)*VelU_L + UL(3)*VelV_L + UL(4)*VelW_L  & 
                                       +UR(2)*VelU_R + UR(3)*VelV_R + UR(4)*VelW_R  )) !   0.5*(p_L + p_R)
aHat   = SQRT(kappa*p1Hat/rhoHat)
HHat   = kappa*p1Hat/(rhoHat*kappaM1) + 0.5*(uHat*uHat + vHat*vHat + wHat*wHat)
! Ducros et. al. skew-symmetric flux
qLR      = 0.5*uHat
Fstar(1) = rhoHat*(qLR+qLR)                            !rhoHat*uHat
Fstar(2) = qLR*( UL(2) + UR(2)) + p1Hat
Fstar(3) = qLR*( UL(3) + UR(3)) 
Fstar(4) = qLR*( UL(4) + UR(4)) 
Fstar(5) = qLR*((UL(5) + UR(5)) + (p1Hat+p1Hat))
END SUBROUTINE DucrosFlux


!==================================================================================================================================
!> Computes the skew-symmetric 3D flux (in the normal direction) of Ducros' flux for the Euler equations
!> Attention 1: Note that normal in this instance is always xHat, yHat, or zHat
!==================================================================================================================================
SUBROUTINE DucrosFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
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
REAL                                :: uHat,vHat,wHat,rhoHat,p1Hat
REAL                                :: rhouHat,rhovHat,rhowHat,rhoEHat,qHat
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5)  )

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat  = 0.5*( rho_L +  rho_R)
uHat    = 0.5*(VelU_L + VelU_R)
vHat    = 0.5*(VelV_L + VelV_R)
wHat    = 0.5*(VelW_L + VelW_R)
rhouHat = 0.5*(rhoU_L + rhoU_R)
rhovHat = 0.5*(rhoV_L + rhoV_R)
rhowHat = 0.5*(rhoW_L + rhoW_R)
rhoeHat = 0.5*(rhoE_L + rhoE_R)
p1Hat   = 0.5*(   p_L +    p_R)

qHat=uHat*metric(1)+vHat*metric(2)+wHat*metric(3)

! Decros et. al. skew-symmetric flux
Fstar(1) =   rhoHat*qHat
Fstar(2) =  rhouHat*qHat +metric(1)*p1Hat
Fstar(3) =  rhovHat*qHat +metric(2)*p1Hat
Fstar(4) =  rhowHat*qHat +metric(3)*p1Hat 
Fstar(5) = (rhoEhat      + p1Hat)*qHat

!! Decros et. al. skew-symmetric flux
!Fstar(1) =   rhoHat*uHat            !=0.5*rhoHat*(velU_L+velU_R)
!Fstar(2) =  rhouHat*uHat + p1Hat
!Fstar(3) =  rhovHat*uHat 
!Fstar(4) =  rhowHat*uHat 
!Fstar(5) = (rhoEhat      + p1Hat)*uHat
!
!Gstar(1) =   rhoHat*vHat
!Gstar(2) =  rhouHat*vHat
!Gstar(3) =  rhovHat*vHat + p1Hat 
!Gstar(4) =  rhowHat*vHat 
!Gstar(5) = (rhoEhat      + p1Hat)*vHat
!
!Hstar(1) =   rhoHat*wHat
!Hstar(2) =  rhouHat*wHat
!Hstar(3) =  rhovHat*wHat
!Hstar(4) =  rhowHat*wHat + p1Hat  
!Hstar(5) = (rhoEhat      + p1Hat)*wHat
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE DucrosFluxVec


!==================================================================================================================================
!> Computes the skew-symmetric 3D flux (in the normal direction) of Morinishi flux for the Euler equations
!> Attention 1: Note that normal in this instance is always xHat, yHat, or zHat
!==================================================================================================================================
SUBROUTINE MorinishiFlux(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappa,kappaM1,sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: qL,qR,rhoqs2_L,rhoqs2_R,vel2_L,vel2_R
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,p_L,p_R,VelV_L,VelV_R,VelW_L,VelW_R
!==================================================================================================================================
! Get the inverse density, velocity, and pressure on left and right
sRho_L = 1./UL(1); VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1); VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

vel2_L = VelU_L*VelU_L+VelV_L*VelV_L+VelW_L*VelW_L
vel2_R = VelU_R*VelU_R+VelV_R*VelV_R+VelW_R*VelW_R
p_L    = kappaM1*(UL(5) - 0.5*UL(1)*vel2_L)
p_R    = kappaM1*(UR(5) - 0.5*UR(1)*vel2_R)

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = 0.5*(UL(1) + UR(1))
uHat   = 0.5*(VelU_L + VelU_R)
vHat   = 0.5*(VelV_L + VelV_R)
wHat   = 0.5*(VelW_L + VelW_R)
p1Hat  = 0.5*(p_L + p_R)
aHat   = SQRT(kappa*p1Hat/rhoHat)
HHat   = kappa*p1Hat/(rhoHat*kappaM1) + 0.5*(uHat*uHat + vHat*vHat + wHat*wHat)
! Morinishi kinetic energy preserving flux
qL       = VelU_L
qR       = VelU_R
rhoqs2_L = 0.5*UL(1)*qL
rhoqs2_R = 0.5*UR(1)*qR
Fstar(1) = (rhoqs2_L + rhoqs2_R)
Fstar(2) = Fstar(1)*uhat + p1Hat
Fstar(3) = Fstar(1)*vHat 
Fstar(4) = Fstar(1)*wHat 
Fstar(5) = 0.5*( kappa*skappaM1*(p_L*qL + p_R*qR) -(rhoqs2_L*vel2_L + rhoqs2_R*vel2_R) )&
           + rhoqs2_L*(VelU_L*uHat + VelV_L*vHat + VelW_L*wHat) &
           + rhoqs2_R*(VelU_R*uHat + VelV_R*vHat + VelW_R*wHat)
END SUBROUTINE MorinishiFlux


!==================================================================================================================================
!> Computes the skew-symmetric flux transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the Euler equations
!> derived by Morinishi
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE MorinishiFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
USE MOD_Equation_Vars,ONLY:kappa,skappaM1
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
REAL                                :: uHat,vHat,wHat,rhoHat,p1Hat
REAL                                :: q_L,q_R,rhoqs2_L,rhoqs2_R
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5), &
           Vel2_L =>UauxL(6), Vel2_R =>UauxR(6)  )

! Convenience variables for the velocity, density, pressure as well as sound speed and enthalpy
rhoHat = 0.5*(rho_L + rho_R)
uHat   = 0.5*(VelU_L + VelU_R)
vHat   = 0.5*(VelV_L + VelV_R)
wHat   = 0.5*(VelW_L + VelW_R)
p1Hat  = 0.5*(p_L + p_R)

q_L=VelU_L*metric(1)+VelV_L*metric(2)+VelW_L*metric(3)
q_R=VelU_R*metric(1)+VelV_R*metric(2)+VelW_R*metric(3)

! Morinishi kinetic energy preserving flux
rhoqs2_L = 0.5*rho_L*q_L
rhoqs2_R = 0.5*rho_R*q_R
Fstar(1) = (rhoqs2_L + rhoqs2_R)
Fstar(2) = Fstar(1)*uhat + metric(1)*p1Hat
Fstar(3) = Fstar(1)*vHat + metric(2)*p1Hat
Fstar(4) = Fstar(1)*wHat + metric(3)*p1Hat
Fstar(5) = 0.5*( kappa*skappaM1*(p_L*q_L + p_R*q_R) -(rhoqs2_L*vel2_L + rhoqs2_R*vel2_R) ) &
           + rhoqs2_L*(VelU_L*uHat + VelV_L*vHat + VelW_L*wHat) &
           + rhoqs2_R*(VelU_R*uHat + VelV_R*vHat + VelW_R*wHat)

!! Morinishi kinetic energy preserving flux
!rhoqs2_L = 0.5*rho_L*VelU_L
!rhoqs2_R = 0.5*rho_R*VelU_R
!Fstar(1) = (rhoqs2_L + rhoqs2_R)
!Fstar(2) = Fstar(1)*uhat + p1Hat
!Fstar(3) = Fstar(1)*vHat 
!Fstar(4) = Fstar(1)*wHat 
!Fstar(5) = 0.5*( kappa*skappaM1*(p_L*velU_L + p_R*velU_R) -(rhoqs2_L*vel2_L + rhoqs2_R*vel2_R) ) &
!           + rhoqs2_L*vvhat_L+ rhoqs2_R*vvhat_R
!
!rhoqs2_L = 0.5*rho_L*VelV_L
!rhoqs2_R = 0.5*rho_R*VelV_R
!Gstar(1) = (rhoqs2_L + rhoqs2_R)
!Gstar(2) = Gstar(1)*uhat
!Gstar(3) = Gstar(1)*vHat + p1Hat 
!Gstar(4) = Gstar(1)*wHat 
!Gstar(5) = 0.5*( kappa*skappaM1*(p_L*velV_L + p_R*velV_R) -(rhoqs2_L*vel2_L + rhoqs2_R*vel2_R) ) &
!           + rhoqs2_L*vvhat_L+ rhoqs2_R*vvhat_R
!
!rhoqs2_L = 0.5*rho_L*VelW_L
!rhoqs2_R = 0.5*rho_R*VelW_R
!Hstar(1) = (rhoqs2_L + rhoqs2_R)
!Hstar(2) = Hstar(1)*uhat
!Hstar(3) = Hstar(1)*vHat 
!Hstar(4) = Hstar(1)*wHat + p1Hat 
!Hstar(5) = 0.5*( kappa*skappaM1*(p_L*velW_L + p_R*velW_R) -(rhoqs2_L*vel2_L + rhoqs2_R*vel2_R) ) &
!           + rhoqs2_L*vvhat_L+ rhoqs2_R*vvhat_R
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE MorinishiFluxVec


!==================================================================================================================================
!> Computes the entropy and kinetic energy conserving numerical 3D flux (in the
!> normal direction) for the Euler equations
!> Attention 1: Note that normal in this instance is always xHat, yHat, or zHat
!> TODO: could make it like EvalFlux3D and do every direction simutaneously
!==================================================================================================================================
SUBROUTINE EntropyAndEnergyConservingFlux(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,skappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: qHat
REAL                                :: rho_MEAN,beta_MEAN,beta_Hat
REAL                                :: beta_R,beta_L
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,p_L,p_R,VelV_L,VelV_R,VelW_L,VelW_R
REAL                                :: Vel2s2_L,Vel2s2_R,Vel2_M
!==================================================================================================================================

! Get the inverse density, velocity, and pressure on left and right
sRho_L = 1./UL(1); VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1); VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

Vel2s2_L = 0.5*(VelU_L*VelU_L+VelV_L*VelV_L+VelW_L*VelW_L)
Vel2s2_R = 0.5*(VelU_R*VelU_R+VelV_R*VelV_R+VelW_R*VelW_R)
p_L    = kappaM1*(UL(5) - UL(1)*Vel2s2_L)
p_R    = kappaM1*(UR(5) - UR(1)*Vel2s2_R)
beta_L = 0.5*UL(1)/p_L
beta_R = 0.5*UR(1)/p_R
aHat   = 0.
HHat   = 0.

! Get the averages for the numerical flux

rho_MEAN  = 0.5*(   UL(1)+UR(1))
rhoHat    = LN_MEAN(UL(1),UR(1))
beta_MEAN = 0.5*(    beta_L+beta_R)
beta_Hat  = LN_MEAN(beta_L,beta_R)
uHat      = 0.5*(    VelU_L+VelU_R)
vHat      = 0.5*(    VelV_L+VelV_R)
wHat      = 0.5*(    VelW_L+VelW_R)
p1Hat     = 0.5*rho_MEAN/beta_MEAN
Vel2_M    = Vel2s2_L+Vel2s2_R

! Entropy conserving and kinetic energy conserving flux
qHat     = uHat
Fstar(1) = rhoHat*qHat
Fstar(2) = Fstar(1)*uHat + p1Hat
Fstar(3) = Fstar(1)*vHat
Fstar(4) = Fstar(1)*wHat
Fstar(5) = Fstar(1)*0.5*(skappaM1/beta_Hat - Vel2_M)  &
           + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)
END SUBROUTINE EntropyAndEnergyConservingFlux


!==================================================================================================================================
!> Computes the entropy and kinetic energy conserving flux transformed with the metrics
!> (fstar=f*metric1+g*metric2+h*metric3 ) for the Euler equations
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE EntropyAndEnergyConservingFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
USE MOD_Equation_Vars,ONLY:sKappaM1
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
REAL                                :: uHat,vHat,wHat,rhoHat,p1Hat
REAL                                :: rho_MEAN,beta_MEAN,beta_Hat,qHat
REAL                                :: beta_R,beta_L
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5), &
           Vel2_L =>UauxL(6), Vel2_R =>UauxR(6)  )
beta_L = 0.5*rho_L/p_L
beta_R = 0.5*rho_R/p_R

! Get the averages for the numerical flux

rho_MEAN  = 0.5*(   rho_L+rho_R)
rhoHat    = LN_MEAN(rho_L,rho_R)
beta_MEAN = 0.5*(    beta_L+beta_R)
beta_Hat  = LN_MEAN(beta_L,beta_R)
uHat      = 0.5*(    VelU_L+VelU_R)
vHat      = 0.5*(    VelV_L+VelV_R)
wHat      = 0.5*(    VelW_L+VelW_R)
p1Hat     = 0.5*rho_MEAN/beta_MEAN

qHat=uHat*metric(1)+vHat*metric(2)+wHat*metric(3)

! Entropy conserving and kinetic energy conserving flux
Fstar(1) = rhoHat*qHat
Fstar(2) = Fstar(1)*uHat + metric(1)*p1Hat
Fstar(3) = Fstar(1)*vHat + metric(2)*p1Hat
Fstar(4) = Fstar(1)*wHat + metric(3)*p1Hat 
Fstar(5) = Fstar(1)*0.5*(skappaM1/beta_Hat - 0.5*(Vel2_L+Vel2_R)) &
           + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)

!! Entropy conserving and kinetic energy conserving flux
!Fstar(1) = rhoHat*uHat
!Fstar(2) = Fstar(1)*uHat + p1Hat
!Fstar(3) = Fstar(1)*vHat 
!Fstar(4) = Fstar(1)*wHat 
!Fstar(5) = Fstar(1)*HHat + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)
!
!Gstar(1) = rhoHat*vHat
!Gstar(2) = Fstar(3)              !rhoHat*vHat*uHat 
!Gstar(3) = Gstar(1)*vHat + p1Hat
!Gstar(4) = Gstar(1)*wHat 
!Gstar(5) = Gstar(1)*HHat + uHat*Gstar(2) + vHat*Gstar(3) + wHat*Gstar(4)
!
!Hstar(1) = rhoHat*wHat
!Hstar(2) = Fstar(4)              !rhoHat*wHat*uHat 
!Hstar(3) = Gstar(4)              !rhoHat*wHat*vHat 
!Hstar(4) = Hstar(1)*wHat + p1Hat
!Hstar(5) = Hstar(1)*HHat + uHat*Hstar(2) + vHat*Hstar(3) + wHat*Hstar(4)
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE EntropyAndEnergyConservingFluxVec


!==================================================================================================================================
!> Computes the entropy and kinetic energy conserving numerical 3D flux (in the
!> normal direction) for the Euler equations
!> Attention 1: Note that normal in this instance is always xHat, yHat, or zHat
!> TODO: could make it like EvalFlux3D and do every direction simutaneously
!==================================================================================================================================
SUBROUTINE EntropyAndEnergyConservingFlux2(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: qHat
REAL                                :: rho_MEAN,sbeta_MEAN
REAL                                :: beta_R,beta_L
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,p_L,p_R,VelV_L,VelV_R,VelW_L,VelW_R
REAL                                :: Vel2s2_L,Vel2s2_R,Vel2_M
!==================================================================================================================================

! Get the inverse density, velocity, and pressure on left and right
sRho_L = 1./UL(1) ; VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1) ; VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

Vel2s2_L=0.5*(VelU_L*VelU_L+VelV_L*VelV_L+VelW_L*VelW_L)
Vel2s2_R=0.5*(VelU_R*VelU_R+VelV_R*VelV_R+VelW_R*VelW_R)
p_L    = kappaM1*(UL(5) - UL(1)*Vel2s2_L)
p_R    = kappaM1*(UR(5) - UR(1)*Vel2s2_R)
beta_L = 0.5*UL(1)/p_L
beta_R = 0.5*UR(1)/p_R
aHat   = 0.
HHat   = 0.

! Get the averages for the numerical flux

rho_MEAN   = 0.5*(   UL(1)+UR(1))
rhoHat     = LN_MEAN(UL(1),UR(1))
sbeta_MEAN = 2./(     beta_L+beta_R)
uHat       = 0.5*(    VelU_L+VelU_R)
vHat       = 0.5*(    VelV_L+VelV_R)
wHat       = 0.5*(    VelW_L+VelW_R)
p1Hat      = 0.5*rho_MEAN*sbeta_MEAN
Vel2_M     = Vel2s2_L+Vel2s2_R !0.5*(v_L^2 + v_R^2)

! Entropy conserving and kinetic energy conserving flux
qHat     = uHat
Fstar(1) = rho_MEAN*qHat
Fstar(2) = Fstar(1)*uHat + p1Hat
Fstar(3) = Fstar(1)*vHat
Fstar(4) = Fstar(1)*wHat
Fstar(5) = Fstar(1)*0.5*(skappaM1*sbeta_MEAN - Vel2_M) &
           + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)
END SUBROUTINE EntropyAndEnergyConservingFlux2


!==================================================================================================================================
!> Computes the entropy and kinetic energy conserving flux transformed with the metrics
!> (fstar=f*metric1+g*metric2+h*metric3 ) for the Euler equations
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE EntropyAndEnergyConservingFluxVec2(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
USE MOD_Equation_Vars,ONLY:sKappaM1
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
REAL                                :: uHat,vHat,wHat,p1Hat,qHat
REAL                                :: rho_MEAN,sbeta_MEAN
REAL                                :: beta_R,beta_L
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5), &
           Vel2_L =>UauxL(6), Vel2_R =>UauxR(6)  )
beta_L = 0.5*rho_L/p_L
beta_R = 0.5*rho_R/p_R

! Get the averages for the numerical flux

rho_MEAN   = 0.5*(   rho_L+ rho_R)
sbeta_MEAN = 2./(   beta_L+beta_R)
uHat       = 0.5*(  VelU_L+VelU_R)
vHat       = 0.5*(  VelV_L+VelV_R)
wHat       = 0.5*(  VelW_L+VelW_R)
p1Hat      = 0.5*rho_MEAN*sbeta_MEAN

qHat=uHat*metric(1)+vHat*metric(2)+wHat*metric(3)

! Entropy conserving and kinetic energy conserving flux
Fstar(1) = rho_MEAN*qHat
Fstar(2) = Fstar(1)*uHat + metric(1)*p1Hat
Fstar(3) = Fstar(1)*vHat + metric(2)*p1Hat
Fstar(4) = Fstar(1)*wHat + metric(3)*p1Hat
Fstar(5) = Fstar(1)*0.5*(skappaM1*sbeta_MEAN - 0.5*(Vel2_L+Vel2_R)) &
           + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)

!! Entropy conserving and kinetic energy conserving flux
!Fstar(1) = rho_MEAN*uHat
!Fstar(2) = Fstar(1)*uHat + p1Hat
!Fstar(3) = Fstar(1)*vHat
!Fstar(4) = Fstar(1)*wHat
!Fstar(5) = Fstar(1)*HHat + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)
!
!Gstar(1) = rho_MEAN*vHat
!Gstar(2) = Fstar(3)               !rho_MEAN*vHat*uHat
!Gstar(3) = Gstar(1)*vHat + p1Hat
!Gstar(4) = Gstar(1)*wHat
!Gstar(5) = Gstar(1)*HHat + uHat*Gstar(2) + vHat*Gstar(3) + wHat*Gstar(4)
!
!Hstar(1) = rho_MEAN*wHat
!Hstar(2) = Fstar(4)               !rho_MEAN*wHat*uHat
!Hstar(3) = Gstar(4)               !rho_MEAN*wHat*vHat
!Hstar(4) = Hstar(1)*wHat + p1Hat
!Hstar(5) = Hstar(1)*HHat + uHat*Hstar(2) + vHat*Hstar(3) + wHat*Hstar(4)
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE EntropyAndEnergyConservingFluxVec2


!==================================================================================================================================
!> Computes the entropy and kinetic energy conserving numerical 3D flux (in the
!> normal direction) for the Euler equations
!> Attention 1: Note that normal in this instance is always xHat, yHat, or zHat
!> TODO: could make it like EvalFlux3D and do every direction simutaneously
!==================================================================================================================================
SUBROUTINE ggflux(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,skappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: qHat
REAL                                :: rho_MEAN,beta_Hat
REAL                                :: beta_R,beta_L
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,p_L,p_R,VelV_L,VelV_R,VelW_L,VelW_R
REAL                                :: Vel2s2_L,Vel2s2_R,Vel2_M
!==================================================================================================================================

! Get the inverse density, velocity, and pressure on left and right
sRho_L = 1./UL(1); VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1); VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

Vel2s2_L=0.5*(VelU_L*VelU_L+VelV_L*VelV_L+VelW_L*VelW_L)
Vel2s2_R=0.5*(VelU_R*VelU_R+VelV_R*VelV_R+VelW_R*VelW_R)
p_L    = (kappaM1)*(UL(5) - UL(1)*Vel2s2_L)
p_R    = (kappaM1)*(UR(5) - UR(1)*Vel2s2_R)
beta_L = 0.5*UL(1)/p_L
beta_R = 0.5*UR(1)/p_R
aHat   = 0.
HHat   = 0.

! Get the averages for the numerical flux

rho_MEAN   = 0.5*(UL(1)+UR(1))
rhoHat     = LN_MEAN(UL(1),UR(1))
beta_Hat   = LN_MEAN(beta_L,beta_R)
uHat       = 0.5*(VelU_L+VelU_R)
vHat       = 0.5*(VelV_L+VelV_R)
wHat       = 0.5*(VelW_L+VelW_R)
p1Hat      = 0.5*(p_L+p_R)!rho_MEAN/beta_MEAN !1./!!!
Vel2_M     = Vel2s2_L+Vel2s2_R

! Entropy conserving and kinetic energy conserving flux
qHat     = uHat
Fstar(1) = rhoHat*qHat
Fstar(2) = Fstar(1)*uHat + p1Hat
Fstar(3) = Fstar(1)*vHat
Fstar(4) = Fstar(1)*wHat
Fstar(5) = Fstar(1)*0.5*(skappaM1/beta_Hat - Vel2_M) &
           + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)
END SUBROUTINE ggFlux


!==================================================================================================================================
!> Computes the entropy and kinetic energy conserving flux transformed with the metrics
!> (fstar=f*metric1+g*metric2+h*metric3 ) for the Euler equations
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE ggfluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
USE MOD_Equation_Vars,ONLY:skappaM1
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
REAL                                :: uHat,vHat,wHat,rhoHat,p1Hat,qHat
REAL                                :: rho_MEAN,beta_Hat
REAL                                :: beta_R,beta_L
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
          !srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5), &
           Vel2_L =>UauxL(6), Vel2_R =>UauxR(6)  )
beta_L = 0.5*rho_L/p_L
beta_R = 0.5*rho_R/p_R

! Get the averages for the numerical flux

rho_MEAN   = 0.5*(rho_L+rho_R)
rhoHat     = LN_MEAN(rho_L,rho_R)
beta_Hat   = LN_MEAN(beta_L,beta_R)
uHat       = 0.5*(VelU_L+VelU_R)
vHat       = 0.5*(VelV_L+VelV_R)
wHat       = 0.5*(VelW_L+VelW_R)
p1Hat      = 0.5*(p_L+p_R)

qHat=uHat*metric(1)+vHat*metric(2)+wHat*metric(3)

! Entropy conserving and kinetic energy conserving flux
Fstar(1) = rho_MEAN*qHat
Fstar(2) = Fstar(1)*uHat + metric(1)*p1Hat
Fstar(3) = Fstar(1)*vHat + metric(2)*p1Hat
Fstar(4) = Fstar(1)*wHat + metric(3)*p1Hat
Fstar(5) = Fstar(1)*0.5*(skappaM1/beta_Hat -0.5*(Vel2_L+Vel2_R))  &
           + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)

!! Entropy conserving and kinetic energy conserving flux
!Fstar(1) = rhoHat*uHat
!Fstar(2) = Fstar(1)*uHat + p1Hat
!Fstar(3) = Fstar(1)*vHat
!Fstar(4) = Fstar(1)*wHat
!Fstar(5) = Fstar(1)*HHat + uHat*Fstar(2) + vHat*Fstar(3) + wHat*Fstar(4)
!
!Gstar(1) = rhoHat*vHat
!Gstar(2) = Fstar(3)                !rhoHat*vHat*uHat
!Gstar(3) = Gstar(1)*vHat + p1Hat
!Gstar(4) = Gstar(1)*wHat
!Gstar(5) = Gstar(1)*HHat + uHat*Gstar(2) + vHat*Gstar(3) + wHat*Gstar(4)
!
!Hstar(1) = rhoHat*wHat
!Hstar(2) = Fstar(4)                !rhoHat*wHat*uHat
!Hstar(3) = Gstar(4)                !rhoHat*wHat*vHat
!Hstar(4) = Hstar(1)*wHat + p1Hat
!Hstar(5) = Hstar(1)*HHat + uHat*Hstar(2) + vHat*Hstar(3) + wHat*Hstar(4)
END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE ggFluxVec


!==================================================================================================================================
!> Computes the GW^2 flux for the Rho variables split form
!> Attention 1: Note that normal in this instance is always xHat, yHat, or zHat
!> TODO: could make it like EvalFlux3D and do every direction simutaneously
!==================================================================================================================================
SUBROUTINE GassnerWintersWalchFlux(Fstar,UL,UR,uHat,vHat,wHat,aHat,HHat,p1Hat,rhoHat)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:kappaM1,sKappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UL      !< left state
REAL,DIMENSION(PP_nVar),INTENT(IN)  :: UR      !< right state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: Fstar   !<  flux in x
REAL                   ,INTENT(OUT) :: uHat,vHat,wHat,aHat,rhoHat,HHat,p1Hat !additional variables for riemann
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: w1L,w1R,w2L,w2R,w3L,w3R,w4L,w4R,w5L,w5R,qHat,h_L,h_R
REAL                                :: sRho_L,sRho_R,VelU_L,VelU_R,p_L,p_R,VelV_L,VelV_R,VelW_L,VelW_R
!==================================================================================================================================

! Get the inverse density, velocity, and pressure on left and right
sRho_L = 1./UL(1); VelU_L = sRho_L*UL(2) ;  VelV_L = sRho_L*UL(3) ; VelW_L = sRho_L*UL(4)
sRho_R = 1./UR(1); VelU_R = sRho_R*UR(2) ;  VelV_R = sRho_R*UR(3) ; VelW_R = sRho_R*UR(4)

p_L    = kappaM1*(UL(5) - 0.5*(UL(2)*VelU_L + UL(3)*VelV_L + UL(4)*VelW_L))
p_R    = kappaM1*(UR(5) - 0.5*(UR(2)*VelU_R + UR(3)*VelV_R + UR(4)*VelW_R))
h_L    = sRho_L*(UL(5)+p_L) ! kappa*UL(1)*p_L/kappaM1 + 0.5*(VelU_L*VelU_L + VelV_L*VelV_L + VelW_L*VelW_L)
h_R    = sRho_R*(UR(5)+p_R) ! kappa*UR(1)*p_R/kappaM1 + 0.5*(VelU_R*VelU_R + VelV_R*VelV_R + VelW_R*VelW_R)
! get the Rho variables on both sides
w1L = SQRT(UL(1))
w1R = SQRT(UR(1))
w2L = w1L*VelU_L
w2R = w1R*VelU_R
w3L = w1L*VelV_L
w3R = w1R*VelV_R
w4L = w1L*VelW_L
w4R = w1R*VelW_R
w5L = w1L*h_L
w5R = w1R*h_R
! compute the average of each Rho variable
rhoHat = 0.5*(w1L+w1R) ! w1 average
uHat   = 0.5*(w2L+w2R) ! w2 average
vHat   = 0.5*(w3L+w3R) ! w3 average
wHat   = 0.5*(w4L+w4R) ! w4 average
Hhat   = 0.5*(w5L+w5R) ! w5 average
p1Hat  = kappaM1*skappa*(rhoHat*Hhat - 0.5*(uHat*uHat + vHat*vHat + wHat*wHat))
aHat   = 0.
! normal "velocity"
qHat   = uHat
! compute the Rho variable split form flux
Fstar(1) = qHat*rhoHat
Fstar(2) = qHat*uHat + p1Hat
Fstar(3) = qHat*vHat
Fstar(4) = qHat*wHat
Fstar(5) = qhat*Hhat

END SUBROUTINE GassnerWintersWalchFlux


!==================================================================================================================================
!> Computes the GW^2 flux for the Rho variables split form
!> transformed with the metrics (fstar=f*metric1+g*metric2+h*metric3 ) for the Euler equations
!> for curved metrics, 1/2(metric_L+metric_R) is taken!
!==================================================================================================================================
SUBROUTINE GassnerWintersWalchFluxVec(UL,UR,UauxL,UauxR,metric_L,metric_R,Fstar)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:nAuxVar
USE MOD_Equation_Vars,ONLY:kappaM1,sKappa
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
REAL                                :: uHat,vHat,wHat,rhoHat,p1Hat,Hhat,qHat
REAL                                :: sqrtRhoL,sqrtRhoR
REAL                                :: metric(3)
!==================================================================================================================================
metric = 0.5*(metric_L+metric_R)
ASSOCIATE(  rho_L =>   UL(1),  rho_R =>   UR(1), &
           rhoU_L =>   UL(2), rhoU_R =>   UR(2), &
           rhoV_L =>   UL(3), rhoV_R =>   UR(3), &
           rhoW_L =>   UL(4), rhoW_R =>   UR(4), &
           rhoE_L =>   UL(5), rhoE_R =>   UR(5), &
           srho_L =>UauxL(1), srho_R =>UauxR(1), &
           VelU_L =>UauxL(2), VelU_R =>UauxR(2), &
           VelV_L =>UauxL(3), VelV_R =>UauxR(3), &
           VelW_L =>UauxL(4), VelW_R =>UauxR(4), &
              p_L =>UauxL(5),    p_R =>UauxR(5)  )
! get the Rho variables on both sides
sqrtRhoL = SQRT(rho_L)
sqrtRhoR = SQRT(rho_R)
! compute the average of each Rho variable
rhoHat = 0.5*(sqrtRhoL        + sqrtRhoR) ! w1 average
uHat   = 0.5*(sqrtRhoL*VelU_L + sqrtRhoR*VelU_R) ! w2 average
vHat   = 0.5*(sqrtRhoL*VelV_L + sqrtRhoR*VelV_R) ! w3 average
wHat   = 0.5*(sqrtRhoL*VelW_L + sqrtRhoR*VelW_R) ! w4 average
Hhat   = 0.5*( sqrtRhoL*sRho_L*(rhoE_L+p_L) &
              +sqrtRhoR*sRho_R*(rhoE_R+p_R)) ! w1L*HL average
p1Hat  = kappaM1*skappa*(rhoHat*Hhat - 0.5*(uHat*uHat + vHat*vHat + wHat*wHat))

qHat=uHat*metric(1)+vHat*metric(2)+wHat*metric(3)

! compute the Rho variable split form flux
Fstar(1) = rhoHat*qHat
Fstar(2) = qHat*uHat + metric(1)*p1Hat
Fstar(3) = qHat*vHat + metric(2)*p1Hat
Fstar(4) = qHat*wHat + metric(3)*p1Hat
Fstar(5) = qhat*Hhat

!! compute the Rho variable split form flux
!Fstar(1) = rhoHat*uHat
!Fstar(2) = uHat*uHat + p1Hat
!Fstar(3) = uHat*vHat
!Fstar(4) = uHat*wHat
!Fstar(5) = uhat*Hhat
!
!Gstar(1) = rhoHat*vHat
!Gstar(2) = Fstar(3) !vHat*uHat
!Gstar(3) = vHat*vHat + p1Hat
!Gstar(4) = vHat*wHat 
!Gstar(5) = vhat*Hhat
!
!Hstar(1) = rhoHat*wHat
!Hstar(2) = Fstar(4) !wHat*uHat
!Hstar(3) = Gstar(4) !wHat*vHat
!Hstar(4) = wHat*wHat + p1Hat
!Hstar(5) = what*Hhat

END ASSOCIATE !rho_L/R,rhov1_L/R,...
END SUBROUTINE GassnerWintersWalchFluxVec

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
REAL FUNCTION LN_MEAN(aL,aR)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL         :: aL  !< left value
REAL         :: aR  !< right value
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
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
