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
!> Contains the routine EvalFluxTilde3D which computes the complete tranformed MHD flux f,g,h for all DOFs in one Element,
!> the routine EvalDiffFluxTilde3D which computes only the diffusive part of the tranformed MHD flux f,g,h  ,
!> the routine EvalFlux1D_Adv which computes the Advection flux f for all DOFs of one Element side: used in Riemann_Adv
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

!INTERFACE EvalFluxTilde3D
!  MODULE PROCEDURE EvalFluxTilde3D
!END INTERFACE

INTERFACE EvalAdvectionFlux1D
  MODULE PROCEDURE EvalAdvectionFlux1D
END INTERFACE

#if PARABOLIC
! no interface because of dimension change
!INTERFACE EvalDiffFlux3D
!  MODULE PROCEDURE EvalDiffFlux3D
!END INTERFACE

!INTERFACE EvalDiffFluxTilde3D
!  MODULE PROCEDURE EvalDiffFluxTilde3D
!END INTERFACE

!INTERFACE EvalLiftingVolumeFlux
!  MODULE PROCEDURE EvalLiftingVolumeFlux
!END INTERFACE

INTERFACE EvalLiftingSurfFlux
  MODULE PROCEDURE EvalLiftingSurfFlux
END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::EvalAdvFluxTilde3D
PUBLIC::EvalFluxTilde3D
PUBLIC::EvalAdvectionFlux1D
#if PARABOLIC
PUBLIC::EvalDiffFlux3D
PUBLIC::EvalDiffFluxTilde3D
PUBLIC::EvalLiftingVolumeFlux
PUBLIC::EvalLiftingSurfFlux
#endif /*PARABOLIC*/

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute MHD transformed advection fluxes using conservative variables and derivatives for every volume Gauss point.
!> directly apply metrics and output the tranformed flux 
!==================================================================================================================================
PURE SUBROUTINE EvalAdvFluxTilde3D(U_in,M_f,M_g,M_h,ftilde,gtilde,htilde)
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
REAL,INTENT(IN )   :: U_in(PP_nVar,1:nTotal_vol) !< solution state (conservative vars)
REAL,INTENT(IN )   :: M_f(       3,1:nTotal_vol) !< metrics for ftilde                 
REAL,INTENT(IN )   :: M_g(       3,1:nTotal_vol) !< metrics for gtilde                 
REAL,INTENT(IN )   :: M_h(       3,1:nTotal_vol) !< metrics for htilde                 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(PP_nVar,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(PP_nVar,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(PP_nVar,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
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
  ASSOCIATE(rho   =>U_in(IRHO1,i), &
            rhov1 =>U_in(IRHOU,i), &
            rhov2 =>U_in(IRHOV,i), &
            rhov3 =>U_in(IRHOW,i), &
#ifdef PP_GLM
            Etotal=>U_in(IRHOE,i)-0.5*smu_0*U_in(IPSI,i)**2, &
#else
            Etotal=>U_in(IRHOE,i), &
#endif /*def PP_GLM*/
            b1    =>U_in(IB1,i), &
            b2    =>U_in(IB2,i), &
            b3    =>U_in(IB3,i)  ) 
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
  
  ! Advection part
  ! Advection fluxes x-direction
  f(IRHO1)=rhov1                     ! rho*u
  f(IRHOU)=rhov1*v1+pt -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
  f(IRHOV)=rhov1*v2    -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
  f(IRHOW)=rhov1*v3    -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
  f(IRHOE)=Ep*v1       -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
  f(IB1)=0.
  f(IB2)=v1*b2-b1*v2
  f(IB3)=v1*b3-b1*v3
  ! Advection fluxes y-direction
  g(IRHO1)=rhov2                     ! rho*v      
  g(IRHOU)=f(IRHOV)                      ! rho*u*v      -1/mu_0*b2*b1
  g(IRHOV)=rhov2*v2+pt -smu_0*b2*b2  ! rho*v²+p     -1/mu_0*b2*b2
  g(IRHOW)=rhov2*v3    -smu_0*b2*b3  ! rho*v*w      -1/mu_0*b2*b3
  g(IRHOE)=Ep*v2       -smu_0*b2*vb  ! (rho*e+p)*v  -1/mu_0*b2*(v dot B)
  g(IB1)=-f(IB2)                     ! (v2*b1-b2*v1)
  g(IB2)=0.
  g(IB3)=v2*b3-b2*v3
  ! Advection fluxes z-direction
  h(IRHO1)=rhov3                     ! rho*v
  h(IRHOU)=f(IRHOW)                      ! rho*u*w      -1/mu_0*b3*b1
  h(IRHOV)=g(IRHOW)                      ! rho*v*w      -1/mu_0*b3*b2
  h(IRHOW)=rhov3*v3+pt -smu_0*b3*b3  ! rho*v²+p     -1/mu_0*b3*b3
  h(IRHOE)=Ep*v3       -smu_0*b3*vb  ! (rho*e+p)*w  -1/mu_0*b3*(v dot B)
  h(IB1)=-f(IB3)                     ! v3*b1-b3*v1 
  h(IB2)=-g(IB3)                     ! v3*b2-b3*v2
  h(IB3)=0.

#ifdef PP_GLM
  f(IRHOE) = f(IRHOE)+smu_0*GLM_ch*b1*U_in(IPSI,i)
  f(IB1) = f(IB1)+GLM_ch*U_in(IPSI,i)
  f(IPSI) = GLM_ch*b1

  g(IRHOE) = g(IRHOE)+smu_0*GLM_ch*b2*U_in(IPSI,i)
  g(IB2) = g(IB2)+GLM_ch*U_in(IPSI,i)
  g(IPSI) = GLM_ch*b2

  h(IRHOE) = h(IRHOE)+smu_0*GLM_ch*b3*U_in(IPSI,i)
  h(IB3) = h(IB3)+GLM_ch*U_in(IPSI,i)
  h(IPSI) = GLM_ch*b3
#endif /* PP_GLM */

  END ASSOCIATE ! rho,rhov1,rhov2,rhov3,Etotal,b1,b2,b3

  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(:,i) =   f(:)*M_f(1,i) + g(:)*M_f(2,i) + h(:)*M_f(3,i)
  gtilde(:,i) =   f(:)*M_g(1,i) + g(:)*M_g(2,i) + h(:)*M_g(3,i)
  htilde(:,i) =   f(:)*M_h(1,i) + g(:)*M_h(2,i) + h(:)*M_h(3,i)
END DO ! i
END SUBROUTINE EvalAdvFluxTilde3D

!==================================================================================================================================
!> Compute transformed MHD fluxes using the conservative variables and gradients for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFluxTilde3D(U_in,M_f,M_g,M_h, &
#if PARABOLIC
                           gradPx_in,gradPy_in,gradPz_in,&
#endif /*PARABOLIC*/
                           ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars      ,ONLY:kappaM1,smu_0,s2mu_0,sKappaM1
#ifdef PP_GLM
USE MOD_Equation_Vars      ,ONLY: GLM_ch
#endif /*PP_GLM*/
#if PARABOLIC
USE MOD_Equation_Vars      ,ONLY:mu,s23,etasmu_0
#ifdef PP_ANISO_HEAT
USE MOD_Equation_vars      ,ONLY:kperp,kpar
#else
USE MOD_Equation_vars      ,ONLY:kappasPr
#endif 
#endif /*PARABOLIC*/
USE MOD_DG_Vars            ,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN )   :: U_in(PP_nVar,1:nTotal_vol) !< state in conservative variables
REAL,INTENT(IN )   :: M_f(       3,1:nTotal_vol) !< metrics for ftilde                 
REAL,INTENT(IN )   :: M_g(       3,1:nTotal_vol) !< metrics for gtilde                 
REAL,INTENT(IN )   :: M_h(       3,1:nTotal_vol) !< metrics for htilde                 
#if PARABOLIC
REAL,INTENT(IN )   :: gradPx_in(PP_nVar,1:nTotal_vol) !< gradient x in primitive variables 
REAL,INTENT(IN )   :: gradPy_in(PP_nVar,1:nTotal_vol) !< gradient y in primitive variables 
REAL,INTENT(IN )   :: gradPz_in(PP_nVar,1:nTotal_vol) !< gradient z in primitive variables 
#endif /*PARABOLIC*/
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(PP_nVar,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(PP_nVar,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(PP_nVar,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:PP_nVar) :: f,g,h                             ! Cartesian fluxes (iVar)
REAL                :: srho                                    ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,p,ptilde                       ! velocity and pressure(including magnetic pressure
REAL                :: bb2,vb                                  ! magnetic field, bb2=|bvec|^2, v dot b
REAL                :: Ep                                      ! E + p
#if PARABOLIC
REAL                :: divv
REAL                :: lambda 
REAL                :: cv_gradTx,cv_gradTy,cv_gradTz
REAL                :: Qx,Qy,Qz
REAL,DIMENSION(IRHOU:PP_nVar) :: f_visc,g_visc,h_visc
#endif /*PARABOLIC*/
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  ASSOCIATE(rho   =>U_in(IRHO1,i), &
            rhov1 =>U_in(IRHOU,i), &
            rhov2 =>U_in(IRHOV,i), &
            rhov3 =>U_in(IRHOW,i), &
            b1    =>U_in(IB1,i), &
            b2    =>U_in(IB2,i), &
            b3    =>U_in(IB3,i), &
#ifdef PP_GLM
            Etotal=>U_in(IRHOE,i)-s2mu_0*(U_in(IPSI,i)**2)  &
#else
            Etotal=>U_in(IRHOE,i)  &
#endif /*def PP_GLM*/
            )
  ! auxiliary variables
  srho = 1. / rho ! 1/rho
  v1   = rhov1*srho 
  v2   = rhov2*srho 
  v3   = rhov3*srho 
  bb2  = (b1*b1+b2*b2+b3*b3)
  vb   = (b1*v1+b2*v2+b3*v3)
  !gas pressure
  p    = kappaM1*(Etotal - 0.5*(rho*(v1*v1+v2*v2+v3*v3)) - s2mu_0*bb2)
  ! ptilde includes magnetic pressure
  ptilde = p + s2mu_0*bb2
  Ep   = (Etotal + ptilde)
  ! Advection part
  ! Advection fluxes x-direction
  f(IRHO1)=rhov1                          ! rho*u
  f(IRHOU)=rhov1*v1+ptilde  -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
  f(IRHOV)=rhov1*v2         -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
  f(IRHOW)=rhov1*v3         -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
  f(IRHOE)=Ep*v1            -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
  f(IB1)=0.
  f(IB2)=v1*b2-b1*v2
  f(IB3)=v1*b3-b1*v3
  ! Advection fluxes y-direction
  g(IRHO1)=rhov2                          ! rho*v      
  g(IRHOU)=f(IRHOV)                           ! rho*u*v      -1/mu_0*b2*b1
  g(IRHOV)=rhov2*v2+ptilde  -smu_0*b2*b2  ! rho*v²+p     -1/mu_0*b2*b2
  g(IRHOW)=rhov2*v3         -smu_0*b2*b3  ! rho*v*w      -1/mu_0*b2*b3
  g(IRHOE)=Ep*v2            -smu_0*b2*vb  ! (rho*e+p)*v  -1/mu_0*b2*(v dot B)
  g(IB1)=-f(IB2)                          ! v2*b1-b2*v1 
  g(IB2)=0.
  g(IB3)=v2*b3-b2*v3
  ! Advection fluxes z-direction
  h(IRHO1)=rhov3                          ! rho*v
  h(IRHOU)=f(IRHOW)                           ! rho*u*w      -1/mu_0*b3*b1
  h(IRHOV)=g(IRHOW)                           ! rho*v*w      -1/mu_0*b3*b2
  h(IRHOW)=rhov3*v3+ptilde  -smu_0*b3*b3  ! rho*v²+p     -1/mu_0*b3*b3
  h(IRHOE)=Ep*v3            -smu_0*b3*vb  ! (rho*e+p)*w  -1/mu_0*b3*(v dot B)
  h(IB1)=-f(IB3)                          ! v3*b1-b3*v1 
  h(IB2)=-g(IB3)                          ! v3*b2-b3*v2
  h(IB3)=0.

#ifdef PP_GLM
  f(IRHOE) = f(IRHOE)+smu_0*GLM_ch*b1*U_in(IPSI,i)
  f(IB1) = f(IB1)+      GLM_ch   *U_in(IPSI,i)
  f(IPSI) =            GLM_ch*b1

  g(IRHOE) = g(IRHOE)+smu_0*GLM_ch*b2*U_in(IPSI,i)
  g(IB2) = g(IB2)+      GLM_ch   *U_in(IPSI,i)
  g(IPSI) =            GLM_ch*b2

  h(IRHOE) = h(IRHOE)+smu_0*GLM_ch*b3*U_in(IPSI,i)
  h(IB3) = h(IB3)+      GLM_ch   *U_in(IPSI,i)
  h(IPSI) =            GLM_ch*b3
#endif /* PP_GLM */

#if PARABOLIC
  ! Viscous part
  ASSOCIATE( gradv1x => gradPx_in(IU,i), gradB1x => gradPx_in(IB1,i), & 
             gradv2x => gradPx_in(IV,i), gradB2x => gradPx_in(IB2,i), & 
             gradv3x => gradPx_in(IW,i), gradB3x => gradPx_in(IB3,i), & 
             gradv1y => gradPy_in(IU,i), gradB1y => gradPy_in(IB1,i), & 
             gradv2y => gradPy_in(IV,i), gradB2y => gradPy_in(IB2,i), & 
             gradv3y => gradPy_in(IW,i), gradB3y => gradPy_in(IB3,i), & 
             gradv1z => gradPz_in(IU,i), gradB1z => gradPz_in(IB1,i), & 
             gradv2z => gradPz_in(IV,i), gradB2z => gradPz_in(IB2,i), & 
             gradv3z => gradPz_in(IW,i), gradB3z => gradPz_in(IB3,i)  )

  divv    = gradv1x+gradv2y+gradv3z
  cv_gradTx  = sKappaM1*sRho*(gradPx_in(IP,i)-srho*p*gradPx_in(IRHO1,i))  ! cv*T_x = 1/(kappa-1) *1/rho *(p_x - p/rho*rho_x)
  cv_gradTy  = sKappaM1*sRho*(gradPy_in(IP,i)-srho*p*gradPy_in(IRHO1,i)) 
  cv_gradTz  = sKappaM1*sRho*(gradPz_in(IP,i)-srho*p*gradPz_in(IRHO1,i)) 

#ifndef PP_ANISO_HEAT
  !isotropic heat flux
  lambda=mu*KappasPr
  Qx=lambda*cv_gradTx  !q=lambda*gradT= (mu*kappa/Pr)*(cv*gradT)
  Qy=lambda*cv_gradTy
  Qz=lambda*cv_gradTz
#else
  !IF(bb2.GT. 0.)THEN! ATTENTION: |B|^2 /= 0 !!!
  lambda=(kpar-kperp)/bb2  
  Qx=kappaM1*(lambda*(b1*b1*cv_gradTx+b1*b2*cv_gradTy+b1*b3*cv_gradTz)+kperp*cv_gradTx)
  Qy=kappaM1*(lambda*(b2*b1*cv_gradTx+b2*b2*cv_gradTy+b2*b3*cv_gradTz)+kperp*cv_gradTy)
  Qz=kappaM1*(lambda*(b3*b1*cv_gradTx+b3*b2*cv_gradTy+b3*b3*cv_gradTz)+kperp*cv_gradTz)
  !ELSE
  !  Qx=kperp*cv_gradTx
  !  Qy=kperp*cv_gradTy
  !  Qz=kperp*cv_gradTz
  !END IF 
#endif /*PP_ANISO_HEAT*/
  ! viscous fluxes in x-direction      
  f_visc(IRHOU)=mu*(2*gradv1x-s23*divv)
  f_visc(IRHOV)=mu*(  gradv2x+gradv1y)   
  f_visc(IRHOW)=mu*(  gradv3x+gradv1z)   
  f_visc(IB1)=0.
  f_visc(IB2)=etasmu_0*(gradB2x-gradB1y)
  f_visc(IB3)=etasmu_0*(gradB3x-gradB1z)
  !energy
  f_visc(IRHOE)= v1*f_visc(IRHOU)+v2*f_visc(IRHOV)+v3*f_visc(IRHOW) +smu_0*(b1*f_visc(IB1)+b2*f_visc(IB2)+b3*f_visc(IB3)) + Qx

  ! viscous fluxes in y-direction      
  g_visc(IRHOU)= f_visc(IRHOV)                  !mu*(  gradv1y+gradv2x)  
  g_visc(IRHOV)=mu*(2*gradv2y-s23*divv)     
  g_visc(IRHOW)=mu*(  gradv3y+gradv2z)      
  g_visc(IB1)=-f_visc(IB2)                  !etasmu_0*(gradB1y-gradB2x)
  g_visc(IB2)=0.
  g_visc(IB3)=etasmu_0*(gradB3y-gradB2z)
  !energy
  g_visc(IRHOE)=v1*g_visc(IRHOU)+v2*g_visc(IRHOV)+v3*g_visc(IRHOW) + smu_0*(b1*g_visc(IB1)+b2*g_visc(IB2)+b3*g_visc(IB3)) + Qy 



  ! viscous fluxes in z-direction      
  h_visc(IRHOU)= f_visc(IRHOW)                       !mu*(  gradv1z+gradv3x)                 
  h_visc(IRHOV)= g_visc(IRHOW)                       !mu*(  gradv2z+gradv3y)                
  h_visc(IRHOW)=mu*(2*gradv3z-s23*divv )             
  h_visc(IB1)=-f_visc(IB3)                       !etasmu_0*(gradB1z-gradB3x)
  h_visc(IB2)=-g_visc(IB3)                       !etasmu_0*(gradB2z-gradB3y)
  h_visc(IB3)=0.
  !energy
  h_visc(IRHOE)=v1*h_visc(IRHOU)+v2*h_visc(IRHOV)+v3*h_visc(IRHOW)+smu_0*(b1*h_visc(IB1)+b2*h_visc(IB2)+b3*h_visc(IB3)) + Qz

  f(IRHOU:IB3)=f(IRHOU:IB3)-f_visc(IRHOU:IB3)
  g(IRHOU:IB3)=g(IRHOU:IB3)-g_visc(IRHOU:IB3)
  h(IRHOU:IB3)=h(IRHOU:IB3)-h_visc(IRHOU:IB3)

END ASSOCIATE ! gradB1x => gradPx(IB1 ...

#endif /*PARABOLIC*/
END ASSOCIATE ! rho,rhov1,rhov2,rhov3,Etotal,b1,b2,b3

  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(:,i) =   f(:)*M_f(1,i) + g(:)*M_f(2,i) + h(:)*M_f(3,i)
  gtilde(:,i) =   f(:)*M_g(1,i) + g(:)*M_g(2,i) + h(:)*M_g(3,i)
  htilde(:,i) =   f(:)*M_h(1,i) + g(:)*M_h(2,i) + h(:)*M_h(3,i)
END DO ! i
END SUBROUTINE EvalFluxTilde3D


!==================================================================================================================================
!> Computes the 1D Advection fluxes using the conservative variables for every volume Gauss point.
!==================================================================================================================================
pure SUBROUTINE EvalAdvectionFlux1D(U_Face,F_Face)
! MODULES
USE MOD_Equation_Vars,ONLY:KappaM1,KappaM2,smu_0,s2mu_0
#ifdef PP_GLM
USE MOD_Equation_vars,ONLY:GLM_ch
#endif 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face(PP_nVar)    !< input conservative state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: F_Face(PP_nVar)    !< F_Face(iVar), Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: srho               ! reciprocal values for density and kappa/Pr
REAL                :: v1,v2,v3,ptilde    ! velocity and pressure(including magnetic pressure
REAL                :: vb                 ! magnetic field, bb2=|bvec|^2
!==================================================================================================================================
! auxiliary variables
srho = 1. / U_Face(IRHO1) ! 1/rho
v1   = U_Face(IRHOU)*srho ! u
v2   = U_Face(IRHOV)*srho ! v
v3   = U_Face(IRHOW)*srho ! w
ASSOCIATE( b1 => U_Face(IB1), &
           b2 => U_Face(IB2), &
           b3 => U_Face(IB3), &
#ifdef PP_GLM
           Etotal => U_Face(IRHOE)-s2mu_0*U_Face(IPSI)**2  &
#else
           Etotal => U_Face(IRHOE)  &
#endif
)
vb   = (b1*v1+b2*v2+b3*v3)
! ptilde includes magnetic pressure
ptilde = kappaM1*(Etotal-0.5*U_Face(IRHO1)*(v1*v1+v2*v2+v3*v3))-kappaM2*s2mu_0*(b1*b1+b2*b2+b3*b3)
! Advection fluxes x-direction
F_Face(IRHO1)= U_Face(IRHOU)          ! rho*u
F_Face(IRHOU)= U_Face(IRHOU)*v1+ptilde    -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
F_Face(IRHOV)= U_Face(IRHOU)*v2           -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
F_Face(IRHOW)= U_Face(IRHOU)*v3           -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
F_Face(IRHOE)=(Etotal + ptilde)*v1    -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
F_Face(IB1)=0.
F_Face(IB2)=v1*b2-b1*v2
F_Face(IB3)=v1*b3-b1*v3
#ifdef PP_GLM
F_Face(IRHOE)=F_Face(IRHOE)+smu_0*GLM_ch*b1*U_Face(IPSI)
F_Face(IB1)=F_Face(IB1)+      GLM_ch   *U_Face(IPSI)
F_Face(IPSI)=                GLM_ch*b1
#endif /* PP_GLM */
END ASSOCIATE ! b1 => UFace(IB1) ...
END SUBROUTINE EvalAdvectionFlux1D


#if PARABOLIC
!==================================================================================================================================
!> Compute the cartesian diffusion flux of the MHD equations, for all face points (PP_N+1)**2
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D(f,g,h,U_Face,gradPx_Face,gradPy_Face,gradPz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:s23,smu_0,s2mu_0,kappaM1,sKappaM1
USE MOD_Equation_Vars,ONLY:mu,etasmu_0
USE MOD_DG_Vars,ONLY:nTotal_face
#ifdef PP_ANISO_HEAT
USE MOD_Equation_vars,ONLY:kperp,kpar
#else
USE MOD_Equation_vars,ONLY:kappasPr
#endif /*PP_ANISO_HEAT*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face(     PP_nVar,nTotal_face)    !< state in conservative variables on surface points 
REAL,INTENT(IN)     :: gradPx_Face(PP_nVar,nTotal_face)    !< x gradient of state 
REAL,INTENT(IN)     :: gradPy_Face(PP_nVar,nTotal_face)    !< y gradient of state 
REAL,INTENT(IN)     :: gradPz_Face(PP_nVar,nTotal_face)    !< z gradient of state 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: f       !< Cartesian diffusion flux in x
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: g       !< Cartesian diffusion flux in y
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: h       !< Cartesian diffusion flux in z
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: srho,p
REAL                :: v1,v2,v3                                  ! auxiliary variables
REAL                :: bb2
REAL                :: divv
REAL                :: lambda
REAL                :: cv_gradTx,cv_gradTy,cv_gradTz
REAL                :: Qx,Qy,Qz
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_face
  ! auxiliary variables
  ASSOCIATE(rho   =>U_Face(IRHO1,i), &
            b1    =>U_Face(IB1,i), &
            b2    =>U_Face(IB2,i), &
            b3    =>U_Face(IB3,i)  )
  srho = 1. / rho
  v1   = U_Face(IRHOU,i)*srho
  v2   = U_Face(IRHOV,i)*srho
  v3   = U_Face(IRHOW,i)*srho
  p    = kappaM1*(U_Face(IRHOE,i) - 0.5*(rho*(v1*v1+v2*v2+v3*v3)) - s2mu_0*SUM(U_Face(IB1:PP_nVar,i)**2) )

  ! Viscous part
  ! ideal gas law
ASSOCIATE( gradv1x => gradPx_Face(IU,i),  gradB1x => gradPx_Face(IB1,i), & 
           gradv2x => gradPx_Face(IV,i),  gradB2x => gradPx_Face(IB2,i), & 
           gradv3x => gradPx_Face(IW,i),  gradB3x => gradPx_Face(IB3,i), & 
           gradv1y => gradPy_Face(IU,i),  gradB1y => gradPy_Face(IB1,i), & 
           gradv2y => gradPy_Face(IV,i),  gradB2y => gradPy_Face(IB2,i), & 
           gradv3y => gradPy_Face(IW,i),  gradB3y => gradPy_Face(IB3,i), & 
           gradv1z => gradPz_Face(IU,i),  gradB1z => gradPz_Face(IB1,i), & 
           gradv2z => gradPz_Face(IV,i),  gradB2z => gradPz_Face(IB2,i), & 
           gradv3z => gradPz_Face(IW,i),  gradB3z => gradPz_Face(IB3,i)  )

  divv    = gradv1x+gradv2y+gradv3z
  cv_gradTx  = sKappaM1*sRho*(gradPx_face(IP,i)-srho*p*gradPx_face(IRHO1,i))  ! cv*T_x = 1/(kappa-1) *1/rho *(p_x - p/rho*rho_x)
  cv_gradTy  = sKappaM1*sRho*(gradPy_face(IP,i)-srho*p*gradPy_face(IRHO1,i)) 
  cv_gradTz  = sKappaM1*sRho*(gradPz_face(IP,i)-srho*p*gradPz_face(IRHO1,i)) 
#ifndef PP_ANISO_HEAT
  !isotropic heat flux
  lambda=mu*KappasPr

  Qx=lambda*cv_gradTx  !q=lambda*gradT= (mu*kappa/Pr)*(cv*gradT)
  Qy=lambda*cv_gradTy
  Qz=lambda*cv_gradTz
#else
  bb2  = (b1*b1+b2*b2+b3*b3)
  !IF(bb2.GT. 0.)THEN! ATTENTION: |B|^2 /= 0 !!!
  lambda=(kpar-kperp)/bb2  
  Qx=kappaM1*(lambda*(b1*b1*cv_gradTx+b1*b2*cv_gradTy+b1*b3*cv_gradTz)+kperp*cv_gradTx)
  Qy=kappaM1*(lambda*(b2*b1*cv_gradTx+b2*b2*cv_gradTy+b2*b3*cv_gradTz)+kperp*cv_gradTy)
  Qz=kappaM1*(lambda*(b3*b1*cv_gradTx+b3*b2*cv_gradTy+b3*b3*cv_gradTz)+kperp*cv_gradTz)
  !ELSE
  !  Qx=kperp*cv_gradTx
  !  Qy=kperp*cv_gradTy
  !  Qz=kperp*cv_gradTz
  !END IF 
#endif /*PP_ANISO_HEAT*/
  ! viscous fluxes in x-direction      
  f(IRHO1,i)=0.
  f(IRHOU,i)=-mu*(2*gradv1x-s23*divv)
  f(IRHOV,i)=-mu*(  gradv2x+gradv1y)   
  f(IRHOW,i)=-mu*(  gradv3x+gradv1z)   
  f(IB1,i)= 0.
  f(IB2,i)=-etasmu_0*(gradB2x-gradB1y)
  f(IB3,i)=-etasmu_0*(gradB3x-gradB1z)
  !energy
  f(IRHOE,i)=v1*f(IRHOU,i)+v2*f(IRHOV,i)+v3*f(IRHOW,i) +smu_0*(b1*f(IB1,i)+b2*f(IB2,i)+b3*f(IB3,i)) - Qx

  ! viscous fluxes in y-direction      
  g(IRHO1,i)=0.
  g(IRHOU,i)= f(IRHOV,i)                  !mu*(  gradv1y+gradv2x)  
  g(IRHOV,i)=-mu*(2*gradv2y-s23*divv)
  g(IRHOW,i)=-mu*(  gradv3y+gradv2z)      
  g(IB1,i)=-f(IB2,i)                  !etasmu_0*(gradB1y-gradB2x)
  g(IB2,i)= 0.
  g(IB3,i)=-etasmu_0*(gradB3y-gradB2z)
  !energy
  g(IRHOE,i)=v1*g(IRHOU,i)+v2*g(IRHOV,i)+v3*g(IRHOW,i) + smu_0*(b1*g(IB1,i)+b2*g(IB2,i)+b3*g(IB3,i)) - Qy

  ! viscous fluxes in z-direction      
  h(IRHO1,i)=0.
  h(IRHOU,i)= f(IRHOW,i)                       !mu*(  gradv1z+gradv3x)                 
  h(IRHOV,i)= g(IRHOW,i)                       !mu*(  gradv2z+gradv3y)                
  h(IRHOW,i)=-mu*(2*gradv3z-s23*divv )        
  h(IB1,i)=-f(IB3,i)                       !etasmu_0*(gradB1z-gradB3x)
  h(IB2,i)=-g(IB3,i)                       !etasmu_0*(gradB2z-gradB3y)
  h(IB3,i)= 0.
  !energy
  h(IRHOE,i)=v1*h(IRHOU,i)+v2*h(IRHOV,i)+v3*h(IRHOW,i) +smu_0*(b1*h(IB1,i)+b2*h(IB2,i)+b3*h(IB3,i)) - Qz
  
#ifdef PP_GLM
  f(IPSI,i) = 0.
  g(IPSI,i) = 0.
  h(IPSI,i) = 0.
#endif /*PP_GLM*/

END ASSOCIATE ! gradv1x => gradPx(2,....
END ASSOCIATE ! b1 => UFace(IB1,i) ...
END DO !i=1,nTotal
END SUBROUTINE EvalDiffFlux3D


!==================================================================================================================================
!> Compute only diffusive part of tranformed MHD fluxes, using conservative variables and derivatives for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFluxTilde3D(U_in,M_f,M_g,M_h,gradPx_in,gradPy_in,gradPz_in,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars      ,ONLY:smu_0,mu,s23,etasmu_0,s2mu_0,kappaM1,sKappaM1
#ifdef PP_ANISO_HEAT
USE MOD_Equation_vars      ,ONLY:kperp,kpar
#else
USE MOD_Equation_vars      ,ONLY:kappasPr
#endif 
USE MOD_DG_Vars            ,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN )   :: U_in(     PP_nVar,1:nTotal_vol) !< state in conservative variables
REAL,INTENT(IN )   :: M_f(            3,1:nTotal_vol) !< metrics for ftilde                 
REAL,INTENT(IN )   :: M_g(            3,1:nTotal_vol) !< metrics for gtilde                 
REAL,INTENT(IN )   :: M_h(            3,1:nTotal_vol) !< metrics for htilde                 
REAL,INTENT(IN )   :: gradPx_in(PP_nVar,1:nTotal_vol) !< gradient x in primitive variables 
REAL,INTENT(IN )   :: gradPy_in(PP_nVar,1:nTotal_vol) !< gradient y in primitive variables 
REAL,INTENT(IN )   :: gradPz_in(PP_nVar,1:nTotal_vol) !< gradient z in primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(PP_nVar,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(PP_nVar,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(PP_nVar,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:PP_nVar) :: f_Visc,g_visc,h_visc              ! Cartesian fluxes (iVar)
REAL                :: srho                                    ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,bb2,p 
REAL                :: divv
REAL                :: lambda 
REAL                :: cv_gradTx,cv_gradTy,cv_gradTz
REAL                :: Qx,Qy,Qz
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  ! auxiliary variables
ASSOCIATE( rho => U_in(IRHO1,i), &
           b1  => U_in(IB1,i), &
           b2  => U_in(IB2,i), &
           b3  => U_in(IB3,i)  )
  srho = 1. / rho
  
  v1   = U_in(IRHOU,i)*srho 
  v2   = U_in(IRHOV,i)*srho 
  v3   = U_in(IRHOW,i)*srho 

  p    = kappaM1*(U_in(IRHOE,i) - 0.5*(rho*(v1*v1+v2*v2+v3*v3)) - s2mu_0*SUM(U_in(IB1:PP_nVar,i)**2) )
  ! Viscous part
ASSOCIATE( gradv1x => gradPx_in(IU,i), gradB1x => gradPx_in(IB1,i), & 
           gradv2x => gradPx_in(IV,i), gradB2x => gradPx_in(IB2,i), & 
           gradv3x => gradPx_in(IW,i), gradB3x => gradPx_in(IB3,i), & 
           gradv1y => gradPy_in(IU,i), gradB1y => gradPy_in(IB1,i), & 
           gradv2y => gradPy_in(IV,i), gradB2y => gradPy_in(IB2,i), & 
           gradv3y => gradPy_in(IW,i), gradB3y => gradPy_in(IB3,i), & 
           gradv1z => gradPz_in(IU,i), gradB1z => gradPz_in(IB1,i), & 
           gradv2z => gradPz_in(IV,i), gradB2z => gradPz_in(IB2,i), & 
           gradv3z => gradPz_in(IW,i), gradB3z => gradPz_in(IB3,i))

  divv    = gradv1x+gradv2y+gradv3z
  cv_gradTx  = sKappaM1*sRho*(gradPx_in(IP,i)-srho*p*gradPx_in(IRHO1,i))  ! cv*T_x = 1/(kappa-1) *1/rho *(p_x - p/rho*rho_x)
  cv_gradTy  = sKappaM1*sRho*(gradPy_in(IP,i)-srho*p*gradPy_in(IRHO1,i)) 
  cv_gradTz  = sKappaM1*sRho*(gradPz_in(IP,i)-srho*p*gradPz_in(IRHO1,i)) 
#ifndef PP_ANISO_HEAT
  !isotropic heat flux
  
  lambda=mu*KappasPr
  
  Qx=lambda*cv_gradTx  !q=lambda*gradT= (mu*kappa/Pr)*(cv*gradT)
  Qy=lambda*cv_gradTy
  Qz=lambda*cv_gradTz
#else
  bb2  = (b1*b1+b2*b2+b3*b3)
  !IF(bb2.GT. 0.)THEN! ATTENTION: |B|^2 /= 0 !!!
  lambda=(kpar-kperp)/bb2  
  Qx=kappaM1*(lambda*(b1*b1*cv_gradTx+b1*b2*cv_gradTy+b1*b3*cv_gradTz)+kperp*cv_gradTx)
  Qy=kappaM1*(lambda*(b2*b1*cv_gradTx+b2*b2*cv_gradTy+b2*b3*cv_gradTz)+kperp*cv_gradTy)
  Qz=kappaM1*(lambda*(b3*b1*cv_gradTx+b3*b2*cv_gradTy+b3*b3*cv_gradTz)+kperp*cv_gradTz)
  !ELSE
  !  Qx=kperp*cv_gradTx
  !  Qy=kperp*cv_gradTy
  !  Qz=kperp*cv_gradTz
  !END IF 
#endif /*PP_ANISO_HEAT*/
  ! viscous fluxes in x-direction      
  f_visc(IRHO1)=0.
  f_visc(IRHOU)=-mu*(2*gradv1x-s23*divv)
  f_visc(IRHOV)=-mu*(  gradv2x+gradv1y)   
  f_visc(IRHOW)=-mu*(  gradv3x+gradv1z)   
  f_visc(IB1)=0.
  f_visc(IB2)=-etasmu_0*(gradB2x-gradB1y)
  f_visc(IB3)=-etasmu_0*(gradB3x-gradB1z)
  !energy
  f_visc(IRHOE)= v1*f_visc(IRHOU)+v2*f_visc(IRHOV)+v3*f_visc(IRHOW) +smu_0*(b1*f_visc(IB1)+b2*f_visc(IB2)+b3*f_visc(IB3)) - Qx

  ! viscous fluxes in y-direction      
  g_visc(IRHO1)=0.
  g_visc(IRHOU)= f_visc(IRHOV)                  !-mu*(  gradv1y+gradv2x)  
  g_visc(IRHOV)=-mu*(2*gradv2y-s23*divv)     
  g_visc(IRHOW)=-mu*(  gradv3y+gradv2z)      
  g_visc(IB1)=-f_visc(IB2)                  !etasmu_0*(gradB1y-gradB2x)
  g_visc(IB2)=0.
  g_visc(IB3)=-etasmu_0*(gradB3y-gradB2z)
  !energy
  g_visc(IRHOE)=v1*g_visc(IRHOU)+v2*g_visc(IRHOV)+v3*g_visc(IRHOW) + smu_0*(b1*g_visc(IB1)+b2*g_visc(IB2)+b3*g_visc(IB3)) - Qy 

  ! viscous fluxes in z-direction      
  h_visc(IRHO1)=0.
  h_visc(IRHOU)= f_visc(IRHOW)                       !-mu*(  gradv1z+gradv3x)                 
  h_visc(IRHOV)= g_visc(IRHOW)                       !-mu*(  gradv2z+gradv3y)                
  h_visc(IRHOW)=-mu*(2*gradv3z-s23*divv )             
  h_visc(IB1)=-f_visc(IB3)                       !etasmu_0*(gradB1z-gradB3x)
  h_visc(IB2)=-g_visc(IB3)                       !etasmu_0*(gradB2z-gradB3y)
  h_visc(IB3)=0.
  !energy
  h_visc(IRHOE)=v1*h_visc(IRHOU)+v2*h_visc(IRHOV)+v3*h_visc(IRHOW)+smu_0*(b1*h_visc(IB1)+b2*h_visc(IB2)+b3*h_visc(IB3)) - Qz

#ifdef PP_GLM
  f_visc(IPSI) = 0. 
  g_visc(IPSI) = 0. 
  h_visc(IPSI) = 0. 
#endif /*PP_GLM*/

END ASSOCIATE ! gradB1x => gradPx(IB1 ...

END ASSOCIATE ! b1 = U(IB1,....

  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(:,i) =   f_visc(:)*M_f(1,i) + g_visc(:)*M_f(2,i) + h_visc(:)*M_f(3,i)
  gtilde(:,i) =   f_visc(:)*M_g(1,i) + g_visc(:)*M_g(2,i) + h_visc(:)*M_g(3,i)
  htilde(:,i) =   f_visc(:)*M_h(1,i) + g_visc(:)*M_h(2,i) + h_visc(:)*M_h(3,i)
END DO ! i
END SUBROUTINE EvalDiffFluxTilde3D


!==================================================================================================================================
!> Compute the lifting flux depending on the variable to be used for the gradient: cons_var / prim_var / entropy_var 
!==================================================================================================================================
SUBROUTINE EvalLiftingVolumeFlux(U_in,Flux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,ONLY:nTotal_vol
#if PP_Lifting_Var==2
USE MOD_Equation_Vars,ONLY: ConsToPrimVec
#elif PP_Lifting_Var==3
USE MOD_Equation_Vars,ONLY: ConsToEntropyVec
#endif /*PP_Lifting_Var**/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN )   :: U_in(PP_nVar,1:nTotal_vol) !< state in conservative variables
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: flux(PP_nVar,1:nTotal_vol) !< lifting flux, depending on lifting_var
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if PP_Lifting_Var==1
  Flux=U_in
#elif PP_Lifting_Var==2
  CALL ConsToPrimVec(nTotal_vol,Flux,U_in) !prim_var
#elif PP_Lifting_Var==3
  CALL ConsToEntropyVec(nTotal_vol,Flux,U_in) !entropy_var
#endif /*PP_Lifting_Var**/

END SUBROUTINE EvalLiftingVolumeFlux


!==================================================================================================================================
!> Compute the average lifting surface flux in strong form U*=1/2(U_s+U_m)  (BR1/BR2),
!> careful, we use the strong form, so that the surface flux becomes: F=U^*-U_m = 1/2(U_s-U_m)
!> depending on the variable to be used for the gradient: 
!> cons_var / prim_var / entropy_var 
!==================================================================================================================================
SUBROUTINE EvalLiftingSurfFlux(SideID,Flux)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars   ,ONLY:nTotal_face
USE MOD_DG_Vars   ,ONLY: U_master,U_slave
USE MOD_Mesh_Vars ,ONLY: SurfElem
#if (PP_Lifting_Var==2)
USE MOD_Equation_Vars,ONLY: ConsToPrimVec
#elif (PP_Lifting_Var==3)
USE MOD_Equation_Vars,ONLY: ConsToEntropyVec
#endif /*PP_Lifting_Var**/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)     :: SideID       !< current side index
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)       :: flux(PP_nVar,0:PP_N,0:PP_N) !< lifting surface flux, depending on lifting_var
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER  :: p,q 
REAL     :: F_m(PP_nVar,0:PP_N,0:PP_N)
REAL     :: F_s(PP_nVar,0:PP_N,0:PP_N)
!==================================================================================================================================
#if (PP_Lifting_Var==1)
  F_m=U_master(:,:,:,SideID)
  F_s=U_slave(:,:,:,SideID)
#elif (PP_Lifting_Var==2)
  CALL ConsToPrimVec(nTotal_face,F_m,U_master(:,:,:,SideID)) !prim_var
  CALL ConsToPrimVec(nTotal_face,F_s,U_slave( :,:,:,SideID)) !prim_var
#elif (PP_Lifting_Var==3)
  CALL ConsToEntropyVec(nTotal_face,F_m,U_master(:,:,:,SideID)) !prim_var
  CALL ConsToEntropyVec(nTotal_face,F_s,U_slave( :,:,:,SideID)) !prim_var
#endif /*PP_Lifting_Var**/
  !strong lifting flux: 
  DO q=0,PP_N; DO p=0,PP_N
      Flux(:,p,q)=0.5*(F_s(:,p,q)-F_m(:,p,q))*SurfElem(p,q,SideID)
  END DO; END DO

END SUBROUTINE EvalLiftingSurfFlux

#endif /*PARABOLIC*/

END MODULE MOD_Flux
