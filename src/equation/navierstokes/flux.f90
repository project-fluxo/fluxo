!==================================================================================================================================
! Copyright (c) 2016 - 2017 Gregor Gassner
! Copyright (c) 2016 - 2017 Florian Hindenlang
! Copyright (c) 2016 - 2017 Andrew Winters
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

!==================================================================================================================================
!> Contains the routine EvalFlux3D which computes the complete flux f,g,h for all DOFs in one Element: used in volume integral
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE EvalFluxTilde3D
  MODULE PROCEDURE EvalFluxTilde3D
END INTERFACE

INTERFACE EvalEulerFlux1D
  MODULE PROCEDURE EvalEulerFlux1D
END INTERFACE


#if PARABOLIC
INTERFACE EvalDiffFlux1D_Outflow
  MODULE PROCEDURE EvalDiffFlux1D_Outflow
END INTERFACE
INTERFACE EvalDiffFlux3D
  MODULE PROCEDURE EvalDiffFlux3D
END INTERFACE

INTERFACE EvalDiffFluxTilde3D
  MODULE PROCEDURE EvalDiffFluxTilde3D
END INTERFACE
#endif /*PARABOLIC*/


PUBLIC::EvalFluxTilde3D
PUBLIC::EvalEulerFlux1D
#if PARABOLIC
PUBLIC::EvalDiffFlux1D_Outflow
PUBLIC::EvalDiffFlux3D
PUBLIC::EvalDiffFluxTilde3D
#endif /*PARABOLIC*/
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute transformed 3D Navier-Stokes fluxes(Euler+diffusion) for every volume Gauss point of element iElem.
!> In comparison to EvalFlux3D, metrics are directly applied 
!==================================================================================================================================
SUBROUTINE EvalFluxTilde3D(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Equation_Vars ,ONLY:kappaM1,R
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
#if PARABOLIC
USE MOD_Lifting_Vars  ,ONLY:gradUx,gradUy,gradUz
USE MOD_Equation_Vars ,ONLY:mu0,KappasPr,s43,s23
#if PP_VISC==1
USE MOD_Equation_Vars ,ONLY:muSuth
#endif
#if PP_VISC==2
USE MOD_Equation_Vars ,ONLY:ExpoSuth
#endif
#endif /*PARABOLIC*/
#ifdef OPTIMIZED
USE MOD_DG_Vars       ,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem  !< current element treated in volint
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(5,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde !< transformed flux f(iVar,i,j,k)
REAL,DIMENSION(5,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: gtilde !< transformed flux g(iVar,i,j,k)
REAL,DIMENSION(5,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: htilde !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: f(5),g(5),h(5)                          ! Cartesian fluxes (iVar)
REAL                :: srho,e                                  ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,p                              ! auxiliary variables
REAL                :: Ep                                      ! E + p
#if PARABOLIC
REAL                :: f_visc(5),g_visc(5),h_visc(5)           ! viscous cartesian fluxes (iVar)
REAL                :: muS,lambda                              ! viscosity,heat coeff.
REAL                :: gradv1x,gradv2x,gradv3x
REAL                :: gradv1y,gradv2y,gradv3y
REAL                :: gradv1z,gradv2z,gradv3z
REAL                :: divv 
REAL                :: gradTx,gradTy,gradTz
#if (PP_VISC == 1) || (PP_VISC == 2) 
REAL                :: T                                       ! temperature
#endif
#endif /*PARABOLIC*/
INTEGER             :: i 
#ifndef OPTIMIZED
INTEGER             :: j,k
#endif
#ifdef CARTESIANFLUX 
REAL                :: X_xi,Y_eta,Z_zeta
#endif 
!==================================================================================================================================
#ifdef CARTESIANFLUX 
X_xi   = Metrics_fTilde(1,0,0,0,iElem)
Y_eta  = Metrics_gTilde(2,0,0,0,iElem)
Z_zeta = Metrics_hTilde(3,0,0,0,iElem)
#endif 
#ifdef OPTIMIZED
DO i=0,nTotal_vol-1
#else /*OPTIMIZED*/
DO k=0,PP_N;  DO j=0,PP_N; DO i=0,PP_N
#endif /*OPTIMIZED*/
  ASSOCIATE(rho   =>U(1,PP_IJK,iElem), &
            rhov1 =>U(2,PP_IJK,iElem), &
            rhov2 =>U(3,PP_IJK,iElem), &
            rhov3 =>U(4,PP_IJK,iElem), &
            rhoE  =>U(5,PP_IJK,iElem)   )
  ! auxiliary variables
  srho = 1./rho
  v1   = rhov1*srho 
  v2   = rhov2*srho 
  v3   = rhov3*srho 
  p    = kappaM1*(rhoE-0.5*(rhov1*v1+rhov2*v2+rhov3*v3))
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

#if PARABOLIC
  ! Viscous part
  ! ideal gas law
#if PP_VISC == 0
  muS=mu0 ! Constant mu
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
  T=p*srho/R ! Calculate temperature
#endif
#if PP_VISC == 1
  muS=muSuth(T) ! compute viscosity with Sutherlands law
#endif
#if PP_VISC == 2
  muS=mu0*T**ExpoSuth  ! mu0=mu0/T0^n: compute vsicosity using the power-law
#endif
  ! Viscous part
  ! ideal gas law
  e=rhoE*srho                              ! e...specific energy
  ASSOCIATE(  rho_x =>gradUx(1,PP_IJK,iElem),   rho_y =>gradUy(1,PP_IJK,iElem),   rho_z =>gradUz(1,PP_IJK,iElem), &
            rhov1_x =>gradUx(2,PP_IJK,iElem), rhov1_y =>gradUy(2,PP_IJK,iElem), rhov1_z =>gradUz(2,PP_IJK,iElem), &
            rhov2_x =>gradUx(3,PP_IJK,iElem), rhov2_y =>gradUy(3,PP_IJK,iElem), rhov2_z =>gradUz(3,PP_IJK,iElem), &
            rhov3_x =>gradUx(4,PP_IJK,iElem), rhov3_y =>gradUy(4,PP_IJK,iElem), rhov3_z =>gradUz(4,PP_IJK,iElem), &
             rhoE_x =>gradUx(5,PP_IJK,iElem),  rhoE_y =>gradUy(5,PP_IJK,iElem),  rhoE_z =>gradUz(5,PP_IJK,iElem)   )
  ! compute derivatives via product rule (a*b)'=a'*b+a*b'
  gradv1x = srho*(rhov1_x - v1*rho_x)
  gradv1y = srho*(rhov1_y - v1*rho_y)
  gradv1z = srho*(rhov1_z - v1*rho_z)
  gradv2x = srho*(rhov2_x - v2*rho_x)
  gradv2y = srho*(rhov2_y - v2*rho_y)
  gradv2z = srho*(rhov2_z - v2*rho_z)
  gradv3x = srho*(rhov3_x - v3*rho_x)
  gradv3y = srho*(rhov3_y - v3*rho_y)
  gradv3z = srho*(rhov3_z - v3*rho_z)
  divv    = gradv1x+gradv2y+gradv3z
  ! rho*e=p/(kappa-1)+1/2*rho*|v|^2 => e = cv*T+1/2*|v|^2 => cv*T_x= e_x - v*v_x
  gradTx  = srho*(rhoE_x -  e*rho_x) -(v1*gradv1x+v2*gradv2x+v3*gradv3x) !=cv*gradT = grad_e-v*gradv
  gradTy  = srho*(rhoE_y -  e*rho_y) -(v1*gradv1y+v2*gradv2y+v3*gradv3y)
  gradTz  = srho*(rhoE_z -  e*rho_z) -(v1*gradv1z+v2*gradv2z+v3*gradv3z)
  END ASSOCIATE !rho_x/y/z,rhov1_x/y/z ...
  !isotropic heat flux
  lambda=muS*KappasPr
  ! viscous fluxes in x-direction      
  f_visc(2)=-muS*(2*gradv1x-s23*divv)
  f_visc(3)=-muS*(  gradv2x+gradv1y)   
  f_visc(4)=-muS*(  gradv3x+gradv1z)   
  !energy
  f_visc(5)= v1*f_visc(2)+v2*f_visc(3)+v3*f_visc(4) -lambda*gradTx
                                                     
  ! viscous fluxes in y-direction      
  g_visc(2)= f_visc(3)                  !muS*(  gradv1y+gradv2x)  
  g_visc(3)=-muS*(2*gradv2y-s23*divv)     
  g_visc(4)=-muS*(  gradv3y+gradv2z)      
  !energy
  g_visc(5)= v1*g_visc(2)+v2*g_visc(3)+v3*g_visc(4) - lambda*gradTy

  ! viscous fluxes in z-direction      
  h_visc(2)= f_visc(4)                  !muS*(  gradv1z+gradv3x)                 
  h_visc(3)= g_visc(4)                  !muS*(  gradv2z+gradv3y)                
  h_visc(4)=-muS*(2*gradv3z-s23*divv )             
  !energy
  h_visc(5)= v1*h_visc(2)+v2*h_visc(3)+v3*h_visc(4)-lambda*gradTz

  ! add viscous flux
  f(2:5)=f(2:5)+f_visc(2:5)            
  g(2:5)=g(2:5)+g_visc(2:5)
  h(2:5)=h(2:5)+h_visc(2:5)
#endif /*PARABOLIC*/

  END ASSOCIATE !rho,rhov1,rhov2,rhov3,rhoE

  !now transform fluxes to reference ftilde,gtilde,htilde
#ifdef CARTESIANFLUX
  !for cartesian meshes, metric tensor is constant and diagonal:
  ftilde(:,PP_IJK) =  f(:)*X_xi
  gtilde(:,PP_IJK) =  g(:)*Y_eta
  htilde(:,PP_IJK) =  h(:)*Z_zeta
#else /* CURVED FLUX*/
  ! general curved metrics
  ftilde(:,PP_IJK) =   f(:)*Metrics_fTilde(1,PP_IJK,iElem)  &
                     + g(:)*Metrics_fTilde(2,PP_IJK,iElem)  &
                     + h(:)*Metrics_fTilde(3,PP_IJK,iElem)
  gtilde(:,PP_IJK) =   f(:)*Metrics_gTilde(1,PP_IJK,iElem)  &
                     + g(:)*Metrics_gTilde(2,PP_IJK,iElem)  &
                     + h(:)*Metrics_gTilde(3,PP_IJK,iElem)
  htilde(:,PP_IJK) =   f(:)*Metrics_hTilde(1,PP_IJK,iElem)  &
                     + g(:)*Metrics_hTilde(2,PP_IJK,iElem)  &
                     + h(:)*Metrics_hTilde(3,PP_IJK,iElem)
#endif /*CARTESIANFLUX*/
#ifdef OPTIMIZED
END DO ! i
#else /*OPTIMIZED*/
END DO; END DO; END DO ! i,j,k
#endif /*OPTIMIZED*/
END SUBROUTINE EvalFluxTilde3D


!==================================================================================================================================
!> Computes the first cartesian Euler flux (F_x) using the conservative variables for every surface point
!==================================================================================================================================
SUBROUTINE EvalEulerFlux1D(U_Face,F_Face)
! MODULES
USE MOD_Equation_Vars,ONLY:KappaM1
USE MOD_PreProc,ONLY:PP_N
#ifdef OPTIMIZED
USE MOD_DG_Vars,ONLY:nTotal_face
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face(5,0:PP_N,0:PP_N)                    !< state on surface points
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: F_Face(5,0:PP_N,0:PP_N)                    !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: srho                                       ! reciprocal values for density and kappa/Pr
REAL                :: v1,v2,v3,p                                 ! auxiliary variables
INTEGER             :: i 
#ifndef OPTIMIZED
INTEGER             :: j
#endif
!==================================================================================================================================
#ifdef OPTIMIZED
DO i=0,nTotal_face-1
#else /*OPTIMIZED*/
DO j=0,PP_N ; DO i=0,PP_N
#endif /*OPTIMIZED*/
    ! auxiliary variables
    srho = 1. / U_Face(1,PP_ij) ! 1/rho
    v1   = U_Face(2,PP_ij)*srho ! u
    v2   = U_Face(3,PP_ij)*srho ! v
    v3   = U_Face(4,PP_ij)*srho ! w
    p    = kappaM1*(U_Face(5,PP_ij)-0.5*U_Face(1,PP_ij)*(v1*v1+v2*v2+v3*v3))
    ! Euler fluxes x-direction
    F_Face(1,PP_ij)= U_Face(2,PP_ij)          ! rho*u
    F_Face(2,PP_ij)= U_Face(2,PP_ij)*v1+p     ! rho*u²+p
    F_Face(3,PP_ij)= U_Face(2,PP_ij)*v2       ! rho*u*v
    F_Face(4,PP_ij)= U_Face(2,PP_ij)*v3       ! rho*u*w
    F_Face(5,PP_ij)=(U_Face(5,PP_ij) + p)*v1 ! (rho*e+p)*u    
#ifdef OPTIMIZED
END DO !i=0,nTotal
#else /*OPTIMIZED*/
END DO ; END DO !i,j
#endif /*OPTIMIZED*/
END SUBROUTINE EvalEulerFlux1D



#if PARABOLIC
!==================================================================================================================================
!> Compute Navier-Stokes fluxes 1D using tau_12=tau_13=q1=0 for outflow, for surface points!
!==================================================================================================================================
SUBROUTINE EvalDiffFlux1D_Outflow(f,U_Face,gradRho,gradVel) !gradUx_Face,gradUy_Face,gradUz_Face)
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_Equation_Vars,ONLY:kappaM1,s43,s23,KappasPr,R
#if PP_VISC==0
USE MOD_Equation_Vars,ONLY:mu0
#endif
#if PP_VISC==1
USE MOD_Equation_Vars,ONLY:muSuth
#endif
#if PP_VISC==2
USE MOD_Equation_Vars,ONLY:ExpoSuth,mu0
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face( 5,0:PP_N,0:PP_N)   !< state on surface points (iVar,i,j)
REAL,INTENT(IN)     :: gradRho(3,0:PP_N,0:PP_N)   !< rho_x, rho_y, rho_z
REAL,INTENT(IN)     :: gradVel(3,0:PP_N,0:PP_N)   !< u_x, v_y, w_z 

!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: f(5,0:PP_N,0:PP_N)                         !< Cartesian flux in x direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Uin(5)
REAL                :: muS                                        ! viscosity and Temperature,
#if (PP_VISC == 1) || (PP_VISC == 2)
REAL                :: T,pres
#endif
REAL                :: srho                                       ! reciprocal values for density and the value of specific energy
REAL                :: v(3)
REAL                :: gradv_diag(3)                              ! diagonal of velocity and energy gradient matrix
INTEGER             :: p,q
!==================================================================================================================================
f=0.
DO q=0,PP_N
  DO p=0,PP_N
    Uin=U_Face(:,p,q)
    ! auxiliary variables
    srho = 1. / Uin(1) ! 1/rho
    v    = Uin(2:4)*srho
    ! Viscous part
    ! ideal gas law
#if PP_VISC == 0
    ! Constant mu
    muS=mu0
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
    ! Calculate temperature for Sutherland or power-law
    pres = kappaM1*(Uin(5)-0.5*Uin(1)*(SUM(v(:)*v(:)))
    T    = pres*srho/R                      ! T=p/(rho*R)
#endif
#if PP_VISC == 1
    ! compute viscosity with Sutherlands law
    muS=muSuth(T)
#endif
#if PP_VISC == 2
    ! compute vsicosity using the power-law
    muS=mu0*T**ExpoSuth  ! mu0=mu0/T0^n
#endif
    ! compute derivatives via product rule (a*b)'=a'*b+a*b' and multiply with viscosity
    gradv_diag = muS*srho*(gradVel(:,p,q) - v*gradRho(:,p,q))
    ! viscous fluxes in x-direction      
    f(2,p,q)= s43*gradv_diag(1) - s23*(gradv_diag(2) + gradv_diag(3)) ! tau_11
    f(5,p,q)=f(2,p,q)*v(1) !tau_11 * v1
  END DO ! p
END DO ! q 
END SUBROUTINE EvalDiffFlux1D_Outflow


!==================================================================================================================================
!> Compute 3D Navier-Stokes diffusion fluxes using the conservative variables and derivatives for every surface Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D(f,g,h,U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:kappaM1,s43,s23,KappasPr,R
#if PP_VISC==0
USE MOD_Equation_Vars ,ONLY:mu0
#endif
#if PP_VISC==1
USE MOD_Equation_Vars ,ONLY:muSuth
#endif
#if PP_VISC==2
USE MOD_Equation_Vars ,ONLY:ExpoSuth,mu0
#endif
#ifdef OPTIMIZED
USE MOD_DG_Vars       ,ONLY:nTotal_face
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face(5,0:PP_N,0:PP_N)                    !< state on surface points
REAL,INTENT(IN)     :: gradUx_Face(5,0:PP_N,0:PP_N)               !< gradient x 
REAL,INTENT(IN)     :: gradUy_Face(5,0:PP_N,0:PP_N)               !< gradient y
REAL,INTENT(IN)     :: gradUz_Face(5,0:PP_N,0:PP_N)               !< gradient z

!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(5,0:PP_N,0:PP_N),INTENT(OUT) :: f,g,h             !< 3D Cartesian fluxes on each surface point 
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: muS,lambda  
#if (PP_VISC == 1) || (PP_VISC == 2)
REAL                :: T
#endif
REAL                :: srho,e
REAL                :: v1,v2,v3                                  ! auxiliary variables
REAL                :: gradv1x,gradv2x,gradv3x
REAL                :: gradv1y,gradv2y,gradv3y
REAL                :: gradv1z,gradv2z,gradv3z
REAL                :: divv 
REAL                :: gradTx,gradTy,gradTz
INTEGER             :: i 
#ifndef OPTIMIZED
INTEGER             :: j
#endif
!==================================================================================================================================
#ifdef OPTIMIZED
DO i=0,nTotal_face-1
#else /*OPTIMIZED*/
DO j=0,PP_N ; DO i=0,PP_N
#endif /*OPTIMIZED*/
  ASSOCIATE(rho   =>U_Face(1,PP_ij),&
            rhov1 =>U_face(2,PP_ij),&
            rhov2 =>U_face(3,PP_ij),&
            rhov3 =>U_face(4,PP_ij),&
            rhoE  =>U_face(5,PP_ij),&
            rho_x =>gradUx_Face(1,PP_ij),   rho_y =>gradUy_Face(1,PP_ij),   rho_z =>gradUz_Face(1,PP_ij), &
          rhov1_x =>gradUx_Face(2,PP_ij), rhov1_y =>gradUy_Face(2,PP_ij), rhov1_z =>gradUz_Face(2,PP_ij), &
          rhov2_x =>gradUx_Face(3,PP_ij), rhov2_y =>gradUy_Face(3,PP_ij), rhov2_z =>gradUz_Face(3,PP_ij), &
          rhov3_x =>gradUx_Face(4,PP_ij), rhov3_y =>gradUy_Face(4,PP_ij), rhov3_z =>gradUz_Face(4,PP_ij), &
           rhoE_x =>gradUx_Face(5,PP_ij),  rhoE_y =>gradUy_Face(5,PP_ij),  rhoE_z =>gradUz_Face(5,PP_ij)   )
  ! auxiliary variables
  srho = 1. / rho 
  v1   = rhov1*srho
  v2   = rhov2*srho
  v3   = rhov3*srho
  ! Viscous part
  ! ideal gas law
#if PP_VISC == 0
  ! Constant mu
  muS=mu0
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
  ! Calculate temperature for Sutherland or power-law
  T=kappaM1*(rhoE-0.5*(rhov1*v1+rhov2*v2+rhov3*v3))*srho/R                      ! T=p/(rho*R)
#endif
#if PP_VISC == 1
  ! compute viscosity with Sutherlands law
  muS=muSuth(T)
#endif
#if PP_VISC == 2
  ! compute vsicosity using the power-law
  muS=mu0*T**ExpoSuth  ! mu0=mu0/T0^n
#endif
  e=rhoE*srho                              ! e...specific energy
  ! compute derivatives via product rule (a*b)'=a'*b+a*b'
  gradv1x = srho*(rhov1_x - v1*rho_x)
  gradv1y = srho*(rhov1_y - v1*rho_y)
  gradv1z = srho*(rhov1_z - v1*rho_z)
  gradv2x = srho*(rhov2_x - v2*rho_x)
  gradv2y = srho*(rhov2_y - v2*rho_y)
  gradv2z = srho*(rhov2_z - v2*rho_z)
  gradv3x = srho*(rhov3_x - v3*rho_x)
  gradv3y = srho*(rhov3_y - v3*rho_y)
  gradv3z = srho*(rhov3_z - v3*rho_z)
  divv    = gradv1x+gradv2y+gradv3z
  ! rho*e=p/(kappa-1)+1/2*rho*|v|^2 => e = cv*T+1/2*|v|^2 => cv*T_x= e_x - v*v_x
  gradTx  = srho*(rhoE_x -  e*rho_x) -(v1*gradv1x+v2*gradv2x+v3*gradv3x) !=cv*gradT = grad_e-v*gradv
  gradTy  = srho*(rhoE_y -  e*rho_y) -(v1*gradv1y+v2*gradv2y+v3*gradv3y)
  gradTz  = srho*(rhoE_z -  e*rho_z) -(v1*gradv1z+v2*gradv2z+v3*gradv3z)
  END ASSOCIATE !rho,rhov1,rhov2,rhov3,rhoE & rho_x/y/z,rhov1_x/y/z...
  lambda=muS*KappasPr
  ! viscous fluxes in x-direction      
  f(1,PP_ij)=0.
  f(2,PP_ij)=-muS*(2*gradv1x-s23*divv)
  f(3,PP_ij)=-muS*(  gradv2x+gradv1y)   
  f(4,PP_ij)=-muS*(  gradv3x+gradv1z)   
  !energy
  f(5,PP_ij)= v1*f(2,PP_ij)+v2*f(3,PP_ij)+v3*f(4,PP_ij) -lambda*gradTx
                                                     
  ! viscous fluxes in y-direction      
  g(1,PP_ij)= 0. 
  g(2,PP_ij)= f(3,PP_ij)                  !muS*(  gradv1y+gradv2x)  
  g(3,PP_ij)=-muS*(2*gradv2y-s23*divv)     
  g(4,PP_ij)=-muS*(  gradv3y+gradv2z)      
  !energy
  g(5,PP_ij)= v1*g(2,PP_ij)+v2*g(3,PP_ij)+v3*g(4,PP_ij) - lambda*gradTy

  ! viscous fluxes in z-direction      
  h(1,PP_ij)= 0. 
  h(2,PP_ij)= f(4,PP_ij)                  !muS*(  gradv1z+gradv3x)                 
  h(3,PP_ij)= g(4,PP_ij)                  !muS*(  gradv2z+gradv3y)                
  h(4,PP_ij)=-muS*(2*gradv3z-s23*divv )             
  !energy
  h(5,PP_ij)= v1*h(2,PP_ij)+v2*h(3,PP_ij)+v3*h(4,PP_ij)-lambda*gradTz
#ifdef OPTIMIZED
END DO !i=0,nTotal
#else /*OPTIMIZED*/
END DO ; END DO !i,j
#endif /*OPTIMIZED*/
END SUBROUTINE EvalDiffFlux3D


!==================================================================================================================================
!> Compute transformed 3D Navier-Stokes diffusion fluxes for every volume Gauss point of element iElem.
!==================================================================================================================================
SUBROUTINE EvalDiffFluxTilde3D(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars       ,ONLY:U
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Lifting_Vars  ,ONLY:gradUx,gradUy,gradUz
USE MOD_Equation_Vars ,ONLY:mu0,KappasPr,s43,s23
#if PP_VISC==1
USE MOD_Equation_Vars ,ONLY:muSuth,kappaM1,R
#endif
#if PP_VISC==2
USE MOD_Equation_Vars ,ONLY:ExpoSuth,KappaM1,R
#endif
#ifdef OPTIMIZED
USE MOD_DG_Vars,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem  !< current element treated in volint
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(5,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde !< transformed diffusion flux f(iVar,i,j,k)
REAL,DIMENSION(5,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: gtilde !< transformed diffusion flux g(iVar,i,j,k)
REAL,DIMENSION(5,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: htilde !< transformed diffusion flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: f_visc(5),g_visc(5),h_visc(5)           !cartesian fluxes 
REAL                :: srho,e                                  ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3                                ! auxiliary variables
REAL                :: muS,lambda                              ! viscosity,heat coeff.
REAL                :: gradv1x,gradv2x,gradv3x
REAL                :: gradv1y,gradv2y,gradv3y
REAL                :: gradv1z,gradv2z,gradv3z
REAL                :: divv 
REAL                :: gradTx,gradTy,gradTz
#if (PP_VISC == 1) || (PP_VISC == 2) 
REAL                :: T                                       ! temperature
#endif
INTEGER             :: i 
#ifndef OPTIMIZED
INTEGER             :: j,k
#endif
#ifdef CARTESIANFLUX 
REAL                :: X_xi,Y_eta,Z_zeta
#endif 
!==================================================================================================================================
#ifdef CARTESIANFLUX 
X_xi   = Metrics_fTilde(1,0,0,0,iElem)
Y_eta  = Metrics_gTilde(2,0,0,0,iElem)
Z_zeta = Metrics_hTilde(3,0,0,0,iElem)
#endif 
#ifdef OPTIMIZED
DO i=0,nTotal_vol-1
#else /*OPTIMIZED*/
DO k=0,PP_N;  DO j=0,PP_N; DO i=0,PP_N
#endif /*OPTIMIZED*/
  ASSOCIATE(rho   =>U(1,PP_IJK,iElem), &
            rhov1 =>U(2,PP_IJK,iElem), &
            rhov2 =>U(3,PP_IJK,iElem), &
            rhov3 =>U(4,PP_IJK,iElem), &
            rhoE  =>U(5,PP_IJK,iElem), & 
            rho_x =>gradUx(1,PP_IJK,iElem),   rho_y =>gradUy(1,PP_IJK,iElem),   rho_z =>gradUz(1,PP_IJK,iElem), &
          rhov1_x =>gradUx(2,PP_IJK,iElem), rhov1_y =>gradUy(2,PP_IJK,iElem), rhov1_z =>gradUz(2,PP_IJK,iElem), &
          rhov2_x =>gradUx(3,PP_IJK,iElem), rhov2_y =>gradUy(3,PP_IJK,iElem), rhov2_z =>gradUz(3,PP_IJK,iElem), &
          rhov3_x =>gradUx(4,PP_IJK,iElem), rhov3_y =>gradUy(4,PP_IJK,iElem), rhov3_z =>gradUz(4,PP_IJK,iElem), &
           rhoE_x =>gradUx(5,PP_IJK,iElem),  rhoE_y =>gradUy(5,PP_IJK,iElem),  rhoE_z =>gradUz(5,PP_IJK,iElem)   )
  ! auxiliary variables
  srho = 1./rho
  v1   = rhov1*srho 
  v2   = rhov2*srho 
  v3   = rhov3*srho 

  ! Viscous part
  ! ideal gas law
#if PP_VISC == 0
  muS=mu0 ! Constant mu
#endif
#if (PP_VISC == 1) || (PP_VISC == 2)
  T=kappaM1*(rhoE-0.5*(rhov1*v1+rhov2*v2+rhov3*v3))*srho/R ! Calculate temperature
#endif
#if PP_VISC == 1
  muS=muSuth(T) ! compute viscosity with Sutherlands law
#endif
#if PP_VISC == 2
  muS=mu0*T**ExpoSuth  ! mu0=mu0/T0^n: compute vsicosity using the power-law
#endif
  ! Viscous part
  ! ideal gas law
  e=rhoE*srho                              ! e...specific energy
  ! compute derivatives via product rule (a*b)'=a'*b+a*b'
  gradv1x = srho*(rhov1_x - v1*rho_x)
  gradv1y = srho*(rhov1_y - v1*rho_y)
  gradv1z = srho*(rhov1_z - v1*rho_z)
  gradv2x = srho*(rhov2_x - v2*rho_x)
  gradv2y = srho*(rhov2_y - v2*rho_y)
  gradv2z = srho*(rhov2_z - v2*rho_z)
  gradv3x = srho*(rhov3_x - v3*rho_x)
  gradv3y = srho*(rhov3_y - v3*rho_y)
  gradv3z = srho*(rhov3_z - v3*rho_z)
  divv    = gradv1x+gradv2y+gradv3z
  ! rho*e=p/(kappa-1)+1/2*rho*|v|^2 => e = cv*T+1/2*|v|^2 => cv*T_x= e_x - v*v_x
  gradTx  = srho*(rhoE_x -  e*rho_x) -(v1*gradv1x+v2*gradv2x+v3*gradv3x) !=cv*gradT = grad_e-v*gradv
  gradTy  = srho*(rhoE_y -  e*rho_y) -(v1*gradv1y+v2*gradv2y+v3*gradv3y)
  gradTz  = srho*(rhoE_z -  e*rho_z) -(v1*gradv1z+v2*gradv2z+v3*gradv3z)
  END ASSOCIATE !rho,rhov1,rhov2,rhov3,rhoE & rho_x/y/z,rhov1_x/y/z ...

  lambda=muS*KappasPr
  ! viscous fluxes in x-direction      
  f_visc(1)= 0.
  f_visc(2)=-muS*(2*gradv1x-s23*divv)
  f_visc(3)=-muS*(  gradv2x+gradv1y)   
  f_visc(4)=-muS*(  gradv3x+gradv1z)   
  !energy
  f_visc(5)= v1*f_visc(2)+v2*f_visc(3)+v3*f_visc(4) -lambda*gradTx
                                                     
  ! viscous fluxes in y-direction      
  g_visc(1)= 0.
  g_visc(2)= f_visc(3)                  !muS*(  gradv1y+gradv2x)  
  g_visc(3)=-muS*(2*gradv2y-s23*divv)     
  g_visc(4)=-muS*(  gradv3y+gradv2z)      
  !energy
  g_visc(5)= v1*g_visc(2)+v2*g_visc(3)+v3*g_visc(4) - lambda*gradTy

  ! viscous fluxes in z-direction      
  h_visc(1)= 0.
  h_visc(2)= f_visc(4)                  !muS*(  gradv1z+gradv3x)                 
  h_visc(3)= g_visc(4)                  !muS*(  gradv2z+gradv3y)                
  h_visc(4)=-muS*(2*gradv3z-s23*divv )             
  !energy
  h_visc(5)= v1*h_visc(2)+v2*h_visc(3)+v3*h_visc(4)-lambda*gradTz

  !now transform fluxes to reference ftilde,gtilde,htilde
#ifdef CARTESIANFLUX
  !for cartesian meshes, metric tensor is constant and diagonal:
  ftilde(:,PP_IJK) =  f_visc(:) *X_xi
  gtilde(:,PP_IJK) =  g_visc(:) *Y_eta
  htilde(:,PP_IJK) =  h_visc(:) *Z_zeta
#else /* CURVED FLUX*/
  ! general curved metrics
  ftilde(:,PP_IJK) =   f_visc(:)*Metrics_fTilde(1,PP_IJK,iElem)  &
                     + g_visc(:)*Metrics_fTilde(2,PP_IJK,iElem)  &
                     + h_visc(:)*Metrics_fTilde(3,PP_IJK,iElem)
  gtilde(:,PP_IJK) =   f_visc(:)*Metrics_gTilde(1,PP_IJK,iElem)  &
                     + g_visc(:)*Metrics_gTilde(2,PP_IJK,iElem)  &
                     + h_visc(:)*Metrics_gTilde(3,PP_IJK,iElem)
  htilde(:,PP_IJK) =   f_visc(:)*Metrics_hTilde(1,PP_IJK,iElem)  &
                     + g_visc(:)*Metrics_hTilde(2,PP_IJK,iElem)  &
                     + h_visc(:)*Metrics_hTilde(3,PP_IJK,iElem)
#endif /*CARTESIANFLUX*/
#ifdef OPTIMIZED
END DO ! i
#else /*OPTIMIZED*/
END DO; END DO; END DO ! i,j,k
#endif /*OPTIMIZED*/
END SUBROUTINE EvalDiffFluxTilde3D
#endif /*PARABOLIC*/

END MODULE MOD_Flux
