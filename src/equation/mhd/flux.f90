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

INTERFACE EvalFluxTilde3D
  MODULE PROCEDURE EvalFluxTilde3D
END INTERFACE

INTERFACE EvalAdvectionFlux1D
  MODULE PROCEDURE EvalAdvectionFlux1D
END INTERFACE

#if PARABOLIC
! no interface because of dimension change
!INTERFACE EvalDiffFlux3D
!  MODULE PROCEDURE EvalDiffFlux3D
!END INTERFACE

INTERFACE EvalDiffFluxTilde3D
  MODULE PROCEDURE EvalDiffFluxTilde3D
END INTERFACE
#endif /*PARABOLIC*/

PUBLIC::EvalFluxTilde3D
PUBLIC::EvalAdvectionFlux1D
#if PARABOLIC
PUBLIC::EvalDiffFlux3D
PUBLIC::EvalDiffFluxTilde3D
#endif /*PARABOLIC*/

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute transformed MHD fluxes using the conservative variables and gradients for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFluxTilde3D(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,ONLY:U
USE MOD_Equation_Vars ,ONLY:kappaM1,kappaM2,smu_0,s2mu_0
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
#ifdef PP_GLM
USE MOD_Equation_Vars,ONLY: GLM_ch
#endif /*PP_GLM*/
#if PARABOLIC
USE MOD_Lifting_Vars  ,ONLY:gradUx,gradUy,gradUz
USE MOD_Equation_Vars ,ONLY:mu,s23,etasmu_0
#ifdef PP_ANISO_HEAT
USE MOD_Equation_vars ,ONLY:kperp,kpar
#else
USE MOD_Equation_vars ,ONLY:kappasPr
#endif 
#endif /*PARABOLIC*/
#ifdef OPTIMIZED
USE MOD_DG_Vars       ,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem  !< Determines the actual element
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde !< transformed flux f(iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: gtilde !< transformed flux g(iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: htilde !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:PP_nVar) :: f,g,h                             ! Cartesian fluxes (iVar)
REAL                :: srho                                    ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,p                              ! velocity and pressure(including magnetic pressure
REAL                :: bb2,vb                                  ! magnetic field, bb2=|bvec|^2, v dot b
REAL                :: Ep                                      ! E + p
#if PARABOLIC
REAL                :: esrho
REAL                :: divv
REAL                :: gradv1x,gradv2x,gradv3x
REAL                :: gradv1y,gradv2y,gradv3y
REAL                :: gradv1z,gradv2z,gradv3z
REAL                :: gradex,gradey,gradez
REAL                :: lambda 
REAL                :: gradTx,gradTy,gradTz
REAL                :: Qx,Qy,Qz
REAL,DIMENSION(2:PP_nVar) :: f_visc,g_visc,h_visc
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
            Etotal=>U(5,PP_IJK,iElem), &
            b1    =>U(6,PP_IJK,iElem), &
            b2    =>U(7,PP_IJK,iElem), &
            b3    =>U(8,PP_IJK,iElem)  &
#ifdef PP_GLM
           ,Psi   =>U(9,PP_IJK,iElem)  &
#endif /*def PP_GLM*/
            )
  ! auxiliary variables
  srho = 1. / rho ! 1/rho
  v1   = rhov1*srho 
  v2   = rhov2*srho 
  v3   = rhov3*srho 
  bb2  = (b1*b1+b2*b2+b3*b3)
  vb   = (b1*v1+b2*v2+b3*v3)
  !p = ptilde (includes magnetic pressure)
#ifdef PP_GLM
  p    = kappaM1*(Etotal-0.5*(rho*(v1*v1+v2*v2+v3*v3)+smu_0*psi*psi))-KappaM2*s2mu_0*bb2
  Ep   = (Etotal-0.5*smu_0*psi*psi + p)
#else
  p    = kappaM1*(Etotal-0.5*(rho*(v1*v1+v2*v2+v3*v3)))-KappaM2*s2mu_0*bb2
  Ep   = (Etotal + p)
#endif /*def PP_GLM*/
  ! Advection part
  ! Advection fluxes x-direction
  f(1)=rhov1                     ! rho*u
  f(2)=rhov1*v1+p  -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
  f(3)=rhov1*v2    -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
  f(4)=rhov1*v3    -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
  f(5)=Ep*v1       -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
  f(6)=0.
  f(7)=v1*b2-b1*v2
  f(8)=v1*b3-b1*v3
  ! Advection fluxes y-direction
  g(1)=rhov2                     ! rho*v      
  g(2)=f(3)                      ! rho*u*v      -1/mu_0*b2*b1
  g(3)=rhov2*v2+p  -smu_0*b2*b2  ! rho*v²+p     -1/mu_0*b2*b2
  g(4)=rhov2*v3    -smu_0*b2*b3  ! rho*v*w      -1/mu_0*b2*b3
  g(5)=Ep*v2       -smu_0*b2*vb  ! (rho*e+p)*v  -1/mu_0*b2*(v dot B)
  g(6)=-f(7)                     ! v2*b1-b2*v1 
  g(7)=0.
  g(8)=v2*b3-b2*v3
  ! Advection fluxes z-direction
  h(1)=rhov3                     ! rho*v
  h(2)=f(4)                      ! rho*u*w      -1/mu_0*b3*b1
  h(3)=g(4)                      ! rho*v*w      -1/mu_0*b3*b2
  h(4)=rhov3*v3+p  -smu_0*b3*b3  ! rho*v²+p     -1/mu_0*b3*b3
  h(5)=Ep*v3       -smu_0*b3*vb  ! (rho*e+p)*w  -1/mu_0*b3*(v dot B)
  h(6)=-f(8)                     ! v3*b1-b3*v1 
  h(7)=-g(8)                     ! v3*b2-b3*v2
  h(8)=0.

#ifdef PP_GLM
  f(5) = f(5)+smu_0*GLM_ch*b1*U(9,PP_IJK,iElem)
  f(6) = f(6)+GLM_ch   *U(9,PP_IJK,iElem)
  f(9) =      GLM_ch*b1

  g(5) = g(5)+smu_0*GLM_ch*b2*U(9,PP_IJK,iElem)
  g(7) = g(7)+GLM_ch   *U(9,PP_IJK,iElem)
  g(9) =      GLM_ch*b2

  h(5) = h(5)+smu_0*GLM_ch*b3*U(9,PP_IJK,iElem)
  h(8) = h(8)+GLM_ch   *U(9,PP_IJK,iElem)
  h(9) =      GLM_ch*b3
#endif /* PP_GLM */

#if PARABOLIC
  ! Viscous part
  ! ideal gas law
  esrho=Etotal*srho                              ! e...specific energy
  ! compute derivatives via product rule (a*b)'=a'*b+a*b'
  gradv1x = srho*(gradUx(2,PP_IJK,iElem) -   v1*gradUx(1,PP_IJK,iElem))    !gradu(1,1)= u_x, gradu(1,2)=u_y, gradu(1,3)=u_z
  gradv2x = srho*(gradUx(3,PP_IJK,iElem) -   v2*gradUx(1,PP_IJK,iElem))    !gradu(2,1)= v_x, gradu(2,2)=v_y, gradu(2,3)=v_z
  gradv3x = srho*(gradUx(4,PP_IJK,iElem) -   v3*gradUx(1,PP_IJK,iElem))    !gradu(3,1)= w_x, gradu(3,2)=w_y, gradu(3,3)=w_z
  gradex  = srho*(gradUx(5,PP_IJK,iElem) -esrho*gradUx(1,PP_IJK,iElem))    !gradex  = e_x, gradey  =e_y, gradez  =e_z
  gradv1y = srho*(gradUy(2,PP_IJK,iElem) -   v1*gradUy(1,PP_IJK,iElem))
  gradv2y = srho*(gradUy(3,PP_IJK,iElem) -   v2*gradUy(1,PP_IJK,iElem))
  gradv3y = srho*(gradUy(4,PP_IJK,iElem) -   v3*gradUy(1,PP_IJK,iElem))
  gradey  = srho*(gradUy(5,PP_IJK,iElem) -esrho*gradUy(1,PP_IJK,iElem))
  gradv1z = srho*(gradUz(2,PP_IJK,iElem) -   v1*gradUz(1,PP_IJK,iElem))
  gradv2z = srho*(gradUz(3,PP_IJK,iElem) -   v2*gradUz(1,PP_IJK,iElem))
  gradv3z = srho*(gradUz(4,PP_IJK,iElem) -   v3*gradUz(1,PP_IJK,iElem))
  gradez  = srho*(gradUz(5,PP_IJK,iElem) -esrho*gradUz(1,PP_IJK,iElem))
  divv    = gradv1x+gradv2y+gradv3z
ASSOCIATE( gradB1x => gradUx(6,PP_IJK,iElem), & 
           gradB2x => gradUx(7,PP_IJK,iElem), & 
           gradB3x => gradUx(8,PP_IJK,iElem), & 
           gradB1y => gradUy(6,PP_IJK,iElem), & 
           gradB2y => gradUy(7,PP_IJK,iElem), & 
           gradB3y => gradUy(8,PP_IJK,iElem), & 
           gradB1z => gradUz(6,PP_IJK,iElem), & 
           gradB2z => gradUz(7,PP_IJK,iElem), & 
           gradB3z => gradUz(8,PP_IJK,iElem))
#ifdef PP_GLM
ASSOCIATE( gradPsix => gradUx(9,PP_IJK,iElem), & 
           gradPsiy => gradUy(9,PP_IJK,iElem), & 
           gradPsiz => gradUz(9,PP_IJK,iElem))
  gradTx  = gradex           -(v1*gradv1x+v2*gradv2x+v3*gradv3x)    &  !=cv*gradT 
                  -srho*smu_0*((b1*gradB1x+b2*gradB2x+b3*gradB3x)+Psi*gradPsix     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUx(1,PP_IJK,iElem))
  gradTy  = gradey           -(v1*gradv1y+v2*gradv2y+v3*gradv3y)    &
                  -srho*smu_0*((b1*gradB1y+b2*gradB2y+b3*gradB3y)+Psi*gradPsiy     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUy(1,PP_IJK,iElem))
  gradTz  = gradez           -(v1*gradv1z+v2*gradv2z+v3*gradv3z)    &
                  -srho*smu_0*((b1*gradB1z+b2*gradB2z+b3*gradB3z)+Psi*gradPsiz     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUz(1,PP_IJK,iElem))
END ASSOCIATE
#else
  gradTx  = gradex           -(v1*gradv1x+v2*gradv2x+v3*gradv3x)    &  !=cv*gradT 
                  -smu_0*srho*(b1*gradB1x+b2*gradB2x+b3*gradB3x     &
                               -0.5*bb2*sRho*gradUx(1,PP_IJK,iElem))
  gradTy  = gradey           -(v1*gradv1y+v2*gradv2y+v3*gradv3y)    &
                  -smu_0*srho*(b1*gradB1y+b2*gradB2y+b3*gradB3y     &
                               -0.5*bb2*sRho*gradUy(1,PP_IJK,iElem))
  gradTz  = gradez           -(v1*gradv1z+v2*gradv2z+v3*gradv3z)    &
                  -smu_0*srho*(b1*gradB1z+b2*gradB2z+b3*gradB3z     &
                               -0.5*bb2*sRho*gradUz(1,PP_IJK,iElem))
#endif /*P_GLM*/

#ifndef PP_ANISO_HEAT
  !isotropic heat flux
  lambda=mu*KappasPr
  Qx=lambda*gradTx  !q=lambda*gradT= (mu*kappa/Pr)*(cv*gradT)
  Qy=lambda*gradTy
  Qz=lambda*gradTz
#else
  !IF(bb2.GT. 0.)THEN! ATTENTION: |B|^2 /= 0 !!!
  lambda=(kpar-kperp)/bb2  
  Qx=kappaM1*(lambda*(b1*b1*gradTx+b1*b2*gradTy+b1*b3*gradTz)+kperp*gradTx)
  Qy=kappaM1*(lambda*(b2*b1*gradTx+b2*b2*gradTy+b2*b3*gradTz)+kperp*gradTy)
  Qz=kappaM1*(lambda*(b3*b1*gradTx+b3*b2*gradTy+b3*b3*gradTz)+kperp*gradTz)
  !ELSE
  !  Qx=kperp*gradTx
  !  Qy=kperp*gradTy
  !  Qz=kperp*gradTz
  !END IF 
#endif /*PP_ANISO_HEAT*/
  ! viscous fluxes in x-direction      
  f_visc(2)=mu*(2*gradv1x-s23*divv)
  f_visc(3)=mu*(  gradv2x+gradv1y)   
  f_visc(4)=mu*(  gradv3x+gradv1z)   
  f_visc(6)=0.
  f_visc(7)=etasmu_0*(gradB2x-gradB1y)
  f_visc(8)=etasmu_0*(gradB3x-gradB1z)
  !energy
  f_visc(5)= v1*f_visc(2)+v2*f_visc(3)+v3*f_visc(4) +smu_0*(b1*f_visc(6)+b2*f_visc(7)+b3*f_visc(8)) + Qx



  ! viscous fluxes in y-direction      
  g_visc(2)= f_visc(3)                  !mu*(  gradv1y+gradv2x)  
  g_visc(3)=mu*(2*gradv2y-s23*divv)     
  g_visc(4)=mu*(  gradv3y+gradv2z)      
  g_visc(6)=-f_visc(7)                  !etasmu_0*(gradB1y-gradB2x)
  g_visc(7)=0.
  g_visc(8)=etasmu_0*(gradB3y-gradB2z)
  !energy
  g_visc(5)=v1*g_visc(2)+v2*g_visc(3)+v3*g_visc(4) + smu_0*(b1*g_visc(6)+b2*g_visc(7)+b3*g_visc(8)) + Qy 



  ! viscous fluxes in z-direction      
  h_visc(2)= f_visc(4)                       !mu*(  gradv1z+gradv3x)                 
  h_visc(3)= g_visc(4)                       !mu*(  gradv2z+gradv3y)                
  h_visc(4)=mu*(2*gradv3z-s23*divv )             
  h_visc(6)=-f_visc(8)                       !etasmu_0*(gradB1z-gradB3x)
  h_visc(7)=-g_visc(8)                       !etasmu_0*(gradB2z-gradB3y)
  h_visc(8)=0.
  !energy
  h_visc(5)=v1*h_visc(2)+v2*h_visc(3)+v3*h_visc(4)+smu_0*(b1*h_visc(6)+b2*h_visc(7)+b3*h_visc(8)) + Qz

  f(2:8)=f(2:8)-f_visc(2:8)
  g(2:8)=g(2:8)-g_visc(2:8)
  h(2:8)=h(2:8)-h_visc(2:8)

END ASSOCIATE ! gradB1x => gradUx(6 ...

#endif /*PARABOLIC*/
END ASSOCIATE ! rho,rhov1,rhov2,rhov3,Etotal,b1,b2,b3

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
!> Computes the 1D Advection fluxes using the conservative variables for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalAdvectionFlux1D(U_Face,F_Face)
! MODULES
USE MOD_Equation_Vars,ONLY:KappaM1,KappaM2,smu_0,s2mu_0
#ifdef PP_GLM
USE MOD_Equation_vars,ONLY:GLM_ch
#endif 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face(PP_nVar)                    !< input conservative state
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: F_Face(PP_nVar)                    !< F_Face(iVar), Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: srho                               ! reciprocal values for density and kappa/Pr
REAL                :: v1,v2,v3,p                         ! velocity and pressure(including magnetic pressure
REAL                :: vb                                 ! magnetic field, bb2=|bvec|^2
!==================================================================================================================================
! auxiliary variables
srho = 1. / U_Face(1) ! 1/rho
v1   = U_Face(2)*srho ! u
v2   = U_Face(3)*srho ! v
v3   = U_Face(4)*srho ! w
ASSOCIATE( b1 => U_Face(6), &
           b2 => U_Face(7), &
           b3 => U_Face(8), &
#ifdef PP_GLM
           Etotal => U_Face(5)-0.5*smu_0*U_Face(9)**2  &
#else
           Etotal => U_Face(5)  &
#endif
)
vb   = (b1*v1+b2*v2+b3*v3)
p    = kappaM1*(Etotal-0.5*U_Face(1)*(v1*v1+v2*v2+v3*v3))-kappaM2*s2mu_0*(b1*b1+b2*b2+b3*b3)
! Advection fluxes x-direction
F_Face(1)= U_Face(2)          ! rho*u
F_Face(2)= U_Face(2)*v1+p    -smu_0*b1*b1  ! rho*u²+p     -1/mu_0*b1*b1
F_Face(3)= U_Face(2)*v2      -smu_0*b1*b2  ! rho*u*v      -1/mu_0*b1*b2
F_Face(4)= U_Face(2)*v3      -smu_0*b1*b3  ! rho*u*w      -1/mu_0*b1*b3
F_Face(5)=(Etotal + p)*v1 -smu_0*b1*vb  ! (rho*e+p)*u  -1/mu_0*b1*(v dot B)
F_Face(6)=0.
F_Face(7)=v1*b2-b1*v2
F_Face(8)=v1*b3-b1*v3
#ifdef PP_GLM
F_Face(5)=F_Face(5)+smu_0*GLM_ch*b1*U_Face(9)
F_Face(6)=F_Face(6)+GLM_ch   *U_Face(9)
F_Face(9)=          GLM_ch*b1
#endif /* PP_GLM */
END ASSOCIATE ! b1 => UFace(6) ...
END SUBROUTINE EvalAdvectionFlux1D


#if PARABOLIC
!==================================================================================================================================
!> Compute the cartesian diffusion flux of the MHD equations, for all face points (PP_N+1)**2
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D(f,g,h,U_Face,gradUx_Face,gradUy_Face,gradUz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:s23,smu_0
USE MOD_Equation_Vars,ONLY:mu,etasmu_0
USE MOD_DG_Vars,ONLY:nTotal_face
#ifdef PP_ANISO_HEAT
USE MOD_Equation_vars,ONLY:kperp,kpar,kappaM1
#else
USE MOD_Equation_vars,ONLY:kappasPr
#endif 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face(PP_nVar,nTotal_face)         !< state variable 
REAL,INTENT(IN)     :: gradUx_Face(PP_nVar,nTotal_face)    !< x gradient of state 
REAL,INTENT(IN)     :: gradUy_Face(PP_nVar,nTotal_face)    !< y gradient of state 
REAL,INTENT(IN)     :: gradUz_Face(PP_nVar,nTotal_face)    !< z gradient of state 

!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: f       !< Cartesian diffusiion flux in x
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: g       !< Cartesian diffusiion flux in y
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: h       !< Cartesian diffusiion flux in z
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: srho,esrho
REAL                :: v1,v2,v3                                  ! auxiliary variables
REAL                :: bb2
REAL                :: divv
REAL                :: gradv1x,gradv2x,gradv3x
REAL                :: gradv1y,gradv2y,gradv3y
REAL                :: gradv1z,gradv2z,gradv3z
REAL                :: gradex,gradey,gradez
REAL                :: lambda 
REAL                :: gradTx,gradTy,gradTz
REAL                :: Qx,Qy,Qz
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_face
  ! auxiliary variables
  srho = 1. / U_Face(1,i) ! 1/rho
  v1   = U_Face(2,i)*srho
  v2   = U_Face(3,i)*srho
  v3   = U_Face(4,i)*srho
ASSOCIATE( b1 => U_Face(6,i), &
           b2 => U_Face(7,i), &
           b3 => U_Face(8,i)  )
  bb2  = (b1*b1+b2*b2+b3*b3)
  ! Viscous part
  ! ideal gas law
  esrho=U_Face(5,i)*srho                              ! e...specific energy
  ! compute derivatives via product rule (a*b)'=a'*b+a*b'
  gradv1x = srho*(gradUx_Face(2,i) -   v1*gradUx_Face(1,i))    !gradu(1,1)= u_x, gradu(1,2)=u_y, gradu(1,3)=u_z
  gradv2x = srho*(gradUx_Face(3,i) -   v2*gradUx_Face(1,i))    !gradu(2,1)= v_x, gradu(2,2)=v_y, gradu(2,3)=v_z
  gradv3x = srho*(gradUx_Face(4,i) -   v3*gradUx_Face(1,i))    !gradu(3,1)= w_x, gradu(3,2)=w_y, gradu(3,3)=w_z
  gradex  = srho*(gradUx_Face(5,i) -esrho*gradUx_Face(1,i))    !gradex  = e_x, gradey  =e_y, gradez  =e_z
  gradv1y = srho*(gradUy_Face(2,i) -   v1*gradUy_Face(1,i))
  gradv2y = srho*(gradUy_Face(3,i) -   v2*gradUy_Face(1,i))
  gradv3y = srho*(gradUy_Face(4,i) -   v3*gradUy_Face(1,i))
  gradey  = srho*(gradUy_Face(5,i) -esrho*gradUy_Face(1,i))
  gradv1z = srho*(gradUz_Face(2,i) -   v1*gradUz_Face(1,i))
  gradv2z = srho*(gradUz_Face(3,i) -   v2*gradUz_Face(1,i))
  gradv3z = srho*(gradUz_Face(4,i) -   v3*gradUz_Face(1,i))
  gradez  = srho*(gradUz_Face(5,i) -esrho*gradUz_Face(1,i))
  divv    = gradv1x+gradv2y+gradv3z
ASSOCIATE( gradB1x => gradUx_Face(6,i), & 
           gradB2x => gradUx_Face(7,i), & 
           gradB3x => gradUx_Face(8,i), & 
           gradB1y => gradUy_Face(6,i), & 
           gradB2y => gradUy_Face(7,i), & 
           gradB3y => gradUy_Face(8,i), & 
           gradB1z => gradUz_Face(6,i), & 
           gradB2z => gradUz_Face(7,i), & 
           gradB3z => gradUz_Face(8,i) )
#ifdef PP_GLM
ASSOCIATE( Psi =>U_Face(9,i), &
           gradPsix => gradUx_Face(9,i), & 
           gradPsiy => gradUy_Face(9,i), & 
           gradPsiz => gradUz_Face(9,i))
  gradTx  = gradex           -(v1*gradv1x+v2*gradv2x+v3*gradv3x)    &  !=cv*gradT 
                  -srho*smu_0*((b1*gradB1x+b2*gradB2x+b3*gradB3x)+Psi*gradPsix     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUx_Face(1,i))
  gradTy  = gradey           -(v1*gradv1y+v2*gradv2y+v3*gradv3y)    &
                  -srho*smu_0*((b1*gradB1y+b2*gradB2y+b3*gradB3y)+Psi*gradPsiy     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUy_Face(1,i))
  gradTz  = gradez           -(v1*gradv1z+v2*gradv2z+v3*gradv3z)    &
                  -srho*smu_0*((b1*gradB1z+b2*gradB2z+b3*gradB3z)+Psi*gradPsiz     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUz_Face(1,i))
END ASSOCIATE
#else
  gradTx  = gradex           -(v1*gradv1x+v2*gradv2x+v3*gradv3x)    &  !=cv*gradT 
                  -smu_0*srho*(b1*gradB1x+b2*gradB2x+b3*gradB3x     &
                               -0.5*bb2*sRho*gradUx_Face(1,i))
  gradTy  = gradey           -(v1*gradv1y+v2*gradv2y+v3*gradv3y)    &
                  -smu_0*srho*(b1*gradB1y+b2*gradB2y+b3*gradB3y     &
                               -0.5*bb2*sRho*gradUy_Face(1,i))
  gradTz  = gradez           -(v1*gradv1z+v2*gradv2z+v3*gradv3z)    &
                  -smu_0*srho*(b1*gradB1z+b2*gradB2z+b3*gradB3z     &
                               -0.5*bb2*sRho*gradUz_Face(1,i))
#endif /*P_GLM*/
#ifndef PP_ANISO_HEAT
  !isotropic heat flux
  lambda=mu*KappasPr
  Qx=lambda*gradTx  !q=lambda*gradT= (mu*kappa/Pr)*(cv*gradT)
  Qy=lambda*gradTy
  Qz=lambda*gradTz
#else
  !IF(bb2.GT. 0.)THEN! ATTENTION: |B|^2 /= 0 !!!
  lambda=(kpar-kperp)/bb2  
  Qx=kappaM1*(lambda*(b1*b1*gradTx+b1*b2*gradTy+b1*b3*gradTz)+kperp*gradTx)
  Qy=kappaM1*(lambda*(b2*b1*gradTx+b2*b2*gradTy+b2*b3*gradTz)+kperp*gradTy)
  Qz=kappaM1*(lambda*(b3*b1*gradTx+b3*b2*gradTy+b3*b3*gradTz)+kperp*gradTz)
  !ELSE
  !  Qx=kperp*gradTx
  !  Qy=kperp*gradTy
  !  Qz=kperp*gradTz
  !END IF 
#endif /*PP_ANISO_HEAT*/
  ! viscous fluxes in x-direction      
  f(2,i)=-mu*(2*gradv1x-s23*divv)
  f(3,i)=-mu*(  gradv2x+gradv1y)   
  f(4,i)=-mu*(  gradv3x+gradv1z)   
  f(6,i)= 0.
  f(7,i)=-etasmu_0*(gradB2x-gradB1y)
  f(8,i)=-etasmu_0*(gradB3x-gradB1z)
  !energy
  f(5,i)=v1*f(2,i)+v2*f(3,i)+v3*f(4,i) +smu_0*(b1*f(6,i)+b2*f(7,i)+b3*f(8,i)) - Qx



  ! viscous fluxes in y-direction      
  g(2,i)= f(3,i)                  !mu*(  gradv1y+gradv2x)  
  g(3,i)=-mu*(2*gradv2y-s23*divv)     
  g(4,i)=-mu*(  gradv3y+gradv2z)      
  g(6,i)=-f(7,i)                  !etasmu_0*(gradB1y-gradB2x)
  g(7,i)= 0.
  g(8,i)=-etasmu_0*(gradB3y-gradB2z)
  !energy
  g(5,i)=v1*g(2,i)+v2*g(3,i)+v3*g(4,i) + smu_0*(b1*g(6,i)+b2*g(7,i)+b3*g(8,i)) - Qy



  ! viscous fluxes in z-direction      
  h(2,i)= f(4,i)                       !mu*(  gradv1z+gradv3x)                 
  h(3,i)= g(4,i)                       !mu*(  gradv2z+gradv3y)                
  h(4,i)=-mu*(2*gradv3z-s23*divv )             
  h(6,i)=-f(8,i)                       !etasmu_0*(gradB1z-gradB3x)
  h(7,i)=-g(8,i)                       !etasmu_0*(gradB2z-gradB3y)
  h(8,i)= 0.
  !energy
  h(5,i)=v1*h(2,i)+v2*h(3,i)+v3*h(4,i) +smu_0*(b1*h(6,i)+b2*h(7,i)+b3*h(8,i)) - Qz
#ifdef PP_GLM
  f(9,i) = 0.
  g(9,i) = 0.
  h(9,i) = 0.
#endif /*PP_GLM*/
END ASSOCIATE ! gradB1x => gradUx(6,....
END ASSOCIATE ! b1 => UFace(6,i) ...
END DO !i=1,nTotal
END SUBROUTINE EvalDiffFlux3D


!==================================================================================================================================
!> Compute only diffusive part of tranformed MHD fluxes, using conservative variables and derivatives for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFluxTilde3D(iElem,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars,ONLY:U
USE MOD_Mesh_Vars     ,ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
USE MOD_Lifting_Vars  ,ONLY:gradUx,gradUy,gradUz
USE MOD_Equation_Vars ,ONLY:smu_0,mu,s23,etasmu_0
#ifdef PP_ANISO_HEAT
USE MOD_Equation_vars,ONLY:kperp,kpar,kappaM1
#else
USE MOD_Equation_vars,ONLY:kappasPr
#endif 
#ifdef OPTIMIZED
USE MOD_DG_Vars       ,ONLY:nTotal_vol
#endif /*OPTIMIZED*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                                 :: iElem  !< Determines the actual element
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: ftilde !< transformed flux f(iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: gtilde !< transformed flux g(iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N),INTENT(OUT) :: htilde !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(1:PP_nVar) :: f_Visc,g_visc,h_visc              ! Cartesian fluxes (iVar)
REAL                :: srho,e                                  ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,bb2 
REAL                :: divv
REAL                :: gradv1x,gradv2x,gradv3x
REAL                :: gradv1y,gradv2y,gradv3y
REAL                :: gradv1z,gradv2z,gradv3z
REAL                :: gradex,gradey,gradez
REAL                :: lambda 
REAL                :: gradTx,gradTy,gradTz
REAL                :: Qx,Qy,Qz
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
  ! auxiliary variables
  srho = 1. / U(1,PP_IJK,iElem) ! 1/rho
  v1   = U(2,PP_IJK,iElem)*srho 
  v2   = U(3,PP_IJK,iElem)*srho 
  v3   = U(4,PP_IJK,iElem)*srho 
ASSOCIATE( b1 => U(6,PP_IJK,iElem), &
           b2 => U(7,PP_IJK,iElem), &
           b3 => U(8,PP_IJK,iElem)  )
  bb2  = (b1*b1+b2*b2+b3*b3)

  ! Viscous part
  ! ideal gas law
  e=U(5,PP_IJK,iElem)*srho                              ! e...specific energy
  ! compute derivatives via product rule (a*b)'=a'*b+a*b'
  gradv1x = srho*(gradUx(2,PP_IJK,iElem) - v1*gradUx(1,PP_IJK,iElem))    !gradu(1,1)= u_x, gradu(1,2)=u_y, gradu(1,3)=u_z
  gradv2x = srho*(gradUx(3,PP_IJK,iElem) - v2*gradUx(1,PP_IJK,iElem))    !gradu(2,1)= v_x, gradu(2,2)=v_y, gradu(2,3)=v_z
  gradv3x = srho*(gradUx(4,PP_IJK,iElem) - v3*gradUx(1,PP_IJK,iElem))    !gradu(3,1)= w_x, gradu(3,2)=w_y, gradu(3,3)=w_z
  gradex  = srho*(gradUx(5,PP_IJK,iElem) -  e*gradUx(1,PP_IJK,iElem))    !gradex  = e_x, gradey  =e_y, gradez  =e_z
  gradv1y = srho*(gradUy(2,PP_IJK,iElem) - v1*gradUy(1,PP_IJK,iElem))
  gradv2y = srho*(gradUy(3,PP_IJK,iElem) - v2*gradUy(1,PP_IJK,iElem))
  gradv3y = srho*(gradUy(4,PP_IJK,iElem) - v3*gradUy(1,PP_IJK,iElem))
  gradey  = srho*(gradUy(5,PP_IJK,iElem) -  e*gradUy(1,PP_IJK,iElem))
  gradv1z = srho*(gradUz(2,PP_IJK,iElem) - v1*gradUz(1,PP_IJK,iElem))
  gradv2z = srho*(gradUz(3,PP_IJK,iElem) - v2*gradUz(1,PP_IJK,iElem))
  gradv3z = srho*(gradUz(4,PP_IJK,iElem) - v3*gradUz(1,PP_IJK,iElem))
  gradez  = srho*(gradUz(5,PP_IJK,iElem) -  e*gradUz(1,PP_IJK,iElem))
  divv    = gradv1x+gradv2y+gradv3z
ASSOCIATE( gradB1x => gradUx(6,PP_IJK,iElem), & 
           gradB2x => gradUx(7,PP_IJK,iElem), & 
           gradB3x => gradUx(8,PP_IJK,iElem), & 
           gradB1y => gradUy(6,PP_IJK,iElem), & 
           gradB2y => gradUy(7,PP_IJK,iElem), & 
           gradB3y => gradUy(8,PP_IJK,iElem), & 
           gradB1z => gradUz(6,PP_IJK,iElem), & 
           gradB2z => gradUz(7,PP_IJK,iElem), & 
           gradB3z => gradUz(8,PP_IJK,iElem))
#ifdef PP_GLM
ASSOCIATE( Psi =>U(9,PP_IJK,iElem), &
           gradPsix => gradUx(9,PP_IJK,iElem), & 
           gradPsiy => gradUy(9,PP_IJK,iElem), & 
           gradPsiz => gradUz(9,PP_IJK,iElem))
  gradTx  = gradex           -(v1*gradv1x+v2*gradv2x+v3*gradv3x)    &  !=cv*gradT 
                  -srho*smu_0*((b1*gradB1x+b2*gradB2x+b3*gradB3x)+Psi*gradPsix     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUx(1,PP_IJK,iElem))
  gradTy  = gradey           -(v1*gradv1y+v2*gradv2y+v3*gradv3y)    &
                  -srho*smu_0*((b1*gradB1y+b2*gradB2y+b3*gradB3y)+Psi*gradPsiy     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUy(1,PP_IJK,iElem))
  gradTz  = gradez           -(v1*gradv1z+v2*gradv2z+v3*gradv3z)    &
                  -srho*smu_0*((b1*gradB1z+b2*gradB2z+b3*gradB3z)+Psi*gradPsiz     &
                               -0.5*(bb2+Psi*Psi)*sRho*gradUz(1,PP_IJK,iElem))
END ASSOCIATE
#else
  gradTx  = gradex           -(v1*gradv1x+v2*gradv2x+v3*gradv3x)    &  !=cv*gradT 
                  -smu_0*srho*(b1*gradB1x+b2*gradB2x+b3*gradB3x     &
                               -0.5*bb2*sRho*gradUx(1,PP_IJK,iElem))
  gradTy  = gradey           -(v1*gradv1y+v2*gradv2y+v3*gradv3y)    &
                  -smu_0*srho*(b1*gradB1y+b2*gradB2y+b3*gradB3y     &
                               -0.5*bb2*sRho*gradUy(1,PP_IJK,iElem))
  gradTz  = gradez           -(v1*gradv1z+v2*gradv2z+v3*gradv3z)    &
                  -smu_0*srho*(b1*gradB1z+b2*gradB2z+b3*gradB3z     &
                               -0.5*bb2*sRho*gradUz(1,PP_IJK,iElem))
#endif /*P_GLM*/
#ifndef PP_ANISO_HEAT
  !isotropic heat flux
  lambda=mu*KappasPr
  Qx=lambda*gradTx  !q=lambda*gradT= (mu*kappa/Pr)*(cv*gradT)
  Qy=lambda*gradTy
  Qz=lambda*gradTz
#else
  !IF(bb2.GT. 0.)THEN! ATTENTION: |B|^2 /= 0 !!!
  lambda=(kpar-kperp)/bb2  
  Qx=kappaM1*(lambda*(b1*b1*gradTx+b1*b2*gradTy+b1*b3*gradTz)+kperp*gradTx)
  Qy=kappaM1*(lambda*(b2*b1*gradTx+b2*b2*gradTy+b2*b3*gradTz)+kperp*gradTy)
  Qz=kappaM1*(lambda*(b3*b1*gradTx+b3*b2*gradTy+b3*b3*gradTz)+kperp*gradTz)
  !ELSE
  !  Qx=kperp*gradTx
  !  Qy=kperp*gradTy
  !  Qz=kperp*gradTz
  !END IF 
#endif /*PP_ANISO_HEAT*/
  ! viscous fluxes in x-direction      
  f_visc(1)=0.
  f_visc(2)=-mu*(2*gradv1x-s23*divv)
  f_visc(3)=-mu*(  gradv2x+gradv1y)   
  f_visc(4)=-mu*(  gradv3x+gradv1z)   
  f_visc(6)=0.
  f_visc(7)=-etasmu_0*(gradB2x-gradB1y)
  f_visc(8)=-etasmu_0*(gradB3x-gradB1z)
  !energy
  f_visc(5)= v1*f_visc(2)+v2*f_visc(3)+v3*f_visc(4) +smu_0*(b1*f_visc(6)+b2*f_visc(7)+b3*f_visc(8)) - Qx



  ! viscous fluxes in y-direction      
  g_visc(1)=0.
  g_visc(2)= f_visc(3)                  !-mu*(  gradv1y+gradv2x)  
  g_visc(3)=-mu*(2*gradv2y-s23*divv)     
  g_visc(4)=-mu*(  gradv3y+gradv2z)      
  g_visc(6)=-f_visc(7)                  !etasmu_0*(gradB1y-gradB2x)
  g_visc(7)=0.
  g_visc(8)=-etasmu_0*(gradB3y-gradB2z)
  !energy
  g_visc(5)=v1*g_visc(2)+v2*g_visc(3)+v3*g_visc(4) + smu_0*(b1*g_visc(6)+b2*g_visc(7)+b3*g_visc(8)) - Qy 



  ! viscous fluxes in z-direction      
  h_visc(1)=0.
  h_visc(2)= f_visc(4)                       !-mu*(  gradv1z+gradv3x)                 
  h_visc(3)= g_visc(4)                       !-mu*(  gradv2z+gradv3y)                
  h_visc(4)=-mu*(2*gradv3z-s23*divv )             
  h_visc(6)=-f_visc(8)                       !etasmu_0*(gradB1z-gradB3x)
  h_visc(7)=-g_visc(8)                       !etasmu_0*(gradB2z-gradB3y)
  h_visc(8)=0.
  !energy
  h_visc(5)=v1*h_visc(2)+v2*h_visc(3)+v3*h_visc(4)+smu_0*(b1*h_visc(6)+b2*h_visc(7)+b3*h_visc(8)) - Qz
#ifdef PP_GLM
  f_visc(9) = 0. 
  g_visc(9) = 0. 
  h_visc(9) = 0. 
#endif /*PP_GLM*/


END ASSOCIATE ! gradB1x => gradUx(6 ...

END ASSOCIATE ! b1 = U(6,....

  !now transform fluxes to reference ftilde,gtilde,htilde
#ifdef CARTESIANFLUX
  !for cartesian meshes, metric tensor is constant and diagonal:
  ftilde(:,PP_IJK) =  f_visc(:)*X_xi
  gtilde(:,PP_IJK) =  g_visc(:)*Y_eta
  htilde(:,PP_IJK) =  h_visc(:)*Z_zeta
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
