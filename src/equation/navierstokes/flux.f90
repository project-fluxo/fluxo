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
!> Contains routines to evaluate the Euler and diffusion part of the Navier-Stokes flux
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------

!INTERFACE EvalFluxTilde3D
!  MODULE PROCEDURE EvalFluxTilde3D
!END INTERFACE

!INTERFACE EvalEulerFlux1D
!  MODULE PROCEDURE EvalEulerFlux1D
!END INTERFACE


#if PARABOLIC
!INTERFACE EvalDiffFlux1D_Outflow
!  MODULE PROCEDURE EvalDiffFlux1D_Outflow
!END INTERFACE

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
PUBLIC::EvalEulerFlux1D
#if PARABOLIC
PUBLIC::EvalDiffFlux1D_Outflow
PUBLIC::EvalDiffFlux3D
PUBLIC::EvalDiffFluxTilde3D
PUBLIC::EvalLiftingVolumeFlux
PUBLIC::EvalLiftingSurfFlux
#endif /*PARABOLIC*/

!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Compute Euler fluxes using the conservative variables and derivatives for every volume Gauss point.
!> directly apply metrics and output the tranformed flux 
!==================================================================================================================================
SUBROUTINE EvalAdvFluxTilde3D(U_in,M_f,M_g,M_h,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:kappaM1
USE MOD_DG_Vars       ,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN )   :: U_in(5,1:nTotal_vol) !< solution state (conservative vars)
REAL,INTENT(IN )   :: M_f( 3,1:nTotal_vol) !< metrics for ftilde                 
REAL,INTENT(IN )   :: M_g( 3,1:nTotal_vol) !< metrics for gtilde                 
REAL,INTENT(IN )   :: M_h( 3,1:nTotal_vol) !< metrics for htilde                 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(5,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(5,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(5,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: f(5),g(5),h(5)                    !cartesi an fluxes 
REAL                :: srho                              ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,v_2,p                    ! auxiliary variables
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  ! auxiliary variables
  srho = 1./U_in(1,i)
  v1   = U_in(2,i)*srho 
  v2   = U_in(3,i)*srho 
  v3   = U_in(4,i)*srho 
  v_2  = v1*v1+v2*v2+v3*v3 
  p    = kappaM1*(U_in(5,i)-0.5*U_in(1,i)*v_2)
!  ! Euler part
  ! Euler fluxes x-direction
  f(1)=U_in(2,i)         
  f(2)=U_in(2,i)*v1+p    
  f(3)=U_in(2,i)*v2      
  f(4)=U_in(2,i)*v3      
  f(5)=(U_in(5,i)+p)*v1         
  ! Euler fluxes y-direction
  g(1)=U_in(3,i)
  g(2)=f(3)                      ! rho*u*v
  g(3)=U_in(2,i)*v2+p  
  g(4)=U_in(2,i)*v3
  g(5)=(U_in(5,i)+p)*v2 
  ! Euler fluxes z-direction
  h(1)=U_in(4,i)
  h(2)=f(4)               ! rho*v1*v3
  h(3)=g(4)               ! rho*v2*v3  
  h(4)=U_in(4,i)*v3+p    
  h(5)=(U_in(5,i)+p)*v3  
  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(:,i) =   f(:)*M_f(1,i) + g(:)*M_f(2,i) + h(:)*M_f(3,i)
  gtilde(:,i) =   f(:)*M_g(1,i) + g(:)*M_g(2,i) + h(:)*M_g(3,i)
  htilde(:,i) =   f(:)*M_h(1,i) + g(:)*M_h(2,i) + h(:)*M_h(3,i)

END DO ! i
END SUBROUTINE EvalAdvFluxTilde3D

!==================================================================================================================================
!> Compute transformed 3D Navier-Stokes fluxes(Euler+diffusion) for every volume Gauss point of element iElem.
!> In comparison to EvalFlux3D, metrics are directly applied 
!==================================================================================================================================
SUBROUTINE EvalFluxTilde3D(iElem,U_in,M_f,M_g,M_h, &
#if PARABOLIC
                           gradPx_in,gradPy_in,gradPz_in,&
#endif /*PARABOLIC*/
                           ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
#if SHOCK_ARTVISC
USE MOD_ShockCapturing_Vars,ONLY:nu
#endif /*SHOCK_ARTVISC*/
USE MOD_Equation_Vars      ,ONLY:kappaM1
#if PARABOLIC
USE MOD_Equation_Vars      ,ONLY:mu0,sKappaM1,KappasPr,s23
#if PP_VISC==1
USE MOD_Equation_Vars      ,ONLY:muSuth,R
#endif
#if PP_VISC==2
USE MOD_Equation_Vars      ,ONLY:ExpoPow,R
#endif
#endif /*PARABOLIC*/
USE MOD_DG_Vars            ,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN ):: iElem                !< element number
REAL,INTENT(IN )   :: U_in(5,1:nTotal_vol) !< state in conservative variables
REAL,INTENT(IN )   :: M_f( 3,1:nTotal_vol) !< metrics for ftilde                 
REAL,INTENT(IN )   :: M_g( 3,1:nTotal_vol) !< metrics for gtilde                 
REAL,INTENT(IN )   :: M_h( 3,1:nTotal_vol) !< metrics for htilde                 
#if PARABOLIC
REAL,INTENT(IN )   :: gradPx_in(5,1:nTotal_vol) !< gradient x in primitive variables 
REAL,INTENT(IN )   :: gradPy_in(5,1:nTotal_vol) !< gradient y in primitive variables 
REAL,INTENT(IN )   :: gradPz_in(5,1:nTotal_vol) !< gradient z in primitive variables 
#endif /*PARABOLIC*/
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(5,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(5,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(5,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: f(5),g(5),h(5)                 ! Cartesian fluxes (iVar)
REAL                :: srho                           ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,p                     ! auxiliary variables
#if PARABOLIC
REAL                :: f_visc(5),g_visc(5),h_visc(5)           ! viscous cartesian fluxes (iVar)
REAL                :: muS,lambda                              ! viscosity,heat coeff.
REAL                :: mu_eff                                  ! effective viscosity can be a mix of physical and artificial
REAL                :: divv 
REAL                :: cv_gradTx,cv_gradTy,cv_gradTz
#if (PP_VISC == 1) || (PP_VISC == 2) 
REAL                :: T                              ! temperature
#endif
#endif /*PARABOLIC*/
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  ASSOCIATE(rho   =>U_in(1,i), &
            rhov1 =>U_in(2,i), &
            rhov2 =>U_in(3,i), &
            rhov3 =>U_in(4,i), &
            rhoE  =>U_in(5,i)   )
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
  muS=mu0*T**ExpoPow  ! mu0=mu0/T0^n: compute vsicosity using the power-law
#endif
#if SHOCK_ARTVISC
  mu_eff = muS+nu(iElem)
#else
  mu_eff = muS
#endif /*SHOCK_ARTVISC*/
  ! Viscous part
  ASSOCIATE( gradv1x => gradPx_in(2,i), & 
             gradv2x => gradPx_in(3,i), & 
             gradv3x => gradPx_in(4,i), & 
             gradv1y => gradPy_in(2,i), & 
             gradv2y => gradPy_in(3,i), & 
             gradv3y => gradPy_in(4,i), & 
             gradv1z => gradPz_in(2,i), & 
             gradv2z => gradPz_in(3,i), & 
             gradv3z => gradPz_in(4,i)  )
  
  divv    = gradv1x+gradv2y+gradv3z
  cv_gradTx  = sKappaM1*sRho*(gradPx_in(5,i)-srho*p*gradPx_in(1,i))  ! cv*T_x = 1/(kappa-1) *1/rho *(p_x - p/rho*rho_x)
  cv_gradTy  = sKappaM1*sRho*(gradPy_in(5,i)-srho*p*gradPy_in(1,i)) 
  cv_gradTz  = sKappaM1*sRho*(gradPz_in(5,i)-srho*p*gradPz_in(1,i)) 
  
  !isotropic heat flux
  lambda=mu_eff*KappasPr
  ! viscous fluxes in x-direction      
  f_visc(2)=-mu_eff*(2*gradv1x-s23*divv)
  f_visc(3)=-mu_eff*(  gradv2x+gradv1y)
  f_visc(4)=-mu_eff*(  gradv3x+gradv1z)
  !energy
  f_visc(5)= v1*f_visc(2)+v2*f_visc(3)+v3*f_visc(4) -lambda*cv_gradTx
                                                     
  ! viscous fluxes in y-direction      
  g_visc(2)= f_visc(3)                  !muS*(  gradv1y+gradv2x)  
  g_visc(3)=-mu_eff*(2*gradv2y-s23*divv)
  g_visc(4)=-mu_eff*(  gradv3y+gradv2z)
  !energy
  g_visc(5)= v1*g_visc(2)+v2*g_visc(3)+v3*g_visc(4) - lambda*cv_gradTy

  ! viscous fluxes in z-direction      
  h_visc(2)= f_visc(4)                  !muS*(  gradv1z+gradv3x)                 
  h_visc(3)= g_visc(4)                  !muS*(  gradv2z+gradv3y)                
  h_visc(4)=-mu_eff*(2*gradv3z-s23*divv )
  !energy
  h_visc(5)= v1*h_visc(2)+v2*h_visc(3)+v3*h_visc(4)-lambda*cv_gradTz

  ! add viscous flux
  f(2:5)=f(2:5)+f_visc(2:5)            
  g(2:5)=g(2:5)+g_visc(2:5)
  h(2:5)=h(2:5)+h_visc(2:5)
  END ASSOCIATE !v_x/y/z...
#endif /*PARABOLIC*/

  END ASSOCIATE !rho,rhov1,rhov2,rhov3,rhoE

  !now transform fluxes to reference ftilde,gtilde,htilde
  ftilde(:,i) =   f(:)*M_f(1,i) + g(:)*M_f(2,i) + h(:)*M_f(3,i)
  gtilde(:,i) =   f(:)*M_g(1,i) + g(:)*M_g(2,i) + h(:)*M_g(3,i)
  htilde(:,i) =   f(:)*M_h(1,i) + g(:)*M_h(2,i) + h(:)*M_h(3,i)
END DO ! i
END SUBROUTINE EvalFluxTilde3D


!==================================================================================================================================
!> Computes the first cartesian Euler flux (F_x) using the conservative variables for every surface point
!==================================================================================================================================
SUBROUTINE EvalEulerFlux1D(U_Face,F_Face)
! MODULES
USE MOD_Equation_Vars,ONLY:KappaM1
USE MOD_DG_Vars,ONLY:nTotal_face
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face(5,nTotal_Face)   !< state (conservative vars) on surface points
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: F_Face(5,nTotal_Face)   !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                :: srho                    ! reciprocal values for density and kappa/Pr
REAL                :: v1,v2,v3,p              ! auxiliary variables
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_face
    ! auxiliary variables
    srho = 1. / U_Face(1,i) ! 1/rho
    v1   = U_Face(2,i)*srho ! u
    v2   = U_Face(3,i)*srho ! v
    v3   = U_Face(4,i)*srho ! w
    p    = kappaM1*(U_Face(5,i)-0.5*(U_Face(2,i)*v1+U_Face(3,i)*v2+U_Face(4,i)*v3))
    ! Euler fluxes x-direction
    F_Face(1,i)= U_Face(2,i)          ! rho*u
    F_Face(2,i)= U_Face(2,i)*v1+p     ! rho*u²+p
    F_Face(3,i)= U_Face(2,i)*v2       ! rho*u*v
    F_Face(4,i)= U_Face(2,i)*v3       ! rho*u*w
    F_Face(5,i)=(U_Face(5,i) + p)*v1 ! (rho*e+p)*u    
END DO !i=1,nTotal
END SUBROUTINE EvalEulerFlux1D



#if PARABOLIC
!==================================================================================================================================
!> Compute Navier-Stokes fluxes 1D using tau_12=tau_13=q1=0 for outflow, for surface points!
!==================================================================================================================================
SUBROUTINE EvalDiffFlux1D_Outflow(f,U_Face,gradVel) !gradPx_Face,gradPy_Face,gradPz_Face)
! MODULES
USE MOD_PreProc ! PP_N
USE MOD_Equation_Vars,ONLY:s23
#if PP_VISC==0
USE MOD_Equation_Vars,ONLY:mu0
#endif
#if PP_VISC==1
USE MOD_Equation_Vars,ONLY:kappaM1,R,muSuth
#endif
#if PP_VISC==2
USE MOD_Equation_Vars,ONLY:kappaM1,R,ExpoPow,mu0
#endif
USE MOD_DG_Vars,ONLY:nTotal_face
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face( 5,nTotal_Face)   !< state (conservative vars) on surface points (iVar,i,j)
REAL,INTENT(IN)     :: gradVel(3,nTotal_Face)   !< u_x, v_y, w_z 

!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: f(5,nTotal_Face)         !< Cartesian flux in x direction
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: Uin(5)
REAL                :: muS                     ! viscosity and Temperature,
#if (PP_VISC == 1) || (PP_VISC == 2)
REAL                :: T
#endif
REAL                :: srho                    ! reciprocal values for density and the value of specific energy
REAL                :: v(3)
REAL                :: gradv_diag(3)           ! diagonal of velocity and energy gradient matrix
INTEGER             :: i
!==================================================================================================================================
f=0.
DO i=1,nTotal_Face
  Uin=U_Face(:,i)
  ! auxiliary variables
  srho = 1. / Uin(1) ! 1/rho
  v    = Uin(2:4)*srho
  ! Viscous part
  ! ideal gas law
#if PP_VISC == 0
  ! Constant mu
  muS=mu0
#elif PP_VISC == 1
  ! compute viscosity with Sutherlands law
  T  = kappaM1*(Uin(5)*srho-0.5*(SUM(v(:)*v(:))))/R             ! T=p/(rho*R)
  muS=muSuth(T)
#elif PP_VISC == 2
  ! compute vsicosity using the power-law
  T  = kappaM1*(Uin(5)*srho-0.5*(SUM(v(:)*v(:))))/R             ! T=p/(rho*R)
  muS=mu0*T**ExpoPow  ! mu0=mu0/T0^n
#endif
  ! compute derivatives via product rule (a*b)'=a'*b+a*b' and multiply with viscosity
  gradv_diag = muS*gradVel(:,i)
  ! viscous fluxes in x-direction      
  f(2,i)= 2.*gradv_diag(1) - s23*(gradv_diag(1)+gradv_diag(2) + gradv_diag(3)) ! tau_11
  f(5,i)=f(2,i)*v(1) !tau_11 * v1
END DO ! i 
END SUBROUTINE EvalDiffFlux1D_Outflow


!==================================================================================================================================
!> Compute 3D Navier-Stokes diffusion fluxes using the conservative variables and derivatives for every surface Gauss point.
!==================================================================================================================================
SUBROUTINE EvalDiffFlux3D(f,g,h,U_Face,gradPx_Face,gradPy_Face,gradPz_Face)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars ,ONLY:s23,KappasPr,kappaM1,skappaM1
#if PP_VISC==0
USE MOD_Equation_Vars ,ONLY:mu0
#endif
#if PP_VISC==1
USE MOD_Equation_Vars ,ONLY:R,muSuth
#endif
#if PP_VISC==2
USE MOD_Equation_Vars ,ONLY:R,ExpoPow,mu0
#endif
USE MOD_DG_Vars       ,ONLY:nTotal_face
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: U_Face(     5,nTotal_face)    !< state in conservative variables on surface points 
REAL,INTENT(IN)     :: gradPx_Face(5,nTotal_face)    !< x gradient of state 
REAL,INTENT(IN)     :: gradPy_Face(5,nTotal_face)    !< y gradient of state 
REAL,INTENT(IN)     :: gradPz_Face(5,nTotal_face)    !< z gradient of state 

!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: f       !< Cartesian diffusion flux in x
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: g       !< Cartesian diffusion flux in y
REAL,DIMENSION(PP_nVar,nTotal_face),INTENT(OUT) :: h       !< Cartesian diffusion flux in z
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: muS,lambda  
#if (PP_VISC == 1) || (PP_VISC == 2)
REAL                :: T
#endif
REAL                :: srho,p
REAL                :: v1,v2,v3                                ! auxiliary variables
REAL                :: divv 
REAL                :: cv_gradTx,cv_gradTy,cv_gradTz
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_face
  ASSOCIATE(rho   =>U_Face(1,i),&
            rhov1 =>U_face(2,i),&
            rhov2 =>U_face(3,i),&
            rhov3 =>U_face(4,i),&
            rhoE  =>U_face(5,i),&
            gradv1x => gradPx_Face(2,i), & 
            gradv2x => gradPx_Face(3,i), & 
            gradv3x => gradPx_Face(4,i), & 
            gradv1y => gradPy_Face(2,i), & 
            gradv2y => gradPy_Face(3,i), & 
            gradv3y => gradPy_Face(4,i), & 
            gradv1z => gradPz_Face(2,i), & 
            gradv2z => gradPz_Face(3,i), & 
            gradv3z => gradPz_Face(4,i)  )
  ! auxiliary variables
  srho = 1. / rho 
  v1 = rhov1*srho
  v2 = rhov2*srho
  v3 = rhov3*srho
  p  = kappaM1*(rhoE-0.5*(rhov1*v1+rhov2*v2+rhov3*v3))
  ! Viscous part
  ! ideal gas law
#if PP_VISC == 0
  ! Constant mu
  muS=mu0
  ! Calculate temperature for Sutherland or power-law
#elif PP_VISC == 1
  ! compute viscosity with Sutherlands law
  T=p*srho/R                      ! T=p/(rho*R)
  muS=muSuth(T)
#elif PP_VISC == 2
  ! compute vsicosity using the power-law
  T=p*srho/R                      ! T=p/(rho*R)
  muS=mu0*T**ExpoPow  ! mu0=mu0/T0^n
#endif
  ! compute derivatives via product rule (a*b)'=a'*b+a*b'
  divv    = gradv1x+gradv2y+gradv3z
  cv_gradTx  = sKappaM1*sRho*(gradPx_Face(5,i)-srho*p*gradPx_Face(1,i))  ! cv*T_x = 1/(kappa-1) *1/rho *(p_x - p/rho*rho_x)
  cv_gradTy  = sKappaM1*sRho*(gradPy_Face(5,i)-srho*p*gradPy_Face(1,i)) 
  cv_gradTz  = sKappaM1*sRho*(gradPz_Face(5,i)-srho*p*gradPz_Face(1,i)) 
  lambda=muS*KappasPr
  ! viscous fluxes in x-direction      
  f(1,i)=0.
  f(2,i)=-muS*(2*gradv1x-s23*divv)
  f(3,i)=-muS*(  gradv2x+gradv1y)   
  f(4,i)=-muS*(  gradv3x+gradv1z)   
  !energy
  f(5,i)= v1*f(2,i)+v2*f(3,i)+v3*f(4,i) -lambda*cv_gradTx
                                                     
  ! viscous fluxes in y-direction      
  g(1,i)= 0. 
  g(2,i)= f(3,i)                  !muS*(  gradv1y+gradv2x)  
  g(3,i)=-muS*(2*gradv2y-s23*divv)     
  g(4,i)=-muS*(  gradv3y+gradv2z)      
  !energy
  g(5,i)= v1*g(2,i)+v2*g(3,i)+v3*g(4,i) - lambda*cv_gradTy

  ! viscous fluxes in z-direction      
  h(1,i)= 0. 
  h(2,i)= f(4,i)                  !muS*(  gradv1z+gradv3x)                 
  h(3,i)= g(4,i)                  !muS*(  gradv2z+gradv3y)                
  h(4,i)=-muS*(2*gradv3z-s23*divv )             
  !energy
  h(5,i)= v1*h(2,i)+v2*h(3,i)+v3*h(4,i)-lambda*cv_gradTz
  END ASSOCIATE !rho,rhov1,rhov2,rhov3,rhoE & v_x/y/z...
END DO !i=1,nTotal
END SUBROUTINE EvalDiffFlux3D


!==================================================================================================================================
!> Compute transformed 3D Navier-Stokes diffusion fluxes for every volume Gauss point of element iElem.
!==================================================================================================================================
SUBROUTINE EvalDiffFluxTilde3D(iElem,U_in,M_f,M_g,M_h,gradPx_in,gradPy_in,gradPz_in,ftilde,gtilde,htilde)
! MODULES
USE MOD_PreProc
#if SHOCK_ARTVISC
USE MOD_ShockCapturing_Vars,ONLY:nu
#endif /*SHOCK_ARTVISC*/
USE MOD_Equation_Vars ,ONLY:mu0,kappaM1,sKappaM1,KappasPr,s23
#if PP_VISC==1
USE MOD_Equation_Vars ,ONLY:muSuth,R
#endif
#if PP_VISC==2
USE MOD_Equation_Vars ,ONLY:ExpoPow,R
#endif
USE MOD_DG_Vars,ONLY:nTotal_vol
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN ):: iElem                     !< element number
REAL,INTENT(IN )   :: U_in(     5,1:nTotal_vol) !< state in conservative variables
REAL,INTENT(IN )   :: M_f(      3,1:nTotal_vol) !< metrics for ftilde  
REAL,INTENT(IN )   :: M_g(      3,1:nTotal_vol) !< metrics for gtilde  
REAL,INTENT(IN )   :: M_h(      3,1:nTotal_vol) !< metrics for htilde  
REAL,INTENT(IN )   :: gradPx_in(5,1:nTotal_vol) !< gradient x in primitive variables 
REAL,INTENT(IN )   :: gradPy_in(5,1:nTotal_vol) !< gradient y in primitive variables 
REAL,INTENT(IN )   :: gradPz_in(5,1:nTotal_vol) !< gradient z in primitive variables 
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: ftilde(5,1:nTotal_vol) !< transformed flux f(iVar,i,j,k)
REAL,INTENT(OUT)   :: gtilde(5,1:nTotal_vol) !< transformed flux g(iVar,i,j,k)
REAL,INTENT(OUT)   :: htilde(5,1:nTotal_vol) !< transformed flux h(iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: f_visc(5),g_visc(5),h_visc(5)           !cartesian fluxes 
REAL                :: srho,e                                  ! reciprocal values for density and the value of specific energy
REAL                :: v1,v2,v3,p                              ! auxiliary variables
REAL                :: muS,lambda                              ! viscosity,heat coeff.
REAL                :: mu_eff                                  ! effective viscosity can be a mix of physical and artificial
REAL                :: divv 
REAL                :: cv_gradTx,cv_gradTy,cv_gradTz
#if (PP_VISC == 1) || (PP_VISC == 2) 
REAL                :: T                                       ! temperature
#endif
INTEGER             :: i 
!==================================================================================================================================
DO i=1,nTotal_vol
  ASSOCIATE(rho   =>U_in(1,i), &
            rhov1 =>U_in(2,i), &
            rhov2 =>U_in(3,i), &
            rhov3 =>U_in(4,i), &
            rhoE  =>U_in(5,i), & 
            gradv1x => gradPx_in(2,i), & 
            gradv2x => gradPx_in(3,i), & 
            gradv3x => gradPx_in(4,i), & 
            gradv1y => gradPy_in(2,i), & 
            gradv2y => gradPy_in(3,i), & 
            gradv3y => gradPy_in(4,i), & 
            gradv1z => gradPz_in(2,i), & 
            gradv2z => gradPz_in(3,i), & 
            gradv3z => gradPz_in(4,i))
  ! auxiliary variables
  srho = 1./rho
  v1   = rhov1*srho 
  v2   = rhov2*srho 
  v3   = rhov3*srho 

  p    = kappaM1*(rhoE - 0.5*(rhov1*v1+rhov2*v2+rhov3*v3)) 
  ! Viscous part
  ! ideal gas law
#if PP_VISC == 0
  muS=mu0 ! Constant mu
#elif PP_VISC == 1
  T=p*srho/R ! Calculate temperature
  muS=muSuth(T) ! compute viscosity with Sutherlands law
#elif PP_VISC == 2
  T=p*srho/R ! Calculate temperature
  muS=mu0*T**ExpoPow  ! mu0=mu0/T0^n: compute vsicosity using the power-law
#endif
#if SHOCK_ARTVISC
  mu_eff = muS+nu(iElem)
#else
  mu_eff = muS
#endif /*SHOCK_ARTVISC*/
  ! Viscous part
  divv       = gradv1x+gradv2y+gradv3z
  cv_gradTx  = sKappaM1*sRho*(gradPx_in(5,i)-srho*p*gradPx_in(1,i))  ! cv*T_x = 1/(kappa-1) *1/rho *(p_x - p/rho*rho_x)
  cv_gradTy  = sKappaM1*sRho*(gradPy_in(5,i)-srho*p*gradPy_in(1,i)) 
  cv_gradTz  = sKappaM1*sRho*(gradPz_in(5,i)-srho*p*gradPz_in(1,i)) 

  lambda=mu_eff*KappasPr
  ! viscous fluxes in x-direction      
  f_visc(1)= 0.
  f_visc(2)=-mu_eff*(2*gradv1x-s23*divv)
  f_visc(3)=-mu_eff*(  gradv2x+gradv1y)
  f_visc(4)=-mu_eff*(  gradv3x+gradv1z)
  !energy
  f_visc(5)= v1*f_visc(2)+v2*f_visc(3)+v3*f_visc(4) -lambda*cv_gradTx
                                                     
  ! viscous fluxes in y-direction      
  g_visc(1)= 0.
  g_visc(2)= f_visc(3)                  !muS*(  gradv1y+gradv2x)  
  g_visc(3)=-mu_eff*(2*gradv2y-s23*divv)
  g_visc(4)=-mu_eff*(  gradv3y+gradv2z)
  !energy
  g_visc(5)= v1*g_visc(2)+v2*g_visc(3)+v3*g_visc(4) - lambda*cv_gradTy

  ! viscous fluxes in z-direction      
  h_visc(1)= 0.
  h_visc(2)= f_visc(4)                  !muS*(  gradv1z+gradv3x)                 
  h_visc(3)= g_visc(4)                  !muS*(  gradv2z+gradv3y)                
  h_visc(4)=-mu_eff*(2*gradv3z-s23*divv )
  !energy
  h_visc(5)= v1*h_visc(2)+v2*h_visc(3)+v3*h_visc(4)-lambda*cv_gradTz

  END ASSOCIATE !rho,rhov1,rhov2,rhov3,rhoE & v_x/y/z,p_x/y/z ...
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
