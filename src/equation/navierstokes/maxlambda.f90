!==================================================================================================================================
! Copyright (c) 2020 - 2025 Andrés Rueda
! Copyright (c) 2016 - 2016 Guermond, J. L., & Popov, B.
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
!
! This module contains the routines to estimate the maximum wave speed from above according to:
!   Guermond, J. L., & Popov, B. (2016). Fast estimation from above of the maximum wave speed in the Riemann problem for the Euler equations. Journal of Computational Physics, 321, 908–926. https://doi.org/10.1016/j.jcp.2016.05.054
!
! Authors: Jean-Luc Guermond and Bojan Popov, Texas A&M, April 5, 2016
! Adapted to FLUXO by Andrés Rueda
!==================================================================================================================================
#include "defines.h"
MODULE MOD_MaxLambda
  USE MOD_Equation_Vars, ONLY: Kappa, KappaM1, KappaP1, sKappaM1, sKappa
  PRIVATE
  PUBLIC :: lambda, MaxLambdaInit
  
  REAL    :: Mgas   ! number of degrees of freedom of the molecules composing the gas (should be an integer, but we use real to allow any Kappa)
  REAL    :: expo
  REAL, PARAMETER :: b=0.0
CONTAINS

!==================================================================================================================================
!> Initialize module to estimate the maximum wave speed from above
!==================================================================================================================================
  SUBROUTINE MaxLambdaInit()
    USE MOD_Globals
    IMPLICIT NONE
    
    Mgas = 2.0*sKappaM1
    expo = 0.5*KappaM1*sKappa
    
  END SUBROUTINE MaxLambdaInit
!==================================================================================================================================
!> Estimate the maximum wave speed from above
!==================================================================================================================================
  pure SUBROUTINE lambda(tol,rhol,ul,pl,rhor,ur,pr,lambda_max,pstar,k)
    IMPLICIT NONE
    REAL, INTENT(IN) :: tol, rhol, ul, pl, rhor, ur, pr
    REAL, INTENT(OUT):: lambda_max, pstar
    INTEGER,      INTENT(OUT):: k
    REAL             :: capAmin, capBmin, acovmin, acovmax, ratio
    REAL             :: lambda_min, phimax, phimin, ptilde
    REAL             :: phi1, phi11, phi12, phi112, phi2, phi22, phi221
    REAL             :: p1, p2, pmin, pmax, rhomin, rhomax, v11, v12, v31, v32
    
    REAL             :: al, capAl, capBl, covl, ar, capAr, capBr, covr
    !===Initialization
    CALL GetAuxVars(rhol,pl,rhor,pr,al, capAl, capBl, covl, ar, capAr, capBr, covr)
    k = 0
    IF (pl.LE.pr) THEN
       pmin   = pl
       rhomin = rhol
       pmax   = pr
       rhomax = rhor
    ELSE
       pmin   = pr
       rhomin = rhor
       pmax   = pl
       rhomax = rhol
    END IF
    capAmin = 2/((KappaP1)*rhomin)
    capBmin = pmin*(KappaM1)/(KappaP1)
    acovmin = SQRT(Kappa*pmin*(1-b*rhomin)/rhomin)
    acovmax = SQRT(Kappa*pmax*(1-b*rhomax)/rhomax)
    ratio = (pmin/pmax)**expo
    phimin = (2.0*sKappaM1)*acovmax*(ratio-1.d0) + ur-ul
    IF (phimin.GE.0) THEN
       pstar = 0.d0
       lambda_max = MAX(MAX(-(ul-al),0.d0), MAX(ur+ar,0.d0))
       RETURN
    END IF
    phimax = (pmax-pmin)*SQRT(capAmin/(pmax+capBmin)) + ur-ul
    ptilde = pmin*((acovmin+acovmax - (ur-ul)*(KappaM1)*0.5)&
         /(acovmin + acovmax*ratio))**(Mgas+2)
    
    IF (phimax < 0.d0) THEN
       p1 = pmax
       p2 = ptilde
    ELSE
       p1=pmin
       p2 = MIN(pmax,ptilde)
    END IF
    !===Check for accuracy after initialization
    v11 = lambdaz(ul,pl,al/covl,p2,-1)
    v12 = lambdaz(ul,pl,al/covl,p1,-1)
    v31 = lambdaz(ur,pr,ar/covr,p1,1)
    v32 = lambdaz(ur,pr,ar/covr,p2,1)
    lambda_max = MAX(MAX(v32,0.d0),MAX(-v11,0.d0))
    lambda_min = MAX(MAX(MAX(v31,0.d0),MAX(-v12,0.d0)),0.d0)
    IF (lambda_min>0.d0) THEN
       IF (lambda_max/lambda_min -1.d0 .LE. tol) THEN
          pstar = p2
          RETURN
       END IF
    END IF
    p1 = MAX(p1,p2-phi(p2,ul,pl,ur,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr)/phi_prime(p2,pl,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr))
    !===Iterations
    DO WHILE(.TRUE.)
       v11 = lambdaz(ul,pl,al/covl,p2,-1)
       v12 = lambdaz(ul,pl,al/covl,p1,-1)
       v31 = lambdaz(ur,pr,ar/covr,p1,1)
       v32 = lambdaz(ur,pr,ar/covr,p2,1)
       lambda_max = MAX(MAX(v32,0.d0),MAX(-v11,0.d0))
       lambda_min = MAX(MAX(MAX(v31,0.d0),MAX(-v12,0.d0)),0.d0)
       IF (lambda_min>0.d0) THEN
          IF (lambda_max/lambda_min -1.d0 .LE. tol) THEN
             pstar = p2
             RETURN
          END IF
       END IF
       phi1 =  phi(p1,ul,pl,ur,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr) ! phi(p1,ul,pl,ur,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr)
       phi11 = phi_prime(p1,pl,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr)
       phi2 =  phi(p2,ul,pl,ur,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr)
       phi22 = phi_prime(p2,pl,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr)
       IF (phi1>0.d0) THEN
          lambda_max = lambda_min
          RETURN
       END IF
       IF (phi2<0.d0) RETURN
       phi12 = (phi2-phi1)/(p2-p1)
       phi112 = (phi12-phi11)/(p2-p1)
       phi221 = (phi22-phi12)/(p2-p1)
       p1 = p1 - 2*phi1/(phi11 + SQRT(phi11**2 - 4*phi1*phi112))
       p2 = p2 - 2*phi2/(phi22 + SQRT(phi22**2 - 4*phi2*phi221))
       k = k+1
    END DO
  END SUBROUTINE lambda
!==================================================================================================================================
!> Compute auxiliary variables to estimate the maximum wave speed from above
!==================================================================================================================================
  pure SUBROUTINE GetAuxVars(rhol,pl,rhor,pr,al, capAl, capBl, covl, ar, capAr, capBr, covr)
    IMPLICIT NONE
    REAL, INTENT(IN) :: rhol, pl, rhor, pr
    REAL, INTENT(OUT):: al, capAl, capBl, covl, ar, capAr, capBr, covr
    al = SQRT(Kappa*pl/rhol)
    capAl = 2/((KappaP1)*rhol)
    capBl = pl*(KappaM1)/(KappaP1)
    covl = SQRT(1-b*rhol)
    ar = SQRT(Kappa*pr/rhor)
    capAr = 2/((KappaP1)*rhor)
    capBr = pr*(KappaM1)/(KappaP1)
    covr = SQRT(1-b*rhor)
  END SUBROUTINE GetAuxVars
!==================================================================================================================================
!> Auxiliary function to estimate the maximum wave speed from above
!==================================================================================================================================
  pure FUNCTION lambdaz(uz,pz,az,pstar,z) RESULT(vv)
    IMPLICIT NONE
    REAL, INTENT(IN) :: uz,pz,az,pstar
    INTEGER,      INTENT(IN) :: z
    REAL             :: vv
    vv = uz + z*az*SQRT(1+MAX((pstar-pz)/pz,0.d0)*(KappaP1)/(2*Kappa))
  END FUNCTION lambdaz
!==================================================================================================================================
!> Auxiliary function to estimate the maximum wave speed from above
!==================================================================================================================================
  pure FUNCTION phi(p,ul,pl,ur,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr) RESULT(vv)
    IMPLICIT NONE
    REAL, INTENT(IN) :: p, ul, pl, ur, pr, capAl, capBl, al, capAr, capBr, ar, covl, covr
    REAL             :: vv, fl, fr
    IF (p>pl) THEN
       fl = (p-pl)*SQRT(capAl/(p+capBl))
    ELSE
       fl = (2*al/(KappaM1))*((p/pl)**expo-1)
    END IF
    IF (p>pr) THEN
       fr = (p-pr)*SQRT(capAr/(p+capBr))
    ELSE
       fr = (2*ar/(KappaM1))*((p/pr)**expo-1)
    END IF
    vv = fl*covl + fr*covr + ur - ul
  END FUNCTION phi
!==================================================================================================================================
!> Auxiliary function to estimate the maximum wave speed from above
!==================================================================================================================================
  pure FUNCTION phi_prime(p,pl,pr,capAl,capBl,al,capAr,capBr,ar,covl,covr) RESULT(vv)
    IMPLICIT NONE
    REAL, INTENT(IN) :: p, pl, pr, capAl, capBl, al, capAr, capBr, ar, covl, covr
    REAL             :: vv, fl, fr
    IF (p>pl) THEN
       fl = SQRT(capAl/(p+capBl))*(1-(p-pl)/(2*(capBl+p)))
    ELSE
       fl = (al/(Kappa*pl))*(p/pl)**(-(KappaP1)/(2*Kappa))
    END IF
    IF (p>pr) THEN
       fr = SQRT(capAr/(p+capBr))*(1-(p-pr)/(2*(capBr+p)))
    ELSE
       fr = (ar/(Kappa*pr))*(p/pr)**(-(KappaP1)/(2*Kappa))
    END IF
    vv = fl*covl + fr*covr
  END FUNCTION phi_prime
END MODULE MOD_MaxLambda
