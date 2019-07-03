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
!> Subroutines defining one specific testcase with all necessary variables
!==================================================================================================================================
MODULE MOD_Testcase_ExactFunc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES

INTERFACE TestcaseExactFunc
  MODULE PROCEDURE TestcaseExactFunc
END INTERFACE


PUBLIC:: TestcaseExactFunc

CONTAINS

!==================================================================================================================================
!> Specifies all the initial conditions, all case numbers  here are >10000.
!==================================================================================================================================
SUBROUTINE TestcaseExactFunc(ExactFunction,tIn,x,resu,resu_t,resu_tt,Apot) 
! MODULES
USE MOD_Globals,ONLY:Abort,CROSS
USE MOD_Preproc
USE MOD_Equation_Vars,ONLY:RefStatePrim,IniRefState
USE MOD_Equation_Vars,ONLY:smu_0,mu_0
USE MOD_Equation_Vars,ONLY:IniFrequency,IniAmplitude,IniHalfwidth,IniWaveNumber
USE MOD_Equation_Vars,ONLY:IniDisturbance
USE MOD_Equation_Vars,ONLY:PrimToCons
USE MOD_Testcase_Vars,ONLY:EvalEquilibrium
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: ExactFunction    !< determines the exact function
REAL,INTENT(IN)                 :: tIn              !< current simulation time
REAL,INTENT(IN)                 :: x(3)             !< position in physical coordinates
REAL,INTENT(OUT)                :: Resu(PP_nVar)    !< exact fuction evaluated at tIn, returning state in conservative variables
REAL,INTENT(OUT)                :: Resu_t(PP_nVar)  !< first time deriv of exact fuction
REAL,INTENT(OUT)                :: Resu_tt(PP_nVar) !< second time deriv of exact fuction
REAL,INTENT(OUT),OPTIONAL       :: Apot(3)          !< optional magnetic vector potential
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                         :: i,j
REAL                            :: Omega,a,xc(3)
REAL                            :: Prim(1:PP_nVar) 
REAL                            :: r,r2(1:16),Bphi,dp
REAL                            :: q0,q1,Lz
REAL                            :: PsiN
REAL                            :: psi,dpsi_dx,dpsi_dy,psi_axis,psi_0
REAL                            :: rR2,zR2,sR0,Fprof,Faxis,dF2,p0
!==================================================================================================================================
Resu=0.
Resu_t=0.
Resu_tt=0.
! Determine the value, the first and the second time derivative
SELECT CASE (ExactFunction)
CASE(10010) ! mhd exact equilibrium, from potential A=(0,0,A3), A3=IniAmplitude*PRODUCT(sin(omega*x(:)))
         !domain should be a cube [0,1]^2, boundary conditions either periodic of perfectly conducting wall
  xc(:)  = x(:) + (/-0.0343, 0.0776, 0./)
  Prim(:)= RefStatePrim(IniRefState,:)
  Prim(2:4)=0.
  Omega=2*PP_Pi*IniFrequency
  a=SQRT(IniAmplitude*Prim(5))/omega !IniAmplitude is related to the change of pressure (IniAmplitude=0.1: 10% change) 
  Prim(6)= a*omega*SIN(Omega*xc(1))*COS(Omega*xc(2))
  Prim(7)=-a*omega*COS(Omega*xc(1))*SIN(Omega*xc(2))
  Prim(5)=Prim(5)*(1+ IniAmplitude*(SIN(Omega*xc(1))*SIN(Omega*xc(2)))**2) !a^2omega^2=p0*IniAmplitude

  CALL PrimToCons(Prim,Resu)
  IF(PRESENT(Apot))THEN
    Apot(1:2) = 0.
    Apot(3)   = a*SIN(Omega*xc(1))*SIN(Omega*xc(2)) 
  END IF
CASE(10013) ! 3D mhd exact equilibrium with constant pressure, 
         ! Potential A_1= (y^3*z-y*z^3), A_2=x^3*z-x*z^3, A_3= x^3*y-x*y^3
         ! B = [3*x*z^2-3*x*y^2, -3*y*z^2+2*y^3-3*x^2*y, 3*x^2*z-3*y^2*z]
         ! J=curl(B) = (0,0,0) -> dp=0, constant pressure
  Prim(:)= RefStatePrim(IniRefState,:)
  Prim(2:4)=0.
  Prim(6)=3.*x(1)*(x(3)**2-x(2)**2)
  Prim(7)=   x(2)*(2.*x(2)**2-3.*(x(3)**2+x(1)**2))
  Prim(8)=3.*x(3)*(x(1)**2-x(2)**2)

  CALL PrimToCons(Prim,Resu)
  IF(PRESENT(Apot))THEN
    Apot(1) = x(2)*x(3)*(x(2)**2-x(3)**2) 
    Apot(2) = x(1)*x(3)*(x(1)**2-x(3)**2) 
    Apot(3) = x(1)*x(2)*(x(1)**2-x(2)**2)
  END IF
CASE(10014) ! 3D mhd exact equilibrium with constant pressure & constant magnetic field
         ! Potential A = 1/2 [B_2*z-B_3*y ,B_3*x-B_1*z,, B_1*y-B_2*x] 
         ! B = [B_1,B_2,B_3] =const
         ! J=curl(B) = (0,0,0) -> dp=0, constant pressure
  Prim(:)= RefStatePrim(IniRefState,:)
  Prim(2:4)=0.
  Prim(6)=0.15
  Prim(7)=0.3
  Prim(8)=1.0

  CALL PrimToCons(Prim,Resu)
  IF(PRESENT(Apot))THEN
    Apot(1) = 0.5*(Prim(5+2)*x(3)-Prim(5+3)*x(2)) 
    Apot(2) = 0.5*(Prim(5+3)*x(1)-Prim(5+1)*x(3)) 
    Apot(3) = 0.5*(Prim(5+1)*x(2)-Prim(5+2)*x(1)) 
  END IF
CASE(10071) !Tearing mode instability, of paper Landi et al. , domain [0,6*pi]x[-pi/2,pi/2]x[0:2Pi]
        ! "Three-dimensional simulations of compressible tearing instability"
        ! rho_0=1, p0=0.5*beta (choose with refstate)
        ! Re_eta=5000, mu=0.,kappa=5/3 1/delta=0.1(=IniHalfwidth)  IniAmplitude=1.0E-04
  Prim=0.
  Prim(6)=TANH((x(2)-0.01)/IniHalfwidth) !tanh(y*delta) delta=10. !NOT FULLY CENTERED
  Prim(8)=SQRT(1-Prim(6)*Prim(6))
  
  Prim(1)=RefStatePrim(IniRefState,1) 
  IF(EvalEquilibrium)THEN
    Prim(3)=0.
  ELSE
    DO j=0,NINT(IniWaveNumber(3))
      DO i=0,NINT(IniWaveNumber(1))
        a=REAL(1+0.8*i+0.9*j)/(1+0.8*IniWaveNumber(1)+0.9*IniWaveNumber(3))
        Prim(3)=Prim(3)+SIN(x(1)/3.*i+ x(3)*j+2*PP_Pi*a)
      END DO
    END DO
    Prim(3)=IniDisturbance*Prim(6)*Prim(8)*Prim(3)
  END IF
  Prim(5)=RefStatePrim(IniRefState,5)
  !Prim(6:8)=sSqrt4pi*Prim(6:8) ! scaling with sqrt(4pi)!?!
  CALL PrimToCons(Prim,Resu)


CASE(10090) !cylindrical equilibrium for ideal MHD for current hole (Czarny, JCP, 2008), current Jz in z direction is given:
         ! cylindrical domain r[0,1], z[0,20] (from q(r=1)=4.4 =2*pi*B0/(Lz*Bphi(r=1)) => B0/Lz=0.364 B0~7.44, L0~20.)
         ! Jz=j1*(1-r^4)-j2*(1-r^2)^8, j1=0.2, j2=0.266
         ! from J=rotB (Br=0) follows
         ! Bphi(r) =mu_0 1/r  \int_0^r r*Jc dr
         ! pressure difference from gradp=J x B:
         ! dp(r) = -smu_0 ( 0.5*Bphi^2 + \int_0^r Bphi^2/r dr)
  Prim(:)=0.
  Prim(1)=1. 
  Prim(5)=0.025/9. ! =pmax, pmin=p(r=1)=pmax-0.0025 and pmax/pmin ~=10
  Prim(8)=SQRT(2*mu_0*Prim(5)*1.0E4) !beta=p/(0.5*smu0 B^2) =1.0E-04
  Prim(5)=Prim(5)+0.01  ! BETA CHANGED!!! ~ 1E-02  
  r2(1)=SUM(x(1:2)**2)
  r=SQRT(r2(1))
  DO i=2,16
    r2(i)=r2(i-1)*r2(1)
  END DO  
  Bphi=-mu_0*(133*r2(8)-1197*r2(7)+4788*r2(6)-11172*r2(5)+16758*r2(4)-16758*r2(3)+11472*r2(2)-4788*r2(1)+297)/9.0E+3 !*r
  dp  =-r2(1)*(      607086480.*r2(16)   -10965499545.*r2(15)   +93572262784.*r2(14)  -501279979200.*r2(13) &
                +1889439921600.*r2(12) -5321922445840.*r2(11)+11614288882560.*r2(10)-20096857548768.*r2( 9) &
               +27971578034400.*r2( 8)-31525903338930.*r2( 7)+28714629715200.*r2( 6)-20859299741280.*r2( 5) &
               +11747978409024.*r2( 4) -4854306857400.*r2( 3) +1285266977280.*r2( 2)  -138278780640.*r2( 1) &
               +5718295440.)/5.25096E+12
  IF(PRESENT(Apot))THEN
    Apot(1:2) = Prim(8)*0.5*(/-x(2),x(1)/)
    Apot(3)   = r2(1)*(    5320.*r2(8)  -53865.*r2(7) +246240.*r2(6)-670320.*r2(5) & 
                       +1206576.*r2(4)-1508220.*r2(3)+1376640.*r2(2)-861840.*r2(1)+106920.)/6.48E+6
  END IF
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  Prim(5)=Prim(5)+smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
  !disturbance of the velocity in xy plane  (scaled with (1-r^2)^2 and IniDisturbance) 
  IF(.NOT.EvalEquilibrium) THEN
    Prim(2)=IniDisturbance*(1.-r2(1))**2
    Prim(3)=x(2)*Prim(2)
  END IF
              
  CALL PrimToCons(Prim,Resu)
CASE(10091) !like 90, BUT WITH REAL BETA!!
         !cylindrical equilibrium for ideal MHD for current hole (Czarny, JCP, 2008), current Jz in z direction is given:
         ! cylindrical domain r[0,1], z[0,20] (from q(r=1)=4.4 =2*pi*B0/(Lz*Bphi(r=1)) => B0/Lz=0.364 B0~7.44, L0~20.)
         ! Jz=j1*(1-r^4)-j2*(1-r^2)^8, j1=0.2, j2=0.266
         ! from J=rotB (Br=0) follows
         ! Bphi(r) =mu_0 1/r  \int_0^r r*Jc dr
         ! pressure difference from gradp=J x B:
         ! dp(r) = -smu_0 ( 0.5*Bphi^2 + \int_0^r Bphi^2/r dr)
  Prim(:)=0.
  Prim(1)=1. 
  Prim(5)=0.025/9. ! =pmax, pmin=p(r=1)=pmax-0.0025 and pmax/pmin ~=10
  Prim(8)=SQRT(2*mu_0*Prim(5)*1.0E4) !beta=p/(0.5*smu0 B^2) =1.0E-04

  r2(1)=SUM(x(1:2)**2)
  r=SQRT(r2(1))
  DO i=2,16
    r2(i)=r2(i-1)*r2(1)
  END DO  
  Bphi=-mu_0*(133.*r2(8)-1197.*r2(7)+4788.*r2(6)-11172.*r2(5)+16758.*r2(4)-16758.*r2(3)+11472.*r2(2)-4788.*r2(1)+297.)/9.0E+3 !*r
  dp  =-r2(1)*(      607086480.*r2(16)   -10965499545.*r2(15)   +93572262784.*r2(14)  -501279979200.*r2(13) &
                +1889439921600.*r2(12) -5321922445840.*r2(11)+11614288882560.*r2(10)-20096857548768.*r2( 9) &
               +27971578034400.*r2( 8)-31525903338930.*r2( 7)+28714629715200.*r2( 6)-20859299741280.*r2( 5) &
               +11747978409024.*r2( 4) -4854306857400.*r2( 3) +1285266977280.*r2( 2)  -138278780640.*r2( 1) &
               +5718295440.)/5.25096E+12
  IF(PRESENT(Apot))THEN
    Apot(1:2) = 0. !Prim(8)*0.5*(/-x(2),x(1)/)
    Apot(3)   = r2(1)*(    5320.*r2(8)  -53865.*r2(7) +246240.*r2(6)-670320.*r2(5) & 
                       +1206576.*r2(4)-1508220.*r2(3)+1376640.*r2(2)-861840.*r2(1)+106920.)/6.48E+6
  END IF
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  Prim(5)=Prim(5)+smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
  !disturbance of the velocity in xy plane  (scaled with (1-r^2)^2 and IniDisturbance) 
  IF(.NOT.EvalEquilibrium) THEN
    Prim(2)=IniDisturbance*(1.-r2(1))**2
    Prim(3)=x(2)*Prim(2)
  END IF
              
  CALL PrimToCons(Prim,Resu)
CASE(10092) !like 90, BUT WITH REAL BETA!!
         !cylindrical equilibrium for ideal MHD for current hole (Czarny, JCP, 2008), current Jz in z direction is given:
         ! cylindrical domain r[0,1], z[0,20] (from q(r=1)=4.4 =2*pi*B0/(Lz*Bphi(r=1)) => B0/Lz=0.364 B0~7.44, L0~20.)
         ! Jz=q0*(1-r^4)-q1*(1-r^2)^8
  q0=0.2d0
  q1=0.d0 !0.266d0
         ! from J=rotB (Br=0) follows
         ! Bphi(r) =mu_0 1/r  \int_0^r r*Jc dr
         ! pressure difference from gradp=J x B:
         ! dp(r) = -smu_0 ( 0.5*Bphi^2 + \int_0^r Bphi^2/r dr)
  Prim(:)=0.d0
  Prim(1)=1.d0 
  dp = (68612544.d0*q0**2-40824446.d0*q0*q1+8842385.d0*q1**2)/4.4108064D+8 !dpmax=-dp(r=1)
  Prim(5)=10.d0*dp/9.d0 ! =pmax, pmin=p(r=1)=pmax-dp and pmax/pmin ~=10
  Prim(8)=SQRT(2.d0*mu_0*Prim(5)*1.0D4) !beta=p/(0.5*smu0 B^2) =1.0E-04

  r2(1)=SUM(x(1:2)**2)
  r=SQRT(r2(1))
  DO i=2,16
    r2(i)=r2(i-1)*r2(1)
  END DO  
  Bphi=-mu_0*(9.d0*(q1-q0)-36.d0*q1*r2(1)+(84.d0*q1+3.d0*q0)*r2(2)-126.d0*q1*r2(3) &
              +126.d0*q1*r2(4)-84.d0*q1*r2(5)+36.d0*q1*r2(6)-9.d0*q1*r2(7)+q1*r2(8))/18.d0 !*r 


  dp=-r2(1)*(       (110270160.d0*q1**2-220540320.d0*q0*q1+110270160.d0*q0**2)+(661620960.d0*q0*q1-661620960.d0*q1**2)*r2( 1) &
                                                      +(2548465920.d0*q1**2-1323241920.d0*q0*q1-49008960.d0*q0**2)*r2( 2) &
         +(1745944200.d0*q0*q1-7075668600.d0*q1**2)*r2( 3)+(15026147136.d0*q1**2-1440863424.d0*q0*q1+7351344.d0*q0**2)*r2( 4) &
         +(600359760.d0*q0*q1-25215109920.d0*q1**2)*r2( 5)                 +(34026220800.d0*q1**2+84015360.d0*q0*q1)*r2( 6) &
      +((-37229962770.d0*q1**2)-261891630.d0*q0*q1)*r2( 7)                +(33094661600.d0*q1**2+149749600.d0*q0*q1)*r2( 8) &
       +((-23828156352.d0*q1**2)-40432392.d0*q0*q1)*r2( 9)                  +(13784883840.d0*q1**2+4455360.d0*q0*q1)*r2(10)&
                             -6318071760.d0*q1**2*r2(11)     +2243102400.d0*q1**2*r2(12)       -595108800.d0*q1**2*r2(13)&
                         +111086976.d0*q1**2*r2(14)       -13018005.d0*q1**2*r2(15)          +720720.d0*q1**2*r2(16))/4.4108064D+8

  IF(PRESENT(Apot))THEN
    Apot(1:2) = Prim(8)*0.5d0*(/-x(2),x(1)/)
    Apot(3)   = r2(1)*((22680.d0*q1-22680.d0*q0)-45360.d0*q1*r2(1)+(70560.d0*q1+2520.*q0)*r2(2)-79380.d0*q1*r2(3) &
                           +63504.d0*q1*r2(4)-35280.d0*q1*r2(5)+12960.d0*q1*r2(6)-2835.d0*q1*r2(7)+280.d0*q1*r2(8))/9.072D+4
  END IF
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  Prim(5)=Prim(5)+smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
  !disturbance of the velocity in xy plane  (scaled with (1-r^2)^2 and IniDisturbance) 
  IF(.NOT.EvalEquilibrium) THEN
    Prim(2)=IniDisturbance*(1.-r2(1))**2
    Prim(3)=x(2)*Prim(2)
  END IF
              
  CALL PrimToCons(Prim,Resu)
CASE(10100) !solov'ev equilibrium (data from iphyton script...)
            !R0=10, a=3 !!! Raxis=10.57057107272790
  sR0=0.1 !R0=10.
  p0=2.0e-04
  Faxis=10.
  psi_axis=-6.178530464221542e-02
  psi_0 = 1.087903535041322e+01
  dp=-1.46250030038e-04
  dF2 =  -1.170000240304823e-01 

  r2(1)=SUM(x(1:2)**2)
  r  =SQRT(r2(1))
  rR2=r2(1)*sR0*sR0  !x^2=(R/R0)^2
  zR2=(x(3)*sR0)**2  !y^2=(Z/R0)^)2
  
  psi = rR2*(rR2*(-0.00323934432272746d0*rR2 + 0.0310403721993646d0*zR2 - &
              0.0143450850199111d0) + zR2*(-0.0112954698582052d0*zR2 + &
              0.0668243748907208d0) - 0.120048805732835d0) + rR2*(rR2*( &
              0.000783175967336492d0*rR2 - 0.00939811160803791d0*zR2 + &
              0.0150926608648206d0) + zR2*(0.00626540773869194d0*zR2 - &
              0.0603706434592823d0) + 0.12811145973562d0)*LOG(rR2) + zR2*(zR2*( &
              -0.000835387698492258d0*zR2 + 0.0201235478197608d0) + &
              0.143777080528761d0) + 0.0766840762309537d0

  psin=(psi-psi_axis)/psi_axis

  dpsi_dx = r*sR0*(rR2*(-0.0178697140016918d0*rR2 + 0.105365265581383d0*zR2 - &
                       0.0271950183500031d0) + zR2*(-0.0100601242390265d0*zR2 + &
                       0.0129074628628771d0) + (rR2*(0.00469905580401895d0*rR2 - &
                       0.0375924464321516d0*zR2 + 0.0603706434592823d0) + zR2*( &
                       0.0125308154773839d0*zR2 - 0.120741286918565d0) + &
                       0.256222919471239d0)*LOG(rR2) + 0.0161253080055692d0)
  
  dpsi_dy = x(3)*sR0*(rR2*(-0.0187962232160758d0*rR2 + 0.0250616309547677d0*zR2 - &
                          0.120741286918565d0)*LOG(rR2) + rR2*(0.0620807443987293d0*rR2 - &
                          0.0451818794328207d0*zR2 + 0.133648749781442d0) + zR2*( &
                          -0.00501232619095355d0*zR2 + 0.080494191279043d0) + &
                          0.287554161057522d0)
  !real psi scaled with psi_0
  dpsi_dx=dpsi_dx*psi_0 
  dpsi_dy=dpsi_dy*psi_0 
  psi    =    psi*psi_0

  Prim(:)=0.
  Prim(1)=1. 
  Prim(5)=p0+dp*psin !linear pressure profile
  Fprof=SQRT(Faxis**2+dF2*psin) !linear current profile

  IF(PRESENT(Apot))THEN
    !A_xy= Ar*r^ +Aphi*phi^ Ar=F/R*Z => Bphi=dAr/dz=F/R
    ! Aphi=psi/R
    Apot(1:2) = (Fprof*x(3)*(/x(1),x(2)/) + psi*(/-x(2),x(1)/) )/r2(1)
    Apot(3)   = 0.    
  END IF
  !BR=-1/R*dpsi/dZ = -1/R*dpsi/dy*(1/R0)
  !BZ=1/R*dpsi/dR = 1/R*dpsi/dx*(1/R0)
  !Bphi=F/r
  Prim(6)=(-x(2)*Fprof- dpsi_dy*x(1)*sR0)/r2(1)
  Prim(7)=( x(1)*Fprof- dpsi_dy*x(2)*sR0)/r2(1)
  Prim(8)= dpsi_dx*sR0/r
 
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
  !disturbance of the velocity in xy plane  (scaled with (1-r^2)^2 and IniDisturbance) 
  IF(.NOT.EvalEquilibrium) THEN
    Prim(2)=IniDisturbance*(1.-psin)**2
    Prim(3)=x(2)*Prim(2)
  END IF
              
  CALL PrimToCons(Prim,Resu)
CASE(92) !cylindrical equilibrium for ideal MHD for internal kink (Jorek paper Huysmans),q profile given, Bz=1
         ! cylindrical domain r[0,1], z[0,100]
         ! q(r)=0.73*(1-r^2)+1.6*r^2 => Bphi(r)=2*pi*r*Bz/(L0*q(r))
         ! pressure difference from gradp=J x B:
         ! dp(r) = smu_0 ( 0.5*Bphi^2 + \int_0^r Bphi^2/r dr)
         ! dp(r) = smu_0*2*pi^2*r^2/Lz^2*((q1-q0)*r^2+2*q0)/(q0*(((q1-q0)*r^2+2*q0)*(q1-q0)*r^2 + q0^2))
         ! p=p0-dp, p0=dp(1)/0.98, dp(1) = smu_0*2*pi^2/Lz^2*(q0+q1)/(q0*q1^2), so that pmax/pmin = 50
         ! small perturbation u_z=1.0E-12*r^2*(1-r^2)*cos(2*pi*z/Lz)
         !   in the first mode of the z velocity z \in[0,100] !! 
  r2(1)=SUM(x(1:2)**2)
  q0=0.73  !q(r=0)
  q1=1.60 !q(r=1)
  Lz=100.
  Prim(:)=0.
  Prim(1)=1.-0.9*r2(1) 
  Prim(4)=IniDisturbance*r2(1)*(1-r2(1))*SIN(2*PP_Pi*x(3)/Lz)
  Prim(5)=smu_0*2.*PP_Pi*PP_Pi*(q0+q1)/(0.98*Lz*Lz*q0*q1*q1)
  Prim(8)=1.  !Bz=1.
  Bphi= 2*PP_Pi/(Lz*((q1-q0)*r2(1)+q0)) !*r !L_0=100 !
  Prim(6)=-x(2)*Bphi  !/r
  Prim(7)= x(1)*Bphi  !/r
  dp  = 2*PP_Pi*PP_Pi*r2(1)*((q1-q0)*r2(1)+2*q0)/(Lz*Lz*q0*(((q1-q0)*r2(1)+2*q0)*(q1-q0)*r2(1)+q0**2))

  Prim(5)=Prim(5)-smu_0*dp
  IF(Prim(5).LT.0.) CALL abort(__STAMP__,  & 
                'negative pressure in exactfunc 90 !',999,prim(5))
              
  CALL PrimToCons(Prim,Resu)
CASE DEFAULT 
  CALL abort(__STAMP__, &
        " exactfunc does not exist in testcase")
END SELECT ! ExactFunction

END SUBROUTINE TestcaseExactFunc


END MODULE MOD_Testcase_ExactFunc
