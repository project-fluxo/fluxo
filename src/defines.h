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
!==================================================================================================================================
! Here, preprocessor variables for different equation systems and abbreviations for specific expressions are defined
!==================================================================================================================================

! Abbrevations
#ifdef SUN
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SUN COMPILER'
#  define IEEE_ISNAN
#elif SX
#  define __DATE__ '__TIME__ and __DATE__ not'
#  define __TIME__ 'available for SX COMPILER'
#elif PGI
#  define NO_ISNAN
#endif
#ifndef __FILENAME__ 
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#define SIZEOF_F(x) INT(STORAGE_SIZE(x)/8)

#ifdef GNU
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  k).OR.x<-HUGE(1_  k))       CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ k).OR.x<-HUGE(1._ k))       CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#elif CRAY
#define CHECKSAFEINT(x,k)  
#define CHECKSAFEREAL(x,k) 
#else
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  ## k).OR.x<-HUGE(1_  ## k)) CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ ## k).OR.x<-HUGE(1._ ## k)) CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#endif

#if MPI
#  define SWRITE IF(MPIRoot) WRITE
#  define IPWRITE(a,b) WRITE(a,b)myRank,
#  define GETTIME(a) a=MPI_WTIME()
#else
#  define SWRITE WRITE
#  define IPWRITE WRITE
#  define GETTIME(a) CALL CPU_TIME(a)
#endif
#define ERRWRITE(a,b) CALL CreateErrFile(); IF(ErrorFiles) WRITE(UNIT_errOut,b) 
#define LOGWRITE(a,b) IF(Logging) WRITE(UNIT_logOut,b)
#define SDEALLOCATE(A) IF(ALLOCATED(A)) DEALLOCATE(A)

! Loop variables
#ifdef OPTIMIZED
#define PP_IJK     i,0,0
#define PP_ij      i,0
#else
#define PP_IJK     i,j,k
#define PP_ij      i,j
#endif

! Predefined "PARAMETER-like" variables
#define XI_MINUS   5
#define XI_PLUS    3
#define ETA_MINUS  2
#define ETA_PLUS   4
#define ZETA_MINUS 1
#define ZETA_PLUS  6

! Entry position in SideToElem
#define S2E_ELEM_ID        1
#define S2E_NB_ELEM_ID     2
#define S2E_LOC_SIDE_ID    3
#define S2E_NB_LOC_SIDE_ID 4
#define S2E_FLIP           5

! Entry position in ElemToSide
#define E2S_SIDE_ID 1
#define E2S_FLIP    2

! Entry position in BC
#define MI_SIDEID 1
#define MI_FLIP   2

! Entry position in BC
#define BC_TYPE  1
#define BC_STATE 2
#define BC_ALPHA 3

! Entry position in BC
#define SEND 1
#define RECV 2

!#define DEBUGMESH

#define KILL(x) SWRITE(*,*) __FILE__,__LINE__,x; stop
