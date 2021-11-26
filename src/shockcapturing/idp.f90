!==================================================================================================================================
! Copyright (c) 2020 - 2020 Andrés Rueda
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
! This module includes IDP correction routines that rely on the native LGL subcell FV discretization
!
!==================================================================================================================================
#include "defines.h"
#define barStates 1
module MOD_IDP
#if NFVSE_CORR
  implicit none
  
  private
  public :: DefineParameters_IDP
  public :: Init_IDP
  public :: Apply_IDP
  public :: Finalize_IDP
  
contains
!===================================================================================================================================
!> Defines the parameters of the IDP module
!===================================================================================================================================
  subroutine DefineParameters_IDP
    use MOD_ReadInTools,  only: prms
    implicit none
    
!   The different IDP limiters
!   --------------------------
    call prms%CreateLogicalOption("IDPMathEntropy",  " IDP correction on mathematical entropy?", "F")
    call prms%CreateLogicalOption("IDPSpecEntropy",  " IDP correction on specific entropy?", "F")
    call prms%CreateLogicalOption( "IDPDensityTVD",  " IDP(TVD) correction on density? (uses a Zalesak limiter with LOCAL_ALPHA=ON)", "F")
    call prms%CreateLogicalOption( "IDPPositivity",  " IDP correction for positivity of density and pressure?", "F")
    
!   Additional options
!   ------------------
    call prms%CreateRealOption(  "PositCorrFactor",  " The correction factor for IDPPositivity=T", "0.1")
    call prms%CreateIntOption(        "IDPMaxIter",  " Maximum number of iterations for positivity limiter", "10")
#if LOCAL_ALPHA
    call prms%CreateRealOption(         "IDPgamma",  " Constant for the subcell limiting of convex (nonlinear) constraints (must be IDPgamma>=2*d, where d is the number of dimensions of the problem)", "6.0")
#endif /*LOCAL_ALPHA*/
    
  end subroutine DefineParameters_IDP
!===================================================================================================================================
!> Initializes the IDP module
!===================================================================================================================================
  subroutine Init_IDP
    use MOD_NFVSE_Vars
    use MOD_IDP_Vars
    use MOD_Globals    , only: MPIRoot, UNIT_stdOut
    use MOD_PreProc    , only: PP_N
    use MOD_Mesh_Vars  , only: nElems
#if !(barStates)
    use MOD_Mesh_Vars  , only: firstSlaveSide, LastSlaveSide, nSides
#endif /*!(barStates)*/
    use MOD_ReadInTools, only: GETINT, GETREAL, GETLOGICAL
    implicit none
    !-local-variables----------------------------------------
    integer :: i
    !--------------------------------------------------------
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' IDP Methods: '
    
!   Get parameters
!   --------------
    
    ! IDP limiters
    IDPPositivity  = GETLOGICAL('IDPPositivity' ,'F')
    IDPDensityTVD  = GETLOGICAL('IDPDensityTVD' ,'F')
    IDPMathEntropy = GETLOGICAL('IDPMathEntropy','F')
    IDPSpecEntropy = GETLOGICAL('IDPSpecEntropy','F')
    
    ! Consistency check
    if (IDPMathEntropy .and. IDPSpecEntropy) then
      stop 'Only one of the two can be selected: IDPMathEntropy/IDPSpecEntropy'
    end if
    
    ! Additional options
    IDPMaxIter      = GETINT    ('IDPMaxIter','10')
    if (IDPPositivity) then
      PositCorrFactor = GETREAL   ('PositCorrFactor','0.1')
    end if
#if LOCAL_ALPHA
    if (IDPSpecEntropy .or. IDPMathEntropy .or. IDPPositivity) then
      IDPgamma = GETREAL   ('IDPgamma','6.0')
    end if
#endif /*LOCAL_ALPHA*/
    
!   Internal definitions (all are .FALSE. by default)
!   -------------------------------------------------
    if (IDPPositivity ) then
      IDPneedsUsafe     = .TRUE.
    end if
    
    if (IDPDensityTVD ) then
#if barStates
#if LOCAL_ALPHA
      IDPneedsUsafe = .TRUE.
#endif /*LOCAL_ALPHA*/
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev     = .TRUE.
      IDPneedsUprev_ext = .TRUE.
#else
      IDPneedsUsafe = .TRUE.
#endif /*barStates*/
    end if
    
    if (IDPMathEntropy) then
      IDPneedsUprev     = .TRUE.
#if barStates
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev_ext = .TRUE.
#endif /*barStates*/
    end if
    
    if (IDPSpecEntropy) then
      IDPneedsUprev     = .TRUE.
#if barStates
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev_ext = .TRUE.
#endif /*barStates*/
#if LOCAL_ALPHA
      IDPneedsUsafe = .TRUE.
#endif /*LOCAL_ALPHA*/
    end if
    
!   Allocate storage
!   ----------------
    ! Alpha before limiting
    allocate ( alpha_old(nElems) )
    alpha_old = 0.
    
    ! Udot_FV-Udot_DG
    allocate ( FFV_m_FDG    (PP_nVar, 0:PP_N  , 0:PP_N  , 0:PP_N  ,nElems) )
    FFV_m_FDG = 0.
    
    ! Usafe = U_FV (with external DOFs)
    if (IDPneedsUsafe) then
      allocate ( Usafe      (PP_nVar,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,nElems) )
      allocate ( p_safe             ( 0:PP_N  , 0:PP_N  , 0:PP_N  ,nElems) )
    end if
    
    ! Solution in the previous step (with external DOFs)
    if (IDPneedsUprev) then
      allocate ( Uprev      (PP_nVar,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,nElems) )
    end if
    
    ! Container for external Uprev
    if (IDPneedsUprev_ext) then
      allocate( Uprev_ext    (PP_nVar, 0:PP_N, 0:PP_N, 6,nElems) )
    end if
    
    ! Allocate bar states if needed
#if barStates
    if (IDPneedsUbar) then
      allocate( Ubar_xi      (PP_nVar,-1:PP_N  , 0:PP_N  , 0:PP_N  ,nElems) )
      allocate( Ubar_eta     (PP_nVar, 0:PP_N  ,-1:PP_N  , 0:PP_N  ,nElems) )
      allocate( Ubar_zeta    (PP_nVar, 0:PP_N  , 0:PP_N  ,-1:PP_N  ,nElems) )
    end if
#else
    ! Containers for previous entropy (with external DOFs)
    if (IDPMathEntropy .or. IDPSpecEntropy) then
      allocate( EntPrev       (     1,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1,nElems) )
      allocate( EntPrev_master(     1, 0:PP_N  , 0:PP_N            ,nSides) )
      allocate( EntPrev_slave (     1, 0:PP_N  , 0:PP_N            ,firstSlaveSide:LastSlaveSide) )
      allocate( EntPrev_ext   (     1, 0:PP_N  , 0:PP_N          ,6,nElems) )
    end if
    if (IDPDensityTVD) then
      allocate( Usafe_ext    (PP_nVar, 0:PP_N  , 0:PP_N          ,6,nElems) )
    end if
#endif /*barStates*/
    
    ! Variables for local alpha
#if LOCAL_ALPHA
    allocate ( alpha_loc(0:PP_N,0:PP_N,0:PP_N,1:nElems) )
    allocate ( ftilde_FV(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N,nElems) )
    allocate ( gtilde_FV(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N,nElems) )
    allocate ( htilde_FV(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N,nElems) )
    allocate ( ftilde_DG(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N,nElems) )
    allocate ( gtilde_DG(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N,nElems) )
    allocate ( htilde_DG(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N,nElems) )
    ftilde_FV = 0.0
    gtilde_FV = 0.0
    htilde_FV = 0.0
    ftilde_DG = 0.0
    gtilde_DG = 0.0
    htilde_DG = 0.0
    
    allocate ( dalpha_loc     (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1) )
#endif /*LOCAL_ALPHA*/
    
    ! Bounds containers
    if (IDPDensityTVD .or. IDPPositivity) then
      allocate ( rho_min     (0:PP_N,0:PP_N,0:PP_N) )
    end if
    if (IDPDensityTVD) then
      allocate ( rho_max     (0:PP_N,0:PP_N,0:PP_N) )
    end if
    if (IDPSpecEntropy) then
      allocate ( s_min       (0:PP_N,0:PP_N,0:PP_N) )
    end if
    if (IDPMathEntropy) then
      allocate ( s_max       (0:PP_N,0:PP_N,0:PP_N) )
    end if
    if (IDPPositivity) then
      allocate ( p_min       (0:PP_N,0:PP_N,0:PP_N) )
    end if
    
    ! Stencil for bounds
    ! ------------------
!#    allocate ( idx_p1(0:PP_N) )
!#    allocate ( idx_m1(0:PP_N) )
!#    !FV inner stencil
!#    idx_p1 = 1
!#    idx_p1(PP_N) = 0
    
!#    idx_m1 = -1
!#    idx_m1(0) = 0
    
!#    !FV stencil
!#    idx_p1 = 1
!#    idx_m1 = -1
    
!#    !DG stencil
!#    idx_p1 = (/(PP_N-i, i=0, PP_N)/)
!#    idx_m1 = (/(-i, i=0, PP_N)/)
!#    idx_p1(PP_N) = 1
!#    idx_m1(0)    =-1
    
!   Initialize variables for the analyze routines
!   ---------------------------------------------
    maximum_alpha = 0.0
    amount_alpha  = 0.0
    amount_alpha_steps = 0
    
!   All done
!   --------
#if barStates
    SWRITE(UNIT_stdOut,'(A)') ' INIT IDP DONE WITH BAR STATES'
#else
    SWRITE(UNIT_stdOut,'(A)') ' INIT IDP DONE WITHOUT BAR STATES'
#endif /*barStates*/
    
  end subroutine Init_IDP
!===================================================================================================================================
!> Applies all the activated IDP limiters.
!> This subroutine has to be called after after every SSP-RK stage
!===================================================================================================================================
  subroutine Apply_IDP(U,Ut,dt)
    use MOD_PreProc     , only: PP_N
    use MOD_Mesh_Vars   , only: nElems
    use MOD_NFVSE_Vars  , only: alpha, alpha_old, maximum_alpha, amount_alpha, amount_alpha_steps
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars  , only: alpha_loc
    use MOD_Analyze_Vars, only: wGPVol
    use MOD_Mesh_Vars   , only: nElems, sJ
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars    , only: IDPSpecEntropy, IDPMathEntropy, IDPDensityTVD, IDPPositivity
#if LOCAL_ALPHA
    use MOD_IDP_Vars    , only: dalpha_loc
#endif /*LOCAL_ALPHA*/
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    !-local-variables----------------------------------------
    real :: sdt
    real :: curr_amount_alpha, dalpha
    integer :: eID, i, j, k
    !--------------------------------------------------------

!   Initialize
!   ----------
    alpha_old = alpha
    sdt = 1./dt
    
!   Get the variables for limiting in the right position
!   ----------------------------------------------------
    call Get_IDP_Variables(U,dt)
    
!   Perform limiting!
!   -----------------
    do eID=1, nElems
#if LOCAL_ALPHA
      dalpha_loc = 0.0
#endif /*LOCAL_ALPHA*/
!     Call all user-defined limiters
!     (Element-wise limiting is performed directly in the routines)
!     -------------------------------------------------------------
      if (IDPDensityTVD)  call IDP_LimitDensityTVD (U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      if (IDPSpecEntropy) call IDP_LimitSpecEntropy(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      if (IDPMathEntropy) call IDP_LimitMathEntropy(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      if (IDPPositivity)  call IDP_LimitPositivity (U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      
#if LOCAL_ALPHA
!     Perform subcell-wise limiting
!     -----------------------------
      dalpha = maxval (dalpha_loc)
      if (dalpha>0.0 .or. any(isnan(dalpha_loc))) then
        call PerformCorrection(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dalpha    ,alpha(eID)          , &
                                      dalpha_loc,alpha_loc(:,:,:,eID), &
                                      dt,sdt,eID)
      end if
#endif /*LOCAL_ALPHA*/
!     Check that we are within bounds
!     -------------------------------
#if DEBUG || IDP_CHECKBOUNDS
      call CheckBounds(U(:,:,:,:,eID))
#endif /*DEBUG || IDP_CHECKBOUNDS*/
    end do
    
!   Update variables for the analyze routines
!   -----------------------------------------
    maximum_alpha = max(maximum_alpha,maxval(alpha-alpha_old))
    
    amount_alpha = amount_alpha*amount_alpha_steps
#if LOCAL_ALPHA
    curr_amount_alpha = 0.0
    do eID=1, nElems
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        curr_amount_alpha = curr_amount_alpha + (alpha_loc(i,j,k,eID)-alpha_old(eID))*wGPVol(i,j,k)/sJ(i,j,k,eID)
      end do       ; end do       ; end do
    end do
    amount_alpha = amount_alpha + curr_amount_alpha
#else
    amount_alpha = amount_alpha + sum(alpha-alpha_old)
#endif /*LOCAL_ALPHA*/
    amount_alpha_steps = amount_alpha_steps+1
    amount_alpha = amount_alpha/amount_alpha_steps
    
  end subroutine Apply_IDP
!===================================================================================================================================
!> Check that all bounds are met
!===================================================================================================================================
  subroutine CheckBounds(U)
    use MOD_Preproc
    use MOD_IDP_Vars      , only: rho_min, rho_max, s_min, s_max, p_min
    use MOD_IDP_Vars      , only: IDPDensityTVD, IDPSpecEntropy, IDPMathEntropy, IDPPositivity
    use MOD_Equation_Vars , only: Get_SpecEntropy, Get_MathEntropy, Get_Pressure
    implicit none
    !-arguments------------------------------------------------------------
    real, intent(in) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N)
    !-local-variables------------------------------------------------------
    integer :: i,j,k
    real    :: p
    real, parameter :: tolerance=1.e-12
    !----------------------------------------------------------------------
    
    
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      if (IDPDensityTVD .or. IDPPositivity) then
        if (U(1,i,j,k) < rho_min(i,j,k) - tolerance) then
          print*, 'WARNING: rho below min (curr/min):', U(1,i,j,k), rho_min(i,j,k)
          stop
        end if
      end if
        
      if (IDPDensityTVD) then
        if (U(1,i,j,k) > rho_max(i,j,k) + tolerance) then
          print*, 'WARNING: rho above max (curr/max):', U(1,i,j,k), rho_max(i,j,k)
          stop
        end if
      end if
      
      if (IDPSpecEntropy) then
        if (Get_SpecEntropy(U(:,i,j,k)) < s_min(i,j,k) - tolerance) then
          print*, 'WARNING: specific entropy below min (curr/min):', Get_SpecEntropy(U(:,i,j,k)), s_min(i,j,k)
          stop
        end if
      end if
      
      if (IDPMathEntropy) then
        if (Get_MathEntropy(U(:,i,j,k)) > s_max(i,j,k) + tolerance) then
          print*, 'WARNING: mathematical entropy above max (curr/max):', Get_MathEntropy(U(:,i,j,k)), s_max(i,j,k)
          stop
        end if
      end if
      
      if (IDPPositivity) then
        call Get_Pressure(U(:,i,j,k),p)
        if (p < p_min(i,j,k) - tolerance) then
          print*, 'WARNING: rho below min (curr/min):', p, p_min(i,j,k)
          stop
        end if
      end if
    end do       ; end do       ; end do ! i,j,k
  
  end subroutine CheckBounds
!===================================================================================================================================
!> Get the IDP variables in the right position to perform limiting
!> ATTENTION: 1) U_master and U_slave need to have the previous solution!
!>            2) If the bar states are *deactivated* and IDPDensityTVD, U_master and U_slave contain the safe solution when this 
!>               subroutine returns
!===================================================================================================================================
  subroutine Get_IDP_Variables(U,dt)
    use MOD_NFVSE_Vars, only: alpha
    use MOD_PreProc       , only: PP_N
    use MOD_IDP_Vars      , only: IDPneedsUprev_ext, IDPneedsUsafe
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: Uprev_ext,Uprev
    use MOD_IDP_Vars      , only: Usafe, p_safe
#if barStates
    use MOD_IDP_Vars      , only: IDPneedsUbar
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta
    use MOD_NFVSE_Vars    , only: SubCellMetrics
#else
    use MOD_ProlongToFace , only: ProlongToFace
    use MOD_IDP_Vars      , only: Usafe_ext
    use MOD_Mesh_Vars     , only: nSides, firstSlaveSide, lastSlaveSide
    use MOD_IDP_Vars      , only: EntPrev, EntPrev_master, EntPrev_slave, EntPrev_ext, IDPDensityTVD, IDPMathEntropy, IDPSpecEntropy
#if MPI
    use MOD_MPI_Vars      , only: MPIRequest_U, DataSizeSide, nNbProcs
    use MOD_NFVSE_Vars    , only: MPIRequest_Umaster
    use MOD_MPI           , only: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
#endif /*MPI*/
#endif /*barStates*/
    use MOD_DG_Vars       , only: U_master,U_slave
    use MOD_Mesh_Vars     , only: nElems, offsetElem
    use MOD_NFVSE_MPI     , only: Get_externalU
    use MOD_Equation_Vars , only: Get_MathEntropy, Get_SpecEntropy, Get_Pressure
    implicit none
    !-arguments------------------------------------------------------------
    real, intent(in) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
    real, intent(in) :: dt
    !-local-variables------------------------------------------------------
    integer :: i,j,k
    integer :: eID, sideID
    !----------------------------------------------------------------------
    
!   Get the previous solution in place if needed!
!   *********************************************
    if (IDPneedsUprev_ext) then
      ! Gather the external Uprev in the right location
      call Get_externalU(PP_nVar,Uprev_ext,Uprev(:,0:PP_N,0:PP_N,0:PP_N,:),U_master,U_slave)
      ! FIll Uprev with info
      do eID=1, nElems
        Uprev(:,    -1,0:PP_N,0:PP_N,eID) = Uprev_ext(:,0:PP_N,0:PP_N,5,eID)
        Uprev(:,PP_N+1,0:PP_N,0:PP_N,eID) = Uprev_ext(:,0:PP_N,0:PP_N,3,eID)
        Uprev(:,0:PP_N,    -1,0:PP_N,eID) = Uprev_ext(:,0:PP_N,0:PP_N,2,eID)
        Uprev(:,0:PP_N,PP_N+1,0:PP_N,eID) = Uprev_ext(:,0:PP_N,0:PP_N,4,eID)
        Uprev(:,0:PP_N,0:PP_N,    -1,eID) = Uprev_ext(:,0:PP_N,0:PP_N,1,eID)
        Uprev(:,0:PP_N,0:PP_N,PP_N+1,eID) = Uprev_ext(:,0:PP_N,0:PP_N,6,eID)
      end do !eID
    end if !allocated(Uprev_ext)
    
!   Get the bar states if needed!
!   *****************************
#if barStates
    if (IDPneedsUbar) then
      do eID=1, nElems
        associate (SCM => SubCellMetrics(eID))
        do i=-1, PP_N ! i goes through the interfaces
          do k=0, PP_N  ; do j=0, PP_N ! j and k go through DOFs
            !xi
            Ubar_xi  (:,i,j,k,eID) = GetBarStates(Uprev(:,i,j,k,eID),Uprev(:,i+1,j,k,eID),SCM % xi   % nv (:,j,k,i),SCM % xi   % t1 (:,j,k,i),SCM % xi   % t2 (:,j,k,i))
            !eta
            Ubar_eta (:,j,i,k,eID) = GetBarStates(Uprev(:,j,i,k,eID),Uprev(:,j,i+1,k,eID),SCM % eta  % nv (:,j,k,i),SCM % eta  % t1 (:,j,k,i),SCM % eta  % t2 (:,j,k,i))
            !zeta
            Ubar_zeta(:,j,k,i,eID) = GetBarStates(Uprev(:,j,k,i,eID),Uprev(:,j,k,i+1,eID),SCM % zeta % nv (:,j,k,i),SCM % zeta % t1 (:,j,k,i),SCM % zeta % t2 (:,j,k,i))
          end do        ; end do  ! j,k
        end do
        end associate
      end do !eID
    end if
#else
!   Otherwise get the previous entropy if needed
!   ********************************************
    ! Mathematical entropy
    if (IDPMathEntropy) then
      do eID=1, nElems
        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          EntPrev(1,i,j,k,eID) = Get_MathEntropy(Uprev(:,i,j,k,eID))
        end do       ; end do       ; end do ! i,j,k
      end do !eID
      do sideID=1, nSides
        do j=0, PP_N ; do i=0, PP_N
          EntPrev_master(1,i,j,sideID) = Get_MathEntropy(U_master(:,i,j,sideID))
        end do       ; end do ! i,j
      end do !sideID
      do sideID=firstSlaveSide, LastSlaveSide
        do j=0, PP_N ; do i=0, PP_N
          EntPrev_slave (1,i,j,sideID) = Get_MathEntropy(U_slave (:,i,j,sideID))
        end do       ; end do ! i,j
      end do !sideID
    ! Specific entropy
    elseif (IDPSpecEntropy) then
      do eID=1, nElems
        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          EntPrev(1,i,j,k,eID) = Get_SpecEntropy(Uprev(:,i,j,k,eID))
        end do       ; end do       ; end do ! i,j,k
      end do !eID
      do sideID=1, nSides
        do j=0, PP_N ; do i=0, PP_N
          EntPrev_master(1,i,j,sideID) = Get_SpecEntropy(U_master(:,i,j,sideID))
        end do       ; end do ! i,j
      end do !sideID
      do sideID=firstSlaveSide, LastSlaveSide
        do j=0, PP_N ; do i=0, PP_N
          EntPrev_slave (1,i,j,sideID) = Get_SpecEntropy(U_slave (:,i,j,sideID))
        end do       ; end do ! i,j
      end do !sideID
    end if
    
    if (IDPMathEntropy .or. IDPSpecEntropy) then
      ! Get neighbor info
      ! *****************
    
      ! Gather the external math entropy in the right location
      call Get_externalU(1,EntPrev_ext,EntPrev(:,0:PP_N,0:PP_N,0:PP_N,:),EntPrev_master,EntPrev_slave)
      ! FIll EntPrev with info
      do eID=1, nElems
        EntPrev(1,    -1,0:PP_N,0:PP_N,eID) = EntPrev_ext(1,0:PP_N,0:PP_N,5,eID)
        EntPrev(1,PP_N+1,0:PP_N,0:PP_N,eID) = EntPrev_ext(1,0:PP_N,0:PP_N,3,eID)
        EntPrev(1,0:PP_N,    -1,0:PP_N,eID) = EntPrev_ext(1,0:PP_N,0:PP_N,2,eID)
        EntPrev(1,0:PP_N,PP_N+1,0:PP_N,eID) = EntPrev_ext(1,0:PP_N,0:PP_N,4,eID)
        EntPrev(1,0:PP_N,0:PP_N,    -1,eID) = EntPrev_ext(1,0:PP_N,0:PP_N,1,eID)
        EntPrev(1,0:PP_N,0:PP_N,PP_N+1,eID) = EntPrev_ext(1,0:PP_N,0:PP_N,6,eID)
      end do !eID
    end if
#endif /*barStates*/
    
!   Get Usafe for each point and check validity
!   *******************************************
    if (IDPneedsUsafe) then
      
      do eID=1, nElems
        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          
          ! Compute Usafe
          Usafe(:,i,j,k,eID) = U(:,i,j,k,eID) + dt * FFV_m_FDG(:,i,j,k,eID) * (1. - alpha(eID))
          
          ! Check if this is a valid state
          call Get_Pressure(Usafe(:,i,j,k,eID),p_safe(i,j,k,eID))
          
          if (p_safe(i,j,k,eID) < 0.) then
            print*, 'ERROR: safe pressure not safe el=', eID+offsetElem, p_safe(i,j,k,eID)
            stop
          end if
          if (Usafe(1,i,j,k,eID) < 0.) then
            print*, 'ERROR: safe dens not safe el=', eID+offsetElem, Usafe(1,i,j,k,eID)
            stop
          end if
        end do       ; end do       ; end do ! i,j,k
        
      end do !eID
      
    end if
    
#if !(barStates)
    if (IDPDensityTVD) then
!     Get the safe (FV) solution in place (if not using bar states)
!     * This overwrites U_master and U_slave
!     * The MPI communication here is blocking!!
!     *************************************************************
#if MPI
      ! receive the slave
      CALL StartReceiveMPIData(U_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide, &
                               MPIRequest_U(:,SEND),SendID=2) ! Receive MINE (sendID=2) 
      
      ! prolong MPI sides and do the mortar on the MPI sides
      CALL ProlongToFace(PP_nVar,Usafe(:,0:PP_N,0:PP_N,0:PP_N,:),U_master,U_slave,doMPISides=.TRUE.)
      
      ! TODO: Mortars are not really working yet!!
      !CALL U_Mortar(U_master,U_slave,doMPISides=.TRUE.)
      
      ! send the slave
      CALL StartSendMPIData(U_slave,DataSizeSide,FirstSlaveSide,LastSlaveSide, &
                            MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR (sendID=2) 
      
      ! receive the master
      call StartReceiveMPIData(U_master, DataSizeSide, 1, nSides, &
                               MPIRequest_Umaster(:,1), SendID=1) ! Receive YOUR  (sendID=1) 
      
#endif /* MPI */
      ! Get external Usafe
      call ProlongToFace(PP_nVar,Usafe(:,0:PP_N,0:PP_N,0:PP_N,:),U_master,U_slave,doMPISides=.FALSE.)
      ! TODO: Add mortars!!
#if MPI
      ! Send the master
      call StartSendMPIData   (U_master, DataSizeSide, 1, nSides, &
                               MPIRequest_Umaster(:,2),SendID=1) 
      
      
      call FinishExchangeMPIData(2*nNbProcs,MPIRequest_U) 
      call FinishExchangeMPIData(2*nNbProcs,MPIRequest_Umaster) 
      
#endif /* MPI */
      ! Gather the external Usafe in the right location
      call Get_externalU(PP_nVar,Usafe_ext,Usafe(:,0:PP_N,0:PP_N,0:PP_N,:),U_master,U_slave)
      ! FIll Usafe with info
      do eID=1, nElems
        Usafe(:,    -1,0:PP_N,0:PP_N,eID) = Usafe_ext(:,0:PP_N,0:PP_N,5,eID)
        Usafe(:,PP_N+1,0:PP_N,0:PP_N,eID) = Usafe_ext(:,0:PP_N,0:PP_N,3,eID)
        Usafe(:,0:PP_N,    -1,0:PP_N,eID) = Usafe_ext(:,0:PP_N,0:PP_N,2,eID)
        Usafe(:,0:PP_N,PP_N+1,0:PP_N,eID) = Usafe_ext(:,0:PP_N,0:PP_N,4,eID)
        Usafe(:,0:PP_N,0:PP_N,    -1,eID) = Usafe_ext(:,0:PP_N,0:PP_N,1,eID)
        Usafe(:,0:PP_N,0:PP_N,PP_N+1,eID) = Usafe_ext(:,0:PP_N,0:PP_N,6,eID)
      end do !eID
    end if !IDPDensityTVD
#endif /*!(barStates)*/
    
  end subroutine Get_IDP_Variables
!===================================================================================================================================
!> Density TVD correction
!===================================================================================================================================
  subroutine IDP_LimitDensityTVD(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_Mesh_Vars     , only: nElems
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_NFVSE_Vars    , only: ftilde_DG, gtilde_DG, htilde_DG, ftilde_FV, gtilde_FV, htilde_FV
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ, offsetElem
    use MOD_IDP_Vars      , only: Usafe, dalpha_loc
#endif /*LOCAL_ALPHA*/
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#else
    use MOD_IDP_Vars      , only: Usafe
#endif /*barStates*/
    use MOD_IDP_Vars      , only: FFV_m_FDG, rho_min, rho_max
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    real    :: dalpha, dalpha1
    real    :: rho_safe
    real    :: a   ! a  = PositCorrFactor * rho_safe - rho
    real    :: Qp, Qm, Pp, Pm
    integer :: i,j,k,l
    real, parameter :: eps = 1.e-10           ! Very small value
    !--------------------------------------------------------
    
      dalpha = -epsilon(1.0) ! Safe initialization
!       Compute correction factors
!       --------------------------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          
        ! Get the limit states
        !*********************
        rho_min(i,j,k) =  huge(1.0)
        rho_max(i,j,k) = -huge(1.0)
        
#if barStates
        ! Previous sol
        rho_min(i,j,k) = min(rho_min(i,j,k), Uprev  (1,i  ,j  ,k  ,eID))
        rho_max(i,j,k) = max(rho_max(i,j,k), Uprev  (1,i  ,j  ,k  ,eID))
        
        !xi
        do l=i-1, i !min(i-1,0), max(i,PP_N-1)
          rho_min(i,j,k) = min(rho_min(i,j,k), Ubar_xi  (1,l  ,j  ,k  ,eID))
          rho_max(i,j,k) = max(rho_max(i,j,k), Ubar_xi  (1,l  ,j  ,k  ,eID))
        end do
        !eta
        do l=j-1, j !l=max(j-1,0), min(j,PP_N-1)
          rho_min(i,j,k) = min(rho_min(i,j,k), Ubar_eta (1,i  ,l  ,k  ,eID))
          rho_max(i,j,k) = max(rho_max(i,j,k), Ubar_eta (1,i  ,l  ,k  ,eID))
        end do
        !zeta
        do l=k-1, k !l=max(k-1,0), min(k,PP_N-1)
          rho_min(i,j,k) = min(rho_min(i,j,k), Ubar_zeta(1,i  ,j  ,l  ,eID))
          rho_max(i,j,k) = max(rho_max(i,j,k), Ubar_zeta(1,i  ,j  ,l  ,eID))
        end do
#else
        ! check stencil in xi
!          do l = i+idx_m1(i), i+idx_p1(i) !no neighbor
        do l = i-1, i+1
          rho_min(i,j,k) = min(rho_min(i,j,k), Usafe(1,l,j,k,eID))
          rho_max(i,j,k) = max(rho_max(i,j,k), Usafe(1,l,j,k,eID))
        end do
        ! check stencil in eta
!          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          rho_min(i,j,k) = min(rho_min(i,j,k), Usafe(1,i,l,k,eID))
          rho_max(i,j,k) = max(rho_max(i,j,k), Usafe(1,i,l,k,eID))
        end do
        ! check stencil in zeta
!          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          rho_min(i,j,k) = min(rho_min(i,j,k), Usafe(1,i,j,l,eID))
          rho_max(i,j,k) = max(rho_max(i,j,k), Usafe(1,i,j,l,eID))
        end do
#endif /*barStates*/
        
#if LOCAL_ALPHA
        ! Real Zalesak type limiter
        ! * Zalesak (1979). "Fully multidimensional flux-corrected transport algorithms for fluids"
        ! * Kuzmin et al. (2010). "Failsafe flux limiting and constrained data projections for equations of gas dynamics"
        ! ATTENTION: 1) The Zalesak limiter has to be computed, even if the state is valid, because the correction is 
        !               for each interface, not each node
        !****************************************************************************************************************
        
        ! Upper/lower bounds for admissible increments
        Qp = max(0.0,(rho_max(i,j,k)-Usafe(1,i,j,k,eID))*sdt)
        Qm = min(0.0,(rho_min(i,j,k)-Usafe(1,i,j,k,eID))*sdt)
        
        ! Check the sign of Qp and Qm... !!!
        if (Qp < 0.0) then
          print*, 'Usafe is not within bounds: 0>Qp=', Qp, i, j, k, eID+offsetElem
          print*, 'rho_max', rho_max(i,j,k)
          print*, 'rho_saf', Usafe(1,i,j,k,eID)
          stop 
        end if
        if (Qm > 0.0)  then
          print*, 'Usafe is not within bounds: 0<Qm=', Qm, i, j, k, eID+offsetElem
          print*, 'rho_min', rho_min(i,j,k)
          print*, 'rho_saf', Usafe(1,i,j,k,eID)
          stop
        end if
        
        ! Positive contributions
        Pp = 0.0
        Pp = Pp + max(0.0, sWGP(i) * (ftilde_DG(1,i-1,j  ,k  ,eID) - ftilde_FV(1,i-1,j  ,k  ,eID)) )
        Pp = Pp + max(0.0,-sWGP(i) * (ftilde_DG(1,i  ,j  ,k  ,eID) - ftilde_FV(1,i  ,j  ,k  ,eID)) )
        Pp = Pp + max(0.0, sWGP(j) * (gtilde_DG(1,i  ,j-1,k  ,eID) - gtilde_FV(1,i  ,j-1,k  ,eID)) )
        Pp = Pp + max(0.0,-sWGP(j) * (gtilde_DG(1,i  ,j  ,k  ,eID) - gtilde_FV(1,i  ,j  ,k  ,eID)) )
        Pp = Pp + max(0.0, sWGP(k) * (htilde_DG(1,i  ,j  ,k-1,eID) - htilde_FV(1,i  ,j  ,k-1,eID)) )
        Pp = Pp + max(0.0,-sWGP(k) * (htilde_DG(1,i  ,j  ,k  ,eID) - htilde_FV(1,i  ,j  ,k  ,eID)) )
        Pp = Pp*sJ(i,j,k,eID)
        
        ! Negative contributions
        Pm = 0.0
        Pm = Pm + min(0.0, sWGP(i) * (ftilde_DG(1,i-1,j  ,k  ,eID) - ftilde_FV(1,i-1,j  ,k  ,eID)) )
        Pm = Pm + min(0.0,-sWGP(i) * (ftilde_DG(1,i  ,j  ,k  ,eID) - ftilde_FV(1,i  ,j  ,k  ,eID)) )
        Pm = Pm + min(0.0, sWGP(j) * (gtilde_DG(1,i  ,j-1,k  ,eID) - gtilde_FV(1,i  ,j-1,k  ,eID)) )
        Pm = Pm + min(0.0,-sWGP(j) * (gtilde_DG(1,i  ,j  ,k  ,eID) - gtilde_FV(1,i  ,j  ,k  ,eID)) )
        Pm = Pm + min(0.0, sWGP(k) * (htilde_DG(1,i  ,j  ,k-1,eID) - htilde_FV(1,i  ,j  ,k-1,eID)) )
        Pm = Pm + min(0.0,-sWGP(k) * (htilde_DG(1,i  ,j  ,k  ,eID) - htilde_FV(1,i  ,j  ,k  ,eID)) )
        Pm = Pm*sJ(i,j,k,eID)
        
        if (Pp==0.0) then
          Qp = 1.0
        else
          Qp = abs(Qp/Pp)
        end if
        
        if (Pm==0.0) then
          Qm = 1.0
        else
          Qm = abs(Qm/Pm)
        end if
        
        ! Compute correction as: (needed_alpha) - current_alpha = (1.0 - min(1.0,Qp,Qm)) - alpha_loc(i,j,k,eID)
        dalpha1 = 1.0 - min(1.0,Qp,Qm) - alpha_loc(i,j,k,eID)
        
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k),dalpha1)
        dalpha = max(dalpha,dalpha1)
        
#else
        ! Simple element-wise limiter
        !****************************
        if ( U(1,i,j,k) < rho_min(i,j,k)) then
          rho_safe = rho_min(i,j,k)
        elseif (U(1,i,j,k) > rho_max(i,j,k)) then
          rho_safe = rho_max(i,j,k)
        else
          cycle !nothing to do here!
        end if
        
        if ( abs(FFV_m_FDG(1,i,j,k,eID)) == 0.0) cycle !nothing to do here!
        
        ! Density correction
        a = (rho_safe - U(1,i,j,k))
        dalpha1 = a*sdt / FFV_m_FDG(1,i,j,k,eID)
        dalpha = max(dalpha,dalpha1)
        
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; end do ! i,j,k

!       Do the correction if needed
!       ---------------------------
#if !(LOCAL_ALPHA)
      if ( dalpha > 0. ) then
        call PerformCorrection(U,Ut,dalpha    ,alpha(eID)          , &
                                    dt,sdt,eID)
      end if
#endif /*!(LOCAL_ALPHA)*/
  end subroutine IDP_LimitDensityTVD
!===================================================================================================================================
!> Specific entropy correction (discrete local minimum principle)
!===================================================================================================================================
  subroutine IDP_LimitSpecEntropy(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: Usafe, dalpha_loc
    use MOD_NFVSE_Vars    , only: ftilde_DG, gtilde_DG, htilde_DG, ftilde_FV, gtilde_FV, htilde_FV
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ
#endif /*LOCAL_ALPHA*/
    use MOD_Mesh_Vars     , only: nElems, offsetElem
    use MOD_Equation_Vars , only: Get_SpecEntropy, ConsToSpecEntropy
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#else
    use MOD_IDP_Vars      , only: EntPrev
#endif /*barStates*/
    use MOD_IDP_Vars      , only: NEWTON_ABSTOL, NEWTON_RELTOL
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: alpha_maxIDP, IDPMaxIter, s_min
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    real    :: dalpha, dalpha1, dalpha_old
    real    :: as, U_curr(PP_nVar), dSdalpha
    integer :: i,j,k, l
    integer :: iter
    logical :: notInIter
    real, parameter :: eps = 1.e-14           ! Very small value
    !--------------------------------------------------------
      
      dalpha = -epsilon(1.0) ! Safe initialization
      
!       Compute correction factors
!       --------------------------
      notInIter = .FALSE.
      
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Get the limit states
        !*********************
        s_min(i,j,k) = huge(1.0)
        
#if barStates
        ! Previous entropy of the node (ubar_ii)
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Uprev    (:,i  ,j  ,k  ,eID)))
        
        ! TODO: Compute them for all interfaces before...
        !xi+
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Ubar_xi  (:,i  ,j  ,k  ,eID)))
        !xi-
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Ubar_xi  (:,i-1,j  ,k  ,eID)))
        !eta+
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Ubar_eta (:,i  ,j  ,k  ,eID)))
        !eta-
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Ubar_eta (:,i  ,j-1,k  ,eID)))
        !zeta+
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Ubar_zeta(:,i  ,j  ,k  ,eID)))
        !zeta-
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Ubar_zeta(:,i  ,j  ,k-1,eID)))
#else
        ! check stencil in xi
!#          do l = i+idx_m1(i), i+idx_p1(i) !no neighbor
        do l = i-1, i+1
          s_min(i,j,k) = min(s_min(i,j,k), EntPrev(1,l,j,k,eID))
        end do
        ! check stencil in eta
!#          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          s_min(i,j,k) = min(s_min(i,j,k), EntPrev(1,i,l,k,eID))
        end do
        ! check stencil in zeta
!#          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          s_min(i,j,k) = min(s_min(i,j,k), EntPrev(1,i,j,l,eID))
        end do
#endif /*barStates*/
        
#if LOCAL_ALPHA
        ! Subcell-wise limiter
        ! ATTENTION: 1) We need to find the minimum alpha for each interface
        !            2) Even if the current state is valid, some limiting might be needed for one/some of the interfaces
        !***************************************************************************************************************
        
        ! Initialization
        notInIter = .FALSE.
        dalpha1 = alpha_loc(i,j,k,eID) + dalpha_loc(i,j,k)  ! Here dalpha1 is a place-holder for the current blending coefficient
        
        ! Find blending coefficient for each neighbor node:
        ! -------------------------------------------------
        
        ! xi-
        call NewtonLoopNeighbor_SpecEntropy( Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                                             sJ(i,j,k,eID) * sWGP(i) * (ftilde_DG(:,i-1,j  ,k  ,eID) - ftilde_FV(:,i-1,j  ,k  ,eID)), & ! Anti-difussive flux in xi-
                                             s_min(i,j,k),dt,sdt,eps, & ! some constants
                                             dalpha1,notInIter)  ! in/out quantities
        ! xi+
        call NewtonLoopNeighbor_SpecEntropy( Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                                            -sJ(i,j,k,eID) * sWGP(i) * (ftilde_DG(:,i  ,j  ,k  ,eID) - ftilde_FV(:,i  ,j  ,k  ,eID)), & ! Anti-difussive flux in xi+
                                             s_min(i,j,k),dt,sdt,eps, & ! some constants
                                             dalpha1,notInIter)  ! in/out quantities
        ! eta-
        call NewtonLoopNeighbor_SpecEntropy( Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                                             sJ(i,j,k,eID) * sWGP(j) * (gtilde_DG(:,i  ,j-1,k  ,eID) - gtilde_FV(:,i  ,j-1,k  ,eID)), & ! Anti-difussive flux in eta-
                                             s_min(i,j,k),dt,sdt,eps, & ! some constants
                                             dalpha1,notInIter)  ! in/out quantities
        ! eta+
        call NewtonLoopNeighbor_SpecEntropy( Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                                            -sJ(i,j,k,eID) * sWGP(j) * (gtilde_DG(:,i  ,j  ,k  ,eID) - gtilde_FV(:,i  ,j  ,k  ,eID)), & ! Anti-difussive flux in eta+
                                             s_min(i,j,k),dt,sdt,eps, & ! some constants
                                             dalpha1,notInIter)  ! in/out quantities
        ! zeta-
        call NewtonLoopNeighbor_SpecEntropy( Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                                             sJ(i,j,k,eID) * sWGP(k) * (htilde_DG(:,i  ,j  ,k-1,eID) - htilde_FV(:,i  ,j  ,k-1,eID)), & ! Anti-difussive flux in zeta-
                                             s_min(i,j,k),dt,sdt,eps, & ! some constants
                                             dalpha1,notInIter)  ! in/out quantities
        ! zeta+
        call NewtonLoopNeighbor_SpecEntropy( Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                                            -sJ(i,j,k,eID) * sWGP(k) * (htilde_DG(:,i  ,j  ,k  ,eID) - htilde_FV(:,i  ,j  ,k  ,eID)), & ! Anti-difussive flux in zeta+
                                             s_min(i,j,k),dt,sdt,eps, & ! some constants
                                             dalpha1,notInIter)  ! in/out quantities
        
        ! Update dalpha_loc and dalpha
        dalpha_loc(i,j,k) = dalpha1 - alpha_loc(i,j,k,eID)
        dalpha = max(dalpha,dalpha1) 
        
#else
        ! Simple element-wise limiter
        !****************************
        
        ! Difference between goal entropy and current entropy (works as cycling criterion ONLY for element-wise limiting)
        as = (s_min(i,j,k) - Get_SpecEntropy(U(:,i,j,k)))
        if (as <= max(eps,abs(s_min(i,j,k))*NEWTON_ABSTOL)) cycle ! this DOF does NOT need entropy correction   

        ! Newton initialization:
        U_curr = U(:,i,j,k)
        dalpha1 = 0.0
        
        ! Perform Newton iterations
        NewtonLoopSpec: do iter=1, IDPMaxIter
          dalpha_old = dalpha1
          
          ! Evaluate dS/d(alpha)
          dSdalpha = dot_product(ConsToSpecEntropy(U_curr),dt*FFV_m_FDG(:,i,j,k,eID))
          
          if ( abs(dSdalpha)<eps) exit NewtonLoopSpec ! Nothing to do here!
          
          ! Update correction
          dalpha1 = dalpha1 + as / dSdalpha
          if (alpha(eID) + dalpha1 > alpha_maxIDP+eps) then
            dalpha1 = alpha_maxIDP - alpha(eID)
          elseif (dalpha1<0.0) then
            dalpha1=0.0
          else
            ! Check relative tolerance
            if ( abs(dalpha_old-dalpha1)<= NEWTON_RELTOL ) exit NewtonLoopSpec
          end if
          
          ! Get new U
          U_curr = U (:,i,j,k) + dalpha1 * dt * FFV_m_FDG(:,i,j,k,eID)
          
          ! Evaluate if goal entropy was achieved (and exit the Newton loop if that's the case)
          as = s_min(i,j,k)-Get_SpecEntropy(U_curr)
          
          ! Check absolute tolerance
          if ( abs(as) < max(eps,abs(s_min(i,j,k))*NEWTON_ABSTOL) ) exit NewtonLoopSpec  
          
        end do NewtonLoopSpec ! iter
        
        as = alpha(eID) + dalpha1 ! using as for new alpha
        if ( (as > alpha_maxIDP) .and. (as <= alpha_maxIDP+eps)) then
          dalpha1 = alpha_maxIDP - alpha(eID)
        end if
        
        if (iter > IDPMaxIter) notInIter =.TRUE.
        
        dalpha = max(dalpha,dalpha1) ! Compute the element-wise maximum correction
        
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; enddo !i,j,k
      
      if (notInIter) then
        write(*,'(A,I0,A,I0)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem
      end if
      
!       Do the correction if needed
!       ---------------------------
#if !(LOCAL_ALPHA)
      if ( dalpha > 0. ) then
        call PerformCorrection(U,Ut,dalpha    ,alpha(eID)          , &
                                    dt,sdt,eID)
      end if
#endif /*!(LOCAL_ALPHA)*/
  end subroutine IDP_LimitSpecEntropy  
#if LOCAL_ALPHA
!===================================================================================================================================
!> Newton loop for each neighbor node's contribution
!> ATTENTION: 1) We solve for beta:=1-alpha (to match FCT literature) and then return fluxo's alpha
!===================================================================================================================================
  subroutine NewtonLoopNeighbor_SpecEntropy(Usafe,Fan,smin,dt,sdt,eps,alpha,notInIter)
    use MOD_Equation_Vars , only: Get_SpecEntropy, ConsToSpecEntropy
    use MOD_IDP_Vars      , only: NEWTON_ABSTOL, NEWTON_RELTOL, IDPMaxIter, IDPgamma
    implicit none
    !-arguments----------------------------------------
    real, intent(in)    :: Usafe(PP_nVar) ! (safe) low-order solution
    real, intent(in)    :: Fan  (PP_nVar) ! Mass-matrix-scaled antidiffusive flux with the right sign. Fan := ± sWGP*sJ*(F_DG - F_FV) ... (± must be consistent with the equation update)
    real, intent(in)    :: smin           ! Goal entropy
    real, intent(in)    :: dt             ! Time-step size
    real, intent(in)    :: sdt            ! inverse of time-step size
    real, intent(in)    :: eps            ! Small value for the entropy deviation
    real, intent(inout) :: alpha          ! Current (nodal) blending coefficient. Gets updated here.
    logical, intent(inout) :: notInIter   ! 
    !-local-variables----------------------------------
    real :: beta           ! FCT blending coefficient (beta:=1-alpha). u^{n+1} = uFV^{n+1} + (beta*dt) \sum_i Fan_i
    real :: beta_old       ! beta from last iteration
    real :: Ucurr(PP_nVar)
    real :: dSdbeta        ! ds/d(beta)
    real :: as             ! Goal function as := smin - s (we find alpha such that as=0)
    real :: new_alpha
    integer :: iter
    !--------------------------------------------------
    
    ! Compute current FCT blending coef
    beta = (1.0 - alpha)
    
    ! Compute current directional update
    Ucurr = Usafe + beta * dt * IDPgamma * Fan
    
    ! Check entropy of this update
    as = smin - Get_SpecEntropy(Ucurr)
    if (as <= max(eps,abs(smin)*NEWTON_ABSTOL)) return ! this Neighbor contribution does NOT need entropy correction
    
    ! Perform Newton iterations
    NewtonLoopSpec: do iter=1, IDPMaxIter
      beta_old = beta
      
      ! Evaluate dS/d(alpha)
      dSdbeta = -dot_product(ConsToSpecEntropy(Ucurr),dt*IDPgamma*Fan)
      
      if ( abs(dSdbeta) == 0.0) then
        beta = 0.0
        exit NewtonLoopSpec ! Nothing to do here!
      end if
      
      ! Update correction
      beta = beta - as / dSdbeta
      
      ! Check bounds
      if ( (1.0 - beta) < alpha) then
        beta = (1.0 - alpha)
      elseif (beta < -eps) then
        beta = 0.0
      else
        ! Check relative tolerance
        if ( abs(beta_old-beta)<= NEWTON_RELTOL ) exit NewtonLoopSpec
      end if
      
      ! Get new U
      Ucurr = Usafe + beta *dt * IDPgamma * Fan
      
      ! Evaluate if goal entropy was achieved (and exit the Newton loop if that's the case)
      as = smin-Get_SpecEntropy(Ucurr)
      
      ! Check absolute tolerance
      if ( abs(as) < max(eps,abs(smin)*NEWTON_ABSTOL) ) exit NewtonLoopSpec  
          
    end do NewtonLoopSpec ! iter
    
    if ( (beta<0.0) .and. (beta>=-eps) ) then
      print*, '0>beta=',beta
      beta = 0.0
    end if
    
    new_alpha = 1.0 - beta
    if (alpha > new_alpha+eps) then
      print*, 'WTF... alpha is getting smaller (old/new)', alpha, new_alpha
      stop
    else
      alpha = new_alpha
    end if
    
    if (iter > IDPMaxIter) notInIter = .TRUE.
    
  end subroutine NewtonLoopNeighbor_SpecEntropy
#endif /*LOCAL_ALPHA*/
!===================================================================================================================================
!> Mathematical entropy correction (discrete local maximum principle)
!===================================================================================================================================
  subroutine IDP_LimitMathEntropy(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: dalpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_Mesh_Vars     , only: nElems, offsetElem
    use MOD_Equation_Vars , only: Get_MathEntropy, ConsToEntropy
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#else
    use MOD_IDP_Vars      , only: EntPrev
#endif /*barStates*/
    use MOD_IDP_Vars      , only: NEWTON_ABSTOL, NEWTON_RELTOL
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: alpha_maxIDP, IDPMaxIter, s_max
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    real    :: dalpha, dalpha1, dalpha_old
    real    :: as, U_curr(PP_nVar), dSdalpha
    integer :: i,j,k, l
    integer :: iter
    logical :: notInIter
    real, parameter :: eps = 1.e-14           ! Very small value
    !--------------------------------------------------------
      
      dalpha = -epsilon(1.0) ! Safe initialization
      
!       Compute correction factors
!       --------------------------
      notInIter = .FALSE.
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Get the limit states
        !*********************
        s_max(i,j,k) = -huge(1.0)
        
#if barStates
        ! Previous entropy of the node (ubar_ii)
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Uprev    (:,i  ,j  ,k  ,eID)))
        
        ! TODO: Compute for all interfaces before the loop!
        !xi+
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Ubar_xi  (:,i  ,j  ,k  ,eID)))
        !xi-
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Ubar_xi  (:,i-1,j  ,k  ,eID)))
        !eta+
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Ubar_eta (:,i  ,j  ,k  ,eID)))
        !eta-
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Ubar_eta (:,i  ,j-1,k  ,eID)))
        !zeta+
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Ubar_zeta(:,i  ,j  ,k  ,eID)))
        !zeta-
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Ubar_zeta(:,i  ,j  ,k-1,eID)))
#else
        ! check stencil in xi
!          do l = i+idx_m1(i), i+idx_p1(i) !no neighbor
        do l = i-1, i+1
          s_max(i,j,k) = max(s_max(i,j,k), EntPrev(1,l,j,k,eID))
        end do
        ! check stencil in eta
!          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          s_max(i,j,k) = max(s_max(i,j,k), EntPrev(1,i,l,k,eID))
        end do
        ! check stencil in zeta
!          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          s_max(i,j,k) = max(s_max(i,j,k), EntPrev(1,i,j,l,eID))
        end do
#endif /*barStates*/
        
#if LOCAL_ALPHA
        ! Subcell-wise limiter
        ! ATTENTION: 1) We need to find the minimum alpha for each interface
        !            2) Even if the current state is valid, some limiting might be needed for one/some of the interfaces
        !***************************************************************************************************************
        stop 'subcell limiter not yet implemented for mathematical entropy!!'
        !TODO!!
        
#else
        ! Simple element-wise limiter
        !****************************
        
        ! Difference between goal entropy and current entropy (works as cycling criterion ONLY for element-wise limiting)
        as = (s_max(i,j,k) - Get_MathEntropy(U(:,i,j,k)))
        if (as >= -max(eps,abs(s_max(i,j,k))*NEWTON_ABSTOL)) cycle ! this DOF does NOT need pressure correction   
      
        ! Newton initialization:
        U_curr = U(:,i,j,k)
        dalpha1 = 0.0
        
        ! Perform Newton iterations
        NewtonLoop: do iter=1, IDPMaxIter
          dalpha_old = dalpha1
          
          ! Evaluate dS/d(alpha)
          dSdalpha = dot_product(ConsToEntropy(U_curr),dt*FFV_m_FDG(:,i,j,k,eID))
          
          if ( abs(dSdalpha)<eps ) exit NewtonLoop ! Nothing to do here!
          
          ! Update correction
          dalpha1 = dalpha1 + as / dSdalpha
          if (alpha(eID) + dalpha1 > alpha_maxIDP) then
            dalpha1 = alpha_maxIDP - alpha(eID)
          elseif (dalpha1<0) then
            dalpha1=0.0
          else
            ! Check relative tolerance
            if ( abs(dalpha_old-dalpha1)<= NEWTON_RELTOL ) exit NewtonLoop
          end if
          
          ! Get new U
          U_curr = U (:,i,j,k) + dalpha1 * dt * FFV_m_FDG(:,i,j,k,eID)
          
          ! Evaluate if goal entropy was achieved (and exit the Newton loop if that's the case)
          as = s_max(i,j,k)-Get_MathEntropy(U_curr)
          
          ! Check absolute tolerance
          if ( abs(as) < max(eps,abs(s_max(i,j,k))*NEWTON_ABSTOL) ) exit NewtonLoop  
          
        end do NewtonLoop ! iter
        
        if (iter > IDPMaxIter) notInIter =.TRUE.
        
#if LOCAL_ALPHA
        ! Not the right way to do local alpha
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k),dalpha1)
#endif /*LOCAL_ALPHA*/
        dalpha = max(dalpha,dalpha1) ! Compute the element-wise maximum correction
        
#endif /*LOCAL_ALPHA*/
      end do       ; end do       ; enddo !i,j,k
      
      if (notInIter) then
        write(*,'(A,I0,A,I0)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem
      end if
      
!       Do the correction if needed
!       ---------------------------
#if !(LOCAL_ALPHA)
      if ( dalpha > 0. ) then
        call PerformCorrection(U,Ut,dalpha    ,alpha(eID)          , &
                                    dt,sdt,eID)
      end if
#endif /*!(LOCAL_ALPHA)*/
    
  end subroutine IDP_LimitMathEntropy
!===================================================================================================================================
!> Density and pressure positivity limiter. As explained in:
!> * Rueda-Ramírez, A. M., & Gassner, G. J. (2021). A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations of the Euler Equations. arXiv preprint arXiv:2102.06017.
!>    ... But with the possibility to be used locally
!===================================================================================================================================
  subroutine IDP_LimitPositivity(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha, PositCorrFactor
    use MOD_Mesh_Vars     , only: nElems, offsetElem
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: dalpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: Usafe, p_safe, rho_min, p_min, IDPDensityTVD
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: alpha_maxIDP, IDPMaxIter, NEWTON_ABSTOL, NEWTON_RELTOL
    use MOD_Equation_Vars , only: Get_Pressure, Get_dpdU
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    real    :: dalpha, dalpha1, dalpha_old
    real    :: a   ! a  = PositCorrFactor * rho_safe - rho
    real    :: ap  ! ap = (PositCorrFactor * p_safe   - p) / (kappa-1)
    real    :: pres
    real    :: alphadiff
    real    :: dpdU(PP_nVar), U_curr(PP_nVar)
    real    :: dp_dalpha
    integer :: i,j,k, iter
    logical :: NotInIter
    real, parameter :: eps = 1.e-14           ! Very small value
    !--------------------------------------------------------
      
!     ----------------------------------
!     Check if it makes sense correcting
!     ----------------------------------
      alphadiff = alpha(eID) - alpha_maxIDP
      if ( abs(alphadiff) < eps ) return ! Not much to do for this element...
      
!     ---------------
!     Correct density
!     ---------------
      
      dalpha = -epsilon(1.0) ! Safe initialization
        
!     Compute correction factors
!     --------------------------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Density correction
        rho_min(i,j,k) = merge (max(rho_min(i,j,k), PositCorrFactor * Usafe(1,i,j,k,eID)), PositCorrFactor * Usafe(1,i,j,k,eID), IDPDensityTVD) ! This writes the more restrictive bound into rho_min
        a = (rho_min(i,j,k) - U(1,i,j,k))
        if (a > 0.) then ! This DOF needs a correction
          if (abs(FFV_m_FDG(1,i,j,k,eID)) < eps) cycle
          dalpha1 = a * sdt / FFV_m_FDG(1,i,j,k,eID)
#if LOCAL_ALPHA
          ! TODO: Use a one-sided Zalesak limiter here!!
          dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k),dalpha1)
#endif /*LOCAL_ALPHA*/
          dalpha = max(dalpha,dalpha1)
        end if
        
      end do       ; end do       ; end do ! i,j,k
        
      
!       Do the correction if needed
!       ---------------------------
#if !(LOCAL_ALPHA)
      if ( dalpha > 0. ) then
        call PerformCorrection(U,Ut,dalpha    ,alpha(eID)          , &
                                    dt,sdt,eID)
      end if
#endif /*!(LOCAL_ALPHA)*/
!     ---------------
!     Correct pressure
!     ---------------
      
      dalpha = -epsilon(1.0) ! Safe initialization
#if LOCAL_ALPHA
      !TODO: Use a neigbor by neighbor Newton method here!!
      stop 'not implemented local alpha for pressure corr'
#endif /*LOCAL_ALPHA*/
      
!     Compute correction factors
!     --------------------------
      notInIter = .FALSE.
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        ! Current pressure and goal
        call Get_Pressure(U(:,i,j,k),pres)
        p_min(i,j,k) = PositCorrFactor * p_safe(i,j,k,eID)
        ap = (p_min(i,j,k) - pres)
        
        if (ap <= 0.) cycle ! this DOF does NOT need pressure correction
        
        ! Newton initialization:
        U_curr = U(:,i,j,k)
        dalpha1 = 0.0
        
        ! Perform Newton iterations
        NewtonLoop: do iter=1, IDPMaxIter
          dalpha_old = dalpha1
          
          ! Evaluate dp/d(alpha)
          call Get_dpdU(U_curr,dpdU)
          dp_dalpha = dot_product(dpdU,dt*FFV_m_FDG(:,i,j,k,eID))
          if ( abs(dp_dalpha)<eps ) exit NewtonLoop ! Nothing to do here!
          
          ! Update correction
          dalpha1 = dalpha1 + ap / dp_dalpha
          if (alpha(eID) + dalpha1 > alpha_maxIDP+eps) then
            dalpha1 = alpha_maxIDP - alpha(eID)
          elseif (dalpha1<0.0) then
            dalpha1=0.0
          else
            ! Check relative tolerance
            if ( abs(dalpha_old-dalpha1)<= NEWTON_RELTOL ) exit NewtonLoop
          end if
          
          ! Get new U and pressure
          U_curr = U (:,i,j,k) + dalpha1 * dt * FFV_m_FDG(:,i,j,k,eID)
          call Get_Pressure(U_curr,pres)
          
          ! Evaluate if goal pressure was achieved (and exit the Newton loop if that's the case)
          ap = p_min(i,j,k)-pres
          if ( (ap <= epsilon(p_min(i,j,k))) .and. (ap > -NEWTON_ABSTOL*p_min(i,j,k)) ) exit NewtonLoop  ! Note that we use an asymmetric tolerance!
        end do NewtonLoop ! iter
        
        if (iter > IDPMaxIter) notInIter =.TRUE.
        
#if LOCAL_ALPHA
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k),dalpha1)
#endif /*LOCAL_ALPHA*/
        dalpha = max(dalpha,dalpha1) ! Compute the element-wise maximum correction
      
      end do       ; end do       ; enddo !i,j,k
      
      if (notInIter) then
        write(*,'(A,I0,A,I0)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem
      end if
      
!       Do the correction if needed
!       ---------------------------
#if !(LOCAL_ALPHA)
      if ( dalpha > 0. ) then
        call PerformCorrection(U,Ut,dalpha    ,alpha(eID)          , &
                                    dt,sdt,eID)
      end if
#endif /*!(LOCAL_ALPHA)*/
  end subroutine IDP_LimitPositivity
!===================================================================================================================================
!> Takes dalpha/dalpha_loc U and Ut, and outputs the corrected U and Ut, and alpha/alpha_loc for visualization
!===================================================================================================================================
  pure subroutine PerformCorrection(U,Ut,dalpha    ,alpha    , &
#if LOCAL_ALPHA
                                         dalpha_loc,alpha_loc, &
#endif /*LOCAL_ALPHA*/
                                                             dt,sdt,eID)
    use MOD_PreProc   , only: PP_N
    use MOD_IDP_Vars  , only: alpha_maxIDP
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars, only: ftilde_DG, gtilde_DG, htilde_DG, ftilde_FV, gtilde_FV, htilde_FV
    use MOD_NFVSE_Vars, only: sWGP
    use MOD_Mesh_Vars , only: sJ
#else
    use MOD_IDP_Vars  , only: FFV_m_FDG
#endif /*LOCAL_ALPHA*/
    use MOD_NFVSE_Vars  , only: alpha_old
    implicit none
    !-arguments--------------------------------------
    real, intent(inout) :: U (PP_nVar, 0:PP_N  , 0:PP_N  , 0:PP_N)
    real, intent(inout) :: Ut(PP_nVar, 0:PP_N  , 0:PP_N  , 0:PP_N)
    real, intent(inout) :: dalpha                             ! Scaled element-wise d_alpha
    real, intent(inout) :: alpha                            ! Current element-wise alpha
#if LOCAL_ALPHA
    real, intent(inout) :: dalpha_loc  (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1) ! Scaled local d_alpha
    real, intent(inout) :: alpha_loc   ( 0:PP_N  , 0:PP_N  , 0:PP_N) ! Current local  alpha
#endif /*LOCAL_ALPHA*/
    real, intent(in)    :: dt,sdt                          ! time step and inverse
    integer, intent(in) :: eID
    !-local-variables--------------------------------
    real :: alphacont_loc ! A local container
    real :: alphacont     ! An element-wise container
    real :: my_corr(PP_nVar)
    integer :: i,j,k
    !------------------------------------------------
    
    ! Change the alpha for output
    alphacont  = alpha
    alpha      = alpha + dalpha
    
    ! Change inconsistent alphas
    if (alpha > alpha_maxIDP) then
      alpha = alpha_maxIDP
      dalpha  = alpha_maxIDP - alphacont
    end if
          
#if LOCAL_ALPHA
    ! Change the alpha for output
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      alphacont_loc    = alpha_loc(i,j,k)
      alpha_loc(i,j,k) = alpha_loc(i,j,k) + dalpha_loc(i,j,k)
      if (alpha_loc(i,j,k) > alpha_maxIDP .or. isnan(alpha_loc(i,j,k))) then
        alpha_loc(i,j,k) = alpha_maxIDP
        dalpha_loc (i,j,k) = alpha_maxIDP - alpha_old(eID)
      end if
    end do       ; end do       ; enddo
    
    ! Correct!
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      ! xi correction
      ! -------------
      ! left
      my_corr=-max(dalpha_loc(i-1,j  ,k  ),dalpha_loc(i  ,j  ,k  )) * sWGP(i) * (ftilde_DG(:,i-1,j,k,eID) - ftilde_FV(:,i-1,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! right
      my_corr=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i+1,j  ,k  )) * sWGP(i) * (ftilde_DG(:,i  ,j,k,eID) - ftilde_FV(:,i  ,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! eta correction
      ! --------------
      ! left
      my_corr=-max(dalpha_loc(i  ,j-1,k  ),dalpha_loc(i  ,j  ,k  )) * sWGP(j) * (gtilde_DG(:,i,j-1,k,eID) - gtilde_FV(:,i,j-1,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! right
      my_corr=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i  ,j+1,k  )) * sWGP(j) * (gtilde_DG(:,i,  j,k,eID) - gtilde_FV(:,i  ,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! zeta correction
      ! ---------------
      ! left
      my_corr=-max(dalpha_loc(i  ,j  ,k-1),dalpha_loc(i  ,j  ,k  )) * sWGP(k) * (htilde_DG(:,i,j,k-1,eID) - htilde_FV(:,i,j,k-1,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! right
      my_corr=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i  ,j  ,k+1)) * sWGP(k) * (htilde_DG(:,i,  j,k,eID) - htilde_FV(:,i  ,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
    end do       ; end do       ; enddo
#else
    ! Element-wise correction!
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      ! Correct U
      U (:,i,j,k) = U (:,i,j,k) + dalpha * dt * FFV_m_FDG(:,i,j,k,eID)
      ! Correct Ut
      Ut(:,i,j,k) = Ut(:,i,j,k) + dalpha * FFV_m_FDG(:,i,j,k,eID)
    end do       ; end do       ; enddo
#endif /*LOCAL_ALPHA*/
  end subroutine PerformCorrection
!===================================================================================================================================
!> Function to get the bar states! (LLF)
!===================================================================================================================================
#if barStates
  pure function GetBarStates(UL,UR,nv,t1,t2) result(Ubar)
    use MOD_Riemann       , only: RotateState, RotateFluxBack
    use MOD_Equation_Vars , only: SoundSpeed2
    USE MOD_Flux          , only: EvalOneEulerFlux1D
    implicit none
    !-arguments----------------------------------------
    real, intent(in) :: UL  (PP_nVar)
    real, intent(in) :: UR  (PP_nVar)
    real, intent(in) :: nv  (3)
    real, intent(in) :: t1  (3)
    real, intent(in) :: t2  (3)
    real             :: Ubar(PP_nVar)
    !-local-variables----------------------------------
    real :: UL_r(PP_nVar) ! Rotated left state
    real :: UR_r(PP_nVar) ! Rotated right state
    real :: FL_r(PP_nVar) ! Rotated left flux
    real :: FR_r(PP_nVar) ! Rotated right flux
    real :: lambdamax
    !--------------------------------------------------
    
    UL_r = RotateState(UL,nv,t1,t2)
    UR_r = RotateState(UR,nv,t1,t2)
    
    lambdamax = MAX(ABS(UL_r(2)/UL_r(1)),ABS(UR_r(2)/UR_r(1))) + SQRT(MAX(SoundSpeed2(UL_r),SoundSpeed2(UR_r)) )
!#    lambdamax = MAX(SQRT(SoundSpeed2(UL_r)) + ABS(UL_r(2)/UL_r(1)), sqrt(SoundSpeed2(UR_r)) + ABS(UR_r(2)/UR_r(1)) )  ! Will's lambda_max
    
    call EvalOneEulerFlux1D(UL_r, FL_r)
    call EvalOneEulerFlux1D(UR_r, FR_r)
    
    Ubar = 0.5*( UL_r + UR_r ) - (0.5/(lambdamax)) * (FR_r-FL_r)
    
    call RotateFluxBack(Ubar,nv,t1,t2)
  end function GetBarStates
#endif /*barStates*/
!===================================================================================================================================
!> Finalizes the IDP module
!===================================================================================================================================
  subroutine Finalize_IDP
    use MOD_NFVSE_Vars
    use MOD_IDP_Vars
    implicit none
    
    SDEALLOCATE (alpha_old)
    SDEALLOCATE (FFV_m_FDG)
    
    SDEALLOCATE (Usafe)
    SDEALLOCATE (p_safe)
    SDEALLOCATE (UPrev)
    
    SDEALLOCATE (EntPrev)
    
    SDEALLOCATE( Usafe_ext )
    SDEALLOCATE( Uprev_ext )
    SDEALLOCATE( EntPrev_master)
    SDEALLOCATE( EntPrev_slave )
    SDEALLOCATE( EntPrev_ext   )
    SDEALLOCATE( Ubar_xi  )
    SDEALLOCATE( Ubar_eta  )
    SDEALLOCATE( Ubar_zeta  )
    
#if LOCAL_ALPHA
    SDEALLOCATE ( alpha_loc )
    SDEALLOCATE ( ftilde_FV )
    SDEALLOCATE ( gtilde_FV )
    SDEALLOCATE ( htilde_FV )
    SDEALLOCATE ( ftilde_DG )
    SDEALLOCATE ( gtilde_DG )
    SDEALLOCATE ( htilde_DG )
    SDEALLOCATE ( dalpha_loc )
#endif /*LOCAL_ALPHA*/
    
    SDEALLOCATE ( rho_min )
    SDEALLOCATE ( rho_max )
    SDEALLOCATE ( s_min )
    SDEALLOCATE ( s_max )
    SDEALLOCATE ( p_min )
    
  end subroutine Finalize_IDP
#endif /*NFVSE_CORR*/  
end module MOD_IDP

