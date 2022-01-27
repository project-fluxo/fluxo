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
    call prms%CreateLogicalOption("IDPPressureTVD",  " IDP(TVD) correction on pressure? (uses a Zalesak limiter with LOCAL_ALPHA=ON)", "F")
    call prms%CreateLogicalOption( "IDPPositivity",  " IDP correction for positivity of density and pressure?", "F")
    
!   Additional options
!   ------------------
    call prms%CreateRealOption(  "PositCorrFactor",  " The correction factor for IDPPositivity=T", "0.1")
    call prms%CreateIntOption(        "IDPMaxIter",  " Maximum number of iterations for positivity limiter", "10")
    
    call prms%CreateLogicalOption( "IDPForce2D"   ,  " Force a 2D solution x-y / xi-eta??", "F")
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
    use MOD_Mesh_Vars  , only: nElems, sJ
#if !(barStates)
    use MOD_Mesh_Vars  , only: firstSlaveSide, LastSlaveSide, nSides
#endif /*!(barStates)*/
    use MOD_ReadInTools, only: GETINT, GETREAL, GETLOGICAL
    implicit none
    !-local-variables----------------------------------------
    integer :: i,j,k,eID
    !--------------------------------------------------------
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' IDP Methods: '
    
!   Get parameters
!   --------------
    
    ! IDP limiters
    IDPPositivity  = GETLOGICAL('IDPPositivity' ,'F')
    IDPDensityTVD  = GETLOGICAL('IDPDensityTVD' ,'F')
    IDPPressureTVD = GETLOGICAL('IDPPressureTVD' ,'F')
    IDPMathEntropy = GETLOGICAL('IDPMathEntropy','F')
    IDPSpecEntropy = GETLOGICAL('IDPSpecEntropy','F')
    
    ! Consistency check
    if (IDPMathEntropy .and. IDPSpecEntropy) then
      stop 'Only one of the two can be selected: IDPMathEntropy/IDPSpecEntropy'
    end if
    
    ! Additional options
    IDPForce2D      = GETLOGICAL('IDPForce2D','F')
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
    
    if (IDPDensityTVD) then
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
    
    if (IDPPressureTVD) then
      IDPneedsUprev     = .TRUE.
#if barStates
#if LOCAL_ALPHA
      IDPneedsUsafe = .TRUE.
#endif /*LOCAL_ALPHA*/
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev_ext = .TRUE.
#else
      IDPneedsUsafe = .TRUE.
#endif /*barStates*/
    end if
    
    if (IDPMathEntropy) then
      IDPneedsUprev     = .TRUE.
      IDPneedsUsafe     = .TRUE.
#if barStates
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev_ext = .TRUE.
#endif /*barStates*/
    end if
    
    if (IDPSpecEntropy) then
      IDPneedsUprev     = .TRUE.
      IDPneedsUsafe     = .TRUE.
#if barStates
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev_ext = .TRUE.
#endif /*barStates*/
    end if
   
! TODO: Add c-preprocessor definition for strict time-step restriction 
    IDPneedsUbar = .TRUE.
    IDPneedsUprev     = .TRUE.
    IDPneedsUprev_ext = .TRUE.
    
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
    if (IDPDensityTVD .or. IDPPressureTVD) then
      allocate( Usafe_ext    (PP_nVar, 0:PP_N  , 0:PP_N          ,6,nElems) )
    end if
#endif /*barStates*/
    
    ! Variables for local alpha
#if LOCAL_ALPHA
    allocate ( alpha_loc(0:PP_N,0:PP_N,0:PP_N,1:nElems) )
    allocate ( ftilde_DG(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N,nElems) )
    allocate ( gtilde_DG(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N,nElems) )
    allocate ( htilde_DG(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N,nElems) )
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
    if (IDPPositivity .or. IDPPressureTVD) then
      allocate ( p_min       (0:PP_N,0:PP_N,0:PP_N) )
    end if
    if (IDPPressureTVD) then
      allocate ( p_max       (0:PP_N,0:PP_N,0:PP_N) )
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
    
    ! Now the analyze bounds
#if DEBUG || IDP_CHECKBOUNDS
    idp_bounds_num = 0
    
    if (IDPDensityTVD .or. IDPPositivity) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'rho_min'
    end if
    
    if (IDPDensityTVD) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'rho_max'
    end if
    
    if (IDPSpecEntropy) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'ent_min'
    end if
    
    if (IDPMathEntropy) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'ent_max'
    end if
    
    if (IDPPositivity .or. IDPPressureTVD) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'p_min'
    end if
    
    if (IDPPressureTVD) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'p_max'
    end if
    idp_bounds_delta = 0.0
#endif /*DEBUG || IDP_CHECKBOUNDS*/
    
!   Finally enforce 2D condition
  if (IDPForce2D) then
    do eID=1, nElems
      do k=1, PP_N ; do j=0, PP_N ; do i=0, PP_N
        sJ(i,j,k,eID) = sJ(i,j,0,eID)
      end do       ; end do       ; end do
    end do
  end if
  
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
  subroutine Apply_IDP(U,Ut,dt,tIn)
    use MOD_PreProc     , only: PP_N
    use MOD_Mesh_Vars   , only: nElems
    use MOD_NFVSE_Vars  , only: alpha, alpha_old, maximum_alpha, amount_alpha, amount_alpha_steps
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars  , only: alpha_loc
    use MOD_Analyze_Vars, only: wGPVol
    use MOD_Mesh_Vars   , only: nElems, sJ
    use MOD_IDP_Vars    , only: dalpha_loc
    use MOD_NFVSE_Vars  , only: ftilde_DG, gtilde_DG, htilde_DG
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars    , only: IDPSpecEntropy, IDPMathEntropy, IDPDensityTVD, IDPPressureTVD, IDPPositivity, dalpha
    use MOD_IDP_Vars    , only: IDPForce2D, FFV_m_FDG
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: tIn                                       !< Current RK time (at the beginning of RK stage)
    !-local-variables----------------------------------------
    real :: sdt
    real :: curr_amount_alpha
    integer :: eID, i, j, k
    !--------------------------------------------------------

!   Initialize
!   ----------
    alpha_old = alpha
    sdt = 1./dt
    
!   Get the variables for limiting in the right position
!   ----------------------------------------------------
    call Get_IDP_Variables(U,dt,tIn)
    
!   Perform limiting!
!   -----------------
    do eID=1, nElems
!     Initialize dalpha
!     -----------------
      dalpha = 0.0
#if LOCAL_ALPHA
      dalpha_loc = 0.0
#endif /*LOCAL_ALPHA*/
!     Enforce 2D condition
!     --------------------
      if (IDPForce2D) then
        do k=1, PP_N ; do j=0, PP_N ; do i=0, PP_N
          FFV_m_FDG(:,i,j,k,eID) = FFV_m_FDG(:,i,j,0,eID)
        end do       ; end do       ; end do
#if LOCAL_ALPHA
        do k=1, PP_N ; do j=0, PP_N ; do i=-1, PP_N
          ftilde_DG(:,i,j,k,eID) = ftilde_DG(:,i,j,0,eID)
        end do       ; end do       ; end do
        do k=1, PP_N ; do j=-1, PP_N ; do i=0, PP_N
          gtilde_DG(:, i,j,k,eID) = gtilde_DG(:, i,j,0,eID)
        end do       ; end do       ; end do
#endif /*LOCAL_ALPHA*/  
      end if
!     Call all user-defined limiters to obtain dalpha
!     -----------------------------------------------
      if (IDPDensityTVD)  call IDP_LimitDensityTVD (U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      if (IDPPressureTVD) call IDP_LimitPressureTVD(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      if (IDPSpecEntropy) call IDP_LimitSpecEntropy(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      if (IDPMathEntropy) call IDP_LimitMathEntropy(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      if (IDPPositivity)  call IDP_LimitPositivity (U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      
!     Enforce 2D condition
!     --------------------
      if (IDPForce2D) then
        
#if LOCAL_ALPHA
        do k=1, PP_N ; do j=0, PP_N ; do i=0, PP_N
          U(:,i,j,k,eID) = U(:,i,j,0,eID)
          Ut(:,i,j,k,eID) = Ut(:,i,j,0,eID)
          dalpha_loc(i,j,k) = maxval(dalpha_loc(i,j,:))
        end do       ; end do       ; end do
#endif /*LOCAL_ALPHA*/  
      end if
!     Perform limiting
!     ----------------
      if ( dalpha>0.0 .or. isnan(dalpha) &
#if LOCAL_ALPHA
           .or. any(isnan(dalpha_loc))   &
#endif /*LOCAL_ALPHA*/
                       ) then
        call PerformCorrection(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dalpha    ,alpha(eID)          , &
#if LOCAL_ALPHA
                                                              dalpha_loc,alpha_loc(:,:,:,eID), &
#endif /*LOCAL_ALPHA*/
                                                                         dt,sdt,eID)
      end if

!     Check that we are within bounds
!     -------------------------------
#if DEBUG || IDP_CHECKBOUNDS
      call CheckBounds(U(:,:,:,:,eID),eID)
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
#if DEBUG || IDP_CHECKBOUNDS
  subroutine CheckBounds(U,eID)
    use MOD_Preproc
    use MOD_Globals
    use MOD_IDP_Vars      , only: rho_min, rho_max, s_min, s_max, p_min, p_max, idp_bounds_delta
    use MOD_IDP_Vars      , only: IDPDensityTVD, IDPSpecEntropy, IDPMathEntropy, IDPPositivity, IDPForce2D, IDPPressureTVD
    use MOD_Equation_Vars , only: Get_SpecEntropy, Get_MathEntropy, Get_Pressure
    use MOD_NFVSE_Vars    , only: alpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
#endif /*LOCAL_ALPHA*/
    implicit none
    !-arguments------------------------------------------------------------
    real   , intent(in) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N)
    integer, intent(in) :: eID
    !-local-variables------------------------------------------------------
    integer :: i,j,k,counter
    real    :: p
    !----------------------------------------------------------------------
    
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      counter=0
      if (IDPDensityTVD .or. IDPPositivity) then
        counter=counter+1
        if (IDPForce2D) rho_min(i,j,k) = rho_min(i,j,0)
        idp_bounds_delta(counter) = max(idp_bounds_delta(counter), rho_min(i,j,k) - U(1,i,j,k))
      end if
        
      if (IDPDensityTVD) then
        counter=counter+1
        if (IDPForce2D) rho_max(i,j,k) = rho_max(i,j,0)  
        idp_bounds_delta(counter) = max(idp_bounds_delta(counter), U(1,i,j,k) - rho_max(i,j,k))
      end if
      
      if (IDPSpecEntropy) then
        counter=counter+1
        if (IDPForce2D) s_min(i,j,k) = s_min(i,j,0)
        idp_bounds_delta(counter) = max(idp_bounds_delta(counter),s_min(i,j,k) - Get_SpecEntropy(U(:,i,j,k)))
      end if
      
      if (IDPMathEntropy) then
        counter=counter+1
        if (IDPForce2D) s_max(i,j,k) = s_max(i,j,0)
        idp_bounds_delta(counter) = max(idp_bounds_delta(counter),Get_MathEntropy(U(:,i,j,k)) - s_max(i,j,k))
      end if
      
      if (IDPPositivity .or. IDPPressureTVD) then
        counter=counter+1
        if (IDPForce2D) p_min(i,j,k) = p_min(i,j,0)
        call Get_Pressure(U(:,i,j,k),p)
        idp_bounds_delta(counter) = max(idp_bounds_delta(counter),p_min(i,j,k) - p)
      end if

      if (IDPPressureTVD) then
        counter=counter+1
        if (IDPForce2D) p_max(i,j,k) = p_max(i,j,0)
        call Get_Pressure(U(:,i,j,k),p)
        idp_bounds_delta(counter) = max(idp_bounds_delta(counter),p - p_max(i,j,k))
      end if
    end do       ; end do       ; end do ! i,j,k
  
  end subroutine CheckBounds
#endif /*DEBUG || IDP_CHECKBOUNDS*/
!===================================================================================================================================
!> Get the IDP variables in the right position to perform limiting
!> ATTENTION: 1) U_master and U_slave need to have the previous solution!
!>            2) If the bar states are *deactivated* and IDPDensityTVD, U_master and U_slave contain the safe solution when this 
!>               subroutine returns
!===================================================================================================================================
  subroutine Get_IDP_Variables(U,dt,tIn)
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
    use MOD_IDP_Vars      , only: EntPrev, EntPrev_master, EntPrev_slave, EntPrev_ext, IDPDensityTVD, IDPPressureTVD, IDPMathEntropy, IDPSpecEntropy
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
    ! For time step
    use MOD_Globals
    use MOD_IDP_Vars, only: maxdt_IDP
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ
    implicit none
    !-arguments------------------------------------------------------------
    real, intent(in) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems)
    real, intent(in) :: dt
    real, intent(in) :: tIn
    !-local-variables------------------------------------------------------
    integer :: i,j,k
    integer :: eID, sideID
    real :: lambdamax_xi  (-1:PP_N, 0:PP_N, 0:PP_N)
    real :: lambdamax_eta ( 0:PP_N,-1:PP_N, 0:PP_N)
    real :: lambdamax_zeta( 0:PP_N, 0:PP_N,-1:PP_N)
    real :: inv_dt
    !----------------------------------------------------------------------
    
!   Get the previous solution in place if needed!
!   *********************************************
    if (IDPneedsUprev_ext) then
      ! Gather the external Uprev in the right location
      call Get_externalU(PP_nVar,Uprev_ext,Uprev(:,0:PP_N,0:PP_N,0:PP_N,:),U_master,U_slave,tIn)
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
      maxdt_IDP = huge(1.0)
      do eID=1, nElems
        associate (SCM => SubCellMetrics(eID))
        do i=-1, PP_N ! i goes through the interfaces
          do k=0, PP_N  ; do j=0, PP_N ! j and k go through DOFs
            !xi
            call GetBarStates(Uprev(:,i,j,k,eID),Uprev(:,i+1,j,k,eID),SCM % xi   % nv (:,j,k,i),SCM % xi   % t1 (:,j,k,i),SCM % xi   % t2 (:,j,k,i), Ubar_xi  (:,i,j,k,eID), lambdamax_xi  (i,j,k))
            !eta
            call GetBarStates(Uprev(:,j,i,k,eID),Uprev(:,j,i+1,k,eID),SCM % eta  % nv (:,j,k,i),SCM % eta  % t1 (:,j,k,i),SCM % eta  % t2 (:,j,k,i), Ubar_eta (:,j,i,k,eID), lambdamax_eta (j,i,k))
            !zeta
            call GetBarStates(Uprev(:,j,k,i,eID),Uprev(:,j,k,i+1,eID),SCM % zeta % nv (:,j,k,i),SCM % zeta % t1 (:,j,k,i),SCM % zeta % t2 (:,j,k,i), Ubar_zeta(:,j,k,i,eID), lambdamax_zeta(j,k,i))
          end do        ; end do  ! j,k
        end do
        
        ! Compute maximum time step
        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          inv_dt = ( sWGP(i) * (lambdamax_xi  (i-1,j  ,k  ) * SCM % xi   % norm(j,k,i-1) + lambdamax_xi  (i,j,k) * SCM % xi   % norm(j,k,i)) + &
                     sWGP(j) * (lambdamax_eta (i  ,j-1,k  ) * SCM % eta  % norm(i,k,j-1) + lambdamax_eta (i,j,k) * SCM % eta  % norm(i,k,j)) + &
                     sWGP(k) * (lambdamax_zeta(i  ,j  ,k-1) * SCM % zeta % norm(i,j,k-1) + lambdamax_zeta(i,j,k) * SCM % zeta % norm(i,j,k)) ) * sJ(i,j,k,eID)
          maxdt_IDP = min (maxdt_IDP, 1./inv_dt)
        end do       ; end do       ; end do
        end associate
      end do !eID
      
#if MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxdt_IDP,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#endif /*MPI*/
      if (dt > maxdt_IDP) then
        SWRITE(*,'(A,2ES21.12)') "MAYDAY, we have a problem with the time step. Consider reducing CFLScale (dt, maxdt_IDP): ", dt, maxdt_IDP
      end if
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
      call Get_externalU(1,EntPrev_ext,EntPrev(:,0:PP_N,0:PP_N,0:PP_N,:),EntPrev_master,EntPrev_slave,tIn)
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
    if (IDPDensityTVD .or. IDPPressureTVD) then
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
      call Get_externalU(PP_nVar,Usafe_ext,Usafe(:,0:PP_N,0:PP_N,0:PP_N,:),U_master,U_slave,tIn)
      ! FIll Usafe with info
      do eID=1, nElems
        Usafe(:,    -1,0:PP_N,0:PP_N,eID) = Usafe_ext(:,0:PP_N,0:PP_N,5,eID)
        Usafe(:,PP_N+1,0:PP_N,0:PP_N,eID) = Usafe_ext(:,0:PP_N,0:PP_N,3,eID)
        Usafe(:,0:PP_N,    -1,0:PP_N,eID) = Usafe_ext(:,0:PP_N,0:PP_N,2,eID)
        Usafe(:,0:PP_N,PP_N+1,0:PP_N,eID) = Usafe_ext(:,0:PP_N,0:PP_N,4,eID)
        Usafe(:,0:PP_N,0:PP_N,    -1,eID) = Usafe_ext(:,0:PP_N,0:PP_N,1,eID)
        Usafe(:,0:PP_N,0:PP_N,PP_N+1,eID) = Usafe_ext(:,0:PP_N,0:PP_N,6,eID)
      end do !eID
    end if !IDPDensityTVD .or. IDPPressureTVD
#endif /*!(barStates)*/
    
  end subroutine Get_IDP_Variables
!===================================================================================================================================
!> Density TVD correction
!===================================================================================================================================
  subroutine IDP_LimitDensityTVD(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_IDP_Vars      , only: dalpha, alpha_maxIDP
    use MOD_Mesh_Vars     , only: nElems
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_NFVSE_Vars    , only: ftilde_DG, gtilde_DG, htilde_DG
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ
    use MOD_IDP_Vars      , only: Usafe, dalpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: IDPForce2D
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#else
#if !(LOCAL_ALPHA)
    use MOD_IDP_Vars      , only: Usafe
#endif /*!(LOCAL_ALPHA)*/
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
    real    :: dalpha1
    real    :: rho_safe
    real    :: a   ! a  = PositCorrFactor * rho_safe - rho
    real    :: Qp, Qm, Pp, Pm
    integer :: i,j,k,l
    !--------------------------------------------------------
    
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
        if (.not. IDPForce2D) then
        !zeta
        do l=k-1, k !l=max(k-1,0), min(k,PP_N-1)
          rho_min(i,j,k) = min(rho_min(i,j,k), Ubar_zeta(1,i  ,j  ,l  ,eID))
          rho_max(i,j,k) = max(rho_max(i,j,k), Ubar_zeta(1,i  ,j  ,l  ,eID))
        end do
        end if
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
        
        ! Positive contributions
        Pp = 0.0
        Pp = Pp + max(0.0, sWGP(i) * ftilde_DG(1,i-1,j  ,k  ,eID))
        Pp = Pp + max(0.0,-sWGP(i) * ftilde_DG(1,i  ,j  ,k  ,eID))
        Pp = Pp + max(0.0, sWGP(j) * gtilde_DG(1,i  ,j-1,k  ,eID))
        Pp = Pp + max(0.0,-sWGP(j) * gtilde_DG(1,i  ,j  ,k  ,eID))
        if (.not. IDPForce2D) then
        Pp = Pp + max(0.0, sWGP(k) * htilde_DG(1,i  ,j  ,k-1,eID))
        Pp = Pp + max(0.0,-sWGP(k) * htilde_DG(1,i  ,j  ,k  ,eID))
        end if
        Pp = Pp*sJ(i,j,k,eID)
        
        ! Negative contributions
        Pm = 0.0
        Pm = Pm + min(0.0, sWGP(i) * ftilde_DG(1,i-1,j  ,k  ,eID))
        Pm = Pm + min(0.0,-sWGP(i) * ftilde_DG(1,i  ,j  ,k  ,eID))
        Pm = Pm + min(0.0, sWGP(j) * gtilde_DG(1,i  ,j-1,k  ,eID))
        Pm = Pm + min(0.0,-sWGP(j) * gtilde_DG(1,i  ,j  ,k  ,eID))
        if (.not. IDPForce2D) then
        Pm = Pm + min(0.0, sWGP(k) * htilde_DG(1,i  ,j  ,k-1,eID))
        Pm = Pm + min(0.0,-sWGP(k) * htilde_DG(1,i  ,j  ,k  ,eID))
        end if
        Pm = Pm*sJ(i,j,k,eID)
        
        if (Pp==0.0) then
          Qp = 1.0
        else
          Qp = Qp/Pp
        end if
        
        if (Pm==0.0) then
          Qm = 1.0
        else
          Qm = Qm/Pm
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
        
        ! Change inconsistent alphas
        if ( (alpha(eID)+dalpha1 > alpha_maxIDP) .or. isnan(dalpha1)) then
          dalpha  = alpha_maxIDP - alpha(eID)
        else
          dalpha = max(dalpha,dalpha1)
        end if
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; end do ! i,j,k
      
  end subroutine IDP_LimitDensityTVD
!===================================================================================================================================
!> Pressure TVD correction
!===================================================================================================================================
  subroutine IDP_LimitPressureTVD(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_IDP_Vars      , only: dalpha, alpha_maxIDP
    use MOD_Mesh_Vars     , only: nElems
    use MOD_Equation_Vars , only: Get_Pressure, KappaM1
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_NFVSE_Vars    , only: ftilde_DG, gtilde_DG, htilde_DG
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ
    use MOD_IDP_Vars      , only: Usafe, dalpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: IDPForce2D
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta
#else
#if !(LOCAL_ALPHA)
    use MOD_IDP_Vars      , only: Usafe
#endif /*!(LOCAL_ALPHA)*/
#endif /*barStates*/
    use MOD_IDP_Vars      , only: FFV_m_FDG, p_min, p_max, p_safe
    use MOD_IDP_Vars      , only: Uprev
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    real    :: dalpha1
    real    :: p_goal, p
    real    :: a   ! a  = PositCorrFactor * p_safe - p
    real    :: Qp, Qm, Pp, Pm
    real    :: fp_xi_plus
    real    :: fp_xi_minus
    real    :: fp_eta_plus
    real    :: fp_eta_minus
    real    :: fp_zeta_plus
    real    :: fp_zeta_minus
    real    :: vel(3), v2s2
    real    :: FFV_m_FDG_p
    integer :: i,j,k,l
    !--------------------------------------------------------
    
!       Compute correction factors
!       --------------------------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          
        ! Get the limit states
        !*********************
        p_min(i,j,k) =  huge(1.0)
        p_max(i,j,k) = -huge(1.0)
        
#if barStates
        ! Previous sol
        call Get_Pressure(Uprev  (:,i  ,j  ,k  ,eID),p)
        p_min(i,j,k) = min(p_min(i,j,k), p)
        p_max(i,j,k) = max(p_max(i,j,k), p)
        
        !xi
        do l=i-1, i !min(i-1,0), max(i,PP_N-1)
          call Get_Pressure(Ubar_xi  (:,l  ,j  ,k  ,eID),p)
          p_min(i,j,k) = min(p_min(i,j,k), p)
          p_max(i,j,k) = max(p_max(i,j,k), p)
        end do
        !eta
        do l=j-1, j !l=max(j-1,0), min(j,PP_N-1)
          call Get_Pressure(Ubar_eta (:,i  ,l  ,k  ,eID),p)
          p_min(i,j,k) = min(p_min(i,j,k), p)
          p_max(i,j,k) = max(p_max(i,j,k), p)
        end do
        if (.not. IDPForce2D) then
        !zeta
        do l=k-1, k !l=max(k-1,0), min(k,PP_N-1)
          call Get_Pressure(Ubar_zeta(:,i  ,j  ,l  ,eID),p)
          p_min(i,j,k) = min(p_min(i,j,k), p)
          p_max(i,j,k) = max(p_max(i,j,k), p)
        end do
        end if
#else
        ! check stencil in xi
!          do l = i+idx_m1(i), i+idx_p1(i) !no neighbor
        do l = i-1, i+1
          call Get_Pressure(Usafe(:,l,j,k,eID),p)
          p_min(i,j,k) = min(p_min(i,j,k), p)
          p_max(i,j,k) = max(p_max(i,j,k), p)
        end do
        ! check stencil in eta
!          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          call Get_Pressure(Usafe(:,i,l,k,eID),p)
          p_min(i,j,k) = min(p_min(i,j,k), p)
          p_max(i,j,k) = max(p_max(i,j,k), p)
        end do
        ! check stencil in zeta
!          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          call Get_Pressure(Usafe(:,i,j,l,eID),p)
          p_min(i,j,k) = min(p_min(i,j,k), p)
          p_max(i,j,k) = max(p_max(i,j,k), p)
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
        Qp = max(0.0,(p_max(i,j,k)-p_safe(i,j,k,eID)*sdt))
        Qm = min(0.0,(p_min(i,j,k)-p_safe(i,j,k,eID)*sdt))
        
        ! Compute pressure antidiffusive fluxes
        vel = Uprev(2:4,i,j,k,eID)/Uprev(1,i,j,k,eID)
        v2s2 = 0.5 * sum(vel**2)
        fp_xi_minus  = KappaM1 * (ftilde_DG(5,i-1,j  ,k  ,eID) + v2s2 * ftilde_DG(1,i-1,j  ,k  ,eID) - sum(vel*ftilde_DG(2:4,i-1,j  ,k  ,eID)) )
        fp_xi_plus   = KappaM1 * (ftilde_DG(5,i  ,j  ,k  ,eID) + v2s2 * ftilde_DG(1,i  ,j  ,k  ,eID) - sum(vel*ftilde_DG(2:4,i  ,j  ,k  ,eID)) )
        fp_eta_minus = KappaM1 * (gtilde_DG(5,i  ,j-1,k  ,eID) + v2s2 * gtilde_DG(1,i  ,j-1,k  ,eID) - sum(vel*gtilde_DG(2:4,i  ,j-1,k  ,eID)) )
        fp_eta_plus  = KappaM1 * (gtilde_DG(5,i  ,j  ,k  ,eID) + v2s2 * gtilde_DG(1,i  ,j  ,k  ,eID) - sum(vel*gtilde_DG(2:4,i  ,j  ,k  ,eID)) )
        fp_zeta_minus= KappaM1 * (htilde_DG(5,i  ,j  ,k-1,eID) + v2s2 * htilde_DG(1,i  ,j  ,k-1,eID) - sum(vel*htilde_DG(2:4,i  ,j  ,k-1,eID)) )
        fp_zeta_plus = KappaM1 * (htilde_DG(5,i  ,j  ,k  ,eID) + v2s2 * htilde_DG(1,i  ,j  ,k  ,eID) - sum(vel*htilde_DG(2:4,i  ,j  ,k  ,eID)) )
        ! Positive contributions
        Pp = 0.0
        Pp = Pp + max(0.0, sWGP(i) * fp_xi_minus  )
        Pp = Pp + max(0.0,-sWGP(i) * fp_xi_plus   )
        Pp = Pp + max(0.0, sWGP(j) * fp_eta_minus )
        Pp = Pp + max(0.0,-sWGP(j) * fp_eta_plus  )
        if (.not. IDPForce2D) then   
        Pp = Pp + max(0.0, sWGP(k) * fp_zeta_minus)
        Pp = Pp + max(0.0,-sWGP(k) * fp_zeta_plus )
        end if
        Pp = Pp*sJ(i,j,k,eID)
        
        ! Negative contributions
        Pm = 0.0
        Pm = Pm + min(0.0, sWGP(i) * fp_xi_minus  )
        Pm = Pm + min(0.0,-sWGP(i) * fp_xi_plus   )
        Pm = Pm + min(0.0, sWGP(j) * fp_eta_minus )
        Pm = Pm + min(0.0,-sWGP(j) * fp_eta_plus  )
        if (.not. IDPForce2D) then   
        Pm = Pm + min(0.0, sWGP(k) * fp_zeta_minus)
        Pm = Pm + min(0.0,-sWGP(k) * fp_zeta_plus )
        end if
        Pm = Pm*sJ(i,j,k,eID)
        
        if (Pp==0.0) then
          Qp = 1.0
        else
          Qp = Qp/Pp
        end if
        
        if (Pm==0.0) then
          Qm = 1.0
        else
          Qm = Qm/Pm
        end if
        
        ! Compute correction as: (needed_alpha) - current_alpha = (1.0 - min(1.0,Qp,Qm)) - alpha_loc(i,j,k,eID)
        dalpha1 = 1.0 - min(1.0,Qp,Qm) - alpha_loc(i,j,k,eID)
        
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k),dalpha1)
        dalpha = max(dalpha,dalpha1)
        
#else
        call Get_Pressure(U(:,i,j,k),p)
        ! Simple element-wise limiter
        !****************************
        if ( p < p_min(i,j,k)) then
          p_goal = p_min(i,j,k)
        elseif (p > p_max(i,j,k)) then
          p_goal = p_max(i,j,k)
        else
          cycle !nothing to do here!
        end if
        
        ! Get pressure correcting flux 
        vel = Uprev(2:4,i,j,k,eID)/Uprev(1,i,j,k,eID)
        v2s2 = 0.5 * sum(vel**2)
        FFV_m_FDG_p = KappaM1 * (FFV_m_FDG(5,i,j,k,eID) + v2s2 * FFV_m_FDG(1,i,j,k,eID) - sum(vel*FFV_m_FDG(2:4,i,j,k,eID)) )
        if ( abs(FFV_m_FDG_p) == 0.0) cycle !nothing to do here!
        
        ! Density correction
        a = (p_goal - p)
        dalpha1 = a*sdt / FFV_m_FDG_p
        
        ! Change inconsistent alphas
        if ( (alpha(eID)+dalpha1 > alpha_maxIDP) .or. isnan(dalpha1)) then
          dalpha  = alpha_maxIDP - alpha(eID)
        else
          dalpha = max(dalpha,dalpha1)
        end if
        
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; end do ! i,j,k
      
  end subroutine IDP_LimitPressureTVD
!===================================================================================================================================
!> Specific entropy correction (discrete local minimum principle)
!===================================================================================================================================
  subroutine IDP_LimitSpecEntropy(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_IDP_Vars      , only: dalpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: dalpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: Usafe, IDPForce2D
    use MOD_IDP_Vars      , only: IDPparam_t
    use MOD_Mesh_Vars     , only: nElems
    use MOD_Equation_Vars , only: Get_SpecEntropy, ConsToSpecEntropy
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#else
    use MOD_IDP_Vars      , only: EntPrev
#endif /*barStates*/
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: IDPMaxIter, s_min
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    integer :: i,j,k, l
    logical :: notInIter
    type(IDPparam_t) :: param ! Parameters for Newton's method
    real             :: new_alpha
    !--------------------------------------------------------
      
!     Compute correction factors
!     --------------------------
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
        if (.not. IDPForce2D) then
        !zeta+
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Ubar_zeta(:,i  ,j  ,k  ,eID)))
        !zeta-
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Ubar_zeta(:,i  ,j  ,k-1,eID)))
        end if
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
        
        ! Compute the needed blending coefficient with a Newton's method
        !***************************************************************
        
        ! Initialization
        param % dt    = dt
        param % bound = s_min(i,j,k)
        
#if LOCAL_ALPHA
        ! Get the current alpha
        new_alpha = alpha_loc(i,j,k,eID) + dalpha_loc(i,j,k)
        ! Perform Newton's method to find the new alpha (the antidiffusive fluxes are specified therein)
        call NewtonLoops_LocalAlpha(param,i,j,k,eID,new_alpha,notInIter,SpecEntropy_Goal,SpecEntropy_dGoal_dbeta,SpecEntropy_InitialCheck,Standard_FinalCheck)
        ! Update dalpha_loc and dalpha
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k), new_alpha - alpha_loc(i,j,k,eID))
        dalpha = max(dalpha,dalpha_loc(i,j,k)) 
#else
        ! Get the current alpha
        new_alpha = alpha(eID) + dalpha
        ! Specify the antidiffusive flux
        param % F_antidiff = -FFV_m_FDG(:,i,j,k,eID)
        ! Perform Newton's method to find the new alpha
        call NewtonLoop(Usafe(:,i,j,k,eID),param,new_alpha,notInIter,SpecEntropy_Goal,SpecEntropy_dGoal_dbeta,SpecEntropy_InitialCheck,Standard_FinalCheck)
        ! Update dalpha
        dalpha = max(dalpha,new_alpha-alpha(eID)) 
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; enddo !i,j,k
      
!~      if (notInIter) then
!~        write(*,'(A,I0,A,I0,A,ES21.12)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem, '. alpha = ', alpha(eID)+dalpha
!~      end if

  contains
!===================================================================================================================================
!>  Goal function for the specific entropy Newton's method
!===================================================================================================================================
    pure function SpecEntropy_Goal(param,Ucurr) result(goal)
      use MOD_Equation_Vars , only: Get_SpecEntropy
      use MOD_IDP_Vars      , only: IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: Ucurr(PP_nVar) ! Current solution
      real                         :: goal
      
      goal = param % bound - Get_SpecEntropy(Ucurr)
      
    end function SpecEntropy_Goal
!===================================================================================================================================
!>  Derivative of goal function with respect to (FCT) blending coefficient (beta:=1-alpha) for the specific entropy Newton's method
!===================================================================================================================================
    pure function SpecEntropy_dGoal_dbeta(param,Ucurr) result(dGoal_dbeta)
      use MOD_Equation_Vars , only: ConsToSpecEntropy
      use MOD_IDP_Vars      , only: IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: Ucurr(PP_nVar) ! Current solution
      real                         :: dGoal_dbeta
      
      dGoal_dbeta = -dot_product(ConsToSpecEntropy(Ucurr),param % dt * param % F_antidiff)
      
    end function SpecEntropy_dGoal_dbeta
!===================================================================================================================================
!>  InitialCheck for the specific entropy Newton's method: Is the current state admissible?
!===================================================================================================================================
    pure function SpecEntropy_InitialCheck(param,goalFunction) result(check)
      use MOD_IDP_Vars, only: NEWTON_ABSTOL, IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: goalFunction ! Current solution
      logical                      :: check
      
      check = (goalFunction <= max(NEWTON_ABSTOL,abs(param % bound)*NEWTON_ABSTOL))
      
    end function SpecEntropy_InitialCheck
  end subroutine IDP_LimitSpecEntropy  
!===================================================================================================================================
!> Mathematical entropy correction (discrete local maximum principle)
!===================================================================================================================================
  subroutine IDP_LimitMathEntropy(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_IDP_Vars      , only: dalpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: dalpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: IDPparam_t, Usafe
    use MOD_Mesh_Vars     , only: nElems
    use MOD_Equation_Vars , only: Get_MathEntropy, ConsToEntropy
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#else
    use MOD_IDP_Vars      , only: EntPrev
#endif /*barStates*/
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: IDPMaxIter, s_max
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    real    :: dalpha1
    integer :: i,j,k, l
    logical :: notInIter
    type(IDPparam_t) :: param ! Parameters for Newton's method
    real             :: new_alpha
    !--------------------------------------------------------
      
      
!     Compute correction factors
!     --------------------------
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
        
        ! Compute the needed blending coefficient with a Newton's method
        !***************************************************************
        
        ! Initialization
        param % dt    = dt
        param % bound = s_max(i,j,k)
        
#if LOCAL_ALPHA
        ! Get the current alpha
        new_alpha = alpha_loc(i,j,k,eID) + dalpha_loc(i,j,k)
        ! Perform Newton's method to find the new alpha (the antidiffusive fluxes are specified therein)
        call NewtonLoops_LocalAlpha(param,i,j,k,eID,new_alpha,notInIter,MathEntropy_Goal,MathEntropy_dGoal_dbeta,MathEntropy_InitialCheck,Standard_FinalCheck)
        ! Update dalpha_loc and dalpha
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k), new_alpha - alpha_loc(i,j,k,eID))
        dalpha = max(dalpha,dalpha_loc(i,j,k)) 
#else
        ! Get the current alpha
        new_alpha = alpha(eID) + dalpha
        ! Specify the antidiffusive flux
        param % F_antidiff = -FFV_m_FDG(:,i,j,k,eID)
        ! Perform Newton's method to find the new alpha
        call NewtonLoop(Usafe(:,i,j,k,eID),param,new_alpha,notInIter,MathEntropy_Goal,MathEntropy_dGoal_dbeta,MathEntropy_InitialCheck,Standard_FinalCheck)
        ! Update dalpha
        dalpha = max(dalpha,new_alpha-alpha(eID)) 
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; enddo !i,j,k
      
!~      if (notInIter) then
!~        write(*,'(A,I0,A,I0,A,ES21.12)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem, '. alpha = ', alpha(eID)+dalpha
!~      end if
      
  contains
!===================================================================================================================================
!>  Goal function for the mathematical entropy Newton's method
!===================================================================================================================================
    pure function MathEntropy_Goal(param,Ucurr) result(goal)
      use MOD_Equation_Vars , only: Get_MathEntropy
      use MOD_IDP_Vars      , only: IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: Ucurr(PP_nVar) ! Current solution
      real                         :: goal
      
      goal = param % bound - Get_MathEntropy(Ucurr)
      
    end function MathEntropy_Goal
!===================================================================================================================================
!>  Derivative of goal function with respect to (FCT) blending coefficient (beta:=1-alpha) for the mathematical entropy Newton's method
!===================================================================================================================================
    pure function MathEntropy_dGoal_dbeta(param,Ucurr) result(dGoal_dbeta)
      use MOD_Equation_Vars , only: ConsToEntropy
      use MOD_IDP_Vars      , only: IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: Ucurr(PP_nVar) ! Current solution
      real                         :: dGoal_dbeta
      
      dGoal_dbeta = -dot_product(ConsToEntropy(Ucurr),param % dt * param % F_antidiff)
      
    end function MathEntropy_dGoal_dbeta
!===================================================================================================================================
!>  InitialCheck for the mathematical entropy Newton's method: Is the current state admissible?
!===================================================================================================================================
    pure function MathEntropy_InitialCheck(param,goalFunction) result(check)
      use MOD_IDP_Vars, only: NEWTON_ABSTOL, IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: goalFunction ! Current solution
      logical                      :: check
      
      check = (goalFunction >= -max(NEWTON_ABSTOL,abs(param % bound)*NEWTON_ABSTOL))
      
    end function MathEntropy_InitialCheck

  end subroutine IDP_LimitMathEntropy
!===================================================================================================================================
!> Density and pressure positivity limiter. As explained in:
!> * Rueda-Ramírez, A. M., & Gassner, G. J. (2021). A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations of the Euler Equations. arXiv preprint arXiv:2102.06017.
!>    ... But with the possibility to be used locally
!===================================================================================================================================
  subroutine IDP_LimitPositivity(U,Ut,dt,sdt,eID)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha, PositCorrFactor
    use MOD_Mesh_Vars     , only: nElems
    use MOD_IDP_Vars      , only: dalpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: dalpha_loc
    use MOD_NFVSE_Vars    , only: ftilde_DG, gtilde_DG, htilde_DG
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: Usafe, p_safe, rho_min, p_min, IDPDensityTVD, IDPPressureTVD, IDPForce2D
    use MOD_IDP_Vars      , only: FFV_m_FDG, IDPparam_t, IDPMaxIter
    use MOD_Equation_Vars , only: Get_Pressure, Get_dpdU
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    real    :: dalpha1
    real    :: a      ! a  = PositCorrFactor * rho_safe - rho
    real    :: Qm, Pm ! Zalesak's limiter's variables
    integer :: i,j,k, iter
    logical :: NotInIter
    type(IDPparam_t) :: param ! Parameters for Newton's method
    real             :: new_alpha
    !--------------------------------------------------------
      
!     ---------------
!     Correct density
!     ---------------
        
!     Compute correction factors
!     --------------------------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Compute density bound
        !**********************
        rho_min(i,j,k) = merge (max(rho_min(i,j,k), PositCorrFactor * Usafe(1,i,j,k,eID)), PositCorrFactor * Usafe(1,i,j,k,eID), IDPDensityTVD) ! This writes the more restrictive bound into rho_min
#if LOCAL_ALPHA
        ! Real one-sided Zalesak-type limiter
        ! * Zalesak (1979). "Fully multidimensional flux-corrected transport algorithms for fluids"
        ! * Kuzmin et al. (2010). "Failsafe flux limiting and constrained data projections for equations of gas dynamics"
        ! ATTENTION: 1) The Zalesak limiter has to be computed, even if the state is valid, because the correction is 
        !               for each interface, not each node
        !****************************************************************************************************************
        
        ! Upper/lower bounds for admissible increments
        Qm = min(0.0,(rho_min(i,j,k)-Usafe(1,i,j,k,eID))*sdt)
        
        ! Negative contributions
        Pm = 0.0
        Pm = Pm + min(0.0,  sWGP(i) * ftilde_DG(1,i-1,j  ,k  ,eID))
        Pm = Pm + min(0.0, -sWGP(i) * ftilde_DG(1,i  ,j  ,k  ,eID))
        Pm = Pm + min(0.0,  sWGP(j) * gtilde_DG(1,i  ,j-1,k  ,eID))
        Pm = Pm + min(0.0, -sWGP(j) * gtilde_DG(1,i  ,j  ,k  ,eID))
        if (.not. IDPForce2D) then
        Pm = Pm + min(0.0,  sWGP(k) * htilde_DG(1,i  ,j  ,k-1,eID))
        Pm = Pm + min(0.0, -sWGP(k) * htilde_DG(1,i  ,j  ,k  ,eID))
        end if
        Pm = Pm*sJ(i,j,k,eID)
        
        if (Pm<0.0) then
          Qm = min(1.0,Qm/Pm)
        else
          Qm = 1.0
        end if
        
        ! Compute correction as: (needed_alpha) - current_alpha = (1.0 - Qm) - alpha_loc(i,j,k,eID)
        dalpha1 = (1.0 - Qm) - alpha_loc(i,j,k,eID)
        
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k),dalpha1)
        dalpha = max(dalpha,dalpha1)
        
#else
        a = (rho_min(i,j,k) - U(1,i,j,k)) * sdt
        if (a > 0.) then ! This DOF needs a correction
          if ( abs(FFV_m_FDG(1,i,j,k,eID)) == 0.0) cycle !nothing to do here!
          dalpha1 = a / FFV_m_FDG(1,i,j,k,eID)
          dalpha = max(dalpha,dalpha1)
        end if

#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; end do ! i,j,k
      
!     ---------------
!     Correct pressure
!     ---------------
      
!     Compute correction factors
!     --------------------------
      notInIter = .FALSE.
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Compute bound
        ! *************
        p_min(i,j,k) = PositCorrFactor * p_safe(i,j,k,eID)
        p_min(i,j,k) = merge (max(p_min(i,j,k), PositCorrFactor * p_safe(i,j,k,eID)), PositCorrFactor * p_safe(i,j,k,eID), IDPPressureTVD) ! This writes the more restrictive bound into p_min
        
        ! Compute the needed blending coefficient with a Newton's method
        ! TODO: Check if an asymmetric tolerance is really needed
        !***************************************************************
        
        ! Initialization
        param % dt    = dt
        param % bound = p_min(i,j,k)
        
#if LOCAL_ALPHA
        ! Get the current alpha
        new_alpha = alpha_loc(i,j,k,eID) + dalpha_loc(i,j,k)
        ! Perform Newton's method to find the new alpha (the antidiffusive fluxes are specified therein)
        call NewtonLoops_LocalAlpha(param,i,j,k,eID,new_alpha,notInIter,Pressure_Goal,Pressure_dGoal_dbeta,Pressure_InitialCheck,Pressure_FinalCheck)
        ! Update dalpha_loc and dalpha
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k), new_alpha - alpha_loc(i,j,k,eID))
        dalpha = max(dalpha,dalpha_loc(i,j,k)) 
#else
        ! Get the current alpha
        new_alpha = alpha(eID) + dalpha
        ! Specify the antidiffusive flux
        param % F_antidiff = -FFV_m_FDG(:,i,j,k,eID)
        ! Perform Newton's method to find the new alpha
        call NewtonLoop(Usafe(:,i,j,k,eID),param,new_alpha,notInIter,Pressure_Goal,Pressure_dGoal_dbeta,Pressure_InitialCheck,Pressure_FinalCheck)
        ! Update dalpha
        dalpha = max(dalpha,new_alpha-alpha(eID)) 
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; enddo !i,j,k
      
!~      if (notInIter) then
!~        write(*,'(A,I0,A,I0,A,ES21.12)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem, '. alpha = ', alpha(eID)+dalpha
!~      end if
      
  contains
!===================================================================================================================================
!>  Goal function for the pressure Newton's method
!===================================================================================================================================
    pure function Pressure_Goal(param,Ucurr) result(goal)
      use MOD_Equation_Vars , only: Get_Pressure
      use MOD_IDP_Vars      , only: IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: Ucurr(PP_nVar) ! Current solution
      real                         :: goal
      !-local-variables--------------------------
      real :: p
      !------------------------------------------
      
      call Get_Pressure(Ucurr,p)
      goal = param % bound - p
      
    end function Pressure_Goal
!===================================================================================================================================
!>  Derivative of goal function with respect to (FCT) blending coefficient (beta:=1-alpha) for the pressure Newton's method
!===================================================================================================================================
    pure function Pressure_dGoal_dbeta(param,Ucurr) result(dGoal_dbeta)
      use MOD_Equation_Vars , only: Get_dpdU
      use MOD_IDP_Vars      , only: IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: Ucurr(PP_nVar) ! Current solution
      real                         :: dGoal_dbeta
      !-local-variables--------------------------
      real :: dpdu(PP_nVar)
      !------------------------------------------
      
      call Get_dpdU(Ucurr,dpdu)
      dGoal_dbeta = -dot_product(dpdu,param % dt * param % F_antidiff)
      
    end function Pressure_dGoal_dbeta
!===================================================================================================================================
!>  InitialCheck for the pressure Newton's method: Is the current state admissible?
!===================================================================================================================================
    pure function Pressure_InitialCheck(param,goalFunction) result(check)
      use MOD_IDP_Vars, only: IDPparam_t
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: goalFunction ! Current solution
      logical                      :: check
      
      check = (goalFunction <= 0.0)
      
    end function Pressure_InitialCheck
!===================================================================================================================================
!>  FinalCheck for the pressure Newton's method: Is the current state admissible? (asymmetric condition)
!===================================================================================================================================
    pure function Pressure_FinalCheck(param,goalFunction) result(check)
      use MOD_IDP_Vars, only: IDPparam_t, NEWTON_ABSTOL
      implicit none
      type(IDPparam_t), intent(in) :: param
      real            , intent(in) :: goalFunction ! Current solution
      logical                      :: check
      
      check = ( goalFunction <= epsilon(1.0) ) .and. (goalFunction > -max(NEWTON_ABSTOL,abs(param % bound)*NEWTON_ABSTOL))
      
    end function Pressure_FinalCheck
  end subroutine IDP_LimitPositivity
#if LOCAL_ALPHA
!===================================================================================================================================
!> Run all Newton loops for a particular subcell
!> ATTENTION: 1) We need to find the minimum alpha for each interface using the correspoinding antidiffusive flux
!>            2) Even if the current state is valid, some limiting might be needed for one/some of the interfaces
!===================================================================================================================================
  subroutine NewtonLoops_LocalAlpha(param,i,j,k,eID,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    use MOD_NFVSE_Vars    , only: ftilde_DG, gtilde_DG, htilde_DG
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ
    use MOD_IDP_Vars      , only: IDPForce2D, Usafe, IDPparam_t, i_sub_Goal, i_sub_InitialCheck, IDPgamma
    implicit none
    !-arguments--------------------------------------
    type(IDPparam_t), intent(inout) :: param
    integer         , intent(in)    :: i,j,k      ! local coordinate index
    integer         , intent(in)    :: eID        ! Element index
    real            , intent(inout) :: alpha      ! 
    logical         , intent(inout) :: notInIter  ! 
    procedure(i_sub_Goal)           :: Goal
    procedure(i_sub_Goal)           :: dGoal_dbeta
    procedure(i_sub_InitialCheck)   :: InitialCheck
    procedure(i_sub_InitialCheck)   :: FinalCheck
    !------------------------------------------------
    
    ! xi-
    param % F_antidiff =  IDPgamma * sJ(i,j,k,eID) * sWGP(i) * (ftilde_DG(:,i-1,j  ,k  ,eID)) ! Anti-difussive flux in xi-
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    
    ! xi+
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(i) * (ftilde_DG(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in xi+
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    ! eta-
    param % F_antidiff =  IDPgamma * sJ(i,j,k,eID) * sWGP(j) * (gtilde_DG(:,i  ,j-1,k  ,eID)) ! Anti-difussive flux in eta-
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    
    ! eta+
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(j) * (gtilde_DG(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in eta+
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    
    if (.not. IDPForce2D) then
    ! zeta-
    param % F_antidiff =  IDPgamma * sJ(i,j,k,eID) * sWGP(k) * (htilde_DG(:,i  ,j  ,k-1,eID)) ! Anti-difussive flux in zeta-
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    
    ! zeta+
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(k) * (htilde_DG(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in zeta+
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    end if
  
  end subroutine NewtonLoops_LocalAlpha
#endif /*LOCAL_ALPHA*/
!===================================================================================================================================
!> General Newton loop (can be used for the element-wise or subcell-wise limiters)
!> ATTENTION: 1) We solve for beta:=1-alpha (to match FCT literature) and then return fluxo's alpha
!>            2) This is a Newton-bisection algorithm. We specify lower and upper bounds for beta and solve with Newton's method,
!>               but we switch to a bisection step if our Newton's method goes out of bounds or d(goal)/d(beta) = 0.
!===================================================================================================================================
  subroutine NewtonLoop(Usafe,param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    use MOD_IDP_Vars      , only: NEWTON_ABSTOL, NEWTON_RELTOL, IDPMaxIter
    use MOD_IDP_Vars      , only: i_sub_Goal, i_sub_InitialCheck, IDPparam_t
    implicit none
    !-arguments----------------------------------------
    real            , intent(in)    :: Usafe(PP_nVar) ! (safe) low-order solution
    type(IDPparam_t), intent(in)    :: param          ! IDP parameters 
    real            , intent(inout) :: alpha          ! Current (nodal) blending coefficient. Gets updated here.
    logical         , intent(inout) :: notInIter      ! Didn't Newton's method converge?
    procedure(i_sub_Goal)           :: Goal           ! Goal function (that we want to bring to <= or >= 0)
    procedure(i_sub_Goal)           :: dGoal_dbeta    ! Derivative of goal function with respect to beta
    procedure(i_sub_InitialCheck)   :: InitialCheck   ! Function to check if goal is fulfilled
    procedure(i_sub_InitialCheck)   :: FinalCheck     ! Function to exit the Newton loop (strict check of goal function)
    !-local-variables----------------------------------
    real :: beta           ! FCT blending coefficient (beta:=1-alpha). u^{n+1} = uFV^{n+1} + (beta*dt) \sum_i Fan_i
    real :: beta_old       ! beta from last iteration
    real :: Ucurr(PP_nVar)
    real :: dSdbeta        ! d(goal)/d(beta)
    real :: as             ! Goal function. For example: as = smin - s (we find alpha such that as=0)
    real :: new_alpha
    real :: beta_L, beta_R ! Left and right bounds for beta
    integer :: iter
    !--------------------------------------------------
    
    ! Compute current FCT blending coef
    beta = (1.0 - alpha)
    
    ! Get initial lower and upper bounds
    beta_L = 0.0   ! (alpha=1)
    beta_R = beta  ! We do not allow a higher beta (lower alpha) than the current one
    
    ! Compute current directional update
    Ucurr = Usafe + beta * param % dt * param % F_antidiff
    
    ! Perform initial check
    as = Goal(param,Ucurr)
    if (InitialCheck(param,as)) return ! this node does NOT need a correction
    
    ! Perform Newton iterations
    TheNewtonLoop: do iter=1, IDPMaxIter
      beta_old = beta
      
      ! Evaluate d(goal)/d(beta)
      dSdbeta = dGoal_dbeta(param,Ucurr)
      
      if (dSdbeta /= 0.0) then
        ! Update beta with Newton's method
        beta = beta - as / dSdbeta
      end if
      
      ! Check bounds
      if ( (beta < beta_L) .or. (beta > beta_R) .or. (dSdbeta == 0.0) .or. isnan(beta) ) then
        ! Out of bounds, do a bisection step
        beta = 0.5 * (beta_L + beta_R)
        ! Get new U
        Ucurr = Usafe + beta * param % dt * param % F_antidiff
        ! Check if new beta fulfills the condition and update bounds
        as = Goal(param,Ucurr)
        if (InitialCheck(param,as)) then
          ! New beta fulfills condition
          beta_L = beta
        else
          ! New beta does not fulfill condition
          beta_R = beta
        end if
      else
        ! Get new U
        Ucurr = Usafe + beta * param % dt * param % F_antidiff
        ! Evaluate goal function
        as = Goal(param,Ucurr)
      end if

      ! Check relative tolerance
      if ( abs(beta_old-beta)<= NEWTON_RELTOL ) exit TheNewtonLoop
      
      ! Check absolute tolerance
      if ( FinalCheck(param,as) ) exit TheNewtonLoop  
          
    end do TheNewtonLoop ! iter
    
    new_alpha = 1.0 - beta
    if (alpha > new_alpha+NEWTON_ABSTOL) then
      print*, 'WTF... alpha is getting smaller (old/new)', alpha, new_alpha
      stop
    else
      alpha = new_alpha
    end if
    
    if (iter > IDPMaxIter) notInIter = .TRUE.
    
  end subroutine NewtonLoop
!===================================================================================================================================
!>  standard FinalCheck for the Newton's method: Is the current state admissible?
!===================================================================================================================================
  pure function Standard_FinalCheck(param,goalFunction) result(check)
    use MOD_IDP_Vars, only: IDPparam_t, NEWTON_ABSTOL
    implicit none
    type(IDPparam_t), intent(in) :: param
    real            , intent(in) :: goalFunction ! Current solution
    logical                      :: check
    
    check = (abs(goalFunction) < max(NEWTON_ABSTOL,abs(param % bound)*NEWTON_ABSTOL))
    
  end function Standard_FinalCheck
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
    use MOD_NFVSE_Vars, only: ftilde_DG, gtilde_DG, htilde_DG
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
    if ( (alpha > alpha_maxIDP) .or. isnan(dalpha)) then
      alpha = alpha_maxIDP
      dalpha  = alpha_maxIDP - alphacont
    end if
          
#if LOCAL_ALPHA
    ! Change the alpha for output
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      alphacont_loc    = alpha_loc(i,j,k)
      alpha_loc(i,j,k) = alpha_loc(i,j,k) + dalpha_loc(i,j,k)
      if ( (alpha_loc(i,j,k) > alpha_maxIDP) .or. isnan(dalpha_loc(i,j,k))) then
        alpha_loc(i,j,k) = alpha_maxIDP
        dalpha_loc (i,j,k) = alpha_maxIDP - alphacont_loc
      end if
    end do       ; end do       ; enddo
    
    ! Correct!
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      ! xi correction
      ! -------------
      ! left
      my_corr=-max(dalpha_loc(i-1,j  ,k  ),dalpha_loc(i  ,j  ,k  )) * sWGP(i) * (ftilde_DG(:,i-1,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! right
      my_corr=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i+1,j  ,k  )) * sWGP(i) * (ftilde_DG(:,i  ,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! eta correction
      ! --------------
      ! left
      my_corr=-max(dalpha_loc(i  ,j-1,k  ),dalpha_loc(i  ,j  ,k  )) * sWGP(j) * (gtilde_DG(:,i,j-1,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! right
      my_corr=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i  ,j+1,k  )) * sWGP(j) * (gtilde_DG(:,i,  j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! zeta correction
      ! ---------------
      ! left
      my_corr=-max(dalpha_loc(i  ,j  ,k-1),dalpha_loc(i  ,j  ,k  )) * sWGP(k) * (htilde_DG(:,i,j,k-1,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr * dt
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr
      
      ! right
      my_corr=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i  ,j  ,k+1)) * sWGP(k) * (htilde_DG(:,i,  j,k,eID))*sJ(i,j,k,eID)
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
  pure subroutine GetBarStates(UL,UR,nv,t1,t2,Ubar,lambdamax)
    use MOD_Riemann       , only: RotateState, RotateFluxBack, MaxEigenvalRiemann
    USE MOD_Flux          , only: EvalAdvectionFlux1D
    implicit none
    !-arguments----------------------------------------
    real, intent(in) :: UL  (PP_nVar)
    real, intent(in) :: UR  (PP_nVar)
    real, intent(in) :: nv  (3)
    real, intent(in) :: t1  (3)
    real, intent(in) :: t2  (3)
    real, intent(out):: Ubar(PP_nVar)
    real, intent(out):: lambdamax
    !-local-variables----------------------------------
    real :: UL_r(PP_nVar) ! Rotated left state
    real :: UR_r(PP_nVar) ! Rotated right state
    real :: FL_r(PP_nVar) ! Rotated left flux
    real :: FR_r(PP_nVar) ! Rotated right flux 
    !--------------------------------------------------
    
    UL_r = RotateState(UL,nv,t1,t2)
    UR_r = RotateState(UR,nv,t1,t2)
    
    lambdamax = MaxEigenvalRiemann(UL_r,UR_r)
    
    call EvalAdvectionFlux1D(UL_r, FL_r)
    call EvalAdvectionFlux1D(UR_r, FR_r)
    
    Ubar = 0.5*( UL_r + UR_r ) - (0.5/(lambdamax)) * (FR_r-FL_r)
    
    call RotateFluxBack(Ubar,nv,t1,t2)
  end subroutine GetBarStates
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

