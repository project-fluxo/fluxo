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
! * If FV_TIMESTEP is active and the bar states are computed, the allowable time-step is updated with the very restrictive
!   LLF CFL condition... Eq (40) of Pazner (2020) "Sparse Invariant Domain Preserving Discontinuous Galerkin Methods With Subcell Convex Limiting".
!==================================================================================================================================
#include "defines.h"
! Use bar states by default, except if the equation has non-conservative terms
!TODO: define  barstates in another manner!!!
#if NONCONS
#define barStates 0
#else
#define barStates 1
#endif /*NONCONS*/
! A switch to define if the output variable alpha_old contains a cumulative alpha (1), or the instant alpha before limiting (0). Must be changed in src/analyze/analyze.f90 as well!
#define cumulativeAlphaOld 1
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
    call prms%CreateLogicalOption(   "IDPStateTVD",  " IDP(TVD) correction on state quantities? (to be used in combination with IDPStateTVDVarsNum, IDPStateTVDVars, and IDPStateTVDeqWise)", "F")
    call prms%CreateLogicalOption( "IDPPositivity",  " IDP correction for positivity of density and pressure?", "F")
    
!   Additional options
!   ------------------
    ! For IDPStateTVD
    call prms%CreateIntOption(   "IDPStateTVDVarsNum",  " Number of state variables to impose TVD correction", "1")
    call prms%CreateIntArrayOption( "IDPStateTVDVars",  " Variables to impose TVD correction", "1")
    call prms%CreateLogicalOption("IDPStateTVDeqWise",  " Perform IDP(TVD) correction equation-wise?", "F")
    
    ! For IDPPositivity
    call prms%CreateIntOption(   "IDPPositiveVarsNum",  " Number of state variables to impose positivity correction", "1")
    call prms%CreateIntArrayOption( "IDPPositiveVars",  " Variables to impose positivity correction", "1")
    
    ! For IDPMathEntropy and IDPSpecEntropy
    call prms%CreateLogicalOption("IDPNonlinearIfState",  " Do IDP correction on nonlinear quantities only if there was limiting on state quantities?", "F")
    
    call prms%CreateLogicalOption( "IDPafterIndicator"," If true, the IDP limiters (except for IDPPositivity) are only used where shock indicator alpha>IDPalpha_min (and the indicator's alpha is neglected)", "F")
    call prms%CreateRealOption(  "IDPalpha_min"   ,  " Parameter for IDPafterIndicator=T (Default: alpha_min)")
    
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
    use MOD_Globals
    use MOD_PreProc
    use MOD_Mesh_Vars  , only: nElems, sJ, MeshIsNonConforming
#if !(barStates)
    use MOD_Mesh_Vars  , only: firstSlaveSide, LastSlaveSide, nSides
#endif /*!(barStates)*/
#if USE_AMR
    use MOD_AMR_Vars           , only: UseAMR
#endif /*USE_AMR*/
    use MOD_ReadInTools, only: GETINT, GETREAL, GETLOGICAL,GETINTARRAY
    use MOD_StringTools, only: REALTOSTR
    implicit none
    !-local-variables----------------------------------------
    integer :: i,j,k,eID
    logical :: MeshNonConforming
    character(len=255) :: VarName
    !--------------------------------------------------------
    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' IDP Methods: '
    
!   Get parameters
!   --------------
    
    ! IDP limiters
    IDPPositivity  = GETLOGICAL('IDPPositivity' ,'F')
    IDPStateTVD    = GETLOGICAL('IDPStateTVD' ,'F')
    IDPMathEntropy = GETLOGICAL('IDPMathEntropy','F')
    IDPSpecEntropy = GETLOGICAL('IDPSpecEntropy','F')
    
    ! For IDPStateTVD
    if (IDPStateTVD) then
      IDPStateTVDeqWise  = GETLOGICAL('IDPStateTVDeqWise' ,'F')
      IDPStateTVDVarsNum = GETINT('IDPStateTVDVarsNum','1')
      IDPStateTVDVars    = GETINTARRAY('IDPStateTVDVars',IDPStateTVDVarsNum,'1')
    end if
    
    ! For IDPPositivity
    if (IDPPositivity) then
      IDPPositiveVarsNum = GETINT('IDPPositiveVarsNum','1')
      IDPPositiveVars    = GETINTARRAY('IDPPositiveVars',IDPPositiveVarsNum,'1')
    end if
    
    ! Specifications for IDPMathEntropy and IDPMathEntropy
    IDPNonlinearIfState = GETLOGICAL('IDPNonlinearIfState','F')
    
    IDPafterIndicator = GETLOGICAL('IDPafterIndicator','F')
    IDPalpha_min      = GETREAL   ('IDPalpha_min',REALTOSTR(alpha_min))
    
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
    
!   Some IDP methods are not implemented -yet- for nonconforming meshes
!   -------------------------------------------------------------------
    
    ! Check if the mesh can be non-conforming
    MeshNonConforming = MeshIsNonConforming
#if USE_AMR
    MeshNonConforming = UseAMR .or. MeshNonConforming
#endif /*USE_AMR*/
    
    if ( MeshNonConforming .and. (IDPStateTVD .or. IDPMathEntropy .or. IDPSpecEntropy) ) then
      CALL abort(__STAMP__,'IDPStateTVD/IDPMathEntropy/IDPSpecEntropy cannot be used with nonconforming meshes!',999,999.)
      RETURN
    end if
    
!   Internal definitions (all are .FALSE. by default)
!   -------------------------------------------------
    ! IDPneedsUprev made .TRUE. by default, such that the analyze step can compute dSdt after IDP. To avoid additional allocations, make false.
    IDPneedsUprev     = .TRUE.
    
    if (IDPPositivity ) then
      IDPneedsUsafe     = .TRUE.
    end if
    
    if (IDPStateTVD) then
#if barStates
#if LOCAL_ALPHA
      IDPneedsUsafe = .TRUE.
#endif /*LOCAL_ALPHA*/
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev     = .TRUE.
      IDPneedsUprev_ext = .TRUE.
#else
      IDPneedsUsafe = .TRUE.
      IDPneedsUsafe_ext = .TRUE.
#endif /*barStates*/
    end if
    
    if (IDPMathEntropy) then
      IDPneedsUsafe     = .TRUE.
#if barStates
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev_ext = .TRUE.
      IDPneedsUprev     = .TRUE.
#else
      IDPneedsUsafe_ext = .TRUE.
#endif /*barStates*/
    end if
    
    if (IDPSpecEntropy) then
      IDPneedsUsafe     = .TRUE.
#if barStates
      IDPneedsUbar      = .TRUE.
      IDPneedsUprev_ext = .TRUE.
      IDPneedsUprev     = .TRUE.
#else
      IDPneedsUsafe_ext = .TRUE.
#endif /*barStates*/
    end if
       
!   Allocate storage
!   ----------------
    ! Alpha before limiting
#if LOCAL_ALPHA
    allocate ( alpha_old(0:PP_N,0:PP_N,0:PP_N,nElems) )
#else
    allocate ( alpha_old(nElems) )
#endif /*LOCAL_ALPHA*/
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
    if (IDPneedsUsafe_ext) then
      allocate( Usafe_ext    (PP_nVar, 0:PP_N  , 0:PP_N          ,6,nElems) )
    end if
#endif /*barStates*/
    
    ! Variables for local alpha
#if LOCAL_ALPHA
    allocate ( alpha_loc(0:PP_N,0:PP_N,0:PP_N,1:nElems) )
    allocate ( f_antidiff(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N,nElems) )
    allocate ( g_antidiff(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N,nElems) )
    allocate ( h_antidiff(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N,nElems) )
    f_antidiff = 0.0
    g_antidiff = 0.0
    h_antidiff = 0.0
#if NONCONS
    allocate ( f_antidiffR(PP_nVar,-1:PP_N, 0:PP_N, 0:PP_N,nElems) )
    allocate ( g_antidiffR(PP_nVar, 0:PP_N,-1:PP_N, 0:PP_N,nElems) )
    allocate ( h_antidiffR(PP_nVar, 0:PP_N, 0:PP_N,-1:PP_N,nElems) )
    f_antidiffR = 0.0
    g_antidiffR = 0.0
    h_antidiffR = 0.0
#endif /*NONCONS*/
    allocate ( dalpha_loc     (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1) )
#endif /*LOCAL_ALPHA*/
    
    ! Bounds containers
    if (IDPStateTVD .or. IDPPositivity) then
      allocate ( state_min   (PP_nVar,0:PP_N,0:PP_N,0:PP_N) )
    end if
    if (IDPStateTVD) then
      allocate ( state_max   (PP_nVar,0:PP_N,0:PP_N,0:PP_N) )
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
    
    ! Now the analyze bounds
    idp_bounds_num = 0
    
    if (IDPStateTVD) then
      do i=1, IDPStateTVDVarsNum
        write(VarName,'(A,I0)') 'u', IDPStateTVDVars(i)
        idp_bounds_num = idp_bounds_num + 1
        idp_bounds_names(idp_bounds_num) = trim(VarName)//'_min'
        idp_bounds_num = idp_bounds_num + 1
        idp_bounds_names(idp_bounds_num) = trim(VarName)//'_max'
      end do
    end if
    
    if (IDPPositivity) then
      do i=1, IDPPositiveVarsNum
        if (.not. any(IDPStateTVDVars==IDPPositiveVars(i))) then
          write(VarName,'(A,I0)') 'u', IDPPositiveVars(i)
          idp_bounds_num = idp_bounds_num + 1
          idp_bounds_names(idp_bounds_num) = trim(VarName)//'_min'
        end if
      end do
    end if
    
    if (IDPSpecEntropy) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'ent_min'
    end if
    
    if (IDPMathEntropy) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'ent_max'
    end if
    
    if (IDPPositivity) then
      idp_bounds_num = idp_bounds_num + 1
      idp_bounds_names(idp_bounds_num) = 'p_min'
    end if
    
    idp_bounds_delta = 0.0
    
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
    use MOD_PreProc
    use MOD_Mesh_Vars   , only: nElems
    use MOD_NFVSE_Vars  , only: alpha, alpha_old, maximum_alpha, amount_alpha, amount_alpha_steps, alpha_min
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars  , only: alpha_loc
    use MOD_Analyze_Vars, only: wGPVol
    use MOD_Mesh_Vars   , only: nElems, sJ
    use MOD_IDP_Vars    , only: dalpha_loc
    use MOD_NFVSE_Vars  , only: f_antidiff, g_antidiff, h_antidiff
#if NONCONS
    use MOD_NFVSE_Vars  , only: f_antidiffR, g_antidiffR, h_antidiffR
#endif /*NONCONS*/
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars    , only: IDPSpecEntropy, IDPMathEntropy, IDPStateTVD, IDPPositivity, dalpha
    use MOD_IDP_Vars    , only: IDPForce2D, FFV_m_FDG, IDPafterIndicator, IDPalpha_min
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
#if LOCAL_ALPHA
    real :: alphaold (0:PP_N,0:PP_N,0:PP_N,nElems)  ! TODO: allocate beforehand!!
#else
    real :: alphaold(nElems)
#endif /*LOCAL_ALPHA*/
    logical :: doIDP
    !--------------------------------------------------------
    
!   Initialize
!   ----------
#if LOCAL_ALPHA
    alphaold = alpha_loc
#else
    alphaold = alpha
#endif /*LOCAL_ALPHA*/
    
    sdt = 1./dt
    doIDP = .not. IDPafterIndicator
    
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
      call ResetBounds()
!     Enforce 2D condition
!     --------------------
      if (IDPForce2D) then
        do k=1, PP_N ; do j=0, PP_N ; do i=0, PP_N
          FFV_m_FDG(:,i,j,k,eID) = FFV_m_FDG(:,i,j,0,eID)
        end do       ; end do       ; end do
#if LOCAL_ALPHA
        do k=1, PP_N ; do j=0, PP_N ; do i=-1, PP_N
          f_antidiff(:,i,j,k,eID) = f_antidiff(:,i,j,0,eID)
#if NONCONS
          f_antidiffR(:,i,j,k,eID) = f_antidiffR(:,i,j,0,eID)
#endif /*NONCONS*/
        end do       ; end do       ; end do
        do k=1, PP_N ; do j=-1, PP_N ; do i=0, PP_N
          g_antidiff(:, i,j,k,eID) = g_antidiff(:, i,j,0,eID)
#if NONCONS
          g_antidiffR(:, i,j,k,eID) = g_antidiffR(:, i,j,0,eID)
#endif /*NONCONS*/
        end do       ; end do       ; end do
#endif /*LOCAL_ALPHA*/  
      end if
!     Apply IDP limiters only where indicator fires?
!     ----------------------------------------------
      if (IDPafterIndicator) then
        if (alpha(eID) >= IDPalpha_min) then
          doIDP = .TRUE.
        else
          doIDP = .FALSE.
        end if
        ! Reverse blending
        U (:,:,:,:,eID) = U (:,:,:,:,eID) - alpha(eID) * FFV_m_FDG(:,:,:,:,eID) * dt
        Ut(:,:,:,:,eID) = Ut(:,:,:,:,eID) - alpha(eID) * FFV_m_FDG(:,:,:,:,eID)
        ! Zero alpha
        alpha(eID) = 0.0
#if LOCAL_ALPHA
        alpha_loc(:,:,:,eID) = 0.0
        alphaold (:,:,:,eID) = 0.0
#else
        alphaold(eID) = 0.0
#endif /*LOCAL_ALPHA*/  
      end if
      
!     Call all user-defined limiters to obtain dalpha and set bounds
!     --------------------------------------------------------------
      if (doIDP) then
        if (IDPStateTVD)    call IDP_LimitStateTVD   (U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
        if (IDPSpecEntropy) call IDP_LimitSpecEntropy(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
        if (IDPMathEntropy) call IDP_LimitMathEntropy(U(:,:,:,:,eID),Ut(:,:,:,:,eID),dt,sdt,eID)
      end if
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
      call CheckBounds(U(:,:,:,:,eID),eID)
    end do
    
!   Update variables for the analyze routines
!   -----------------------------------------
#if LOCAL_ALPHA
    maximum_alpha = max(maximum_alpha,maxval(alpha_loc-alphaold))
#else
    maximum_alpha = max(maximum_alpha,maxval(alpha-alphaold))
#endif /*LOCAL_ALPHA*/
    
    amount_alpha = amount_alpha*amount_alpha_steps
#if LOCAL_ALPHA
    curr_amount_alpha = 0.0
    do eID=1, nElems
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        curr_amount_alpha = curr_amount_alpha + (alpha_loc(i,j,k,eID)-alphaold(i,j,k,eID))*wGPVol(i,j,k)/sJ(i,j,k,eID)
      end do       ; end do       ; end do
    end do
    amount_alpha = amount_alpha + curr_amount_alpha
#else
    amount_alpha = amount_alpha + sum(alpha-alphaold)
#endif /*LOCAL_ALPHA*/
    amount_alpha_steps = amount_alpha_steps+1
    amount_alpha = amount_alpha/amount_alpha_steps
    
#if cumulativeAlphaOld
    do eID=1, nElems
#if LOCAL_ALPHA
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        alpha_old(i,j,k,eID) = max(alpha_old(i,j,k,eID),alpha_loc(i,j,k,eID))
      end do       ; end do       ; end do
#else
      alpha_old(eID) = max(alpha_old(eID),alpha(eID))
#endif /*LOCAL_ALPHA*/
    end do
#else
    alpha_old = alphaold
#endif /*cumulativeAlphaOld*/
    
  end subroutine Apply_IDP
!===================================================================================================================================
!> Check that all bounds are met
!===================================================================================================================================
  subroutine ResetBounds
    use MOD_Preproc
    use MOD_Globals
    use MOD_IDP_Vars      , only: state_min, state_max, s_min, s_max, p_min, p_max
    use MOD_IDP_Vars      , only: IDPStateTVD, IDPSpecEntropy, IDPMathEntropy, IDPPositivity, IDPForce2D
    implicit none
    !-arguments------------------------------------------------------------
    !----------------------------------------------------------------------
    
    if (IDPStateTVD .or. IDPPositivity) then
      state_min =-huge(1.0)
    end if
      
    if (IDPStateTVD) then
      state_max = huge(1.0)
    end if
    
    if (IDPSpecEntropy) then
      s_min =-huge(1.0)
    end if
    
    if (IDPMathEntropy) then
      s_max = huge(1.0)
    end if
    
    if (IDPPositivity) then
      p_min =0.0
    end if
  
  end subroutine ResetBounds
!===================================================================================================================================
!> Check that all bounds are met
!===================================================================================================================================
  subroutine CheckBounds(U,eID)
    use MOD_Preproc
    use MOD_Globals
    use MOD_IDP_Vars      , only: state_min, state_max, s_min, s_max, p_min, p_max, idp_bounds_delta
    use MOD_IDP_Vars      , only: IDPStateTVD, IDPSpecEntropy, IDPMathEntropy, IDPPositivity, IDPForce2D
    use MOD_IDP_Vars      , only: IDPStateTVDVars, IDPStateTVDVarsNum, IDPPositiveVars, IDPPositiveVarsNum
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
    integer :: i,j,k,counter,var
    real    :: p
    !----------------------------------------------------------------------
    
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      counter=0
      if (IDPStateTVD) then
        do var=1, IDPStateTVDVarsNum
          counter=counter+1
          if (IDPForce2D) state_min(IDPStateTVDVars(var),i,j,k) = state_min(IDPStateTVDVars(var),i,j,0)
          idp_bounds_delta(counter) = max(idp_bounds_delta(counter), state_min(IDPStateTVDVars(var),i,j,k) - U(IDPStateTVDVars(var),i,j,k))
          counter=counter+1
          if (IDPForce2D) state_max(IDPStateTVDVars(var),i,j,k) = state_max(IDPStateTVDVars(var),i,j,0)  
          idp_bounds_delta(counter) = max(idp_bounds_delta(counter), U(IDPStateTVDVars(var),i,j,k) - state_max(IDPStateTVDVars(var),i,j,k))
        end do
      end if
        
      if (IDPPositivity) then
        do var=1, IDPPositiveVarsNum
          if (.not. any(IDPStateTVDVars==IDPPositiveVars(var))) then
            counter=counter+1
            if (IDPForce2D) state_min(IDPPositiveVars(var),i,j,k) = state_min(IDPPositiveVars(var),i,j,0)
            idp_bounds_delta(counter) = max(idp_bounds_delta(counter), state_min(IDPPositiveVars(var),i,j,k) - U(IDPPositiveVars(var),i,j,k))
          end if
        end do
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
      
      if (IDPPositivity) then
        counter=counter+1
        if (IDPForce2D) p_min(i,j,k) = p_min(i,j,0)
        call Get_Pressure(U(:,i,j,k),p)
        idp_bounds_delta(counter) = max(idp_bounds_delta(counter),p_min(i,j,k) - p)
      end if
      
    end do       ; end do       ; end do ! i,j,k
  
  end subroutine CheckBounds
!===================================================================================================================================
!> Get the IDP variables in the right position to perform limiting
!> ATTENTION: 1) U_master and U_slave need to have the previous solution!
!>            2) If the bar states are *deactivated* and IDPneedsUsafe_ext, U_master and U_slave contain the safe solution when this 
!>               subroutine returns
!===================================================================================================================================
  subroutine Get_IDP_Variables(U,dt,tIn)
    use MOD_NFVSE_Vars, only: alpha
    use MOD_PreProc
    use MOD_IDP_Vars      , only: IDPneedsUprev_ext, IDPneedsUsafe, IDPneedsUsafe_ext, IDPalpha_min, IDPafterIndicator
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
    use MOD_IDP_Vars      , only: IDPStateTVD, IDPMathEntropy, IDPSpecEntropy
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
        if (IDPafterIndicator .and. (alpha(eID) < IDPalpha_min) ) cycle
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
        
#if FV_TIMESTEP
        ! Compute maximum time step
        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          inv_dt = ( sWGP(i) * (lambdamax_xi  (i-1,j  ,k  ) * SCM % xi   % norm(j,k,i-1) + lambdamax_xi  (i,j,k) * SCM % xi   % norm(j,k,i)) + &
                     sWGP(j) * (lambdamax_eta (i  ,j-1,k  ) * SCM % eta  % norm(i,k,j-1) + lambdamax_eta (i,j,k) * SCM % eta  % norm(i,k,j)) + &
                     sWGP(k) * (lambdamax_zeta(i  ,j  ,k-1) * SCM % zeta % norm(i,j,k-1) + lambdamax_zeta(i,j,k) * SCM % zeta % norm(i,j,k)) ) * sJ(i,j,k,eID)
          maxdt_IDP = min (maxdt_IDP, 1./inv_dt)
        end do       ; end do       ; end do
#endif /*FV_TIMESTEP*/
        end associate
      end do !eID
      
#if FV_TIMESTEP
#if MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,maxdt_IDP,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,iError)
#endif /*MPI*/
      if (dt > maxdt_IDP) then
        SWRITE(*,'(A,2ES21.12)') "MAYDAY, we have a problem with the time step. Consider reducing CFLScale (dt, maxdt_IDP): ", dt, maxdt_IDP
      end if
#endif /*FV_TIMESTEP*/
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
    if (IDPneedsUsafe_ext) then
!     Get the safe (FV) solution in place (if not using bar states)
!     * This overwrites U_master and U_slave (TODO: Change this?)
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
    end if !IDPneedsUsafe_ext
#endif /*!(barStates)*/
  end subroutine Get_IDP_Variables
!===================================================================================================================================
!> Density TVD correction
!===================================================================================================================================
  subroutine IDP_LimitStateTVD(U,Ut,dt,sdt,eID)
    use MOD_PreProc
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_IDP_Vars      , only: dalpha, alpha_maxIDP, IDPStateTVDeqWise, IDPStateTVDactive
    use MOD_Mesh_Vars     , only: nElems
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_NFVSE_Vars    , only: f_antidiff, g_antidiff, h_antidiff
#if NONCONS
    use MOD_NFVSE_Vars    , only: f_antidiffR, g_antidiffR, h_antidiffR
#endif /*NONCONS*/
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
    use MOD_IDP_Vars      , only: FFV_m_FDG, state_min, state_max, IDPStateTVDVarsNum, IDPStateTVDVars
    use MOD_DG_Vars       , only: Source
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    integer,intent(in) :: eID
    !-local-variables----------------------------------------
    real    :: dalpha1
    real    :: u_safe
    real    :: a   ! a  = PositCorrFactor * u_safe - rho
    real    :: Qp, Qm, Pp, Pm
#if LOCAL_ALPHA
    real    :: dalpha_locState(-1:PP_N+1,-1:PP_N+1,-1:PP_N+1)  ! Local copy of dalpha_loc
    real    :: alpha_locState(0:PP_N,0:PP_N,0:PP_N)   ! Local copy of alpha_loc
#endif /*LOCAL_ALPHA*/
    real    :: dalphaState                            ! Local copy of dalpha
    real    :: alphaState                             ! Local copy of alpha
    integer :: i,j,k,l,var
    !--------------------------------------------------------
    
    IDPStateTVDactive = .FALSE.
    
#if LOCAL_ALPHA
    dalpha_locState = dalpha_loc
#endif /*LOCAL_ALPHA*/
    dalphaState = dalpha
    
!   Compute bounds and correction factors for each variable
!   -------------------------------------------------------
    do var=1, IDPStateTVDVarsNum
      associate (ivar => IDPStateTVDVars(var))
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          
        ! Get the limit states
        !*********************
        state_min(ivar,i,j,k) =  huge(1.0)
        state_max(ivar,i,j,k) = -huge(1.0)
        
#if barStates
        ! Previous sol
        state_min(ivar,i,j,k) = min(state_min(ivar,i,j,k), Uprev  (ivar,i  ,j  ,k  ,eID))
        state_max(ivar,i,j,k) = max(state_max(ivar,i,j,k), Uprev  (ivar,i  ,j  ,k  ,eID))
        
        ! Source term
        state_min(ivar,i,j,k) = min(state_min(ivar,i,j,k), Uprev  (ivar,i  ,j  ,k  ,eID) + 2.0 * dt * Source(ivar, i, j, k, eID))
        state_max(ivar,i,j,k) = max(state_max(ivar,i,j,k), Uprev  (ivar,i  ,j  ,k  ,eID) + 2.0 * dt * Source(ivar, i, j, k, eID))
        
        !xi
        do l=i-1, i !min(i-1,0), max(i,PP_N-1)
          state_min(ivar,i,j,k) = min(state_min(ivar,i,j,k), Ubar_xi  (ivar,l  ,j  ,k  ,eID))
          state_max(ivar,i,j,k) = max(state_max(ivar,i,j,k), Ubar_xi  (ivar,l  ,j  ,k  ,eID))
        end do
        !eta
        do l=j-1, j !l=max(j-1,0), min(j,PP_N-1)
          state_min(ivar,i,j,k) = min(state_min(ivar,i,j,k), Ubar_eta (ivar,i  ,l  ,k  ,eID))
          state_max(ivar,i,j,k) = max(state_max(ivar,i,j,k), Ubar_eta (ivar,i  ,l  ,k  ,eID))
        end do
        if (.not. IDPForce2D) then
        !zeta
        do l=k-1, k !l=max(k-1,0), min(k,PP_N-1)
          state_min(ivar,i,j,k) = min(state_min(ivar,i,j,k), Ubar_zeta(ivar,i  ,j  ,l  ,eID))
          state_max(ivar,i,j,k) = max(state_max(ivar,i,j,k), Ubar_zeta(ivar,i  ,j  ,l  ,eID))
        end do
        end if
#else
        ! check stencil in xi
!          do l = i+idx_m1(i), i+idx_p1(i) !no neighbor
        do l = i-1, i+1
          state_min(ivar,i,j,k) = min(state_min(ivar,i,j,k), Usafe(ivar,l,j,k,eID))
          state_max(ivar,i,j,k) = max(state_max(ivar,i,j,k), Usafe(ivar,l,j,k,eID))
        end do
        ! check stencil in eta
!          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          state_min(ivar,i,j,k) = min(state_min(ivar,i,j,k), Usafe(ivar,i,l,k,eID))
          state_max(ivar,i,j,k) = max(state_max(ivar,i,j,k), Usafe(ivar,i,l,k,eID))
        end do
        ! check stencil in zeta
!          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          state_min(ivar,i,j,k) = min(state_min(ivar,i,j,k), Usafe(ivar,i,j,l,eID))
          state_max(ivar,i,j,k) = max(state_max(ivar,i,j,k), Usafe(ivar,i,j,l,eID))
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
        Qp = max(0.0,(state_max(ivar,i,j,k)-Usafe(ivar,i,j,k,eID))*sdt)
        Qm = min(0.0,(state_min(ivar,i,j,k)-Usafe(ivar,i,j,k,eID))*sdt)
        
        ! Positive contributions
        Pp = 0.0
        Pp = Pp + max(0.0, sWGP(i) * f_antidiff(ivar,i-1,j  ,k  ,eID))
        Pp = Pp + max(0.0, sWGP(j) * g_antidiff(ivar,i  ,j-1,k  ,eID))
#if NONCONS
        Pp = Pp + max(0.0,-sWGP(i) * f_antidiffR(ivar,i  ,j  ,k  ,eID))
        Pp = Pp + max(0.0,-sWGP(j) * g_antidiffR(ivar,i  ,j  ,k  ,eID))
#else
        Pp = Pp + max(0.0,-sWGP(i) * f_antidiff(ivar,i  ,j  ,k  ,eID))
        Pp = Pp + max(0.0,-sWGP(j) * g_antidiff(ivar,i  ,j  ,k  ,eID))
#endif /*NONCONS*/
        if (.not. IDPForce2D) then
        Pp = Pp + max(0.0, sWGP(k) * h_antidiff(ivar,i  ,j  ,k-1,eID))
#if NONCONS
        Pp = Pp + max(0.0,-sWGP(k) * h_antidiffR(ivar,i  ,j  ,k  ,eID))
#else
        Pp = Pp + max(0.0,-sWGP(k) * h_antidiff(ivar,i  ,j  ,k  ,eID))
#endif /*NONCONS*/
        end if
        Pp = Pp*sJ(i,j,k,eID)
        
        ! Negative contributions
        Pm = 0.0
        Pm = Pm + min(0.0, sWGP(i) * f_antidiff(ivar,i-1,j  ,k  ,eID))
        Pm = Pm + min(0.0, sWGP(j) * g_antidiff(ivar,i  ,j-1,k  ,eID))
#if NONCONS
        Pm = Pm + min(0.0,-sWGP(i) * f_antidiffR(ivar,i  ,j  ,k  ,eID))
        Pm = Pm + min(0.0,-sWGP(j) * g_antidiffR(ivar,i  ,j  ,k  ,eID))
#else
        Pm = Pm + min(0.0,-sWGP(i) * f_antidiff(ivar,i  ,j  ,k  ,eID))
        Pm = Pm + min(0.0,-sWGP(j) * g_antidiff(ivar,i  ,j  ,k  ,eID))
#endif /*NONCONS*/
        if (.not. IDPForce2D) then
        Pm = Pm + min(0.0, sWGP(k) * h_antidiff(ivar,i  ,j  ,k-1,eID))
#if NONCONS
        Pm = Pm + min(0.0,-sWGP(k) * h_antidiffR(ivar,i  ,j  ,k  ,eID))
#else
        Pm = Pm + min(0.0,-sWGP(k) * h_antidiff(ivar,i  ,j  ,k  ,eID))
#endif /*NONCONS*/
        end if
        Pm = Pm*sJ(i,j,k,eID)
        
        ! Compute blending coefficient avoiding division by zero
        ! (as in paper of [Guermond, Nazarov, Popov, Thomas] (4.8))
        Qp = (abs(Qp))/(abs(Pp) + epsilon(1.0)*100.*abs(state_max(ivar,i,j,k)))
        Qm = (abs(Qm))/(abs(Pm) + epsilon(1.0)*100.*abs(state_max(ivar,i,j,k)))
        
        ! Compute correction as: (needed_alpha) - current_alpha = (1.0 - min(1.0,Qp,Qm)) - alpha_loc(i,j,k,eID)
        dalpha1 = 1.0 - min(1.0,Qp,Qm) - alpha_loc(i,j,k,eID)
        
        dalpha_locState(i,j,k) = max(dalpha_locState(i,j,k),dalpha1)
        dalphaState = max(dalphaState,dalpha1)
        
#else
        ! Simple element-wise limiter
        !****************************
        if ( U(ivar,i,j,k) < state_min(ivar,i,j,k)) then
          u_safe = state_min(ivar,i,j,k)
        elseif (U(ivar,i,j,k) > state_max(ivar,i,j,k)) then
          u_safe = state_max(ivar,i,j,k)
        else
          cycle !nothing to do here!
        end if
        
        if ( abs(FFV_m_FDG(ivar,i,j,k,eID)) == 0.0) cycle !nothing to do here!
        
        ! Density correction
        a = (u_safe - U(ivar,i,j,k))
        dalpha1 = a*sdt / FFV_m_FDG(ivar,i,j,k,eID)
        
        ! Change inconsistent alphas
        if ( (alpha(eID)+dalpha1 > alpha_maxIDP) .or. isnan(dalpha1)) then
          dalphaState  = alpha_maxIDP - alpha(eID)
        else
          dalphaState = max(dalphaState,dalpha1)
        end if
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; end do ! i,j,k
      
      ! Check if TVD limiter had to act
      if (dalphaState > 0.0) IDPStateTVDactive = .TRUE.
      
      ! Perform correction equation wise if needed
      if (IDPStateTVDeqWise) then
        if ( dalphaState>0.0 .or. isnan(dalphaState) &
#if LOCAL_ALPHA
                             .or. any(isnan(dalpha_locState))   &
#endif /*LOCAL_ALPHA*/
                        ) then
          
          ! Fill in containers for alpha and alpha_loc
          alphaState = alpha(eID)
#if LOCAL_ALPHA
          alpha_locState = alpha_loc(:,:,:,eID)
#endif /*LOCAL_ALPHA*/
          
          ! Perform correction
          call PerformCorrection(U,Ut,dalphaState    ,alphaState    , &
#if LOCAL_ALPHA 
                                      dalpha_locState,alpha_locState, &
#endif /*LOCAL_ALPHA*/
                                      dt,sdt,eID,ivar,ivar)
          
          ! Recompute antidiffusive fluxes
#if LOCAL_ALPHA 
          do i=0, PP_N-1 ! i goes through the inner interfaces
            do k=0, PP_N  ; do j=0, PP_N ! j and k go through DOFs
              f_antidiff(ivar,i,j,k,eID) = (1.0 - max(dalpha_locState(i,j,k),dalpha_locState(i+1,j,k))) * f_antidiff(ivar,i,j,k,eID)
              g_antidiff(ivar,j,i,k,eID) = (1.0 - max(dalpha_locState(j,i,k),dalpha_locState(j,i+1,k))) * g_antidiff(ivar,j,i,k,eID)
              h_antidiff(ivar,j,k,i,eID) = (1.0 - max(dalpha_locState(j,k,i),dalpha_locState(j,k,i+1))) * h_antidiff(ivar,j,k,i,eID)
#if NONCONS
              f_antidiffR(ivar,i,j,k,eID) = (1.0 - max(dalpha_locState(i,j,k),dalpha_locState(i+1,j,k))) * f_antidiffR(ivar,i,j,k,eID)
              g_antidiffR(ivar,j,i,k,eID) = (1.0 - max(dalpha_locState(j,i,k),dalpha_locState(j,i+1,k))) * g_antidiffR(ivar,j,i,k,eID)
              h_antidiffR(ivar,j,k,i,eID) = (1.0 - max(dalpha_locState(j,k,i),dalpha_locState(j,k,i+1))) * h_antidiffR(ivar,j,k,i,eID)
#endif /*NONCONS*/
            end do        ; end do
          end do
#else
          FFV_m_FDG(ivar,:,:,:,eID) = (1.0 - dalphaState) * FFV_m_FDG(ivar,:,:,:,eID)
#endif /*LOCAL_ALPHA*/
          
          ! Restore dalphaState and dalpha_locState
#if LOCAL_ALPHA
          dalpha_locState = dalpha_loc
#endif /*LOCAL_ALPHA*/
          dalphaState = dalpha
        end if
      end if
      
      end associate
    end do !var 
    
    ! Modify variables for element
#if LOCAL_ALPHA
    dalpha_loc = dalpha_locState
#endif /*LOCAL_ALPHA*/
    dalpha = dalphaState
     
  end subroutine IDP_LimitStateTVD
!===================================================================================================================================
!> Specific entropy correction (discrete local minimum principle)
!===================================================================================================================================
  subroutine IDP_LimitSpecEntropy(U,Ut,dt,sdt,eID)
    use MOD_PreProc
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_IDP_Vars      , only: dalpha
    use MOD_IDP_Vars      , only: IDPNonlinearIfState, IDPStateTVDactive
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: dalpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: Usafe, IDPForce2D
    use MOD_IDP_Vars      , only: IDPparam_t
    use MOD_Mesh_Vars     , only: nElems
    use MOD_Equation_Vars , only: Get_SpecEntropy
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#endif /*barStates*/
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: IDPMaxIter, s_min
    use MOD_DG_Vars       , only: Source
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
    
    if (IDPNonlinearIfState .and. (.not. IDPStateTVDactive) ) return
    
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
        
        ! Source term
        s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Uprev    (:,i  ,j  ,k  ,eID) + 2.0 * dt * Source(:, i, j, k, eID)) )
        
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
          s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Usafe(:,l,j,k,eID)))
        end do
        ! check stencil in eta
!#          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Usafe(:,i,l,k,eID)))
        end do
        ! check stencil in zeta
!#          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          s_min(i,j,k) = min(s_min(i,j,k), Get_SpecEntropy(Usafe(:,i,j,l,eID)))
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
    use MOD_PreProc
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_IDP_Vars      , only: dalpha
    use MOD_IDP_Vars      , only: IDPNonlinearIfState, IDPStateTVDactive
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: dalpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: IDPparam_t, Usafe
    use MOD_Mesh_Vars     , only: nElems
    use MOD_Equation_Vars , only: Get_MathEntropy, ConsToEntropy
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#endif /*barStates*/
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: IDPMaxIter, s_max
    use MOD_DG_Vars       , only: Source
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
    
    if (IDPNonlinearIfState .and. (.not. IDPStateTVDactive) ) return
    
!     Compute correction factors
!     --------------------------
      notInIter = .FALSE.
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Get the limit states
        !*********************
        s_max(i,j,k) = -huge(1.0)
        
#if barStates
        ! Previous entropy of the node (ubar_ii)
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Uprev  (:,i  ,j  ,k  ,eID)))
        
        ! Source term
        s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Uprev  (:,i  ,j  ,k  ,eID) + 2.0 * dt * Source(:, i, j, k, eID)))
        
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
          s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Usafe(:,l,j,k,eID)))
        end do
        ! check stencil in eta
!          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Usafe(:,i,l,k,eID)))
        end do
        ! check stencil in zeta
!          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          s_max(i,j,k) = max(s_max(i,j,k), Get_MathEntropy(Usafe(:,i,j,l,eID)))
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
    use MOD_PreProc
    use MOD_NFVSE_Vars    , only: alpha, PositCorrFactor
    use MOD_Mesh_Vars     , only: nElems
    use MOD_IDP_Vars      , only: dalpha
    use MOD_IDP_Vars      , only: IDPPositiveVars, IDPPositiveVarsNum
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_IDP_Vars      , only: dalpha_loc
    use MOD_NFVSE_Vars    , only: f_antidiff, g_antidiff, h_antidiff
#if NONCONS
    use MOD_NFVSE_Vars    , only: f_antidiffR, g_antidiffR, h_antidiffR
#endif /*NONCONS*/
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: Usafe, p_safe, state_min, p_min, IDPStateTVD, IDPStateTVDVars, IDPForce2D
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
    real    :: a      ! a  = PositCorrFactor * state_safe - rho
    real    :: Qm, Pm ! Zalesak's limiter's variables
    integer :: i,j,k, iter
    logical :: NotInIter
    type(IDPparam_t) :: param ! Parameters for Newton's method
    real             :: new_alpha
    integer :: var
    !--------------------------------------------------------
      
!     ---------------
!     Correct density
!     ---------------
        
!     Compute correction factors
!     --------------------------
    do var=1, IDPPositiveVarsNum
      associate (ivar => IDPPositiveVars(var))
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Compute density bound
        !********************** 
        ! This writes the more restrictive bound into state_min
        if (IDPStateTVD .and. any(IDPStateTVDVars==ivar)) then
          state_min(ivar,i,j,k) = max(state_min(ivar,i,j,k), PositCorrFactor * Usafe(ivar,i,j,k,eID))
        else
          state_min(ivar,i,j,k) = PositCorrFactor * Usafe(ivar,i,j,k,eID)
        end if
#if LOCAL_ALPHA
        ! Real one-sided Zalesak-type limiter
        ! * Zalesak (1979). "Fully multidimensional flux-corrected transport algorithms for fluids"
        ! * Kuzmin et al. (2010). "Failsafe flux limiting and constrained data projections for equations of gas dynamics"
        ! ATTENTION: 1) The Zalesak limiter has to be computed, even if the state is valid, because the correction is 
        !               for each interface, not each node
        !****************************************************************************************************************
        
        ! Upper/lower bounds for admissible increments
        Qm = min(0.0,(state_min(ivar,i,j,k)-Usafe(ivar,i,j,k,eID))*sdt)
        
        ! Negative contributions
        Pm = 0.0
        Pm = Pm + min(0.0, sWGP(i) * f_antidiff(ivar,i-1,j  ,k  ,eID))
        Pm = Pm + min(0.0, sWGP(j) * g_antidiff(ivar,i  ,j-1,k  ,eID))
#if NONCONS
        Pm = Pm + min(0.0,-sWGP(i) * f_antidiffR(ivar,i  ,j  ,k  ,eID))
        Pm = Pm + min(0.0,-sWGP(j) * g_antidiffR(ivar,i  ,j  ,k  ,eID))
#else
        Pm = Pm + min(0.0,-sWGP(i) * f_antidiff(ivar,i  ,j  ,k  ,eID))
        Pm = Pm + min(0.0,-sWGP(j) * g_antidiff(ivar,i  ,j  ,k  ,eID))
#endif /*NONCONS*/
        if (.not. IDPForce2D) then
        Pm = Pm + min(0.0, sWGP(k) * h_antidiff(ivar,i  ,j  ,k-1,eID))
#if NONCONS
        Pm = Pm + min(0.0,-sWGP(k) * h_antidiffR(ivar,i  ,j  ,k  ,eID))
#else
        Pm = Pm + min(0.0,-sWGP(k) * h_antidiff(ivar,i  ,j  ,k  ,eID))
#endif /*NONCONS*/
        end if
        Pm = Pm*sJ(i,j,k,eID)
        
        ! Compute blending coefficient avoiding division by zero
        ! (modification of technique of [Guermond, Nazarov, Popov, Thomas] (4.8))
        Qm = (abs(Qm))/(abs(Pm) + epsilon(1.0)*100.)
        
        ! Compute correction as: (needed_alpha) - current_alpha = (1.0 - Qm) - alpha_loc(i,j,k,eID)
        dalpha1 = (1.0 - Qm) - alpha_loc(i,j,k,eID)
        
        dalpha_loc(i,j,k) = max(dalpha_loc(i,j,k),dalpha1)
        dalpha = max(dalpha,dalpha1)
        
#else
        a = (state_min(ivar,i,j,k) - U(ivar,i,j,k)) * sdt
        if (a > 0.) then ! This DOF needs a correction
          if ( abs(FFV_m_FDG(ivar,i,j,k,eID)) == 0.0) cycle !nothing to do here!
          dalpha1 = a / FFV_m_FDG(ivar,i,j,k,eID)
          dalpha = max(dalpha,dalpha1)
        end if

#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; end do ! i,j,k
      end associate
    end do
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
    use MOD_NFVSE_Vars    , only: f_antidiff, g_antidiff, h_antidiff
#if NONCONS
    use MOD_NFVSE_Vars    , only: f_antidiffR, g_antidiffR, h_antidiffR
#endif /*NONCONS*/
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
    param % F_antidiff =  IDPgamma * sJ(i,j,k,eID) * sWGP(i) * (f_antidiff(:,i-1,j  ,k  ,eID)) ! Anti-difussive flux in xi-
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    
    ! xi+
#if NONCONS
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(i) * (f_antidiffR(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in xi+
#else
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(i) * (f_antidiff(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in xi+
#endif /*NONCONS*/
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    ! eta-
    param % F_antidiff =  IDPgamma * sJ(i,j,k,eID) * sWGP(j) * (g_antidiff(:,i  ,j-1,k  ,eID)) ! Anti-difussive flux in eta-
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    
    ! eta+
#if NONCONS
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(j) * (g_antidiffR(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in eta+
#else
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(j) * (g_antidiff(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in eta+
#endif /*NONCONS*/
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    
    if (.not. IDPForce2D) then
    ! zeta-
    param % F_antidiff =  IDPgamma * sJ(i,j,k,eID) * sWGP(k) * (h_antidiff(:,i  ,j  ,k-1,eID)) ! Anti-difussive flux in zeta-
    call NewtonLoop(Usafe(:,i,j,k,eID), &  ! Safe (FV) solution
                              param,alpha,notInIter,Goal,dGoal_dbeta,InitialCheck,FinalCheck)
    
    ! zeta+
#if NONCONS
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(k) * (h_antidiffR(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in zeta+
#else
    param % F_antidiff = -IDPgamma * sJ(i,j,k,eID) * sWGP(k) * (h_antidiff(:,i  ,j  ,k  ,eID)) ! Anti-difussive flux in zeta+
#endif /*NONCONS*/
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
    use MOD_Equation_Vars , only: StateIsValid
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
    
    ! If the state is valid, perform initial check and return if correction is not needed
    if (StateIsValid(Ucurr)) then
      as = Goal(param,Ucurr)
      if (InitialCheck(param,as)) return
    end if
    
    ! Perform Newton iterations
    TheNewtonLoop: do iter=1, IDPMaxIter
      beta_old = beta
      
      ! If the state is valid, evaluate d(goal)/d(beta)
      if (StateIsValid(Ucurr)) then
        dSdbeta = dGoal_dbeta(param,Ucurr)
      else ! Otherwise, perform a bisection step
        dSdbeta= 0.0
      end if
      
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
        
        ! If the state is invalid, finish bisection step without checking tolerance and iterate further
        if (.not. StateIsValid(Ucurr)) then
          beta_R = beta
          cycle
        end if
        
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
        ! If the state is invalid, redefine right bound without checking tolerance and iterate further
        if (.not. StateIsValid(Ucurr)) then
          beta_R = beta
          cycle
        end if !debug
        
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
                                                             dt,sdt,eID,vs_in,ve_in)
    use MOD_PreProc
    use MOD_IDP_Vars  , only: alpha_maxIDP
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars, only: f_antidiff, g_antidiff, h_antidiff
#if NONCONS
    use MOD_NFVSE_Vars, only: f_antidiffR, g_antidiffR, h_antidiffR
#endif /*NONCONS*/
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
    integer, intent(in), optional :: vs_in,ve_in           ! Start and end variables to perform correction
    !-local-variables--------------------------------
    real :: alphacont_loc ! A local container
    real :: alphacont     ! An element-wise container
    real :: my_corr(PP_nVar)
    integer :: i,j,k
    integer :: vs, ve
    !------------------------------------------------
    
    if ( present(vs_in) .and. present(ve_in)) then
      vs = vs_in
      ve = ve_in
    else
      vs = 1
      ve = PP_nVar
    end if
    
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
      my_corr(vs:ve)=-max(dalpha_loc(i-1,j  ,k  ),dalpha_loc(i  ,j  ,k  )) * sWGP(i) * (f_antidiff(vs:ve,i-1,j,k,eID))*sJ(i,j,k,eID)
      U (vs:ve,i,j,k) = U (vs:ve,i,j,k) + my_corr(vs:ve) * dt
      Ut(vs:ve,i,j,k) = Ut(vs:ve,i,j,k) + my_corr(vs:ve)
      
      ! right
#if NONCONS
      my_corr(vs:ve)=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i+1,j  ,k  )) * sWGP(i) * (f_antidiffR(vs:ve,i  ,j,k,eID))*sJ(i,j,k,eID)
#else
      my_corr(vs:ve)=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i+1,j  ,k  )) * sWGP(i) * (f_antidiff (vs:ve,i  ,j,k,eID))*sJ(i,j,k,eID)
#endif /*NONCONS*/
      U (vs:ve,i,j,k) = U (vs:ve,i,j,k) + my_corr(vs:ve) * dt
      Ut(vs:ve,i,j,k) = Ut(vs:ve,i,j,k) + my_corr(vs:ve)
      
      ! eta correction
      ! --------------
      ! left
      my_corr(vs:ve)=-max(dalpha_loc(i  ,j-1,k  ),dalpha_loc(i  ,j  ,k  )) * sWGP(j) * (g_antidiff(vs:ve,i,j-1,k,eID))*sJ(i,j,k,eID)
      U (vs:ve,i,j,k) = U (vs:ve,i,j,k) + my_corr(vs:ve) * dt
      Ut(vs:ve,i,j,k) = Ut(vs:ve,i,j,k) + my_corr(vs:ve)
      
      ! right
#if NONCONS
      my_corr(vs:ve)=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i  ,j+1,k  )) * sWGP(j) * (g_antidiffR(vs:ve,i,  j,k,eID))*sJ(i,j,k,eID)
#else
      my_corr(vs:ve)=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i  ,j+1,k  )) * sWGP(j) * (g_antidiff (vs:ve,i,  j,k,eID))*sJ(i,j,k,eID)
#endif /*NONCONS*/
      U (vs:ve,i,j,k) = U (vs:ve,i,j,k) + my_corr(vs:ve) * dt
      Ut(vs:ve,i,j,k) = Ut(vs:ve,i,j,k) + my_corr(vs:ve)
      
      ! zeta correction
      ! ---------------
      ! left
      my_corr(vs:ve)=-max(dalpha_loc(i  ,j  ,k-1),dalpha_loc(i  ,j  ,k  )) * sWGP(k) * (h_antidiff(vs:ve,i,j,k-1,eID))*sJ(i,j,k,eID)
      U (vs:ve,i,j,k) = U (vs:ve,i,j,k) + my_corr(vs:ve) * dt
      Ut(vs:ve,i,j,k) = Ut(vs:ve,i,j,k) + my_corr(vs:ve)
      
      ! right
#if NONCONS
      my_corr(vs:ve)=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i  ,j  ,k+1)) * sWGP(k) * (h_antidiffR(vs:ve,i,  j,k,eID))*sJ(i,j,k,eID)
#else
      my_corr(vs:ve)=max(dalpha_loc(i  ,j  ,k  ),dalpha_loc(i  ,j  ,k+1)) * sWGP(k) * (h_antidiff (vs:ve,i,  j,k,eID))*sJ(i,j,k,eID)
#endif /*NONCONS*/
      U (vs:ve,i,j,k) = U (vs:ve,i,j,k) + my_corr(vs:ve) * dt
      Ut(vs:ve,i,j,k) = Ut(vs:ve,i,j,k) + my_corr(vs:ve)
    end do       ; end do       ; enddo
#else
    ! Element-wise correction!
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      ! Correct U
      U (vs:ve,i,j,k) = U (vs:ve,i,j,k) + dalpha * dt * FFV_m_FDG(vs:ve,i,j,k,eID)
      ! Correct Ut
      Ut(vs:ve,i,j,k) = Ut(vs:ve,i,j,k) + dalpha * FFV_m_FDG(vs:ve,i,j,k,eID)
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
    
    SDEALLOCATE( Usafe_ext )
    SDEALLOCATE( Uprev_ext )
    SDEALLOCATE( Ubar_xi  )
    SDEALLOCATE( Ubar_eta  )
    SDEALLOCATE( Ubar_zeta  )
    
#if LOCAL_ALPHA
    SDEALLOCATE ( alpha_loc )
    SDEALLOCATE ( f_antidiff )
    SDEALLOCATE ( g_antidiff )
    SDEALLOCATE ( h_antidiff )
#if NONCONS
    SDEALLOCATE ( f_antidiffR )
    SDEALLOCATE ( g_antidiffR )
    SDEALLOCATE ( h_antidiffR )
#endif /*NONCONS*/
    SDEALLOCATE ( dalpha_loc )
#endif /*LOCAL_ALPHA*/
    
    SDEALLOCATE ( state_min )
    SDEALLOCATE ( state_max )
    SDEALLOCATE ( s_min )
    SDEALLOCATE ( s_max )
    SDEALLOCATE ( p_min )
    SDEALLOCATE ( p_max )
    
  end subroutine Finalize_IDP
#endif /*NFVSE_CORR*/  
end module MOD_IDP

