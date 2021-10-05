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
    call prms%CreateLogicalOption("IDPSemiDiscEnt",  " IDP correction on semi-discrete entropy balance?", "F")
    call prms%CreateLogicalOption( "IDPDensityTVD",  " IDP(TVD) correction on density? (uses a Zalesak limiter with LOCAL_ALPHA=ON)", "F")
    call prms%CreateLogicalOption( "IDPPositivity",  " IDP correction for positivity of density and pressure?", "F")
    
!   Additional options
!   ------------------
    call prms%CreateRealOption(  "PositCorrFactor",  " The correction factor for IDPPositivity=T", "0.1")
    call prms%CreateIntOption(        "IDPMaxIter",  " Maximum number of iterations for positivity limiter", "10")
    
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
    IDPSemiDiscEnt = GETLOGICAL('IDPSemiDiscEnt','F')
    
    ! Consistency check
#if !defined(LOCAL_ALPHA)
    if (IDPSemiDiscEnt) stop 'IDPSemiDiscEnt needs LOCAL_ALPHA'
#endif
    if (IDPMathEntropy .and. IDPSpecEntropy) then
      stop 'Only one of the two can be selected: IDPMathEntropy/IDPSpecEntropy'
    end if
    
    ! Additional options
    IDPMaxIter      = GETINT    ('IDPMaxIter','10')
    if (IDPPositivity) then
      PositCorrFactor = GETREAL   ('PositCorrFactor','0.1')
    end if
    
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
    end if
    
    if (IDPSemiDiscEnt) then
      IDPneedsUprev     = .TRUE.
      IDPneedsUprev_ext = .TRUE.
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
      Uprev = 0.0
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
    if (IDPSemiDiscEnt) then
      allocate( Flux_ext     (PP_nVar, 0:PP_N  , 0:PP_N          ,6,nElems) )
    end if
    
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
    
    allocate ( rf_DG(-1:PP_N, 0:PP_N, 0:PP_N,1:nElems) )
    allocate ( rg_DG( 0:PP_N,-1:PP_N, 0:PP_N,1:nElems) )
    allocate ( rh_DG( 0:PP_N, 0:PP_N,-1:PP_N,1:nElems) )
    rf_DG = 0.0
    rg_DG = 0.0
    rh_DG = 0.0
#endif /*LOCAL_ALPHA*/
    
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
    use MOD_IDP_Vars    , only: IDPSpecEntropy, IDPMathEntropy, IDPDensityTVD, IDPSemiDiscEnt, IDPPositivity
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
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
    call Get_IDP_Variables(U,dt)
    
!   Perform limiting!
!   -----------------
    if (IDPSemiDiscEnt) call IDP_LimitSemiDiscEnt(U,Ut,dt,sdt)
    if (IDPDensityTVD)  call IDP_LimitDensityTVD (U,Ut,dt,sdt)
    if (IDPSpecEntropy) call IDP_LimitSpecEntropy(U,Ut,dt,sdt)
    if (IDPMathEntropy) call IDP_LimitMathEntropy(U,Ut,dt,sdt)
    if (IDPPositivity)  call IDP_LimitPositivity (U,Ut,dt,sdt)
    
!   Update variables for the analyze routines
!   -----------------------------------------
    maximum_alpha = max(maximum_alpha,maxval(alpha-alpha_old))
    
    amount_alpha = amount_alpha*amount_alpha_steps
    ! TODO: how to do the local alpha amount??... This is not working:
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
  subroutine IDP_LimitDensityTVD(U,Ut,dt,sdt)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
    use MOD_Mesh_Vars     , only: nElems
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
    use MOD_NFVSE_Vars    , only: ftilde_DG, gtilde_DG, htilde_DG, ftilde_FV, gtilde_FV, htilde_FV
    use MOD_NFVSE_Vars    , only: sWGP
    use MOD_Mesh_Vars     , only: sJ
    use MOD_IDP_Vars      , only: Usafe
#endif /*LOCAL_ALPHA*/
#if barStates
    use MOD_IDP_Vars      , only: Ubar_xi, Ubar_eta, Ubar_zeta, Uprev
#else
    use MOD_IDP_Vars      , only: Usafe
#endif /*barStates*/
    use MOD_IDP_Vars      , only: FFV_m_FDG
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    !-local-variables----------------------------------------
    real    :: corr, corr1
#if LOCAL_ALPHA
    real    :: corr_loc     (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1)
#endif /*LOCAL_ALPHA*/
    real    :: rho_min, rho_max, rho_safe
    real    :: a   ! a  = PositCorrFactor * rho_safe - rho
    real    :: Qp, Qm, Pp, Pm
    integer :: eID
    integer :: i,j,k,l
    real, parameter :: eps = 1.e-10           ! Very small value
    !--------------------------------------------------------
    
    do eID=1, nElems
      
      corr = -epsilon(1.0) ! Safe initialization
#if LOCAL_ALPHA
      corr_loc = 0.0
#endif /*LOCAL_ALPHA*/
!       Compute correction factors
!       --------------------------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
          
        ! Get the limit states
        !*********************
        rho_min =  huge(1.0)
        rho_max = -huge(1.0)
        
#if barStates
        ! Previous sol
        rho_min = min(rho_min, Uprev  (1,i  ,j  ,k  ,eID))
        rho_max = max(rho_max, Uprev  (1,i  ,j  ,k  ,eID))
        
        !xi
        do l=i-1, i
          rho_min = min(rho_min, Ubar_xi  (1,l  ,j  ,k  ,eID))
          rho_max = max(rho_max, Ubar_xi  (1,l  ,j  ,k  ,eID))
        end do
        !eta
        do l=j-1, j
          rho_min = min(rho_min, Ubar_eta (1,i  ,l  ,k  ,eID))
          rho_max = max(rho_max, Ubar_eta (1,i  ,l  ,k  ,eID))
        end do
        !zeta
        do l=k-1, k
          rho_min = min(rho_min, Ubar_zeta(1,i  ,j  ,l  ,eID))
          rho_max = max(rho_max, Ubar_zeta(1,i  ,j  ,l  ,eID))
        end do
#else
        ! check stencil in xi
!          do l = i+idx_m1(i), i+idx_p1(i) !no neighbor
        do l = i-1, i+1
          rho_min = min(rho_min, Usafe(1,l,j,k,eID))
          rho_max = max(rho_max, Usafe(1,l,j,k,eID))
        end do
        ! check stencil in eta
!          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          rho_min = min(rho_min, Usafe(1,i,l,k,eID))
          rho_max = max(rho_max, Usafe(1,i,l,k,eID))
        end do
        ! check stencil in zeta
!          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          rho_min = min(rho_min, Usafe(1,i,j,l,eID))
          rho_max = max(rho_max, Usafe(1,i,j,l,eID))
        end do
#endif /*barStates*/
        
#if LOCAL_ALPHA
        ! Real Zalesak type limiter
        ! * Zalesak (1979). "Fully multidimensional flux-corrected transport algorithms for fluids"
        ! * Kuzmin et al. (2010). "Failsafe flux limiting and constrained data projections for equations of gas dynamics"
        !****************************************************************************************************************
        if ( (U(1,i,j,k,eID) >= rho_min) .and. (U(1,i,j,k,eID) <= rho_max) ) cycle
        
        ! Upper/lower bounds for admissible increments
        Qp = max(0.0,(rho_max-Usafe(1,i,j,k,eID))*sdt)
        Qm = min(0.0,(rho_min-Usafe(1,i,j,k,eID))*sdt)
        
        ! Check the sign of Qp and Qm... !!!
        if (Qp < 0) then
          print*, '0>Qp=', Qp, i, j, k, eID
          print*, 'rho_max', rho_max
          print*, 'rho_saf', Usafe(1,i,j,k,eID)
          read(*,*)
!          stop 
        end if
        if (Qm > 0)  then
          print*, '0<Qm=', Qm, i, j, k, eID
          print*, 'rho_min', rho_min
          print*, 'rho_saf', Usafe(1,i,j,k,eID)
          read(*,*)
!          stop
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
        
        if (Pp==0) then
          Qp = 1.0
        else
          Qp = Qp/Pp
        end if
        
        if (Pm==0) then
          Qm = 1.0
        else
          Qm = Qm/Pm
        end if
        
        corr1 = 1.0 - min(1.0,Qp,Qm)
        
        corr_loc(i,j,k) = max(corr_loc(i,j,k),corr1)
        corr = max(corr,corr1)
        
#else
        ! Naive limiter that gets out of bounds
        !********************************************
        if ( U(1,i,j,k,eID) < rho_min) then
          rho_safe = rho_min
        elseif (U(1,i,j,k,eID) > rho_max) then
          rho_safe = rho_max
        else
          cycle
        end if
        
        if ( abs(FFV_m_FDG(1,i,j,k,eID)) < eps) cycle !nothing to do here!
        
        ! Density correction
        a = (rho_safe - U(1,i,j,k,eID))
        corr1 = a / FFV_m_FDG(1,i,j,k,eID)
#if LOCAL_ALPHA
        corr_loc(i,j,k) = max(corr_loc(i,j,k),corr1)
#endif /*LOCAL_ALPHA*/
        corr = max(corr,corr1)
        
#endif /*LOCAL_ALPHA*/
        
      end do       ; end do       ; end do ! i,j,k
      
!       Do the correction if needed
!       ---------------------------
      
      if ( corr > 0. ) then
        call PerformCorrection(U(:,:,:,:,eID),Ut(:,:,:,:,eID),corr    ,alpha(eID)          , &
#if LOCAL_ALPHA
                                                              corr_loc,alpha_loc(:,:,:,eID), &
#endif /*LOCAL_ALPHA*/
                                                              dt,sdt,eID)
      end if
      
    end do !eID
    
  end subroutine IDP_LimitDensityTVD
!===================================================================================================================================
!> Semi-discrete entropy correction
!> ATTENTION: 1) Needs debugging
!>            2) Does not work with MPI
!>            3) TODO: Check if the subcell limiting is valid
!>            4) Only works for local alpha and disc2
!===================================================================================================================================
  subroutine IDP_LimitSemiDiscEnt(U,Ut,dt,sdt)
    use MOD_PreProc       , only: PP_N
    use MOD_Mesh_Vars     , only: nElems
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_Equation_Vars , only: ConsToEntropy, GetEntropyPot
    use MOD_IDP_Vars      , only: Flux_ext, Uprev
    use MOD_Globals       , only: nProcessors
    use MOD_NFVSE_MPI     , only: Get_externalU
    use MOD_DG_Vars       , only: Flux_master, Flux_slave
    use MOD_NFVSE_Vars    , only: SubCellMetrics
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: rf_DG, rg_DG, rh_DG, ftilde_FV, gtilde_FV, htilde_FV
#endif /*LOCAL_ALPHA*/
#endif /*LOCAL_ALPHA*/
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    !-local-variables----------------------------------------
    real    :: corr, corr1
#if LOCAL_ALPHA
    real    :: corr_loc     (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1)
#endif /*LOCAL_ALPHA*/
    integer :: eID
    integer :: i,j,k
    
    real :: Ent_Jump(PP_nVar), Psi_Jump, r, rFV
    real :: entVar(PP_nVar,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1)
    real :: entPot(3      ,-1:PP_N+1,-1:PP_N+1,-1:PP_N+1)
    
    real :: rf_FV(-1:PP_N, 0:PP_N, 0:PP_N)
    real :: rg_FV( 0:PP_N,-1:PP_N, 0:PP_N)
    real :: rh_FV( 0:PP_N, 0:PP_N,-1:PP_N)
    
    real, parameter :: eps = 1.e-10           ! Very small value
    !--------------------------------------------------------
    
!   Some definitions
!   ****************
    rf_FV = 0.0
    rg_FV = 0.0
    rh_FV = 0.0
    
#if LOCAL_ALPHA
    ! ATTENTION: 1) we are supposing that alpha=0 before limiting
    
    ! Get the boundary entropy production (same for FV and DG)
    ! ********************************************************
    
    
    ! TODO: Send the F_master across MPI interfaces
#if MPI
    if (nProcessors > 1) then
      stop 'IDPSemiDiscEnt not working in parallel.. Send Flux_master!!'
    end if
#endif /*MPI*/
    
    ! Get the external fluxes in place (only working without BCs!!! TODO..)
    call Get_externalU(PP_nVar,Flux_ext,Uprev(:,0:PP_N,0:PP_N,0:PP_N,:),Flux_master,Flux_slave,fluxSign=.TRUE.)
    
    do eID=1, nElems
      
      ! Get entropy vars and potential
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        entVar(:,i,j,k) = ConsToEntropy(Uprev(:,i,j,k,eID))
        entPot(:,i,j,k) = GetEntropyPot(Uprev(:,i,j,k,eID),entVar(:,i,j,k))
      end do       ; end do       ; end do !i,j,k
      
      ! Get entropy vars and potential (boundaries)
      do k=0, PP_N ; do j=0, PP_N
        !xi-
        entVar(:,-1,j,k) = ConsToEntropy(Uprev(:,-1,j,k,eID))
        entPot(:,-1,j,k) = GetEntropyPot(Uprev(:,-1,j,k,eID),entVar(:,-1,j,k))
        !xi+
        entVar(:,PP_N+1,j,k) = ConsToEntropy(Uprev(:,PP_N+1,j,k,eID))
        entPot(:,PP_N+1,j,k) = GetEntropyPot(Uprev(:,PP_N+1,j,k,eID),entVar(:,PP_N+1,j,k))
        !eta-
        entVar(:,j,-1,k) = ConsToEntropy(Uprev(:,j,-1,k,eID))
        entPot(:,j,-1,k) = GetEntropyPot(Uprev(:,j,-1,k,eID),entVar(:,j,-1,k))
        !eta+
        entVar(:,j,PP_N+1,k) = ConsToEntropy(Uprev(:,j,PP_N+1,k,eID))
        entPot(:,j,PP_N+1,k) = GetEntropyPot(Uprev(:,j,PP_N+1,k,eID),entVar(:,j,PP_N+1,k))
        !zeta-
        entVar(:,j,k,-1) = ConsToEntropy(Uprev(:,j,k,-1,eID))
        entPot(:,j,k,-1) = GetEntropyPot(Uprev(:,j,k,-1,eID),entVar(:,j,k,-1))
        !zeta+
        entVar(:,j,k,PP_N+1) = ConsToEntropy(Uprev(:,j,k,PP_N+1,eID))
        entPot(:,j,k,PP_N+1) = GetEntropyPot(Uprev(:,j,k,PP_N+1,eID),entVar(:,j,k,PP_N+1))
      end do          ; end do !j,k
      
      ! Compute boundary entropy production
      ! -----------------------------------
      do k=0, PP_N  ; do j=0, PP_N
        !xi- (minus flux)
        Ent_Jump = entVar(:,0,j,k) - entVar(:,-1,j,k)
        Psi_Jump = dot_product(SubCellMetrics(eID) % xi % nv(:,j,k,-1  )*SubCellMetrics(eID) % xi % norm(j,k,-1  ), entPot(:,0     ,j,k) - entPot(:,-1  ,j,k) )
        rf_DG(-1,j,k,eID) = -dot_product(Ent_Jump,Flux_ext(:,j,k,5,eID)) - Psi_Jump
        rf_FV(-1,j,k)     = rf_DG(-1,j,k,eID)
        !xi+
        Ent_Jump = entVar(:,PP_N+1,j,k) - entVar(:,PP_N,j,k)
        Psi_Jump = dot_product(SubCellMetrics(eID) % xi % nv(:,j,k,PP_N)*SubCellMetrics(eID) % xi % norm(j,k,PP_N), entPot(:,PP_N+1,j,k) - entPot(:,PP_N,j,k) )
        rf_DG(PP_N,j,k,eID) = dot_product(Ent_Jump,Flux_ext(:,j,k,3,eID)) - Psi_Jump
        rf_FV(PP_N,j,k)     = rf_DG(PP_N,j,k,eID)
        !eta- (minus flux)
        Ent_Jump = entVar(:,j,0,k) - entVar(:,j,-1,k)
        Psi_Jump = dot_product(SubCellMetrics(eID) % eta % nv(:,j,k,-1  )*SubCellMetrics(eID) % eta % norm(j,k,-1  ), entPot(:,j,0     ,k) - entPot(:,j,-1  ,k) )
        rg_DG(j,-1,k,eID) = -dot_product(Ent_Jump,Flux_ext(:,j,k,2,eID)) - Psi_Jump
        rg_FV(j,-1,k)     = rg_DG(j,-1,k,eID)
        !eta+
        Ent_Jump = entVar(:,j,PP_N+1,k) - entVar(:,j,PP_N,k)
        Psi_Jump = dot_product(SubCellMetrics(eID) % eta % nv(:,j,k,PP_N)*SubCellMetrics(eID) % eta % norm(j,k,PP_N), entPot(:,j,PP_N+1,k) - entPot(:,j,PP_N,k) )
        rg_DG(j,PP_N,k,eID) = dot_product(Ent_Jump,Flux_ext(:,j,k,4,eID)) - Psi_Jump
        rg_FV(j,PP_N,k)     = rg_DG(j,PP_N,k,eID)
        !zeta- (minus flux)
        Ent_Jump = entVar(:,j,k,0) - entVar(:,j,k,-1)
        Psi_Jump = dot_product(SubCellMetrics(eID) % zeta % nv(:,j,k,-1  )*SubCellMetrics(eID) % zeta % norm(j,k,-1  ), entPot(:,j,k,0     ) - entPot(:,j,k,-1  ) )
        rh_DG(j,k,-1,eID) = -dot_product(Ent_Jump,Flux_ext(:,j,k,1,eID)) - Psi_Jump
        rh_FV(j,k,-1)     = rh_DG(j,k,-1,eID)
        !zeta+
        Ent_Jump = entVar(:,j,k,PP_N+1) - entVar(:,j,k,PP_N)
        Psi_Jump = dot_product(SubCellMetrics(eID) % zeta % nv(:,j,k,PP_N)*SubCellMetrics(eID) % zeta % norm(j,k,PP_N), entPot(:,j,k,PP_N+1) - entPot(:,j,k,PP_N) )
        rh_DG(j,k,PP_N,eID) = dot_product(Ent_Jump,Flux_ext(:,j,k,6,eID)) - Psi_Jump
        rh_FV(j,k,PP_N)     = rh_DG(j,k,PP_N,eID)
      end do        ; end do
      
      ! compute entropy production of FV
      ! --------------------------------
      !xi
      do i=0, PP_N-1
        do k=0, PP_N ; do j=0, PP_N
          ! Jumps
          Ent_Jump = entVar(:,i+1,j,k) - entVar(:,i,j,k)
          Psi_Jump = dot_product(SubCellMetrics(eID) % xi % nv(:,j,k,i)*SubCellMetrics(eID) % xi % norm(j,k,i), entPot(:,i+1,j,k) - entPot(:,i,j,k) )
          ! Compute entropy production
          rf_FV(i,j,k) = dot_product(Ent_Jump, ftilde_FV(:,i,j,k,eID)) - Psi_Jump
        end do       ; end do !j,k
      end do !i
      !Eta
      do j=0, PP_N-1
        do k=0, PP_N ; do i=0, PP_N
          ! Jumps
          Ent_Jump = entVar(:,i,j+1,k) - entVar(:,i,j,k)
          Psi_Jump = dot_product(SubCellMetrics(eID) % eta % nv(:,i,k,j)*SubCellMetrics(eID) % eta % norm(i,k,j), entPot(:,i,j+1,k) - entPot(:,i,j,k) )
          ! Compute entropy production
          rg_FV(i,j,k) = dot_product(Ent_Jump, gtilde_FV(:,i,j,k,eID)) - Psi_Jump
        end do       ; end do !j,k
      end do !i
      !Zeta
      do k=0, PP_N-1
        do j=0, PP_N ; do i=0, PP_N
          ! Jumps
          Ent_Jump = entVar(:,i,j,k+1) - entVar(:,i,j,k)
          Psi_Jump = dot_product(SubCellMetrics(eID) % zeta % nv(:,i,j,k)*SubCellMetrics(eID) % zeta % norm(i,j,k), entPot(:,i,j,k+1) - entPot(:,i,j,k) )
          ! Compute entropy production
          rh_FV(i,j,k) = dot_product(Ent_Jump, htilde_FV(:,i,j,k,eID)) - Psi_Jump
        end do       ; end do !j,k
      end do !i
      
      ! Get corr factors!!!!
      corr = -epsilon(1.0) ! Safe initialization
#if LOCAL_ALPHA
      corr_loc = 0.0
#endif /*LOCAL_ALPHA*/

      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Get the total DG production
        !****************************
        r = rf_DG(i,j,k,eID) + rf_DG(i-1,j,k,eID) + &
            rg_DG(i,j,k,eID) + rg_DG(i,j-1,k,eID) + &
            rh_DG(i,j,k,eID) + rh_DG(i,j,k-1,eID) 
        
        if (r < eps) cycle
        
        rFV = rf_FV(i,j,k) + rf_FV(i-1,j,k) + &
              rg_FV(i,j,k) + rg_FV(i,j-1,k) + &
              rh_FV(i,j,k) + rh_FV(i,j,k-1) 
        
        if (rFV > eps) then
          print*, 'WARNING: Low-order method is not entropy dissipative!!', eID, i,j,k
          stop
        end if
        
        rFV = r - rFV
        if (abs(rFV) < eps) cycle !nothing to do here
        
        ! correction
        corr1 = r*dt/rFV
#if LOCAL_ALPHA
        corr_loc(i,j,k) = max(corr_loc(i,j,k),corr1)
#endif /*LOCAL_ALPHA*/
        corr = max(corr,corr1)
        
      end do       ; end do       ; end do ! i,j,k
      
!       Do the correction if needed
!       ---------------------------
      
      if ( corr > 0. ) then
        ! TODO: Check if this holds at the ubcell level
        call PerformCorrection(U(:,:,:,:,eID),Ut(:,:,:,:,eID),corr    ,alpha(eID)          , &
#if LOCAL_ALPHA
                                                              corr_loc,alpha_loc(:,:,:,eID), &
#endif /*LOCAL_ALPHA*/
                                                              dt,sdt,eID)
      end if
      
    end do !nElems
#else
    stop 'IDPSemiDiscEnt only for local alpha'
#endif /*LOCAL_ALPHA*/
    
  end subroutine IDP_LimitSemiDiscEnt
!===================================================================================================================================
!> Specific entropy correction (discrete local minimum principle)
!===================================================================================================================================
  subroutine IDP_LimitSpecEntropy(U,Ut,dt,sdt)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
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
    use MOD_IDP_Vars      , only: alpha_maxIDP, IDPMaxIter
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    !-local-variables----------------------------------------
    real    :: corr, corr1, corr_old
#if LOCAL_ALPHA
    real    :: corr_loc     (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1)
#endif /*LOCAL_ALPHA*/
    real    :: s_min, as, U_curr(PP_nVar), dSdalpha
    integer :: eID
    integer :: i,j,k, l
    integer :: iter
    logical :: notInIter
    real, parameter :: eps = 1.e-10           ! Very small value
    !--------------------------------------------------------
    
    do eID=1, nElems
      
      corr = -epsilon(1.0) ! Safe initialization
#if LOCAL_ALPHA
      corr_loc = 0.0
#endif /*LOCAL_ALPHA*/
!       Compute correction factors
!       --------------------------
      notInIter = .FALSE.
      
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Get the limit states
        !*********************
        s_min = huge(1.0)
        
#if barStates
        ! Previous entropy of the node (ubar_ii)
        s_min = min(s_min, Get_SpecEntropy(Uprev    (:,i  ,j  ,k  ,eID)))
        
        ! TODO: Compute them for all interfaces before...
        !xi+
        s_min = min(s_min, Get_SpecEntropy(Ubar_xi  (:,i  ,j  ,k  ,eID)))
        !xi-
        s_min = min(s_min, Get_SpecEntropy(Ubar_xi  (:,i-1,j  ,k  ,eID)))
        !eta+
        s_min = min(s_min, Get_SpecEntropy(Ubar_eta (:,i  ,j  ,k  ,eID)))
        !eta-
        s_min = min(s_min, Get_SpecEntropy(Ubar_eta (:,i  ,j-1,k  ,eID)))
        !zeta+
        s_min = min(s_min, Get_SpecEntropy(Ubar_zeta(:,i  ,j  ,k  ,eID)))
        !zeta-
        s_min = min(s_min, Get_SpecEntropy(Ubar_zeta(:,i  ,j  ,k-1,eID)))
#else
        ! check stencil in xi
!#          do l = i+idx_m1(i), i+idx_p1(i) !no neighbor
        do l = i-1, i+1
          s_min = min(s_min, EntPrev(1,l,j,k,eID))
        end do
        ! check stencil in eta
!#          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          s_min = min(s_min, EntPrev(1,i,l,k,eID))
        end do
        ! check stencil in zeta
!#          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          s_min = min(s_min, EntPrev(1,i,j,l,eID))
        end do
#endif /*barStates*/
        
        ! Difference between goal entropy and current entropy
        as = (s_min - Get_SpecEntropy(U(:,i,j,k,eID)))

        if (as <= max(eps,abs(s_min)*NEWTON_ABSTOL)) cycle ! this DOF does NOT need entropy correction   
      
        ! Newton initialization:
        U_curr = U(:,i,j,k,eID)
        corr1 = 0.0
        
        ! Perform Newton iterations
        NewtonLoopSpec: do iter=1, IDPMaxIter
          corr_old = corr1
          
          ! Evaluate dS/d(alpha)
          dSdalpha = dot_product(ConsToSpecEntropy(U_curr),FFV_m_FDG(:,i,j,k,eID))
          
          if ( abs(dSdalpha)<eps) exit NewtonLoopSpec ! Nothing to do here!
          
          ! Update correction
          corr1 = corr1 + as / dSdalpha
          if (alpha(eID) + corr1 * sdt > alpha_maxIDP) then
            corr1 = (alpha_maxIDP - alpha(eID) ) * dt
          end if
          if (corr1<0) corr1=0.0
          
          ! Check relative tolerance
          if ( abs(corr_old-corr1)<= NEWTON_RELTOL ) exit NewtonLoopSpec
          
          ! Get new U
          U_curr = U (:,i,j,k,eID) + corr1 * FFV_m_FDG(:,i,j,k,eID)
          
          ! Evaluate if goal entropy was achieved (and exit the Newton loop if that's the case)
          as = s_min-Get_SpecEntropy(U_curr)
          
          ! Check absolute tolerance
          if ( abs(as) < max(eps,abs(s_min)*NEWTON_ABSTOL) ) exit NewtonLoopSpec  
          
        end do NewtonLoopSpec ! iter
        
        if (iter > IDPMaxIter) notInIter =.TRUE.
        
#if LOCAL_ALPHA
        corr_loc(i,j,k) = max(corr_loc(i,j,k),corr1)
#endif /*LOCAL_ALPHA*/
        corr = max(corr,corr1) ! Compute the element-wise maximum correction
      end do       ; end do       ; enddo !i,j,k
      
      if (notInIter) then
        write(*,'(A,I0,A,I0)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem
      end if
      
!       Do the correction if needed
!       ---------------------------
      if ( corr > 0. ) then
        call PerformCorrection(U(:,:,:,:,eID),Ut(:,:,:,:,eID),corr    ,alpha(eID)          , &
#if LOCAL_ALPHA
                                                              corr_loc,alpha_loc(:,:,:,eID), &
#endif /*LOCAL_ALPHA*/
                                                              dt,sdt,eID)
      end if
    
    end do !eID
    
  end subroutine IDP_LimitSpecEntropy  
!===================================================================================================================================
!> Mathematical entropy correction (discrete local maximum principle)
!===================================================================================================================================
  subroutine IDP_LimitMathEntropy(U,Ut,dt,sdt)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
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
    use MOD_IDP_Vars      , only: alpha_maxIDP, IDPMaxIter
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    !-local-variables----------------------------------------
    real    :: corr, corr1, corr_old
#if LOCAL_ALPHA
    real    :: corr_loc     (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1)
#endif /*LOCAL_ALPHA*/
    real    :: s_max, as, U_curr(PP_nVar), dSdalpha
    integer :: eID
    integer :: i,j,k, l
    integer :: iter
    logical :: notInIter
    real, parameter :: eps = 1.e-10           ! Very small value
    !--------------------------------------------------------
    
    do eID=1, nElems
      
      corr = -epsilon(1.0) ! Safe initialization
#if LOCAL_ALPHA
      corr_loc = 0.0
#endif /*LOCAL_ALPHA*/
!       Compute correction factors
!       --------------------------
      notInIter = .FALSE.
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Get the limit states
        !*********************
        s_max = -huge(1.0)
        
#if barStates
        ! Previous entropy of the node (ubar_ii)
        s_max = max(s_max, Get_MathEntropy(Uprev    (:,i  ,j  ,k  ,eID)))
        
        ! TODO: Compute for all interfaces before the loop!
        !xi+
        s_max = max(s_max, Get_MathEntropy(Ubar_xi  (:,i  ,j  ,k  ,eID)))
        !xi-
        s_max = max(s_max, Get_MathEntropy(Ubar_xi  (:,i-1,j  ,k  ,eID)))
        !eta+
        s_max = max(s_max, Get_MathEntropy(Ubar_eta (:,i  ,j  ,k  ,eID)))
        !eta-
        s_max = max(s_max, Get_MathEntropy(Ubar_eta (:,i  ,j-1,k  ,eID)))
        !zeta+
        s_max = max(s_max, Get_MathEntropy(Ubar_zeta(:,i  ,j  ,k  ,eID)))
        !zeta-
        s_max = max(s_max, Get_MathEntropy(Ubar_zeta(:,i  ,j  ,k-1,eID)))
#else
        ! check stencil in xi
!          do l = i+idx_m1(i), i+idx_p1(i) !no neighbor
        do l = i-1, i+1
          s_max = max(s_max, EntPrev(1,l,j,k,eID))
        end do
        ! check stencil in eta
!          do l = j+idx_m1(j), j+idx_p1(j) !no neighbor
        do l = j-1, j+1
          s_max = max(s_max, EntPrev(1,i,l,k,eID))
        end do
        ! check stencil in zeta
!          do l = k+idx_m1(k), k+idx_p1(k) !no neighbor
        do l = k-1, k+1
          s_max = max(s_max, EntPrev(1,i,j,l,eID))
        end do
#endif /*barStates*/
        
        ! Difference between goal entropy and current entropy
        as = (s_max - Get_MathEntropy(U(:,i,j,k,eID)))

        if (as >= -max(eps,abs(s_max)*NEWTON_ABSTOL)) cycle ! this DOF does NOT need pressure correction   
      
        ! Newton initialization:
        U_curr = U(:,i,j,k,eID)
        corr1 = 0.0
        
        ! Perform Newton iterations
        NewtonLoop: do iter=1, IDPMaxIter
          corr_old = corr1
          
          ! Evaluate dS/d(alpha)
          dSdalpha = dot_product(ConsToEntropy(U_curr),FFV_m_FDG(:,i,j,k,eID))
          
          if ( abs(dSdalpha)<eps ) exit NewtonLoop ! Nothing to do here!
          
          ! Update correction
          corr1 = corr1 + as / dSdalpha
          if (alpha(eID) + corr1 * sdt > alpha_maxIDP) then
            corr1 = (alpha_maxIDP - alpha(eID) ) * dt
          end if
          if (corr1<0) corr1=0.0
          
          ! Check relative tolerance
          if ( abs(corr_old-corr1)<= NEWTON_RELTOL ) exit NewtonLoop
          
          ! Get new U
          U_curr = U (:,i,j,k,eID) + corr1 * FFV_m_FDG(:,i,j,k,eID)
          
          ! Evaluate if goal entropy was achieved (and exit the Newton loop if that's the case)
          as = s_max-Get_MathEntropy(U_curr)
          
          ! Check absolute tolerance
          if ( abs(as) < max(eps,abs(s_max)*NEWTON_ABSTOL) ) exit NewtonLoop  
          
        end do NewtonLoop ! iter
        
        if (iter > IDPMaxIter) notInIter =.TRUE.
        
#if LOCAL_ALPHA
        corr_loc(i,j,k) = max(corr_loc(i,j,k),corr1)
#endif /*LOCAL_ALPHA*/
        corr = max(corr,corr1) ! Compute the element-wise maximum correction
      end do       ; end do       ; enddo !i,j,k
      
      if (notInIter) then
        write(*,'(A,I0,A,I0)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem
      end if
      
!       Do the correction if needed
!       ---------------------------
      if ( corr > 0. ) then
        call PerformCorrection(U(:,:,:,:,eID),Ut(:,:,:,:,eID),corr    ,alpha(eID)          , &
#if LOCAL_ALPHA
                                                              corr_loc,alpha_loc(:,:,:,eID), &
#endif /*LOCAL_ALPHA*/
                                                              dt,sdt,eID)
      end if
    
    end do !eID
    
  end subroutine IDP_LimitMathEntropy
!===================================================================================================================================
!> Density and pressure positivity limiter. As explained in:
!> * Rueda-Ramírez, A. M., & Gassner, G. J. (2021). A Subcell Finite Volume Positivity-Preserving Limiter for DGSEM Discretizations of the Euler Equations. arXiv preprint arXiv:2102.06017.
!>    ... But with the possibility to be used locally
!===================================================================================================================================
  subroutine IDP_LimitPositivity(U,Ut,dt,sdt)
    use MOD_PreProc       , only: PP_N
    use MOD_NFVSE_Vars    , only: alpha, PositCorrFactor
    use MOD_Mesh_Vars     , only: nElems, offsetElem
#if LOCAL_ALPHA
    use MOD_NFVSE_Vars    , only: alpha_loc
#endif /*LOCAL_ALPHA*/
    use MOD_IDP_Vars      , only: Usafe, p_safe
    use MOD_IDP_Vars      , only: FFV_m_FDG
    use MOD_IDP_Vars      , only: alpha_maxIDP, IDPMaxIter, NEWTON_ABSTOL
    use MOD_Equation_Vars , only: Get_Pressure, Get_dpdU
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(inout) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current Ut (in RK stage)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    real,intent(in)    :: sdt                                       !< Inverse of current RK time-step size (in RK stage)
    !-local-variables----------------------------------------
    real    :: corr, corr1
#if LOCAL_ALPHA
    real    :: corr_loc     (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1)
#endif /*LOCAL_ALPHA*/
    real    :: a   ! a  = PositCorrFactor * rho_safe - rho
    real    :: ap  ! ap = (PositCorrFactor * p_safe   - p) / (kappa-1)
    real    :: pres
    real    :: alphadiff
    real    :: dpdU(PP_nVar), U_curr(PP_nVar), p_goal
    real    :: dp_dalpha
    integer :: eID
    integer :: i,j,k, iter
    logical :: NotInIter
    real, parameter :: eps = 1.e-14           ! Very small value
    !--------------------------------------------------------
    
    do eID=1, nElems
      
!     ----------------------------------
!     Check if it makes sense correcting
!     ----------------------------------
      alphadiff = alpha(eID) - alpha_maxIDP
      if ( abs(alphadiff) < eps ) cycle ! Not much to do for this element...
      
!     ---------------
!     Correct density
!     ---------------
      
      corr = -epsilon(1.0) ! Safe initialization
#if LOCAL_ALPHA
      corr_loc = 0.0
#endif /*LOCAL_ALPHA*/
        
!     Compute correction factors
!     --------------------------
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        
        ! Density correction
        a = (PositCorrFactor * Usafe(1,i,j,k,eID) - U(1,i,j,k,eID))
        if (a > 0.) then ! This DOF needs a correction
          if (abs(FFV_m_FDG(1,i,j,k,eID)) < eps) cycle
          corr1 = a / FFV_m_FDG(1,i,j,k,eID)
#if LOCAL_ALPHA
          corr_loc(i,j,k) = max(corr_loc(i,j,k),corr1)
#endif /*LOCAL_ALPHA*/
          corr = max(corr,corr1)
        end if
        
      end do       ; end do       ; end do ! i,j,k
        
      
!       Do the correction if needed
!       ---------------------------
      if ( corr > 0. ) then
        call PerformCorrection(U(:,:,:,:,eID),Ut(:,:,:,:,eID),corr    ,alpha(eID)          , &
#if LOCAL_ALPHA
                                                              corr_loc,alpha_loc(:,:,:,eID), &
#endif /*LOCAL_ALPHA*/
                                                              dt,sdt,eID)
      end if
      
!     ---------------
!     Correct pressure
!     ---------------
      
      corr = -epsilon(1.0) ! Safe initialization
#if LOCAL_ALPHA
      corr_loc = 0.0
#endif /*LOCAL_ALPHA*/
      
!     Compute correction factors
!     --------------------------
      notInIter = .FALSE.
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        ! Current pressure and goal
        call Get_Pressure(U(:,i,j,k,eID),pres)
        p_goal = PositCorrFactor * p_safe(i,j,k,eID)
        ap = (p_goal - pres)
        
        if (ap <= 0.) cycle ! this DOF does NOT need pressure correction
        
        ! Newton initialization:
        U_curr = U(:,i,j,k,eID)
        corr1 = 0.0
        
        ! Perform Newton iterations
        NewtonLoop: do iter=1, IDPMaxIter
          ! Evaluate dp/d(alpha)
          call Get_dpdU(U_curr,dpdU)
          dp_dalpha = dot_product(dpdU,FFV_m_FDG(:,i,j,k,eID))
          if ( abs(dp_dalpha)<eps ) exit NewtonLoop ! Nothing to do here!
          
          ! Update correction
          corr1 = corr1 + ap / dp_dalpha
          
          ! Get new U and pressure
          U_curr = U (:,i,j,k,eID) + corr1 * FFV_m_FDG(:,i,j,k,eID)
          call Get_Pressure(U_curr,pres)
          
          ! Evaluate if goal pressure was achieved (and exit the Newton loop if that's the case)
          ap = p_goal-pres
          if ( (ap <= epsilon(p_goal)) .and. (ap > -NEWTON_ABSTOL*p_goal) ) exit NewtonLoop  ! Note that we use an asymmetric tolerance!
        end do NewtonLoop ! iter
        
        if (iter > IDPMaxIter) notInIter =.TRUE.
        
#if LOCAL_ALPHA
        corr_loc(i,j,k) = max(corr_loc(i,j,k),corr1)
#endif /*LOCAL_ALPHA*/
        corr = max(corr,corr1) ! Compute the element-wise maximum correction
      
      end do       ; end do       ; enddo !i,j,k
      
      if (notInIter) then
        write(*,'(A,I0,A,I0)') 'WARNING: Not able to perform NFVSE correction within ', IDPMaxIter, ' Newton iterations. Elem: ', eID + offsetElem
      end if
      
!       Do the correction if needed
!       ---------------------------
      if ( corr > 0. ) then
        call PerformCorrection(U(:,:,:,:,eID),Ut(:,:,:,:,eID),corr    ,alpha(eID)          , &
#if LOCAL_ALPHA
                                                              corr_loc,alpha_loc(:,:,:,eID), &
#endif /*LOCAL_ALPHA*/
                                                              dt,sdt,eID)
      end if
      
    end do !eID
    
  end subroutine IDP_LimitPositivity
!===================================================================================================================================
!> Takes corr/corr_loc U and Ut, and outputs the corrected U and Ut, and alpha/alpha_loc for visualization
!===================================================================================================================================
  pure subroutine PerformCorrection(U,Ut,corr    ,alpha    , &
#if LOCAL_ALPHA
                                         corr_loc,alpha_loc, &
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
    implicit none
    !-arguments--------------------------------------
    real, intent(inout) :: U (PP_nVar, 0:PP_N  , 0:PP_N  , 0:PP_N)
    real, intent(inout) :: Ut(PP_nVar, 0:PP_N  , 0:PP_N  , 0:PP_N)
    real, intent(inout) :: corr                             ! Scaled element-wise d_alpha
    real, intent(inout) :: alpha                            ! Current element-wise alpha
#if LOCAL_ALPHA
    real, intent(inout) :: corr_loc  (-1:PP_N+1,-1:PP_N+1,-1:PP_N+1) ! Scaled local d_alpha
    real, intent(inout) :: alpha_loc ( 0:PP_N  , 0:PP_N  , 0:PP_N) ! Current local  alpha
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
    alpha      = alpha + corr * sdt
    
    ! Change inconsistent alphas
    if (alpha > alpha_maxIDP) then
      alpha = alpha_maxIDP
      corr  = (alpha_maxIDP - alphacont ) * dt
    end if
          
#if LOCAL_ALPHA
    ! Change the alpha for output
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      alphacont_loc    = alpha_loc(i,j,k)
      alpha_loc(i,j,k) = alpha_loc(i,j,k) + corr_loc(i,j,k)* sdt
      if (alpha_loc(i,j,k) > alpha_maxIDP) then
        alpha_loc(i,j,k) = alpha_maxIDP
        corr_loc (i,j,k) = (alpha_maxIDP - alphacont_loc ) * dt
      end if
    end do       ; end do       ; enddo
    
    ! Correct!
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      ! xi correction
      ! -------------
      ! left
      my_corr=-max(corr_loc(i-1,j  ,k  ),corr_loc(i  ,j  ,k  )) * sWGP(i) * (ftilde_DG(:,i-1,j,k,eID) - ftilde_FV(:,i-1,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr 
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr * sdt
      
      ! right
      my_corr=max(corr_loc(i  ,j  ,k  ),corr_loc(i+1,j  ,k  )) * sWGP(i) * (ftilde_DG(:,i  ,j,k,eID) - ftilde_FV(:,i  ,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr 
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr * sdt
      
      ! eta correction
      ! --------------
      ! left
      my_corr=-max(corr_loc(i  ,j-1,k  ),corr_loc(i  ,j  ,k  )) * sWGP(j) * (gtilde_DG(:,i,j-1,k,eID) - gtilde_FV(:,i,j-1,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr * sdt
      
      ! right
      my_corr=max(corr_loc(i  ,j  ,k  ),corr_loc(i  ,j+1,k  )) * sWGP(j) * (gtilde_DG(:,i,  j,k,eID) - gtilde_FV(:,i  ,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr * sdt
      
      ! zeta correction
      ! ---------------
      ! left
      my_corr=-max(corr_loc(i  ,j  ,k-1),corr_loc(i  ,j  ,k  )) * sWGP(k) * (htilde_DG(:,i,j,k-1,eID) - htilde_FV(:,i,j,k-1,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr * sdt
      
      ! right
      my_corr=max(corr_loc(i  ,j  ,k  ),corr_loc(i  ,j  ,k+1)) * sWGP(k) * (htilde_DG(:,i,  j,k,eID) - htilde_FV(:,i  ,j,k,eID))*sJ(i,j,k,eID)
      U (:,i,j,k) = U (:,i,j,k) + my_corr
      Ut(:,i,j,k) = Ut(:,i,j,k) + my_corr * sdt
      
    end do       ; end do       ; enddo
#else
    ! Element-wise correction!
    do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
      ! Correct U
      U (:,i,j,k) = U (:,i,j,k) + corr * FFV_m_FDG(:,i,j,k,eID)
      ! Correct Ut
      Ut(:,i,j,k) = Ut(:,i,j,k) + (alpha-alphacont) * FFV_m_FDG(:,i,j,k,eID)
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
    SDEALLOCATE( Flux_ext  )
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
    SDEALLOCATE ( rf_DG )
    SDEALLOCATE ( rg_DG )
    SDEALLOCATE ( rh_DG )
#endif /*LOCAL_ALPHA*/
    
  end subroutine Finalize_IDP
#endif /*NFVSE_CORR*/  
end module MOD_IDP

