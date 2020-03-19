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
!> Computes the spatial operator using Native Finite Volume Sub-Elements (NFVSE)
!> Based on: Hennemann and Gassner (2020). "Entropy stable shock capturing for the discontinuous galerkin spectral element
!>                                          method with native finite volume sub elements"
!> Attention 1: The current strategy is to perform the operations on "inner faces". This can be improved, but the Riemann solver 
!>              routines have to be adapted
!==================================================================================================================================
module MOD_NFVSE
#if SHOCK_NFVSE
  use MOD_PreProc
  implicit none
  
  private
  public :: VolInt_NFVSE, InitNFVSE, FinalizeNFVSE, Apply_NFVSE_Correction
  
contains
!===================================================================================================================================
!> Initializes the NFVSE module
!===================================================================================================================================
  subroutine InitNFVSE()
    use MOD_NFVSE_Vars         , only: SubCellMetrics, sWGP, MPIRequest_alpha, Fsafe, Fblen
    use MOD_MPI_Vars           , only: nNbProcs
    use MOD_Mesh_Vars          , only: nElems, Metrics_fTilde, Metrics_gTilde, Metrics_hTilde
    use MOD_Interpolation_Vars , only: wGP
    use MOD_ShockCapturing_Vars, only: alpha_old
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES 
    integer :: iElem
    integer :: i,j,k      !DOF counters
    real, parameter :: half = 0.5d0
    !--------------------------------------------------------------------------------------------------------------------------------
    
    ! Allocate storage
    allocate ( SubCellMetrics(nElems) )
    call SubCellMetrics % construct(PP_N)
    
    allocate ( sWGP(0:PP_N) )
    
    ! Compute inner normal and tangent vectors (the storage access to Metrics_*Tilde is not optimized... but this is done only once!)
    do iElem=1, nElems
      
      do i=0, PP_N-1
        ! Compute vectors
        SubCellMetrics(iElem) % xi   % nv(:,:,:,i) = half * ( Metrics_fTilde(:,i,:,:,iElem) + Metrics_fTilde(:,i+1,:,:,iElem) )
        SubCellMetrics(iElem) % xi   % t1(:,:,:,i) = half * ( Metrics_gTilde(:,i,:,:,iElem) + Metrics_gTilde(:,i+1,:,:,iElem) )
        SubCellMetrics(iElem) % xi   % t2(:,:,:,i) = half * ( Metrics_hTilde(:,i,:,:,iElem) + Metrics_hTilde(:,i+1,:,:,iElem) )
        
        SubCellMetrics(iElem) % eta  % nv(:,:,:,i) = half * ( Metrics_gTilde(:,:,i,:,iElem) + Metrics_gTilde(:,:,i+1,:,iElem) )
        SubCellMetrics(iElem) % eta  % t1(:,:,:,i) = half * ( Metrics_hTilde(:,:,i,:,iElem) + Metrics_hTilde(:,:,i+1,:,iElem) )
        SubCellMetrics(iElem) % eta  % t2(:,:,:,i) = half * ( Metrics_fTilde(:,:,i,:,iElem) + Metrics_fTilde(:,:,i+1,:,iElem) )
        
        SubCellMetrics(iElem) % zeta % nv(:,:,:,i) = half * ( Metrics_hTilde(:,:,:,i,iElem) + Metrics_hTilde(:,:,:,i+1,iElem) )
        SubCellMetrics(iElem) % zeta % t1(:,:,:,i) = half * ( Metrics_fTilde(:,:,:,i,iElem) + Metrics_fTilde(:,:,:,i+1,iElem) )
        SubCellMetrics(iElem) % zeta % t2(:,:,:,i) = half * ( Metrics_gTilde(:,:,:,i,iElem) + Metrics_gTilde(:,:,:,i+1,iElem) )
        
        ! Normalize each
        do k=0, PP_N ; do j=0, PP_N
          SubCellMetrics(iElem) % xi   % norm(j,k,i) = norm2 (SubCellMetrics(iElem) % xi   % nv(:,j,k,i))
          SubCellMetrics(iElem) % xi   % nv(:,j,k,i) = SubCellMetrics(iElem) % xi   % nv(:,j,k,i) / SubCellMetrics(iElem) % xi   % norm(j,k,i)
          SubCellMetrics(iElem) % xi   % t1(:,j,k,i) = SubCellMetrics(iElem) % xi   % t1(:,j,k,i) / norm2 (SubCellMetrics(iElem) % xi   % t1(:,j,k,i))
          SubCellMetrics(iElem) % xi   % t2(:,j,k,i) = SubCellMetrics(iElem) % xi   % t2(:,j,k,i) / norm2 (SubCellMetrics(iElem) % xi   % t2(:,j,k,i))
          
          SubCellMetrics(iElem) % eta  % norm(j,k,i) = norm2 (SubCellMetrics(iElem) % eta  % nv(:,j,k,i))
          SubCellMetrics(iElem) % eta  % nv(:,j,k,i) = SubCellMetrics(iElem) % eta  % nv(:,j,k,i) / SubCellMetrics(iElem) % eta  % norm(j,k,i)
          SubCellMetrics(iElem) % eta  % t1(:,j,k,i) = SubCellMetrics(iElem) % eta  % t1(:,j,k,i) / norm2 (SubCellMetrics(iElem) % eta  % t1(:,j,k,i))
          SubCellMetrics(iElem) % eta  % t2(:,j,k,i) = SubCellMetrics(iElem) % eta  % t2(:,j,k,i) / norm2 (SubCellMetrics(iElem) % eta  % t2(:,j,k,i))
          
          SubCellMetrics(iElem) % zeta % norm(j,k,i) = norm2 (SubCellMetrics(iElem) % zeta % nv(:,j,k,i))
          SubCellMetrics(iElem) % zeta % nv(:,j,k,i) = SubCellMetrics(iElem) % zeta % nv(:,j,k,i) / SubCellMetrics(iElem) % zeta % norm(j,k,i)
          SubCellMetrics(iElem) % zeta % t1(:,j,k,i) = SubCellMetrics(iElem) % zeta % t1(:,j,k,i) / norm2 (SubCellMetrics(iElem) % zeta % t1(:,j,k,i))
          SubCellMetrics(iElem) % zeta % t2(:,j,k,i) = SubCellMetrics(iElem) % zeta % t2(:,j,k,i) / norm2 (SubCellMetrics(iElem) % zeta % t2(:,j,k,i))
        end do       ; end do
      end do
    end do
    
    ! Compute the inverse of the quadrature weights (sub-cell dimensions)
    sWGP = 1.d0 / wGP
    
#if NFVSE_CORR
    allocate ( Fsafe(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) )
    allocate ( Fblen(1:PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) )
    allocate ( alpha_old(nElems) )
    alpha_old = 0.
#endif /*NFVSE_CORR*/
    
#if MPI
    allocate(MPIRequest_alpha(nNbProcs,4)    ) ! 1: send slave, 2: send master, 3: receive slave, 4, receive master
#endif
  end subroutine InitNFVSE
  
!===================================================================================================================================
!> Computes the "volume integral": Spatial contribution to Ut by the subcell finite volumes
!> Attention 1: 1/J(i,j,k) is not yet accounted for
!> Attention 2: input Ut=0. and is updated with the volume flux derivatives
!> Attention 3: The element boundaries are neglected in this routine
!===================================================================================================================================
  subroutine VolInt_NFVSE(Ut)
  !----------------------------------------------------------------------------------------------------------------------------------
    ! Modules
    use MOD_PreProc
    use MOD_DG_Vars            , only: U
    use MOD_Mesh_Vars          , only: nElems
    use MOD_NFVSE_Vars         , only: SubCellMetrics, sWGP, Fsafe, Fblen
    use MOD_ShockCapturing_Vars, only: alpha, alpha_max
    use MOD_Basis              , only: ALMOSTEQUAL
    use MOD_NFVSE_MPI          , only: PropagateBlendingCoeff
    ! IMPLICIT VARIABLE HANDLING
    implicit none
    !-------------------------------------------------------------------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    real,intent(inout)                                 :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems)
    !< Adds volume contribution to time derivative Ut contained in MOD_DG_Vars 
    !-------------------------------------------------------------------------------------------------------------------------------
    ! LOCAL VARIABLES
    real,dimension(PP_nVar,-1:PP_N,-1:PP_N,-1:PP_N) :: ftilde, gtilde, htilde ! transformed flux inter-subcell fluxes (with ghost cells)
    real,dimension(PP_nVar)                         :: F_FV
    integer                                         :: i,j,k,iElem
    !===============================================================================================================================
    
    call PropagateBlendingCoeff()
    
    ftilde = 0.d0
    gtilde = 0.d0
    htilde = 0.d0
    
    do iElem=1,nElems
#if !defined(NFVSE_CORR)
      if ( ALMOSTEQUAL(alpha(iElem),0.d0) ) cycle
#endif /*NFVSE_CORR*/
      !compute inner Riemann solutions (TODO: Not optimal.. Improve)
      call Compute_VolFluxes( U(:,:,:,:,iElem), SubCellMetrics(iElem), ftilde(:,0:PP_N-1, 0:PP_N  , 0:PP_N  ), &
                                                                       gtilde(:,0:PP_N  , 0:PP_N-1, 0:PP_N  ), &
                                                                       htilde(:,0:PP_N  , 0:PP_N  , 0:PP_N-1)  )
      
      ! Update the inner nodes
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        ! Get Finite Volume Ut
        F_FV =   sWGP(i) * ( ftilde(:,i,j,k) - ftilde(:,i-1,j  ,k  ) ) &
               + sWGP(j) * ( gtilde(:,i,j,k) - gtilde(:,i  ,j-1,k  ) ) &
               + sWGP(k) * ( htilde(:,i,j,k) - htilde(:,i  ,j  ,k-1) )
#if NFVSE_CORR
        Fsafe(:,i,j,k,iElem) = (1. - alpha_max   ) * Ut(:,i,j,k,iElem) +  alpha_max    * F_FV
        Fblen(:,i,j,k,iElem) = (1. - alpha(iElem)) * Ut(:,i,j,k,iElem) +  alpha(iElem) * F_FV
        Ut(:,i,j,k,iElem)    = Fblen(:,i,j,k,iElem)
#else
        Ut   (:,i,j,k,iElem) = (1. - alpha(iElem)) * Ut(:,i,j,k,iElem) +  alpha(iElem) * F_FV
#endif /*NFVSE_CORR*/
      end do       ; end do       ; end do ! i,j,k
      
    end do ! iElem
  end subroutine VolInt_NFVSE
  
!============================================================================================================================
!> Get Pressure
!============================================================================================================================
  pure subroutine GetPressure(U,p)
#ifdef mhd
    use MOD_Equation_Vars      , only: s2mu_0
#endif /*mhd*/
    use MOD_Equation_Vars      , only: KappaM1
    implicit none
    real, intent(in)  :: U(PP_nVar)
    real, intent(out) :: p
    
    p = KappaM1*(U(5)-0.5*(SUM(U(2:4)*U(2:4))/U(1)))
#ifdef mhd
    p = p - KappaM1*s2mu_0*SUM(U(6:8)*U(6:8))
#ifdef PP_GLM
    p = p - 0.5*KappaM1*U(9)*U(9)
#endif /*PP_GLM*/
#endif /*mhd*/
  end subroutine GetPressure
!===================================================================================================================================
!> Corrects the Solution after the Runge-Kutta stage
!===================================================================================================================================
  subroutine Apply_NFVSE_Correction(U,t,dt)
    use MOD_ShockCapturing_Vars, only: alpha, alpha_max, beta, alpha_old
    use MOD_NFVSE_Vars         , only: Fsafe, Fblen
    use MOD_Mesh_Vars          , only: nElems, offsetElem
    use MOD_Basis              , only: ALMOSTEQUAL
    use MOD_Equation_Vars      , only: KappaM1, Kappa
    use MOD_Mesh_Vars          , only: sJ
    use MOD_NFVSE_MPI
    USE MOD_Globals
    implicit none
    !-arguments----------------------------------------------
    real,intent(inout) :: U (PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< Current solution (in RK stage)
    real,intent(in)    :: t                                         !< Current time (in time step!)
    real,intent(in)    :: dt                                        !< Current RK time-step size (in RK stage)
    !-local-variables----------------------------------------
    real    :: Usafe(PP_nVar,0:PP_N,0:PP_N,0:PP_N)
    real    :: a   ! a  = beta * rho_safe - rho
    real    :: ap  ! ap = (beta * p_safe   - p) / (kappa-1)
    real    :: p   (nElems), pres, lowest_psafe, lowest_pres, lowest_rhosafe, lowest_rho
    real    :: corr, corr1
    real    :: p_safe(0:PP_N,0:PP_N,0:PP_N) ! pressure obtained with alpha_max for an element
    real    :: alphadiff
    real    :: a_max, sdt
    real :: pprev, rhoprev
    real :: max_pdev, pdev
    real :: minS, minSm, maxS, maxSm, ent
    real, parameter :: eps = 1.e-8
    integer :: a_loc(3), ijkl(4)
    integer :: corrElems(nElems)
    integer :: numCorrElems
    integer :: eID
    integer :: i,j,k
    integer :: fID
    integer :: iter
    integer, parameter :: MAX_NEWTON_ITER = 1
    !--------------------------------------------------------
    
    numCorrElems = 0
    corrElems = 0
    alpha_old = alpha 
    sdt = 1./dt
    
    !<debug
!#    p = huge(1.)
!#    minS  = huge(1.)
!#    minSm = huge(1.)
!#    maxS  =-huge(1.)
!#    maxSm =-huge(1.)
!#    do eID=1, nElems
!#      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
!#        call GetPressure(U(:,i,j,k,eID),pres)
!#        p(eID) = min(p(eID),pres)
!#        ent = log(pres) - Kappa * log(U(1,i,j,k,eID))
!#        minS = min (minS, ent)
!#        maxS = max (maxS, ent)
!#        ent = -ent * U(1,i,j,k,eID) / KappaM1
!#        minSm = min (minSm, ent)
!#        maxSm = max (maxSm, ent)
!#      end do       ; end do       ; end do ! i,j,k
!#    end do
!#    pprev = minval(p)
!#    rhoprev = minval(U(1,:,:,:,:))
!~     print*, '####################'
!~     print*, '### Gonna correct?' !, rhoprev, pprev, minS, minSm, maxS, maxSm
    !debug>
!
!   First do the easy density correction (with the hope that it solves the problem)
!   *******************************************************************************

    do eID=1, nElems
      ! Check if it makes sense correcting
      alphadiff = alpha(eID) - alpha_max
      if ( abs(alphadiff) < eps ) cycle ! Not much to do for this element...
      
      ! Compute a (density correction)
      alphadiff = 1./alphadiff
      corr = 0.
      
!#      !<debug
!#      lowest_psafe = huge(1.)
!#      lowest_pres = huge(1.)
!#      lowest_rhosafe = huge(1.)
!#      !debug>
      
      ! Compute corrections
      do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        ! Density correction
        Usafe(:,i,j,k) = U(:,i,j,k,eID) + dt * ( Fsafe(:,i,j,k,eID) - Fblen(:,i,j,k,eID) ) * -sJ(i,j,k,eID)
        
        a = (beta * Usafe(1,i,j,k) - U(1,i,j,k,eID))
        if (a > 0.) then ! This DOF doesn't need correction
          corr1 = (Fblen(1,i,j,k,eID) - Fsafe(1,i,j,k,eID)) * alphadiff * -sJ(i,j,k,eID)
          corr1 = a / corr1
          if (corr1 > corr) then
            corr = corr1
            a_loc = [i,j,k]
          end if
!~           corr = max(corr,corr1)
        end if
        
        !<<<Initial pressure correction
        call GetPressure(U(:,i,j,k,eID),pres)
        call GetPressure(Usafe(:,i,j,k),p_safe(i,j,k))
        if (p_safe(i,j,k) < 0.) then
          print*, 'ERROR: safe pressure not safe el=', eID+offsetElem, p_safe(i,j,k)
          stop
        end if
        if (Usafe(1,i,j,k) < 0.) then
          print*, 'ERROR: safe dens not safe el=', eID+offsetElem, Usafe(1,i,j,k)
          stop
        end if
        ap = (beta * p_safe(i,j,k) - pres) / KappaM1
        if (ap > 0.) then
          corr1 = (Fblen(5,i,j,k,eID) - Fsafe(5,i,j,k,eID)) * alphadiff * -sJ(i,j,k,eID)
          corr1 = ap / corr1
          if (corr1 > corr) then
            corr = corr1
            a_loc = [i,j,k]
          end if
        end if
        !Initial pressure correction>>>
        
!#        !<debug
!#!~         call GetPressure(U(:,i,j,k,eID),pres)
!#        lowest_pres = min ( lowest_pres, pres)
!#        lowest_psafe = min (lowest_psafe,p_safe(i,j,k))
!#        lowest_rhosafe = min ( lowest_rhosafe, Usafe(1,i,j,k))
!#        !debug>
      end do       ; end do       ; end do ! i,j,k
      lowest_rho = minval (U(1,:,:,:,eID))
      
      ! Check if the density went below the allowed value and correct
      if ( corr > 0. ) then
!#        !<debug
!#        print*, '####################'
!#        print*, '### CORRECTING', rhoprev, pprev
!#        call GetPressure(U(:,a_loc(1),a_loc(2),a_loc(3),eID),pres)
!#        print*, '---'
!#        print*, '* s', p_safe(a_loc(1),a_loc(2),a_loc(3)), Usafe(1,a_loc(1),a_loc(2),a_loc(3))
!#        print*, '* _', pres, U(1,a_loc(1),a_loc(2),a_loc(3),eID)
!#        !debug>
        
        alpha_old(eID) = alpha(eID)
        
        ! Change the alpha for output
        alpha(eID) = alpha_old(eID) + corr * sdt
        
        ! Change inconsistent alphas
        if (alpha(eID) > alpha_max) then
          alpha(eID) = alpha_max
          corr = (alpha_max - alpha_old(eID) ) * dt
        end if
        
        ! Correct!
        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
        U(:,i,j,k,eID) = U(:,i,j,k,eID) + corr * (Fblen(:,i,j,k,eID) - Fsafe(:,i,j,k,eID)) * -sJ(i,j,k,eID) / (alpha_old(eID) - alpha_max)
        end do       ; end do       ; enddo
!#        !<debug
!#        call GetPressure(U(:,a_loc(1),a_loc(2),a_loc(3),eID),pres)
!#        print*, '*af', pres, U(1,a_loc(1),a_loc(2),a_loc(3),eID)
!#        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
!#          call GetPressure(U(:,i,j,k,eID),pres)
!#          if(pres           +eps < beta *  p_safe(i,j,k)) print*, 'WARNING: p  is too low', i, j, k, pres,  beta * p_safe(i,j,k)
!#          if(U(1,i,j,k,eID) +eps < beta * Usafe(1,i,j,k)) print*, 'WARNING:rho is too low', i, j, k, U(1,i,j,k,eID), beta * Usafe(1,i,j,k)
!#        end do       ; end do       ; end do ! i,j,k
!#        !debug>
        
        numCorrElems = numCorrElems + 1
        corrElems(numCorrElems) = eID + offsetElem
        
!#        print*, 'el ', eID + offsetElem, alpha_old(eID), alpha(eID)
!#        print*, '  s', lowest_psafe, lowest_rhosafe
!#        print*, '  _', lowest_pres, lowest_rho
        
!#        lowest_pres = huge(1.)
!#        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
!#          call GetPressure(U(:,i,j,k,eID),pres)
!#          lowest_pres = min ( lowest_pres, pres)
!#        end do       ; end do       ; end do ! i,j,k
!#        print*, ' af', lowest_pres, minval(U(1,:,:,:,eID))
        
      end if
      
    end do !eID
    
!#    ! Debug
!#    if (numCorrElems > 0) then
!#      p = huge(1.)
!#      do eID=1, nElems
!#        do k=0, PP_N ; do j=0, PP_N ; do i=0, PP_N
!#          call GetPressure(U(:,i,j,k,eID),pres)
!#          p(eID) = min(p(eID),pres)
!#        end do       ; end do       ; end do ! i,j,k
!#        if (eID == 3254) then
!#          print*, 'el 3254. rho, p:', minval(U(1,:,:,:,3254)),p(3254)
!#        end if
!#      end do
!#      WRITE(UNIT_StdOut,'(A,ES16.7,a,ES16.7)') '-->Corrected alphas!! at t=', t, 'dt_rk=', dt
!#      print*, 'new: rho_min, p_min=', minval(U(1,:,:,:,:)), minval(p)
!#      print*, 'old: rho_min, p_min=', rhoprev, pprev
!#      print*, 'minloc:  rho, p     ', minloc(U(1,:,:,:,:)), minloc(p)
!#      print*, corrElems(1:numCorrElems)
!#    end if
  end subroutine Apply_NFVSE_Correction

!===================================================================================================================================
!> Solves the inner Riemann problems and outputs a FV consistent flux
!> Attention 1: It's not optimal to reshape the vectors (TODO: Improve)
!===================================================================================================================================
  subroutine Compute_VolFluxes(U,sCM,F,G,H)
    use MOD_PreProc
    use MOD_NFVSE_Vars, only: SubCellMetrics_t
    use MOD_Riemann   , only: Riemann
    implicit none
    !
    real,dimension(PP_nVar,0:PP_N  ,0:PP_N  ,0:PP_N  ), intent(in)    :: U         !< The element solution
    type(SubCellMetrics_t)                            , intent(in)    :: sCM       !< Sub-cell metric terms
    real,dimension(PP_nVar,0:PP_N-1,0:PP_N  ,0:PP_N  ), intent(inout) :: F
    real,dimension(PP_nVar,0:PP_N  ,0:PP_N-1,0:PP_N  ), intent(inout) :: G
    real,dimension(PP_nVar,0:PP_N  ,0:PP_N  ,0:PP_N-1), intent(inout) :: H
    !
    integer  :: i,j,k
    real :: U_(PP_nVar,0:PP_N,0:PP_N,0:PP_N  )
    real :: F_(PP_nVar,0:PP_N,0:PP_N,0:PP_N-1)
    
    !--
    
!   Xi-planes
!   ---------
    U_ = reshape(U , shape(U_), order = [1,4,2,3])
    do i=0, PP_N-1
      
      call Riemann(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1), &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R, &
#endif
                   sCM % xi   % nv(:,:,:,i),sCM % xi   % t1(:,:,:,i), sCM % xi   % t2(:,:,:,i))
      
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_(:,j,k,i) = F_(:,j,k,i) * sCM % xi   % norm(j,k,i)
      end do       ; end do
    end do
    F  = reshape(F_, shape(F ), order = [1,3,4,2])
    
!   Eta-planes
!   ----------
    U_ = reshape(U , shape(U_), order = [1,2,4,3])
    do i=0, PP_N-1
      
      call Riemann(F_(:,:,:,i),U_(:,:,:,i),U_(:,:,:,i+1), &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R, &
#endif
                   sCM % eta  % nv(:,:,:,i),sCM % eta  % t1(:,:,:,i), sCM % eta  % t2(:,:,:,i))
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        F_(:,j,k,i) = F_(:,j,k,i) * sCM % eta  % norm(j,k,i)
      end do       ; end do
    end do
    G  = reshape(F_, shape(G ), order = [1,2,4,3])
    
!   Zeta-planes
!   -----------
    do i=0, PP_N-1
      
      call Riemann(H(:,:,:,i),U(:,:,:,i),U(:,:,:,i+1), &
#if PARABOLIC
                   gradUx_L,gradUx_R,gradUy_L,gradUy_R,gradUz_L,gradUz_R, &
#endif
                   sCM % zeta % nv(:,:,:,i),sCM % zeta % t1(:,:,:,i), sCM % zeta % t2(:,:,:,i))
      ! Scale flux
      do k=0, PP_N ; do j=0, PP_N
        H (:,j,k,i) = H (:,j,k,i) * sCM % zeta % norm(j,k,i)
      end do       ; end do
    end do
    
  end subroutine Compute_VolFluxes
  
!===================================================================================================================================
!> Finalizes the NFVSE module
!===================================================================================================================================
  subroutine FinalizeNFVSE()
    use MOD_NFVSE_Vars, only: SubCellMetrics, sWGP, MPIRequest_alpha, Fsafe
    implicit none
    
    SDEALLOCATE (SubCellMetrics)
    SDEALLOCATE (sWGP)
    SDEALLOCATE (MPIRequest_alpha)
    SDEALLOCATE (Fsafe)
    
  end subroutine FinalizeNFVSE
#endif /*SHOCK_NFVSE*/
end module MOD_NFVSE

