!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_adc
!
!> \brief MPAS ocean vertical mixing interface to new ADC closure
!> \author Luke Van Roekel
!> \date   10 April 2020
!> \details
!>  This module contains the routines for calls into the new ADC library
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_adc

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_log
   use mpas_vector_reconstruction
   use mpas_timekeeping

   use ocn_adc_mixing_fused
   use ocn_surface_bulk_forcing
   use ocn_surface_land_ice_fluxes
   use ocn_turbulence
   use ocn_constants
   use ocn_config

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vmix_adc_mix, &
             ocn_vmix_adc_init, &
             ocn_vmix_adc_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: adcOn
   integer :: local_nCells
   integer, parameter :: num_water_types = 5
   real (kind=RKIND), dimension(num_water_types) ::                       &
      rfac   = (/ 0.58_RKIND, 0.62_RKIND, 0.67_RKIND, 0.77_RKIND, 0.78_RKIND /), &
      depth1 = (/ 0.35_RKIND, 0.60_RKIND, 1.00_RKIND, 1.50_RKIND, 1.40_RKIND /), &
      depth2 = (/ 23.0_RKIND, 20.0_RKIND, 17.0_RKIND, 14.0_RKIND, 7.90_RKIND /)

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_adc_mix
!
!> \brief   Computes mixing tendency using the assumed distribution closure
!> \author  Luke Van Roekel
!> \date    10 April 2020
!> \details
!>  This routine computes the vertical mixing tendencies for momentum
!>  and tracers by calling ADC mixing library
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_adc_mix(meshPool, statePool, forcingPool, diagnosticsPool, &
      mixedLayerDepthsAMPool, dt, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, intent(in), optional :: timeLevelIn !< Input: time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool         !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         diagnosticsPool   !< Input/Output: diagnostic information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool   !< Input/Output: forcing information

      type (mpas_pool_type), intent(in) :: &
         mixedLayerDepthsAMPool

      type (mpas_pool_type), pointer :: tracersPool, tracersSurfaceFluxPool

      real(KIND=RKIND), intent(in) :: dt
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer, dimension(:), pointer :: &
        maxLevelCell, nEdgesOnCell

      real (kind=RKIND), dimension(:), pointer :: &
        stressZonal, stressMeridional, fCell, &
        surfaceStress, surfaceStressMagnitude, bottomDepth,  &
        penetrativeTemperatureFlux, angleEdge, surfaceThicknessFlux, &
        dThreshMLD, surfaceThicknessFluxRunoff, boundaryLayerDepth

      real (kind=RKIND), dimension(:,:), pointer :: &
        normalVelocity, layerThickness, velocityZonal, velocityMeridional, &
        tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
        tracerGroupSurfaceFluxRemoved, BruntVaisalaFreqTop,  &
        inSituThermalExpansionCoeff, inSituSalineContractionCoeff

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer, pointer :: indexTemperature, indexSalinity

      integer :: k, i, iCell, jCell, iNeighbor, iter, timeLevel, kIndexOBL,  &
            ntracers, kav, iEdge, nCells, nEdges, cell1, cell2
      integer, pointer :: nVertLevels, nVertLevelsP1
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), allocatable :: uw_sfc, vw_sfc, ws_sfc, wt_sfc

      real (kind=RKIND), dimension(:,:), allocatable :: utend, vtend, ttend, stend
      real (kind=RKIND) :: fac, dep1, dep2, fracAbsorbed, fracAbsorbedRunoff, &
                           utendEdge, vtendEdge

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing mixing-related fields
      ! note that the user can choose multiple options and the
      !   mixing fields have to be added/merged together
      !
      !-----------------------------------------------------------------

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      !
      ! only build up viscosity/diffusivity if CVMix is turned on
      !
      if ( .not. adcOn ) return

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      !
      ! set pointers for fields related to position on sphere
      !
      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFlux", tracerGroupSurfaceFlux)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFluxRunoff", tracerGroupSurfaceFluxRunoff)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFluxRemoved", tracerGroupSurfaceFluxRemoved)
      call mpas_pool_get_array(forcingPool, "surfaceThicknessFlux", surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, "surfaceThicknessFluxRunoff", surfaceThicknessFluxRunoff)

      call mpas_pool_get_array(mixedLayerDepthsAMPool, 'dThreshMLD', dThreshMLD)
      call mpas_pool_get_array(diagnosticsPool, "boundaryLayerDepth", boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'inSituThermalExpansionCoeff', inSituThermalExpansionCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'inSituSalineContractionCoeff', inSituSalineContractionCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
      !
      call mpas_pool_get_array(meshPool, 'fCell', fCell)
      !
      ! set pointers for fields related to vertical mesh
      !
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      !
      ! set pointers for fields related ocean state
      !
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop',BruntVaisalaFreqTop)

      ! need to use RBFs to get stress at cell centers
      call mpas_pool_get_array(forcingPool, 'windStressZonal', stressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', stressMeridional)

      call mpas_pool_get_array(forcingPool, 'surfaceStress', surfaceStress)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)

      nCells = nCellsArray( 4 )
      nEdges = nEdgesArray( 3 )

      allocate(uw_sfc(nCells), vw_sfc(nCells), ws_sfc(nCells), wt_sfc(nCells))
      allocate(utend(nVertLevels,nCells), vtend(nVertLevels,nCells))
      allocate(ttend(nVertLevels,nCells), stend(nVertLevels,nCells))

      !$omp do schedule(runtime)
      do iEdge=1,nEdges
         surfaceStress(iEdge) = 0.0_RKIND
      end do
      !$omp end do

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         surfaceThicknessFlux(iCell) = 0.0_RKIND
         surfaceThicknessFluxRunoff(iCell) = 0.0_RKIND
         surfaceStressMagnitude(iCell) = 0.0_RKIND
         tracerGroupSurfaceFluxRunoff(:, iCell) = 0.0_RKIND
         tracerGroupSurfaceFlux(:, iCell) = 0.0_RKIND
         tracerGroupSurfaceFluxRemoved(:, iCell) = 0.0_RKIND
      end do
      !$omp end do

      call ocn_surface_bulk_forcing_vel(meshPool, forcingPool, surfaceStress, surfaceStressMagnitude, err)

      call ocn_surface_land_ice_fluxes_vel(meshPool, diagnosticsPool, surfaceStress, surfaceStressMagnitude, err)

      if (config_use_activeTracers_surface_bulk_forcing) then
         call ocn_surface_bulk_forcing_tracers(meshPool, "activeTracers", forcingPool, activeTracers, &
                                               tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
                                               tracerGroupSurfaceFluxRemoved, dt, layerThickness, err)
      end if

      call ocn_surface_land_ice_fluxes_tracers(meshPool, "activeTracers", forcingPool, tracerGroupSurfaceFlux, err)

      call ocn_surface_bulk_forcing_thick(meshPool, forcingPool, surfaceThicknessFlux, surfaceThicknessFluxRunoff, err)

      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)

      ! water type parameters
      fac = rfac(config_jerlov_water_type)
      dep1 = depth1(config_jerlov_water_type)
      dep2 = depth2(config_jerlov_water_type)

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         ! set surface forcing
         fracAbsorbed = 1.0_RKIND - exp( max(-layerThickness(1, iCell) / config_flux_attenuation_coefficient, -100.0_RKIND) )
         fracAbsorbedRunoff = 1.0_RKIND - exp( max(-layerThickness(1, iCell) / config_flux_attenuation_coefficient_runoff, -100.0_RKIND) )

         wt_sfc(iCell) = tracerGroupSurfaceFlux(indexTemperature, iCell) &
                       - fracAbsorbed * surfaceThicknessFlux(iCell) * activeTracers(indexTemperature,1,iCell) &
                       + fracAbsorbedRunoff * tracerGroupSurfaceFluxRunoff(indexTemperature,iCell)

         ws_sfc(iCell) = tracerGroupSurfaceFlux(indexSalinity, iCell) &
                       - fracAbsorbed * surfaceThicknessFlux(iCell) * activeTracers(indexSalinity,1,iCell)

         uw_sfc(iCell) = stressZonal(iCell) / rho_sw

         vw_sfc(iCell) = stressMeridional(iCell) / rho_sw

         ! set vertical grid
         ze(1,iCell) = 0.0_RKIND
         zm(1,iCell) = -layerThickness(1,iCell)/2.0_RKIND
         do k=2,maxLevelCell(iCell)
            ze(k,iCell) = ze(k-1,iCell) - layerThickness(k-1,iCell)
            zm(k,iCell) = ze(k, iCell) - layerThickness(k,iCell)/2.0_RKIND
         end do
         k = maxLevelCell(iCell)+1
         ze(k,iCell) = ze(k-1,iCell) - layerThickness(k-1,iCell)
         do k = maxLevelCell(iCell) + 1, nVertLevels
            ze(k+1,iCell) = ze(maxLevelCell(iCell)+1,iCell)
            zm(k,iCell) = ze(maxLevelCell(iCell)+1,iCell)
         end do
      end do
      !$omp end do

      !this will update T in place, so verify it sends time level n+1

      ntracers = size(activeTracers, dim=1)

      if ( config_adc_use_single_column ) then
         local_nCells = 1
      else
         local_nCells = nCells
      end if

      ! use xtend arrays to store velocity and tracers before ADC
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         do k = 1, nVertLevels
            utend(k,iCell) = velocityZonal(k,iCell)
            vtend(k,iCell) = velocityMeridional(k,iCell)
            ttend(k,iCell) = activeTracers(indexTemperature,k,iCell)
            stend(k,iCell) = activeTracers(indexSalinity,k,iCell)
         end do
      end do
      !$omp end do

      ! do ADC
      call compute_ADC_tends(local_nCells, nVertLevels, ntracers, dt, &
                             activeTracers, velocityZonal, velocityMeridional, BruntVaisalaFreqTop, &
                             uw_sfc, vw_sfc, -wt_sfc, ws_sfc, &
                             inSituThermalExpansionCoeff, inSituSalineContractionCoeff, fCell, dThreshMLD)

      ! compute the tendencies
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         do k = 1, nVertLevels
            utend(k,iCell) = (velocityZonal(k,iCell)      - utend(k,iCell)) / dt
            vtend(k,iCell) = (velocityMeridional(k,iCell) - vtend(k,iCell)) / dt
         end do
      end do
      !$omp end do

      !need to add interpolation back to N. Velocity!
      if ( config_adc_use_single_column ) then
         !$omp do schedule(runtime)
         do iCell = 2, nCells
            do k = 1, nVertLevels
               activeTracers(indexTemperature,k,iCell) = activeTracers(indexTemperature,k,1)
               activeTracers(indexSalinity,k,iCell) = activeTracers(indexSalinity,k,1)
               utend(k,iCell) = utend(k,1)
               vtend(k,iCell) = vtend(k,1)
            end do
         end do
         !$omp end do
      end if

      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k = 1, nVertLevels
            utendEdge = 0.5_RKIND*(utend(k,cell1) + utend(k,cell2))
            vtendEdge = 0.5_RKIND*(vtend(k,cell1) + vtend(k,cell2))
            normalVelocity(k,iEdge) = normalVelocity(k,iEdge) &
                                    + dt * (cos(angleEdge(iEdge)) * utendEdge + sin(angleEdge(iEdge)) * vtendEdge)
         end do
      end do
      !$omp end do

      deallocate(uw_sfc, vw_sfc, ws_sfc, wt_sfc)
      deallocate(utend,vtend,ttend,stend)

   end subroutine ocn_vmix_adc_mix!}}}

!***********************************************************************
!
!  routine ocn_vmix_adc_init
!
!> \brief   Initializes ocean vertical mixing quantities by using
!> \ get and puts into ADC library
!> \author  Luke Van Roekel
!> \date    10 April 2020
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical mixing in the ocean. Parameters are set by calling into ADC library
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_adc_init(domain,err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      type (mpas_pool_type), pointer :: meshPool, statePool
      type (block_type), pointer :: block
      type (MPAS_TimeInterval_type) :: timeStep
      integer, intent(out) :: err !< Output: error flag

      integer, pointer :: nVertLevels
      integer :: k, iCell, nCells, err_tmp, nEdges
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      real(KIND=RKIND),dimension(:,:),pointer :: layerThickness

      real(KIND=RKIND) :: dt, numiterations

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      adcOn = config_use_adc

      !
      ! only initialize if Assumed Distribution Closure is turned on
      !
      if (.not.config_use_adc) return

      if (config_use_cvmix) then
         call mpas_log_write( "Cannot use ADC mixing and CVMix at the same time:" // &
               "set either config_use_cvmix or config_use_adc to .false.", MPAS_LOG_CRIT)
      end if

      if (.not. config_AM_mixedLayerDepths_enable) then
         call mpas_log_write("ADC mixing requires the MLD calculation to be enabled!", MPAS_LOG_CRIT)
      end if

      !
      ! pull nVertLevels out of the mesh structure
      !

      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      call mpas_get_timeInterval(timeStep, dt=dt)

      if ( mod(dt, config_adc_timestep) .ne. 0 ) then
         numiterations = ceiling(dt/config_adc_timestep)
         config_adc_timestep = dt / numiterations
         call mpas_log_write( "WARNING: timestep for ADC scheme was not an even multiple of" // &
               "config_dt, config_adc_timestep has been reset appropriately", MPAS_LOG_WARN)
      end if

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      call ocn_turbulenceCreate(domain)

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
         nCells = nCellsArray( size(nCellsArray) )
         do iCell = 1, nCells
            do k = 1, nVertLevels
               w3(:,k,iCell) = 0.0_RKIND
            end do
            do k = 1, nVertLevels+1
               KspsU(:,k,iCell) = epsilon
               KspsD(:,k,iCell) = epsilon
               w2(:,k,iCell) = 0.0_RKIND
               areaFraction(k,iCell) = 0.5_RKIND
               tumd(k,iCell) = 0.0_RKIND
               sumd(k,iCell) = 0.0_RKIND
               wumd(k,iCell) = 0.0_RKIND
               uw(:,k,iCell) = 0.0_RKIND
               vw(:,k,iCell) = 0.0_RKIND
               u2(:,k,iCell) = 0.0_RKIND
               v2(:,k,iCell) = 0.0_RKIND
               uv(:,k,iCell) = 0.0_RKIND
               ut(:,k,iCell) = 0.0_RKIND
               vt(:,k,iCell) = 0.0_RKIND
               wt(:,k,iCell) = 0.0_RKIND
               us(:,k,iCell) = 0.0_RKIND
               vs(:,k,iCell) = 0.0_RKIND
               ws(:,k,iCell) = 0.0_RKIND
               eps(:,k,iCell) = 0.0_RKIND
            end do
            length(1,iCell) = 0.5_RKIND*layerThickness(1,iCell)
            do k=2,nVertLevels
               length(k,iCell) = 0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))
            end do
            length(nVertLevels+1,iCell) = 0.5_RKIND*layerThickness(nVertLevels,iCell)
         end do
         block => block % next
      end do

   end subroutine ocn_vmix_adc_init!}}}

!***********************************************************************

   subroutine ocn_vmix_adc_finalize(ierr)

     integer, intent(out) :: ierr

     deallocate(KspsU, KspsD, ze, zm,          &
                eps, length, lenspsD, lenspsU, &
                KhU, KhD, KmU, KmD,            &
                wt_spsU, wt_spsD, ws_spsU, ws_spsD, &
                uw2, vw2, u2w, v2w, w2t,  &
                w2s, wts, uvw, uwt, vwt,  &
                ws2, wt2, areaFraction, Entrainment, &
                Detrainment, tumd, sumd, wumd,       &
                w2tend1, w2tend2, w2tend3, w2tend4,  &
                w2tend5, wttend1, wttend2, wttend3,  &
                wttend4, wttend5, wstend1, wstend2,  &
                wstend3, wstend4, wstend5, uwtend1,  &
                uwtend2, uwtend3, uwtend4, uwtend5,  &
                vwtend1, vwtend2, vwtend3, vwtend4,  &
                vwtend5, u2tend1, u2tend2, u2tend3,  &
                u2tend4, u2tend5, v2tend1, v2tend2,  &
                v2tend3, v2tend4, v2tend5, u2cliptend, &
                v2cliptend, w2cliptend, u2, v2,      &
                w2, wt, ws, uw, vw,        &
                w3, uv, t2, s2, ut,        &
                vt, us, vs, ts)

   end subroutine ocn_vmix_adc_finalize

end module ocn_vmix_adc

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
