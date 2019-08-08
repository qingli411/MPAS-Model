!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_LES
!
!> \brief MPAS ocean vertical mixing interface to PALM les
!> \author Luke Van Roekel
!> \date   09 January 2019
!> \details
!>  This module contains the routines for calls into LES
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_LES

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_log
   use mpas_vector_reconstruction

   use palm_mod
   use ocn_constants
   use ocn_surface_bulk_forcing
   use ocn_surface_land_ice_fluxes

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !-------------------------------------------------------------------

   logical :: first = .true.
   integer, parameter :: num_water_types = 5
   real (kind=RKIND), dimension(num_water_types) ::                       &
      rfac   = (/ 0.58_RKIND, 0.62_RKIND, 0.67_RKIND, 0.77_RKIND, 0.78_RKIND /), &
      depth1 = (/ 0.35_RKIND, 0.60_RKIND, 1.00_RKIND, 1.50_RKIND, 1.40_RKIND /), &
      depth2 = (/ 23.0_RKIND, 20.0_RKIND, 17.0_RKIND, 14.0_RKIND, 7.90_RKIND /)

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vmix_LES_interface, &
             ocn_vmix_LES_apply_tend

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_LES_interface
!
!> \brief   Compute mixing tendency due to LES
!> \author  Luke Van Roekel
!> \date    09 January 2019
!> \details
!>  This routine computes the sub grid mixing tendencies
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_LES_interface(meshPool, statePool, forcingPool, diagnosticsPool, mixedLayerDepthsAMPool, dt, err, timeLevelIn)!{{{

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

      type (mpas_pool_type), intent(in) :: mixedLayerDepthsAMPool

      type (mpas_pool_type), intent(inout) :: &
         diagnosticsPool   !< Input/Output: diagnostic information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool   !< Input/Output: forcing information

      type (mpas_pool_type), pointer :: tracersPool, tracersSurfaceFluxPool

      real(KIND=RKIND), intent(in) :: dt

!      logical, intent(in) :: first

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

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:), pointer :: &
        stressX, stressY, stressZ, stressZonal, stressMeridional, fCell, &
        surfaceStress, surfaceStressMagnitude,   &
        penetrativeTemperatureFlux, angleEdge, surfaceThicknessFlux, &
        surfaceThicknessFluxRunoff, dThreshMLD

      real (kind=RKIND), dimension(:,:), pointer :: &
        layerThickness, velocityZonal, velocityMeridional, &
        normalVelocity, uIncrementLES, vIncrementLES, tIncrementLES, &
        sIncrementLES, normalVelocityIncrementLES,     &
        tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
        tracerGroupSurfaceFluxRemoved, temperatureLES, zLES, &
        salinityLES, velocityZonalLES, velocityMeridionalLES, &
        zMid

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer, pointer :: indexTemperature, indexSalinity, config_jerlov_water_type
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer :: k, i, iCell, timeLevel, iEdge, nCells, nEdges
      integer :: nLevs, edgeCount, cell1, cell2
      integer, pointer :: nVertLevels, config_les_disturb_nblocks
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      character(len=strKind), pointer :: config_dt

      logical, pointer :: config_use_les_for_vmix, config_use_ecosysTracers,    &
                          config_use_activeTracers_surface_bulk_forcing,        &
                          config_AM_mixedLayerDepths_enable,                    &
                          constant_dz, mixed_layer_refine

      real (kind=RKIND), pointer :: config_les_dz, config_les_dx, config_les_dy,               &
                                    config_les_dt_data_output,   &
                                    config_les_dt_disturb, config_les_dt_data_output_av,       &
                                    config_les_dt_dopr, config_les_end_time,                   &
                                    config_les_dt_dots, config_flux_attenuation_coefficient,   &
                                    config_flux_attenuation_coefficient_runoff,                &
                                    config_les_disturb_top, config_les_disturb_bottom,         &
                                    config_les_disturb_amplitude, config_les_bottom,           &
                                    config_les_disturb_energy_limit, config_les_time_avg,      &
                                    config_les_restoring_strength

      integer, pointer :: config_les_nz, config_les_nx, config_les_ny

      real (kind=RKIND), dimension(:), allocatable :: zonalAverage, meridionalAverage, uw_sfc, &
                           vw_sfc, latVal, ws_sfc, wt_sfc, les_disturb_bottom, les_bottom

      real (kind=RKIND) :: fac, dep1, dep2, fracAbsorbed, fracAbsorbedRunoff,  &
                           etime

      err=0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 2
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_dt', config_dt)
      call mpas_pool_get_config(ocnConfigs, 'config_use_les_for_vmix', config_use_les_for_vmix)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers', config_use_ecosysTracers)
      call mpas_pool_get_config(ocnConfigs, 'config_use_activeTracers_surface_bulk_forcing', config_use_activeTracers_surface_bulk_forcing)
      call mpas_pool_get_config(ocnConfigs, 'config_AM_mixedLayerDepths_enable', config_AM_mixedLayerDepths_enable)

      if ( .not. config_use_les_for_vmix ) return

      if (config_use_ecosysTracers) then
         call mpas_log_write('Error, superparameterization is not configured for use with BGC',MPAS_LOG_CRIT)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_les_restoring_strength', &
                                config_les_restoring_strength)
      call mpas_pool_get_config(ocnConfigs, 'config_flux_attenuation_coefficient',  &
                               config_flux_attenuation_coefficient)
      call mpas_pool_get_config(ocnConfigs, 'config_flux_attenuation_coefficient_runoff',  &
                               config_flux_attenuation_coefficient_runoff)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dx', config_les_dx)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dy', config_les_dy)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dz', config_les_dz)
      call mpas_pool_get_config(ocnConfigs, 'config_les_nx', config_les_nx)
      call mpas_pool_get_config(ocnConfigs, 'config_les_ny', config_les_ny)
      call mpas_pool_get_config(ocnConfigs, 'config_les_nz', config_les_nz)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dt_data_output', config_les_dt_data_output)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dt_disturb', config_les_dt_disturb)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dt_data_output_av', config_les_dt_data_output_av)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dt_dopr', config_les_dt_dopr)
      call mpas_pool_get_config(ocnConfigs, 'config_les_end_time', config_les_end_time)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dt_dots', config_les_dt_dots)
      call mpas_pool_get_config(ocnConfigs, 'config_jerlov_water_type', config_jerlov_water_type)
      call mpas_pool_get_config(ocnConfigs, 'config_les_disturb_energy_limit', config_les_disturb_energy_limit)
      call mpas_pool_get_config(ocnConfigs, 'config_les_disturb_amplitude', config_les_disturb_amplitude)
      call mpas_pool_get_config(ocnConfigs, 'config_les_disturb_bottom', config_les_disturb_bottom)
      call mpas_pool_get_config(ocnConfigs, 'config_les_bottom', config_les_bottom)
      call mpas_pool_get_config(ocnConfigs, 'config_les_disturb_top', config_les_disturb_top)
      call mpas_pool_get_config(ocnConfigs, 'config_les_disturb_nblocks', config_les_disturb_nblocks)
      call mpas_pool_get_config(ocnConfigs, 'config_les_time_avg', config_les_time_avg)
      call mpas_pool_get_config(ocnConfigs, 'config_les_constant_dz', constant_dz)
      call mpas_pool_get_config(ocnConfigs, 'config_les_mixed_layer_refine', mixed_layer_refine)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)

      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFlux", tracerGroupSurfaceFlux)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFluxRunoff", tracerGroupSurfaceFluxRunoff)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFluxRemoved", tracerGroupSurfaceFluxRemoved)
      call mpas_pool_get_array(forcingPool, "surfaceThicknessFlux", surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, "surfaceThicknessFluxRunoff", surfaceThicknessFluxRunoff)
      call mpas_pool_get_array(meshPool, 'fCell', fCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'tIncrementLES', tIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'sIncrementLES', sIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'uIncrementLES', uIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'vIncrementLES', vIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityIncrementLES', normalVelocityIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

      call mpas_pool_get_array(mixedLayerDepthsAMPool, 'dThreshMLD', dThreshMLD)
      call mpas_pool_get_array(diagnosticsPool, 'zLES', zLES)
      call mpas_pool_get_array(diagnosticsPool, 'temperatureLES', temperatureLES)
      call mpas_pool_get_array(diagnosticsPool, 'salinityLES', salinityLES)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonalLES', velocityZonalLES)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridionalLES', velocityMeridionalLES)

      ! need to use RBFs to get stress at cell centers
      call mpas_pool_get_array(forcingPool, 'windStressZonal', stressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', stressMeridional)
      call mpas_pool_get_array(diagnosticsPool, 'stressX', stressX)
      call mpas_pool_get_array(diagnosticsPool, 'stressY', stressY)
      call mpas_pool_get_array(diagnosticsPool, 'stressZ', stressZ)

      call mpas_pool_get_array(forcingPool, 'surfaceStress', surfaceStress)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)
      etime = config_les_end_time

      if(constant_dz .and. mixed_layer_refine) then
        print *, 'Cannot set constant_dz and mixed_layer_refine, setting constant_dz to false'
        constant_dz = .false.
      endif

      allocate(zonalAverage(nVertLevels),meridionalAverage(nVertLevels))

      nCells = nCellsArray( 2 )
      nEdges = nEdgesArray( 1 )

      allocate(uw_sfc(nCells), vw_sfc(nCells), ws_sfc(nCells), wt_sfc(nCells))
      allocate(les_disturb_bottom(nCells), les_bottom(nCells),latVal(nCells))
      !$omp do schedule(runtime)
      do iEdge=1,nEdges
         surfaceStress(iEdge) = 0.0_RKIND
      enddo
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

      if(config_use_activeTracers_surface_bulk_forcing) then
        call ocn_surface_bulk_forcing_tracers(meshPool, "activeTracers", forcingPool, activeTracers,  &
                                            tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
                                            tracerGroupSurfaceFluxRemoved, dt, layerThickness, err)
      endif

      call ocn_surface_land_ice_fluxes_tracers(meshPool, "activeTracers", forcingPool, tracerGroupSurfaceFlux, err)

      call ocn_surface_bulk_forcing_thick(meshPool, forcingPool, surfaceThicknessFlux, surfaceThicknessFluxRunoff, err)

      !need surface thickness flux too

      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)

      !interpolate edge stresses to cell centers to get X and Y directions
!      call mpas_reconstruct(meshPool, surfaceStress, stressX, stressY, stressZ,  &
!                            stressZonal, stressMeridional, includeHalos = .true.)

      fac = rfac(config_jerlov_water_type)
      dep1 = depth1(config_jerlov_water_type)
      dep2 = depth2(config_jerlov_water_type)

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         les_bottom(iCell) = max(zMid(maxLevelCell(iCell),iCell), config_les_bottom)

         if(config_AM_mixedLayerDepths_enable) then
            les_disturb_bottom(iCell) = -2.0_RKIND*dThreshMLD(iCell)
         else
            les_disturb_bottom(iCell) = config_les_disturb_bottom
         endif

         fracAbsorbed = 1.0_RKIND - exp( max(-layerThickness(1, iCell) / config_flux_attenuation_coefficient, -100.0_RKIND) )
         fracAbsorbedRunoff = 1.0_RKIND - exp( max(-layerThickness(1, iCell) / config_flux_attenuation_coefficient_runoff, &
                              -100.0_RKIND) )

         wt_sfc(iCell) = tracerGroupSurfaceFlux(indexTemperature, iCell) - fracAbsorbed * surfaceThicknessFlux(iCell) * &
                        activeTracers(indexTemperature,1,iCell) + tracerGroupSurfaceFluxRunoff(indexTemperature,iCell) * &
                        fracAbsorbedRunoff

         ws_sfc(iCell) = tracerGroupSurfaceFlux(indexSalinity, iCell) - fracAbsorbed * surfaceThicknessFlux(iCell) *  &
                        activeTracers(indexSalinity,1,iCell)

         uw_sfc(iCell) = stressZonal(iCell) / rho_sw

         if(abs(stressMeridional(iCell)) < 1e-6_RKIND) then
                 stressMeridional(iCell) = 0.0_RKIND
         endif

         vw_sfc(iCell) = stressMeridional(iCell) / rho_sw

         latVal(iCell) = ASIN(fCell(iCell) / (2.0_RKIND*7.292E-5_RKIND))

     enddo
     !$omp end do

      if(first) then

        call palm_init(nCells, nVertLevels, activeTracers(1,:,:), activeTracers(2,:,:),              &
             velocityZonal, velocityMeridional, layerThickness, latVal,                       &
             maxLevelCell, -wt_sfc, penetrativeTemperatureFlux, ws_sfc,                       &
             -uw_sfc, -vw_sfc, fac, dep1, dep2, config_les_dx, config_les_dy, config_les_dz,  &
             config_les_nx, config_les_ny, config_les_nz,                 &
             etime, config_les_dt_disturb, tIncrementLES, sIncrementLES, uIncrementLES,       &
             vIncrementLES, temperatureLES, salinityLES, velocityZonalLES,                    &
             velocityMeridionalLES, dt, zLES, config_les_disturb_energy_limit,                &
             config_les_disturb_amplitude, les_disturb_bottom, config_les_disturb_top,        &
             config_les_disturb_nblocks, les_bottom, config_les_time_avg, constant_dz,        &
             dThreshMLD, mixed_layer_refine,0.0_RKIND)

      endif

      etime = dt

      call palm_main(nCells, nVertLevels, activeTracers(1,:,:), activeTracers(2,:,:),              &
             velocityZonal, velocityMeridional, layerThickness, latVal,                       &
             maxLevelCell, -wt_sfc, penetrativeTemperatureFlux, ws_sfc,                       &
             -uw_sfc, -vw_sfc, fac, dep1, dep2, config_les_dx, config_les_dy, config_les_dz,  &
             config_les_nx, config_les_ny, config_les_nz,                 &
             etime, config_les_dt_disturb, tIncrementLES, sIncrementLES, uIncrementLES,       &
             vIncrementLES, temperatureLES, salinityLES, velocityZonalLES,                    &
             velocityMeridionalLES, dt, zLES, config_les_disturb_energy_limit,                &
             config_les_disturb_amplitude, les_disturb_bottom, config_les_disturb_top,        &
             config_les_disturb_nblocks, les_bottom, config_les_time_avg, constant_dz,        &
             dThreshMLD, mixed_layer_refine, config_les_restoring_strength)

      normalVelocityIncrementLES(:,:) = 0.0_RKIND
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         zonalAverage = 0.5_RKIND * (uIncrementLES(:,cell1) + uIncrementLES(:,cell2))
         meridionalAverage = 0.5_RKIND * (vIncrementLES(:,cell1) + vIncrementLES(:,cell2))

         normalVelocityIncrementLES(:,iEdge) = cos(angleEdge(iEdge)) * zonalAverage  &
                                               + sin(angleEdge(iEdge)) * meridionalAverage

      end do
      !$omp end do

      if(first) first = .false.

      deallocate(uw_sfc, vw_sfc, ws_sfc, wt_sfc)
      deallocate(les_disturb_bottom, les_bottom,latVal)

      deallocate(zonalAverage)
      deallocate(meridionalAverage)

   end subroutine ocn_vmix_LES_interface!}}}

!***********************************************************************
!
!  routine ocn_vmix_LES_apply_tend
!
!> \brief   Updates fields due to tends from LES
!> \author  Luke Van Roekel
!> \date    09 January 2019
!> \details
!>  This routine updates state fields due to LES tendencies.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_LES_apply_tend(dt, meshPool, statePool, forcingPool, diagnosticsPool, timeLevelIn)!{{{

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

      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND) :: dt
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), pointer :: &
              normalVelocity, tIncrementLES, sIncrementLES, normalVelocityIncrementLES

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      integer, pointer :: indexTemperature, indexSalinity

      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdge

      integer :: iCell, iEdge, k, nEdges, nCells

      logical,pointer :: config_use_les_for_vmix, config_use_ecosysTracers

      call mpas_pool_get_config(ocnConfigs, 'config_use_les_for_vmix', config_use_les_for_vmix)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers', config_use_ecosysTracers)
      if ( .not. config_use_les_for_vmix ) return
      if (config_use_ecosysTracers) then
         call mpas_log_write('Error, superparameterization is not configured for use with BGC',MPAS_LOG_CRIT)
      endif

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevelIn)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdge)

      call mpas_pool_get_array(diagnosticsPool, 'tIncrementLES', tIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'sIncrementLES', sIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityIncrementLES', normalVelocityIncrementLES)

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevelIn)

      nEdges = nEdgesArray(1)
      nCells = nCellsArray(1)

      !$omp do schedule(runtime)
      do iCell=1,nCells
         do k=1,maxLevelCell(iCell)
             activeTracers(indexTemperature,k,iCell) = activeTracers(indexTemperature,k,iCell) + &
                                                      dt*tIncrementLES(k,iCell)
            activeTracers(indexSalinity,k,iCell) = activeTracers(indexSalinity,k,iCell) +       &
                                                   dt*sIncrementLES(k,iCell)
         end do
      end do
      !$omp end do


      !$omp do schedule(runtime)
      do iEdge=1,nEdges
         do k=1,maxLevelEdge(iEdge)
            normalVelocity(k,iEdge) = normalVelocity(k,iEdge) + dt*normalVelocityIncrementLES(k,iEdge)
         end do
      end do
      !$omp end do

   end subroutine ocn_vmix_LES_apply_tend!}}}

end module ocn_vmix_LES

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
