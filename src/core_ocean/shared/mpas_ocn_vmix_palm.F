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

   integer, parameter :: num_water_types = 5
   real (kind=RKIND), dimension(num_water_types) ::                       &
      rfac   = (/ 0.58_RKIND, 0.62_RKIND, 0.67_RKIND, 0.77_RKIND, 0.78_RKIND /), &
      depth1 = (/ 0.35_RKIND, 0.60_RKIND, 1.00_RKIND, 1.50_RKIND, 1.40_RKIND /), &
      depth2 = (/ 23.0_RKIND, 20.0_RKIND, 17.0_RKIND, 14.0_RKIND, 7.90_RKIND /)
   logical, public :: run_palm_finalize = .false.

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_vmix_LES_interface, &
             ocn_vmix_LES_tracer_tend, &
             ocn_vmix_LES_vel_tend

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_LES_interface
!
!> \brief   Compute mixing tendency due to LES
!> \author  Luke Van Roekel, Qing Li
!> \date    09 January 2019
!> \details
!>  This routine computes the sub grid mixing tendencies
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_LES_interface(meshPool, statePool, forcingPool, diagnosticsPool, dt, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

      logical, save :: first_palm = .true.
      integer, save :: nzMPAS
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:), pointer :: lesCell, indexToCellID

      real (kind=RKIND), dimension(:), pointer :: &
        stressZonal, stressMeridional, fCell, &
        surfaceStress, surfaceStressMagnitude,   &
        penetrativeTemperatureFlux, angleEdge, surfaceThicknessFlux, &
        surfaceThicknessFluxRunoff, dThreshMLD

      real (kind=RKIND), dimension(:,:), pointer :: &
        layerThickness, velocityZonal, velocityMeridional, &
        normalVelocityCur, normalVelocityNew, &
        uIncrementLES, vIncrementLES, tIncrementLES, &
        sIncrementLES, normalVelocityIncrementLES,     &
        tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
        tracerGroupSurfaceFluxRemoved, temperatureLES, zLES, &
        salinityLES, velocityZonalLES, velocityMeridionalLES, &
        zMid, zTop, &
        velocityX, velocityY, velocityZ

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersCur, activeTracersNew
      integer, pointer :: indexTemperature, indexSalinity, config_jerlov_water_type
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer :: k, i, iCell, iEdge, nCells, nEdges
      integer :: nLevs, edgeCount, cell1, cell2
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      character(len=strKind), pointer :: config_dt

      logical, pointer :: config_use_ecosysTracers,                             &
                          config_use_activeTracers_surface_bulk_forcing,        &
                          config_les_idealized_diurnal

      real (kind=RKIND), pointer :: config_les_dx, config_les_dy,                              &
                                    config_les_spinup_time,                                    &
                                    config_les_bottom, config_les_ls_factor_tracer,            &
                                    config_les_ls_factor_vel,                                  &
                                    config_flux_attenuation_coefficient,                       &
                                    config_flux_attenuation_coefficient_runoff

      integer, pointer :: config_les_nz, config_les_nx, config_les_ny, config_les_ls_mode

      real (kind=RKIND), dimension(:), allocatable :: zonalAverage, meridionalAverage
      real (kind=RKIND), dimension(:), allocatable :: tLSForcing, sLSForcing, uLSForcing, vLSForcing
      real (kind=RKIND), dimension(:,:), allocatable :: uCur, uNew, vCur, vNew

      real (kind=RKIND) :: fac, dep1, dep2, fracAbsorbed, fracAbsorbedRunoff,  &
                           spinuptime, les_bottom, uw_sfc, vw_sfc, wt_sfc, ws_sfc

      err=0

      call mpas_pool_get_config(ocnConfigs, 'config_dt', config_dt)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers', config_use_ecosysTracers)
      call mpas_pool_get_config(ocnConfigs, 'config_use_activeTracers_surface_bulk_forcing', config_use_activeTracers_surface_bulk_forcing)

      if (config_use_ecosysTracers) then
         call mpas_log_write('Error, superparameterization is not configured for use with BGC',MPAS_LOG_CRIT)
      endif

      call mpas_pool_get_config(ocnConfigs, 'config_flux_attenuation_coefficient',  &
                               config_flux_attenuation_coefficient)
      call mpas_pool_get_config(ocnConfigs, 'config_flux_attenuation_coefficient_runoff',  &
                               config_flux_attenuation_coefficient_runoff)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dx', config_les_dx)
      call mpas_pool_get_config(ocnConfigs, 'config_les_dy', config_les_dy)
      call mpas_pool_get_config(ocnConfigs, 'config_les_nx', config_les_nx)
      call mpas_pool_get_config(ocnConfigs, 'config_les_ny', config_les_ny)
      call mpas_pool_get_config(ocnConfigs, 'config_les_nz', config_les_nz)
      call mpas_pool_get_config(ocnConfigs, 'config_les_ls_mode', config_les_ls_mode)
      call mpas_pool_get_config(ocnConfigs, 'config_les_ls_factor_tracer', config_les_ls_factor_tracer)
      call mpas_pool_get_config(ocnConfigs, 'config_les_ls_factor_vel', config_les_ls_factor_vel)
      call mpas_pool_get_config(ocnConfigs, 'config_les_spinup_time', config_les_spinup_time)
      call mpas_pool_get_config(ocnConfigs, 'config_les_bottom', config_les_bottom)
      call mpas_pool_get_config(ocnConfigs, 'config_les_idealized_diurnal', config_les_idealized_diurnal)
      call mpas_pool_get_config(ocnConfigs, 'config_jerlov_water_type', config_jerlov_water_type)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocityNew, 2)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
      call mpas_pool_get_array(meshPool, 'lesCell', lesCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersCur, 1)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersNew, 2)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFlux", tracerGroupSurfaceFlux)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFluxRunoff", tracerGroupSurfaceFluxRunoff)
      call mpas_pool_get_array(tracersSurfaceFluxPool, "activeTracersSurfaceFluxRemoved", tracerGroupSurfaceFluxRemoved)
      call mpas_pool_get_array(forcingPool, "surfaceThicknessFlux", surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, "surfaceThicknessFluxRunoff", surfaceThicknessFluxRunoff)
      call mpas_pool_get_array(meshPool, 'fCell', fCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

      call mpas_pool_get_array(diagnosticsPool, 'tIncrementLES', tIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'sIncrementLES', sIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'uIncrementLES', uIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'vIncrementLES', vIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityIncrementLES', normalVelocityIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
      call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
      call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)

      call mpas_pool_get_array(diagnosticsPool, 'zLES', zLES)
      call mpas_pool_get_array(diagnosticsPool, 'temperatureLES', temperatureLES)
      call mpas_pool_get_array(diagnosticsPool, 'salinityLES', salinityLES)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonalLES', velocityZonalLES)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridionalLES', velocityMeridionalLES)

      ! need to use RBFs to get stress at cell centers
      call mpas_pool_get_array(forcingPool, 'windStressZonal', stressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', stressMeridional)
      call mpas_pool_get_array(forcingPool, 'surfaceStress', surfaceStress)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)
      spinuptime = config_les_spinup_time

      nCells = nCellsArray( 1 )
      nEdges = nEdgesArray( 1 )

      allocate(zonalAverage(nVertLevels),meridionalAverage(nVertLevels))
      allocate(tLSForcing(nVertLevels), sLSForcing(nVertLevels))
      allocate(uLSForcing(nVertLevels), vLSForcing(nVertLevels))
      allocate(uCur(nVertLevels,nCells), vCur(nVertLevels,nCells))
      allocate(uNew(nVertLevels,nCells), vNew(nVertLevels,nCells))

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
         call ocn_surface_bulk_forcing_tracers(meshPool, "activeTracers", forcingPool, activeTracersCur,  &
                                               tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
                                               tracerGroupSurfaceFluxRemoved, dt, layerThickness, err)
      endif

      call ocn_surface_land_ice_fluxes_tracers(meshPool, "activeTracers", forcingPool, tracerGroupSurfaceFlux, err)

      call ocn_surface_bulk_forcing_thick(meshPool, forcingPool, surfaceThicknessFlux, surfaceThicknessFluxRunoff, err)

      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)

      ! reconstruct zonal and meridional velocity
      call mpas_reconstruct(meshPool, normalVelocityCur,  &
                       velocityX, velocityY, velocityZ,   &
                       uCur, vCur, includeHalos = .true.)
      call mpas_reconstruct(meshPool, normalVelocityNew,  &
                       velocityX, velocityY, velocityZ,   &
                       uNew, vNew, includeHalos = .true.)

      ! water type parameters
      fac = rfac(config_jerlov_water_type)
      dep1 = depth1(config_jerlov_water_type)
      dep2 = depth2(config_jerlov_water_type)

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         ! only run LES where lesCell == 1
         if ( lesCell(iCell) == 1 ) then
            ! set surface forcing
            fracAbsorbed = 1.0_RKIND - exp( max(-layerThickness(1, iCell) / config_flux_attenuation_coefficient, -100.0_RKIND) )
            fracAbsorbedRunoff = 1.0_RKIND - exp( max(-layerThickness(1, iCell) / config_flux_attenuation_coefficient_runoff, &
                              -100.0_RKIND) )

            wt_sfc = tracerGroupSurfaceFlux(indexTemperature, iCell) - fracAbsorbed * surfaceThicknessFlux(iCell) * &
                     activeTracersCur(indexTemperature,1,iCell) + tracerGroupSurfaceFluxRunoff(indexTemperature,iCell) * &
                     fracAbsorbedRunoff

            ws_sfc = tracerGroupSurfaceFlux(indexSalinity, iCell) - fracAbsorbed * surfaceThicknessFlux(iCell) *  &
                     activeTracersCur(indexSalinity,1,iCell)

            uw_sfc = stressZonal(iCell) / rho_sw

            vw_sfc = stressMeridional(iCell) / rho_sw

            ! run PALM
            if ( first_palm ) then

               ! find the MPAS level that corresponds to the LES bottom
               les_bottom = max(zMid(maxLevelCell(iCell),iCell), config_les_bottom)
               do k = 1, maxLevelCell(iCell)
                  if ( zMid(k,iCell) - les_bottom <= 0.0_RKIND ) then
                     nzMPAS = k - 1
                     exit
                  endif
               enddo

               ! initial step of PALM, longer initialization
               call palm_init(dt, nzMPAS, zTop(1:nzMPAS+1,iCell),                                 &
                    activeTracersCur(indexTemperature,1:nzMPAS,iCell)+273.15_RKIND,               &
                    activeTracersCur(indexSalinity,1:nzMPAS,iCell),                               &
                    uCur(1:nzMPAS,iCell), vCur(1:nzMPAS,iCell),                                   &
                    -wt_sfc, ws_sfc, -uw_sfc, -vw_sfc, -penetrativeTemperatureFlux(iCell),        &
                    config_les_nz, config_les_nx, config_les_ny, config_les_dx, config_les_dy,    &
                    config_les_idealized_diurnal, fac, dep1, dep2,                                &
                    spinuptime, config_les_ls_mode, config_les_ls_factor_tracer,                  &
                    config_les_ls_factor_vel, fCell(iCell), indexToCellID(iCell))

               ! flag to finalize palm
               run_palm_finalize = .true.

            endif

            tLSForcing = activeTracersNew(indexTemperature,:,iCell) - activeTracersCur(indexTemperature,:,iCell)
            sLSForcing = activeTracersNew(indexSalinity,:,iCell) - activeTracersCur(indexSalinity,:,iCell)
            ! uLSForcing = uNew(:,iCell) - uCur(:,iCell) + fCell(iCell) * vCur(:,iCell) * dt
            ! vLSForcing = vNew(:,iCell) - vCur(:,iCell) - fCell(iCell) * uCur(:,iCell) * dt
            uLSForcing = 0.0_RKIND
            vLSForcing = 0.0_RKIND

            ! continue PALM
            call palm_main(nzMPAS, zTop(1:nzMPAS+1,iCell),                                         &
                 activeTracersCur(indexTemperature,1:nzMPAS,iCell)+273.15_RKIND,                   &
                 activeTracersCur(indexSalinity,1:nzMPAS,iCell),                                   &
                 uCur(1:nzMPAS,iCell), vCur(1:nzMPAS,iCell),                                       &
                 -wt_sfc, ws_sfc, -uw_sfc, -vw_sfc, -penetrativeTemperatureFlux(iCell),            &
                 tLSForcing(1:nzMPAS), sLSForcing(1:nzMPAS),                                       &
                 uLSForcing(1:nzMPAS), vLSForcing(1:nzMPAS),                                       &
                 tIncrementLES(1:nzMPAS,iCell), sIncrementLES(1:nzMPAS,iCell),                     &
                 uIncrementLES(1:nzMPAS,iCell), vIncrementLES(1:nzMPAS,iCell),                     &
                 zLES(:,iCell), temperatureLES(:,iCell), salinityLES(:,iCell),                     &
                 velocityZonalLES(:,iCell), velocityMeridionalLES(:,iCell))

         else
            tIncrementLES(:,iCell) = 0.0_RKIND
            sIncrementLES(:,iCell) = 0.0_RKIND
            uIncrementLES(:,iCell) = 0.0_RKIND
            vIncrementLES(:,iCell) = 0.0_RKIND
            temperatureLES(:,iCell) = 0.0_RKIND
            salinityLES(:,iCell) = 0.0_RKIND
            velocityZonalLES(:,iCell) = 0.0_RKIND
            velocityMeridionalLES(:,iCell) = 0.0_RKIND
            zLES(:,iCell) = 0.0_RKIND
         endif

      enddo
      !$omp end do

      normalVelocityIncrementLES(:,:) = 0.0_RKIND
      !$omp do schedule(runtime) private(cell1, cell2)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)

         zonalAverage = 0.5_RKIND * (uIncrementLES(:,cell1) + uIncrementLES(:,cell2))
         meridionalAverage = 0.5_RKIND * (vIncrementLES(:,cell1) + vIncrementLES(:,cell2))

         normalVelocityIncrementLES(:,iEdge) = cos(angleEdge(iEdge)) * zonalAverage  &
                                               + sin(angleEdge(iEdge)) * meridionalAverage

      end do
      !$omp end do

      if ( first_palm ) first_palm = .false.

      deallocate(zonalAverage)
      deallocate(meridionalAverage)
      deallocate(tLSForcing, sLSForcing, uLSForcing, vLSForcing)
      deallocate(uCur, uNew, vCur, vNew)

   end subroutine ocn_vmix_LES_interface!}}}

!***********************************************************************
!
!  routine ocn_vmix_LES_tracer_tend
!
!> \brief   Updates fields due to tends from LES for tracers
!> \author  Luke Van Roekel, Qing Li
!> \date    09 January 2019
!> \details
!>  This routine updates tracer state fields due to LES tendencies.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_LES_tracer_tend(dt, meshPool, statePool, diagnosticsPool, timeLevelIn)!{{{

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

      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND) :: dt
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), pointer :: &
              tIncrementLES, sIncrementLES

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      integer, pointer :: indexTemperature, indexSalinity

      integer, dimension(:), pointer :: nCellsArray

      integer, dimension(:), pointer :: maxLevelCell

      integer :: iCell, k, nCells

      logical,pointer :: config_use_ecosysTracers, config_run_les

      call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers', config_use_ecosysTracers)
      if (config_use_ecosysTracers) then
         call mpas_log_write('Error, superparameterization is not configured for use with BGC',MPAS_LOG_CRIT)
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_run_les', config_run_les)
      if (.not. config_run_les) then
         call mpas_log_write('Error, config_run_les = .true. required if config_use_les_tracer_vmix = .true.',MPAS_LOG_CRIT)
      end if

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevelIn)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call mpas_pool_get_array(diagnosticsPool, 'tIncrementLES', tIncrementLES)
      call mpas_pool_get_array(diagnosticsPool, 'sIncrementLES', sIncrementLES)

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

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

   end subroutine ocn_vmix_LES_tracer_tend!}}}


!***********************************************************************
!
!  routine ocn_vmix_LES_vel_tend
!
!> \brief   Updates fields due to tends from LES for the velocity
!> \author  Luke Van Roekel, Qing Li
!> \date    09 January 2019
!> \details
!>  This routine updates the velocity due to LES tendencies.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_LES_vel_tend(dt, meshPool, statePool, diagnosticsPool, timeLevelIn)!{{{

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

      real (kind=RKIND) :: dt
      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), pointer :: &
              normalVelocity, normalVelocityIncrementLES

      integer, dimension(:), pointer :: nEdgesArray

      integer, dimension(:), pointer :: maxLevelEdge

      integer :: iEdge, k, nEdges

      logical,pointer :: config_run_les

      call mpas_pool_get_config(ocnConfigs, 'config_run_les', config_run_les)
      if (.not. config_run_les) then
         call mpas_log_write('Error, config_run_les = .true. required if config_use_les_tracer_vmix = .true.',MPAS_LOG_CRIT)
      end if

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdge)

      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityIncrementLES', normalVelocityIncrementLES)

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevelIn)

      nEdges = nEdgesArray(1)

      !$omp do schedule(runtime)
      do iEdge=1,nEdges
         do k=1,maxLevelEdge(iEdge)
            normalVelocity(k,iEdge) = normalVelocity(k,iEdge) + dt*normalVelocityIncrementLES(k,iEdge)
         end do
      end do
      !$omp end do

   end subroutine ocn_vmix_LES_vel_tend!}}}

end module ocn_vmix_LES

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
