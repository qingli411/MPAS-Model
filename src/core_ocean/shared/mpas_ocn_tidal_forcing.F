! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tidal_forcing
!
!> \brief MPAS ocean tidal forcing module
!> \author Phillip J. Wolfram
!> \date   04/05/2019
!> \details
!>  This module contains routines for the tidal forcing.
!
!-----------------------------------------------------------------------

module ocn_tidal_forcing

   use mpas_kind_types
   use mpas_constants
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
   use mpas_log
   use ocn_constants
   use ocn_config
   use ocn_equation_of_state
   use ocn_diagnostics_variables

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

   public :: ocn_tidal_forcing_build_array, &
             ocn_tidal_forcing_layer_thickness, &
             ocn_tidal_forcing_tracers, &
             ocn_tidal_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: tidalFluxOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tidal_forcing_layer_thickness
!
!> \brief   Add tendency due to tidal forcing
!> \author  Phillip J. Wolfram
!> \date    04/05/2019
!> \details
!>  This routine adds a tendency to layer thickness due to tidal forcing
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_forcing_layer_thickness(meshPool, forcingPool, layerThicknessTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(inout), dimension(:,:) :: layerThicknessTend

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k, nCells
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), pointer :: minLevelCell, maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: tidalLayerThicknessTendency

      err = 0

      if ( .not. tidalFluxOn ) return

      call mpas_timer_start("tidal thickness tendency")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(forcingPool, 'tidalLayerThicknessTendency', &
                                             tidalLayerThicknessTendency)

      ! frazil fields are needed only over 0 and 1 halos
      nCells = nCellsArray( 2 )

      ! Build surface fluxes at cell centers
      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells
        do k = minLevelCell(iCell), maxLevelCell(iCell)
          layerThicknessTend(k,iCell) = layerThicknessTend(k,iCell) + &
                                        tidalLayerThicknessTendency(k,iCell)

        end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("tidal thickness tendency")

   end subroutine ocn_tidal_forcing_layer_thickness !}}}

!***********************************************************************
!
!  routine ocn_tidal_forcing_tracers
!
!> \brief   Add the tracer tendency due to tidal forcing
!> \author  Qing Li
!> \date    27 April 2021
!> \details
!>
!>  This routine adds to the tracer tendency arrays the tidal forcing
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_forcing_tracers(meshPool, tracersPool, groupName, forcingPool, tracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      character (len=*) :: &
         groupName         !< Input: Name of tracer group
      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      type (mpas_pool_type), intent(in) :: &
         tracersPool       !< Input: tracer tendency pool
      type (mpas_pool_type), intent(in) :: &
         forcingPool       !< Input: forcing pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tracersTend   !< Input/Output: tracers tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

      if ( .not. tidalFluxOn ) return

      call mpas_timer_start("tidal tracer tendency")

      if ( trim(groupName) == 'activeTracers' ) then
         call ocn_tidal_forcing_active_tracers(meshPool, tracersPool, forcingPool, tracersTend, err)
      end if

      call mpas_timer_stop("tidal tracer tendency")

   end subroutine ocn_tidal_forcing_tracers!}}}

!***********************************************************************
!
!  routine ocn_tidal_forcing_active_tracers
!
!> \brief   Add the active tracer tendency due to tidal forcing
!> \author  Qing Li
!> \date    27 April 2021
!> \details
!>  This routine adds the active tracer tendency due to tidal forcing.
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_forcing_active_tracers(meshPool, tracersPool, forcingPool, activeTracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      type (mpas_pool_type), intent(in) :: &
         tracersPool       !< Input:
      type (mpas_pool_type), intent(in) :: &
         forcingPool       !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         activeTracersTend          !< Input/Output: tracers tendency

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

      integer :: k, iCell, nCells
      real (kind=RKIND) :: inflow, outflow
      integer, pointer :: indexTemperature
      integer, pointer :: indexSalinity
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), pointer :: maxLevelCell, minLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: tidalLayerThicknessTendency
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if (.not. tidalFluxOn) return

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)

      call mpas_pool_get_array(forcingPool, 'tidalLayerThicknessTendency', &
                                             tidalLayerThicknessTendency)

      nCells = nCellsArray( 2 )

      !$omp parallel
      !$omp do schedule(runtime) private(k, iCell, inflow, outflow)
      do iCell = 1, nCells
        do k = minLevelCell(iCell), maxLevelCell(iCell)
           inflow  = max(0.0_RKIND, tidalLayerThicknessTendency(k,iCell))
           outflow = min(0.0_RKIND, tidalLayerThicknessTendency(k,iCell))
           activeTracersTend(indexTemperature,k,iCell) = activeTracersTend(indexTemperature,k,iCell) &
               + inflow * config_tidal_forcing_inflow_temperature &
               + outflow * activeTracers(indexTemperature,k,iCell)
           activeTracersTend(indexSalinity,k,iCell) = activeTracersTend(indexSalinity,k,iCell) &
               + inflow * config_tidal_forcing_inflow_salinity &
               + outflow * activeTracers(indexSalinity,k,iCell)
        end do
      end do
      !$omp end do
      !$omp end parallel

   !--------------------------------------------------------------------

   end subroutine ocn_tidal_forcing_active_tracers!}}}

!***********************************************************************
!
!  routine ocn_tidal_forcing_build_array
!
!> \brief   Computes tidal forcing functions.
!> \author  Phillip J. Wolfram
!> \date    04/05/2019
!> \details
!>   ocn_tidal_forcing_build_array computes the tendencies to layer thickness
!>   due to tidal forcing of layer thickness
!>
!>   this routine should be call at the beginning of whatever time stepping method is utilized
!>      and the tendencies should be retieved when building up the RHS of the thickess, temperature
!>      and salinity equations.
!>
!>   this routine is only applicable to the thickness fields
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_forcing_build_array(domain, meshPool, forcingPool, statePool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer, intent(inout) :: statePool !< Input: State information
      type (mpas_pool_type), pointer, intent(inout) :: forcingPool !< Input: Forcing information
      integer, intent(inout) :: err !< Error flag

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), pointer :: tidalLayerThicknessTendency
      real (kind=RKIND), dimension(:), pointer :: tidalInputMask
      real (kind=RKIND), dimension(:), pointer :: tidalBCValue
      real (kind=RKIND), dimension(:), pointer :: tidalInputFactor
      real (kind=RKIND), dimension(:), pointer :: bottomDepth

      integer :: iCell, k, nCells
      integer, dimension(:), pointer :: nCellsArray
      integer, pointer, dimension(:) :: minLevelCell, maxLevelCell
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: dt

      type (MPAS_time_type) :: currentTime
      real (kind=RKIND), pointer, dimension(:)     :: ssh
      real (kind=RKIND), pointer, dimension(:,:)   :: zMid
      real (kind=RKIND), pointer, dimension(:,:)   :: layerThickness

      ! local variables
      real (kind=RKIND) :: totalDepth, tidalHeight

      character (len=StrKIND), pointer :: simulationStartTime, xtime
      type (MPAS_Time_type) :: startTime, xtime_timeType, simulationStartTime_timeType

      if ( .not. tidalFluxOn ) return

      call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(forcingPool, 'tidalLayerThicknessTendency', &
                                             tidalLayerThicknessTendency)
      call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
      call mpas_pool_get_array(forcingPool, 'tidalBCValue', tidalBCValue)
      call mpas_pool_get_array(forcingPool, 'tidalInputFactor', tidalInputFactor)

      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      ! tidal fields are needed only over 0 and 1 halos
      nCells = nCellsArray( 2 )

      if (trim(config_tidal_forcing_type) == 'boundary_flux') then

        ! loop over all columns
        do iCell = 1, nCells
          ! compute total depth for relative thickness contribution
          totalDepth = 0.0_RKIND
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            totalDepth = totalDepth + layerThickness(k,iCell)
          end do

          ! compute layer thickness
          tidalLayerThicknessTendency(:,iCell) = 0.0_RKIND
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            ! compute the tidal layer thickness tendency from the tidal velocity
            tidalLayerThicknessTendency(k,iCell) = tidalInputFactor(iCell) &
               * config_tidal_forcing_monochromatic_velocity &
               * (layerThickness(k,iCell)/totalDepth) &
               * sin(2.0_RKIND * pii / config_tidal_forcing_monochromatic_period &
               * daysSinceStartOfSim - pii * config_tidal_forcing_monochromatic_phaseLag / 180.0_RKIND)
          end do
        end do   ! do iCell = 1, nCells

      else

        ! loop over all columns
        do iCell = 1, nCells
          ! compute the tidalHeight
          tidalHeight = config_tidal_forcing_monochromatic_amp * &
            sin(2.0_RKIND*pii/config_tidal_forcing_monochromatic_period * daysSinceStartOfSim - &
            pii*config_tidal_forcing_monochromatic_phaseLag/180.0_RKIND) - &
            config_tidal_forcing_monochromatic_baseline

          ! ensure tidalHeight can't push flow below anticipated minimum
          if (config_use_wetting_drying .and. tidalInputMask(iCell) == 1.0_RKIND) then
            ! ensure that tidal height can't force below total minimum thickness
            ! condition wrong to ensure that there isn't any drying according to criteria
            tidalHeight = max(-bottomDepth(iCell) + (real(maxLevelCell(iCell), RKIND)+1.0_RKIND)*config_drying_min_cell_height, tidalHeight)
          end if

          ! compute total depth for relative thickness contribution
          totalDepth = 0.0_RKIND
          do k = minLevelCell(iCell), maxLevelCell(iCell)
            totalDepth = totalDepth + layerThickness(k,iCell)
          end do

          tidalLayerThicknessTendency(:,iCell) = 0.0_RKIND
          if (trim(config_tidal_forcing_type) == 'thickness_source') then
            ! distribute tidal forcing tendency fractionally over water column
            do k = minLevelCell(iCell), maxLevelCell(iCell)
              tidalLayerThicknessTendency(k,iCell) = tidalInputMask(iCell) / config_use_tidal_forcing_tau &
                * (layerThickness(k,iCell)/totalDepth) * (tidalHeight - ssh(iCell))
            end do
            ! else if
            ! using the 'direct' option to immediately force the free surface as a prescribe boundary
            ! in src/core_ocean/mode_forward/mpas_ocn_time_integration_rk4.F
          end if

          ! store tidal boundary cell values (e.g., for config_tidal_forcing_type == 'direct')
          tidalBCValue(iCell) = tidalHeight

        end do   ! do iCell = 1, nCells

     end if

   end subroutine ocn_tidal_forcing_build_array !}}}


!***********************************************************************
!
!  routine ocn_tidal_forcing_init
!
!> \brief   Initializes ocean tidal forcing module.
!> \author  Phillip J. Wolfram
!> \date    04/05/2019
!> \details
!>  This routine initializes the ocean tidal forcing module and variables.
!
!-----------------------------------------------------------------------

   subroutine ocn_tidal_forcing_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

      tidalFluxOn = .false.

      if(config_use_tidal_forcing) then
        tidalFluxOn = .true.
        if (trim(config_tidal_forcing_model) .ne. 'monochromatic') then
          call mpas_log_write( &
             "Unknown value for config_tidal_forcing_model, supported values are:" // &
             "         monochromatic", &
             MPAS_LOG_CRIT)
          err = 1
          return
        end if
        if (trim(config_tidal_forcing_type) .ne. 'direct' .and. &
            trim(config_tidal_forcing_type) .ne. 'thickness_source' .and. &
            trim(config_tidal_forcing_type) .ne. 'boundary_flux') then
          call mpas_log_write( &
             "Unknown value for config_tidal_forcing_type, supported values are:" // &
             "         direct, thickness_source, boundary_flux", &
             MPAS_LOG_CRIT)
          err = 1
          return
        end if
      end if

   end subroutine ocn_tidal_forcing_init!}}}


!***********************************************************************

end module ocn_tidal_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
