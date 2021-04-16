! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_idealized_estuary
!
!> \brief MPAS ocean initialize case -- zonally periodic idealized estuary
!> \author Qing Li
!> \date   09/17/2020
!> \details
!>  This module contains the routines for initializing the
!>  idealized estuary initial condition
!
!-----------------------------------------------------------------------

module ocn_init_idealized_estuary

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

   use ocn_constants
   use ocn_config
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

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

   public :: ocn_init_setup_idealized_estuary, &
             ocn_init_validate_idealized_estuary

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_idealized_estuary
!
!> \brief   Setup for this initial condition
!> \author  Qing Li
!> \date    09/17/2020
!> \details
!>  This routine sets up the initial conditions for this case.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_idealized_estuary(domain, iErr)!{{{

    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: forcingPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: tracersInteriorRestoringFieldsPool

    ! local variables
    integer :: iCell, k, idx
    real (kind=RKIND) :: yMin, yMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal
    real (kind=RKIND) :: salinity, invAreaCell

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell, minLevelCell
    real (kind=RKIND), dimension(:), pointer :: tidalInputMask
    real (kind=RKIND), dimension(:), pointer :: yCell,refBottomDepth, refZMid, &
                                                vertCoordMovementWeights, bottomDepth, &
                                                fCell, fEdge, fVertex, dcEdge, areaCell
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:), pointer :: boundaryThickTend
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersInteriorRestoringValue, &
                                                    activeTracersInteriorRestoringRate
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

    ! Define local interfaceLocations variable
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations

    iErr = 0

    if(config_init_configuration .ne. trim('idealized_estuary')) return

    ! check configuration parameters
    if (config_idealized_estuary_depth_ocean .lt. config_idealized_estuary_depth_river) then
       call mpas_log_write('Ocean boundary must be deeper than River boundary!', MPAS_LOG_CRIT)
    end if

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! you may restrict your case geometry as follows:
    if ( on_a_sphere ) call mpas_log_write('The idealized_estuary configuration can only be applied ' &
        // 'to a planar mesh. Exiting...', MPAS_LOG_CRIT)

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

    ! Initalize min/max values to large positive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

    ! Determine local min and max values.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

       yMin = min( yMin, minval(yCell(1:nCellsSolve)))
       yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
       dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))

       block_ptr => block_ptr % next
    end do

    ! Determine global min and max values.
    call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

    yMidGlobal = ( yMinGlobal + yMaxGlobal ) * 0.5_RKIND

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
       call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
       call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
       call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
       call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
       call mpas_pool_get_array(meshPool, 'fCell', fCell)
       call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
       call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

       call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

       call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
       call mpas_pool_get_array(forcingPool, 'boundaryThickTend', boundaryThickTend)

       call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
       call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

       call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, &
                 'activeTracersInteriorRestoringRate', activeTracersInteriorRestoringRate, 1)
       call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, &
                 'activeTracersInteriorRestoringValue', activeTracersInteriorRestoringValue, 1)

       call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
       call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

       tidalInputMask(:) = 0.0_RKIND
       boundaryThickTend(:,:) = 0.0_RKIND
       activeTracersInteriorRestoringRate(:,:,:) = 0.0_RKIND
       activeTracersInteriorRestoringValue(:,:,:) = 0.0_RKIND

       ! Set refBottomDepth and refZMid
       do k = 1, nVertLevels
          refBottomDepth(k) = config_idealized_estuary_depth_ocean * interfaceLocations(k+1)
          refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_idealized_estuary_depth_ocean
       end do

       ! Set vertCoordMovementWeights
       vertCoordMovementWeights(:) = 1.0_RKIND

       do iCell = 1, nCellsSolve

          ! Set temperature
          idx = index_temperature
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = config_idealized_estuary_initial_temperature
          end do

          ! Set salinity
          salinity = config_idealized_estuary_initial_salinity_river &
                   + 0.5_RKIND * ( tanh( 2.0_RKIND * (yCell(iCell) - yMidGlobal) / &
                       config_idealized_estuary_initial_salinity_width ) + 1.0_RKIND ) &
                   * ( config_idealized_estuary_initial_salinity_ocean - &
                       config_idealized_estuary_initial_salinity_river )
          idx = index_salinity
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = salinity
          end do

          ! Set layerThickness and restingThickness
          do k = 1, nVertLevels
              layerThickness(k, iCell) = config_idealized_estuary_depth_ocean * ( interfaceLocations(k+1) &
                                       - interfaceLocations(k) )
              restingThickness(k, iCell) = config_idealized_estuary_depth_ocean * ( interfaceLocations(k+1) &
                                         - interfaceLocations(k) )
          end do

          ! Set bottomDepth
          bottomDepth(iCell) = config_idealized_estuary_depth_river &
             + (yCell(iCell) - yMinGlobal) / (yMaxGlobal - yMinGlobal) &
             * (config_idealized_estuary_depth_ocean - config_idealized_estuary_depth_river)

          ! Set maxLevelCell
          maxLevelCell(iCell) = nVertLevels
          if (nVertLevels > 1) then
             do k = 2, nVertLevels
                if (bottomDepth(iCell) <= -refZMid(k) .and. bottomDepth(iCell) > refBottomDepth(k-1)) then
                   maxLevelCell(iCell) = k-1
                   bottomDepth(iCell) = refBottomDepth(k-1)
                   exit
                end if
                if (bottomDepth(iCell) > -refZMid(k) .and. bottomDepth(iCell) < refBottomDepth(k)) then
                   maxLevelCell(iCell) = k
                   bottomDepth(iCell) = refBottomDepth(k)
                   exit
                end if
             end do
          end if

          ! Set minLevelCell
          minLevelCell(iCell) = 1

          ! Set tidal boundary mask
          if (yMaxGlobal-yCell(iCell) <= 3.0_RKIND*config_idealized_estuary_tidal_boundary_width) then
             ! tidal forcing exponentially decay away from the boundary at the ocean side
             tidalInputMask(iCell) = exp( (yCell(iCell) - yMaxGlobal) / &
                                           config_idealized_estuary_tidal_boundary_width )
          end if

          ! Set restoring at the ocean side boundary
          if (yMaxGlobal-yCell(iCell) <= 3.0_RKIND*config_idealized_estuary_restoring_salinity_width_ocean) then
             do k = minLevelCell(iCell), maxLevelCell(iCell)
                ! no restoring for temperature
                idx = index_temperature
                activeTracersInteriorRestoringValue(idx, k, iCell) = config_idealized_estuary_initial_temperature
                activeTracersInteriorRestoringRate(idx, k, iCell) = 0.0_RKIND
                idx = index_salinity
                activeTracersInteriorRestoringValue(idx, k, iCell) = config_idealized_estuary_initial_salinity_ocean
                activeTracersInteriorRestoringRate(idx, k, iCell) = &
                                 exp( (yCell(iCell) - yMaxGlobal) / &
                                 config_idealized_estuary_restoring_salinity_width_ocean ) * &
                                 ( 1.0_RKIND / (config_idealized_estuary_restoring_salinity_tau*86400.0_RKIND))
             end do
          end if

          ! Set restoring at the river side boundary
          if (yCell(iCell)-yMinGlobal <= 3.0_RKIND*config_idealized_estuary_restoring_salinity_width_river) then
             do k = minLevelCell(iCell), maxLevelCell(iCell)
                ! no restoring for temperature
                idx = index_temperature
                activeTracersInteriorRestoringValue(idx, k, iCell) = config_idealized_estuary_initial_temperature
                activeTracersInteriorRestoringRate(idx, k, iCell) = 0.0_RKIND
                idx = index_salinity
                activeTracersInteriorRestoringValue(idx, k, iCell) = config_idealized_estuary_initial_salinity_river
                activeTracersInteriorRestoringRate(idx, k, iCell) = &
                                 exp( (yMidGlobal - yCell(iCell)) / &
                                 config_idealized_estuary_restoring_salinity_width_river ) * &
                                 ( 1.0_RKIND / (config_idealized_estuary_restoring_salinity_tau*86400.0_RKIND))
             end do
          end if

          ! Set inflow at river side boundary
          if (yCell(iCell)-yMinGlobal <= 0.5_RKIND*dcEdgeMinGlobal) then
             invAreaCell = 1.0_RKIND / areaCell(iCell)
             do k = minLevelCell(iCell), maxLevelCell(iCell)
               boundaryThickTend(k, iCell) = config_idealized_estuary_velocity_river * &
                                             dcEdgeMinGlobal * layerThickness(k, iCell) * invAreaCell
             end do
          end if

          ! Set outflow at ocean side boundary
          if (yMaxGlobal-yCell(iCell) <= 0.5_RKIND*dcEdgeMinGlobal) then
             invAreaCell = 1.0_RKIND / areaCell(iCell)
             do k = minLevelCell(iCell), maxLevelCell(iCell)
               boundaryThickTend(k, iCell) = - config_idealized_estuary_velocity_river * &
                                             config_idealized_estuary_depth_river / &
                                             config_idealized_estuary_depth_ocean * &
                                             dcEdgeMinGlobal * layerThickness(k, iCell) * invAreaCell
             end do
          end if

          ! Set Coriolis parameters, if other than zero
          fCell(iCell) = config_idealized_estuary_coriolis_parameter
          fEdge(iCell) = config_idealized_estuary_coriolis_parameter
          fVertex(iCell) = config_idealized_estuary_coriolis_parameter

       end do

       block_ptr => block_ptr % next
    end do

    deallocate(interfaceLocations)
    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_idealized_estuary!}}}

!***********************************************************************
!
!  routine ocn_init_validate_idealized_estuary
!
!> \brief   Validation for this initial condition
!> \author  Qing Li
!> \date    09/17/2020
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_idealized_estuary(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_idealized_estuary_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('idealized_estuary')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_idealized_estuary_vert_levels', config_idealized_estuary_vert_levels)

      if(config_vert_levels <= 0 .and. config_idealized_estuary_vert_levels > 0) then
         config_vert_levels = config_idealized_estuary_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for idealized_estuary. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_idealized_estuary!}}}


!***********************************************************************

end module ocn_init_idealized_estuary

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
