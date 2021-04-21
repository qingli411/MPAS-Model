! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_lateral_boundary_forcing
!
!> \brief MPAS ocean lateral flux for thickness
!> \author Qing Li
!> \date   15 April 2021
!> \details
!>  This module contains the routine for computing
!>  tendencies for thickness from lateral flux at boundaries
!
!-----------------------------------------------------------------------

module ocn_lateral_boundary_forcing

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
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

   public :: ocn_lateral_boundary_forcing_layer_thickness, &
             ocn_lateral_boundary_forcing_tracers, &
             ocn_lateral_boundary_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: lateralBoundaryForcingOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_lateral_boundary_forcing_tracers
!
!> \brief   Add the tracer tendency due to lateral boundary forcing
!> \author  Qing Li
!> \date    19 April 2021
!> \details
!>
!>  This routine adds to the tracer tendency arrays the lateral
!>  boundary forcing
!
!-----------------------------------------------------------------------

   subroutine ocn_lateral_boundary_forcing_tracers(meshPool, tracersPool, groupName, forcingPool, tracersTend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      character (len=*) :: groupName !< Input: Name of tracer group

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: meshPool !< Input/Output: mesh information
      type (mpas_pool_type), intent(inout) :: tracersPool !< Input/Output: tracer tendency pool
      type (mpas_pool_type), intent(inout) :: forcingPool  !< Input/Output: forcing pool holding lateral boundary forcing tendencies
      real (kind=RKIND), dimension(:,:,:), intent(inout) ::  tracersTend
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

      if ( .not. lateralBoundaryForcingOn ) return

      call mpas_timer_start("lateral boundary forcing tracers")

      if ( trim(groupName) == 'activeTracers' ) then
         call ocn_lateral_boundary_forcing_active_tracers(meshPool, tracersPool, forcingPool, tracersTend, err)
      end if

      call mpas_timer_stop("lateral boundary forcing tracers")

   end subroutine ocn_lateral_boundary_forcing_tracers!}}}

!***********************************************************************
!
!  routine ocn_lateral_boundary_forcing_layer_thickness
!
!> \brief   Add layer thickness tendenty due to lateral boundary forcing
!> \author  Qing Li
!> \date    15 April 2021
!> \details
!>  This routine adds a tendency to layer thickness due to lateral
!>  boundary forcing
!
!-----------------------------------------------------------------------

   subroutine ocn_lateral_boundary_forcing_layer_thickness(meshPool, forcingPool, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information
      type (mpas_pool_type), intent(inout) :: &
         forcingPool       !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: layer thickness tendency

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
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), pointer :: maxLevelCell, minLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: boundaryLayerThicknessTendency

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if (.not.lateralBoundaryForcingOn) return

      call mpas_timer_start("lateral boundary forcing layer thickness")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)

      call mpas_pool_get_array(forcingPool, 'boundaryLayerThicknessTendency', boundaryLayerThicknessTendency)

      nCells = nCellsArray( 2 )

      !$omp parallel
      !$omp do schedule(runtime) private(k, iCell)
      do iCell = 1, nCells
        do k = minLevelCell(iCell), maxLevelCell(iCell)
           tend(k, iCell) = tend(k, iCell) + boundaryLayerThicknessTendency(k, iCell)
        end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("lateral boundary forcing layer thickness")

   !--------------------------------------------------------------------

   end subroutine ocn_lateral_boundary_forcing_layer_thickness!}}}

!***********************************************************************
!
!  routine ocn_lateral_boundary_forcing_active_tracers
!
!> \brief   Add the active tracers forcing due to lateral boundary forcing
!> \author  Qing Li
!> \date    19 April 2021
!> \details
!>  This routine adds the active tracers forcing due to lateral
!>  boundary forcing.
!
!-----------------------------------------------------------------------

   subroutine ocn_lateral_boundary_forcing_active_tracers(meshPool, tracersPool, forcingPool, activeTracersTend, err)!{{{

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
      integer, pointer :: indexTemperature
      integer, pointer :: indexSalinity
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), pointer :: maxLevelCell, minLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: boundaryTemperatureTendency, &
                                                    boundarySalinityTendency

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if (.not.lateralBoundaryForcingOn) return

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)

      call mpas_pool_get_array(forcingPool, 'boundaryTemperatureTendency', boundaryTemperatureTendency)
      call mpas_pool_get_array(forcingPool, 'boundarySalinityTendency', boundarySalinityTendency)

      nCells = nCellsArray( 2 )

      !$omp parallel
      !$omp do schedule(runtime) private(k, iCell)
      do iCell = 1, nCells
        do k = minLevelCell(iCell), maxLevelCell(iCell)
           activeTracersTend(indexTemperature,k,iCell) = activeTracersTend(indexTemperature,k,iCell) &
                                                       + boundaryTemperatureTendency(k,iCell)
           activeTracersTend(indexSalinity,k,iCell) = activeTracersTend(indexSalinity,k,iCell) &
                                                       + boundarySalinityTendency(k,iCell)
        end do
      end do
      !$omp end do
      !$omp end parallel

   !--------------------------------------------------------------------

   end subroutine ocn_lateral_boundary_forcing_active_tracers!}}}

!***********************************************************************
!
!  routine ocn_lateral_boundary_forcing_init
!
!> \brief   Initializes ocean lateral boundary thickness forcing
!> \author  Qing Li
!> \date    15 April 2021
!> \details
!>  This routine initializes quantities related to lateral boundary
!>  thickness forcing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_lateral_boundary_forcing_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      lateralBoundaryForcingOn = .false.

      if (config_use_lateral_boundary_forcing) lateralBoundaryForcingOn = .true.

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_lateral_boundary_forcing_init!}}}

!***********************************************************************

end module ocn_lateral_boundary_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
