! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_lateral_boundary_flux
!
!> \brief MPAS ocean lateral flux for thickness
!> \author Qing Li
!> \date   15 April 2021
!> \details
!>  This module contains the routine for computing
!>  tendencies for thickness from lateral flux at boundaries
!
!-----------------------------------------------------------------------

module ocn_thick_lateral_boundary_flux

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

   public :: ocn_thick_lateral_boundary_flux_tend, &
             ocn_thick_lateral_boundary_flux_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: thickBoundaryForcingOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_thick_lateral_boundary_flux_tend
!
!> \brief   Computes tendency term from lateral boundary flux of thickness
!> \author  Qing Li
!> \date    15 April 2021
!> \details
!>  This routine computes the thickness tendency from lateral
!>  boundary fluxes
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_lateral_boundary_flux_tend(meshPool, forcingPool, tend, err)!{{{

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
         tend          !< Input/Output: velocity tendency

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
      real (kind=RKIND), dimension(:,:), pointer :: boundaryThickTend

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if (.not.thickBoundaryForcingOn) return

      call mpas_timer_start("thick lateral boundary forcing")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)

      call mpas_pool_get_array(forcingPool, 'boundaryThickTend', boundaryThickTend)

      nCells = nCellsArray( 1 )

      !$omp parallel
      !$omp do schedule(runtime) private(k, iCell)
      do iCell = 1, nCells
        do k = minLevelCell(iCell), maxLevelCell(iCell)
           tend(k, iCell) = tend(k, iCell) + boundaryThickTend(k, iCell)
        end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("thick lateral boundary forcing")

   !--------------------------------------------------------------------

   end subroutine ocn_thick_lateral_boundary_flux_tend!}}}

!***********************************************************************
!
!  routine ocn_thick_lateral_boundary_flux_init
!
!> \brief   Initializes ocean lateral boundary thickness forcing
!> \author  Qing Li
!> \date    15 April 2021
!> \details
!>  This routine initializes quantities related to lateral boundary
!>  thickness forcing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_lateral_boundary_flux_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      thickBoundaryForcingOn = .false.

      if (config_use_lateral_boundary_forcing_thickness) thickBoundaryForcingOn = .true.

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_thick_lateral_boundary_flux_init!}}}

!***********************************************************************

end module ocn_thick_lateral_boundary_flux

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
