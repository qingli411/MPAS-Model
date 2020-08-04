!> @file init_grid.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2018 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
!
!
! Former revisions:
! -----------------
! $Id: init_grid.f90 3068 2018-06-12 14:49:41Z Giersch $
! New warning message concerning grid stretching has been introduced
!
! 3066 2018-06-12 08:55:55Z Giersch
! Bugfix in IF statement before error message
!
! 3065 2018-06-12 07:03:02Z Giersch
! New vertical stretching mechanism introduced
!
! 3051 2018-05-30 17:43:55Z suehring
! Minor bugfix concerning mapping 3D buildings on top of terrain
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
!
! 2968 2018-04-13 11:52:24Z suehring
! Bugfix in initialization in case of elevated model surface. Introduce
! index for minimum topography-top.
!
! 2955 2018-04-09 15:14:01Z suehring
! Improve topography filter routine and add ghost-point exchange for building
! ID and building type.
!
! 2927 2018-03-23 15:13:00Z suehring
! Bugfix, setting boundary conditions for topography index array.
!
! 2918 2018-03-21 15:52:14Z gronemeier
! Moved init_mixing_length to turbulence_closure_mod.f90
!
! 2897 2018-03-15 11:47:16Z suehring
! Relax restrictions for topography input, terrain and building heights can be
! input separately and are not mandatory any more.
!
! 2893 2018-03-14 16:20:52Z suehring
! Revise informative message concerning filtered topography (1 grid-point
! holes).
!
! 2892 2018-03-14 15:06:29Z suehring
! Bugfix, uninitialized array in case of single_building.
!
! 2867 2018-03-09 09:40:23Z suehring
! Revise mapping of 3D buildings onto onto orography.
!
! 2823 2018-02-20 15:31:45Z Giersch
! Set boundary conditions for 3D topography in case of non-cyclic boundary
! conditions
!
! 2796 2018-02-08 12:25:39Z suehring
! Bugfix in 3D building initialization
!
! 2747 2018-01-15 12:44:17Z suehring
! Bugfix, topography height is rounded to the nearest discrete grid level
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2701 2017-12-15 15:40:50Z suehring
! Changes from last commit documented
!
! 2698 2017-12-14 18:46:24Z suehring
! Bugfix in get_topography_top_index
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Revised topography input
! Set nzb_max not for the entire nest domain, only for boundary PEs
! Re-organize routine, split-up into several subroutines
! Modularize poismg_noopt
! Remove setting bit 26, 27, 28 in wall_flags_0, indicating former '_outer'
! arrays (not required any more).
! Bugfix in generic tunnel setup (MS)
!
! 2550 2017-10-16 17:12:01Z boeske
! Set lateral boundary conditions for topography on all three ghost layers
!
! 2478 2017-09-18 13:37:24Z suehring
! Bugfix, correct flag for use_top
!
! 2365 2017-08-21 14:59:59Z kanani
! Vertical nesting implemented (SadiqHuq)
!
! 2319 2017-07-20 17:33:17Z suehring
! Remove print statements
!
! 2318 2017-07-20 17:27:44Z suehring
! Get topography top index via Function call
!
! 2317 2017-07-20 17:27:19Z suehring
! Bugfixes in reading 3D topography from file
!
! 2274 2017-06-09 13:27:48Z Giersch
! Changed error messages
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! - Adjustments according to new topography representation
! - Bugfix: Move determination of nzb_max behind topography modification in
!   cell-edge case
! - Get rid off global arrays required for topography output
! - Enable topography input via netcdf
! - Generic tunnel set-up added
!
! 2200 2017-04-11 11:37:51Z suehring
! monotonic_adjustment removed
!
! 2169 2017-03-06 18:16:35Z suehring
! Bugfix, move setting for topography grid convention to init_grid, else, if no
! value is set, the simulation may abort in case of restarts
!
! 2128 2017-01-23 15:00:03Z suehring
! Bugfix in setting topography from file in case of ocean simulations
!
! 2088 2016-12-19 16:30:25Z suehring
! Bugfix in generic topography in case of ocean simulations
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2021 2016-10-07 14:08:57Z suehring
! Bugfix: setting Neumann boundary conditions for topography required for
! topography flags in multigrid_noopt solver
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1994 2016-08-15 09:52:21Z suehring
! Bugfix in definition of generic topography
!
! 1982 2016-08-01 11:04:48Z suehring
! Bugfix concering consistency check for topography
!
! 1968 2016-07-18 12:01:49Z suehring
! Changed: PE-wise reading of topography file in order to avoid global definition
! of arrays nzb_local and nzb_tmp. Thereby, topography definition for single
! buildings and street canyons has changed, as well as flag setting for
! multigrid scheme.
!
! Bugfix in checking l_grid anisotropy.
! Simplify initial computation of lwall and vertical_influence, i.e. remove
! nzb_s_inner as it is still zero at this point.
!
! 1942 2016-06-14 12:18:18Z suehring
! Topography filter implemented to fill holes resolved by only one grid point.
! Initialization of flags for ws-scheme moved to advec_ws.
!
! 1931 2016-06-10 12:06:59Z suehring
! Rename multigrid into multigrid_noopt and multigrid_fast into multigrid
!
! 1910 2016-05-26 06:49:46Z raasch
! Bugfix: if topography is read from file, Neumann conditions are used for the
! nzb_local array (instead of cyclic conditions) in case that non-cyclic
! boundary conditions are switched on for the run
!
! 1902 2016-05-09 11:18:56Z suehring
! Set topography flags for multigrid solver only (not for multigrid_fast)
!
! 1886 2016-04-21 11:20:47Z suehring
! Bugfix: setting advection flags near walls
! reformulated index values for nzb_v_inner
! variable discriptions added in declaration block
!
! 1845 2016-04-08 08:29:13Z raasch
! nzb_2d removed
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1779 2016-03-03 08:01:28Z raasch
! coupling_char is trimmed at every place it occurs, because it can have
! different length now
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1743 2016-01-13 10:23:51Z raasch
! Bugfix for calculation of nzb_s_outer and nzb_u_outer at north boundary of
! total domain
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1677 2015-10-02 13:25:23Z boeske
! Bugfix: Ghost points are included in wall_flags_0 and wall_flags_00
!
! 1675 2015-10-02 08:28:59Z gronemeier
! Bugfix: Definition of topography grid levels
!
! 1660 2015-09-21 08:15:16Z gronemeier
! Bugfix: Definition of topography grid levels if vertical grid stretching
!         starts below the maximum topography height.
!
! 1580 2015-04-10 13:43:49Z suehring
! Bugfix: setting flags for 5th order scheme near buildings
!
! 1575 2015-03-27 09:56:27Z raasch
! adjustments for psolver-queries
!
! 1557 2015-03-05 16:43:04Z suehring
! Adjustment for monotoinic limiter
!
! 1418 2014-06-06 13:05:08Z fricke
! Bugfix: Change if-condition for stretched grid in the ocean, with the old
!          condition and a negative value for dz_stretch_level the condition
!          was always true for the whole model domain
!
! 1409 2014-05-23 12:11:32Z suehring
! Bugfix: set wall_flags_0 at inflow and outflow boundary also for i <= nxlu
! j <= nysv
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1221 2013-09-10 08:59:13Z raasch
! wall_flags_00 introduced to hold bits 32-63,
! additional 3D-flag arrays for replacing the 2D-index array nzb_s_inner in
! loops optimized for openACC (pres + flow_statistics)
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1069 2012-11-28 16:18:43Z maronga
! bugfix: added coupling_char to TOPOGRAPHY_DATA to allow topography in the
!         ocean model in case of coupled runs
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! lower index for calculating wall_flags_0 set to nzb_w_inner instead of
! nzb_w_inner+1
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! Bugfix: nzb_max is set to nzt at non-cyclic lateral boundaries
! Bugfix: Set wall_flags_0 for inflow boundary
!
! 927 2012-06-06 19:15:04Z raasch
! Wall flags are not set for multigrid method in case of masking method
!
! 864 2012-03-27 15:10:33Z gryschka
! In case of ocean and Dirichlet bottom bc for u and v dzu_mg and ddzu_pres
! were not correctly defined for k=1.
!
! 861 2012-03-26 14:18:34Z suehring
! Set wall_flags_0. The array is needed for degradation in ws-scheme near walls,
! inflow and outflow boundaries as well as near the bottom and the top of the
! model domain.!
! Initialization of nzb_s_inner and nzb_w_inner.
! gls has to be at least nbgp to do not exceed the array bounds of nzb_local
! while setting wall_flags_0
!
! 843 2012-02-29 15:16:21Z gryschka
! In case of ocean and dirichlet bc for u and v at the bottom
! the first u-level ist defined at same height as the first w-level
!
! 818 2012-02-08 16:11:23Z maronga
! Bugfix: topo_height is only required if topography is used. It is thus now
! allocated in the topography branch
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  1997/08/11 06:17:45  raasch
! Initial revision (Testversion)
!
!
! Description:
! -----------------------------------------------------------------------------!
!> Creating grid depending constants
!> @todo: Rearrange topo flag list
!> @todo: reference 3D buildings on top of orography is not tested and may need
!>        further improvement for steep slopes
!> @todo: Use more advanced setting of building type at filled holes
!------------------------------------------------------------------------------!
 SUBROUTINE init_grid

    USE arrays_3d,                                                             &
        ONLY:  dd2zu, ddzu, ddzu_pres, ddzw, dzu, dzw, zu, zw, ddzuw

    USE control_parameters,                                                    &
        ONLY:  dp_level_ind_b, dz, dz_max, dz_stretch_factor,                  &
               dz_stretch_factor_array, dz_stretch_level, dz_stretch_level_end,&
               dz_stretch_level_end_index, dz_stretch_level_start_index,       &
               dz_stretch_level_start, ibc_uv_b, message_string,               &
               number_stretch_level_end, number_stretch_level_start

    USE grid_variables,                                                        &
        ONLY:  ddx, ddx2, ddy, ddy2, dx, dx2, dy, dy2

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nz,   &
               nzb, nzt

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i                           !< index variable along x
    INTEGER(iwp) ::  j                           !< index variable along y
    INTEGER(iwp) ::  k                           !< index variable along z
    INTEGER(iwp) ::  n                           !< loop variable for stretching
    INTEGER(iwp) ::  number_dz                   !< number of user-specified dz values

    REAL(wp) ::  dz_level_end  !< distance between calculated height level for u/v-grid and user-specified end level for stretching
    REAL(wp) ::  dz_stretched  !< stretched vertical grid spacing

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  min_dz_stretch_level_end !< Array that contains all minimum heights where the stretching can end


!
!-- Calculation of horizontal array bounds including ghost layers
    nxlg = nxl - nbgp
    nxrg = nxr + nbgp
    nysg = nys - nbgp
    nyng = nyn + nbgp

!
!-- Allocate grid arrays
    ALLOCATE( ddzu(1:nzt+1), ddzw(1:nzt+1), dd2zu(1:nzt), dzu(1:nzt+1),        &
              zu(nzb:nzt+1), dzw(1:nzt+1), zw(nzb:nzt+1) )

!
!-- Compute height of u-levels from constant grid length and dz stretch factors
    IF ( dz(1) == -1.0_wp )  THEN
       message_string = 'missing dz'
       CALL message( 'init_grid', 'PA0200', 1, 2, 0, 6, 0 )
    ELSEIF ( dz(1) <= 0.0_wp )  THEN
       WRITE( message_string, * ) 'dz=',dz(1),' <= 0.0'
       CALL message( 'init_grid', 'PA0201', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize dz_stretch_level_start with the value of dz_stretch_level
!-- if it was set by the user
    IF ( dz_stretch_level /= -9999999.9_wp ) THEN
       dz_stretch_level_start(1) = dz_stretch_level
    ENDIF

!
!-- Determine number of dz values and stretching levels specified by the
!-- user to allow right controlling of the stretching mechanism and to
!-- perform error checks
    number_dz = COUNT( dz /= -1.0_wp )
    number_stretch_level_start = COUNT( dz_stretch_level_start /=              &
                                       -9999999.9_wp )
    number_stretch_level_end = COUNT( dz_stretch_level_end /=                  &
                                      9999999.9_wp )

!
!-- The number of specified end levels +1 has to be the same than the number
!-- of specified dz values
    IF ( number_dz /= number_stretch_level_end + 1 ) THEN
       WRITE( message_string, * ) 'The number of values for dz = ',         &
                                   number_dz, 'has to be the same than& ',  &
                                   'the number of values for ',             &
                                   'dz_stretch_level_end + 1 = ',           &
                                   number_stretch_level_end+1
          CALL message( 'init_grid', 'PA0156', 1, 2, 0, 6, 0 )
    ENDIF

!
!--    The number of specified start levels has to be the same or one less than
!--    the number of specified dz values
    IF ( number_dz /= number_stretch_level_start + 1 .AND.                  &
         number_dz /= number_stretch_level_start ) THEN
       WRITE( message_string, * ) 'The number of values for dz = ',         &
                                   number_dz, 'has to be the same or one ', &
                                   'more than& the number of values for ',  &
                                   'dz_stretch_level_start = ',             &
                                   number_stretch_level_start
          CALL message( 'init_grid', 'PA0211', 1, 2, 0, 6, 0 )
    ENDIF

!--    The number of specified start levels has to be the same or one more than
!--    the number of specified end levels
    IF ( number_stretch_level_start /= number_stretch_level_end + 1 .AND.   &
         number_stretch_level_start /= number_stretch_level_end ) THEN
       WRITE( message_string, * ) 'The number of values for ',              &
                                  'dz_stretch_level_start = ',              &
                                   dz_stretch_level_start, 'has to be the ',&
                                   'same or one more than& the number of ', &
                                   'values for dz_stretch_level_end = ',    &
                                   number_stretch_level_end
          CALL message( 'init_grid', 'PA0216', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Initialize dz for the free atmosphere with the value of dz_max
    IF ( dz(number_stretch_level_start+1) == -1.0_wp .AND.                     &
         number_stretch_level_start /= 0 ) THEN
       dz(number_stretch_level_start+1) = dz_max
    ENDIF

!
!-- Initialize the stretching factor if (infinitely) stretching in the free
!-- atmosphere is desired (dz_stretch_level_end was not specified for the
!-- free atmosphere)
    IF ( number_stretch_level_start == number_stretch_level_end + 1 ) THEN
       dz_stretch_factor_array(number_stretch_level_start) =                   &
       dz_stretch_factor
    ENDIF

!
!-- Allocation of arrays for stretching
    ALLOCATE( min_dz_stretch_level_end(number_stretch_level_start) )

!
!-- Define the vertical grid levels
!
!-- The stretching region has to be large enough to allow for a smooth
!-- transition between two different grid spacings
    DO n = 1, number_stretch_level_start
       min_dz_stretch_level_end(n) = dz_stretch_level_start(n) -            &
                                     4 * MAX( dz(n),dz(n+1) )
    ENDDO

    IF ( ANY( min_dz_stretch_level_end (1:number_stretch_level_start) <     &
              dz_stretch_level_end(1:number_stretch_level_start) ) ) THEN
          message_string= 'Eeach dz_stretch_level_end has to be less ' //   &
                          'than its corresponding value for &' //           &
                          'dz_stretch_level_start - 4*MAX(dz(n),dz(n+1)) '//&
                          'to allow for smooth grid stretching'
          CALL message( 'init_grid', 'PA0224', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Stretching must not be applied close to the surface (last two grid
!-- points). For the default case dz_stretch_level_start is negative.
    IF ( ANY( dz_stretch_level_start > - dz(1) * 1.5_wp ) ) THEN
       WRITE( message_string, * ) 'Eeach dz_stretch_level_start has to be ',&
                                  'less than ', dz(1) * 1.5
          CALL message( 'init_grid', 'PA0226', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- The stretching has to start and end on a grid level. Therefore
!-- user-specified values have to ''interpolate'' to the next highest level
    IF ( number_stretch_level_start /= 0 ) THEN
       dz_stretch_level_start(1) = INT( (dz_stretch_level_start(1) +        &
                                         dz(1)/2.0) / dz(1) )               &
                                   * dz(1) - dz(1)/2.0
    ENDIF

    IF ( number_stretch_level_start > 1 ) THEN
       DO n = 2, number_stretch_level_start
          dz_stretch_level_start(n) = INT( dz_stretch_level_start(n) /      &
                                           dz(n) ) * dz(n)
       ENDDO
    ENDIF

    IF ( number_stretch_level_end /= 0 ) THEN
       DO n = 1, number_stretch_level_end
          dz_stretch_level_end(n) = INT( dz_stretch_level_end(n) /          &
                                         dz(n+1) ) * dz(n+1)
       ENDDO
    ENDIF

!
!-- Determine stretching factor if necessary
    IF ( number_stretch_level_end >= 1 ) THEN
       CALL calculate_stretching_factor( number_stretch_level_end )
    ENDIF

!
!-- Grid for ocean with free water surface is at k=nzt (w-grid).
!-- In case of neumann bc at the ground the first first u-level (k=0) lies
!-- below the first w-level (k=0). In case of dirichlet bc the first u- and
!-- w-level are defined at same height, but staggered from the second level.
!-- The second u-level (k=1) corresponds to the top of the Prandtl-layer.
!-- z values are negative starting from z=0 (surface)
    zu(nzt+1) =   dz(1) * 0.5_wp
    zu(nzt)   = - dz(1) * 0.5_wp

!
!-- Determine u and v height levels considering the possibility of grid
!-- stretching in several heights.
    n = 1
    dz_stretch_level_start_index = 0
    dz_stretch_level_end_index = 0
    dz_stretched = dz(1)

    DO  k = nzt-1, 0, -1

       IF ( dz_stretch_level_start(n) >= zu(k+1) ) THEN
          dz_stretched = dz_stretched * dz_stretch_factor_array(n)

          IF ( dz(n) > dz(n+1) ) THEN
             dz_stretched = MAX( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (higher) dz
          ELSE
             dz_stretched = MIN( dz_stretched, dz(n+1) ) !Restrict dz_stretched to the user-specified (lower) dz
          ENDIF

          IF ( dz_stretch_level_start_index(n) == 0 )                             &
          dz_stretch_level_start_index(n) = k+1

       ENDIF

       zu(k) = zu(k+1) - dz_stretched

!
!--    Make sure that the stretching ends exactly at dz_stretch_level_end
       dz_level_end = ABS( zu(k) - dz_stretch_level_end(n) )

       IF ( dz_level_end  < dz(n+1)/3.0 ) THEN
          zu(k) = dz_stretch_level_end(n)
          dz_stretched = dz(n+1)
          dz_stretch_level_end_index(n) = k
          n = n + 1
       ENDIF
    ENDDO


!
!-- Compute the w-levels. They are always staggered half-way between the
!-- corresponding u-levels, except in case of dirichlet bc for u and v
!-- at the ground. In this case the first u- and w-level are defined at
!-- same height. The top w-level (nzt+1) is not used but set for
!-- consistency, since w and all scalar variables are defined up tp nzt+1.
    zw(nzt+1) = dz(1)
    zw(nzt)   = 0.0_wp
    DO  k = 0, nzt
       zw(k) = ( zu(k) + zu(k+1) ) * 0.5_wp
    ENDDO

!
!-- In case of dirichlet bc for u and v the first u- and w-level are defined
!-- at same height.
    IF ( ibc_uv_b == 0 ) THEN
       zu(0) = zw(0)
    ENDIF

!
!-- Compute grid lengths.
    DO  k = 1, nzt+1
       dzu(k)  = zu(k) - zu(k-1)
       ddzu(k) = 1.0_wp / dzu(k)
       dzw(k)  = zw(k) - zw(k-1)
       ddzw(k) = 1.0_wp / dzw(k)
    ENDDO

    DO  k = 1, nzt
       dd2zu(k) = 1.0_wp / ( dzu(k) + dzu(k+1) )
    ENDDO

!
!-- The FFT- SOR-pressure solvers assume grid spacings of a staggered grid
!-- everywhere. For the actual grid, the grid spacing at the lowest level
!-- is only dz/2, but should be dz. Therefore, an additional array
!-- containing with appropriate grid information is created for these
!-- solvers.
    ALLOCATE( ddzu_pres(1:nzt+1),  ddzuw(0:nzt-1,3) )
    ddzu_pres = ddzu
    ddzu_pres(1) = ddzu_pres(2)  ! change for lowest level

    DO  k = 0, nzt-1
       ddzuw(k,1) = ddzu_pres(k+1) * ddzw(k+1)
       ddzuw(k,2) = ddzu_pres(k+2) * ddzw(k+1)
       ddzuw(k,3) = -1.0_wp * &
                    ( ddzu_pres(k+2) * ddzw(k+1) +                          &
                      ddzu_pres(k+1) * ddzw(k+1) )
    ENDDO

!
!-- Compute the reciprocal values of the horizontal grid lengths.
    ddx = 1.0_wp / dx
    ddy = 1.0_wp / dy
    dx2 = dx * dx
    dy2 = dy * dy
    ddx2 = 1.0_wp / dx2
    ddy2 = 1.0_wp / dy2

!
 END SUBROUTINE init_grid


! Description:
! -----------------------------------------------------------------------------!
!> Calculation of the stretching factor through an iterative method. Ideas were
!> taken from the paper "Regional stretched grid generation and its application
!> to the NCAR RegCM (1999)". Normally, no analytic solution exists because the
!> system of equations has two variables (r,l) but four requirements
!> (l=integer, r=[0,88;1,2], Eq(6), Eq(5) starting from index j=1) which
!> results into an overdetermined system.
!------------------------------------------------------------------------------!
 SUBROUTINE calculate_stretching_factor( number_end )

    USE control_parameters,                                                    &
        ONLY:  dz, dz_stretch_factor, dz_stretch_factor_array,                 &
               dz_stretch_level_end, dz_stretch_level_start, message_string

    USE kinds

    IMPLICIT NONE

    INTEGER(iwp) ::  iterations  !< number of iterations until stretch_factor_lower/upper_limit is reached
    INTEGER(iwp) ::  l_rounded   !< after l_rounded grid levels dz(n) is strechted to dz(n+1) with stretch_factor_2
    INTEGER(iwp) ::  n           !< loop variable for stretching

    INTEGER(iwp), INTENT(IN) ::  number_end !< number of user-specified end levels for stretching

    REAL(wp) ::  delta_l               !< absolute difference between l and l_rounded
    REAL(wp) ::  delta_stretch_factor  !< absolute difference between stretch_factor_1 and stretch_factor_2
    REAL(wp) ::  delta_total_new       !< sum of delta_l and delta_stretch_factor for the next iteration (should be as small as possible)
    REAL(wp) ::  delta_total_old       !< sum of delta_l and delta_stretch_factor for the last iteration
    REAL(wp) ::  distance              !< distance between dz_stretch_level_start and dz_stretch_level_end (stretching region)
    REAL(wp) ::  l                     !< value that fulfil Eq. (5) in the paper mentioned above together with stretch_factor_1 exactly
    REAL(wp) ::  numerator             !< numerator of the quotient
    REAL(wp) ::  stretch_factor_1      !< stretching factor that fulfil Eq. (5) togehter with l exactly
    REAL(wp) ::  stretch_factor_2      !< stretching factor that fulfil Eq. (6) togehter with l_rounded exactly

    REAL(wp) ::  dz_stretch_factor_array_2(9) = 1.08_wp  !< Array that contains all stretch_factor_2 that belongs to stretch_factor_1

    REAL(wp), PARAMETER ::  stretch_factor_interval = 1.0E-06  !< interval for sampling possible stretching factors
    REAL(wp), PARAMETER ::  stretch_factor_lower_limit = 0.88  !< lowest possible stretching factor
    REAL(wp), PARAMETER ::  stretch_factor_upper_limit = 1.12  !< highest possible stretching factor


    l = 0
    DO  n = 1, number_end

       iterations = 1
       stretch_factor_1 = 1.0
       stretch_factor_2 = 1.0
       delta_total_old = 1.0

       IF ( dz(n) > dz(n+1) ) THEN
          DO WHILE ( stretch_factor_1 >= stretch_factor_lower_limit )

             stretch_factor_1 = 1.0 - iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) -                   &
                        dz_stretch_level_start(n) )
             numerator = distance*stretch_factor_1/dz(n) +               &
                         stretch_factor_1 - distance/dz(n)

             IF ( numerator > 0.0 ) THEN
                l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0
                l_rounded = NINT( l )
                delta_l = ABS( l_rounded - l ) / l
             ENDIF

             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )

             delta_stretch_factor = ABS( stretch_factor_1 -              &
                                         stretch_factor_2 ) /            &
                                    stretch_factor_2

             delta_total_new = delta_l + delta_stretch_factor

!
!--                stretch_factor_1 is taken to guarantee that the stretching
!--                procedure ends as close as possible to dz_stretch_level_end.
!--                stretch_factor_2 would guarantee that the stretched dz(n) is
!--                equal to dz(n+1) after l_rounded grid levels.
             IF (delta_total_new < delta_total_old) THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF

             iterations = iterations + 1

          ENDDO

       ELSEIF ( dz(n) < dz(n+1) ) THEN
          DO WHILE ( stretch_factor_1 <= stretch_factor_upper_limit )

             stretch_factor_1 = 1.0 + iterations * stretch_factor_interval
             distance = ABS( dz_stretch_level_end(n) -                      &
                        dz_stretch_level_start(n) )
             numerator = distance*stretch_factor_1/dz(n) +                  &
                         stretch_factor_1 - distance/dz(n)

             l = LOG( numerator ) / LOG( stretch_factor_1 ) - 1.0
             l_rounded = NINT( l )
             delta_l = ABS( l_rounded - l ) / l

             stretch_factor_2 = EXP( LOG( dz(n+1)/dz(n) ) / (l_rounded) )

             delta_stretch_factor = ABS( stretch_factor_1 -                 &
                                        stretch_factor_2 ) /                &
                                        stretch_factor_2

             delta_total_new = delta_l + delta_stretch_factor

!
!--                stretch_factor_1 is taken to guarantee that the stretching
!--                procedure ends as close as possible to dz_stretch_level_end.
!--                stretch_factor_2 would guarantee that the stretched dz(n) is
!--                equal to dz(n+1) after l_rounded grid levels.
             IF (delta_total_new < delta_total_old) THEN
                dz_stretch_factor_array(n) = stretch_factor_1
                dz_stretch_factor_array_2(n) = stretch_factor_2
                delta_total_old = delta_total_new
             ENDIF

             iterations = iterations + 1
          ENDDO

       ELSE
          message_string= 'Two adjacent values of dz must be different'
          CALL message( 'init_grid', 'PA0228', 1, 2, 0, 6, 0 )

       ENDIF

!
!--    Check if also the second stretching factor fits into the allowed
!--    interval. If not, print a warning for the user.
       IF ( dz_stretch_factor_array_2(n) < stretch_factor_lower_limit .OR.     &
            dz_stretch_factor_array_2(n) > stretch_factor_upper_limit ) THEN
          WRITE( message_string, * ) 'stretch_factor_2 = ',                    &
                                     dz_stretch_factor_array_2(n), ' which is',&
                                     ' responsible for exactly reaching& dz =',&
                                      dz(n+1), 'after a specific amount of',   &
                                     ' grid levels& exceeds the upper',        &
                                     ' limit =', stretch_factor_upper_limit,   &
                                     ' &or lower limit = ',                    &
                                     stretch_factor_lower_limit
          CALL message( 'init_grid', 'PA0499', 0, 1, 0, 6, 0 )

       ENDIF
    ENDDO

 END SUBROUTINE calculate_stretching_factor


