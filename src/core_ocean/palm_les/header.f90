!> @file header.f90
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
! 2018-10-25 cbegeman
! Add dirichlet bottom boundary conditions for salinity
!
! Former revisions:
! -----------------
! $Id: header.f90 3083 2018-06-19 14:03:12Z gronemeier $
! Print RANS-mode constants
!
! 3065 2018-06-12 07:03:02Z Giersch
! Header output concerning stretching revised
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
!
! 2967 2018-04-13 11:22:08Z raasch
! bugfix: missing parallel cpp-directives added
!
! 2883 2018-03-14 08:29:10Z Giersch
! Format of the output of dt_dopr_listing (325) has been changed
!
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
!
! 2776 2018-01-31 10:44:42Z Giersch
! Variable synthetic_turbulence_generator has been abbreviated
!
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
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
! Print information about turbulence closure (TG)
! Print information about inifor initialization (MS)
!
! 2575 2017-10-24 09:57:58Z maronga
! Added output for complex terrain simulations
!
! 2544 2017-10-13 18:09:32Z maronga
! Moved initial day of year and time to inipar.
!
! 2339 2017-08-07 13:55:26Z gronemeier
! corrected timestamp in header
!
! 2338 2017-08-07 12:15:38Z gronemeier
! Modularize 1D model
!
! 2320 2017-07-21 12:47:43Z suehring
! Modularize large-scale forcing and nudging
!
! 2300 2017-06-29 13:31:14Z raasch
! host-specific code removed
!
! 2299 2017-06-29 10:14:38Z maronga
! Modified output for spinups
!
! 2298 2017-06-29 09:28:18Z raasch
! MPI2 related parts removed
!
! 2270 2017-06-09 12:18:47Z maronga
! Renamed Prandtl layer to constant flux layer
!
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2258 2017-06-08 07:55:13Z suehring
! Bugfix, add pre-preprocessor directives to enable non-parrallel mode
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
! Generic tunnel setup added
!
! 2200 2017-04-11 11:37:51Z suehring
! monotonic_adjustment removed
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC relatec code removed
!
! 2073 2016-11-30 14:34:05Z raasch
! small bugfix concerning output of scalar profiles
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1992 2016-08-12 15:14:59Z suehring
! Adapted for top_scalarflux
!
! 1960 2016-07-12 16:34:24Z suehring
! Treat humidity and passive scalar separately.
! Modify misleading information concerning humidity.
! Bugfix, change unit for humidity flux.
!
! 1957 2016-07-07 10:43:48Z suehring
! flight module added
!
! 1931 2016-06-10 12:06:59Z suehring
! Rename multigrid into multigrid_noopt
!
! 1902 2016-05-09 11:18:56Z suehring
! Write information about masking_method only for multigrid solver
!
! 1849 2016-04-08 11:33:18Z hoffmann
! Adapted for modularization of microphysics
!
! 1833 2016-04-07 14:23:03Z raasch
! spectrum renamed spectra_mod, output of spectra related quantities moved to
! spectra_mod
!
! 1831 2016-04-07 13:15:51Z hoffmann
! turbulence renamed collision_turbulence,
! drizzle renamed cloud_water_sedimentation
!
! 1826 2016-04-07 12:01:39Z maronga
! Moved radiation model header output to the respective module.
! Moved canopy model header output to the respective module.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Tails removed. icloud_scheme replaced by microphysics_*
!
! 1817 2016-04-06 15:44:20Z maronga
! Moved land_surface_model header output to the respective module.
!
! 1808 2016-04-05 19:44:00Z raasch
! routine local_flush replaced by FORTRAN statement
!
! 1797 2016-03-21 16:50:28Z raasch
! output of nesting datatransfer mode
!
! 1791 2016-03-11 10:41:25Z raasch
! output of nesting informations of all domains
!
! 1788 2016-03-10 11:01:04Z maronga
! Parameter dewfall removed
!
! 1786 2016-03-08 05:49:27Z raasch
! cpp-direktives for spectra removed
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf module and variable names changed, output of netcdf_deflate
!
! 1764 2016-02-28 12:45:19Z raasch
! output of nesting informations
!
! 1697 2015-10-28 17:14:10Z raasch
! small E- and F-FORMAT changes to avoid informative compiler messages about
! insufficient field width
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed prandtl_layer to constant_flux_layer, renames rif_min/rif_max to
! zeta_min/zeta_max.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1675 2015-10-02 08:28:59Z gronemeier
! Bugfix: Definition of topography grid levels
!
! 1660 2015-09-21 08:15:16Z gronemeier
! Bugfix: Definition of building/street canyon height if vertical grid stretching
!         starts below the maximum topography height.
!
! 1590 2015-05-08 13:56:27Z maronga
! Bugfix: Added TRIM statements for character strings for LSM and radiation code
!
! 1585 2015-04-30 07:05:52Z maronga
! Further output for radiation model(s).
!
! 1575 2015-03-27 09:56:27Z raasch
! adjustments for psolver-queries, output of seed_follows_topography
!
! 1560 2015-03-06 10:48:54Z keck
! output for recycling y shift
!
! 1557 2015-03-05 16:43:04Z suehring
! output for monotonic limiter
!
! 1551 2015-03-03 14:18:16Z maronga
! Added informal output for land surface model and radiation model. Removed typo.
!
! 1496 2014-12-02 17:25:50Z maronga
! Renamed: "radiation -> "cloud_top_radiation"
!
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   module plant_canopy_model_mod and output for new canopy model parameters
!   (alpha_lad, beta_lad, lai_beta,...) added,
!   drag_coefficient, leaf_surface_concentration and scalar_exchange_coefficient
!   renamed to canopy_drag_coeff, leaf_surface_conc and leaf_scalar_exch_coeff,
!   learde renamed leaf_area_density.
! Bugfix: DO-WHILE-loop for lad header information additionally restricted
! by maximum number of gradient levels (currently 10)
!
! 1482 2014-10-18 12:34:45Z raasch
! information about calculated or predefined virtual processor topology adjusted
!
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores
!
! 1429 2014-07-15 12:53:45Z knoop
! header exended to provide ensemble_member_nr if specified
!
! 1376 2014-04-26 11:21:22Z boeske
! Correction of typos
!
! 1365 2014-04-22 15:03:56Z boeske
! New section 'Large scale forcing and nudging':
! output of large scale forcing and nudging information,
! new section for initial profiles created
!
! 1359 2014-04-11 17:15:14Z hoffmann
! dt_sort_particles removed
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1327 2014-03-21 11:00:16Z raasch
! parts concerning iso2d and avs output removed,
! -netcdf output queries
!
! 1324 2014-03-21 09:13:16Z suehring
! Bugfix: module spectrum added
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL functions provided with KIND-attribute,
! some REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1308 2014-03-13 14:58:42Z fricke
! output of the fixed number of output time levels
! output_format adjusted for masked data if netcdf_data_format > 5
!
! 1299 2014-03-06 13:15:21Z heinze
! output for using large_scale subsidence in combination
! with large_scale_forcing
! reformatting, more detailed explanations
!
! 1241 2013-10-30 11:36:58Z heinze
! output for nudging + large scale forcing from external file
!
! 1216 2013-08-26 09:31:42Z raasch
! output for transpose_compute_overlap
!
! 1212 2013-08-15 08:46:27Z raasch
! output for poisfft_hybrid removed
!
! 1179 2013-06-14 05:57:58Z raasch
! output of reference_state, use_reference renamed use_single_reference_value
!
! 1159 2013-05-21 11:58:22Z fricke
! +use_cmax
!
! 1115 2013-03-26 18:16:16Z hoffmann
! descriptions for Seifert-Beheng-cloud-physics-scheme added
!
! 1111 2013-03-08 23:54:10Z raasch
! output of accelerator board information
! ibc_p_b = 2 removed
!
! 1108 2013-03-05 07:03:32Z raasch
! bugfix for r1106
!
! 1106 2013-03-04 05:31:38Z raasch
! some format changes for coupled runs
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031 2012-10-19 14:35:30Z raasch
! output of netCDF data format modified
!
! 1015 2012-09-27 09:23:24Z raasch
! output of Adjustment of mixing length to the Prandtl mixing length at first
! grid point above ground removed
!
! 1003 2012-09-14 14:35:53Z raasch
! output of information about equal/unequal subdomain size removed
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog- and upstream-spline-scheme removed
!
! 978 2012-08-09 08:28:32Z fricke
! -km_damp_max, outflow_damping_width
! +pt_damping_factor, pt_damping_width
! +z0h
!
! 964 2012-07-26 09:14:24Z raasch
! output of profil-related quantities removed
!
! 940 2012-07-09 14:31:00Z raasch
! Output in case of simulations for pure neutral stratification (no pt-equation
! solved)
!
! 927 2012-06-06 19:15:04Z raasch
! output of masking_method for mg-solver
!
! 868 2012-03-28 12:21:07Z raasch
! translation velocity in Galilean transformation changed to 0.6 * ug
!
! 833 2012-02-22 08:55:55Z maronga
! Adjusted format for leaf area density
!
! 828 2012-02-21 12:00:36Z raasch
! output of dissipation_classes + radius_classes
!
! 825 2012-02-19 03:03:44Z raasch
! Output of cloud physics parameters/quantities complemented and restructured
!
! Revision 1.1  1997/08/11 06:17:20  raasch
! Initial revision
!
!
! Description:
! ------------
!> Writing a header with all important information about the current run.
!> This subroutine is called three times, two times at the beginning
!> (writing information on files RUN_CONTROL and HEADER) and one time at the
!> end of the run, then writing additional information about CPU-usage on file
!> header.
!-----------------------------------------------------------------------------!
 SUBROUTINE header


    USE arrays_3d,                                                             &
        ONLY:  pt_init, sa_init, zu, zw

    USE control_parameters

    USE cpulog,                                                                &
        ONLY:  log_point_s

    USE date_and_time_mod,                                                     &
        ONLY:  day_of_year_init, time_utc_init

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nnx, nny, nnz, nx, ny, nzt

    USE kinds

    USE netcdf_interface,                                                      &
        ONLY:  netcdf_data_format, netcdf_data_format_string, netcdf_deflate

    USE pegrid

    IMPLICIT NONE

    CHARACTER (LEN=1)  ::  prec                !<

    CHARACTER (LEN=5)  ::  section_chr         !<

    CHARACTER (LEN=10) ::  coor_chr            !<
    CHARACTER (LEN=10) ::  host_chr            !<

    CHARACTER (LEN=16) ::  begin_chr           !<

    CHARACTER (LEN=26) ::  ver_rev             !<

    CHARACTER (LEN=32) ::  cpl_name            !<

    CHARACTER (LEN=40) ::  output_format       !<

    CHARACTER (LEN=70) ::  char1               !<
    CHARACTER (LEN=70) ::  char2               !<
    CHARACTER (LEN=70) ::  dopr_chr            !<
    CHARACTER (LEN=70) ::  do3d_chr            !<
    CHARACTER (LEN=70) ::  domask_chr          !<
    CHARACTER (LEN=70) ::  run_classification  !<

    CHARACTER (LEN=85) ::  r_upper             !<
    CHARACTER (LEN=85) ::  r_lower             !<

    CHARACTER (LEN=86) ::  coordinates         !<
    CHARACTER (LEN=86) ::  gradients           !<
    CHARACTER (LEN=86) ::  slices              !<
    CHARACTER (LEN=86) ::  temperatures        !<
    CHARACTER (LEN=86) ::  salinity            !<

    CHARACTER (LEN=1), DIMENSION(1:3) ::  dir = (/ 'x', 'y', 'z' /)  !<

    INTEGER(iwp) ::  av             !<
    INTEGER(iwp) ::  bh             !<
    INTEGER(iwp) ::  blx            !<
    INTEGER(iwp) ::  bly            !<
    INTEGER(iwp) ::  bxl            !<
    INTEGER(iwp) ::  bxr            !<
    INTEGER(iwp) ::  byn            !<
    INTEGER(iwp) ::  bys            !<
    INTEGER(iwp) ::  ch             !<
    INTEGER(iwp) ::  count          !<
    INTEGER(iwp) ::  cpl_parent_id  !<
    INTEGER(iwp) ::  cwx            !<
    INTEGER(iwp) ::  cwy            !<
    INTEGER(iwp) ::  cxl            !<
    INTEGER(iwp) ::  cxr            !<
    INTEGER(iwp) ::  cyn            !<
    INTEGER(iwp) ::  cys            !<
    INTEGER(iwp) ::  dim            !<
    INTEGER(iwp) ::  i              !<
    INTEGER(iwp) ::  io             !<
    INTEGER(iwp) ::  j              !<
    INTEGER(iwp) ::  k              !<
    INTEGER(iwp) ::  l              !<
    INTEGER(iwp) ::  ll             !<
    INTEGER(iwp) ::  my_cpl_id      !<
    INTEGER(iwp) ::  n              !<
    INTEGER(iwp) ::  ncpl           !<
    INTEGER(iwp) ::  npe_total      !<


    REAL(wp) ::  cpuseconds_per_simulated_second  !<
    REAL(wp) ::  lower_left_coord_x               !< x-coordinate of nest domain
    REAL(wp) ::  lower_left_coord_y               !< y-coordinate of nest domain

!
!-- Open the output file. At the end of the simulation, output is directed
!-- to unit 19.
    IF ( ( runnr == 0 .OR. force_print_header )  .AND. &
         .NOT. simulated_time_at_begin /= simulated_time )  THEN
       io = 15   !  header output on file RUN_CONTROL
    ELSE
       io = 19   !  header output on file HEADER
    ENDIF
    CALL check_open( io )

!
!-- At the end of the run, output file (HEADER) will be rewritten with
!-- new information
    IF ( io == 19 .AND. simulated_time_at_begin /= simulated_time ) REWIND( 19 )

!
!-- Determine kind of model run
    IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
       run_classification = 'restart run'
    ELSEIF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0 )  THEN
       run_classification = 'run without 1D - prerun'
    ELSE
       message_string = ' unknown action(s): ' // TRIM( initializing_actions )
       CALL message( 'header', 'PA0191', 0, 0, 0, 6, 0 )
    ENDIF
       run_classification = 'ocean - ' // run_classification

!
!-- Run-identification, date, time, host
    host_chr = host(1:10)
    ver_rev = TRIM( version ) // '  ' // TRIM( revision )
    WRITE ( io, 100 )  ver_rev, TRIM( run_classification )
#if defined( __parallel )
    IF ( npex == -1  .AND.  npey == -1 )  THEN
       char1 = 'calculated'
    ELSE
       char1 = 'predefined'
    ENDIF
    IF ( threads_per_task == 1 )  THEN
       WRITE ( io, 103 )  numprocs, pdims(1), pdims(2), TRIM( char1 )
    ELSE
       WRITE ( io, 104 )  numprocs*threads_per_task, numprocs, &
                          threads_per_task, pdims(1), pdims(2), TRIM( char1 )
    ENDIF

    IF ( pdims(2) == 1 )  THEN
       WRITE ( io, 107 )  'x'
    ELSEIF ( pdims(1) == 1 )  THEN
       WRITE ( io, 107 )  'y'
    ENDIF
    IF ( numprocs /= maximum_parallel_io_streams )  THEN
       WRITE ( io, 108 )  maximum_parallel_io_streams
    ENDIF
#endif

!
!-- Numerical schemes
    WRITE ( io, 110 )
    IF ( transpose_compute_overlap )  WRITE( io, 115 )
    IF ( call_psolver_at_all_substeps ) &
    THEN
       WRITE ( io, 142 )
    ENDIF

    WRITE ( io, 503 )
    WRITE ( io, 504 )
    WRITE ( io, 122 )
    IF ( rayleigh_damping_factor /= 0.0_wp )  THEN
       WRITE ( io, 123 )  'below', rayleigh_damping_height, &
            rayleigh_damping_factor
    ENDIF
    IF ( dp_external )  THEN
       IF ( dp_smooth )  THEN
          WRITE ( io, 152 )  dpdxy, dp_level_b, ', vertically smoothed.'
       ELSE
          WRITE ( io, 152 )  dpdxy, dp_level_b, '.'
       ENDIF
    ENDIF
    WRITE ( io, 99 )

!
!-- Runtime and timestep information
    WRITE ( io, 200 )
    IF ( .NOT. dt_fixed )  THEN
       WRITE ( io, 201 )  dt_max, cfl_factor
    ELSE
       WRITE ( io, 202 )  dt
    ENDIF
    WRITE ( io, 203 )  simulated_time_at_begin, end_time

    IF ( time_restart /= 9999999.9_wp  .AND. &
         simulated_time_at_begin == simulated_time )  THEN
       IF ( dt_restart == 9999999.9_wp )  THEN
          WRITE ( io, 204 )  ' Restart at:       ',time_restart
       ELSE
          WRITE ( io, 205 )  ' Restart at:       ',time_restart, dt_restart
       ENDIF
    ENDIF

    IF ( simulated_time_at_begin /= simulated_time )  THEN
       i = MAX ( log_point_s(10)%counts, 1 )
       IF ( ( simulated_time - simulated_time_at_begin ) == 0.0_wp )  THEN
          cpuseconds_per_simulated_second = 0.0_wp
       ELSE
          cpuseconds_per_simulated_second = log_point_s(10)%sum / &
                                            ( simulated_time -    &
                                              simulated_time_at_begin )
       ENDIF
       WRITE ( io, 206 )  simulated_time, log_point_s(10)%sum,      &
                          log_point_s(10)%sum / REAL( i, KIND=wp ), &
                          cpuseconds_per_simulated_second
       IF ( time_restart /= 9999999.9_wp  .AND.  time_restart < end_time )  THEN
          IF ( dt_restart == 9999999.9_wp )  THEN
             WRITE ( io, 204 )  ' Next restart at:     ',time_restart
          ELSE
             WRITE ( io, 205 )  ' Next restart at:     ',time_restart, dt_restart
          ENDIF
       ENDIF
    ENDIF


!
!-- Start time for coupled runs, if independent precursor runs for atmosphere
!-- and ocean are used or have been used. In this case, coupling_start_time
!-- defines the time when the coupling is switched on.
    WRITE ( io, 250 )  dx, dy
    DO i = 1, number_stretch_level_start+1
       WRITE ( io, 253 )  i, dz(i)
    ENDDO

    WRITE ( io, 251 ) (nx+1)*dx, (ny+1)*dy, zu(0)

    IF ( ANY( dz_stretch_level_start_index > 0 ) )  THEN
       WRITE( io, '(A)', advance='no') ' Vertical stretching starts at height:'
       DO i = 1, number_stretch_level_start
          WRITE ( io, '(F10.1,A3)', advance='no' )  dz_stretch_level_start(i), ' m,'
       ENDDO
       WRITE( io, '(/,A)', advance='no') ' Vertical stretching starts at index: '
       DO i = 1, number_stretch_level_start
          WRITE ( io, '(I12,A1)', advance='no' )  dz_stretch_level_start_index(i), ','
       ENDDO
       WRITE( io, '(/,A)', advance='no') ' Vertical stretching ends at height:  '
       DO i = 1, number_stretch_level_start
          WRITE ( io, '(F10.1,A3)', advance='no' )  dz_stretch_level_end(i), ' m,'
       ENDDO
       WRITE( io, '(/,A)', advance='no') ' Vertical stretching ends at index:   '
       DO i = 1, number_stretch_level_start
          WRITE ( io, '(I12,A1)', advance='no' )  dz_stretch_level_end_index(i), ','
       ENDDO
       WRITE( io, '(/,A)', advance='no') ' Factor used for stretching:          '
       DO i = 1, number_stretch_level_start
          WRITE ( io, '(F12.3,A1)', advance='no' )  dz_stretch_factor_array(i), ','
       ENDDO
    ENDIF
    WRITE ( io, 254 )  nx, ny, nzt+1, MIN( nnx, nx+1 ), MIN( nny, ny+1 ),      &
                       MIN( nnz+2, nzt+2 )

!
!-- Boundary conditions
    IF ( ibc_p_b == 0 )  THEN
       r_lower = 'p(0)     = 0      |'
    ELSEIF ( ibc_p_b == 1 )  THEN
       r_lower = 'p(0)     = p(1)   |'
    ENDIF
    IF ( ibc_p_t == 0 )  THEN
       r_upper  = 'p(nzt+1) = 0      |'
    ELSE
       r_upper  = 'p(nzt+1) = p(nzt) |'
    ENDIF

    IF ( ibc_uv_b == 0 )  THEN
       r_lower = TRIM( r_lower ) // ' uv(0)     = -uv(1)                |'
    ELSE
       r_lower = TRIM( r_lower ) // ' uv(0)     = uv(1)                 |'
    ENDIF
    IF ( TRIM( bc_uv_t ) == 'dirichlet_0' )  THEN
       r_upper  = TRIM( r_upper  ) // ' uv(nzt+1) = 0                     |'
    ELSEIF ( ibc_uv_t == 0 )  THEN
       r_upper  = TRIM( r_upper  ) // ' uv(nzt+1) = ug(nzt+1), vg(nzt+1)  |'
    ELSE
       r_upper  = TRIM( r_upper  ) // ' uv(nzt+1) = uv(nzt)               |'
    ENDIF

    IF ( ibc_pt_b == 0 )  THEN
          r_lower = TRIM( r_lower ) // ' pt(0)     = pt_surface'
    ELSEIF ( ibc_pt_b == 1 )  THEN
       r_lower = TRIM( r_lower ) // ' pt(0)     = pt(1)'
    ELSEIF ( ibc_pt_b == 2 )  THEN
       r_lower = TRIM( r_lower ) // ' pt(0)     = from coupled model'
    ENDIF
    IF ( ibc_pt_t == 0 )  THEN
       r_upper  = TRIM( r_upper  ) // ' pt(nzt+1) = pt_top'
    ELSEIF( ibc_pt_t == 1 )  THEN
       r_upper  = TRIM( r_upper  ) // ' pt(nzt+1) = pt(nzt)'
    ELSEIF( ibc_pt_t == 2 )  THEN
       r_upper  = TRIM( r_upper  ) // ' pt(nzt+1) = pt(nzt) + dpt/dz_ini'

    ENDIF

    WRITE ( io, 300 )  r_lower, r_upper

    IF ( ibc_e_b == 1 )  THEN
       r_lower = 'e(0)     = e(1)'
    ELSE
       r_lower = 'e(0)     = e(1) = (u*/0.1)**2'
    ENDIF
    r_upper = 'e(nzt+1) = e(nzt) = e(nzt-1)'

    WRITE ( io, 301 )  'e', r_lower, r_upper


    IF ( ibc_sa_b == 0 ) THEN
       r_lower = 'sa(0)    = sa_surface'
    ELSE
       r_lower = 'sa(0)    = sa(1)'
    ENDIF
    IF ( ibc_sa_t == 0 )  THEN
       r_upper =  'sa(nzt+1) = sa_top'
    ELSE
       r_upper =  'sa(nzt+1) = sa(nzt)'
    ENDIF
    WRITE ( io, 301 ) 'sa', r_lower, r_upper

    WRITE ( io, 304 )
    WRITE ( io, 320 )  top_momentumflux_u, top_momentumflux_v
    IF ( constant_top_heatflux )  THEN
       WRITE ( io, 306 )  top_heatflux
    ENDIF
    IF ( constant_top_salinityflux )                          &
       WRITE ( io, 309 )  top_salinityflux

    WRITE ( io, 317 )

!
!-- Initial Profiles
    WRITE ( io, 321 )
!
!-- Initial wind profiles
    IF ( u_profile(1) /= 9999999.9_wp )  WRITE ( io, 427 )

!
!-- Initial temperature profile
!-- Building output strings, starting with surface temperature
    WRITE ( temperatures, '(F6.2)' )  pt_surface
    gradients = '------'
    slices = '     0'
    coordinates = '   0.0'
    i = 1
    DO  WHILE ( pt_vertical_gradient_level_ind(i) /= -9999 )

       WRITE (coor_chr,'(F7.2)')  pt_init(pt_vertical_gradient_level_ind(i))
       temperatures = TRIM( temperatures ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F7.2)')  pt_vertical_gradient(i)
       gradients = TRIM( gradients ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(I7)')  pt_vertical_gradient_level_ind(i)
       slices = TRIM( slices ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F7.1)')  pt_vertical_gradient_level(i)
       coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

       IF ( i == 10 )  THEN
          EXIT
       ELSE
          i = i + 1
       ENDIF

    ENDDO

    WRITE ( io, 420 )  TRIM( coordinates ), TRIM( temperatures ), &
                       TRIM( gradients ), TRIM( slices )
!
!-- Initial salinity profile
!-- Building output strings, starting with surface salinity
    WRITE ( salinity, '(F6.2)' )  sa_surface
    gradients = '------'
    slices = '     0'
    coordinates = '   0.0'
    i = 1
    DO  WHILE ( sa_vertical_gradient_level_ind(i) /= -9999 )

       WRITE (coor_chr,'(F7.2)')  sa_init(sa_vertical_gradient_level_ind(i))
       salinity = TRIM( salinity ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F7.2)')  sa_vertical_gradient(i)
       gradients = TRIM( gradients ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(I7)')  sa_vertical_gradient_level_ind(i)
       slices = TRIM( slices ) // ' ' // TRIM( coor_chr )

       WRITE (coor_chr,'(F7.1)')  sa_vertical_gradient_level(i)
       coordinates = TRIM( coordinates ) // ' '  // TRIM( coor_chr )

       IF ( i == 10 )  THEN
          EXIT
       ELSE
          i = i + 1
       ENDIF

    ENDDO

    WRITE ( io, 425 )  TRIM( coordinates ), TRIM( salinity ), &
                       TRIM( gradients ), TRIM( slices )

!-- Reference potential density
    WRITE ( io, 426) rho_reference

!
!-- Listing of 1D-profiles
    WRITE ( io, 325 )  dt_dopr_listing
    IF ( averaging_interval_pr /= 0.0_wp )  THEN
       WRITE ( io, 326 )  averaging_interval_pr, dt_averaging_input_pr
    ENDIF

!
!-- DATA output
    WRITE ( io, 330 )
    IF ( averaging_interval_pr /= 0.0_wp )  THEN
       WRITE ( io, 326 )  averaging_interval_pr, dt_averaging_input_pr
    ENDIF

!
!-- 1D-profiles
    dopr_chr = 'Profile:'
    IF ( dopr_n /= 0 )  THEN
       WRITE ( io, 331 )

       output_format = ''
       output_format = netcdf_data_format_string
       IF ( netcdf_deflate == 0 )  THEN
          WRITE ( io, 344 )  output_format
       ELSE
          WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
       ENDIF

       DO  i = 1, dopr_n
          dopr_chr = TRIM( dopr_chr ) // ' ' // TRIM( data_output_pr(i) ) // ','
          IF ( LEN_TRIM( dopr_chr ) >= 60 )  THEN
             WRITE ( io, 332 )  dopr_chr
             dopr_chr = '       :'
          ENDIF
       ENDDO

       IF ( dopr_chr /= '' )  THEN
          WRITE ( io, 332 )  dopr_chr
       ENDIF
       WRITE ( io, 333 )  dt_dopr, averaging_interval_pr, dt_averaging_input_pr
       IF ( skip_time_dopr /= 0.0_wp )  WRITE ( io, 339 )  skip_time_dopr
    ENDIF

!
!-- 3d-arrays
    DO  av = 0, 1

       i = 1
       do3d_chr = ''
       DO  WHILE ( do3d(av,i) /= ' ' )

          do3d_chr = TRIM( do3d_chr ) // ' ' // TRIM( do3d(av,i) ) // ','
          i = i + 1

       ENDDO

       IF ( do3d_chr /= '' )  THEN
          IF ( av == 0 )  THEN
             WRITE ( io, 336 )  ''
          ELSE
             WRITE ( io, 336 )  '(time-averaged)'
          ENDIF

          output_format = netcdf_data_format_string
          IF ( netcdf_deflate == 0 )  THEN
             WRITE ( io, 344 )  output_format
          ELSE
             WRITE ( io, 354 )  TRIM( output_format ), netcdf_deflate
          ENDIF

          IF ( do3d_at_begin )  THEN
             begin_chr = 'and at the start'
          ELSE
             begin_chr = ''
          ENDIF
          IF ( av == 0 )  THEN
             WRITE ( io, 337 )  do3d_chr, dt_do3d, TRIM( begin_chr ), &
                                zu(nz_do3d), nz_do3d
          ELSE
             WRITE ( io, 343 )  do3d_chr, dt_data_output_av,           &
                                TRIM( begin_chr ), averaging_interval, &
                                dt_averaging_input, zu(nz_do3d), nz_do3d
          ENDIF

          IF ( netcdf_data_format > 4 )  THEN
             WRITE ( io, 352 )  ntdim_3d(av)
          ELSE
             WRITE ( io, 353 )
          ENDIF

          IF ( av == 0 )  THEN
             IF ( skip_time_do3d /= 0.0_wp )  THEN
                WRITE ( io, 339 )  skip_time_do3d
             ENDIF
          ELSE
             IF ( skip_time_data_output_av /= 0.0_wp )  THEN
                WRITE ( io, 339 )  skip_time_data_output_av
             ENDIF
          ENDIF

       ENDIF

    ENDDO

!
    WRITE ( io, 99 )

!
!-- Physical quantities
    WRITE ( io, 400 )

!
!-- Geostrophic parameters
    WRITE ( io, 410 )  latitude, longitude, omega, f, fs

 !
!-- Geostrophic parameters
    WRITE ( io, 456 )  day_of_year_init, time_utc_init

!
!-- Other quantities
    WRITE ( io, 411 )  g


!-- LES / turbulence parameters
    WRITE ( io, 450 )

!--
! ... LES-constants used must still be added here
!--
    IF ( e_init > 0.0_wp )  WRITE ( io, 455 )  e_init
    IF ( e_min > 0.0_wp )  WRITE ( io, 454 )  e_min
    IF ( wall_adjustment )  WRITE ( io, 453 )  wall_adjustment_factor
!
!-- Special actions during the run
    WRITE ( io, 470 )
    IF ( create_disturbances )  THEN
       WRITE ( io, 471 )  dt_disturb, disturbance_amplitude,                   &
                          zu(disturbance_level_ind_b), disturbance_level_ind_b,&
                          zu(disturbance_level_ind_t), disturbance_level_ind_t
       WRITE ( io, 473 )  disturbance_energy_limit
       WRITE ( io, 474 )  TRIM( random_generator )
    ENDIF
    IF ( pt_surface_initial_change /= 0.0_wp )  THEN
       WRITE ( io, 475 )  pt_surface_initial_change
    ENDIF
    WRITE ( io, 99 )

!
!-- Write buffer contents to disc immediately
    FLUSH( io )

!
!-- Here the FORMATs start

 99 FORMAT (1X,78('-'))
100 FORMAT (/1X,'******************************',4X,44('-')/        &
            1X,'* ',A,' *',4X,A/                               &
            1X,'******************************',4X,44('-'))
#if defined( __parallel )
103 FORMAT (' Number of PEs:',10X,I6,4X,'Processor grid (x,y): (',I4,',',I4, &
              ')',1X,A)
104 FORMAT (' Number of PEs:',10X,I6,4X,'Tasks:',I4,'   threads per task:',I4/ &
              35X,'Processor grid (x,y): (',I4,',',I4,')',1X,A)
107 FORMAT (35X,'A 1d-decomposition along ',A,' is used')
108 FORMAT (35X,'Max. # of parallel I/O streams is ',I5)
#endif
110 FORMAT (/' Numerical Schemes:'/ &
             ' -----------------'/)
115 FORMAT ('     FFT and transpositions are overlapping')
122 FORMAT (' --> Time differencing scheme: runge-kutta-3')
123 FORMAT (' --> Rayleigh-Damping active, starts ',A,' z = ',F8.2,' m'/ &
            '     maximum damping coefficient:',F6.3, ' 1/s')
142 FORMAT (' --> Perturbation pressure is calculated at every Runge-Kutta ', &
                  'step')
152 FORMAT (' --> External pressure gradient directly prescribed by the user:',&
           /'     ',2(1X,E12.5),'Pa/m in x/y direction', &
           /'     starting from dp_level_b =', F8.3, 'm', A /)
200 FORMAT (//' Run time and time step information:'/ &
             ' ----------------------------------'/)
201 FORMAT ( ' Timestep:             variable     maximum value: ',F6.3,' s', &
             '    CFL-factor:',F5.2)
202 FORMAT ( ' Timestep:          dt = ',F6.3,' s'/)
203 FORMAT ( ' Start time:          ',F9.3,' s'/ &
             ' End time:            ',F9.3,' s')
204 FORMAT ( A,F9.3,' s')
205 FORMAT ( A,F9.3,' s',5X,'restart every',17X,F9.3,' s')
206 FORMAT (/' Time reached:        ',F9.3,' s'/ &
             ' CPU-time used:       ',F9.3,' s     per timestep:               ', &
               '  ',F9.3,' s'/                                                    &
             '                                      per second of simulated tim', &
               'e: ',F9.3,' s')
250 FORMAT (//' Computational grid and domain size:'/ &
              ' ----------------------------------'// &
              ' Grid length:      dx =    ',F8.3,' m    dy =    ',F8.3, ' m')
251 FORMAT (  /' Domain size:       x = ',F10.3,' m     y = ',F10.3, &
              ' m  z(u) = ',F10.3,' m'/)
253 FORMAT ( '                dz(',I1,') =    ', F8.3, ' m')
254 FORMAT (//' Number of gridpoints (x,y,z):  (0:',I4,', 0:',I4,', 0:',I4,')'/ &
            ' Subdomain size (x,y,z):        (  ',I4,',   ',I4,',   ',I4,')'/)
300 FORMAT (//' Boundary conditions:'/ &
             ' -------------------'// &
             '                     p                    uv             ', &
             '                     pt'// &
             ' B. bound.: ',A/ &
             ' T. bound.: ',A)
301 FORMAT (/'                     ',A// &
             ' B. bound.: ',A/ &
             ' T. bound.: ',A)
304 FORMAT (/' Top surface fluxes are used in diffusion terms at k=nzt')
306 FORMAT ('       Predefined constant heatflux:   ',F9.6,' K m/s')
309 FORMAT ('       Predefined constant salinityflux:   ',F9.6,' psu m/s')
317 FORMAT (//' Lateral boundaries:'/ &
            '       left/right:  cyclic'/    &
            '       north/south: cyclic')
320 FORMAT ('       Predefined constant momentumflux:  u: ',F9.6,' m**2/s**2'/ &
            '                                          v: ',F9.6,' m**2/s**2')
321 FORMAT (//' Initial profiles:'/ &
              ' ----------------')
325 FORMAT (//' List output:'/ &
             ' -----------'//  &
            '    1D-Profiles:'/    &
            '       Output every             ',F10.2,' s')
326 FORMAT ('       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s')
330 FORMAT (//' Data output:'/ &
             ' -----------'/)
331 FORMAT (/'    1D-Profiles:')
332 FORMAT (/'       ',A)
333 FORMAT ('       Output every             ',F8.2,' s',/ &
            '       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s')
334 FORMAT (/'    2D-Arrays',A,':')
335 FORMAT (/'       ',A2,'-cross-section  Arrays: ',A/ &
            '       Output every             ',F8.2,' s  ',A/ &
            '       Cross sections at ',A1,' = ',A/ &
            '       scalar-coordinates:   ',A,' m'/)
336 FORMAT (/'    3D-Arrays',A,':')
337 FORMAT (/'       Arrays: ',A/ &
            '       Output every             ',F8.2,' s  ',A/ &
            '       Upper output limit at    ',F8.2,' m  (GP ',I4,')'/)
339 FORMAT ('       No output during initial ',F8.2,' s')
342 FORMAT (/'       ',A2,'-cross-section  Arrays: ',A/ &
            '       Output every             ',F8.2,' s  ',A/ &
            '       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s'/ &
            '       Cross sections at ',A1,' = ',A/ &
            '       scalar-coordinates:   ',A,' m'/)
343 FORMAT (/'       Arrays: ',A/ &
            '       Output every             ',F8.2,' s  ',A/ &
            '       Time averaged over       ',F8.2,' s'/ &
            '       Averaging input every    ',F8.2,' s'/ &
            '       Upper output limit at    ',F8.2,' m  (GP ',I4,')'/)
344 FORMAT ('       Output format: ',A/)
352 FORMAT  (/'       Number of output time levels allowed: ',I3 /)
353 FORMAT  (/'       Number of output time levels allowed: unlimited' /)
354 FORMAT ('       Output format: ',A, '   compressed with level: ',I1/)
400 FORMAT (//' Physical quantities:'/ &
              ' -------------------'/)
410 FORMAT ('    Geograph. latitude  :   latitude  = ',F4.1,' degr'/   &
            '    Geograph. longitude :   longitude = ',F4.1,' degr'/   &
            '    Angular velocity    :   omega  =',E10.3,' rad/s'/  &
            '    Coriolis parameter  :   f      = ',F9.6,' 1/s'/    &
            '                            f*     = ',F9.6,' 1/s')
411 FORMAT (/'    Gravity             :   g      = ',F4.1,' m/s**2')
420 FORMAT (/'    Characteristic levels of the initial temperature profile:'// &
            '       Height:        ',A,'  m'/ &
            '       Temperature:   ',A,'  K'/ &
            '       Gradient:      ',A,'  K/100m'/ &
            '       Gridpoint:     ',A)
425 FORMAT (/'    Characteristic levels of the initial salinity profile:'// &
            '       Height:     ',A,'  m'/ &
            '       Salinity:   ',A,'  psu'/ &
            '       Gradient:   ',A,'  psu/100m'/ &
            '       Gridpoint:  ',A)
426 FORMAT (/'    Reference potential density for buoyancy:',F8.2,' kg/m**3')
427 FORMAT (/'    Initial wind profiles (u,v) are interpolated from given'// &
                  ' profiles')
450 FORMAT (//' LES / Turbulence quantities:'/ &
              ' ---------------------------'/)
453 FORMAT ('    Mixing length is limited to',F5.2,' * z')
454 FORMAT ('    TKE is not allowed to fall below ',E9.2,' (m/s)**2')
455 FORMAT ('    initial TKE is prescribed as ',E9.2,' (m/s)**2')
456 FORMAT ('    Day of the year at model start :   day_init = ',I3             &
            /'    UTC time at model start        :   time_utc_init = ',F7.1' s')
470 FORMAT (//' Actions during the simulation:'/ &
              ' -----------------------------'/)
471 FORMAT ('    Disturbance impulse (u,v) every :   ',F6.2,' s'/            &
            '    Disturbance amplitude           :    ',F5.2, ' m/s'/       &
            '    Lower disturbance level         : ',F8.2,' m (GP ',I4,')'/  &
            '    Upper disturbance level         : ',F8.2,' m (GP ',I4,')')
473 FORMAT ('    Disturbances cease as soon as the disturbance energy exceeds',&
                 F6.3, ' m**2/s**2')
474 FORMAT ('    Random number generator used    : ',A/)
475 FORMAT ('    The surface temperature is increased (or decreased, ', &
                 'respectively, if'/ &
            '    the value is negative) by ',F5.2,' K at the beginning of the',&
                 ' 3D-simulation'/)
503 FORMAT (' --> Momentum advection via Wicker-Skamarock-Scheme 5th order')
504 FORMAT (' --> Scalar advection via Wicker-Skamarock-Scheme 5th order')

 END SUBROUTINE header
