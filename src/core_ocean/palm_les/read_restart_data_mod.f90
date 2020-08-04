!> @file read_restart_data_mod.f90
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
! $Id: read_restart_data_mod.f90 3065 2018-06-12 07:03:02Z Giersch $
! New parameters concerning vertical grid stretching have been added
!
! 3056 2018-06-04 07:49:35Z Giersch
! found variable has to be set to false inside overlap loop
!
! 3049 2018-05-29 13:52:36Z Giersch
! Error messages revised
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
!
! 3004 2018-04-27 12:33:25Z Giersch
! precipitation_rate_av removed
!
! 3003 2018-04-23 10:22:58Z Giersch
! z_i is also read to use the last known inversion height from the
! initial run as the first inversion height which is written into the
! run control file
!
! 2956 2018-04-10 11:01:03Z Giersch
! spectrum_x and spectrum_y have been moved to global data
!
! 2921 2018-03-22 15:05:23Z Giersch
! spinup_time, day_of_year_init and time_utc_init are also read now
!
! 2912 2018-03-20 13:00:05Z knoop
! Added gust module interface calls
!
! 2894 2018-03-15 09:17:58Z Giersch
! Initial revision
!
!
! Description:
! ------------
!> Reads restart data from restart-file(s) (binary format).
!------------------------------------------------------------------------------!
 MODULE read_restart_data_mod


    USE control_parameters


    IMPLICIT NONE


    INTERFACE rrd_global
       MODULE PROCEDURE rrd_global
    END INTERFACE rrd_global

    INTERFACE rrd_read_parts_of_global
       MODULE PROCEDURE rrd_read_parts_of_global
    END INTERFACE rrd_read_parts_of_global

    INTERFACE rrd_local
       MODULE PROCEDURE rrd_local
    END INTERFACE rrd_local

    INTERFACE rrd_skip_global
       MODULE PROCEDURE rrd_skip_global
    END INTERFACE rrd_skip_global


    PUBLIC rrd_global, rrd_read_parts_of_global, rrd_local,      &
           rrd_skip_global


 CONTAINS

! Description:
! ------------
!> Reads values of global control variables from restart-file (binary format)
!> created by PE0 of the previous run
!------------------------------------------------------------------------------!
    SUBROUTINE rrd_global


       USE arrays_3d,                                                          &
           ONLY:  pt_init, sa_init, u_init, v_init

       USE date_and_time_mod,                                                  &
           ONLY:  day_of_year_init, time_utc_init

           USE grid_variables,                                                     &
           ONLY:  dx, dy

           USE indices,                                                            &
           ONLY:  nz, nx, nx_on_file, ny, ny_on_file

           USE netcdf_interface,                                                   &
           ONLY:  netcdf_precision, output_for_t0

       USE pegrid

       USE statistics,                                                         &
           ONLY:  hom, hom_sum, pr_palm, u_max, u_max_ijk,  &
                  v_max, v_max_ijk, w_max, w_max_ijk

       IMPLICIT NONE

       CHARACTER (LEN=10) ::  binary_version_global, version_on_file

       LOGICAL ::  found


       CALL check_open( 13 )

!
!-- Make version number check first
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       READ ( 13 )  version_on_file
!
!-- Read number of PEs and horizontal index bounds of all PEs used in previous
!-- run
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( TRIM( restart_string(1:length) ) /= 'numprocs' )  THEN
          WRITE( message_string, * ) 'numprocs not found in data from prior ', &
                                     'run on PE ', myid
          CALL message( 'rrd_global', 'PA0297', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  numprocs_previous_run

       IF ( .NOT. ALLOCATED( hor_index_bounds_previous_run ) )  THEN
          ALLOCATE( hor_index_bounds_previous_run(4,0:numprocs_previous_run-1) )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'hor_index_bounds' )  THEN
          WRITE( message_string, * ) 'hor_index_bounds not found in data ',    &
                                     'from prior run on PE ', myid
          CALL message( 'rrd_global', 'PA0298', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  hor_index_bounds_previous_run

!
!-- Read vertical number of gridpoints and number of different areas used
!-- for computing statistics. Allocate arrays depending on these values,
!-- which are needed for the following read instructions.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'nz' )  THEN
          WRITE( message_string, * ) 'nz not found in data from prior run ',   &
                                     'on PE ', myid
          CALL message( 'rrd_global', 'PA0299', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  nz

       IF ( .NOT. ALLOCATED( u_init ) )  THEN
          ALLOCATE( u_init(0:nz+1), v_init(0:nz+1),                &
                    pt_init(0:nz+1), sa_init(0:nz+1),              &
                    hom(0:nz+1,2,pr_palm),     &
                    hom_sum(0:nz+1,pr_palm) )
       ENDIF

!
!-- Now read all control parameters:
!-- Caution: When the following read instructions have been changed, the
!-- -------  version number stored in the variable binary_version_global has to
!--          be increased. The same changes must also be done in
!--          wrd_write_global.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       DO WHILE ( restart_string(1:length) /= 'binary_version_local' )

          found = .FALSE.

          SELECT CASE ( restart_string(1:length) )

             CASE ( 'average_count_pr' )
                READ ( 13 )  average_count_pr
             CASE ( 'average_count_3d' )
                READ ( 13 )  average_count_3d
             CASE ( 'bc_e_b' )
                READ ( 13 )  bc_e_b
             CASE ( 'bc_p_b' )
                READ ( 13 )  bc_p_b
             CASE ( 'bc_p_t' )
                READ ( 13 )  bc_p_t
             CASE ( 'bc_pt_b' )
                READ ( 13 )  bc_pt_b
             CASE ( 'bc_pt_t' )
                READ ( 13 )  bc_pt_t
             CASE ( 'bc_pt_t_val' )
                READ ( 13 )  bc_pt_t_val
             CASE ( 'bc_sa_b' )
                READ ( 13 )  bc_sa_b
             CASE ( 'bc_sa_t' )
                READ ( 13 )  bc_sa_t
             CASE ( 'bc_uv_b' )
                READ ( 13 )  bc_uv_b
             CASE ( 'bc_uv_t' )
                READ ( 13 )  bc_uv_t
             CASE ( 'call_psolver_at_all_substeps' )
                READ ( 13 )  call_psolver_at_all_substeps
             CASE ( 'cfl_factor' )
                READ ( 13 )  cfl_factor
             CASE ( 'collective_wait' )
                READ ( 13 )  collective_wait
             CASE ( 'current_timestep_number' )
                READ ( 13 )  current_timestep_number
             CASE ( 'day_of_year_init' )
                READ ( 13 )  day_of_year_init
             CASE ( 'do3d_time_count' )
                READ ( 13 )  do3d_time_count
             CASE ( 'dp_external' )
                READ ( 13 )  dp_external
             CASE ( 'dp_level_b' )
                READ ( 13 )  dp_level_b
             CASE ( 'dp_smooth' )
                READ ( 13 )  dp_smooth
             CASE ( 'dpdxy' )
                READ ( 13 )  dpdxy
             CASE ( 'dt_3d' )
                READ ( 13 )  dt_3d
             CASE ( 'dx' )
                READ ( 13 )  dx
             CASE ( 'dy' )
                READ ( 13 )  dy
             CASE ( 'dz' )
                READ ( 13 )  dz
             CASE ( 'dz_max' )
                READ ( 13 )  dz_max
             CASE ( 'dz_stretch_factor' )
                READ ( 13 )  dz_stretch_factor
             CASE ( 'dz_stretch_factor_array' )
                READ ( 13 )  dz_stretch_factor_array
             CASE ( 'dz_stretch_level' )
                READ ( 13 )  dz_stretch_level
             CASE ( 'dz_stretch_level_end' )
                READ ( 13 )  dz_stretch_level_end
             CASE ( 'dz_stretch_level_start' )
                READ ( 13 )  dz_stretch_level_start
             CASE ( 'e_min' )
                READ ( 13 )  e_min
             CASE ( 'hom' )
                READ ( 13 )  hom
             CASE ( 'hom_sum' )
                READ ( 13 )  hom_sum
             CASE ( 'latitude' )
                READ ( 13 )  latitude
             CASE ( 'longitude' )
                READ ( 13 )  longitude
             CASE ( 'most_method' )
                READ ( 13 )  most_method
             CASE ( 'netcdf_precision' )
                READ ( 13 )  netcdf_precision
             CASE ( 'nx' )
                READ ( 13 )  nx
                nx_on_file = nx
             CASE ( 'ny' )
                READ ( 13 )  ny
                ny_on_file = ny
             CASE ( 'old_dt' )
                READ ( 13 )  old_dt
             CASE ( 'omega' )
                READ ( 13 )  omega
             CASE ( 'output_for_t0' )
                READ (13)    output_for_t0
             CASE ( 'prandtl_number' )
                READ ( 13 )  prandtl_number
             CASE ( 'pt_init' )
                READ ( 13 )  pt_init
             CASE ( 'pt_vertical_gradient' )
                READ ( 13 )  pt_vertical_gradient
             CASE ( 'pt_vertical_gradient_level' )
                READ ( 13 )  pt_vertical_gradient_level
             CASE ( 'pt_vertical_gradient_level_ind' )
                READ ( 13 )  pt_vertical_gradient_level_ind
             CASE ( 'random_generator' )
                READ ( 13 )  random_generator
             CASE ( 'rayleigh_damping_factor' )
                READ ( 13 )  rayleigh_damping_factor
             CASE ( 'rayleigh_damping_height' )
                READ ( 13 )  rayleigh_damping_height
             CASE ( 'runnr' )
                READ ( 13 )  runnr
             CASE ( 'sa_init' )
                READ ( 13 )  sa_init
             CASE ( 'sa_surface' )
                READ ( 13 )  sa_surface
             CASE ( 'sa_vertical_gradient' )
                READ ( 13 )  sa_vertical_gradient
             CASE ( 'sa_vertical_gradient_level' )
                READ ( 13 )  sa_vertical_gradient_level
             CASE ( 'simulated_time' )
                READ ( 13 )  simulated_time
             CASE ( 'surface_pressure' )
                READ ( 13 )  surface_pressure
             CASE ( 'time_disturb' )
                READ ( 13 )  time_disturb
             CASE ( 'time_do3d' )
                READ ( 13 )  time_do3d
             CASE ( 'time_do_av' )
                READ ( 13 )  time_do_av
             CASE ( 'time_do_sla' )
                READ ( 13 )  time_do_sla
             CASE ( 'time_dopr' )
                READ ( 13 )  time_dopr
             CASE ( 'time_dopr_av' )
                READ ( 13 )  time_dopr_av
             CASE ( 'time_dopr_listing' )
                READ ( 13 )  time_dopr_listing
             CASE ( 'time_dosp' )
                READ ( 13 )  time_dosp
             CASE ( 'time_restart' )
                READ ( 13 )  time_restart
             CASE ( 'time_run_control' )
                READ ( 13 )  time_run_control
             CASE ( 'time_since_reference_point' )
                READ ( 13 )  time_since_reference_point
             CASE ( 'time_utc_init' )
                READ ( 13 )  time_utc_init
             CASE ( 'top_heatflux' )
                READ ( 13 )  top_heatflux
             CASE ( 'top_momentumflux_u' )
                READ ( 13 )  top_momentumflux_u
             CASE ( 'top_momentumflux_v' )
                READ ( 13 )  top_momentumflux_v
             CASE ( 'top_salinityflux' )
                READ ( 13 )  top_salinityflux
             CASE ( 'tsc' )
                READ ( 13 )  tsc
             CASE ( 'turbulence_closure' )
                READ ( 13 )  turbulence_closure
             CASE ( 'u_init' )
                READ ( 13 )  u_init
             CASE ( 'u_max' )
                READ ( 13 )  u_max
             CASE ( 'u_max_ijk' )
                READ ( 13 )  u_max_ijk
             CASE ( 'v_init' )
                READ ( 13 )  v_init
             CASE ( 'v_max' )
                READ ( 13 )  v_max
             CASE ( 'v_max_ijk' )
                READ ( 13 )  v_max_ijk
             CASE ( 'w_max' )
                READ ( 13 )  w_max
             CASE ( 'w_max_ijk' )
                READ ( 13 )  w_max_ijk
             CASE ( 'wall_adjustment' )
                READ ( 13 )  wall_adjustment


             CASE DEFAULT

                IF ( .NOT. found )  THEN
                   WRITE( message_string, * ) 'unknown variable named "',      &
                                           restart_string(1:length),           &
                                          '" found in global data from ',      &
                                          'prior run on PE ', myid
                CALL message( 'rrd_global', 'PA0302', 1, 2, 0, 6, 0 )

                ENDIF

          END SELECT
!
!--       Read next string
          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)

       ENDDO


    CALL close_file( 13 )


    END SUBROUTINE rrd_global



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Skipping the global control variables from restart-file (binary format)
!> except some information needed when reading restart data from a previous
!> run which used a smaller total domain or/and a different domain decomposition
!> (initializing_actions  == 'cyclic_fill').
!------------------------------------------------------------------------------!

    SUBROUTINE rrd_read_parts_of_global


       USE indices,                                                               &
           ONLY:  nz, nx, nx_on_file, ny, ny_on_file

       USE kinds

       USE pegrid

       USE statistics,                                                            &
           ONLY:  hom, hom_sum, pr_palm, u_max, u_max_ijk,     &
                  v_max, v_max_ijk, w_max, w_max_ijk

       IMPLICIT NONE

       CHARACTER (LEN=10) ::  version_on_file
       CHARACTER (LEN=1)  ::  cdum

       INTEGER(iwp) ::  nz_on_file
       INTEGER(iwp) ::  tmp_mpru

       REAL(wp), DIMENSION(:,:),   ALLOCATABLE ::  hom_sum_on_file
       REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  hom_on_file


       CALL check_open( 13 )

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       READ ( 13 )  version_on_file

!
!-- Read number of PEs and horizontal index bounds of all PEs used in previous
!-- run
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'numprocs' )  THEN
          WRITE( message_string, * ) 'numprocs not found in data from prior ', &
                                     'run on PE ', myid
          CALL message( 'rrd_read_parts_of_global', 'PA0297', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  numprocs_previous_run

       IF ( .NOT. ALLOCATED( hor_index_bounds_previous_run ) )  THEN
          ALLOCATE( hor_index_bounds_previous_run(4,0:numprocs_previous_run-1) )
       ENDIF

       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'hor_index_bounds' )  THEN
          WRITE( message_string, * ) 'hor_index_bounds not found in data ',    &
                                     'from prior run on PE ', myid
          CALL message( 'rrd_read_parts_of_global', 'PA0298', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  hor_index_bounds_previous_run

!
!-- Read vertical number of gridpoints and number of different areas used
!-- for computing statistics. Allocate arrays depending on these values,
!-- which are needed for the following read instructions.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       IF ( restart_string(1:length) /= 'nz' )  THEN
          message_string = 'nz not found in restart data file'
          CALL message( 'rrd_read_parts_of_global', 'PA0303', 1, 2, 0, 6, 0 )
       ENDIF
       READ ( 13 )  nz_on_file
       IF ( nz_on_file /= nz )  THEN
          WRITE( message_string, * ) 'mismatch concerning number of ',         &
                                     'gridpoints along z:',                    &
                                     '&nz on file    = "', nz_on_file, '"',    &
                                     '&nz from run   = "', nz, '"'
          CALL message( 'rrd_read_parts_of_global', 'PA0304', 1, 2, 0, 6, 0 )
       ENDIF

!
!-- Now read and check some control parameters and skip the rest
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       DO  WHILE ( restart_string(1:length) /= 'binary_version_local' )

          SELECT CASE ( restart_string(1:length) )

             CASE ( 'average_count_pr' )
                READ ( 13 )  average_count_pr
                IF ( average_count_pr /= 0 )  THEN
                   WRITE( message_string, * ) 'inflow profiles not ',          &
                                  'temporally averaged. &Averaging will be ',  &
                                  'done now using', average_count_pr,          &
                                  ' samples.'
                   CALL message( 'rrd_read_parts_of_global', 'PA0309',         &
                                 0, 1, 0, 6, 0 )
                ENDIF

             CASE ( 'hom' )
                ALLOCATE( hom_on_file(0:nz+1,2,pr_palm) )
                READ ( 13 )  hom_on_file
                hom(:,:,1:pr_palm+tmp_mpru) =                         &
                             hom_on_file(:,:,1:pr_palm+tmp_mpru)
                DEALLOCATE( hom_on_file )

             CASE ( 'hom_sum' )
                ALLOCATE( hom_sum_on_file(0:nz+1,pr_palm) )
                READ ( 13 )  hom_sum_on_file
                hom_sum(:,1:pr_palm+tmp_mpru) =                       &
                             hom_sum_on_file(:,1:pr_palm+tmp_mpru)
                DEALLOCATE( hom_sum_on_file )

             CASE ( 'nx' )
                READ ( 13 )  nx_on_file

             CASE ( 'ny' )
                READ ( 13 )  ny_on_file

             CASE DEFAULT

                READ ( 13 )  cdum

          END SELECT

          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)

       ENDDO

!
!-- Calculate the temporal average of vertical profiles, if neccessary
    IF ( average_count_pr /= 0 )  THEN
       hom_sum = hom_sum / REAL( average_count_pr, KIND=wp )
    ENDIF


    CALL close_file( 13 )


    END SUBROUTINE rrd_read_parts_of_global


! Description:
! ------------
!> Reads processor specific data of variables and arrays from restart file
!> (binary format).
!------------------------------------------------------------------------------!
 SUBROUTINE rrd_local


    USE arrays_3d,                                                             &
        ONLY:  e, kh, km, p, pt, sa, u, v, w

    USE averaging

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point_s

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, nx_on_file, ny, nys, nysg, nyn, &
               nyng, ny_on_file, nzb, nzt

    USE kinds

    USE pegrid

    USE random_function_mod,                                                   &
        ONLY:  random_iv, random_iy

    USE random_generator_parallel,                                             &
        ONLY:  id_random_array, seq_random_array

    IMPLICIT NONE

    CHARACTER (LEN=7)  ::  myid_char_save
    CHARACTER (LEN=10) ::  binary_version_local
    CHARACTER (LEN=10) ::  version_on_file

    INTEGER(iwp) ::  files_to_be_opened  !<
    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  myid_on_file        !<
    INTEGER(iwp) ::  numprocs_on_file    !<
    INTEGER(iwp) ::  nxlc                !<
    INTEGER(iwp) ::  nxlf                !<
    INTEGER(iwp) ::  nxlpr               !<
    INTEGER(iwp) ::  nxl_on_file         !<
    INTEGER(iwp) ::  nxrc                !<
    INTEGER(iwp) ::  nxrf                !<
    INTEGER(iwp) ::  nxrpr               !<
    INTEGER(iwp) ::  nxr_on_file         !<
    INTEGER(iwp) ::  nync                !<
    INTEGER(iwp) ::  nynf                !<
    INTEGER(iwp) ::  nynpr               !<
    INTEGER(iwp) ::  nyn_on_file         !<
    INTEGER(iwp) ::  nysc                !<
    INTEGER(iwp) ::  nysf                !<
    INTEGER(iwp) ::  nyspr               !<
    INTEGER(iwp) ::  nys_on_file         !<
    INTEGER(iwp) ::  nzb_on_file         !<
    INTEGER(iwp) ::  nzt_on_file         !<
    INTEGER(iwp) ::  offset_x            !<
    INTEGER(iwp) ::  offset_y            !<
    INTEGER(iwp) ::  shift_x             !<
    INTEGER(iwp) ::  shift_y             !<

    INTEGER(iwp), DIMENSION(numprocs_previous_run) ::  file_list       !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run) ::  overlap_count   !<

    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxlfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nxrfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nynfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  nysfa      !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_xa  !<
    INTEGER(iwp), DIMENSION(numprocs_previous_run,1000) ::  offset_ya  !<

    INTEGER(isp), DIMENSION(:,:),   ALLOCATABLE ::  tmp_2d_id_random   !< temporary array for storing random generator data
    INTEGER(isp), DIMENSION(:,:,:), ALLOCATABLE ::  tmp_2d_seq_random  !< temporary array for storing random generator data

    LOGICAL ::  found

    REAL(wp) ::  rdummy

    REAL(wp), DIMENSION(:,:),   ALLOCATABLE   ::  tmp_2d      !< temporary array for storing 2D data
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE   ::  tmp_3d      !< temporary array for storing 3D data


!
!-- Read data from previous model run.
    CALL cpu_log( log_point_s(14), 'rrd_local', 'start' )

!
!-- Check which of the restart files contain data needed for the subdomain
!-- of this PE
    files_to_be_opened = 0

    DO  i = 1, numprocs_previous_run
!
!--    Store array bounds of the previous run ("pr") in temporary scalars
       nxlpr = hor_index_bounds_previous_run(1,i-1)
       nxrpr = hor_index_bounds_previous_run(2,i-1)
       nyspr = hor_index_bounds_previous_run(3,i-1)
       nynpr = hor_index_bounds_previous_run(4,i-1)

!
!--    Determine the offsets. They may be non-zero in case that the total domain
!--    on file is smaller than the current total domain.
       offset_x = ( nxl / ( nx_on_file + 1 ) ) * ( nx_on_file + 1 )
       offset_y = ( nys / ( ny_on_file + 1 ) ) * ( ny_on_file + 1 )

!
!--    Start with this offset and then check, if the subdomain on file
!--    matches another time(s) in the current subdomain by shifting it
!--    for nx_on_file+1, ny_on_file+1 respectively

       shift_y = 0
       j       = 0
       DO WHILE (  nyspr+shift_y <= nyn-offset_y )

          IF ( nynpr+shift_y >= nys-offset_y ) THEN

             shift_x = 0
             DO WHILE ( nxlpr+shift_x <= nxr-offset_x )

                IF ( nxrpr+shift_x >= nxl-offset_x ) THEN
                   j = j +1
                   IF ( j > 1000 )  THEN
!
!--                   Array bound exceeded
                      message_string = 'data from subdomain of previous' //    &
                                       ' run mapped more than 1000 times'
                      CALL message( 'rrd_local', 'PA0284', 2, 2, -1,           &
                                       6, 1 )
                   ENDIF

                   IF ( j == 1 )  THEN
                      files_to_be_opened = files_to_be_opened + 1
                      file_list(files_to_be_opened) = i-1
                   ENDIF

                   offset_xa(files_to_be_opened,j) = offset_x + shift_x
                   offset_ya(files_to_be_opened,j) = offset_y + shift_y
!
!--                Index bounds of overlapping data
                   nxlfa(files_to_be_opened,j) = MAX( nxl-offset_x-shift_x,    &
                                                      nxlpr )
                   nxrfa(files_to_be_opened,j) = MIN( nxr-offset_x-shift_x,    &
                                                      nxrpr )
                   nysfa(files_to_be_opened,j) = MAX( nys-offset_y-shift_y,    &
                                                      nyspr )
                   nynfa(files_to_be_opened,j) = MIN( nyn-offset_y-shift_y,    &
                                                      nynpr )

                ENDIF

                shift_x = shift_x + ( nx_on_file + 1 )
             ENDDO

          ENDIF

          shift_y = shift_y + ( ny_on_file + 1 )
       ENDDO

       IF ( j > 0 )  overlap_count(files_to_be_opened) = j

    ENDDO

!
!-- Save the id-string of the current process, since myid_char may now be used
!-- to open files created by PEs with other id.
    myid_char_save = myid_char

    IF ( files_to_be_opened /= 1  .OR.  numprocs /= numprocs_previous_run )    &
    THEN
       WRITE( message_string, * ) 'number of PEs or virtual PE-grid changed ', &
                        'in restart run & PE', myid, ' will read from files ', &
                         file_list(1:files_to_be_opened)
       CALL message( 'rrd_local', 'PA0285', 0, 0, 0, 6, 0 )
    ENDIF

!
!-- Read data from all restart files determined above
    DO  i = 1, files_to_be_opened

       j = file_list(i)
!
!--    Set the filename (underscore followed by four digit processor id)
       WRITE (myid_char,'(''_'',I6.6)')  j

!
!--    Open the restart file. If this file has been created by PE0 (_000000),
!--    the global variables at the beginning of the file have to be skipped
!--    first.
       CALL check_open( 13 )
       IF ( j == 0 )  CALL rrd_skip_global

!
!--    First compare the version numbers
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)
       READ ( 13 )  version_on_file
!
!--    Read number of processors, processor-id, and array ranges.
!--    Compare the array ranges with those stored in the index bound array.
       READ ( 13 )  numprocs_on_file, myid_on_file, nxl_on_file, nxr_on_file,  &
                    nys_on_file, nyn_on_file, nzb_on_file, nzt_on_file

       IF ( nxl_on_file /= hor_index_bounds_previous_run(1,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nxl on ',       &
                            'restart file "', myid_char, '"',                  &
                            '&nxl = ', nxl_on_file, ' but it should be',       &
                            '&= ', hor_index_bounds_previous_run(1,j),         &
                            '&from the index bound information array'
          CALL message( 'rrd_local', 'PA0287', 2, 2, -1, 6, 1 )
       ENDIF

       IF ( nxr_on_file /= hor_index_bounds_previous_run(2,j) )  THEN
           WRITE( message_string, * ) 'problem with index bound nxr on ',      &
                               'restart file "', myid_char, '"'  ,             &
                               ' nxr = ', nxr_on_file, ' but it should be',    &
                               ' = ', hor_index_bounds_previous_run(2,j),      &
                               ' from the index bound information array'
          CALL message( 'rrd_local', 'PA0288', 2, 2, -1, 6, 1 )

       ENDIF

       IF ( nys_on_file /= hor_index_bounds_previous_run(3,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nys on ',       &
                                 'restart file "', myid_char, '"',             &
                                 '&nys = ', nys_on_file, ' but it should be',  &
                                 '&= ', hor_index_bounds_previous_run(3,j),    &
                                 '&from the index bound information array'
          CALL message( 'rrd_local', 'PA0289', 2, 2, -1, 6, 1 )
       ENDIF

       IF ( nyn_on_file /= hor_index_bounds_previous_run(4,j) )  THEN
          WRITE( message_string, * ) 'problem with index bound nyn on ',       &
                               'restart file "', myid_char, '"',               &
                               '&nyn = ', nyn_on_file, ' but it should be',    &
                               '&= ', hor_index_bounds_previous_run(4,j),      &
                               '&from the index bound information array'
          CALL message( 'rrd_local', 'PA0290', 2, 2, -1, 6, 1 )
       ENDIF

       IF ( nzb_on_file /= nzb )  THEN
          WRITE( message_string, * ) 'mismatch between actual data and data ', &
                                     'from prior run on PE ', myid,            &
                                     '&nzb on file = ', nzb_on_file,           &
                                     '&nzb         = ', nzb
          CALL message( 'rrd_local', 'PA0291', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( nzt_on_file /= nzt )  THEN
          WRITE( message_string, * ) 'mismatch between actual data and data ', &
                                     'from prior run on PE ', myid,            &
                                     '&nzt on file = ', nzt_on_file,           &
                                     '&nzt         = ', nzt
          CALL message( 'rrd_local', 'PA0292', 1, 2, 0, 6, 0 )
       ENDIF

!
!--    Allocate temporary arrays sized as the arrays on the restart file
       ALLOCATE( tmp_2d(nys_on_file-nbgp:nyn_on_file+nbgp,                     &
                        nxl_on_file-nbgp:nxr_on_file+nbgp),                    &
                 tmp_3d(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,           &
                        nxl_on_file-nbgp:nxr_on_file+nbgp) )

!
!--    Read arrays
!--    ATTENTION: If the following read commands have been altered, the
!--    ---------- version number of the variable binary_version_local must
!--               be altered, too. Furthermore, the output list of arrays in
!--               wrd_write_local must also be altered
!--               accordingly.
       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)


!
!--    Loop over processor specific field data
       DO  WHILE ( restart_string(1:length) /= '*** end ***' )

!
!--       Map data on file as often as needed (data are read only for k=1)
          DO  k = 1, overlap_count(i)

             found = .FALSE.

!
!--          Get the index range of the subdomain on file which overlap with
!--          the current subdomain
             nxlf = nxlfa(i,k)
             nxlc = nxlfa(i,k) + offset_xa(i,k)
             nxrf = nxrfa(i,k)
             nxrc = nxrfa(i,k) + offset_xa(i,k)
             nysf = nysfa(i,k)
             nysc = nysfa(i,k) + offset_ya(i,k)
             nynf = nynfa(i,k)
             nync = nynfa(i,k) + offset_ya(i,k)


             SELECT CASE ( restart_string(1:length) )

                CASE ( 'e' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   e(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'e_av' )
                   IF ( .NOT. ALLOCATED( e_av ) )  THEN
                      ALLOCATE( e_av(nzb:nzt+1,nys-nbgp:nyn+nbgp,              &
                                     nxl-nbgp:nxr+nbgp) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   e_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'iran' ) ! matching random numbers is still unresolved
                                ! issue
                   IF ( k == 1 )  READ ( 13 )  iran

                CASE ( 'kh' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   kh(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'kh_av' )
                   IF ( .NOT. ALLOCATED( kh_av ) )  THEN
                      ALLOCATE( kh_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   kh_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'km' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   km(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'km_av' )
                   IF ( .NOT. ALLOCATED( km_av ) )  THEN
                      ALLOCATE( km_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg ))
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   km_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'p' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   p(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'p_av' )
                   IF ( .NOT. ALLOCATED( p_av ) )  THEN
                      ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   p_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pt' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pt(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'pt_av' )
                   IF ( .NOT. ALLOCATED( pt_av ) )  THEN
                      ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   pt_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'random_iv' )  ! still unresolved issue
                   IF ( k == 1 )  READ ( 13 )  random_iv
                   IF ( k == 1 )  READ ( 13 )  random_iy

                CASE ( 'seq_random_array' )
                   ALLOCATE( tmp_2d_id_random(nys_on_file:nyn_on_file,         &
                                              nxl_on_file:nxr_on_file) )
                   ALLOCATE( tmp_2d_seq_random(5,nys_on_file:nyn_on_file,      &
                                                 nxl_on_file:nxr_on_file) )
                   IF ( .NOT. ALLOCATED( id_random_array ) )  THEN
                      ALLOCATE( id_random_array(nys:nyn,nxl:nxr) )
                   ENDIF
                   IF ( .NOT. ALLOCATED( seq_random_array ) )  THEN
                      ALLOCATE( seq_random_array(5,nys:nyn,nxl:nxr) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d_id_random
                   IF ( k == 1 )  READ ( 13 )  tmp_2d_seq_random
                   id_random_array(nysc:nync,nxlc:nxrc) =                      &
                      tmp_2d_id_random(nysf:nynf,nxlf:nxrf)
                   seq_random_array(:,nysc:nync,nxlc:nxrc) =                   &
                      tmp_2d_seq_random(:,nysf:nynf,nxlf:nxrf)
                   DEALLOCATE( tmp_2d_id_random, tmp_2d_seq_random )

                CASE ( 'solar3d_av' )
                   IF ( .NOT. ALLOCATED( solar3d_av ) )  THEN
                      ALLOCATE( solar3d_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   solar3d_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)


                CASE ( 'rho_ocean_av' )
                   IF ( .NOT. ALLOCATED( rho_ocean_av ) )  THEN
                      ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   rho_ocean_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'alpha_T_av' )
                   IF ( .NOT. ALLOCATED( alpha_T_av ) )  THEN
                      ALLOCATE( alpha_T_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   alpha_T_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'beta_S_av' )
                   IF ( .NOT. ALLOCATED( beta_S_av ) )  THEN
                      ALLOCATE( beta_S_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   beta_S_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =   &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'sa' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   sa(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =             &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'sa_av' )
                   IF ( .NOT. ALLOCATED( sa_av ) )  THEN
                      ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   sa_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =          &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'shf_av' )
                   IF ( .NOT. ALLOCATED( shf_av ) )  THEN
                      ALLOCATE( shf_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   shf_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =          &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'shf_sol_av' )
                   IF ( .NOT. ALLOCATED( shf_sol_av ) )  THEN
                      ALLOCATE( shf_sol_av(nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_2d
                   shf_sol_av(nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp)  =          &
                      tmp_2d(nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'u' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   u(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'u_av' )
                   IF ( .NOT. ALLOCATED( u_av ) )  THEN
                      ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   u_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'v' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   v(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'v_av' )
                   IF ( .NOT. ALLOCATED( v_av ) )  THEN
                      ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   v_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'w' )
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   w(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =              &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE ( 'w_av' )
                   IF ( .NOT. ALLOCATED( w_av ) )  THEN
                      ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                   ENDIF
                   IF ( k == 1 )  READ ( 13 )  tmp_3d
                   w_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =           &
                      tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

                CASE DEFAULT

!
                   IF ( .NOT. found )  THEN
                      WRITE( message_string, * ) 'unknown variable named "',   &
                                                 restart_string(1:length),     &
                                                '" found in subdomain data ',  &
                                                'from prior run on PE ', myid
                      CALL message( 'rrd_local', 'PA0302', 1, 2, 0, 6, 0 )

                   ENDIF

             END SELECT

          ENDDO ! overlaploop

!
!--       Read next character string
          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)

       ENDDO ! dataloop

!
!--    Close the restart file
       CALL close_file( 13 )

       DEALLOCATE( tmp_2d, tmp_3d )

    ENDDO  ! loop over restart files


!
!-- Restore the original filename for the restart file to be written
    myid_char = myid_char_save

!
!-- End of time measuring for reading binary data
    CALL cpu_log( log_point_s(14), 'rrd_local', 'stop' )

 END SUBROUTINE rrd_local


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Skipping the global control variables from restart-file (binary format)
!------------------------------------------------------------------------------!

    SUBROUTINE rrd_skip_global


       IMPLICIT NONE

       CHARACTER (LEN=10) ::  version_on_file

       CHARACTER (LEN=1) ::  cdum


       READ ( 13 )  length
       READ ( 13 )  restart_string(1:length)

       DO  WHILE ( restart_string(1:length) /= 'binary_version_local' )

          READ ( 13 )  cdum
          READ ( 13 )  length
          READ ( 13 )  restart_string(1:length)

       ENDDO

       BACKSPACE ( 13 )
       BACKSPACE ( 13 )


    END SUBROUTINE rrd_skip_global


 END MODULE read_restart_data_mod
