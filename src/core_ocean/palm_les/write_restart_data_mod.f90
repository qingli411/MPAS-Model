!> @file write_restart_data_mod.f90
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
!
! Former revisions:
! -----------------
! $Id: write_restart_data_mod.f90 3065 2018-06-12 07:03:02Z Giersch $
! New parameters concerning vertical grid stretching have been added
!
! 3004 2018-04-27 12:33:25Z Giersch
! precipitation_rate_av removed
!
! 3003 2018-04-23 10:22:58Z Giersch
! z_i is also written out to use the last known inversion height from the
! initial run as the first inversion height which is written into the
! run control file
!
! 2956 2018-04-10 11:01:03Z Giersch
! spectrum_x and spectrum_y have been moved to global data
!
! 2921 2018-03-22 15:05:23Z Giersch
! spinup_time, day_of_year_init and time_utc_init are also written out now
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
!> Writes restart data into binary file(s) for restart runs.
!------------------------------------------------------------------------------!
 MODULE write_restart_data_mod


    USE control_parameters

    USE kinds

    USE pegrid,                                                                &
        ONLY:  myid, numprocs


    IMPLICIT NONE


    INTERFACE wrd_global
       MODULE PROCEDURE wrd_global
    END INTERFACE wrd_global

    INTERFACE wrd_local
       MODULE PROCEDURE wrd_local
    END INTERFACE wrd_local


    PUBLIC wrd_local, wrd_global


 CONTAINS


! Description:
! ------------
!> Global data of control variables and arrays is written out for
!> restarts (binary format).
!> This information is only written to the file opened by PE0.
!------------------------------------------------------------------------------!
    SUBROUTINE wrd_global


       USE arrays_3d,                                                          &
           ONLY:  pt_init, sa_init, u_init, v_init

       USE date_and_time_mod,                                                  &
           ONLY:  day_of_year_init, time_utc_init

       USE grid_variables,                                                     &
           ONLY:  dx, dy

       USE indices,                                                            &
           ONLY:  nx, ny, nz

       USE netcdf_interface,                                                   &
           ONLY:  netcdf_precision, output_for_t0

       USE pegrid,                                                             &
           ONLY:  hor_index_bounds, collective_wait

       USE statistics,                                                         &
           ONLY:  hom, hom_sum, u_max, u_max_ijk, v_max,    &
                  v_max_ijk, w_max, w_max_ijk

       IMPLICIT NONE

       CHARACTER (LEN=10)  ::  binary_version_global   !<


       binary_version_global = '4.7'

       CALL wrd_write_string( 'binary_version_global' )
       WRITE ( 14 )  binary_version_global

       CALL wrd_write_string( 'numprocs' )
       WRITE ( 14 )  numprocs

       CALL wrd_write_string( 'hor_index_bounds' )
       WRITE ( 14 )  hor_index_bounds

       CALL wrd_write_string( 'nz' )
       WRITE ( 14 )  nz

!
!-- Caution: After changes in the following parameter-list, the
!-- -------  version number stored in the variable binary_version_global has to
!--          be increased. The same changes must also be done in the parameter-
!--          list in rrd_global.

       CALL wrd_write_string( 'average_count_pr' )
       WRITE ( 14 )  average_count_pr

       CALL wrd_write_string( 'average_count_3d' )
       WRITE ( 14 )  average_count_3d

       CALL wrd_write_string( 'bc_e_b' )
       WRITE ( 14 )  bc_e_b

       CALL wrd_write_string( 'bc_p_b' )
       WRITE ( 14 )  bc_p_b

       CALL wrd_write_string( 'bc_p_t' )
       WRITE ( 14 )  bc_p_t

       CALL wrd_write_string( 'bc_pt_b' )
       WRITE ( 14 )  bc_pt_b

       CALL wrd_write_string( 'bc_pt_t' )
       WRITE ( 14 )  bc_pt_t

       CALL wrd_write_string( 'bc_pt_t_val' )
       WRITE ( 14 )  bc_pt_t_val

       CALL wrd_write_string( 'bc_sa_b' )
       WRITE ( 14 )  bc_sa_b

       CALL wrd_write_string( 'bc_sa_t' )
       WRITE ( 14 )  bc_sa_t

       CALL wrd_write_string( 'bc_uv_b' )
       WRITE ( 14 )  bc_uv_b

       CALL wrd_write_string( 'bc_uv_t' )
       WRITE ( 14 )  bc_uv_t

       CALL wrd_write_string( 'call_psolver_at_all_substeps' )
       WRITE ( 14 )  call_psolver_at_all_substeps

       CALL wrd_write_string( 'cfl_factor' )
       WRITE ( 14 )  cfl_factor

       CALL wrd_write_string( 'collective_wait' )
       WRITE ( 14 )  collective_wait

       CALL wrd_write_string( 'current_timestep_number' )
       WRITE ( 14 )  current_timestep_number

       CALL wrd_write_string( 'day_of_year_init' )
       WRITE ( 14 )  day_of_year_init

       CALL wrd_write_string( 'do3d_time_count' )
       WRITE ( 14 )  do3d_time_count

       CALL wrd_write_string( 'dp_external' )
       WRITE ( 14 )  dp_external

       CALL wrd_write_string( 'dp_level_b' )
       WRITE ( 14 )  dp_level_b

       CALL wrd_write_string( 'dp_smooth' )
       WRITE ( 14 )  dp_smooth

       CALL wrd_write_string( 'dpdxy' )
       WRITE ( 14 )  dpdxy

       CALL wrd_write_string( 'dt_3d' )
       WRITE ( 14 )  dt_3d

       CALL wrd_write_string( 'dx' )
       WRITE ( 14 )  dx

       CALL wrd_write_string( 'dy' )
       WRITE ( 14 )  dy

       CALL wrd_write_string( 'dz' )
       WRITE ( 14 )  dz

       CALL wrd_write_string( 'dz_max' )
       WRITE ( 14 )  dz_max

       CALL wrd_write_string( 'dz_stretch_factor' )
       WRITE ( 14 )  dz_stretch_factor

       CALL wrd_write_string( 'dz_stretch_factor_array' )
       WRITE ( 14 )  dz_stretch_factor_array

       CALL wrd_write_string( 'dz_stretch_level' )
       WRITE ( 14 )  dz_stretch_level

       CALL wrd_write_string( 'dz_stretch_level_end' )
       WRITE ( 14 )  dz_stretch_level_end

       CALL wrd_write_string( 'dz_stretch_level_start' )
       WRITE ( 14 )  dz_stretch_level_start

       CALL wrd_write_string( 'e_min' )
       WRITE ( 14 )  e_min

       CALL wrd_write_string( 'hom' )
       WRITE ( 14 )  hom

       CALL wrd_write_string( 'hom_sum' )
       WRITE ( 14 )  hom_sum

       CALL wrd_write_string( 'latitude' )
       WRITE ( 14 )  latitude

       CALL wrd_write_string( 'longitude' )
       WRITE ( 14 )  longitude

       CALL wrd_write_string( 'most_method' )
       WRITE ( 14 )  most_method

       CALL wrd_write_string( 'netcdf_precision' )
       WRITE ( 14 )  netcdf_precision

       CALL wrd_write_string( 'nx' )
       WRITE ( 14 )  nx

       CALL wrd_write_string( 'ny' )
       WRITE ( 14 )  ny

       CALL wrd_write_string( 'old_dt' )
       WRITE ( 14 )  old_dt

       CALL wrd_write_string( 'omega' )
       WRITE ( 14 )  omega

       CALL wrd_write_string( 'output_for_t0' )
       WRITE ( 14 )  output_for_t0

       CALL wrd_write_string( 'prandtl_number' )
       WRITE ( 14 )  prandtl_number

       CALL wrd_write_string( 'pt_init' )
       WRITE ( 14 )  pt_init

       CALL wrd_write_string( 'pt_vertical_gradient' )
       WRITE ( 14 )  pt_vertical_gradient

       CALL wrd_write_string( 'pt_vertical_gradient_level' )
       WRITE ( 14 )  pt_vertical_gradient_level

       CALL wrd_write_string( 'pt_vertical_gradient_level_ind' )
       WRITE ( 14 )  pt_vertical_gradient_level_ind

       CALL wrd_write_string( 'random_generator' )
       WRITE ( 14 )  random_generator

       CALL wrd_write_string( 'rayleigh_damping_factor' )
       WRITE ( 14 )  rayleigh_damping_factor

       CALL wrd_write_string( 'rayleigh_damping_height' )
       WRITE ( 14 )  rayleigh_damping_height

       CALL wrd_write_string( 'runnr' )
       WRITE ( 14 )  runnr

       CALL wrd_write_string( 'sa_init' )
       WRITE ( 14 )  sa_init

       CALL wrd_write_string( 'sa_surface' )
       WRITE ( 14 )  sa_surface

       CALL wrd_write_string( 'sa_vertical_gradient' )
       WRITE ( 14 )  sa_vertical_gradient

       CALL wrd_write_string( 'sa_vertical_gradient_level' )
       WRITE ( 14 )  sa_vertical_gradient_level

       CALL wrd_write_string( 'simulated_time' )
       WRITE ( 14 )  simulated_time

       CALL wrd_write_string( 'surface_pressure' )
       WRITE ( 14 )  surface_pressure

       CALL wrd_write_string( 'time_disturb' )
       WRITE ( 14 )  time_disturb

       CALL wrd_write_string( 'time_do3d' )
       WRITE ( 14 )  time_do3d

       CALL wrd_write_string( 'time_do_av' )
       WRITE ( 14 )  time_do_av

       CALL wrd_write_string( 'time_do_sla' )
       WRITE ( 14 )  time_do_sla

       CALL wrd_write_string( 'time_dopr' )
       WRITE ( 14 )  time_dopr

       CALL wrd_write_string( 'time_dopr_av' )
       WRITE ( 14 )  time_dopr_av

       CALL wrd_write_string( 'time_dopr_listing' )
       WRITE ( 14 )  time_dopr_listing

       CALL wrd_write_string( 'time_dosp' )
       WRITE ( 14 )  time_dosp

       CALL wrd_write_string( 'time_restart' )
       WRITE ( 14 )  time_restart

       CALL wrd_write_string( 'time_run_control' )
       WRITE ( 14 )  time_run_control

       CALL wrd_write_string( 'time_since_reference_point' )
       WRITE ( 14 )  time_since_reference_point

       CALL wrd_write_string( 'time_utc_init' )
       WRITE ( 14 )  time_utc_init

       CALL wrd_write_string( 'top_heatflux' )
       WRITE ( 14 )  top_heatflux

       CALL wrd_write_string( 'top_momentumflux_u' )
       WRITE ( 14 )  top_momentumflux_u

       CALL wrd_write_string( 'top_momentumflux_v' )
       WRITE ( 14 )  top_momentumflux_v

       CALL wrd_write_string( 'top_salinityflux' )
       WRITE ( 14 )  top_salinityflux

       CALL wrd_write_string( 'tsc' )
       WRITE ( 14 )  tsc

       CALL wrd_write_string( 'turbulence_closure' )
       WRITE ( 14 )  turbulence_closure

       CALL wrd_write_string( 'u_init' )
       WRITE ( 14 )  u_init

       CALL wrd_write_string( 'u_max' )
       WRITE ( 14 )  u_max

       CALL wrd_write_string( 'u_max_ijk' )
       WRITE ( 14 )  u_max_ijk

       CALL wrd_write_string( 'v_init' )
       WRITE ( 14 )  v_init

       CALL wrd_write_string( 'v_max' )
       WRITE ( 14 )  v_max

       CALL wrd_write_string( 'v_max_ijk' )
       WRITE ( 14 )  v_max_ijk

       CALL wrd_write_string( 'w_max' )
       WRITE ( 14 )  w_max

       CALL wrd_write_string( 'w_max_ijk' )
       WRITE ( 14 )  w_max_ijk

       CALL wrd_write_string( 'wall_adjustment' )
       WRITE ( 14 )  wall_adjustment


    END SUBROUTINE wrd_global


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Processor specific data of variables and arrays is written out for
!> restarts (binary format).
!> This information is written to the file opened by each PE.
!------------------------------------------------------------------------------!
    SUBROUTINE wrd_local


       USE arrays_3d,                                                          &
           ONLY:  e, kh, km, p, pt, sa, u, v, w

       USE averaging

      USE indices,                                                            &
           ONLY:  nxl, nxr, nys, nyn, nzb, nzt

       USE random_function_mod,                                                &
           ONLY:  random_iv, random_iy

       USE random_generator_parallel,                                          &
           ONLY:  id_random_array, seq_random_array

       IMPLICIT NONE

       CHARACTER (LEN=10) ::  binary_version_local   !<


!
!-- Write arrays.
       binary_version_local = '4.7'

       CALL wrd_write_string( 'binary_version_local' )
       WRITE ( 14 )  binary_version_local

       WRITE ( 14 )  numprocs, myid, nxl, nxr, nys, nyn, nzb, nzt

!
!-- Attention: After changes to the following output commands the version number
!-- ---------  of the variable binary_version_local must be changed!
!--            Also, the list of arrays to be read in rrd_local must be
!--            adjusted accordingly.
       CALL wrd_write_string( 'e' )
       WRITE ( 14 )  e

       IF ( ALLOCATED( e_av ) )  THEN
          CALL wrd_write_string( 'e_av' )
          WRITE ( 14 )  e_av
       ENDIF

       CALL wrd_write_string( 'iran' )
       WRITE ( 14 )  iran

       CALL wrd_write_string( 'kh' )
       WRITE ( 14 )  kh


       IF ( ALLOCATED( kh_av ) )  THEN
          CALL wrd_write_string( 'kh_av' )
          WRITE ( 14 )  kh_av
       ENDIF

       CALL wrd_write_string( 'km' )
       WRITE ( 14 )  km

       IF ( ALLOCATED( km_av ) )  THEN
          CALL wrd_write_string( 'km_av' )
          WRITE ( 14 )  km_av
       ENDIF

       CALL wrd_write_string( 'p' )
       WRITE ( 14 )  p

       IF ( ALLOCATED( p_av ) )  THEN
          CALL wrd_write_string( 'p_av' )
          WRITE ( 14 )  p_av
       ENDIF

       CALL wrd_write_string( 'pt' )
       WRITE ( 14 )  pt

       IF ( ALLOCATED( pt_av ) )  THEN
          CALL wrd_write_string( 'pt_av' )
          WRITE ( 14 )  pt_av
       ENDIF

       IF ( ALLOCATED( rho_ocean_av ) )  THEN
          CALL wrd_write_string( 'rho_ocean_av' )
          WRITE ( 14 )  rho_ocean_av
       ENDIF

       IF ( ALLOCATED( solar3d_av ) )  THEN
          CALL wrd_write_string( 'solar3d_av' )
          WRITE ( 14 )  solar3d_av
       ENDIF

       IF ( ALLOCATED( alpha_T_av ) )  THEN
          CALL wrd_write_string( 'alpha_T_av' )
          WRITE ( 14 )  alpha_T_av
       ENDIF

       IF ( ALLOCATED( beta_S_av ) )  THEN
          CALL wrd_write_string( 'beta_S_av' )
          WRITE ( 14 )  beta_S_av
       ENDIF

       CALL wrd_write_string( 'sa' )
       WRITE ( 14 )  sa

       IF ( ALLOCATED( sa_av ) )  THEN
          CALL wrd_write_string( 'sa_av' )
          WRITE ( 14 )  sa_av
       ENDIF

       CALL wrd_write_string( 'random_iv' )
       WRITE ( 14 )  random_iv
       WRITE ( 14 )  random_iy

       IF ( ALLOCATED( seq_random_array ) )  THEN
          CALL wrd_write_string( 'seq_random_array' )
          WRITE ( 14 )  id_random_array
          WRITE ( 14 )  seq_random_array
       ENDIF

       IF ( ALLOCATED( shf_av ) )  THEN
          CALL wrd_write_string( 'shf_av' )
          WRITE ( 14 )  shf_av
       ENDIF

       IF ( ALLOCATED( shf_sol_av ) )  THEN
          CALL wrd_write_string( 'shf_sol_av' )
          WRITE ( 14 )  shf_sol_av
       ENDIF

       CALL wrd_write_string( 'u' )
       WRITE ( 14 )  u

       IF ( ALLOCATED( u_av ) )  THEN
          CALL wrd_write_string( 'u_av' )
          WRITE ( 14 )  u_av
       ENDIF

       CALL wrd_write_string( 'v' )
       WRITE ( 14 )  v

       CALL wrd_write_string( 'w' )
       WRITE ( 14 )  w

!
!--    Write end label.
       CALL wrd_write_string( '*** end ***' )

    END SUBROUTINE wrd_local


 END MODULE write_restart_data_mod
