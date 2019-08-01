!> @file time_integration.f90
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
! Description:
! ------------
!> Integration in time of the model equations, statistical analysis and graphic
!> output
!------------------------------------------------------------------------------!
 SUBROUTINE time_integration


    USE advec_ws,                                                              &
        ONLY:  ws_statistics, ws_finalize

    USE arrays_3d,                                                             &
        ONLY:  p, hyp, dzu, e, e_p, nc, nc_p, nr, nr_p, prho, pt, pt_p, pt_init, sa_init, &
               q_init, q, qc, qc_p, ql, ql_c, ql_v, ql_vp, qr, qr_p, q_p,      &
               ref_state, rho_ocean, s, s_p, sa_p, tend, u, u_p, v,            &
               v_p, w, w_p, alpha_T, beta_S, solar3d, sa, &
               ddzu, ddzw, dzw, dd2zu, drho_air, drho_air_zw,                  &
               rho_air, rho_air_zw, kh, km,                                    &
               te_m, tu_m, tv_m, tw_m, tpt_m, tsa_m,     &
               u_stk, v_stk, ug, vg, u_init, v_init, rdf, rdf_sc,              &
               ptdf_x, ptdf_y

    USE calc_mean_profile_mod,                                                 &
        ONLY:  calc_mean_profile

    USE control_parameters,                                                    &
        ONLY:  advected_distance_x, advected_distance_y,        &
               average_count_3d, averaging_interval, averaging_interval_pr,    &
               bc_lr_cyc, bc_ns_cyc, bc_pt_t_val,                              &
               bc_q_t_val, call_psolver_at_all_substeps,       &
               constant_flux_layer, constant_heatflux,          &
               create_disturbances, dopr_n, coupling_mode, &
               coupling_start_time, current_timestep_number,                   &
               disturbance_created, disturbance_energy_limit, dist_range,      &
               do_sum, old_dt, dt_3d, dt_averaging_input, dt_averaging_input_pr,       &
               dt_coupling, dt_data_output_av, dt_disturb, dt_do2d_xy,         &
               dt_do2d_xz, dt_do2d_yz, dt_do3d, dt_domask,dt_dopts, dt_dopr,   &
               dt_dopr_listing, dt_dots, dt_run_control, end_time,    &
               forcing, timestep_count,g,dp_smooth_factor,       &
               intermediate_timestep_count, intermediate_timestep_count_max,   &
               masks,                   &
               mid, average_count_meanpr,  &
               neutral, nr_timesteps_this_run, nudging,                        &
               passive_scalar, pt_reference,            &
               pt_slope_offset, random_heatflux,                    &
               run_coupled, simulated_time, simulated_time_chr,    &
               skip_time_do2d_xy, skip_time_do2d_xz, skip_time_do2d_yz,        &
               skip_time_do3d, skip_time_domask, skip_time_dopr,               &
               skip_time_data_output_av, sloping_surface,                      &
               stop_dt, terminate_coupled, terminate_run, timestep_scheme,     &
               time_coupling, time_do2d_xy, time_do2d_xz, time_do2d_yz,        &
               time_do3d, time_domask, time_dopr, time_dopr_av,                &
               time_dopr_listing, time_dopts, time_dosp, time_dosp_av,         &
               time_dots, time_do_av, time_do_sla, time_disturb, time_dvrp,    &
               time_run_control, time_since_reference_point, tsc,              &
               turbulent_inflow, turbulent_outflow, urban_surface,             &
               use_initial_profile_as_reference, dt_avg, time_avg,             &
               use_single_reference_value, uv_exposure, u_gtrans, v_gtrans,    &
               virtual_flight, wind_turbine, ws_scheme_mom, ws_scheme_sca,     &
               stokes_force, disturbFactor, uProfileInit, vProfileInit,        &
               tProfileInit, sProfileInit, dt_LS, dpdxy

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s

    USE indices,                                                               &
        ONLY:  nbgp, nx, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt, &
               wall_flags_0, advc_flags_1, advc_flags_2, ngp_2dh_outer

    USE interfaces

    USE kinds

    USE netcdf_data_input_mod,                                                 &
        ONLY:  force, netcdf_data_input_lsf

    USE pegrid

    USE progress_bar,                                                          &
        ONLY:  finish_progress_bar, output_progress_bar

    USE prognostic_equations_mod,                                              &
        ONLY:  prognostic_equations_vector

    USE statistics,                                                            &
        ONLY:  flow_statistics_called, hom, pr_palm, sums_ls_l, u_max,         &
               u_max_ijk, v_max, v_max_ijk, w_max, w_max_ijk,                  &
               statistic_regions, meanFields_avg, rmask

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  surface_layer_fluxes

    USE surface_mod,                                                           &
        ONLY:  surf_def_h, bc_h

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_diffusivities, production_e_init, &
               l_grid, l_wall

    USE stokes_force_mod,                                                      &
        ONLY:  stokes_pressure_head

    IMPLICIT NONE

    CHARACTER (LEN=9) ::  time_to_string          !<
    INTEGER(iwp)      ::  it
    INTEGER(iwp)      ::  lsp
    INTEGER(iwp)      ::  n

!    #define MY_DEBUG print *,"DEBUG",__LINE__,__FILE__

    REAL(wp) ::  dt_3d_old  !< temporary storage of timestep to be used for
                            !< steering of run control output interval
    REAL(wp) ::  tsrp_org   !< original value of time_since_reference_point

    logical :: first, firstav
!$acc data copyin( g ) &
!$acc      copyin( drho_air ) &
!$acc      copyin( drho_air_zw ) &
!$acc      copyin( rho_air ) &
!$acc      copyin( rho_air_zw ) &
!$acc      copyin( ref_state ) &
!$acc      copyin( dd2zu ) &
!$acc      copyin( ddzu ) &
!$acc      copyin( ddzw ) &
!$acc      copyin( dzw ) &
!$acc      copyin( bc_h ) &
!$acc      copyin( l_grid ) &
!$acc      copyin( l_wall ) &
!$acc      copyin( surf_def_h ) &
!!$acc      copyin( rmask ) &
!$acc      copyin( ngp_2dh_outer ) &
!$acc      copyin( wall_flags_0 ) &
!$acc      copyin( advc_flags_1, advc_flags_2 ) &
!$acc      copyin( dp_smooth_factor, dpdxy ) &
!$acc      copyin( ptdf_x, ptdf_y ) &
!$acc      copyin( hyp ) &
!$acc      copyin( tsc ) &
!$acc       copyin( rdf, rdf_sc ) &
!$acc      copyin( u_stk, v_stk ) &
!$acc      copyin( u_init, v_init ) &
!$acc      copyin( pt_init, sa_init ) &
!$acc      copyin( u, u_p, tu_m ) &
!$acc      copyin( v, v_p, tv_m ) &
!$acc      copyin( w, w_p, tw_m ) &
!$acc      copyin( e, e_p, te_m ) &
!$acc      copyin( pt, pt_p, tpt_m ) &
!$acc      copyin( sa, sa_p, tsa_m ) &
!$acc      copyin( kh, km ) &
!$acc      copyin( prho, rho_ocean ) &
!$acc      copyin( alpha_T, beta_S, solar3d ) &
!$acc      copyin( ug, vg )

!
!-- At beginning determine the first time step
    CALL timestep
!
!
!-- Determine and print out the run control quantities before the first time
!-- step of this run. For the initial run, some statistics (e.g. divergence)
!-- need to be determined first --> CALL flow_statistics at the beginning of
!-- run_control
    CALL run_control
!
!    CALL location_message( 'starting timestep-sequence', .TRUE. )

    disturbFactor = 1.0_wp

    firstav = .true.
    first = .true.
    !
!-- Start of the time loop
    DO  WHILE ( simulated_time < end_time  .AND.  .NOT. stop_dt  .AND. &
                .NOT. terminate_run )

!
       CALL cpu_log( log_point_s(10), 'timesteps', 'start' )
!
!--    Start of intermediate step loop
       intermediate_timestep_count = 0
       DO  WHILE ( intermediate_timestep_count < &
                   intermediate_timestep_count_max )

          intermediate_timestep_count = intermediate_timestep_count + 1

!
!--       Set the steering factors for the prognostic equations which depend
!--       on the timestep scheme
          CALL timestep_scheme_steering
          !$acc update device(tsc)

!--          Horizontally averaged profiles to be used as reference state in
!--          buoyancy terms (WARNING: only the respective last call of
!--          calc_mean_profile defines the reference state!)
                CALL calc_mean_profile( rho_ocean, 14 )

                ref_state(:)  = hom(:,1,14,0)
!--          Assure that ref_state does not become zero at any level
!--          ( might be the case if a vertical level is completely occupied
!--            with topography ).
             ref_state = MERGE( MAXVAL(ref_state), ref_state,                  &
                                ref_state == 0.0_wp )

!          CALL production_e_init
          IF ( ( ws_scheme_mom .OR. ws_scheme_sca )  .AND.  &
               intermediate_timestep_count == 1 )  CALL ws_statistics
!
            CALL prognostic_equations_vector
            !
!i

!--       Exchange of ghost points (lateral boundary conditions)
          CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'start' )

          CALL exchange_horiz( u_p, nbgp )
          CALL exchange_horiz( v_p, nbgp )
          CALL exchange_horiz( w_p, nbgp )
          CALL exchange_horiz( pt_p, nbgp )
          CALL exchange_horiz( e_p, nbgp )
             CALL exchange_horiz( sa_p, nbgp )
             CALL exchange_horiz( rho_ocean, nbgp )
             CALL exchange_horiz( prho, nbgp )
             CALL exchange_horiz( alpha_T, nbgp )
             CALL exchange_horiz( beta_S, nbgp )
             call exchange_horiz( solar3d, nbgp )
          CALL cpu_log( log_point(26), 'exchange-horiz-progn', 'stop' )

          !
!--       Boundary conditions for the prognostic quantities (except of the
!--       velocities at the outflow in case of a non-cyclic lateral wall)
          CALL boundary_conds
          !
!--       Swap the time levels in preparation for the next time step.
          CALL swap_timelevel

#ifndef __GPU
!
!--       Temperature offset must be imposed at cyclic boundaries in x-direction
!--       when a sloping surface is used
          IF ( sloping_surface )  THEN
             IF ( nxl ==  0 )  pt(:,:,nxlg:nxl-1) = pt(:,:,nxlg:nxl-1) - &
                                                    pt_slope_offset
             IF ( nxr == nx )  pt(:,:,nxr+1:nxrg) = pt(:,:,nxr+1:nxrg) + &
                                                    pt_slope_offset
          ENDIF
#endif
          !
!--       Impose a random perturbation on the horizontal velocity field
          IF ( create_disturbances  .AND.                                      &
               ( call_psolver_at_all_substeps  .AND.                           &
               intermediate_timestep_count == intermediate_timestep_count_max )&
          .OR. ( .NOT. call_psolver_at_all_substeps  .AND.                     &
               intermediate_timestep_count == 1 ) )                            &
          THEN
  time_disturb = time_disturb + dt_3d

             IF ( time_disturb < dt_disturb ) then
      !       if (dt_3d_old == dt_3d .and. first) then
                  CALL disturb_field( 'u', tend, u)
              CALL disturb_field( 'v', tend, v)
              !    call disturb_field('pt', tend, pt)
                 ! call disturb_field('sa', tend, sa, sLSforcing,  &
                 !        hom(:,1,23,statistic_regions) )
         else
                 first = .false.
            ENDIF
          ENDIF
          !$acc update self( pt, sa )

!--       Reduce the velocity divergence via the equation for perturbation
!--       pressure.
          CALL pres
          !$acc update self( u, v, w )
          !
!--       Compute the diffusion quantities

#ifndef __GPU
!--          First the vertical (and horizontal) fluxes in the surface
!--          (constant flux) layer are computed
             IF ( constant_flux_layer )  THEN
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'start' )
                CALL surface_layer_fluxes
                CALL cpu_log( log_point(19), 'surface_layer_fluxes', 'stop' )
             ENDIF
#endif

!--          Compute the diffusion coefficients
             CALL cpu_log( log_point(17), 'diffusivities', 'start' )
          CALL tcm_diffusivities
          !$acc update self(km, kh, e)
             CALL cpu_log( log_point(17), 'diffusivities', 'stop' )
!
       ENDDO   ! Intermediate step loop
       !
!--    Update perturbation pressure to account for the Stokes pressure
!--    head, if required
       IF ( stokes_force ) THEN
          CALL stokes_pressure_head
       ENDIF

!
!--    Increase simulation time and output times
       nr_timesteps_this_run      = nr_timesteps_this_run + 1
       current_timestep_number    = current_timestep_number + 1
       simulated_time             = simulated_time   + dt_3d
       time_since_reference_point = simulated_time - coupling_start_time
       simulated_time_chr         = time_to_string( time_since_reference_point )



       time_run_control   = time_run_control + dt_3d

!
!--    Check, if restart is necessary (because cpu-time is expiring or
!--    because it is forced by user) and set stop flag
!--    This call is skipped if the remote model has already initiated a restart.
!       IF ( .NOT. terminate_run )  CALL check_for_restart
!
!--    Set a flag indicating that so far no statistics have been created
!--    for this time step
       flow_statistics_called = .FALSE.

       CALL flow_statistics

       time_avg = time_avg + dt_3d
       If ( time_avg <= dt_avg ) THEN
         if(firstav) meanFields_avg = 0.0_wp

         meanFields_avg(:,1) = meanFields_avg(:,1) + hom(:,1,1,0)
         meanFields_avg(:,2) = meanFields_avg(:,2) + hom(:,1,2,0)
         meanFields_avg(:,3) = meanFields_avg(:,3) + hom(:,1,4,0)
         meanFields_avg(:,4) = meanFields_avg(:,4) + hom(:,1,5,0)

         average_count_meanpr = average_count_meanpr + 1
         firstav = .false.
       else

         meanFields_avg(:,1) = meanFields_avg(:,1) / average_count_meanpr
         meanFields_avg(:,2) = meanFields_avg(:,2) / average_count_meanpr
         meanFields_avg(:,3) = meanFields_avg(:,3) / average_count_meanpr
         meanFields_avg(:,4) = meanFields_avg(:,4) / average_count_meanpr

         average_count_meanpr = 0
         time_avg = 0.0_wp
         firstav = .true.
       endif

!--    Determine size of next time step. Save timestep dt_3d because it is
!--    newly calculated in routine timestep, but required further below for
!--    steering the run control output interval
       dt_3d_old = dt_3d
       CALL timestep

       if(dt_3d_old .ne. dt_3d .and. first) THEN
         ! disturbFactor = 0.0_wp
!          first = .false.
          uProfileInit(:) =  hom(nzb:nzt,1,1,0)
          vProfileInit(:) =  hom(nzb:nzt,1,2,0)
          tProfileInit(:) =  hom(nzb:nzt,1,4,0)
          sProfileInit(:) =  hom(nzb:nzt,1,5,0)
          end_time = simulated_time + dt_LS
       endif

!--    Computation and output of run control parameters.
!--    This is also done whenever perturbations have been imposed
       IF ( time_run_control >= dt_run_control  .OR.                     &
            timestep_scheme(1:5) /= 'runge'  .OR.  disturbance_created ) &
       THEN
          CALL run_control
          IF ( time_run_control >= dt_run_control )  THEN
             time_run_control = MOD( time_run_control, &
                                     MAX( dt_run_control, dt_3d_old ) )
          ENDIF
       ENDIF

!
!--    Output elapsed simulated time in form of a progress bar on stdout
!       IF ( myid == 0 )  CALL output_progress_bar

       CALL cpu_log( log_point_s(10), 'timesteps', 'stop' )


    ENDDO   ! time loop

!$acc end data

   IF ( myid == 0 )  CALL finish_progress_bar
   ! CALL location_message( 'finished time-stepping', .TRUE. )

 END SUBROUTINE time_integration
