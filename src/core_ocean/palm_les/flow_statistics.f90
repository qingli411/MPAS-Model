!> @file flow_statistics.f90
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
!
!
! Description:
! ------------
!> Compute average profiles and further average flow quantities for the different
!> user-defined (sub-)regions. The region indexed 0 is the total model domain.
!>
!> @note For simplicity, nzb_s_inner and nzb_diff_s_inner are being used as a
!>       lower vertical index for k-loops for all variables, although strictly
!>       speaking the k-loops would have to be split up according to the staggered
!>       grid. However, this implies no error since staggered velocity components
!>       are zero at the walls and inside buildings.
!------------------------------------------------------------------------------!
 SUBROUTINE flow_statistics


    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw, e, heatflux_output_conversion, hyp, km, kh,         &
               momentumflux_output_conversion, nc, nr, p, prho, prr, pt, q,    &
               qc, ql, qr, rho_air, rho_air_zw, rho_ocean, s,                  &
               sa, salinityflux_output_conversion,                             &
               u, ug, v, vg, vpt, w, w_subs,                                   &
               zw, alpha_T, beta_S, solar3d, u_stk, v_stk, u_stk_zw, v_stk_zw

    USE control_parameters,                                                    &
        ONLY:   average_count_pr, cloud_droplets, cloud_physics, do_sum,       &
                dt_3d, g, humidity, initializing_actions, kappa, land_surface, &
                large_scale_forcing, large_scale_subsidence, max_pr_user,      &
                message_string, neutral, microphysics_morrison,                &
                microphysics_seifert, ocean, passive_scalar, simulated_time,   &
                simulated_time_at_begin, stokes_force,                         &
                use_subsidence_tendencies,            &
                use_surface_fluxes, use_top_fluxes, ws_scheme_mom,             &
                ws_scheme_sca, idealized_diurnal, top_momentumflux_u,          &
                top_momentumflux_v, top_heatflux, top_salinityflux

    USE cpulog,                                                                &
        ONLY:   cpu_log, log_point

    USE grid_variables,                                                        &
        ONLY:   ddx, ddy

    USE indices,                                                               &
        ONLY:   ngp_2dh, ngp_2dh_s_inner, ngp_3d, ngp_3d_inner, ngp_sums,      &
                ngp_sums_ls, nxl, nxr, nyn, nys, nzb, nzt, topo_min_level,     &
                wall_flags_0

    USE kinds

    USE netcdf_interface,                                                      &
        ONLY:  dots_rad, dots_soil, dots_max

    USE pegrid
    USE statistics

       USE surface_mod,                                                        &
          ONLY :  surf_def_h


    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  ki                  !<
    INTEGER(iwp) ::  k_surface_level     !<
    INTEGER(iwp) ::  m                   !< loop variable over all horizontal wall elements
    INTEGER(iwp) ::  l                   !< loop variable over surface facing -- up- or downward-facing
    INTEGER(iwp) ::  nt                  !<
    INTEGER(iwp) ::  omp_get_thread_num  !<
    INTEGER(iwp) ::  sr                  !<
    INTEGER(iwp) ::  tn                  !<
    INTEGER(iwp) :: top_bottom_flag

    LOGICAL ::  first  !<

    REAL(wp) ::  dptdz_threshold  !<
    REAL(wp) ::  fac              !<
    REAL(wp) ::  flag             !<
    REAL(wp) ::  height           !<
    REAL(wp) ::  pts              !<
    REAL(wp) ::  ust              !<
    REAL(wp) ::  ust2             !<
    REAL(wp) ::  u2               !<
    REAL(wp) ::  vst              !<
    REAL(wp) ::  vst2             !<
    REAL(wp) ::  v2               !<
    REAL(wp) ::  w2               !<

    REAL(wp) ::  dptdz(nzb+1:nzt+1)    !<
    REAL(wp) ::  sums_ll(nzb:nzt+1,2)  !<

    CALL cpu_log( log_point(10), 'flow_statistics', 'start' )


!
!-- To be on the safe side, check whether flow_statistics has already been
!-- called once after the current time step
    IF ( flow_statistics_called )  THEN

       message_string = 'flow_statistics is called two times within one ' // &
                        'timestep'
       CALL message( 'flow_statistics', 'PA0190', 1, 2, 0, 6, 0 )

    ENDIF
!
!
!--    Initialize (local) summation array
    sums_l = 0.0_wp
!changing sums_l
! 1 = u
! 2 = v
! 3 = w
! 4 = pt
! 5 = sa
! 6 = uw_res
! 7 = vw_res
! 8 = wpt_res
! 9 = wsa_res
!10 = uw_sgs
!11 = vw_sgs
!12 = wpt_sgs
!13 = wsa_sgs

    IF ( ws_scheme_mom )  THEN
       ! TODO: sums_us2_ws_l etc are currently not updated in advec_ws <20190814, Qing Li> !

!
!--    According to the Neumann bc for the horizontal velocity components,
!--    the corresponding fluxes has to satisfiy the same bc.
       sums_us2_ws_l(nzt+1,:) = sums_us2_ws_l(nzt,:)
       sums_vs2_ws_l(nzt+1,:) = sums_vs2_ws_l(nzt,:)

       DO  i = 0, threads_per_task-1
!
!--       Swap the turbulent quantities evaluated in advec_ws.
          sums_l(:,10,i) = sums_wsus_ws_l(:,i)                              &
                           * momentumflux_output_conversion ! w*u*
          sums_l(:,11,i) = sums_wsvs_ws_l(:,i)                              &
                           * momentumflux_output_conversion ! w*v*
       ENDDO

    ENDIF
    IF ( ws_scheme_sca )  THEN
       ! TODO: sums_wspts_ws_l etc are currently not updated in advec_ws <20190814, Qing Li> !

       DO i = 0, threads_per_task-1
          sums_l(:,12,i) = sums_wspts_ws_l(:,i)      &
                           * heatflux_output_conversion  ! w*pt*
          sums_l(:,13,i) = sums_wssas_ws_l(:,i)      &
                           * salinityflux_output_conversion ! w*sa*
       ENDDO

    ENDIF
!
!-- Horizontally averaged profiles of horizontal velocities and temperature.
!-- They must have been computed before, because they are already required
!-- for other horizontal averages.
    tn = 0
    !$OMP PARALLEL PRIVATE( i, j, k, tn, flag )
    !$ tn = omp_get_thread_num()
    !$OMP DO
    DO  i = nxl, nxr
       DO  j =  nys, nyn
          DO  k = nzb, nzt+1
             flag = MERGE( 1.0_wp, 0.0_wp, BTEST( wall_flags_0(k,j,i), 22 ) )
             sums_l(k,1,tn)  = sums_l(k,1,tn)  + u(k,j,i)  * rmask(j,i,0)  &
                                                           * flag
             sums_l(k,2,tn)  = sums_l(k,2,tn)  + v(k,j,i)  * rmask(j,i,0)  &
                                                           * flag
             sums_l(k,3,tn)  = sums_l(k,3,tn)  + w(k,j,i)  * rmask(j,i,0)  &
                                                           * flag
             sums_l(k,4,tn)  = sums_l(k,4,tn)  + pt(k,j,i) * rmask(j,i,0)  &
                                                           * flag
             sums_l(k,5,tn)  = sums_l(k,5,tn)  + sa(k,j,i) * rmask(j,i,0)  &
                                                           * flag
          ENDDO
       ENDDO
    ENDDO
!
!-- Summation of thread sums
    IF ( threads_per_task > 1 )  THEN
       DO  i = 1, threads_per_task-1
          sums_l(:,1,0) = sums_l(:,1,0) + sums_l(:,1,i)
          sums_l(:,2,0) = sums_l(:,2,0) + sums_l(:,2,i)
          sums_l(:,3,0) = sums_l(:,3,0) + sums_l(:,2,i)
          sums_l(:,4,0) = sums_l(:,4,0) + sums_l(:,4,i)
          sums_l(:,5,0) = sums_l(:,5,0) + sums_l(:,5,i)
       ENDDO
    ENDIF
!
!-- Horizontally averaged profiles of the vertical fluxes

    !$OMP DO
    DO  i = nxl, nxr
       DO  j = nys, nyn
!
!--       Subgridscale fluxes (without Prandtl layer from k=nzb,
!--       oterwise from k=nzb+1)
!--       NOTE: for simplicity, nzb_diff_s_inner is used below, although
!--       ----  strictly speaking the following k-loop would have to be
!--             split up according to the staggered grid.
!--             However, this implies no error since staggered velocity
!--             components are zero at the walls and inside buildings.
!--       Flag 23 is used to mask surface fluxes as well as model-top fluxes,
!--       which are added further below.
          DO  k = nzb, nzt
             flag = MERGE( 1.0_wp, 0.0_wp,                                  &
                           BTEST( wall_flags_0(k,j,i), 23 ) ) *             &
                    MERGE( 1.0_wp, 0.0_wp,                                  &
                           BTEST( wall_flags_0(k,j,i), 9  ) )
!
!--          Momentum flux w"u"
             sums_l(k,10,tn) = sums_l(k,10,tn) - 0.25_wp * (                &
                            km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) &
                                                        ) * (               &
                                ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)     &
                              + ( w(k,j,i)   - w(k,j,i-1) ) * ddx           &
                                                        ) * rmask(j,i,0)   &
                                      * rho_air_zw(k)                       &
                                      * momentumflux_output_conversion(k)   &
                                      * flag
!
!--          Momentum flux w"v"
             sums_l(k,11,tn) = sums_l(k,11,tn) - 0.25_wp * (                &
                            km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) &
                                                        ) * (               &
                                ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)     &
                              + ( w(k,j,i)   - w(k,j-1,i) ) * ddy           &
                                                        ) * rmask(j,i,0)   &
                                      * rho_air_zw(k)                       &
                                      * momentumflux_output_conversion(k)   &
                                      * flag
!
!--          Heat flux w"pt"
             sums_l(k,12,tn) = sums_l(k,12,tn)                              &
                                      - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )&
                                            * ( pt(k+1,j,i) - pt(k,j,i) )   &
                                            * rho_air_zw(k)                 &
                                            * heatflux_output_conversion(k) &
                                            * ddzu(k+1) * rmask(j,i,0)     &
                                            * flag


!
!--          Salinity flux w"sa"
             sums_l(k,13,tn) = sums_l(k,13,tn)                           &
                                      - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )&
                                            * ( sa(k+1,j,i) - sa(k,j,i) )   &
                                            * rho_air_zw(k)                 &
                                            * salinityflux_output_conversion(k) &
                                            * ddzu(k+1) * rmask(j,i,0)     &
                                            * flag

             sums_l(k,14,tn) = sums_l(k,14,tn) + rho_ocean(k,j,i) *      &
                                                    rmask(j,i,0) * flag
          ENDDO

!--       Subgridscale fluxes at the top surface
          IF ( use_top_fluxes )  THEN
!             m = surf_def_h(2)%start_index(j,i)
             sums_l(nzt:nzt+1,10,tn) = sums_l(nzt:nzt+1,10,tn) + &
                                 momentumflux_output_conversion(nzt:nzt+1) * &
                                 top_momentumflux_u * rmask(j,i,0)    ! w"u"
             sums_l(nzt:nzt+1,11,tn) = sums_l(nzt:nzt+1,11,tn) + &
                                 momentumflux_output_conversion(nzt:nzt+1) * &
                                 top_momentumflux_v * rmask(j,i,0)    ! w"v"
             sums_l(nzt:nzt+1,12,tn) = sums_l(nzt:nzt+1,12,tn) + &
                                 heatflux_output_conversion(nzt:nzt+1) * &
                                 top_heatflux  * rmask(j,i,0)   ! w"pt"
             sums_l(nzt:nzt+1,13,tn) = sums_l(nzt:nzt+1,13,tn) + &
                                 salinityflux_output_conversion(nzt:nzt+1) * &
                                 top_salinityflux * rmask(j,i,0)  ! w"sa"
          ENDIF

       ENDDO
    ENDDO

    sums(:,:) = sums_l(:,:,0) / ngp_2dh(0)
    hom(:,1,:,0) = sums(:,:)

!
!-- If required, sum up horizontal averages for subsequent time averaging.
!-- Do not sum, if flow statistics is called before the first initial time step.
    IF ( do_sum  .AND.  simulated_time /= 0.0_wp )  THEN
       IF ( average_count_pr == 0 )  hom_sum = 0.0_wp
       hom_sum = hom_sum + hom(:,1,:,:)
       average_count_pr = average_count_pr + 1
       do_sum = .FALSE.
    ENDIF
!
!-- Set flag for other UPs (e.g. output routines, but also buoyancy).
!-- This flag is reset after each time step in time_integration.
    flow_statistics_called = .TRUE.

    CALL cpu_log( log_point(10), 'flow_statistics', 'stop' )


 END SUBROUTINE flow_statistics
