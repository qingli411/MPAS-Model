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
! Current revisions:
! -----------------
!
!
! Former revisions:
! -----------------
! $Id: flow_statistics.f90 3042 2018-05-25 10:44:37Z schwenkel $
! Changed the name specific humidity to mixing ratio
!
! 3040 2018-05-25 10:22:08Z schwenkel
! Comments related to the calculation of the inversion height expanded
!
! 3003 2018-04-23 10:22:58Z Giersch
! The inversion height will not be calcuated before the first timestep in
! case of restarts.
!
! 2968 2018-04-13 11:52:24Z suehring
! Bugfix in output of timeseries quantities in case of elevated model surfaces.
!
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
!
! 2773 2018-01-30 14:12:54Z suehring
! Timeseries output of surface temperature.
!
! 2753 2018-01-16 14:16:49Z suehring
! Tile approach for spectral albedo implemented.
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Bugfix in evaluation of surface quantities in case different surface types
! are used (MS)
!
! 2674 2017-12-07 11:49:21Z suehring
! Bugfix in output conversion of resolved-scale momentum fluxes in case of
! PW advections scheme.
!
! 2320 2017-07-21 12:47:43Z suehring
! Modularize large-scale forcing and nudging
!
! 2296 2017-06-28 07:53:56Z maronga
! Enabled output of radiation quantities for radiation_scheme = 'constant'
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2270 2017-06-09 12:18:47Z maronga
! Revised numbering (removed 2 timeseries)
!
! 2252 2017-06-07 09:35:37Z knoop
! perturbation pressure now depending on flux_output_mode
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC version of subroutine removed
!
! 2073 2016-11-30 14:34:05Z raasch
! openmp bugfix: large scale forcing calculations cannot be executed thread
! parallel
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
!
! 2026 2016-10-18 10:27:02Z suehring
! Bugfix, enable output of s*2.
! Change, calculation of domain-averaged perturbation energy.
! Some formatting adjustments.
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1976 2016-07-27 13:28:04Z maronga
! Removed some unneeded __rrtmg preprocessor directives
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
!
! 1918 2016-05-27 14:35:57Z raasch
! in case of Wicker-Skamarock scheme, calculate disturbance kinetic energy here,
! if flow_statistics is called before the first initial time step
!
! 1853 2016-04-11 09:00:35Z maronga
! Adjusted for use with radiation_scheme = constant
!
! 1849 2016-04-08 11:33:18Z hoffmann
! prr moved to arrays_3d
!
! 1822 2016-04-07 07:49:42Z hoffmann
! Output of bulk microphysics simplified.
!
! 1815 2016-04-06 13:49:59Z raasch
! cpp-directives for intel openmp bug removed
!
! 1783 2016-03-06 18:36:17Z raasch
! +module netcdf_interface
!
! 1747 2016-02-08 12:25:53Z raasch
! small bugfixes for accelerator version
!
! 1738 2015-12-18 13:56:05Z raasch
! bugfixes for calculations in statistical regions which do not contain grid
! points in the lowest vertical levels, mean surface level height considered
! in the calculation of the characteristic vertical velocity,
! old upstream parts removed
!
! 1709 2015-11-04 14:47:01Z maronga
! Updated output of Obukhov length
!
! 1691 2015-10-26 16:17:44Z maronga
! Revised calculation of Obukhov length. Added output of radiative heating >
! rates for RRTMG.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1658 2015-09-18 10:52:53Z raasch
! bugfix: temporary reduction variables in the openacc branch are now
! initialized to zero
!
! 1654 2015-09-17 09:20:17Z raasch
! FORTRAN bugfix of r1652
!
! 1652 2015-09-17 08:12:24Z raasch
! bugfix in calculation of energy production by turbulent transport of TKE
!
! 1593 2015-05-16 13:58:02Z raasch
! FORTRAN errors removed from openacc branch
!
! 1585 2015-04-30 07:05:52Z maronga
! Added output of timeseries and profiles for RRTMG
!
! 1571 2015-03-12 16:12:49Z maronga
! Bugfix: output of rad_net and rad_sw_in
!
! 1567 2015-03-10 17:57:55Z suehring
! Reverse modifications made for monotonic limiter.
!
! 1557 2015-03-05 16:43:04Z suehring
! Adjustments for monotonic limiter
!
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s.
!
! 1551 2015-03-03 14:18:16Z maronga
! Added suppport for land surface model and radiation model output.
!
! 1498 2014-12-03 14:09:51Z suehring
! Comments added
!
! 1482 2014-10-18 12:34:45Z raasch
! missing ngp_sums_ls added in accelerator version
!
! 1450 2014-08-21 07:31:51Z heinze
! bugfix: calculate fac only for simulated_time >= 0.0
!
! 1396 2014-05-06 13:37:41Z raasch
! bugfix: "copyin" replaced by "update device" in openacc-branch
!
! 1386 2014-05-05 13:55:30Z boeske
! bugfix: simulated time before the last timestep is needed for the correct
! calculation of the profiles of large scale forcing tendencies
!
! 1382 2014-04-30 12:15:41Z boeske
! Renamed variables which store large scale forcing tendencies
! pt_lsa -> td_lsa_lpt, pt_subs -> td_sub_lpt,
! q_lsa  -> td_lsa_q,   q_subs  -> td_sub_q,
! added Neumann boundary conditions for profile data output of large scale
! advection and subsidence terms at nzt+1
!
! 1374 2014-04-25 12:55:07Z raasch
! bugfix: syntax errors removed from openacc-branch
! missing variables added to ONLY-lists
!
! 1365 2014-04-22 15:03:56Z boeske
! Output of large scale advection, large scale subsidence and nudging tendencies
! +sums_ls_l, ngp_sums_ls, use_subsidence_tendencies
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
! 1299 2014-03-06 13:15:21Z heinze
! Output of large scale vertical velocity w_subs
!
! 1257 2013-11-08 15:18:40Z raasch
! openacc "end parallel" replaced by "end parallel loop"
!
! 1241 2013-10-30 11:36:58Z heinze
! Output of ug and vg
!
! 1221 2013-09-10 08:59:13Z raasch
! ported for openACC in separate #else branch
!
! 1179 2013-06-14 05:57:58Z raasch
! comment for profile 77 added
!
! 1115 2013-03-26 18:16:16Z hoffmann
! ql is calculated by calc_liquid_water_content
!
! 1111 2013-03-08 23:54:10Z raasch
! openACC directive added
!
! 1053 2012-11-13 17:11:03Z hoffmann
! additions for two-moment cloud physics scheme:
! +nr, qr, qc, prr
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1007 2012-09-19 14:30:36Z franke
! Calculation of buoyancy flux for humidity in case of WS-scheme is now using
! turbulent fluxes of WS-scheme
! Bugfix: Calculation of subgridscale buoyancy flux for humidity and cloud
! droplets at nzb and nzt added
!
! 801 2012-01-10 17:30:36Z suehring
! Calculation of turbulent fluxes in advec_ws is now thread-safe.
!
! Revision 1.1  1997/08/11 06:15:17  raasch
! Initial revision
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
        ONLY:  alpha_T, beta_S, dd2zu, ddzu, ddzw, e, km, kh, p, prho, pt,     &
               rho_ocean, sa, solar3d, u, v, w,                                &
               u_stk, v_stk, u_stk_zw, v_stk_zw, sgs_diss

    USE control_parameters,                                                    &
        ONLY:   average_count_pr, do_sum, message_string, simulated_time,      &
                stokes_force, top_heatflux,                                    &
                top_salinityflux, top_momentumflux_u, top_momentumflux_v

    USE cpulog,                                                                &
        ONLY:   cpu_log, log_point

    USE grid_variables,                                                        &
        ONLY:   ddx, ddy

    USE indices,                                                               &
        ONLY:   ngp_2dh, nxl, nxr, nyn, nys, nzb, nzt

    USE kinds

    USE pegrid

    USE statistics


    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !<
    INTEGER(iwp) ::  j                   !<
    INTEGER(iwp) ::  k                   !<
    INTEGER(iwp) ::  nt                  !<
    INTEGER(iwp) ::  omp_get_thread_num  !<
    INTEGER(iwp) ::  tn                  !<

    REAL(wp) ::  ust              !<
    REAL(wp) ::  vst              !<
    REAL(wp) ::  wpup             !<
    REAL(wp) ::  wpvp             !<

    REAL(wp) ::  dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz     !<
    REAL(wp) ::  def              !<

#ifdef __GPU
    REAL(wp) ::  tmp             !<
#endif


    CALL cpu_log( log_point(10), 'flow_statistics', 'start' )
!
!-- To be on the safe side, check whether flow_statistics has already been
!-- called once after the current time step
    IF ( flow_statistics_called )  THEN

       message_string = 'flow_statistics is called two times within one ' // &
                        'timestep'
       CALL message( 'flow_statistics', 'PA0190', 1, 2, 0, 6, 0 )

    ENDIF

    !$acc data create( sums_l )
!
!-- Initialize (local) summation array
    !$acc parallel loop collapse(3)
    DO k = nzb, nzt+1
       DO j = 1, pr_palm
          DO i = 0, threads_per_task-1
             sums_l(k,j,i) = 0.0_wp
          ENDDO
       ENDDO
    ENDDO
    !$acc end parallel

!
!-- Horizontally averaged profiles of horizontal velocities, temperature
!-- and salinity.
!-- They must have been computed before, because they are already required
!-- for other horizontal averages.
    tn = 0

#ifdef __GPU

    !$acc parallel present( sums_l, u, v, w, pt, sa, p, e )
    !$acc loop
    DO  k = nzb, nzt+1
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + u(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,1,tn) = tmp
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + v(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,2,tn) = tmp
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + w(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,3,tn) = tmp
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + pt(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,4,tn) = tmp
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + sa(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,5,tn) = tmp
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + p(k,j,i) - 2.0_wp * e(k,j,i) / 3.0_wp
          ENDDO
       ENDDO
       sums_l(k,6,tn) = tmp

    ENDDO
    !$acc end parallel

#else

    !$OMP PARALLEL PRIVATE( i, j, k, tn )
    !$ tn = omp_get_thread_num()
    !$OMP DO
    DO  i = nxl, nxr
       DO  j =  nys, nyn
          DO  k = nzb, nzt+1
             ! u
             sums_l(k,1,tn) = sums_l(k,1,tn) + u(k,j,i)
             ! v
             sums_l(k,2,tn) = sums_l(k,2,tn) + v(k,j,i)
             ! w
             sums_l(k,3,tn) = sums_l(k,3,tn) + w(k,j,i)
             ! pt
             sums_l(k,4,tn) = sums_l(k,4,tn) + pt(k,j,i)
             ! sa
             sums_l(k,5,tn) = sums_l(k,5,tn) + sa(k,j,i)
             ! p
             sums_l(k,6,tn) = sums_l(k,6,tn) + p(k,j,i) - 2.0_wp * e(k,j,i) / 3.0_wp
          ENDDO
       ENDDO
    ENDDO

#endif

!
!-- Summation of thread sums
    IF ( threads_per_task > 1 )  THEN
       !$acc parallel present( sums_l )
       !$acc loop
       DO  i = 1, threads_per_task-1
          !$acc loop collapse(2)
          DO k = nzb, nzt+1
             DO j = 1, 6
                sums_l(k,j,0) = sums_l(k,j,0) + sums_l(k,j,i)
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel
    ENDIF

#if defined( __parallel )
!
!-- Compute total sum from local sums
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( sums_l(nzb,1,0), sums(nzb,1), (nzt+2-nzb)*6, MPI_REAL,  &
                        MPI_SUM, comm2d, ierr )
#else

    !$acc parallel present( sums, sums_l )
    !$acc loop collapse(2)
    DO k = nzb, nzt+1
       DO j = 1, 6
          sums(k,j) = sums_l(k,j,0)
       ENDDO
    ENDDO
    !$acc end parallel

#endif

!
!-- Final values are obtained by division by the total number of grid points
!-- used for summation. After that store profiles.
    !$acc parallel present ( sums, hom, ngp_2dh )
    !$acc loop collapse(2)
    DO k = nzb, nzt+1
       DO j = 1, 6
          sums(k,j) = sums(k,j) / ngp_2dh
          hom(k,1,j) = sums(k,j)
       ENDDO
    ENDDO
    !$acc end parallel

!
!-- When calcuating horizontally-averaged total (resolved- plus subgrid-
!-- scale) vertical fluxes and velocity variances by using commonly-
!-- applied Reynolds-based methods ( e.g. <w'pt'> = (w-<w>)*(pt-<pt>) )
!-- in combination with the 5th order advection scheme, pronounced
!-- artificial kinks could be observed in the vertical profiles near the
!-- surface. Please note: these kinks were not related to the model truth,
!-- i.e. these kinks are just related to an evaluation problem.
!-- In order avoid these kinks, vertical fluxes and horizontal as well
!-- vertical velocity variances are calculated directly within the advection
!-- routines, according to the numerical discretization, to evaluate the
!-- statistical quantities as they will appear within the prognostic
!-- equations.
!-- Copy the turbulent quantities, evaluated in the advection routines to
!-- the local array sums_l() for further computations.

!
!-- According to the Neumann bc for the horizontal velocity components,
!-- the corresponding fluxes has to satisfiy the same bc.
    !$acc parallel present( sums_l, sums_wsus_ws_l, sums_wsvs_ws_l ) &
    !$acc present( sums_wspts_ws_l, sums_wssas_ws_l ) &
    !$acc present( sums_usvs_ws_l, sums_vsus_ws_l ) &
    !$acc present( sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l )
    !$acc loop
    DO  i = 0, threads_per_task-1
       sums_us2_ws_l(nzt+1,i) = sums_us2_ws_l(nzt,i)
       sums_vs2_ws_l(nzt+1,i) = sums_vs2_ws_l(nzt,i)
       !$acc loop
       DO k = nzb, nzt+1
          sums_l(k,11,i) = sums_wsus_ws_l(k,i)       ! w*u*
          sums_l(k,13,i) = sums_wsvs_ws_l(k,i)       ! w*v*
          sums_l(k,15,i) = sums_wspts_ws_l(k,i)      ! w*pt*
          sums_l(k,17,i) = sums_wssas_ws_l(k,i)      ! w*sa*
          sums_l(k,18,i) = sums_us2_ws_l(k,i)        ! u*2
          sums_l(k,19,i) = sums_vs2_ws_l(k,i)        ! v*2
          sums_l(k,20,i) = sums_ws2_ws_l(k,i)        ! w*2
          sums_l(k,35,i) = sums_usvs_ws_l(k,i)       ! u*v*
          sums_l(k,36,i) = sums_vsus_ws_l(k,i)       ! v*u*
       ENDDO
    ENDDO
    !$acc end parallel

!
!-- Horizontally averaged profiles of the remaining prognostic variables
!-- and variances
    tn = 0

#ifdef __GPU

    !$acc parallel present( sums_l, hom ) &
    !$acc present( u, v, w, pt, sa, p, e, km, kh, ddzu, dd2zu, ddzw ) &
    !$acc present( alpha_T, beta_S, solar3d, rho_ocean, prho )
    !$acc loop
    DO  k = nzb, nzt+1
       ! e
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + e(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,7,tn) = tmp

       ! km
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + km(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,8,tn) = tmp

       ! kh
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + kh(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,9,tn) = tmp

       ! pt*2
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + ( pt(k,j,i)-hom(k,1,4) )**2
          ENDDO
       ENDDO
       sums_l(k,21,tn) = tmp

       ! sa*2
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + ( sa(k,j,i)-hom(k,1,5) )**2
          ENDDO
       ENDDO
       sums_l(k,22,tn) = tmp

       ! w*3
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + w(k,j,i)**3
          ENDDO
       ENDDO
       sums_l(k,23,tn) = tmp

       ! rho_ocean
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + rho_ocean(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,26,tn) = tmp

       ! alpha_T
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + alpha_T(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,27,tn) = tmp

       ! beta_S
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + beta_S(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,28,tn) = tmp

       ! solar3d
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + solar3d(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,29,tn) = tmp

       ! prho
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + prho(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,30,tn) = tmp

       ! sgs dissipation
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + sgs_diss(k,j,i)
          ENDDO
       ENDDO
       sums_l(k,37,tn) = tmp

    ENDDO

    !$acc loop
    DO k = nzb+1, nzt
       ! w"u"
       tmp = 0.0_wp
       !$acc loop reduction(+:tmp)
       DO  j =  nys, nyn
          !$acc loop reduction(+:tmp)
          DO  i = nxl, nxr
             tmp = tmp - 0.25_wp * (                                           &
                            km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1)    &
                                                           ) * (               &
                                ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)        &
                              + ( w(k,j,i)   - w(k,j,i-1) ) * ddx              &
                                   )
          ENDDO
       ENDDO
       sums_l(k,10,tn) = tmp

       ! w"v"
       tmp = 0.0_wp
       !$acc loop reduction(+:tmp)
       DO  j =  nys, nyn
          !$acc loop reduction(+:tmp)
          DO  i = nxl, nxr
             tmp = tmp - 0.25_wp * (                                           &
                            km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i)    &
                                                           ) * (               &
                                ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)        &
                              + ( w(k,j,i)   - w(k,j-1,i) ) * ddy              &
                                   )
          ENDDO
       ENDDO
       sums_l(k,12,tn) = tmp

       ! w"pt"
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )      &
                                * ( pt(k+1,j,i) - pt(k,j,i) )      &
                                * ddzu(k+1)
          ENDDO
       ENDDO
       sums_l(k,14,tn) = tmp
       ! w"sa"
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )      &
                                * ( sa(k+1,j,i) - sa(k,j,i) )      &
                                * ddzu(k+1)
          ENDDO
       ENDDO
       sums_l(k,16,tn) = tmp

       ! transport of TKE due to pressure fluctuation (w*p*)
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + 0.5_wp * ( w(k,j,i) + w(k-1,j,i) -            &
                                    hom(k,1,3) - hom(k-1,1,3) ) *      &
                                  ( p(k,j,i) - hom(k,1,6) -            &
                                    2.0_wp * e(k,j,i) / 3.0_wp )
          ENDDO
       ENDDO
       sums_l(k,24,tn) = tmp

       ! transport of TKE (w*e*)
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + 0.25_wp * ( w(k,j,i) + w(k-1,j,i) -              &
                                     hom(k,1,3) - hom(k-1,1,3) ) *        &
                       0.25_wp * ( ( u(k,j,i) + u(k,j,i+1) -              &
                                     2.0_wp * hom(k,1,1) )**2 +           &
                                   ( v(k,j,i) + v(k,j+1,i) -              &
                                     2.0_wp * hom(k,1,2) )**2 +           &
                                   ( w(k,j,i) + w(k-1,j,i) -              &
                                     hom(k,1,3) - hom(k-1,1,3) )**2 )
          ENDDO
       ENDDO
       sums_l(k,25,tn) = tmp

       ! sgs flux of pressure and e
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + 2.0_wp * km(k,j,i) * ( e(k+1,j,i) - e(k-1,j,i) ) * dd2zu(k)
          ENDDO
       ENDDO
       sums_l(k,38,tn) = tmp

       ! transport of sgs fluxes (w"u"u*)
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             dudz = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -               &
                                u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)
             dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -               &
                                w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             dvdz = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -               &
                                v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)
             dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -               &
                                w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             tmp = tmp + km(k,j,i) * ( dudz + dwdx ) *                    &
                   ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - hom(k,1,1) ) +  &
                         km(k,j,i) * ( dvdz + dwdy ) *                    &
                   ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - hom(k,1,2) )
          ENDDO
       ENDDO
       sums_l(k,39,tn) = tmp

       ! dissipation of resolved TKE
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             dudx  = 0.5_wp   * ( u(k,j,i+1) + u(k+1,j,i+1) -           &
                                  u(k,j,i) - u(k+1,j,i) ) * ddx
             dudy  = 0.125_wp * ( u(k,j+1,i) + u(k+1,j+1,i) +           &
                                  u(k,j+1,i+1) + u(k+1,j+1,i+1) -       &
                                  u(k,j-1,i) - u(k+1,j-1,i) -           &
                                  u(k,j-1,i+1) - u(k+1,j-1,i+1) ) * ddy
             dudz  = 0.5_wp   * ( u(k+1,j,i) + u(k+1,j,i+1) -           &
                                  u(k,j,i) - u(k,j,i+1) ) * ddzu(k)
             dvdx  = 0.125_wp * ( v(k,j,i+1) + v(k+1,j,i+1) +           &
                                  v(k,j+1,i+1) + v(k+1,j+1,i+1) -       &
                                  v(k,j,i-1) - v(k+1,j,i-1) -           &
                                  v(k,j+1,i-1) - v(k+1,j+1,i-1) ) * ddx
             dvdy  = 0.5_wp   * ( v(k,j+1,i) + v(k+1,j+1,i) -           &
                                  v(k,j,i) - v(k+1,j,i) ) * ddy
             dvdz  = 0.5_wp   * ( v(k+1,j,i) + v(k+1,j+1,i) -           &
                                  v(k,j,i) - v(k,j+1,i) ) * ddzu(k)
             dwdx  = 0.5_wp   * ( w(k,j,i+1) - w(k,j,i-1) ) * ddx
             dwdy  = 0.5_wp   * ( w(k,j+1,i) - w(k,j-1,i) ) * ddy
             dwdz  = 0.5_wp   * ( w(k+1,j,i) - w(k-1,j,i) ) * ddzw(k)

             def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +           &
                              dudy**2 + dvdx**2 + dwdx**2 +             &
                              dwdy**2 + dudz**2 + dvdz**2 +             &
                   2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )
             tmp = tmp + 0.5_wp * ( km(k,j,i) + km(k+1,j,i) ) * def
          ENDDO
       ENDDO
       sums_l(k,40,tn) = tmp

       ! e*
       tmp = 0.0_wp
       !$acc loop collapse(2) reduction(+:tmp)
       DO  i = nxl, nxr
          DO  j =  nys, nyn
             tmp = tmp + 0.5_wp * (                                            &
                  ( 0.25_wp * ( u(k,j,i)+u(k+1,j,i)+u(k,j,i+1)+u(k+1,j,i+1) )  &
                            - 0.5_wp * ( hom(k,1,1) + hom(k+1,1,1) ) )**2      &
                + ( 0.25_wp * ( v(k,j,i)+v(k+1,j,i)+v(k,j+1,i)+v(k+1,j+1,i) )  &
                            - 0.5_wp * ( hom(k,1,2) + hom(k+1,1,2) ) )**2      &
                + ( w(k,j,i) - hom(k,1,3) )**2                )
          ENDDO
       ENDDO
       sums_l(k,41,tn) = tmp

    ENDDO
    !$acc end parallel

#else

    !$OMP PARALLEL PRIVATE( i, j, k, tn )
    !$ tn = omp_get_thread_num()
    !$OMP DO
    DO  i = nxl, nxr
       DO  j =  nys, nyn
          DO  k = nzb, nzt+1
             ! e
             sums_l(k,7,tn)  = sums_l(k,7,tn)  + e(k,j,i)
             ! km
             sums_l(k,8,tn)  = sums_l(k,8,tn)  + km(k,j,i)
             ! kh
             sums_l(k,9,tn)  = sums_l(k,9,tn)  + kh(k,j,i)
             ! pt*2
             sums_l(k,21,tn) = sums_l(k,21,tn) + ( pt(k,j,i)-hom(k,1,4) )**2
             ! sa*2
             sums_l(k,22,tn) = sums_l(k,22,tn) + ( sa(k,j,i)-hom(k,1,5) )**2
             ! w*3
             sums_l(k,23,tn) = sums_l(k,23,tn) + w(k,j,i)**3
             ! rho_ocean
             sums_l(k,26,tn) = sums_l(k,26,tn) + rho_ocean(k,j,i)
             ! alpha_T
             sums_l(k,27,tn) = sums_l(k,27,tn) + alpha_T(k,j,i)
             ! beta_S
             sums_l(k,28,tn) = sums_l(k,28,tn) + beta_S(k,j,i)
             ! solar3d
             sums_l(k,29,tn) = sums_l(k,29,tn) + solar3d(k,j,i)
             ! prho
             sums_l(k,30,tn) = sums_l(k,30,tn) + prho(k,j,i)
             ! sgs dissipation
             sums_l(k,37,tn) = sums_l(k,37,tn) + sgs_diss(k,j,i)
          ENDDO
       ENDDO
    ENDDO

    !$OMP PARALLEL PRIVATE( i, j, k, tn, ust, vst, wpup, wpvp )
    !$ tn = omp_get_thread_num()
    !$OMP DO
    DO  i = nxl, nxr
       DO  j =  nys, nyn
          DO  k = nzb+1, nzt
             ! w"u"
             wpup = - 0.25_wp * (                                              &
                            km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1)    &
                                ) * (                                          &
                                ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)        &
                              + ( w(k,j,i)   - w(k,j,i-1) ) * ddx              &
                                    )
             sums_l(k,10,tn) = sums_l(k,10,tn) + wpup
             ! w"v"
             wpvp = - 0.25_wp * (                                              &
                            km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i)    &
                                ) * (                                          &
                                ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)        &
                              + ( w(k,j,i)   - w(k,j-1,i) ) * ddy              &
                                    )
             sums_l(k,12,tn) = sums_l(k,12,tn) + wpvp
             ! w"pt"
             sums_l(k,14,tn) = sums_l(k,14,tn)                                 &
                                      - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )   &
                                            * ( pt(k+1,j,i) - pt(k,j,i) )      &
                                            * ddzu(k+1)
             ! w"sa"
             sums_l(k,16,tn) = sums_l(k,16,tn)                                 &
                                      - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )   &
                                            * ( sa(k+1,j,i) - sa(k,j,i) )      &
                                            * ddzu(k+1)

             ! transport of TKE due to pressure fluctuation (w*p*)
             sums_l(k,24,tn) = sums_l(k,24,tn) +                               &
                               0.5_wp * ( w(k,j,i) + w(k-1,j,i) -              &
                                          hom(k,1,3) - hom(k-1,1,3) ) *        &
                                        ( p(k,j,i) - hom(k,1,6) -              &
                                          2.0_wp * e(k,j,i) / 3.0_wp )

             ! transport of TKE (w*e*)
             sums_l(k,25,tn) = sums_l(k,25,tn) +                               &
                              0.25_wp * ( w(k,j,i) + w(k-1,j,i) -              &
                                          hom(k,1,3) - hom(k-1,1,3) ) *        &
                            0.25_wp * ( ( u(k,j,i) + u(k,j,i+1) -              &
                                          2.0_wp * hom(k,1,1) )**2 +           &
                                        ( v(k,j,i) + v(k,j+1,i) -              &
                                          2.0_wp * hom(k,1,2) )**2 +           &
                                        ( w(k,j,i) + w(k-1,j,i) -              &
                                          hom(k,1,3) - hom(k-1,1,3) )**2 )

             ! sgs transport of pressure and e fluxes
             sums_l(k,38,tn) = sums_l(k,38,tn) +                               &
                    2.0_wp * km(k,j,i) * ( e(k+1,j,i) - e(k-1,j,i) ) * dd2zu(k)

             ! transport of sgs fluxes
             dudz = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -                    &
                                u(k-1,j,i) - u(k-1,j,i+1) ) * dd2zu(k)
             dwdx = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -                    &
                                w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             dvdz = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -                    &
                                v(k-1,j,i) - v(k-1,j+1,i) ) * dd2zu(k)
             dwdy = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -                    &
                                w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             sums_l(k,39,tn) = sums_l(k,39,tn) +                               &
                               km(k,j,i) * ( dudz + dwdx ) *                   &
                        ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - hom(k,1,1) ) +  &
                               km(k,j,i) * ( dvdz + dwdy ) *                   &
                        ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - hom(k,1,2) )

             ! dissipation of resolved TKE
             dudx  = 0.5_wp   * ( u(k,j,i+1) + u(k+1,j,i+1) -                  &
                                  u(k,j,i) - u(k+1,j,i) ) * ddx
             dudy  = 0.125_wp * ( u(k,j+1,i) + u(k+1,j+1,i) +                  &
                                  u(k,j+1,i+1) + u(k+1,j+1,i+1) -              &
                                  u(k,j-1,i) - u(k+1,j-1,i) -                  &
                                  u(k,j-1,i+1) - u(k+1,j-1,i+1) ) * ddy
             dudz  = 0.5_wp   * ( u(k+1,j,i) + u(k+1,j,i+1) -                  &
                                  u(k,j,i) - u(k,j,i+1) ) * ddzu(k)
             dvdx  = 0.125_wp * ( v(k,j,i+1) + v(k+1,j,i+1) +                  &
                                  v(k,j+1,i+1) + v(k+1,j+1,i+1) -              &
                                  v(k,j,i-1) - v(k+1,j,i-1) -                  &
                                  v(k,j+1,i-1) - v(k+1,j+1,i-1) ) * ddx
             dvdy  = 0.5_wp   * ( v(k,j+1,i) + v(k+1,j+1,i) -                  &
                                  v(k,j,i) - v(k+1,j,i) ) * ddy
             dvdz  = 0.5_wp   * ( v(k+1,j,i) + v(k+1,j+1,i) -                  &
                                  v(k,j,i) - v(k,j+1,i) ) * ddzu(k)
             dwdx  = 0.5_wp   * ( w(k,j,i+1) - w(k,j,i-1) ) * ddx
             dwdy  = 0.5_wp   * ( w(k,j+1,i) - w(k,j-1,i) ) * ddy
             dwdz  = 0.5_wp   * ( w(k+1,j,i) - w(k-1,j,i) ) * ddzw(k)

             def = 2.0_wp * ( dudx**2 + dvdy**2 + dwdz**2 ) +                  &
                              dudy**2 + dvdx**2 + dwdx**2 +                    &
                              dwdy**2 + dudz**2 + dvdz**2 +                    &
                   2.0_wp * ( dvdx*dudy + dwdx*dudz + dwdy*dvdz )
             sums_l(k,40,tn) = sums_l(k,40,tn) +                               &
                               0.5_wp * ( km(k,j,i) + km(k+1,j,i) ) * def

             ! e*
             sums_l(k,41,tn) = sums_l(k,41,tn)                                 &
                                            + 0.5_wp * (                       &
                  ( 0.25_wp * ( u(k,j,i)+u(k+1,j,i)+u(k,j,i+1)+u(k+1,j,i+1) )  &
                            - 0.5_wp * ( hom(k,1,1) + hom(k+1,1,1) ) )**2      &
                + ( 0.25_wp * ( v(k,j,i)+v(k+1,j,i)+v(k,j+1,i)+v(k+1,j+1,i) )  &
                            - 0.5_wp * ( hom(k,1,2) + hom(k+1,1,2) ) )**2      &
                + ( w(k,j,i) - hom(k,1,3) )**2         )

          ENDDO
       ENDDO
    ENDDO

#endif

!
!-- Summation of thread sums
    IF ( threads_per_task > 1 )  THEN
       !$acc parallel present( sums_l )
       !$acc loop
       DO  i = 1, threads_per_task-1
          !$acc loop collapse(2)
          DO k = nzb, nzt+1
             DO j = 7, pr_palm
                sums_l(k,j,0) = sums_l(k,j,0) + sums_l(k,j,i)
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel
    ENDIF

#if defined( __parallel )
!
!-- Compute total sum from local sums
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( sums_l(nzb,7,0), sums(nzb,7), (nzt+2-nzb)*(pr_palm-6),     &
                        MPI_REAL, MPI_SUM, comm2d, ierr )
#else
    !$acc parallel present( sums, sums_l )
    !$acc loop collapse(2)
    DO k = nzb, nzt+1
       DO j = 7, pr_palm
          sums(k,j) = sums_l(k,j,0)
       ENDDO
    ENDDO
    !$acc end parallel
#endif

!
!-- Final values are obtained by division by the total number of grid points
!-- used for summation. After that store profiles.
!-- Check, if statistical regions do contain at least one grid point at the
!-- respective k-level, otherwise division by zero will lead to undefined
!-- values, which may cause e.g. problems with NetCDF output
!-- Profiles:
    !$acc parallel present( hom, sums, ngp_2dh )
    !$acc loop collapse(2)
    DO k = nzb, nzt+1
       DO j = 7, pr_palm
          sums(k,j) = sums(k,j) / ngp_2dh
          hom(k,1,j) = sums(k,j)
       ENDDO
    ENDDO
    !$acc end parallel

!
!-- Subgridscale fluxes at the top surface
    !$acc kernels present( hom )
    ! w"u"
    hom(nzt:nzt+1,1,10) = hom(nzt:nzt+1,1,10) + top_momentumflux_u
    ! w"v"
    hom(nzt:nzt+1,1,12) = hom(nzt:nzt+1,1,12) + top_momentumflux_v
    ! w"pt"
    hom(nzt:nzt+1,1,14) = hom(nzt:nzt+1,1,14) + top_heatflux
    ! w"sa"
    hom(nzt:nzt+1,1,16) = hom(nzt:nzt+1,1,16) + top_salinityflux
    !$acc end kernels
!
    IF ( stokes_force ) THEN
       ! save stokes drift
       !$acc parallel present( hom, u_stk, v_stk, u_stk_zw, v_stk_zw )
       !$acc loop
       DO k = nzb, nzt+1
          hom(k,1,31) = u_stk(k)
          hom(k,1,32) = v_stk(k)
          hom(k,1,33) = u_stk_zw(k)
          hom(k,1,34) = v_stk_zw(k)
       ENDDO
       !$acc end parallel
    ENDIF

!
!-- If required, sum up horizontal averages for subsequent time averaging.
!-- Do not sum, if flow statistics is called before the first initial time step.
    IF ( do_sum  .AND.  simulated_time /= 0.0_wp )  THEN
       IF ( average_count_pr == 0 ) THEN
          !$acc parallel present( hom_sum )
          !$acc loop collapse(2)
          DO k = nzb, nzt+1
             DO j = 1, pr_palm
                hom_sum(k,j) = 0.0_wp
             ENDDO
          ENDDO
          !$acc end parallel
       ENDIF
       !$acc parallel present( hom, hom_sum )
       !$acc loop collapse(2)
       DO k = nzb, nzt+1
          DO j = 1, pr_palm
             hom_sum(k,j) = hom_sum(k,j) + hom(k,1,j)
          ENDDO
       ENDDO
       !$acc end parallel
       average_count_pr = average_count_pr + 1
       do_sum = .FALSE.
    ENDIF

    !$acc end data

!
!-- Set flag for other UPs (e.g. output routines, but also buoyancy).
!-- This flag is reset after each time step in time_integration.
    flow_statistics_called = .TRUE.

    CALL cpu_log( log_point(10), 'flow_statistics', 'stop' )


 END SUBROUTINE flow_statistics
