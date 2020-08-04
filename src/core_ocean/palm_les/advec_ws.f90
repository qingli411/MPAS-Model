!> @file advec_ws.f90
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
! ------------------
!
!
! Former revisions:
! -----------------
! $Id: advec_ws.f90 3022 2018-05-18 11:12:35Z suehring $
! Bugfix in calculation of left-sided fluxes for u-component in OpenMP case.
!
! 2731 2018-01-09 17:44:02Z suehring
! Enable loop vectorization by splitting the k-loop
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implement advection for TKE-dissipation in case of RANS-TKE-e closure (TG)
! Allocate advc_flags_1/2 within ws_init_flags instead of init_grid
! Change argument list for exchange_horiz_2d_int (MS)
!
! 2329 2017-08-03 14:24:56Z knoop
! Bugfix concerning density in divergence correction close to buildings
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Rename wall_flags_0 and wall_flags_00 into advc_flags_1 and advc_flags_2,
! respectively.
! Set advc_flags_1/2 on basis of wall_flags_0/00 instead of nzb_s/u/v/w_inner.
! Setting advc_flags_1/2 also for downward-facing walls
!
! 2200 2017-04-11 11:37:51Z suehring
! monotonic_adjustment removed
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC version of subroutines removed
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1996 2016-08-18 11:42:29Z suehring
! Bugfix concerning calculation of turbulent of turbulent fluxes
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
!
! 1942 2016-06-14 12:18:18Z suehring
! Initialization of flags for ws-scheme moved from init_grid.
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
!
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme removed, microphysics_seifert added
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1630 2015-08-26 16:57:23Z suehring
!
!
! 1629 2015-08-26 16:56:11Z suehring
! Bugfix concerning wall_flags at left and south PE boundaries
!
! 1581 2015-04-10 13:45:59Z suehring
!
!
! 1580 2015-04-10 13:43:49Z suehring
! Bugfix: statistical evaluation of scalar fluxes in case of monotonic limiter
!
! 1567 2015-03-10 17:57:55Z suehring
! Bugfixes in monotonic limiter.
!
! 2015-03-09 13:10:37Z heinze
! Bugfix: REAL constants provided with KIND-attribute in call of
! intrinsic functions like MAX and MIN
!
! 1557 2015-03-05 16:43:04Z suehring
! Enable monotone advection for scalars using monotonic limiter
!
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
!
! 1361 2014-04-16 15:17:48Z hoffmann
! accelerator and vector version for qr and nr added
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute,
! module kinds added
! some formatting adjustments
!
! 1322 2014-03-20 16:38:49Z raasch
! REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1257 2013-11-08 15:18:40Z raasch
! accelerator loop directives removed
!
! 1221 2013-09-10 08:59:13Z raasch
! wall_flags_00 introduced, which holds bits 32-...
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1115 2013-03-26 18:16:16Z hoffmann
! calculation of qr and nr is restricted to precipitation
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary expansions according to the two new prognostic equations (nr, qr)
! of the two-moment cloud physics scheme:
! +flux_l_*, flux_s_*, diss_l_*, diss_s_*, sums_ws*s_ws_l
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1027 2012-10-15 17:18:39Z suehring
! Bugfix in calculation indices k_mm, k_pp in accelerator version
!
! 1019 2012-09-28 06:46:45Z raasch
! small change in comment lines
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator versions (*_acc) added
!
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 888 2012-04-20 15:03:46Z suehring
! Number of IBITS() calls with identical arguments is reduced.
!
! 862 2012-03-26 14:21:38Z suehring
! ws-scheme also work with topography in combination with vector version.
! ws-scheme also work with outflow boundaries in combination with
! vector version.
! Degradation of the applied order of scheme is now steered by multiplying with
! Integer advc_flags_1. 2nd order scheme, WS3 and WS5 are calculated on each
! grid point and mulitplied with the appropriate flag.
! 2nd order numerical dissipation term changed. Now the appropriate 2nd order
! term derived according to the 4th and 6th order terms is applied. It turns
! out that diss_2nd does not provide sufficient dissipation near walls.
! Therefore, the function diss_2nd is removed.
! Near walls a divergence correction is necessary to overcome numerical
! instabilities due to too less divergence reduction of the velocity field.
! boundary_flags and logicals steering the degradation are removed.
! Empty SUBROUTINE local_diss removed.
! Further formatting adjustments.
!
! 801 2012-01-10 17:30:36Z suehring
! Bugfix concerning OpenMP parallelization. Summation of sums_wsus_ws_l,
! sums_wsvs_ws_l, sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l, sums_wspts_ws_l,
! sums_wsqs_ws_l, sums_wssas_ws_l is now thread-safe by adding an additional
! dimension.
!
! Initial revision
!
! 411 2009-12-11 12:31:43 Z suehring
!
! Description:
! ------------
!> Advection scheme for scalars and momentum using the flux formulation of
!> Wicker and Skamarock 5th order. Additionally the module contains of a
!> routine using for initialisation and steering of the statical evaluation.
!> The computation of turbulent fluxes takes place inside the advection
!> routines.
!> Near non-cyclic boundaries the order of the applied advection scheme is
!> degraded.
!> A divergence correction is applied. It is necessary for topography, since
!> the divergence is not sufficiently reduced, resulting in erroneous fluxes and
!> partly numerical instabilities.
!-----------------------------------------------------------------------------!
 MODULE advec_ws



    PRIVATE
    PUBLIC   advec_s_ws, advec_u_ws, advec_v_ws, advec_w_ws, ws_init,          &
             ws_finalize, ws_init_flags, ws_statistics

    INTERFACE ws_init
       MODULE PROCEDURE ws_init
    END INTERFACE ws_init

    INTERFACE ws_finalize
       MODULE PROCEDURE ws_finalize
    END INTERFACE ws_finalize

    INTERFACE ws_init_flags
       MODULE PROCEDURE ws_init_flags
    END INTERFACE ws_init_flags

    INTERFACE ws_statistics
       MODULE PROCEDURE ws_statistics
    END INTERFACE ws_statistics

    INTERFACE advec_s_ws
       MODULE PROCEDURE advec_s_ws
    END INTERFACE advec_s_ws

    INTERFACE advec_u_ws
       MODULE PROCEDURE advec_u_ws
    END INTERFACE advec_u_ws

    INTERFACE advec_v_ws
       MODULE PROCEDURE advec_v_ws
    END INTERFACE advec_v_ws

    INTERFACE advec_w_ws
       MODULE PROCEDURE advec_w_ws
    END INTERFACE advec_w_ws

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of WS-scheme
!------------------------------------------------------------------------------!
    SUBROUTINE ws_init

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5, adv_sca_1, adv_sca_3,       &
                  adv_sca_5

       USE indices,                                                            &
           ONLY:  nyn, nys, nzb, nzt

       USE kinds

       USE pegrid

       USE statistics,                                                         &
           ONLY:  sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l,                 &
                  sums_usvs_ws_l, sums_vsus_ws_l,                              &
                  sums_wspts_ws_l, sums_wssas_ws_l,                            &
                  sums_wsus_ws_l, sums_wsvs_ws_l


!
!--    Set the appropriate factors for scalar and momentum advection.
       adv_sca_5 = 1.0_wp /  60.0_wp
       adv_sca_3 = 1.0_wp /  12.0_wp
       adv_sca_1 = 1.0_wp /   2.0_wp
       adv_mom_5 = 1.0_wp / 120.0_wp
       adv_mom_3 = 1.0_wp /  24.0_wp
       adv_mom_1 = 1.0_wp /   4.0_wp
!
!--    Arrays needed for statical evaluation of fluxes.
       ALLOCATE( sums_wsus_ws_l(nzb:nzt+1,0:threads_per_task-1),            &
                 sums_wsvs_ws_l(nzb:nzt+1,0:threads_per_task-1),            &
                 sums_usvs_ws_l(nzb:nzt+1,0:threads_per_task-1),            &
                 sums_vsus_ws_l(nzb:nzt+1,0:threads_per_task-1),            &
                 sums_us2_ws_l(nzb:nzt+1,0:threads_per_task-1),             &
                 sums_vs2_ws_l(nzb:nzt+1,0:threads_per_task-1),             &
                 sums_ws2_ws_l(nzb:nzt+1,0:threads_per_task-1) )

       sums_wsus_ws_l = 0.0_wp
       sums_wsvs_ws_l = 0.0_wp
       sums_usvs_ws_l  = 0.0_wp
       sums_vsus_ws_l  = 0.0_wp
       sums_us2_ws_l  = 0.0_wp
       sums_vs2_ws_l  = 0.0_wp
       sums_ws2_ws_l  = 0.0_wp

       ALLOCATE( sums_wspts_ws_l(nzb:nzt+1,0:threads_per_task-1) )
       sums_wspts_ws_l = 0.0_wp

       ALLOCATE( sums_wssas_ws_l(nzb:nzt+1,0:threads_per_task-1) )
       sums_wssas_ws_l = 0.0_wp

    END SUBROUTINE ws_init

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of flags for WS-scheme used to degrade the order of the scheme
!> near walls.
!------------------------------------------------------------------------------!
    SUBROUTINE ws_init_flags

       USE indices,                                                            &
           ONLY:  advc_flags_1, advc_flags_2, nbgp, nxl, nxlg, nxr, nxrg,      &
                  nyn, nyng, nys, nysg, nysv, nzb, nzt

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< index variable along x
       INTEGER(iwp) ::  j     !< index variable along y
       INTEGER(iwp) ::  k     !< index variable along z

       ALLOCATE( advc_flags_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                  &
                 advc_flags_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       advc_flags_1 = 0
       advc_flags_2 = 0
!
!--    Set flags to steer the degradation of the advection scheme in advec_ws
!--    near the surface and bottom.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             scalar - x-direction
!--             WS1 (0), WS3 (1), WS5 (2)
                advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 2 )
!
!--             scalar - y-direction
!--             WS1 (3), WS3 (4), WS5 (5)
                advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 5 )
!
!--             scalar - z-direction
!--             WS1 (6), WS3 (7), WS5 (8)
                IF ( k == nzt .OR. k == nzb+1 )                                &
                THEN
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 6 )
                ELSEIF ( k == nzt-1 .OR. k == nzb+2 )                          &
                THEN
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 7 )
                ELSE
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 8 )
                ENDIF

             ENDDO
          ENDDO
       ENDDO

!
!--    Set advc_flags_1 to steer the degradation of the advection scheme in advec_ws
!--    near the surface and bottom.
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             u component - x-direction
!--             WS1 (9), WS3 (10), WS5 (11)
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 11 )
!
!--             u component - y-direction
!--             WS1 (12), WS3 (13), WS5 (14)
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 14 )
!
!--             u component - z-direction
!--             WS1 (15), WS3 (16), WS5 (17)
                IF ( k == nzt .OR. k == nzb+1 )                                &
                THEN
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 15 )
                ELSEIF ( k == nzt-1  .OR. k == nzb+2 )                         &
                THEN
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 16 )
                ELSE
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 17 )
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             v component - x-direction
!--             WS1 (18), WS3 (19), WS5 (20)
                advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 20 )
!
!--             v component - y-direction
!--             WS1 (21), WS3 (22), WS5 (23)
                advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 23 )
!
!--             v component - z-direction
!--             WS1 (24), WS3 (25), WS5 (26)
                IF ( k == nzt .OR. k == nzb+1 )                                &
                THEN
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 24 )
                ELSEIF ( k == nzt-1 .OR. k == nzb+2 )                          &
                THEN
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 25 )
                ELSE
                   advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i), 26 )
                ENDIF

             ENDDO
          ENDDO
       ENDDO

       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
!
!--             w component - x-direction
!--             WS1 (27), WS3 (28), WS5 (29)
                advc_flags_1(k,j,i) = IBSET( advc_flags_1(k,j,i),29 )
!
!--             w component - y-direction
!--             WS1 (30), WS3 (31), WS5 (32)
                advc_flags_2(k,j,i) = IBSET( advc_flags_2(k,j,i), 0 )
!
!--             w component - z-direction
!--             WS1 (33), WS3 (34), WS5 (35)
                IF ( k == nzt .OR. k == nzb+1 )                                &
                THEN
                   advc_flags_2(k,j,i) = IBSET( advc_flags_2(k,j,i), 1 )
                ELSEIF ( k == nzt-1 .OR. k == nzb+2 )                          &
                THEN
                   advc_flags_2(k,j,i) = IBSET( advc_flags_2(k,j,i), 2 )
                ELSE
                   advc_flags_2(k,j,i) = IBSET( advc_flags_2(k,j,i), 3 )
                ENDIF

             ENDDO
          ENDDO
       ENDDO

!
!--    Exchange ghost points for advection flags
       CALL exchange_horiz_int( advc_flags_1, nys, nyn, nxl, nxr, nzt, nbgp )
       CALL exchange_horiz_int( advc_flags_2, nys, nyn, nxl, nxr, nzt, nbgp )


    END SUBROUTINE ws_init_flags


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize variables used for storing statistic quantities (fluxes, variances)
!------------------------------------------------------------------------------!
    SUBROUTINE ws_statistics

       USE indices,                                                            &
           ONLY: nzb, nzt

       USE kinds

       USE pegrid

       USE statistics,                                                         &
           ONLY:  sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l,                 &
                  sums_usvs_ws_l, sums_vsus_ws_l,                              &
                  sums_wspts_ws_l, sums_wssas_ws_l,                            &
                  sums_wsus_ws_l, sums_wsvs_ws_l


       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  k      !<
!
!--    The arrays needed for statistical evaluation are set to to 0 at the
!--    beginning of prognostic_equations.
       !$acc parallel present( sums_wsus_ws_l, sums_wsvs_ws_l ) &
       !$acc present( sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l ) &
       !$acc present( sums_usvs_ws_l, sums_vsus_ws_l ) &
       !$acc present( sums_wspts_ws_l, sums_wssas_ws_l )
       !$acc loop collapse(2)
       DO k = nzb, nzt+1
          DO i = 0,threads_per_task-1
             sums_wsus_ws_l(k,i) = 0.0_wp
             sums_wsvs_ws_l(k,i) = 0.0_wp
             sums_usvs_ws_l(k,i)  = 0.0_wp
             sums_vsus_ws_l(k,i)  = 0.0_wp
             sums_us2_ws_l(k,i)  = 0.0_wp
             sums_vs2_ws_l(k,i)  = 0.0_wp
             sums_ws2_ws_l(k,i)  = 0.0_wp
             sums_wspts_ws_l(k,i) = 0.0_wp
             sums_wssas_ws_l(k,i) = 0.0_wp
          ENDDO
       ENDDO
       !$acc end parallel

    END SUBROUTINE ws_statistics


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Scalar advection - Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_s_ws( sk, sk_char )

       USE arrays_3d,                                                         &
           ONLY:  ddzw, tend, u, v, w

       USE constants,                                                         &
           ONLY:  adv_sca_1, adv_sca_3, adv_sca_5

       USE control_parameters,                                                &
           ONLY:  intermediate_timestep_count

       USE grid_variables,                                                    &
           ONLY:  ddx, ddy

       USE indices,                                                           &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb,            &
                  nzt, advc_flags_1

       USE kinds

       USE statistics,                                                        &
           ONLY:  hom, sums_wspts_ws_l, sums_wssas_ws_l, weight_substep

       IMPLICIT NONE

       CHARACTER (LEN = *), INTENT(IN)    ::  sk_char !<

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit0  !<
       INTEGER(iwp) ::  ibit1  !<
       INTEGER(iwp) ::  ibit2  !<
       INTEGER(iwp) ::  ibit3  !<
       INTEGER(iwp) ::  ibit4  !<
       INTEGER(iwp) ::  ibit5  !<
       INTEGER(iwp) ::  ibit6  !<
       INTEGER(iwp) ::  ibit7  !<
       INTEGER(iwp) ::  ibit8  !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_mmm  !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  sk !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  sk !<
#endif

       REAL(wp) ::  div    !<
       REAL(wp) ::  u_comp !<
       REAL(wp) ::  v_comp !<
       REAL(wp) ::  flux_n, flux_s, flux_r, flux_l, flux_t, flux_d !<
       REAL(wp) ::  diss_n, diss_s, diss_r, diss_l, diss_t, diss_d !<

       REAL(wp) ::  tmp    !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  diss_t_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  flux_t_arr !<

       !$acc data present( tend ) &
       !$acc present( u, v, w ) &
       !$acc present( sk ) &
       !$acc create( flux_t_arr, diss_t_arr )

       !$acc parallel present ( advc_flags_1 ) &
       !$acc present( ddzw )
       !$acc loop
       DO  i = nxl, nxr
          !$acc loop
          DO j = nys, nyn
             !$acc loop
             DO  k = nzb+1, nzt

                ! left
                ibit2 = IBITS(advc_flags_1(k,j,i-1),2,1)
                ibit1 = IBITS(advc_flags_1(k,j,i-1),1,1)
                ibit0 = IBITS(advc_flags_1(k,j,i-1),0,1)

                u_comp = u(k,j,i)
                flux_l = u_comp * (                                           &
                                             ( 37.0_wp * ibit2 * adv_sca_5    &
                                          +     7.0_wp * ibit1 * adv_sca_3    &
                                          +              ibit0 * adv_sca_1    &
                                             ) *                              &
                                          ( sk(k,j,i)   + sk(k,j,i-1)    )    &
                                      -      (  8.0_wp * ibit2 * adv_sca_5    &
                                          +              ibit1 * adv_sca_3    &
                                             ) *                              &
                                          ( sk(k,j,i+1) + sk(k,j,i-2)    )    &
                                      +      (           ibit2 * adv_sca_5    &
                                             ) *                              &
                                          ( sk(k,j,i+2) + sk(k,j,i-3)    )    &
                                  )

                 diss_l = -ABS( u_comp ) * (                                  &
                                             ( 10.0_wp * ibit2 * adv_sca_5    &
                                          +     3.0_wp * ibit1 * adv_sca_3    &
                                          +              ibit0 * adv_sca_1    &
                                             ) *                              &
                                          ( sk(k,j,i)   - sk(k,j,i-1) )       &
                                      -      (  5.0_wp * ibit2 * adv_sca_5    &
                                          +              ibit1 * adv_sca_3    &
                                             ) *                              &
                                         ( sk(k,j,i+1) - sk(k,j,i-2)  )       &
                                      +      (           ibit2 * adv_sca_5    &
                                             ) *                              &
                                          ( sk(k,j,i+2) - sk(k,j,i-3) )       &
                                           )

                ! right
                ibit2 = IBITS(advc_flags_1(k,j,i),2,1)
                ibit1 = IBITS(advc_flags_1(k,j,i),1,1)
                ibit0 = IBITS(advc_flags_1(k,j,i),0,1)

                u_comp    = u(k,j,i+1)
                flux_r = u_comp * (                                           &
                          ( 37.0_wp * ibit2 * adv_sca_5                       &
                      +      7.0_wp * ibit1 * adv_sca_3                       &
                      +               ibit0 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j,i+1) + sk(k,j,i)   )                    &
                   -      (  8.0_wp * ibit2 * adv_sca_5                       &
                       +              ibit1 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j,i+2) + sk(k,j,i-1) )                    &
                   +      (           ibit2 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j,i+3) + sk(k,j,i-2) )                    &
                                  )

                diss_r = -ABS( u_comp ) * (                                   &
                          ( 10.0_wp * ibit2 * adv_sca_5                       &
                       +     3.0_wp * ibit1 * adv_sca_3                       &
                       +              ibit0 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j,i+1) - sk(k,j,i)   )                    &
                   -      (  5.0_wp * ibit2 * adv_sca_5                       &
                       +              ibit1 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j,i+2) - sk(k,j,i-1) )                    &
                   +      (           ibit2 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j,i+3) - sk(k,j,i-2) )                    &
                                          )

                ! south
                ibit5 = IBITS(advc_flags_1(k,j-1,i),5,1)
                ibit4 = IBITS(advc_flags_1(k,j-1,i),4,1)
                ibit3 = IBITS(advc_flags_1(k,j-1,i),3,1)

                v_comp = v(k,j,i)
                flux_s = v_comp * (                                           &
                                             ( 37.0_wp * ibit5 * adv_sca_5    &
                                          +     7.0_wp * ibit4 * adv_sca_3    &
                                          +              ibit3 * adv_sca_1    &
                                             ) *                              &
                                         ( sk(k,j,i)  + sk(k,j-1,i)     )     &
                                       -     (  8.0_wp * ibit5 * adv_sca_5    &
                                          +              ibit4 * adv_sca_3    &
                                              ) *                             &
                                         ( sk(k,j+1,i) + sk(k,j-2,i)    )     &
                                      +      (           ibit5 * adv_sca_5    &
                                             ) *                              &
                                        ( sk(k,j+2,i) + sk(k,j-3,i)     )     &
                                  )

                diss_s = -ABS( v_comp ) * (                                   &
                                             ( 10.0_wp * ibit5 * adv_sca_5    &
                                          +     3.0_wp * ibit4 * adv_sca_3    &
                                          +              ibit3 * adv_sca_1    &
                                             ) *                              &
                                          ( sk(k,j,i)   - sk(k,j-1,i)    )    &
                                      -      (  5.0_wp * ibit5 * adv_sca_5    &
                                          +              ibit4 * adv_sca_3    &
                                             ) *                              &
                                          ( sk(k,j+1,i) - sk(k,j-2,i)    )    &
                                      +      (           ibit5 * adv_sca_5    &
                                             ) *                              &
                                          ( sk(k,j+2,i) - sk(k,j-3,i)    )    &
                                          )

                ! north
                ibit5 = IBITS(advc_flags_1(k,j,i),5,1)
                ibit4 = IBITS(advc_flags_1(k,j,i),4,1)
                ibit3 = IBITS(advc_flags_1(k,j,i),3,1)

                v_comp    = v(k,j+1,i)
                flux_n = v_comp * (                                           &
                          ( 37.0_wp * ibit5 * adv_sca_5                       &
                       +     7.0_wp * ibit4 * adv_sca_3                       &
                       +              ibit3 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j+1,i) + sk(k,j,i)   )                    &
                   -      (  8.0_wp * ibit5 * adv_sca_5                       &
                       +              ibit4 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j+2,i) + sk(k,j-1,i) )                    &
                   +      (           ibit5 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j+3,i) + sk(k,j-2,i) )                    &
                                  )

                diss_n = -ABS( v_comp ) * (                                   &
                          ( 10.0_wp * ibit5 * adv_sca_5                       &
                       +     3.0_wp * ibit4 * adv_sca_3                       &
                       +              ibit3 * adv_sca_1                       &
                          ) *                                                 &
                             ( sk(k,j+1,i) - sk(k,j,i)    )                   &
                   -      (  5.0_wp * ibit5 * adv_sca_5                       &
                       +              ibit4 * adv_sca_3                       &
                          ) *                                                 &
                             ( sk(k,j+2,i) - sk(k,j-1,i)  )                   &
                   +      (           ibit5 * adv_sca_5                       &
                          ) *                                                 &
                             ( sk(k,j+3,i) - sk(k,j-2,i) )                    &
                                          )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ! bottom
                ibit8 = IBITS(advc_flags_1(k-1,j,i),8,1)
                ibit7 = IBITS(advc_flags_1(k-1,j,i),7,1)
                ibit6 = IBITS(advc_flags_1(k-1,j,i),6,1)

                k_pp  = k + 2 * ibit8
                k_mm  = k - 2 * ( ibit7 + ibit8 )
                k_mmm = k - 3 * ibit8

                flux_d = w(k-1,j,i) * (                                       &
                           ( 37.0_wp * ibit8 * adv_sca_5                      &
                        +     7.0_wp * ibit7 * adv_sca_3                      &
                        +           ibit6 * adv_sca_1                         &
                           ) *                                                &
                                   ( sk(k,j,i)  + sk(k-1,j,i)    )            &
                    -      (  8.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k+1,j,i) + sk(k_mm,j,i)  )            &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *     ( sk(k_pp,j,i)+ sk(k_mmm,j,i) )            &
                                      )

                diss_d = -ABS( w(k-1,j,i) )  * (                              &
                           ( 10.0_wp * ibit8 * adv_sca_5                      &
                        +     3.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k,j,i)   - sk(k-1,j,i)    )           &
                    -      (  5.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k+1,j,i)  - sk(k_mm,j,i)  )           &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i) - sk(k_mmm,j,i) )           &
                                               )

                ! top
                ibit8 = IBITS(advc_flags_1(k,j,i),8,1)
                ibit7 = IBITS(advc_flags_1(k,j,i),7,1)
                ibit6 = IBITS(advc_flags_1(k,j,i),6,1)

                k_ppp = k + 3 * ibit8
                k_pp  = k + 2 * ( 1 - ibit6  )
                k_mm  = k - 2 * ibit8


                flux_t = w(k,j,i) * (                                         &
                           ( 37.0_wp * ibit8 * adv_sca_5                      &
                        +     7.0_wp * ibit7 * adv_sca_3                      &
                        +           ibit6 * adv_sca_1                         &
                           ) *                                                &
                                   ( sk(k+1,j,i)  + sk(k,j,i)    )            &
                    -      (  8.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i) + sk(k-1,j,i)  )            &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *     ( sk(k_ppp,j,i)+ sk(k_mm,j,i) )            &
                                    )

                diss_t = -ABS( w(k,j,i) ) * (                                 &
                           ( 10.0_wp * ibit8 * adv_sca_5                      &
                        +     3.0_wp * ibit7 * adv_sca_3                      &
                        +              ibit6 * adv_sca_1                      &
                           ) *                                                &
                                   ( sk(k+1,j,i)   - sk(k,j,i)    )           &
                    -      (  5.0_wp * ibit8 * adv_sca_5                      &
                        +              ibit7 * adv_sca_3                      &
                           ) *                                                &
                                   ( sk(k_pp,j,i)  - sk(k-1,j,i)  )           &
                    +      (           ibit8 * adv_sca_5                      &
                           ) *                                                &
                                   ( sk(k_ppp,j,i) - sk(k_mm,j,i) )           &
                                            )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div   =   ( u(k,j,i+1) * ( ibit0 + ibit1 + ibit2 )            &
                          - u(k,j,i)   * ( IBITS(advc_flags_1(k,j,i-1),0,1)   &
                                         + IBITS(advc_flags_1(k,j,i-1),1,1)   &
                                         + IBITS(advc_flags_1(k,j,i-1),2,1)   &
                                         )                                    &
                          ) * ddx                                             &
                        + ( v(k,j+1,i) * ( ibit3 + ibit4 + ibit5 )            &
                          - v(k,j,i)   * ( IBITS(advc_flags_1(k,j-1,i),3,1)   &
                                         + IBITS(advc_flags_1(k,j-1,i),4,1)   &
                                         + IBITS(advc_flags_1(k,j-1,i),5,1)   &
                                         )                                    &
                          ) * ddy                                             &
                        + ( w(k,j,i)   * ( ibit6 + ibit7 + ibit8 )            &
                          - w(k-1,j,i) * ( IBITS(advc_flags_1(k-1,j,i),6,1)   &
                                         + IBITS(advc_flags_1(k-1,j,i),7,1)   &
                                         + IBITS(advc_flags_1(k-1,j,i),8,1)   &
                                         )                                    &
                          ) * ddzw(k)


                tend(k,j,i) = tend(k,j,i) - (                                 &
                        ( flux_r + diss_r - flux_l - diss_l ) * ddx           &
                      + ( flux_n + diss_n - flux_s - diss_s ) * ddy           &
                      + ( flux_t + diss_t - flux_d - diss_d ) * ddzw(k)       &
                                            ) + sk(k,j,i) * div
                flux_t_arr(k,j,i) = flux_t
                diss_t_arr(k,j,i) = diss_t

             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

!--   Evaluation of statistics.

       SELECT CASE ( sk_char )
          CASE ( 'pt' )
             !$acc parallel present( sums_wspts_ws_l ) &
             !$acc present( hom, weight_substep )
             !$acc loop
             DO  k = nzb+1, nzt
                tmp = 0.0_wp
                !$acc loop collapse(2) reduction(+:tmp)
                DO  i = nxl, nxr
                   DO j = nys, nyn
                      tmp = tmp                                                &
                          + ( flux_t_arr(k,j,i)                                &
                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3)                   )  &
                            + diss_t_arr(k,j,i)                                &
                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                               *   ABS(w(k,j,i) - hom(k,1,3)               )   &
                            ) * weight_substep(intermediate_timestep_count)
                   ENDDO
                ENDDO
                sums_wspts_ws_l(k,tn) = sums_wspts_ws_l(k,tn) + tmp
             ENDDO
             !$acc end parallel
          CASE ( 'sa' )
             !$acc parallel present( sums_wssas_ws_l ) &
             !$acc present( hom, weight_substep )
             !$acc loop
             DO  k = nzb+1, nzt
                tmp = 0.0_wp
                !$acc loop collapse(2) reduction(+:tmp)
                DO  i = nxl, nxr
                   DO j = nys, nyn
                      tmp = tmp                                                &
                          + ( flux_t_arr(k,j,i)                                &
                                / ( w(k,j,i) + SIGN( 1.0E-20_wp, w(k,j,i) ) )  &
                                * ( w(k,j,i) - hom(k,1,3)                   )  &
                            + diss_t_arr(k,j,i)                                &
                                / ( ABS(w(k,j,i)) + 1.0E-20_wp              )  &
                                *   ABS(w(k,j,i) - hom(k,1,3)               )  &
                            ) * weight_substep(intermediate_timestep_count)
                   ENDDO
                ENDDO
                sums_wssas_ws_l(k,tn) = sums_wssas_ws_l(k,tn) + tmp
             ENDDO
             !$acc end parallel
       END SELECT
       !$acc end data

    END SUBROUTINE advec_s_ws


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of u - Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_u_ws

       USE arrays_3d,                                                          &
           ONLY:  ddzw, tend, u, v, w

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nzb, nzt, advc_flags_1

       USE kinds

       USE statistics,                                                         &
           ONLY:  hom, sums_us2_ws_l, sums_vsus_ws_l, sums_wsus_ws_l,          &
                  weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit9  !<
       INTEGER(iwp) ::  ibit10 !<
       INTEGER(iwp) ::  ibit11 !<
       INTEGER(iwp) ::  ibit12 !<
       INTEGER(iwp) ::  ibit13 !<
       INTEGER(iwp) ::  ibit14 !<
       INTEGER(iwp) ::  ibit15 !<
       INTEGER(iwp) ::  ibit16 !<
       INTEGER(iwp) ::  ibit17 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_mmm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<

       REAL(wp) ::  div    !<
       REAL(wp) ::  u_comp !<
       REAL(wp) ::  v_comp !<
       REAL(wp) ::  w_comp !<
       REAL(wp) ::  flux_n, flux_s, flux_r, flux_l, flux_t, flux_d !<
       REAL(wp) ::  diss_n, diss_s, diss_r, diss_l, diss_t, diss_d !<

       REAL(wp) ::  tmp    !<
       REAL(wp) ::  tmp2   !<
       REAL(wp) ::  tmp3   !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  diss_r_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  flux_r_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  diss_n_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  flux_n_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  diss_t_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  flux_t_arr !<


!--    Computation of interior fluxes and tendency terms

       !$acc data present( tend ) &
       !$acc present( u, v, w ) &
       !$acc create( flux_r_arr, diss_r_arr ) &
       !$acc create( flux_n_arr, diss_n_arr ) &
       !$acc create( flux_t_arr, diss_t_arr )

       !$acc parallel present ( advc_flags_1 ) &
       !$acc present( tend ) &
       !$acc present( ddzw )
       !$acc loop
       DO i = nxlu, nxr
          !$acc loop
          DO  j = nys, nyn
             !$acc loop
             DO  k = nzb+1, nzt

                ! left
                ibit11 = IBITS(advc_flags_1(k,j,i-1),11,1)
                ibit10 = IBITS(advc_flags_1(k,j,i-1),10,1)
                ibit9  = IBITS(advc_flags_1(k,j,i-1),9,1)

                u_comp = u(k,j,i) + u(k,j,i-1)
                flux_l = u_comp * (                                            &
                          ( 37.0_wp * ibit11 * adv_mom_5                       &
                       +     7.0_wp * ibit10 * adv_mom_3                       &
                       +              ibit9  * adv_mom_1                       &
                          ) *                                                  &
                        ( u(k,j,i)   + u(k,j,i-1) )                            &
                   -      (  8.0_wp * ibit11 * adv_mom_5                       &
                       +              ibit10 * adv_mom_3                       &
                          ) *                                                  &
                        ( u(k,j,i+1) + u(k,j,i-2) )                            &
                   +      (           ibit11 * adv_mom_5                       &
                          ) *                                                  &
                        ( u(k,j,i+2) + u(k,j,i-3) )                            &
                                  )

                diss_l = - ABS( u_comp ) * (                                   &
                          ( 10.0_wp * ibit11 * adv_mom_5                       &
                       +     3.0_wp * ibit10 * adv_mom_3                       &
                       +              ibit9  * adv_mom_1                       &
                          ) *                                                  &
                        ( u(k,j,i)   - u(k,j,i-1) )                            &
                   -      (  5.0_wp * ibit11 * adv_mom_5                       &
                       +              ibit10 * adv_mom_3                       &
                          ) *                                                  &
                        ( u(k,j,i+1) - u(k,j,i-2) )                            &
                   +      (           ibit11 * adv_mom_5                       &
                          ) *                                                  &
                        ( u(k,j,i+2) - u(k,j,i-3) )                            &
                                           )

                ! right
                ibit11 = IBITS(advc_flags_1(k,j,i),11,1)
                ibit10 = IBITS(advc_flags_1(k,j,i),10,1)
                ibit9  = IBITS(advc_flags_1(k,j,i),9,1)

                u_comp = u(k,j,i+1) + u(k,j,i)
                flux_r = u_comp * (                                            &
                          ( 37.0_wp * ibit11 * adv_mom_5                       &
                       +     7.0_wp * ibit10 * adv_mom_3                       &
                       +              ibit9  * adv_mom_1                       &
                          ) *                                                  &
                                 ( u(k,j,i+1) + u(k,j,i)   )                   &
                   -      (  8.0_wp * ibit11 * adv_mom_5                       &
                       +              ibit10 * adv_mom_3                       &
                          ) *                                                  &
                                 ( u(k,j,i+2) + u(k,j,i-1) )                   &
                   +      (           ibit11 * adv_mom_5                       &
                          ) *                                                  &
                                 ( u(k,j,i+3) + u(k,j,i-2) )                   &
                                 )

                diss_r = - ABS( u_comp ) * (                                   &
                          ( 10.0_wp * ibit11 * adv_mom_5                       &
                       +     3.0_wp * ibit10 * adv_mom_3                       &
                       +              ibit9  * adv_mom_1                       &
                          ) *                                                  &
                                 ( u(k,j,i+1) - u(k,j,i)  )                    &
                   -      (  5.0_wp * ibit11 * adv_mom_5                       &
                       +              ibit10 * adv_mom_3                       &
                          ) *                                                  &
                                 ( u(k,j,i+2) - u(k,j,i-1) )                   &
                   +      (           ibit11 * adv_mom_5                       &
                          ) *                                                  &
                                 ( u(k,j,i+3) - u(k,j,i-2) )                   &
                                           )
                ! south
                ibit14 = IBITS(advc_flags_1(k,j-1,i),14,1)
                ibit13 = IBITS(advc_flags_1(k,j-1,i),13,1)
                ibit12 = IBITS(advc_flags_1(k,j-1,i),12,1)

                v_comp = v(k,j,i) + v(k,j,i-1)
                flux_s = v_comp * (                                            &
                              ( 37.0_wp * ibit14 * adv_mom_5                   &
                           +     7.0_wp * ibit13 * adv_mom_3                   &
                           +              ibit12 * adv_mom_1                   &
                              ) *                                              &
                                ( u(k,j,i)   + u(k,j-1,i) )                    &
                       -      (  8.0_wp * ibit14 * adv_mom_5                   &
                       +                  ibit13 * adv_mom_3                   &
                              ) *                                              &
                                ( u(k,j+1,i) + u(k,j-2,i) )                    &
                   +      (               ibit14 * adv_mom_5                   &
                          ) *                                                  &
                                ( u(k,j+2,i) + u(k,j-3,i) )                    &
                                  )

                diss_s = - ABS ( v_comp ) * (                                  &
                          ( 10.0_wp * ibit14 * adv_mom_5                       &
                       +     3.0_wp * ibit13 * adv_mom_3                       &
                       +              ibit12 * adv_mom_1                       &
                          ) *                                                  &
                            ( u(k,j,i)   - u(k,j-1,i) )                        &
                   -      (  5.0_wp * ibit14 * adv_mom_5                       &
                       +              ibit13 * adv_mom_3                       &
                          ) *                                                  &
                            ( u(k,j+1,i) - u(k,j-2,i) )                        &
                   +      (           ibit14 * adv_mom_5                       &
                          ) *                                                  &
                            ( u(k,j+2,i) - u(k,j-3,i) )                        &
                                            )

                ! north
                ibit14 = IBITS(advc_flags_1(k,j,i),14,1)
                ibit13 = IBITS(advc_flags_1(k,j,i),13,1)
                ibit12 = IBITS(advc_flags_1(k,j,i),12,1)

                v_comp    = v(k,j+1,i) + v(k,j+1,i-1)
                flux_n = v_comp * (                                            &
                          ( 37.0_wp * ibit14 * adv_mom_5                       &
                       +     7.0_wp * ibit13 * adv_mom_3                       &
                       +              ibit12 * adv_mom_1                       &
                          ) *                                                  &
                                 ( u(k,j+1,i) + u(k,j,i)   )                   &
                   -      (  8.0_wp * ibit14 * adv_mom_5                       &
                       +              ibit13 * adv_mom_3                       &
                          ) *                                                  &
                                 ( u(k,j+2,i) + u(k,j-1,i) )                   &
                   +      (           ibit14 * adv_mom_5                       &
                          ) *                                                  &
                                 ( u(k,j+3,i) + u(k,j-2,i) )                   &
                                  )

                diss_n = - ABS ( v_comp ) * (                                  &
                          ( 10.0_wp * ibit14 * adv_mom_5                       &
                       +     3.0_wp * ibit13 * adv_mom_3                       &
                       +              ibit12 * adv_mom_1                       &
                          ) *                                                  &
                                 ( u(k,j+1,i) - u(k,j,i)  )                    &
                   -      (  5.0_wp * ibit14 * adv_mom_5                       &
                       +              ibit13 * adv_mom_3                       &
                          ) *                                                  &
                                 ( u(k,j+2,i) - u(k,j-1,i) )                   &
                   +      (           ibit14 * adv_mom_5                       &
                          ) *                                                  &
                                 ( u(k,j+3,i) - u(k,j-2,i) )                   &
                                            )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ! bottom
                ibit17 = IBITS(advc_flags_1(k-1,j,i),17,1)
                ibit16 = IBITS(advc_flags_1(k-1,j,i),16,1)
                ibit15 = IBITS(advc_flags_1(k-1,j,i),15,1)

                k_pp  = k + 2 * ibit17
                k_mm  = k - 2 * ( ibit16 + ibit17 )
                k_mmm = k - 3 * ibit17

                w_comp    = w(k-1,j,i) + w(k-1,j,i-1)
                flux_d = w_comp * (                                            &
                          ( 37.0_wp * ibit17 * adv_mom_5                       &
                       +     7.0_wp * ibit16 * adv_mom_3                       &
                       +              ibit15 * adv_mom_1                       &
                          ) *                                                  &
                             ( u(k,j,i)  + u(k-1,j,i)     )                    &
                   -      (  8.0_wp * ibit17 * adv_mom_5                       &
                       +              ibit16 * adv_mom_3                       &
                          ) *                                                  &
                             ( u(k+1,j,i) + u(k_mm,j,i)   )                    &
                   +      (           ibit17 * adv_mom_5                       &
                          ) *                                                  &
                             ( u(k_pp,j,i) + u(k_mmm,j,i) )                    &
                                  )

                diss_d = - ABS( w_comp ) * (                                   &
                          ( 10.0_wp * ibit17 * adv_mom_5                       &
                       +     3.0_wp * ibit16 * adv_mom_3                       &
                       +              ibit15 * adv_mom_1                       &
                          ) *                                                  &
                             ( u(k,j,i)   - u(k-1,j,i)    )                    &
                   -      (  5.0_wp * ibit17 * adv_mom_5                       &
                       +              ibit16 * adv_mom_3                       &
                          ) *                                                  &
                             ( u(k+1,j,i)  - u(k_mm,j,i)  )                    &
                   +      (           ibit17 * adv_mom_5                       &
                           ) *                                                 &
                             ( u(k_pp,j,i) - u(k_mmm,j,i) )                    &
                                           )

                ! top
                ibit17 = IBITS(advc_flags_1(k,j,i),17,1)
                ibit16 = IBITS(advc_flags_1(k,j,i),16,1)
                ibit15 = IBITS(advc_flags_1(k,j,i),15,1)

                k_ppp = k + 3 * ibit17
                k_pp  = k + 2 * ( 1 - ibit15  )
                k_mm  = k - 2 * ibit17

                w_comp    = w(k,j,i) + w(k,j,i-1)
                flux_t = w_comp * (                                            &
                          ( 37.0_wp * ibit17 * adv_mom_5                       &
                       +     7.0_wp * ibit16 * adv_mom_3                       &
                       +              ibit15 * adv_mom_1                       &
                          ) *                                                  &
                             ( u(k+1,j,i)  + u(k,j,i)     )                    &
                   -      (  8.0_wp * ibit17 * adv_mom_5                       &
                       +              ibit16 * adv_mom_3                       &
                          ) *                                                  &
                             ( u(k_pp,j,i) + u(k-1,j,i)   )                    &
                   +      (           ibit17 * adv_mom_5                       &
                          ) *                                                  &
                             ( u(k_ppp,j,i) + u(k_mm,j,i) )                    &
                                  )

                diss_t = - ABS( w_comp ) * (                                   &
                          ( 10.0_wp * ibit17 * adv_mom_5                       &
                       +     3.0_wp * ibit16 * adv_mom_3                       &
                       +              ibit15 * adv_mom_1                       &
                          ) *                                                  &
                             ( u(k+1,j,i)   - u(k,j,i)    )                    &
                   -      (  5.0_wp * ibit17 * adv_mom_5                       &
                       +              ibit16 * adv_mom_3                       &
                          ) *                                                  &
                             ( u(k_pp,j,i)  - u(k-1,j,i)  )                    &
                   +      (           ibit17 * adv_mom_5                       &
                           ) *                                                 &
                             ( u(k_ppp,j,i) - u(k_mm,j,i) )                    &
                                           )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( u_comp    * ( ibit9 + ibit10 + ibit11 )              &
                - ( u(k,j,i)   + u(k,j,i-1)   )                                &
                                    * ( IBITS(advc_flags_1(k,j,i-1),9,1)       &
                                      + IBITS(advc_flags_1(k,j,i-1),10,1)      &
                                      + IBITS(advc_flags_1(k,j,i-1),11,1)      &
                                      )                                        &
                  ) * ddx                                                      &
               +  ( v_comp          * ( ibit12 + ibit13 + ibit14 )             &
                  - ( v(k,j,i)   + v(k,j,i-1 )  )                              &
                                    * ( IBITS(advc_flags_1(k,j-1,i),12,1)      &
                                      + IBITS(advc_flags_1(k,j-1,i),13,1)      &
                                      + IBITS(advc_flags_1(k,j-1,i),14,1)      &
                                      )                                        &
                  ) * ddy                                                      &
               +  ( w_comp          * ( ibit15 + ibit16 + ibit17 )             &
                - ( w(k-1,j,i) + w(k-1,j,i-1) )                                &
                                    * ( IBITS(advc_flags_1(k-1,j,i),15,1)      &
                                      + IBITS(advc_flags_1(k-1,j,i),16,1)      &
                                      + IBITS(advc_flags_1(k-1,j,i),17,1)      &
                                      )                                        &
                  ) * ddzw(k)                                                  &
                ) * 0.5_wp



                tend(k,j,i) = tend(k,j,i) - (                                  &
                       ( flux_r + diss_r - flux_l - diss_l ) * ddx             &
                     + ( flux_n + diss_n - flux_s - diss_s ) * ddy             &
                     + ( flux_t + diss_t - flux_d - diss_d ) * ddzw(k)         &
                                            ) + div * u(k,j,i)

                flux_r_arr(k,j,i) = flux_r
                diss_r_arr(k,j,i) = diss_r
                flux_n_arr(k,j,i) = flux_n
                diss_n_arr(k,j,i) = diss_n
                flux_t_arr(k,j,i) = flux_t
                diss_t_arr(k,j,i) = diss_t

             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

       !$acc parallel present( sums_us2_ws_l, sums_vsus_ws_l, sums_wsus_ws_l ) &
       !$acc present( hom, weight_substep )
       !$acc loop
       DO  k = nzb+1, nzt
          tmp = 0.0_wp
          tmp2 = 0.0_wp
          tmp3 = 0.0_wp
          !$acc loop collapse(2) reduction(+:tmp,tmp2,tmp3)
          DO  i = nxl, nxr
             DO j = nys, nyn
                u_comp = u(k,j,i) + u(k,j,i+1)
                v_comp = v(k,j+1,i) + v(k,j+1,i-1)
                w_comp = w(k,j,i) + w(k,j,i-1)
!
!--             Statistical Evaluation of u'u'.
                tmp = tmp                                                      &
                    + ( flux_r_arr(k,j,i)                                      &
                      * ( u_comp - 2.0_wp * hom(k,1,1)             )           &
                      / ( u_comp + SIGN( 1.0E-20_wp, u_comp ) )                &
                    + diss_r_arr(k,j,i)                                        &
                      *   ABS( u_comp - 2.0_wp * hom(k,1,1)        )           &
                      / ( ABS( u_comp ) + 1.0E-20_wp               )           &
                      ) *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of v'u'.
                tmp2 = tmp2                                                    &
                    + ( flux_n_arr(k,j,i)                                      &
                      * ( v_comp - 2.0_wp * hom(k,1,2)             )           &
                      / ( v_comp + SIGN( 1.0E-20_wp, v_comp )      )           &
                    + diss_n_arr(k,j,i)                                        &
                      *   ABS( v_comp - 2.0_wp * hom(k,1,2)        )           &
                      / ( ABS( v_comp ) + 1.0E-20_wp               )           &
                      ) *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of w'u'.
                tmp3 = tmp3                                                    &
                    + ( flux_t_arr(k,j,i)                                      &
                      * ( w_comp - 2.0_wp * hom(k,1,3)             )           &
                      / ( w_comp + SIGN( 1.0E-20_wp, w_comp )      )           &
                    + diss_t_arr(k,j,i)                                        &
                      *   ABS( w_comp - 2.0_wp * hom(k,1,3)        )           &
                      / ( ABS( w_comp ) + 1.0E-20_wp               )           &
                      ) *   weight_substep(intermediate_timestep_count)
             ENDDO
          ENDDO
          sums_us2_ws_l(k,tn)  = sums_us2_ws_l(k,tn)  + tmp
          sums_vsus_ws_l(k,tn) = sums_vsus_ws_l(k,tn) + tmp2
          sums_wsus_ws_l(k,tn) = sums_wsus_ws_l(k,tn) + tmp3
       ENDDO
       sums_us2_ws_l(nzb,tn) = sums_us2_ws_l(nzb+1,tn)
       !$acc end parallel
       !$acc end data


    END SUBROUTINE advec_u_ws


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of v - Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_v_ws

       USE arrays_3d,                                                          &
           ONLY:  ddzw, tend, u, v, w

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nysv, nzb, nzt, advc_flags_1

       USE kinds

       USE statistics,                                                         &
           ONLY:  hom, sums_vs2_ws_l, sums_usvs_ws_l, sums_wsvs_ws_l,          &
                  weight_substep

       IMPLICIT NONE


       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit18 !<
       INTEGER(iwp) ::  ibit19 !<
       INTEGER(iwp) ::  ibit20 !<
       INTEGER(iwp) ::  ibit21 !<
       INTEGER(iwp) ::  ibit22 !<
       INTEGER(iwp) ::  ibit23 !<
       INTEGER(iwp) ::  ibit24 !<
       INTEGER(iwp) ::  ibit25 !<
       INTEGER(iwp) ::  ibit26 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_mmm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<

       REAL(wp) ::  div    !<
       REAL(wp) ::  u_comp !<
       REAL(wp) ::  v_comp !<
       REAL(wp) ::  w_comp !<
       REAL(wp) ::  flux_n, flux_s, flux_r, flux_l, flux_t, flux_d !<
       REAL(wp) ::  diss_n, diss_s, diss_r, diss_l, diss_t, diss_d !<

       REAL(wp) ::  tmp    !<
       REAL(wp) ::  tmp2   !<
       REAL(wp) ::  tmp3   !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  diss_n_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  flux_n_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  diss_r_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  flux_r_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  diss_t_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  flux_t_arr !<

!
!--    Computation of interior fluxes and tendency terms

       !$acc data present( tend ) &
       !$acc present( u, v, w ) &
       !$acc create( flux_r_arr, diss_r_arr ) &
       !$acc create( flux_n_arr, diss_n_arr ) &
       !$acc create( flux_t_arr, diss_t_arr )

       !$acc parallel present ( advc_flags_1 ) &
       !$acc present( ddzw )
       !$acc loop
       DO i = nxl, nxr
          !$acc loop
          DO  j = nysv, nyn
             !$acc loop
             DO  k = nzb+1, nzt

                ! left
                ibit20 = IBITS(advc_flags_1(k,j,i-1),20,1)
                ibit19 = IBITS(advc_flags_1(k,j,i-1),19,1)
                ibit18 = IBITS(advc_flags_1(k,j,i-1),18,1)

                u_comp = u(k,j-1,i) + u(k,j,i)
                flux_l = u_comp * (                                            &
                          ( 37.0_wp * ibit20 * adv_mom_5                       &
                       +     7.0_wp * ibit19 * adv_mom_3                       &
                       +              ibit18 * adv_mom_1                       &
                          ) *                                                  &
                         ( v(k,j,i)   + v(k,j,i-1) )                           &
                   -      (  8.0_wp * ibit20 * adv_mom_5                       &
                       +              ibit19 * adv_mom_3                       &
                          ) *                                                  &
                         ( v(k,j,i+1) + v(k,j,i-2) )                           &
                   +      (           ibit20 * adv_mom_5                       &
                          ) *                                                  &
                         ( v(k,j,i+2) + v(k,j,i-3) )                           &
                                  )

                diss_l = - ABS( u_comp ) * (                                   &
                          ( 10.0_wp * ibit20 * adv_mom_5                       &
                       +     3.0_wp * ibit19 * adv_mom_3                       &
                       +              ibit18 * adv_mom_1                       &
                          ) *                                                  &
                         ( v(k,j,i)   - v(k,j,i-1) )                           &
                   -      (  5.0_wp * ibit20 * adv_mom_5                       &
                       +              ibit19 * adv_mom_3                       &
                          ) *                                                  &
                         ( v(k,j,i+1) - v(k,j,i-2) )                           &
                   +      (           ibit20 * adv_mom_5                       &
                          ) *                                                  &
                         ( v(k,j,i+2) - v(k,j,i-3) )                           &
                                           )

                ! right
                ibit20 = IBITS(advc_flags_1(k,j,i),20,1)
                ibit19 = IBITS(advc_flags_1(k,j,i),19,1)
                ibit18 = IBITS(advc_flags_1(k,j,i),18,1)

                u_comp    = u(k,j-1,i+1) + u(k,j,i+1)
                flux_r = u_comp * (                                            &
                          ( 37.0_wp * ibit20 * adv_mom_5                       &
                       +     7.0_wp * ibit19 * adv_mom_3                       &
                       +              ibit18 * adv_mom_1                       &
                          ) *                                                  &
                                 ( v(k,j,i+1) + v(k,j,i)   )                   &
                   -      (  8.0_wp * ibit20 * adv_mom_5                       &
                       +              ibit19 * adv_mom_3                       &
                          ) *                                                  &
                                 ( v(k,j,i+2) + v(k,j,i-1) )                   &
                   +      (           ibit20 * adv_mom_5                       &
                          ) *                                                  &
                                 ( v(k,j,i+3) + v(k,j,i-2) )                   &
                                  )

                diss_r = - ABS( u_comp ) * (                                   &
                          ( 10.0_wp * ibit20 * adv_mom_5                       &
                       +     3.0_wp * ibit19 * adv_mom_3                       &
                       +              ibit18 * adv_mom_1                       &
                          ) *                                                  &
                                 ( v(k,j,i+1) - v(k,j,i)  )                    &
                   -      (  5.0_wp * ibit20 * adv_mom_5                       &
                       +              ibit19 * adv_mom_3                       &
                          ) *                                                  &
                                 ( v(k,j,i+2) - v(k,j,i-1) )                   &
                   +      (           ibit20 * adv_mom_5                       &
                          ) *                                                  &
                                 ( v(k,j,i+3) - v(k,j,i-2) )                   &
                                           )

                ! south
                ibit23 = IBITS(advc_flags_1(k,j-1,i),23,1)
                ibit22 = IBITS(advc_flags_1(k,j-1,i),22,1)
                ibit21 = IBITS(advc_flags_1(k,j-1,i),21,1)

                v_comp = v(k,j,i) + v(k,j-1,i)
                flux_s = v_comp * (                                            &
                                   ( 37.0_wp * ibit23 * adv_mom_5              &
                                +     7.0_wp * ibit22 * adv_mom_3              &
                                +              ibit21 * adv_mom_1              &
                                   ) *                                         &
                                     ( v(k,j,i)   + v(k,j-1,i) )               &
                            -      (  8.0_wp * ibit23 * adv_mom_5              &
                                +              ibit22 * adv_mom_3              &
                                   ) *                                         &
                                     ( v(k,j+1,i) + v(k,j-2,i) )               &
                            +      (           ibit23 * adv_mom_5              &
                                   ) *                                         &
                                     ( v(k,j+2,i) + v(k,j-3,i) )               &
                                  )

                diss_s = - ABS( v_comp ) * (                  &
                                   ( 10.0_wp * ibit23 * adv_mom_5              &
                                +     3.0_wp * ibit22 * adv_mom_3              &
                                +              ibit21 * adv_mom_1              &
                                   ) *                                         &
                                     ( v(k,j,i)   - v(k,j-1,i) )               &
                            -      (  5.0_wp * ibit23 * adv_mom_5              &
                                +              ibit22 * adv_mom_3              &
                                   ) *                                         &
                                     ( v(k,j+1,i) - v(k,j-2,i) )               &
                            +      (           ibit23 * adv_mom_5              &
                                   ) *                                         &
                                     ( v(k,j+2,i) - v(k,j-3,i) )               &
                                           )

                ! north
                ibit23 = IBITS(advc_flags_1(k,j,i),23,1)
                ibit22 = IBITS(advc_flags_1(k,j,i),22,1)
                ibit21 = IBITS(advc_flags_1(k,j,i),21,1)

                v_comp = v(k,j+1,i) + v(k,j,i)
                flux_n = v_comp * (                                            &
                          ( 37.0_wp * ibit23 * adv_mom_5                       &
                       +     7.0_wp * ibit22 * adv_mom_3                       &
                       +              ibit21 * adv_mom_1                       &
                          ) *                                                  &
                                 ( v(k,j+1,i) + v(k,j,i)   )                   &
                   -      (  8.0_wp * ibit23 * adv_mom_5                       &
                       +              ibit22 * adv_mom_3                       &
                          ) *                                                  &
                                 ( v(k,j+2,i) + v(k,j-1,i) )                   &
                   +      (           ibit23 * adv_mom_5                       &
                          ) *                                                  &
                                 ( v(k,j+3,i) + v(k,j-2,i) )                   &
                                  )

                diss_n = - ABS( v_comp ) * (                                   &
                          ( 10.0_wp * ibit23 * adv_mom_5                       &
                       +     3.0_wp * ibit22 * adv_mom_3                       &
                       +              ibit21 * adv_mom_1                       &
                          ) *                                                  &
                                 ( v(k,j+1,i) - v(k,j,i)  )                    &
                   -      (  5.0_wp * ibit23 * adv_mom_5                       &
                       +              ibit22 * adv_mom_3                       &
                          ) *                                                  &
                                 ( v(k,j+2,i) - v(k,j-1,i) )                   &
                   +      (           ibit23 * adv_mom_5                       &
                          ) *                                                  &
                                 ( v(k,j+3,i) - v(k,j-2,i) )                   &
                                           )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ! bottom
                ibit26 = IBITS(advc_flags_1(k-1,j,i),26,1)
                ibit25 = IBITS(advc_flags_1(k-1,j,i),25,1)
                ibit24 = IBITS(advc_flags_1(k-1,j,i),24,1)

                k_pp  = k + 2 * ibit26
                k_mm  = k - 2 * ( ibit25 + ibit26 )
                k_mmm = k - 3 * ibit26

                w_comp    = w(k-1,j-1,i) + w(k-1,j,i)
                flux_d = w_comp * (                                            &
                          ( 37.0_wp * ibit26 * adv_mom_5                       &
                       +     7.0_wp * ibit25 * adv_mom_3                       &
                       +              ibit24 * adv_mom_1                       &
                          ) *                                                  &
                             ( v(k,j,i)   + v(k-1,j,i)    )                    &
                   -      (  8.0_wp * ibit26 * adv_mom_5                       &
                       +              ibit25 * adv_mom_3                       &
                          ) *                                                  &
                             ( v(k+1,j,i)  + v(k_mm,j,i)  )                    &
                   +      (           ibit26 * adv_mom_5                       &
                          ) *                                                  &
                             ( v(k_pp,j,i) + v(k_mmm,j,i) )                    &
                                  )

                diss_d = - ABS( w_comp ) * (                                   &
                          ( 10.0_wp * ibit26 * adv_mom_5                       &
                       +     3.0_wp * ibit25 * adv_mom_3                       &
                       +              ibit24 * adv_mom_1                       &
                          ) *                                                  &
                             ( v(k,j,i)   - v(k-1,j,i)    )                    &
                   -      (  5.0_wp * ibit26 * adv_mom_5                       &
                       +              ibit25 * adv_mom_3                       &
                          ) *                                                  &
                             ( v(k+1,j,i)  - v(k_mm,j,i)  )                    &
                   +      (           ibit26 * adv_mom_5                       &
                          ) *                                                  &
                             ( v(k_pp,j,i) - v(k_mmm,j,i) )                    &
                                           )
                ! top
                ibit26 = IBITS(advc_flags_1(k,j,i),26,1)
                ibit25 = IBITS(advc_flags_1(k,j,i),25,1)
                ibit24 = IBITS(advc_flags_1(k,j,i),24,1)

                k_ppp = k + 3 * ibit26
                k_pp  = k + 2 * ( 1 - ibit24  )
                k_mm  = k - 2 * ibit26

                w_comp    = w(k,j-1,i) + w(k,j,i)
                flux_t = w_comp * (                                            &
                          ( 37.0_wp * ibit26 * adv_mom_5                       &
                       +     7.0_wp * ibit25 * adv_mom_3                       &
                       +              ibit24 * adv_mom_1                       &
                          ) *                                                  &
                             ( v(k+1,j,i)   + v(k,j,i)    )                    &
                   -      (  8.0_wp * ibit26 * adv_mom_5                       &
                       +              ibit25 * adv_mom_3                       &
                          ) *                                                  &
                             ( v(k_pp,j,i)  + v(k-1,j,i)  )                    &
                   +      (           ibit26 * adv_mom_5                       &
                          ) *                                                  &
                             ( v(k_ppp,j,i) + v(k_mm,j,i) )                    &
                                  )

                diss_t = - ABS( w_comp ) * (                                   &
                          ( 10.0_wp * ibit26 * adv_mom_5                       &
                       +     3.0_wp * ibit25 * adv_mom_3                       &
                       +              ibit24 * adv_mom_1                       &
                          ) *                                                  &
                             ( v(k+1,j,i)   - v(k,j,i)    )                    &
                   -      (  5.0_wp * ibit26 * adv_mom_5                       &
                       +              ibit25 * adv_mom_3                       &
                          ) *                                                  &
                             ( v(k_pp,j,i)  - v(k-1,j,i)  )                    &
                   +      (           ibit26 * adv_mom_5                       &
                          ) *                                                  &
                             ( v(k_ppp,j,i) - v(k_mm,j,i) )                    &
                                           )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( u_comp       * ( ibit18 + ibit19 + ibit20 )          &
                - ( u(k,j-1,i)   + u(k,j,i) )                                  &
                                       * ( IBITS(advc_flags_1(k,j,i-1),18,1)   &
                                         + IBITS(advc_flags_1(k,j,i-1),19,1)   &
                                         + IBITS(advc_flags_1(k,j,i-1),20,1)   &
                                         )                                     &
                  ) * ddx                                                      &
               +  ( v_comp             * ( ibit21 + ibit22 + ibit23 )          &
                - ( v(k,j,i)     + v(k,j-1,i) )                                &
                                       * ( IBITS(advc_flags_1(k,j-1,i),21,1)   &
                                         + IBITS(advc_flags_1(k,j-1,i),22,1)   &
                                         + IBITS(advc_flags_1(k,j-1,i),23,1)   &
                                         )                                     &
                  ) * ddy                                                      &
               +  ( w_comp             * ( ibit24 + ibit25 + ibit26 )          &
                - ( w(k-1,j-1,i) + w(k-1,j,i) )                                &
                                       * ( IBITS(advc_flags_1(k-1,j,i),24,1)   &
                                         + IBITS(advc_flags_1(k-1,j,i),25,1)   &
                                         + IBITS(advc_flags_1(k-1,j,i),26,1)   &
                                         )                                     &
                   ) * ddzw(k)                                                 &
                ) * 0.5_wp


                tend(k,j,i) = tend(k,j,i) - (                                  &
                       ( flux_r + diss_r - flux_l - diss_l ) * ddx             &
                     + ( flux_n + diss_n - flux_s - diss_s ) * ddy             &
                     + ( flux_t + diss_t - flux_d - diss_d ) * ddzw(k)         &
                                            ) + v(k,j,i) * div

                flux_n_arr(k,j,i) = flux_n
                diss_n_arr(k,j,i) = diss_n
                flux_r_arr(k,j,i) = flux_r
                diss_r_arr(k,j,i) = diss_r
                flux_t_arr(k,j,i) = flux_t
                diss_t_arr(k,j,i) = diss_t

             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

       !$acc parallel present( sums_usvs_ws_l, sums_vs2_ws_l, sums_wsvs_ws_l ) &
       !$acc present( hom, weight_substep )
       !$acc loop
       DO  k = nzb+1, nzt
          tmp = 0.0_wp
          tmp2 = 0.0_wp
          tmp3 = 0.0_wp
          !$acc loop collapse(2) reduction(+:tmp,tmp2,tmp3)
          DO  i = nxl, nxr
             DO j = nys, nyn
                u_comp = u(k,j-1,i+1) + u(k,j,i+1)
                v_comp = v(k,j+1,i) + v(k,j,i)
                w_comp = w(k,j,i) + w(k,j-1,i)
!
!--             Statistical Evaluation of u'v'.
                tmp = tmp                                                      &
                     + ( flux_r_arr(k,j,i)                                     &
                       * ( u_comp - 2.0_wp * hom(k,1,1)            )           &
                       / ( u_comp + SIGN( 1.0E-20_wp, u_comp )     )           &
                     +   diss_r_arr(k,j,i)                                     &
                       *   ABS( u_comp - 2.0_wp * hom(k,1,1)       )           &
                       / ( ABS( u_comp ) + 1.0E-20_wp              )           &
                     ) *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of v'v'.
                tmp2 = tmp2                                                    &
                    + ( flux_n_arr(k,j,i)                                      &
                      * ( v_comp - 2.0_wp * hom(k,1,2)             )           &
                      / ( v_comp + SIGN( 1.0E-20_wp, v_comp ) )                &
                    +   diss_n_arr(k,j,i)                                      &
                      *   ABS( v_comp - 2.0_wp * hom(k,1,2)        )           &
                      / ( ABS( v_comp ) + 1.0E-20_wp               )           &
                      ) *   weight_substep(intermediate_timestep_count)
!
!--             Statistical Evaluation of w'v'.
                tmp3 = tmp3                                                    &
                     + ( flux_t_arr(k,j,i)                                     &
                       * ( w_comp - 2.0_wp * hom(k,1,3)            )           &
                       / ( w_comp + SIGN( 1.0E-20_wp, w_comp )     )           &
                     +   diss_t_arr(k,j,i)                                     &
                       *   ABS( w_comp - 2.0_wp * hom(k,1,3)       )           &
                       / ( ABS( w_comp ) + 1.0E-20_wp              )           &
                     ) *   weight_substep(intermediate_timestep_count)
             ENDDO
          ENDDO
          sums_usvs_ws_l(k,tn) = sums_usvs_ws_l(k,tn) + tmp
          sums_vs2_ws_l(k,tn)  = sums_vs2_ws_l(k,tn) + tmp2
          sums_wsvs_ws_l(k,tn) = sums_wsvs_ws_l(k,tn) + tmp3
       ENDDO
       sums_vs2_ws_l(nzb,tn) = sums_vs2_ws_l(nzb+1,tn)
       !$acc end parallel
       !$acc end data

    END SUBROUTINE advec_v_ws


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Advection of w - Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE advec_w_ws

       USE arrays_3d,                                                          &
           ONLY:  ddzu, tend, u, v, w

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5

       USE control_parameters,                                                 &
           ONLY:  intermediate_timestep_count

       USE grid_variables,                                                     &
           ONLY:  ddx, ddy

       USE indices,                                                            &
           ONLY:  nxl, nxr, nyn, nys, nzb, nzt, advc_flags_1,         &
                  advc_flags_2

       USE kinds

       USE statistics,                                                         &
           ONLY:  hom, sums_ws2_ws_l, weight_substep

       IMPLICIT NONE

       INTEGER(iwp) ::  i      !<
       INTEGER(iwp) ::  ibit27 !<
       INTEGER(iwp) ::  ibit28 !<
       INTEGER(iwp) ::  ibit29 !<
       INTEGER(iwp) ::  ibit30 !<
       INTEGER(iwp) ::  ibit31 !<
       INTEGER(iwp) ::  ibit32 !<
       INTEGER(iwp) ::  ibit33 !<
       INTEGER(iwp) ::  ibit34 !<
       INTEGER(iwp) ::  ibit35 !<
       INTEGER(iwp) ::  j      !<
       INTEGER(iwp) ::  k      !<
       INTEGER(iwp) ::  k_mm   !<
       INTEGER(iwp) ::  k_mmm   !<
       INTEGER(iwp) ::  k_pp   !<
       INTEGER(iwp) ::  k_ppp  !<
       INTEGER(iwp) ::  tn = 0 !<

       REAL(wp) ::  div    !<
       REAL(wp) ::  u_comp !<
       REAL(wp) ::  v_comp !<
       REAL(wp) ::  w_comp !<
       REAL(wp) ::  flux_n, flux_s, flux_r, flux_l, flux_t, flux_d !<
       REAL(wp) ::  diss_n, diss_s, diss_r, diss_l, diss_t, diss_d !<

       REAL(wp) ::  tmp    !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  diss_t_arr !<
       REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr)   ::  flux_t_arr !<

!
!--    Computation of interior fluxes and tendency terms

       !$acc data present( tend ) &
       !$acc present( u, v, w ) &
       !$acc create( flux_t_arr, diss_t_arr )

       !$acc parallel present ( advc_flags_1, advc_flags_2 ) &
       !$acc present( ddzu )
       !$acc loop
       DO i = nxl, nxr
          !$acc loop
          DO  j = nys, nyn
             !$acc loop
             DO  k = nzb+1, nzt

                ! left
                ibit29 = IBITS(advc_flags_1(k,j,i-1),29,1)
                ibit28 = IBITS(advc_flags_1(k,j,i-1),28,1)
                ibit27 = IBITS(advc_flags_1(k,j,i-1),27,1)

                u_comp = u(k+1,j,i) + u(k,j,i)
                flux_l = u_comp * (                                            &
                          ( 37.0_wp * ibit29 * adv_mom_5                       &
                       +     7.0_wp * ibit28 * adv_mom_3                       &
                       +              ibit27 * adv_mom_1                       &
                          ) *                                                  &
                         ( w(k,j,i)   + w(k,j,i-1) )                           &
                   -      (  8.0_wp * ibit29 * adv_mom_5                       &
                       +              ibit28 * adv_mom_3                       &
                          ) *                                                  &
                         ( w(k,j,i+1) + w(k,j,i-2) )                           &
                   +      (           ibit29 * adv_mom_5                       &
                          ) *                                                  &
                         ( w(k,j,i+2) + w(k,j,i-3) )                           &
                                  )

                diss_l = - ABS( u_comp ) * (                                   &
                          ( 10.0_wp * ibit29 * adv_mom_5                       &
                       +     3.0_wp * ibit28 * adv_mom_3                       &
                       +              ibit27 * adv_mom_1                       &
                          ) *                                                  &
                       ( w(k,j,i)   - w(k,j,i-1) )                             &
                   -      (  5.0_wp * ibit29 * adv_mom_5                       &
                       +              ibit28 * adv_mom_3                       &
                          ) *                                                  &
                       ( w(k,j,i+1) - w(k,j,i-2) )                             &
                   +      (           ibit29 * adv_mom_5                       &
                          ) *                                                  &
                       ( w(k,j,i+2) - w(k,j,i-3) )                             &
                                           )

                ! right
                ibit29 = IBITS(advc_flags_1(k,j,i),29,1)
                ibit28 = IBITS(advc_flags_1(k,j,i),28,1)
                ibit27 = IBITS(advc_flags_1(k,j,i),27,1)

                u_comp    = u(k+1,j,i+1) + u(k,j,i+1)
                flux_r = u_comp * (                                            &
                          ( 37.0_wp * ibit29 * adv_mom_5                       &
                       +     7.0_wp * ibit28 * adv_mom_3                       &
                       +              ibit27 * adv_mom_1                       &
                          ) *                                                  &
                                 ( w(k,j,i+1) + w(k,j,i)   )                   &
                   -      (  8.0_wp * ibit29 * adv_mom_5                       &
                       +              ibit28 * adv_mom_3                       &
                          ) *                                                  &
                                 ( w(k,j,i+2) + w(k,j,i-1) )                   &
                   +      (           ibit29 * adv_mom_5                       &
                          ) *                                                  &
                                 ( w(k,j,i+3) + w(k,j,i-2) )                   &
                                  )

                diss_r = - ABS( u_comp ) * (                                   &
                          ( 10.0_wp * ibit29 * adv_mom_5                       &
                       +     3.0_wp * ibit28 * adv_mom_3                       &
                       +              ibit27 * adv_mom_1                       &
                          ) *                                                  &
                                 ( w(k,j,i+1) - w(k,j,i)  )                    &
                   -      (  5.0_wp * ibit29 * adv_mom_5                       &
                       +              ibit28 * adv_mom_3                       &
                          ) *                                                  &
                                 ( w(k,j,i+2) - w(k,j,i-1) )                   &
                   +      (           ibit29 * adv_mom_5                       &
                          ) *                                                  &
                                 ( w(k,j,i+3) - w(k,j,i-2) )                   &
                                           )

                ! south
                ibit32 = IBITS(advc_flags_2(k,j-1,i),0,1)
                ibit31 = IBITS(advc_flags_1(k,j-1,i),31,1)
                ibit30 = IBITS(advc_flags_1(k,j-1,i),30,1)

                v_comp = v(k+1,j,i) + v(k,j,i)
                flux_s = v_comp * (                                            &
                          ( 37.0_wp * ibit32 * adv_mom_5                       &
                       +     7.0_wp * ibit31 * adv_mom_3                       &
                       +              ibit30 * adv_mom_1                       &
                          ) *                                                  &
                           ( w(k,j,i)   + w(k,j-1,i) )                         &
                   -      (  8.0_wp * ibit32 * adv_mom_5                       &
                       +              ibit31 * adv_mom_3                       &
                          ) *                                                  &
                           ( w(k,j+1,i) + w(k,j-2,i) )                         &
                   +      (           ibit32 * adv_mom_5                       &
                          ) *                                                  &
                           ( w(k,j+2,i) + w(k,j-3,i) )                         &
                                  )

                diss_s = - ABS( v_comp ) * (                                   &
                          ( 10.0_wp * ibit32 * adv_mom_5                       &
                       +     3.0_wp * ibit31 * adv_mom_3                       &
                       +              ibit30 * adv_mom_1                       &
                          ) *                                                  &
                           ( w(k,j,i)   - w(k,j-1,i) )                         &
                   -      (  5.0_wp * ibit32 * adv_mom_5                       &
                       +              ibit31 * adv_mom_3                       &
                          ) *                                                  &
                           ( w(k,j+1,i) - w(k,j-2,i) )                         &
                   +      (           ibit32 * adv_mom_5                       &
                          ) *                                                  &
                           ( w(k,j+2,i) - w(k,j-3,i) )                         &
                                           )

                ! north
                ibit32 = IBITS(advc_flags_2(k,j,i),0,1)
                ibit31 = IBITS(advc_flags_1(k,j,i),31,1)
                ibit30 = IBITS(advc_flags_1(k,j,i),30,1)

                v_comp    = v(k+1,j+1,i) + v(k,j+1,i)
                flux_n = v_comp * (                                            &
                          ( 37.0_wp * ibit32 * adv_mom_5                       &
                       +     7.0_wp * ibit31 * adv_mom_3                       &
                       +              ibit30 * adv_mom_1                       &
                          ) *                                                  &
                                 ( w(k,j+1,i) + w(k,j,i)   )                   &
                   -      (  8.0_wp * ibit32 * adv_mom_5                       &
                       +              ibit31 * adv_mom_3                       &
                          ) *                                                  &
                                 ( w(k,j+2,i) + w(k,j-1,i) )                   &
                   +      (           ibit32 * adv_mom_5                       &
                          ) *                                                  &
                                 ( w(k,j+3,i) + w(k,j-2,i) )                   &
                                  )

                diss_n = - ABS( v_comp ) * (                                   &
                          ( 10.0_wp * ibit32 * adv_mom_5                       &
                       +     3.0_wp * ibit31 * adv_mom_3                       &
                       +              ibit30 * adv_mom_1                       &
                          ) *                                                  &
                                 ( w(k,j+1,i) - w(k,j,i)  )                    &
                   -      (  5.0_wp * ibit32 * adv_mom_5                       &
                       +              ibit31 * adv_mom_3                       &
                          ) *                                                  &
                                 ( w(k,j+2,i) - w(k,j-1,i) )                   &
                   +      (           ibit32 * adv_mom_5                       &
                          ) *                                                  &
                                 ( w(k,j+3,i) - w(k,j-2,i) )                   &
                                           )
!
!--             k index has to be modified near bottom and top, else array
!--             subscripts will be exceeded.
                ! bottom
                ibit35 = IBITS(advc_flags_2(k-1,j,i),3,1)
                ibit34 = IBITS(advc_flags_2(k-1,j,i),2,1)
                ibit33 = IBITS(advc_flags_2(k-1,j,i),1,1)

                k_pp  = k + 2 * ibit35
                k_mm  = k - 2 * ( ibit34 + ibit35 )
                k_mmm = k - 3 * ibit35

                w_comp    = w(k,j,i) + w(k-1,j,i)
                flux_d = w_comp * (                                            &
                          ( 37.0_wp * ibit35 * adv_mom_5                       &
                       +     7.0_wp * ibit34 * adv_mom_3                       &
                       +              ibit33 * adv_mom_1                       &
                          ) *                                                  &
                             ( w(k,j,i)  + w(k-1,j,i)     )                    &
                   -      (  8.0_wp * ibit35 * adv_mom_5                       &
                       +              ibit34 * adv_mom_3                       &
                          ) *                                                  &
                             ( w(k+1,j,i)  + w(k_mm,j,i)  )                    &
                   +      (           ibit35 * adv_mom_5                       &
                          ) *                                                  &
                             ( w(k_pp,j,i) + w(k_mmm,j,i) )                    &
                                  )

                diss_d = - ABS( w_comp ) * (                                   &
                          ( 10.0_wp * ibit35 * adv_mom_5                       &
                       +     3.0_wp * ibit34 * adv_mom_3                       &
                       +              ibit33 * adv_mom_1                       &
                          ) *                                                  &
                             ( w(k,j,i)   - w(k-1,j,i)    )                    &
                   -      (  5.0_wp * ibit35 * adv_mom_5                       &
                       +              ibit34 * adv_mom_3                       &
                          ) *                                                  &
                             ( w(k+1,j,i)  - w(k_mm,j,i)  )                    &
                   +      (           ibit35 * adv_mom_5                       &
                          ) *                                                  &
                             ( w(k_pp,j,i) - w(k_mmm,j,i) )                    &
                                           )

                ! top
                ibit35 = IBITS(advc_flags_2(k,j,i),3,1)
                ibit34 = IBITS(advc_flags_2(k,j,i),2,1)
                ibit33 = IBITS(advc_flags_2(k,j,i),1,1)

                k_ppp = k + 3 * ibit35
                k_pp  = k + 2 * ( 1 - ibit33  )
                k_mm  = k - 2 * ibit35

                w_comp    = w(k+1,j,i) + w(k,j,i)
                flux_t = w_comp * (                                            &
                          ( 37.0_wp * ibit35 * adv_mom_5                       &
                       +     7.0_wp * ibit34 * adv_mom_3                       &
                       +              ibit33 * adv_mom_1                       &
                          ) *                                                  &
                             ( w(k+1,j,i)  + w(k,j,i)     )                    &
                   -      (  8.0_wp * ibit35 * adv_mom_5                       &
                       +              ibit34 * adv_mom_3                       &
                          ) *                                                  &
                             ( w(k_pp,j,i)  + w(k-1,j,i)  )                    &
                   +      (           ibit35 * adv_mom_5                       &
                          ) *                                                  &
                             ( w(k_ppp,j,i) + w(k_mm,j,i) )                    &
                                  )

                diss_t = - ABS( w_comp ) * (                                   &
                          ( 10.0_wp * ibit35 * adv_mom_5                       &
                       +     3.0_wp * ibit34 * adv_mom_3                       &
                       +              ibit33 * adv_mom_1                       &
                          ) *                                                  &
                             ( w(k+1,j,i)   - w(k,j,i)    )                    &
                   -      (  5.0_wp * ibit35 * adv_mom_5                       &
                       +              ibit34 * adv_mom_3                       &
                          ) *                                                  &
                             ( w(k_pp,j,i)  - w(k-1,j,i)  )                    &
                   +      (           ibit35 * adv_mom_5                       &
                          ) *                                                  &
                             ( w(k_ppp,j,i) - w(k_mm,j,i) )                    &
                                           )
!
!--             Calculate the divergence of the velocity field. A respective
!--             correction is needed to overcome numerical instabilities caused
!--             by a not sufficient reduction of divergences near topography.
                div = ( ( u_comp    * ( ibit27 + ibit28 + ibit29 )            &
                  - ( u(k+1,j,i) + u(k,j,i)   )                               &
                                    * ( IBITS(advc_flags_1(k,j,i-1),27,1)     &
                                      + IBITS(advc_flags_1(k,j,i-1),28,1)     &
                                      + IBITS(advc_flags_1(k,j,i-1),29,1)     &
                                      )                                       &
                  ) * ddx                                                     &
              +   ( v_comp          * ( ibit30 + ibit31 + ibit32 )            &
                  - ( v(k+1,j,i) + v(k,j,i)   )                               &
                                    * ( IBITS(advc_flags_1(k,j-1,i),30,1)     &
                                      + IBITS(advc_flags_1(k,j-1,i),31,1)     &
                                      + IBITS(advc_flags_2(k,j-1,i),0,1)      &
                                      )                                       &
                  ) * ddy                                                     &
              +   ( w_comp          * ( ibit33 + ibit34 + ibit35 )            &
                - ( w(k,j,i)   + w(k-1,j,i)   )                               &
                                    * ( IBITS(advc_flags_2(k-1,j,i),1,1)      &
                                      + IBITS(advc_flags_2(k-1,j,i),2,1)      &
                                      + IBITS(advc_flags_2(k-1,j,i),3,1)      &
                                      )                                       &
                  )  * ddzu(k+1)                                              &
                ) * 0.5_wp



                tend(k,j,i) = tend(k,j,i) - (                                 &
                      ( flux_r + diss_r - flux_l - diss_l ) * ddx             &
                    + ( flux_n + diss_n - flux_s - diss_s ) * ddy             &
                    + ( flux_t + diss_t - flux_d - diss_d ) * ddzu(k+1)       &
                                            )  + div * w(k,j,i)

                flux_t_arr(k,j,i) = flux_t
                diss_t_arr(k,j,i) = diss_t

             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

       !$acc parallel present( sums_ws2_ws_l ) &
       !$acc present( hom, weight_substep )
       !$acc loop
       DO  k = nzb+1, nzt
          tmp = 0.0_wp
          !$acc loop collapse(2) reduction(+:tmp)
          DO  i = nxl, nxr
             DO j = nys, nyn
                w_comp = w(k,j,i) + w(k+1,j,i)
!
!--             Statistical Evaluation of w'w'.
                tmp = tmp                                                      &
                     + ( flux_t_arr(k,j,i)                                     &
                       * ( w_comp - 2.0_wp * hom(k,1,3)          )             &
                       / ( w_comp + SIGN( 1.0E-20_wp, w_comp )     )           &
                     +   diss_t_arr(k,j,i)                                     &
                       *   ABS( w_comp - 2.0_wp * hom(k,1,3)     )             &
                       / ( ABS( w_comp ) + 1.0E-20_wp              )           &
                     ) *   weight_substep(intermediate_timestep_count)
             ENDDO
          ENDDO
          sums_ws2_ws_l(k,tn) = sums_ws2_ws_l(k,tn) + tmp
       ENDDO
       !$acc end parallel
       !$acc end data

    END SUBROUTINE advec_w_ws


    SUBROUTINE ws_finalize

       USE constants,                                                          &
           ONLY:  adv_mom_1, adv_mom_3, adv_mom_5, adv_sca_1, adv_sca_3,       &
                  adv_sca_5

       USE indices,                                                            &
           ONLY:  nyn, nys, nzb, nzt

       USE kinds

       USE pegrid

       USE statistics,                                                         &
           ONLY:  sums_us2_ws_l, sums_vs2_ws_l, sums_ws2_ws_l,                 &
                  sums_wspts_ws_l, sums_wssas_ws_l,                            &
                  sums_wsus_ws_l, sums_wsvs_ws_l

!--    Arrays needed for statical evaluation of fluxes.

       DEALLOCATE( sums_wsus_ws_l, sums_wsvs_ws_l,                             &
                   sums_us2_ws_l, sums_vs2_ws_l,                               &
                   sums_ws2_ws_l )
       DEALLOCATE( sums_wspts_ws_l )
       DEALLOCATE( sums_wssas_ws_l )

    END SUBROUTINE ws_finalize


 END MODULE advec_ws
