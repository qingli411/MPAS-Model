!> @file sum_up_3d_data.f90
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
! $Id: sum_up_3d_data.f90 3004 2018-04-27 12:33:25Z Giersch $
! prr field added to ONLY-list, prr* case/pr* case/precipitation_rate_av
! removed, further allocation checks implemented
!
! 2963 2018-04-12 14:47:44Z suehring
! Introduce index for vegetation/wall, pavement/green-wall and water/window
! surfaces, for clearer access of surface fraction, albedo, emissivity, etc. .
!
! 2894 2018-03-15 09:17:58Z Giersch
! Changed comment
!
! 2817 2018-02-19 16:32:21Z suehring
! Preliminary gust module interface implemented
!
! 2798 2018-02-09 17:16:39Z suehring
! Consider also default-type surfaces for surface temperature output.
!
! 2797 2018-02-08 13:24:35Z suehring
! Enable output of ground-heat flux also at urban surfaces.
!
! 2790 2018-02-06 11:57:19Z suehring
! Bugfix in summation of surface sensible and latent heat flux
!
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
!
! 2743 2018-01-12 16:03:39Z suehring
! In case of natural- and urban-type surfaces output surfaces fluxes in W/m2.
!
! 2742 2018-01-12 14:59:47Z suehring
! Enable output of surface temperature
!
! 2735 2018-01-11 12:01:27Z suehring
! output of r_a moved from land-surface to consider also urban-type surfaces
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! - Change in file header (GPL part)
! - Implementation of uv exposure model (FK)
! - output of diss_av, kh_av, km_av (turbulence_closure_mod) (TG)
! - Implementation of chemistry module (FK)
! - Workaround for sum-up usm arrays in case of restart runs, to avoid program
!   crash (MS)
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new surface concept
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
!
! 2024 2016-10-12 16:42:37Z kanani
! Added missing CASE for ssws*
!
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! changed prefix for urban surface model output to "usm_",
! introduced control parameter varnamelength for LEN of trimvar.
!
! 2007 2016-08-24 15:47:17Z kanani
! Added support for new urban surface model (temporary modifications of
! SELECT CASE ( ) necessary, see variable trimvar),
! added comments in variable declaration section
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1992 2016-08-12 15:14:59Z suehring
! Bugfix in summation of passive scalar
!
! 1976 2016-07-27 13:28:04Z maronga
! Radiation actions are now done directly in the respective module
!
! 1972 2016-07-26 07:52:02Z maronga
! Land surface actions are now done directly in the respective module
!
! 1960 2016-07-12 16:34:24Z suehring
! Scalar surface flux added
!
! 1949 2016-06-17 07:19:16Z maronga
! Bugfix: calculation of lai_av, c_veg_av and c_liq_av.
!
! 1849 2016-04-08 11:33:18Z hoffmann
! precipitation_rate moved to arrays_3d
!
! 1788 2016-03-10 11:01:04Z maronga
! Added z0q and z0q_av
!
! 1693 2015-10-27 08:35:45Z maronga
! Last revision text corrected
!
! 1691 2015-10-26 16:17:44Z maronga
! Added output of Obukhov length and radiative heating rates for RRTMG.
! Corrected output of liquid water path.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1585 2015-04-30 07:05:52Z maronga
! Adapted for RRTMG
!
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s
!
! 1551 2015-03-03 14:18:16Z maronga
! Added support for land surface model and radiation model data.
!
! 1359 2014-04-11 17:15:14Z hoffmann
! New particle structure integrated.
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
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
! 1318 2014-03-17 13:35:16Z raasch
! barrier argument removed from cpu_log,
! module interfaces removed
!
! 1115 2013-03-26 18:16:16Z hoffmann
! ql is calculated by calc_liquid_water_content
!
! 1053 2012-11-13 17:11:03Z hoffmann
! +nr, prr, qr
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1007 2012-09-19 14:30:36Z franke
! Bugfix in calculation of ql_vp
!
! 978 2012-08-09 08:28:32Z fricke
! +z0h*
!
! Revision 1.1  2006/02/23 12:55:23  raasch
! Initial revision
!
!
! Description:
! ------------
!> Sum-up the values of 3d-arrays. The real averaging is later done in routine
!> average_3d_data.
!------------------------------------------------------------------------------!
 SUBROUTINE sum_up_3d_data


    USE arrays_3d,                                                             &
        ONLY:  dzw, e, p, pt,         &
               rho_ocean, sa, u, v, w,      &
               alpha_T, beta_S, solar3d

    USE averaging,                                                             &
        ONLY:  e_av, kh_av, km_av,      &
               p_av, pt_av,    &
               rho_ocean_av, sa_av,    &
               u_av, v_av, w_av,          &
               alpha_T_av, beta_S_av, solar3d_av

    USE control_parameters,                                                    &
        ONLY:  average_count_3d, doav, doav_n,   &
               varnamelength

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

        USE indices,                                                               &
        ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

    USE kinds


    IMPLICIT NONE

    INTEGER(iwp) ::  i   !< grid index x direction
    INTEGER(iwp) ::  ii  !< running index
    INTEGER(iwp) ::  j   !< grid index y direction
    INTEGER(iwp) ::  k   !< grid index x direction
    INTEGER(iwp) ::  m   !< running index surface type
    INTEGER(iwp) ::  n   !<

    REAL(wp)     ::  mean_r !<
    REAL(wp)     ::  s_r2   !<
    REAL(wp)     ::  s_r3   !<

    CHARACTER (LEN=varnamelength) ::  trimvar  !< TRIM of output-variable string


    CALL cpu_log (log_point(34),'sum_up_3d_data','start')

!
!-- Allocate and initialize the summation arrays if called for the very first
!-- time or the first time after average_3d_data has been called
!-- (some or all of the arrays may have been already allocated
!-- in rrd_local)
    IF ( average_count_3d == 0 )  THEN

       DO  ii = 1, doav_n
!
!--       Temporary solution to account for data output within the new urban
!--       surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar )
          trimvar = TRIM( doav(ii) )

          SELECT CASE ( trimvar )

             CASE ( 'e' )
                IF ( .NOT. ALLOCATED( e_av ) )  THEN
                   ALLOCATE( e_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                e_av = 0.0_wp

             CASE ( 'p' )
                IF ( .NOT. ALLOCATED( p_av ) )  THEN
                   ALLOCATE( p_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                p_av = 0.0_wp

             CASE ( 'pt' )
                IF ( .NOT. ALLOCATED( pt_av ) )  THEN
                   ALLOCATE( pt_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                pt_av = 0.0_wp

            CASE ( 'solar3d' )
                IF ( .NOT. ALLOCATED( solar3d_av ) )  THEN
                   ALLOCATE( solar3d_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                solar3d_av = 0.0_wp

             CASE ( 'rho_ocean' )
                IF ( .NOT. ALLOCATED( rho_ocean_av ) )  THEN
                   ALLOCATE( rho_ocean_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                rho_ocean_av = 0.0_wp

             CASE ( 'alpha_T' )
                IF ( .NOT. ALLOCATED( alpha_T_av ) )  THEN
                   ALLOCATE( alpha_T_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                alpha_T_av = 0.0_wp

             CASE ( 'beta_S' )
                IF ( .NOT. ALLOCATED( beta_S_av ) )  THEN
                   ALLOCATE( beta_S_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                beta_S_av = 0.0_wp

             CASE ( 'sa' )
                IF ( .NOT. ALLOCATED( sa_av ) )  THEN
                   ALLOCATE( sa_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                sa_av = 0.0_wp

             CASE ( 'u' )
                IF ( .NOT. ALLOCATED( u_av ) )  THEN
                   ALLOCATE( u_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                u_av = 0.0_wp

             CASE ( 'v' )
                IF ( .NOT. ALLOCATED( v_av ) )  THEN
                   ALLOCATE( v_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                v_av = 0.0_wp

             CASE ( 'w' )
                IF ( .NOT. ALLOCATED( w_av ) )  THEN
                   ALLOCATE( w_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                ENDIF
                w_av = 0.0_wp

!
!--          Block of urban surface model outputs

             CASE DEFAULT

                CONTINUE

          END SELECT

       ENDDO

    ENDIF

!
!-- Loop of all variables to be averaged.
    DO  ii = 1, doav_n
!
!--       Temporary solution to account for data output within the new urban
!--       surface model (urban_surface_mod.f90), see also SELECT CASE ( trimvar )
          trimvar = TRIM( doav(ii) )
!
!--    Store the array chosen on the temporary array.
       SELECT CASE ( trimvar )
          CASE ( 'e' )
             IF ( ALLOCATED( e_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         e_av(k,j,i) = e_av(k,j,i) + e(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'p' )
             IF ( ALLOCATED( p_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         p_av(k,j,i) = p_av(k,j,i) + p(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'pt' )
             IF ( ALLOCATED( pt_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                            pt_av(k,j,i) = pt_av(k,j,i) + pt(k,j,i)
                         ENDDO
                      ENDDO
                   ENDDO
             ENDIF

          CASE ( 'solar3d' )
             IF ( ALLOCATED( solar3d_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         solar3d_av(k,j,i) = solar3d_av(k,j,i) + solar3d(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF


          CASE ( 'rho_ocean' )
             IF ( ALLOCATED( rho_ocean_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         rho_ocean_av(k,j,i) = rho_ocean_av(k,j,i) + rho_ocean(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'alpha_T' )
             IF ( ALLOCATED( alpha_T_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         alpha_T_av(k,j,i) = alpha_T_av(k,j,i) + alpha_T(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'beta_S' )
             IF ( ALLOCATED( beta_S_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         beta_S_av(k,j,i) = beta_S_av(k,j,i) + beta_S(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'sa' )
             IF ( ALLOCATED( sa_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         sa_av(k,j,i) = sa_av(k,j,i) + sa(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'u' )
             IF ( ALLOCATED( u_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         u_av(k,j,i) = u_av(k,j,i) + u(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'v' )
             IF ( ALLOCATED( v_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         v_av(k,j,i) = v_av(k,j,i) + v(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'w' )
             IF ( ALLOCATED( w_av ) ) THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         w_av(k,j,i) = w_av(k,j,i) + w(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE DEFAULT
!
             CONTINUE

       END SELECT

    ENDDO

    CALL cpu_log( log_point(34), 'sum_up_3d_data', 'stop' )


 END SUBROUTINE sum_up_3d_data
