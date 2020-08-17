!> @file turbulence_closure_mod.f90
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
! Copyright 2017-2018 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
!
!
! Former revisions:
! -----------------
! $Id: turbulence_closure_mod.f90 3086 2018-06-25 09:08:04Z gronemeier $
! bugfix: set rans_const_sigma(1) = 1.3
!
! 3083 2018-06-19 14:03:12Z gronemeier
! - set limits of diss at the end of prognostic equations
! - call production_e to calculate production term of diss
! - limit change of diss to -90% to +100%
! - remove factor 0.5 from diffusion_diss_ij
! - rename c_m into c_0, and c_h into c_4
! - add rans_const_c and rans_const_sigma as namelist parameters
! - add calculation of mixing length for profile output in case of rans_tke_e
! - changed format of annotations to comply with doxygen standards
! - calculate and save dissipation rate during rans_tke_l mode
! - set bc at vertical walls for e, diss, km, kh
! - bugfix: set l_wall = 0.0 within buildings
! - set l_wall at bottom and top boundary (rans-mode)
! - bugfix in production term for dissipation rate
! - bugfix in diffusion of dissipation rate
! - disable check for 1D model if rans_tke_e is used
! - bugfixes for initialization (rans-mode):
!    - correction of dissipation-rate formula
!    - calculate km based on l_wall
!    - initialize diss if 1D model is not used
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error message revised
!
! 3014 2018-05-09 08:42:38Z maronga
! Bugfix: nzb_do and nzt_do were not used for 3d data output
!
! 3004 2018-04-27 12:33:25Z Giersch
! Further allocation checks implemented
!
! 2938 2018-03-27 15:52:42Z suehring
! Further todo's
!
! 2936 2018-03-27 14:49:27Z gronemeier
! - defined l_grid only within this module
! - Moved l_wall definition from modules.f90
! - Get level of highest topography, used to limit upward distance calculation
! - Consider cyclic boundary conditions for mixing length calculation
! - Moved copy of wall_flags into subarray to subroutine
! - Implemented l_wall calculation in case of RANS simulation
! - Moved init of l_black to tcm_init_mixing_length
! - Moved init_mixing_length from init_grid.f90 and
!   renamed it to tcm_init_mixing_length
!
! 2764 2018-01-22 09:25:36Z gronemeier
! Bugfix: remove duplicate SAVE statements
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
! Initial revision
!
!
!
!
! Authors:
! --------
! @author Tobias Gronemeier
!
!
! Description:
! ------------
!> This module contains the available turbulence closures for PALM.
!>
!>
!> @todo test initialization for all possibilities
!>       add OpenMP directives whereever possible
!>       remove debug output variables (dummy1, dummy2, dummy3)
!> @todo Check for random disturbances
!> @note <Enter notes on the module>
!> @bug  TKE-e closure still crashes due to too small dt
!------------------------------------------------------------------------------!
 MODULE turbulence_closure_mod


#if defined( __nopointer )
    USE arrays_3d,                                                             &
        ONLY:  dd2zu, dzu, e, e_p, kh, km, prho, te_m, tend, u, v, w,          &
               u_stk, v_stk, sgs_diss
#else
    USE arrays_3d,                                                             &
        ONLY:  dd2zu, dzu, e, e_1, e_2, e_3, e_p, kh, km, prho,                &
               te_m, tend, u, v, w, u_stk, v_stk, sgs_diss
#endif

    USE control_parameters,                                                    &
        ONLY:  dt_3d, e_init, g, initializing_actions,                         &
               intermediate_timestep_count, intermediate_timestep_count_max,   &
               kappa, les_mw, prandtl_number, simulated_time, stokes_force,    &
               turbulence_closure, wall_adjustment, wall_adjustment_factor

    USE advec_ws,                                                              &
        ONLY:  advec_s_ws

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

    USE kinds

    USE pegrid

    USE stokes_force_mod,                                                      &
        ONLY:  stokes_force_s, stokes_production_e

    IMPLICIT NONE


    REAL(wp) ::  c_0                !< constant used for diffusion coefficient and dissipation (dependent on mode RANS/LES)
    REAL(wp) ::  dsig_e = 1.0_wp    !< factor to calculate Ke from Km (1/sigma_e)

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  l_grid     !< geometric mean of grid sizes dx, dy, dz

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  l_wall !< near-wall mixing length
    REAL(wp)     :: l               !< mixing length
    REAL(wp)     :: l_stable        !< mixing length according to stratification
    REAL(wp)     :: ll              !< adjusted l_grid


    PUBLIC c_0

!
!-- PALM interfaces:
!-- Input parameter checks to be done in check_parameters
    INTERFACE tcm_check_parameters
       MODULE PROCEDURE tcm_check_parameters
    END INTERFACE tcm_check_parameters

!
!-- Initialization actions
    INTERFACE tcm_init
       MODULE PROCEDURE tcm_init
    END INTERFACE tcm_init

!
!-- deallocate arrays
    INTERFACE tcm_deallocate_arrays
      MODULE PROCEDURE tcm_deallocate_arrays
    END INTERFACE tcm_deallocate_arrays

!
!-- Initialization of arrays
    INTERFACE tcm_init_arrays
       MODULE PROCEDURE tcm_init_arrays
    END INTERFACE tcm_init_arrays

!
!-- Prognostic equations for TKE and TKE dissipation rate
    INTERFACE tcm_prognostic
       MODULE PROCEDURE tcm_prognostic
    END INTERFACE tcm_prognostic

!
!-- Calculate diffusivities
    INTERFACE tcm_diffusivities
       MODULE PROCEDURE tcm_diffusivities
    END INTERFACE tcm_diffusivities

    SAVE

    PRIVATE
!
!-- Add INTERFACES that must be available to other modules (alphabetical order)
    PUBLIC tcm_check_parameters, tcm_deallocate_arrays,                        &
           tcm_diffusivities, tcm_init, tcm_init_arrays, tcm_prognostic,       &
           l_grid, l_wall


 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check parameters routine for turbulence closure module.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_check_parameters

    USE control_parameters,                                                    &
        ONLY:  message_string

    IMPLICIT NONE

!
!-- Define which turbulence closure is going to be used

    c_0 = 0.1_wp !according to Lilly (1967) and Deardorff (1980)

    dsig_e = 1.0_wp !assure to use K_m to calculate TKE instead
                    !of K_e which is used in RANS mode

    SELECT CASE ( TRIM( turbulence_closure ) )

       CASE ( 'Moeng_Wyngaard' )
          les_mw = .TRUE.

       CASE DEFAULT
          !> @todo rework this part so that only one call of this error exists
          message_string = 'Unknown turbulence closure: ' //                &
                           TRIM( turbulence_closure )
          CALL message( 'tcm_check_parameters', 'PA0500', 1, 2, 0, 6, 0 )

    END SELECT


 END SUBROUTINE tcm_check_parameters

 SUBROUTINE tcm_deallocate_arrays
    IMPLICIT NONE

    deallocate(kh, km, sgs_diss)
    deallocate(l_grid, l_wall)
#if defined( __nopointer )
    deallocate(e,e_p,te_m)
#else
    deallocate(e_1,e_2,e_3)
#endif

 END SUBROUTINE tcm_deallocate_arrays

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate arrays and assign pointers.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init_arrays
    IMPLICIT NONE

    ALLOCATE( kh(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( km(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( l_grid(1:nzt) )
    ALLOCATE( l_wall(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( sgs_diss(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

#if defined( __nopointer )
    ALLOCATE( e(nzb:nzt+1,nysg:nyng,nxlg:nxrg)    )
    ALLOCATE( e_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  )
    ALLOCATE( te_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

#else
    ALLOCATE( e_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( e_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( e_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
#endif
!
!-- Allocate arrays required for dissipation.
!-- Please note, if it is a nested run, arrays need to be allocated even if
!-- they do not necessarily need to be transferred, which is attributed to
!-- the design of the model coupler which allocates memory for each variable.
#if ! defined( __nopointer )
!
!-- Initial assignment of pointers
    e  => e_1;   e_p  => e_2;   te_m  => e_3

#endif

 END SUBROUTINE tcm_init_arrays


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of turbulence closure module.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init

    IMPLICIT NONE

    INTEGER(iwp) :: i            !< loop index
    INTEGER(iwp) :: j            !< loop index
    INTEGER(iwp) :: k            !< loop index

!
!-- Initialize mixing length
    CALL tcm_init_mixing_length

!
!-- Actions for initial runs
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN

       IF ( INDEX(initializing_actions, 'set_constant_profiles') /= 0 .OR. &
                INDEX( initializing_actions, 'inifor' ) /= 0 )  THEN

          IF ( e_init > 0.0_wp )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb+1, nzt
                      km(k,j,i) = c_0 * l_wall(k,j,i) * SQRT( e_init )
                   ENDDO
                ENDDO
             ENDDO
             km(nzb,:,:)   = km(nzb+1,:,:)
             km(nzt+1,:,:) = km(nzt,:,:)
             kh = km / prandtl_number
             e  = e_init
          ELSE
             ! there must exist an initial diffusion, because
             ! otherwise no TKE would be produced by the
             ! production terms, as long as not yet
             ! e = (u*/cm)**2 at k=nzb+1
             kh   = 0.00001_wp
             km   = 0.00001_wp
             e    = 0.0_wp
          ENDIF

       ENDIF
!
!--    Initialize old and new time levels.
       te_m = 0.0_wp
       e_p = e

    ELSEIF ( TRIM( initializing_actions ) == 'read_restart_data'  .OR.         &
             TRIM( initializing_actions ) == 'cyclic_fill' )                   &
    THEN

!--    Initialize new time levels (only done in order to set boundary values
!--    including ghost points)
       e_p = e
!
!--    Allthough tendency arrays are set in prognostic_equations, they have
!--    to be predefined here because there they are used (but multiplied with 0)
!--    before they are set.
       te_m = 0.0_wp

    ENDIF

 END SUBROUTINE tcm_init


! Description:
! -----------------------------------------------------------------------------!
!> Pre-computation of grid-dependent and near-wall mixing length.
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_init_mixing_length

    USE arrays_3d,                                                             &
        ONLY:  dzw

    USE control_parameters,                                                    &
        ONLY:  message_string, wall_adjustment_factor

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nzb, nzt

    USE kinds


    IMPLICIT NONE

    INTEGER(iwp) :: k              !< index variable along z

!
!-- Initialize the mixing length in case of an LES-simulation
!
!-- Compute the grid-dependent mixing length.
    DO  k = 1, nzt
       l_grid(k)  = ( dx * dy * dzw(k) )**0.33333333333333_wp
    ENDDO
!
!-- Initialize near-wall mixing length l_wall only in the vertical direction
!-- for the moment, multiplication with wall_adjustment_factor further below
    l_wall(nzb,:,:)   = l_grid(1)
    DO  k = nzb+1, nzt
       l_wall(k,:,:)  = l_grid(k)
    ENDDO
    l_wall(nzt+1,:,:) = l_grid(nzt)

    DO  k = 1, nzt
       IF ( l_grid(k) > 1.5_wp * dx * wall_adjustment_factor .OR.            &
            l_grid(k) > 1.5_wp * dy * wall_adjustment_factor )  THEN
          WRITE( message_string, * ) 'grid anisotropy exceeds ',             &
                                     'threshold given by only local',        &
                                     ' &horizontal reduction of near_wall ', &
                                     'mixing length l_wall',                 &
                                     ' &starting from height level k = ', k, &
                                     '.'
          CALL message( 'init_grid', 'PA0202', 0, 1, 0, 6, 0 )
          EXIT
       ENDIF
    ENDDO
!
!-- Set lateral boundary conditions for l_wall
    !$acc data copy( l_wall )
    CALL exchange_horiz( l_wall, nbgp )
    !$acc end data

 END SUBROUTINE tcm_init_mixing_length


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Prognostic equation for subgrid-scale TKE and TKE dissipation rate.
!> Vector-optimized version
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_prognostic

    USE arrays_3d,                                                             &
        ONLY:  ddzu, ddzw

    USE control_parameters,                                                    &
        ONLY:  tsc, g, prho_reference, wall_adjustment, wall_adjustment_factor

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy, ddx2, ddy2

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< loop index
    INTEGER(iwp) ::  j       !< loop index
    INTEGER(iwp) ::  k       !< loop index

    REAL(wp)     ::  dvar_dz        !< vertical gradient of var
    REAL(wp)     ::  l              !< mixing length
    REAL(wp)     ::  ll             !< adjusted l
    REAL(wp)     ::  l_stable       !< mixing length according to stratification
    REAL(wp)     ::  def

    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr) ::  dudx, dudy, dudz
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr) ::  dvdx, dvdy, dvdz
    REAL(wp), DIMENSION(nzb+1:nzt,nys:nyn,nxl:nxr) ::  dwdx, dwdy, dwdz

!
!-- If required, compute prognostic equation for turbulent kinetic
!-- energy (TKE)

    CALL cpu_log( log_point(16), 'tke-equation', 'start' )

    !$acc data create( tend )

    !$acc parallel present( tend )
    !$acc loop collapse(3)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb, nzt
             tend(k,j,i) = 0.0_wp
          ENDDO
       ENDDO
    ENDDO
    !$acc end parallel

    CALL advec_s_ws( e, 'e' )

    ! Compute Stokes-advection if required
    IF ( stokes_force ) THEN
       CALL stokes_force_s( e )
    ENDIF

!-- TKE production
!   Inline subroutine production_e()

    !$acc parallel present( g ) &
    !$acc present( tend ) &
    !$acc present( e, e_p ) &
    !$acc present( u, v, w ) &
    !$acc present( dd2zu, ddzw ) &
    !$acc present( km, kh, prho ) &
    !$acc create( dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz )
    !$acc loop
    DO  i = nxl, nxr
       !$acc loop
       DO  j = nys, nyn
          !$acc loop
          DO  k = nzb+1, nzt
             !-- Calculate TKE production by shear. Here, no additional
             !-- wall-bounded code is considered.

             dudx(k,j,i)  =           ( u(k,j,i+1) - u(k,j,i)     ) * ddx
             dudy(k,j,i)  = 0.25_wp * ( u(k,j+1,i) + u(k,j+1,i+1) -           &
                                      u(k,j-1,i) - u(k,j-1,i+1) ) * ddy
             dudz(k,j,i)  = 0.5_wp  * ( u(k+1,j,i) + u(k+1,j,i+1) -           &
                                      u(k-1,j,i) - u(k-1,j,i+1) ) *           &
                                                                  dd2zu(k)

             dvdx(k,j,i)  = 0.25_wp * ( v(k,j,i+1) + v(k,j+1,i+1) -           &
                                      v(k,j,i-1) - v(k,j+1,i-1) ) * ddx
             dvdy(k,j,i)  =           ( v(k,j+1,i) - v(k,j,i)     ) * ddy
             dvdz(k,j,i)  = 0.5_wp  * ( v(k+1,j,i) + v(k+1,j+1,i) -           &
                                      v(k-1,j,i) - v(k-1,j+1,i) ) *           &
                                                                  dd2zu(k)

             dwdx(k,j,i)  = 0.25_wp * ( w(k,j,i+1) + w(k-1,j,i+1) -           &
                                      w(k,j,i-1) - w(k-1,j,i-1) ) * ddx
             dwdy(k,j,i)  = 0.25_wp * ( w(k,j+1,i) + w(k-1,j+1,i) -           &
                                      w(k,j-1,i) - w(k-1,j-1,i) ) * ddy
             dwdz(k,j,i)  =           ( w(k,j,i)   - w(k-1,j,i)   ) *         &
                                                                  ddzw(k)

             def = 2.0_wp * (                                                 &
                          dudx(k,j,i)**2 + dvdy(k,j,i)**2 + dwdz(k,j,i)**2    &
                            ) +                                               &
                          dudy(k,j,i)**2 + dvdx(k,j,i)**2 + dwdx(k,j,i)**2 +  &
                          dwdy(k,j,i)**2 + dudz(k,j,i)**2 + dvdz(k,j,i)**2 +  &
                   2.0_wp * (                                                 &
                          dvdx(k,j,i)*dudy(k,j,i) + dwdx(k,j,i)*dudz(k,j,i) + &
                          dwdy(k,j,i)*dvdz(k,j,i)                             &
                            )

             IF ( def < 0.0_wp )  def = 0.0_wp

             tend(k,j,i) = tend(k,j,i) + km(k,j,i) * def

             !-- TKE production by buoyancy
             tend(k,j,i) = tend(k,j,i) +                                       &
                           kh(k,j,i) * g / prho_reference *                    &
                           ( prho(k+1,j,i) - prho(k-1,j,i) ) *                 &
                           dd2zu(k)
          ENDDO
       ENDDO
    ENDDO
    !$acc end parallel

!-- end inline subroutine production_e()


!-- Compute Stokes production term in e equation
    IF ( stokes_force ) THEN
       CALL stokes_production_e
    ENDIF

!
!-- Calculate the tendency terms due to diffusion
!   Inline subroutine diffusion_e()

    !$acc parallel present( g ) &
    !$acc present( tend ) &
    !$acc present( e ) &
    !$acc present( dd2zu, ddzu, ddzw ) &
    !$acc present( l_grid, l_wall) &
    !$acc present( km, sgs_diss, prho )

    !$acc loop
    DO  i = nxl, nxr
       !$acc loop
       DO  j = nys, nyn
          !$acc loop
          DO  k = nzb+1, nzt
    !
    !-- Determine the mixing length for LES closure
    !   Inline subroutine mixing_length_les()
             dvar_dz = -1.0_wp * (prho(k+1,j,i) - prho(k-1,j,i) ) * dd2zu(k)
             IF ( dvar_dz > 0.0_wp ) THEN
                l_stable = 0.76_wp * SQRT( e(k,j,i) )                                &
                                   / SQRT( g / prho_reference * dvar_dz ) + 1E-5_wp
             ELSE
                l_stable = l_grid(k)
             ENDIF
         !
         !-- Adjustment of the mixing length
             IF ( wall_adjustment )  THEN
                l  = MIN( wall_adjustment_factor * l_wall(k,j,i), l_grid(k), l_stable )
                ll = MIN( wall_adjustment_factor * l_wall(k,j,i), l_grid(k) )
             ELSE
                l  = MIN( l_grid(k), l_stable )
                ll = l_grid(k)
             ENDIF
    !-- end of inline subroutine mixing_length_les()

             ! SGS dissipation
             sgs_diss(k,j,i) = ( 0.19_wp + 0.74_wp * l / ll )                  &
                               * e(k,j,i) * SQRT( e(k,j,i) ) / l
             tend(k,j,i) = tend(k,j,i) + (                                     &
                                           (                                   &
                       ( km(k,j,i)+km(k,j,i+1) ) * ( e(k,j,i+1)-e(k,j,i) )     &
                     - ( km(k,j,i)+km(k,j,i-1) ) * ( e(k,j,i)-e(k,j,i-1) )     &
                                           ) * ddx2                            &
                                         + (                                   &
                       ( km(k,j,i)+km(k,j+1,i) ) * ( e(k,j+1,i)-e(k,j,i) )     &
                     - ( km(k,j,i)+km(k,j-1,i) ) * ( e(k,j,i)-e(k,j-1,i) )     &
                                           ) * ddy2                            &
                                         + (                                   &
            ( km(k,j,i)+km(k+1,j,i) ) * ( e(k+1,j,i)-e(k,j,i) ) * ddzu(k+1)    &
          - ( km(k,j,i)+km(k-1,j,i) ) * ( e(k,j,i)-e(k-1,j,i) ) * ddzu(k)      &
                                           ) * ddzw(k)                         &
                                         ) * dsig_e                            &
                                       ! dissipation
                                       - sgs_diss(k,j,i)

          ENDDO
       ENDDO
    ENDDO
    !$acc end parallel

!-- end of inline subroutine diffusion_e()

!
!-- Prognostic equation for TKE.
!-- Eliminate negative TKE values, which can occur due to numerical
!-- reasons in the course of the integration. In such cases the old TKE
!-- value is reduced by 90%.
    !$acc parallel present( tend ) &
    !$acc present(te_m, tsc) &
    !$acc present( e, e_p )

    !$acc loop collapse(3)
    DO  i = nxl, nxr
       DO  j = nys, nyn
          DO  k = nzb+1, nzt

             e_p(k,j,i) = e(k,j,i) + ( dt_3d * ( tsc(2) * tend(k,j,i) +        &
                                                 tsc(3) * te_m(k,j,i) )        &
                                     )
             IF ( e_p(k,j,i) < 0.0_wp )  e_p(k,j,i) = 0.1_wp * e(k,j,i)

          ENDDO
       ENDDO
    ENDDO
    !$acc end parallel

!
!-- Calculate tendencies for the next Runge-Kutta step
    !$acc parallel present( te_m, tend )
    IF ( intermediate_timestep_count == 1 )  THEN
       !$acc loop collapse(3)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                te_m(k,j,i) = tend(k,j,i)
             ENDDO
          ENDDO
       ENDDO
    ELSEIF ( intermediate_timestep_count < &
             intermediate_timestep_count_max )  THEN
       !$acc loop collapse(3)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                te_m(k,j,i) =   -9.5625_wp * tend(k,j,i)                 &
                               + 5.3125_wp * te_m(k,j,i)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel
    !$acc update self(te_m)
    !$acc end data

    CALL cpu_log( log_point(16), 'tke-equation', 'stop' )


 END SUBROUTINE tcm_prognostic


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Computation of the turbulent diffusion coefficients for momentum and heat
!> according to Prandtl-Kolmogorov.
!> @todo consider non-default surfaces
!------------------------------------------------------------------------------!
 SUBROUTINE tcm_diffusivities


    USE arrays_3d,                                                             &
        ONLY:  dd2zu, prho

    USE control_parameters,                                                    &
        ONLY:  e_min, g, prho_reference, wall_adjustment, wall_adjustment_factor

    USE grid_variables,                                                        &
        ONLY:  dx, dy

    IMPLICIT NONE

    INTEGER(iwp) ::  i                   !< loop index
    INTEGER(iwp) ::  j                   !< loop index
    INTEGER(iwp) ::  k                   !< loop index
    INTEGER(iwp) ::  omp_get_thread_num  !< opemmp function to get thread number
    INTEGER(iwp) ::  tn                  !< thread number

    REAL(wp)     ::  dvar_dz             !< vertical gradient of var
    REAL(wp)     ::  l                   !< mixing length
    REAL(wp)     ::  ll                  !< adjusted mixing length
    REAL(wp)     ::  l_stable            !< mixing length according to stratification

!
!-- Default thread number in case of one thread
    tn = 0

!
!-- Compute the turbulent diffusion coefficient for momentum
    !$OMP PARALLEL PRIVATE (i,j,k,l,ll,tn)
!$  tn = omp_get_thread_num()

!
!-- Introduce an optional minimum tke
    IF ( e_min > 0.0_wp )  THEN
       !$OMP DO
       !$acc parallel present(e)
       !$acc loop collapse(3)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             DO  k = nzb+1, nzt
                e(k,j,i) = MAX( e(k,j,i), e_min )
            ENDDO
          ENDDO
       ENDDO
       !$acc end parallel
    ENDIF

    !$OMP DO
    !$acc data present( kh, km, e, prho ) &
    !$acc present( dd2zu, l_grid ) &
    !$acc present( l_wall)

    !$acc parallel
    !$acc loop collapse(2)
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          !$acc loop
          DO  k = nzb+1, nzt
!
!-- Determine the mixing length for LES closure
! inline subroutine mixing_length_les()

             dvar_dz = -1.0_wp * ( prho(k+1,j,i) - prho(k-1,j,i) ) * dd2zu(k)
             IF ( dvar_dz > 0.0_wp ) THEN
                l_stable = 0.76_wp * SQRT( e(k,j,i) )                                &
                                   / SQRT( g / prho_reference * dvar_dz ) + 1E-5_wp
             ELSE
                l_stable = l_grid(k)
             ENDIF
!
!-- Adjustment of the mixing length
             IF ( wall_adjustment )  THEN
                l  = MIN( wall_adjustment_factor * l_wall(k,j,i), l_grid(k), l_stable )
                ll = MIN( wall_adjustment_factor * l_wall(k,j,i), l_grid(k) )
             ELSE
                l  = MIN( l_grid(k), l_stable )
                ll = l_grid(k)
             ENDIF
!-- end of inline subroutine mixing_length_les()
!
!--          Compute diffusion coefficients for momentum and heat
             km(k,j,i) = c_0 * l * SQRT( e(k,j,i) )
             kh(k,j,i) = ( 1.0_wp + 2.0_wp * l / ll ) * km(k,j,i)

          ENDDO
       ENDDO
    ENDDO
    !$acc end parallel
    !$acc end data

!$OMP END PARALLEL

!
!-- Set vertical boundary values (Neumann conditions both at top and bottom)
!
    !$OMP PARALLEL DO
    !$acc parallel present(km, kh)
    !$acc loop collapse(2)
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          !-- model bottom
          km(nzb,j,i) = km(nzb+1,j,i)
          kh(nzb,j,i) = kh(nzb+1,j,i)
          !-- model top
          km(nzt+1,j,i) = km(nzt,j,i)
          kh(nzt+1,j,i) = kh(nzt,j,i)
       ENDDO
    ENDDO
    !$acc end parallel

 END SUBROUTINE tcm_diffusivities


 END MODULE turbulence_closure_mod
