!> @file init_3d_model.f90
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
! Description:
! ------------
!> Allocation of arrays and initialization of the 3D model via
!> a) pre-run the 1D model
!> or
!> b) pre-set constant linear profiles
!> or
!> c) read values of a previous run
!------------------------------------------------------------------------------!
 MODULE configure_3D_MODEL

    USE advec_ws

    USE arrays_3d

    USE constants,                                                             &
        ONLY:  pi
    
    USE control_parameters
    
    USE grid_variables,                                                        &
        ONLY:  dx, dy, ddx2_mg, ddy2_mg

    USE indices

    USE kinds

    USE netcdf_interface,                                                      &
        ONLY:  dots_max, dots_num, dots_unit, dots_label

    USE netcdf_data_input_mod,                                                 &
        ONLY:  init_3d, netcdf_data_input_interpolate, netcdf_data_input_init_3d
    
    USE pegrid
    
    USE random_function_mod 
    
    USE random_generator_parallel,                                             &
        ONLY:  init_parallel_random_generator

    USE read_restart_data_mod,                                                 &
        ONLY:  rrd_read_parts_of_global, rrd_local                                      
    
    USE statistics,                                                            &
        ONLY:  hom, hom_sum, mean_surface_level_height, pr_palm, rmask,        &
               statistic_regions, sums, sums_divnew_l, sums_divold_l, sums_l,  &
               sums_l_l, sums_wsts_bc_l, ts_value,                             &
               weight_pres, weight_substep

    USE stokes_drift_mod,                                                      &   
            ONLY:  init_stokes_drift

    USE surface_layer_fluxes_mod,                                              &
        ONLY:  init_surface_layer_fluxes

    USE surface_mod,                                                           &
        ONLY :  init_surface_arrays, init_surfaces, surf_def_h,     &
                get_topography_top_index_ji, vertical_surfaces_exist, bc_h
   
    USE transpose_indices

    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_init_arrays, tcm_init

    IMPLICIT NONE

    INTEGER(iwp) ::  i             !<
    INTEGER(iwp) ::  ind_array(1)  !<
    INTEGER(iwp) ::  j             !<
    INTEGER(iwp) ::  k             !<
    INTEGER(iwp) ::  k_surf        !< surface level index
    INTEGER(iwp) ::  m             !< index of surface element in surface data type
    INTEGER(iwp) ::  sr            !< index of statistic region

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE   ::  ngp_2dh_l  !<

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  ngp_2dh_outer_l    !<
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  ngp_2dh_s_inner_l  !<

    REAL(wp)     ::  t_surface !< air temperature at the surface

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  p_hydrostatic !< hydrostatic pressure

    INTEGER(iwp) ::  l       !< loop variable
    INTEGER(iwp) ::  nzt_l   !< index of top PE boundary for multigrid level
    REAL(wp) ::  dx_l !< grid spacing along x on different multigrid level
    REAL(wp) ::  dy_l !< grid spacing along y on different multigrid level

    REAL(wp), DIMENSION(1:3) ::  volume_flow_area_l     !<
    REAL(wp), DIMENSION(1:3) ::  volume_flow_initial_l  !<

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  mean_surface_level_height_l    !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ngp_3d_inner_l    !<
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ngp_3d_inner_tmp  !<

    INTEGER(iwp) ::  nz_u_shift   !< 
    INTEGER(iwp) ::  nz_v_shift   !<
    INTEGER(iwp) ::  nz_w_shift   !<
    INTEGER(iwp) ::  nz_s_shift   !<
    INTEGER(iwp) ::  nz_u_shift_l !<
    INTEGER(iwp) ::  nz_v_shift_l !<
    INTEGER(iwp) ::  nz_w_shift_l !<
    INTEGER(iwp) ::  nz_s_shift_l !< 

    public allocate_3d_arrays, init_3d_model, deallocate_3d_variables
    contains

   subroutine deallocate_3d_variables

   Deallocate(mean_surface_level_height_l, ngp_2dh_l, ngp_3d_inner_l)
       DEALLOCATE(ngp_2dh_outer_l, ngp_2dh_s_inner_l)
       DEALLOCATE( mean_surface_level_height,         &
              ngp_2dh, ngp_3d, ngp_3d_inner,          &
              ngp_3d_inner_tmp, sums_divold_l, sums_divnew_l) 
    DEALLOCATE( dp_smooth_factor, rdf, rdf_sc )

    DEALLOCATE( ngp_2dh_outer, ngp_2dh_s_inner,                 &
              rmask,sums, sums_wsts_bc_l,   &
              ts_value )
    DEALLOCATE( ptdf_x, ptdf_y, weight_pres, weight_substep )
    DEALLOCATE( d )
    deallocate(p)
    deallocate(tend)
!    deallocate(sums_l)
!    deallocate(sums_l_l )
#if defined( __nopointer )
    DEALLOCATE( pt, pt_p, u, u_p, v, v_p, w, w_p, tpt_m, tu_m, tv_m, tw_m)
#else
    DEALLOCATE( pt_1, pt_2, pt_3, u_1, u_2, u_3, v_1,v_2, v_3, w_1, w_2, w_3)
#endif
!
!-- Array for storing constant coeffficients of the tridiagonal solver
    IF ( psolver == 'poisfft' )  THEN
       DEALLOCATE( tri, tric)
    ENDIF
#if defined( __nopointer )
       DEALLOCATE( prho, rho_ocean, alpha_T, beta_S, solar3d,                  &
                 sa, sa_p, tsa_m )
#else
       DEALLOCATE( prho_1,rho_1,alpha_T_1, beta_S_1, solar3d_1, sa_1, sa_2, sa_3 )
#endif
!
!-- Allocation of anelastic and Boussinesq approximation specific arrays
    DEALLOCATE( p_hydrostatic )
    DEALLOCATE( rho_air )
    DEALLOCATE( rho_air_zw )
    DEALLOCATE( drho_air )
    DEALLOCATE( drho_air_zw )
    deallocate(uLSforcing,tLSforcing,sLSforcing,vLSforcing)
    DEALLOCATE(uProfileInit, vProfileInit, tProfileInit, sProfileInit)
!
!-- Allocation of flux conversion arrays
    DEALLOCATE( heatflux_input_conversion )
    DEALLOCATE( waterflux_input_conversion )
    DEALLOCATE( momentumflux_input_conversion )
    DEALLOCATE( heatflux_output_conversion )
    DEALLOCATE( waterflux_output_conversion )
    DEALLOCATE( momentumflux_output_conversion )

   end subroutine deallocate_3d_variables

   subroutine allocate_3d_arrays(nCells)

   integer(iwp) :: nCells

!
!-- Allocate arrays
    ALLOCATE( mean_surface_level_height(0:statistic_regions),                  &
              mean_surface_level_height_l(0:statistic_regions),                &
              ngp_2dh(0:statistic_regions), ngp_2dh_l(0:statistic_regions),    &
              ngp_3d(0:statistic_regions),                                     &
              ngp_3d_inner(0:statistic_regions),                               &
              ngp_3d_inner_l(0:statistic_regions),                             &
              ngp_3d_inner_tmp(0:statistic_regions),                           &
              sums_divnew_l(0:statistic_regions),                              &
              sums_divold_l(0:statistic_regions) )
    ALLOCATE( dp_smooth_factor(nzb:nzt), rdf(nzb+1:nzt), rdf_sc(nzb+1:nzt) )

    ALLOCATE( uLSforcing(nzb:nzt), vLSforcing(nzb:nzt), tLSforcing(nzb:nzt) )
    ALLOCATE( sLSforcing(nzb:nzt) )

    ALLOCATE( uProfileInit(nzb:nzt), vProfileInit(nzb:nzt) )
    ALLOCATE( tProfileInit(nzb:nzt), sProfileInit(nzb:nzt) )

    ALLOCATE( ngp_2dh_outer(nzb:nzt+1,0:statistic_regions),                    &
              ngp_2dh_outer_l(nzb:nzt+1,0:statistic_regions),                  &
              ngp_2dh_s_inner(nzb:nzt+1,0:statistic_regions),                  &
              ngp_2dh_s_inner_l(nzb:nzt+1,0:statistic_regions),                &
              rmask(nysg:nyng,nxlg:nxrg,0:statistic_regions),                  &
              sums(nzb:nzt+1,14),                             &
!              sums_l(nzb:nzt+1,pr_palm+max_pr_user,0:threads_per_task-1),      &
              sums_l(nzb:nzt+1,14,0:threads_per_task-1),      &
              sums_l_l(nzb:nzt+1,0:statistic_regions,0:threads_per_task-1),    &
              sums_wsts_bc_l(nzb:nzt+1,0:statistic_regions),                   &
              ts_value(dots_max,0:statistic_regions) )
    ALLOCATE( ptdf_x(nxlg:nxrg), ptdf_y(nysg:nyng) )

    ALLOCATE( d(nzb+1:nzt,nys:nyn,nxl:nxr),                                    &
              p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                &
              tend(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

    ALLOCATE( u_restart(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nCells),                 &
              v_restart(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nCells),                 &
              w_restart(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nCells),                 &
              pt_restart(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nCells),                 &
              sa_restart(nzb:nzt+1,nysg:nyng,nxlg:nxrg,nCells))

#if defined( __nopointer )
    ALLOCATE( pt(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                               &
              pt_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              u(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                &
              u_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                &
              v_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                                &
              w_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              tpt_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                            &
              tu_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              tv_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              tw_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
#else
    ALLOCATE( pt_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              pt_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                             &
              u_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              u_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              u_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              v_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                              &
              w_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    IF (  .NOT.  neutral )  THEN
       ALLOCATE( pt_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ENDIF
#endif

!
!-- Array for storing constant coeffficients of the tridiagonal solver
    IF ( psolver == 'poisfft' )  THEN
       ALLOCATE( tri(nxl_z:nxr_z,nys_z:nyn_z,0:nz-1,2) )
       ALLOCATE( tric(nxl_z:nxr_z,nys_z:nyn_z,0:nz-1) )
    ENDIF

#if defined( __nopointer )
       ALLOCATE( prho(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                          &
                 rho_ocean(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                     &
                 alpha_T(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                 beta_S(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                        &
                 solar3d(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                 sa(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                            &
                 sa_p(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                          &
                 tsa_m(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
#else
       ALLOCATE( prho_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                        &
                 rho_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                         &
                 alpha_T_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                 beta_S_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                        &
                 solar3d_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                       &
                 sa_1(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                          &
                 sa_2(nzb:nzt+1,nysg:nyng,nxlg:nxrg),                          &
                 sa_3(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       prho => prho_1
       rho_ocean  => rho_1  ! routines calc_mean_profile and diffusion_e require
                      ! density to be apointer
       alpha_T => alpha_T_1
       beta_S => beta_S_1
       solar3d => solar3d_1
#endif

!
!-- Allocation of anelastic and Boussinesq approximation specific arrays
    ALLOCATE( p_hydrostatic(nzb:nzt+1) )
    ALLOCATE( rho_air(nzb:nzt+1) )
    ALLOCATE( rho_air_zw(nzb:nzt+1) )
    ALLOCATE( drho_air(nzb:nzt+1) )
    ALLOCATE( drho_air_zw(nzb:nzt+1) )

   rho_air(:) = 1.0_wp
   rho_air_zw(:) = 1.0_wp
!
!-- compute the inverse density array in order to avoid expencive divisions
    drho_air    = 1.0_wp / rho_air
    drho_air_zw = 1.0_wp / rho_air_zw

!
!-- Allocation of flux conversion arrays
    ALLOCATE( heatflux_input_conversion(nzb:nzt+1) )
    ALLOCATE( waterflux_input_conversion(nzb:nzt+1) )
    ALLOCATE( momentumflux_input_conversion(nzb:nzt+1) )
    ALLOCATE( heatflux_output_conversion(nzb:nzt+1) )
    ALLOCATE( waterflux_output_conversion(nzb:nzt+1) )
    ALLOCATE( momentumflux_output_conversion(nzb:nzt+1) )

!
!-- calculate flux conversion factors according to approximation and in-/output mode
    DO  k = nzb, nzt+1

        IF ( TRIM( flux_input_mode ) == 'kinematic' )  THEN
            heatflux_input_conversion(k)      = rho_air_zw(k)
            waterflux_input_conversion(k)     = rho_air_zw(k)
            momentumflux_input_conversion(k)  = rho_air_zw(k)
        ENDIF

        IF ( TRIM( flux_output_mode ) == 'kinematic' )  THEN
            heatflux_output_conversion(k)     = drho_air_zw(k)
            waterflux_output_conversion(k)    = drho_air_zw(k)
            momentumflux_output_conversion(k) = drho_air_zw(k)
        ENDIF

        IF ( .NOT. humidity ) THEN
            waterflux_input_conversion(k)  = 1.0_wp
            waterflux_output_conversion(k) = 1.0_wp
        ENDIF

    ENDDO
!
!-- 1D-array for large scale subsidence velocity
    IF ( .NOT. ALLOCATED( w_subs ) )  THEN
       ALLOCATE ( w_subs(nzb:nzt+1) )
       w_subs = 0.0_wp
    ENDIF

#if ! defined( __nopointer )
!
!-- Initial assignment of the pointers
    IF ( .NOT. neutral )  THEN
       pt => pt_1;  pt_p => pt_2;  tpt_m => pt_3
    ELSE
       pt => pt_1;  pt_p => pt_1;  tpt_m => pt_3
    ENDIF
    u  => u_1;   u_p  => u_2;   tu_m  => u_3
    v  => v_1;   v_p  => v_2;   tv_m  => v_3
    w  => w_1;   w_p  => w_2;   tw_m  => w_3
    sa => sa_1;  sa_p => sa_2;  tsa_m => sa_3
       prho => prho_1
       rho_ocean  => rho_1  ! routines calc_mean_profile and diffusion_e require
                      ! density to be apointer
       alpha_T => alpha_T_1
       beta_S => beta_S_1
       solar3d => solar3d_1

#endif
!
!-- Initialize arrays for turbulence closure
    CALL tcm_init_arrays(nCells)
!
!-- Initialize surface arrays
    CALL init_surface_arrays
!
!-- Allocate arrays containing the RK coefficient for calculation of 
!-- perturbation pressure and turbulent fluxes. At this point values are
!-- set for pressure calculation during initialization (where no timestep
!-- is done). Further below the values needed within the timestep scheme
!-- will be set.
    ALLOCATE( weight_substep(1:intermediate_timestep_count_max),               &
              weight_pres(1:intermediate_timestep_count_max) )
    weight_substep = 1.0_wp
    weight_pres    = 1.0_wp
    intermediate_timestep_count = 0  ! needed when simulated_time = 0.0
       
!    CALL location_message( 'finished', .TRUE. )

!
!-- Initialize time series
    ts_value = 0.0_wp

!
!-- Initialize local summation arrays for routine flow_statistics.
!-- This is necessary because they may not yet have been initialized when they
!-- are called from flow_statistics (or - depending on the chosen model run -
!-- are never initialized)
    sums_divnew_l      = 0.0_wp
    sums_divold_l      = 0.0_wp
    sums_l_l           = 0.0_wp
    sums_wsts_bc_l     = 0.0_wp

   IF ( ws_scheme_sca .OR. ws_scheme_mom )  CALL ws_init        
 !
!--    Initialize the random number generators (from numerical recipes)
       CALL random_function_ini

       IF ( random_generator == 'random-parallel' )  THEN
          CALL init_parallel_random_generator(nx, ny, nys, nyn, nxl, nxr)
       ENDIF  

       CALL init_stokes_drift
end subroutine allocate_3d_arrays

subroutine init_3d_model

!#if ! defined( __nopointer )
!!
!-- Initial assignment of the pointers
!    IF ( .NOT. neutral )  THEN
!       pt => pt_1;  pt_p => pt_2;  tpt_m => pt_3
!    ELSE
!       pt => pt_1;  pt_p => pt_1;  tpt_m => pt_3
!    ENDIF
!    u  => u_1;   u_p  => u_2;   tu_m  => u_3
!    v  => v_1;   v_p  => v_2;   tv_m  => v_3
!    w  => w_1;   w_p  => w_2;   tw_m  => w_3
!    sa => sa_1;  sa_p => sa_2;  tsa_m => sa_3
!#endif

!
!-- Initialize model variables
    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' .OR.                    &
         TRIM( initializing_actions ) /= 'SP_run_continue' )  THEN

       IF ( INDEX(initializing_actions, 'set_constant_profiles') /= 0 )    &
       THEN

!          CALL location_message( 'initializing with constant profiles', .FALSE. )
          !
!--       Use constructed initial profiles (velocity constant with height,
!--       temperature profile with constant gradient)
!          DO  i = nxlg, nxrg
!             DO  j = nysg, nyng
!                pt(:,j,i) = pt_init
!                u(:,j,i)  = u_init
!                v(:,j,i)  = v_init
!             ENDDO
!          ENDDO
!
!--       Mask topography
          u = MERGE( u, 0.0_wp, BTEST( wall_flags_0, 1 ) )
          v = MERGE( v, 0.0_wp, BTEST( wall_flags_0, 2 ) )
!
!--       Set initial horizontal velocities at the lowest computational grid 
!--       levels to zero in order to avoid too small time steps caused by the 
!--       diffusion limit in the initial phase of a run (at k=1, dz/2 occurs
!--       in the limiting formula!). 
!--       Please note, in case land- or urban-surface model is used and a 
!--       spinup is applied, masking the lowest computational level is not 
!--       possible as MOST as well as energy-balance parametrizations will not 
!--       work with zero wind velocity. 
          IF ( ibc_uv_b /= 1  .AND.  .NOT.  spinup )  THEN
             DO  i = nxlg, nxrg
                DO  j = nysg, nyng
                   DO  k = nzb, nzt
                      u(k,j,i) = MERGE( u(k,j,i), 0.0_wp,                      &
                                        BTEST( wall_flags_0(k,j,i), 20 ) )
                      v(k,j,i) = MERGE( v(k,j,i), 0.0_wp,                      &
                                        BTEST( wall_flags_0(k,j,i), 21 ) )
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
!             DO  i = nxlg, nxrg
!                DO  j = nysg, nyng
!                   sa(:,j,i) = sa_init
!                ENDDO
!             ENDDO
!
!--       Compute initial temperature field and other constants used in case
!--       of a sloping surface
          IF ( sloping_surface )  CALL init_slope
!
!--       Initialize surface variables, e.g. friction velocity, momentum 
!--       fluxes, etc. 
          CALL init_surfaces
       ENDIF

!       CALL location_message( 'initializing statistics, boundary conditions, etc.', &
!                              .FALSE. )

!
!--    Bottom boundary
       IF ( ibc_uv_b == 0 .OR. ibc_uv_b == 2  )  THEN
          u(nzb,:,:) = 0.0_wp
          v(nzb,:,:) = 0.0_wp
       ENDIF

!
!--    Apply channel flow boundary condition
       IF ( TRIM( bc_uv_t ) == 'dirichlet_0' )  THEN
          u(nzt+1,:,:) = 0.0_wp
          v(nzt+1,:,:) = 0.0_wp
       ENDIF

!
!--    Set the reference state to be used in the buoyancy terms (for ocean runs
!--    the reference state will be set (overwritten) in init_ocean)
       IF ( use_single_reference_value )  THEN
             ref_state(:) = pt_reference
       ELSE
             ref_state(:) = pt_init(:)
       ENDIF

!
!--    For the moment, vertical velocity is zero
       w = 0.0_wp

!
!--    Initialize array sums (must be defined in first call of pres)
       sums = 0.0_wp

!--    If required, change the surface temperature at the start of the 3D run
       IF ( pt_surface_initial_change /= 0.0_wp )  THEN
          pt(nzb,:,:) = pt(nzb,:,:) + pt_surface_initial_change
       ENDIF

!
!--    Initialize old and new time levels.
       tpt_m = 0.0_wp; tu_m = 0.0_wp; tv_m = 0.0_wp; tw_m = 0.0_wp
       pt_p = pt; u_p = u; v_p = v; w_p = w
       tsa_m = 0.0_wp
       sa_p  = sa
       
!       CALL location_message( 'finished', .TRUE. )

    ELSEIF ( TRIM( initializing_actions ) == 'read_restart_data'  .OR.         &
             TRIM( initializing_actions ) == 'cyclic_fill' )                   &
    THEN

!       CALL location_message( 'initializing in case of restart / cyclic_fill', &
!                              .FALSE. )
!
!--    Initialize surface elements and its attributes, e.g. heat- and 
!--    momentumfluxes, roughness, scaling parameters. As number of surface 
!--    elements might be different between runs, e.g. in case of cyclic fill, 
!--    and not all surface elements are read, surface elements need to be 
!--    initialized before.     
       CALL init_surfaces
!
!--    Read processor specific binary data from restart file
       DO  i = 0, io_blocks-1
          IF ( i == io_group )  THEN
             CALL rrd_local
          ENDIF
#if defined( __parallel )
          CALL MPI_BARRIER( comm2d, ierr )
#endif
       ENDDO
!
!
!--    Calculate initial temperature field and other constants used in case
!--    of a sloping surface
       IF ( sloping_surface )  CALL init_slope

!
!--    Initialize new time levels (only done in order to set boundary values
!--    including ghost points)
       pt_p = pt; u_p = u; v_p = v; w_p = w
       sa_p = sa

!
!--    Allthough tendency arrays are set in prognostic_equations, they have
!--    have to be predefined here because they are used (but multiplied with 0)
!--    there before they are set. 
       tpt_m = 0.0_wp; tu_m = 0.0_wp; tv_m = 0.0_wp; tw_m = 0.0_wp
       tsa_m = 0.0_wp
!
!       CALL location_message( 'finished', .TRUE. )

    ELSE
!
! do nothing
!--    Actually this part of the programm should not be reached
!!       message_string = 'unknown initializing problem'
!       CALL message( 'init_3d_model', 'PA0193', 1, 2, 0, 6, 0 )
    ENDIF

    CALL tcm_init

!
!-- Before initializing further modules, compute total sum of active mask 
!-- grid points and the mean surface level height for each statistic region.
!-- ngp_2dh: number of grid points of a horizontal cross section through the
!--          total domain
!-- ngp_3d:  number of grid points of the total domain
    ngp_2dh_outer_l   = 0
    ngp_2dh_outer     = 0
    ngp_2dh_s_inner_l = 0
    ngp_2dh_s_inner   = 0
    ngp_2dh_l         = 0
    ngp_2dh           = 0
    ngp_3d_inner_l    = 0.0_wp
    ngp_3d_inner      = 0
    ngp_3d            = 0
    ngp_sums          = ( nz + 2 ) * ( pr_palm + max_pr_user )

    mean_surface_level_height   = 0.0_wp
    mean_surface_level_height_l = 0.0_wp
!
!-- Pre-set masks for regional statistics. Default is the total model domain.
!-- Ghost points are excluded because counting values at the ghost boundaries
!-- would bias the statistics
    rmask = 1.0_wp
    rmask(:,nxlg:nxl-1,:) = 0.0_wp;  rmask(:,nxr+1:nxrg,:) = 0.0_wp
    rmask(nysg:nys-1,:,:) = 0.0_wp;  rmask(nyn+1:nyng,:,:) = 0.0_wp
!
!
!-- To do: New concept for these non-topography grid points!
    DO  sr = 0, statistic_regions
       DO  i = nxl, nxr
          DO  j = nys, nyn
             IF ( rmask(j,i,sr) == 1.0_wp )  THEN
!
!--             All xy-grid points
                ngp_2dh_l(sr) = ngp_2dh_l(sr) + 1
!
!--             Determine mean surface-level height. In case of downward-
!--             facing walls are present, more than one surface level exist.
!--             In this case, use the lowest surface-level height. 
                IF ( surf_def_h(0)%start_index(j,i) <=                         &
                     surf_def_h(0)%end_index(j,i) )  THEN
                   m = surf_def_h(0)%start_index(j,i)
                   k = surf_def_h(0)%k(m)
                   mean_surface_level_height_l(sr) =                           &
                                       mean_surface_level_height_l(sr) + zw(k-1)
                ENDIF
                k_surf = k - 1

                DO  k = nzb, nzt+1
!
!--                xy-grid points above topography
                   ngp_2dh_outer_l(k,sr) = ngp_2dh_outer_l(k,sr)     +         &
                                  MERGE( 1, 0, BTEST( wall_flags_0(k,j,i), 24 ) )

                   ngp_2dh_s_inner_l(k,sr) = ngp_2dh_s_inner_l(k,sr) +         &
                                  MERGE( 1, 0, BTEST( wall_flags_0(k,j,i), 22 ) )

                ENDDO
!
!--             All grid points of the total domain above topography
                ngp_3d_inner_l(sr) = ngp_3d_inner_l(sr) + ( nz - k_surf + 2 )



             ENDIF
          ENDDO
       ENDDO
    ENDDO
!
!-- Initialize arrays encompassing number of grid-points in inner and outer
!-- domains, statistic regions, etc. Mainly used for horizontal averaging
!-- of turbulence statistics. Please note, user_init must be called before
!-- doing this.   
    sr = statistic_regions + 1
#if defined( __parallel )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_l(0), ngp_2dh(0), sr, MPI_INTEGER, MPI_SUM,    &
                        comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_outer_l(0,0), ngp_2dh_outer(0,0), (nz+2)*sr,   &
                        MPI_INTEGER, MPI_SUM, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_2dh_s_inner_l(0,0), ngp_2dh_s_inner(0,0),          &
                        (nz+2)*sr, MPI_INTEGER, MPI_SUM, comm2d, ierr )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( ngp_3d_inner_l(0), ngp_3d_inner_tmp(0), sr, MPI_REAL,  &
                        MPI_SUM, comm2d, ierr )
    ngp_3d_inner = INT( ngp_3d_inner_tmp, KIND = SELECTED_INT_KIND( 18 ) )
    IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
    CALL MPI_ALLREDUCE( mean_surface_level_height_l(0),                        &
                        mean_surface_level_height(0), sr, MPI_REAL,            &
                        MPI_SUM, comm2d, ierr )
    mean_surface_level_height = mean_surface_level_height / REAL( ngp_2dh )
#else
    ngp_2dh         = ngp_2dh_l
    ngp_2dh_outer   = ngp_2dh_outer_l
    ngp_2dh_s_inner = ngp_2dh_s_inner_l
    ngp_3d_inner    = INT( ngp_3d_inner_l, KIND = SELECTED_INT_KIND( 18 ) )
    mean_surface_level_height = mean_surface_level_height_l / REAL( ngp_2dh_l )
#endif
    ngp_3d = INT ( ngp_2dh, KIND = SELECTED_INT_KIND( 18 ) ) * &
             INT ( (nz + 2 ), KIND = SELECTED_INT_KIND( 18 ) )

!
!-- Set a lower limit of 1 in order to avoid zero divisions in flow_statistics,
!-- buoyancy, etc. A zero value will occur for cases where all grid points of
!-- the respective subdomain lie below the surface topography
    ngp_2dh_outer   = MAX( 1, ngp_2dh_outer(:,:)   ) 
    ngp_3d_inner    = MAX( INT(1, KIND = SELECTED_INT_KIND( 18 )),             &
                           ngp_3d_inner(:) )
    ngp_2dh_s_inner = MAX( 1, ngp_2dh_s_inner(:,:) ) 

!-- Initialize quantities for special advections schemes
!    CALL init_advec

!-- Impose random perturbation on the horizontal velocity field and then
!-- remove the divergences from the velocity field at the initial stage
    IF ( create_disturbances  .AND.  dt_disturb > 0.0_wp  .AND. &
         TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN

         !$acc data copy( u, v )
!       CALL location_message( 'creating initial disturbances', .FALSE. )
       CALL disturb_field( 'u', tend, u)
       CALL disturb_field( 'v', tend, v )
 !      CALL disturb_field( 'pt', tend, pt )
!       call disturb_field( 'sa', tend, sa )

!$acc end data

!       CALL location_message( 'finished', .TRUE. )

!       CALL location_message( 'calling pressure solver', .FALSE. )
       n_sor = nsor_ini

       !$acc data copyin( rho_air, rho_air_zw, ddzw, ddzu, wall_flags_0, ngp_2dh_outer, bc_h ) &
       !$acc copy( u, v, w )
       CALL pres
       !$acc end data
       n_sor = nsor
!       CALL location_message( 'finished', .TRUE. )
      
    ENDIF

!--    Initialize quantities needed for the ocean model
       CALL init_ocean
#ifndef __GPU
!-- Initialize surface layer (done after LSM as roughness length are required
!-- for initialization
    IF ( constant_flux_layer )  THEN
       CALL location_message( 'initializing surface layer', .FALSE. )
       CALL init_surface_layer_fluxes
       CALL location_message( 'finished', .TRUE. )
    ENDIF
#endif
    !
!-- Initialize the ws-scheme.    

!
!-- Setting weighting factors for calculation of perturbation pressure 
!-- and turbulent quantities from the RK substeps
    IF ( TRIM(timestep_scheme) == 'runge-kutta-3' )  THEN      ! for RK3-method

       weight_substep(1) = 1._wp/6._wp
       weight_substep(2) = 3._wp/10._wp
       weight_substep(3) = 8._wp/15._wp

       weight_pres(1)    = 1._wp/3._wp
       weight_pres(2)    = 5._wp/12._wp
       weight_pres(3)    = 1._wp/4._wp

    ELSEIF ( TRIM(timestep_scheme) == 'runge-kutta-2' )  THEN  ! for RK2-method

       weight_substep(1) = 1._wp/2._wp
       weight_substep(2) = 1._wp/2._wp
          
       weight_pres(1)    = 1._wp/2._wp
       weight_pres(2)    = 1._wp/2._wp        

    ELSE                                     ! for Euler-method

       weight_substep(1) = 1.0_wp      
       weight_pres(1)    = 1.0_wp                   

    ENDIF

!
!-- Initialize Rayleigh damping factors
    rdf    = 0.0_wp
    rdf_sc = 0.0_wp
    IF ( rayleigh_damping_factor /= 0.0_wp )  THEN
       IF (  .NOT.  ocean )  THEN
          DO  k = nzb+1, nzt
             IF ( zu(k) >= rayleigh_damping_height )  THEN
                rdf(k) = rayleigh_damping_factor *                             &
                      ( SIN( pi * 0.5_wp * ( zu(k) - rayleigh_damping_height ) &
                             / ( zu(nzt) - rayleigh_damping_height ) )         &
                      )**2
             ENDIF
          ENDDO
       ELSE
          DO  k = nzt, nzb+1, -1
             IF ( zu(k) <= rayleigh_damping_height )  THEN
                rdf(k) = rayleigh_damping_factor *                             &
                      ( SIN( pi * 0.5_wp * ( rayleigh_damping_height - zu(k) ) &
                             / ( rayleigh_damping_height - zu(nzb+1) ) )       &
                      )**2
             ENDIF
          ENDDO
       ENDIF
    ENDIF
    IF ( scalar_rayleigh_damping )  rdf_sc = rdf

!
!-- Initialize the starting level and the vertical smoothing factor used for 
!-- the external pressure gradient
    dp_smooth_factor = 1.0_wp
    IF ( dp_external )  THEN
!
!--    Set the starting level dp_level_ind_b only if it has not been set before
!--    (e.g. in init_grid).
       IF ( dp_level_ind_b == 0 )  THEN
          ind_array = MINLOC( ABS( dp_level_b - zu ) )
          dp_level_ind_b = ind_array(1) - 1 + nzb 
                                        ! MINLOC uses lower array bound 1
       ENDIF
       IF ( dp_smooth )  THEN
          dp_smooth_factor(:dp_level_ind_b) = 0.0_wp
          DO  k = dp_level_ind_b+1, nzt
             dp_smooth_factor(k) = 0.5_wp * ( 1.0_wp + SIN( pi *               &
                        ( REAL( k - dp_level_ind_b, KIND=wp ) /                &
                          REAL( nzt - dp_level_ind_b, KIND=wp ) - 0.5_wp ) ) )
          ENDDO
       ENDIF
    ENDIF

!
!-- Initialize damping zone for the potential temperature in case of
!-- non-cyclic lateral boundaries. The damping zone has the maximum value 
!-- at the inflow boundary and decreases to zero at pt_damping_width.
    ptdf_x = 0.0_wp
    ptdf_y = 0.0_wp
    IF ( bc_lr_dirrad )  THEN
       DO  i = nxl, nxr
          IF ( ( i * dx ) < pt_damping_width )  THEN
             ptdf_x(i) = pt_damping_factor * ( SIN( pi * 0.5_wp *              &
                            REAL( pt_damping_width - i * dx, KIND=wp ) / (     &
                            REAL( pt_damping_width, KIND=wp ) ) ) )**2 
          ENDIF
       ENDDO
    ELSEIF ( bc_lr_raddir )  THEN
       DO  i = nxl, nxr
          IF ( ( i * dx ) > ( nx * dx - pt_damping_width ) )  THEN
             ptdf_x(i) = pt_damping_factor *                                   &
                         SIN( pi * 0.5_wp *                                    &
                                 ( ( i - nx ) * dx + pt_damping_width ) /      &
                                 REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO 
    ELSEIF ( bc_ns_dirrad )  THEN
       DO  j = nys, nyn
          IF ( ( j * dy ) > ( ny * dy - pt_damping_width ) )  THEN
             ptdf_y(j) = pt_damping_factor *                                   &
                         SIN( pi * 0.5_wp *                                    &
                                 ( ( j - ny ) * dy + pt_damping_width ) /      &
                                 REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO 
    ELSEIF ( bc_ns_raddir )  THEN
       DO  j = nys, nyn
          IF ( ( j * dy ) < pt_damping_width )  THEN
             ptdf_y(j) = pt_damping_factor *                                   &
                         SIN( pi * 0.5_wp *                                    &
                                ( pt_damping_width - j * dy ) /                &
                                REAL( pt_damping_width, KIND=wp ) )**2
          ENDIF
       ENDDO
    ENDIF
!
!-- Check if maximum number of allowed timeseries is exceeded
    IF ( dots_num > dots_max )  THEN
       WRITE( message_string, * ) 'number of time series quantities exceeds',  &
                                  ' its maximum of dots_max = ', dots_max,     &
                                  '&Please increase dots_max in modules.f90.'
       CALL message( 'init_3d_model', 'PA0194', 1, 2, 0, 6, 0 )    
    ENDIF

!
!-- Input binary data file is not needed anymore. This line must be placed
!-- after call of user_init!
    CALL close_file( 13 )
!
!-- In case of nesting, put an barrier to assure that all parent and child 
!-- domains finished initialization. 

!    CALL location_message( 'leaving init_3d_model', .TRUE. )

 END SUBROUTINE init_3d_model

 END MODULE configure_3D_MODEL
