!> @file check_parameters.f90
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
!> Check control parameters and deduce further quantities.
!------------------------------------------------------------------------------!
 SUBROUTINE check_parameters


    USE arrays_3d
    USE constants
    USE control_parameters
    USE grid_variables
    USE indices
    USE kinds

    USE netcdf_data_input_mod,                                                 &
        ONLY:  init_model, input_pids_static, netcdf_data_input_check_dynamic, &
               netcdf_data_input_check_static

    USE netcdf_interface,                                                      &
        ONLY:  dopr_unit, do3d_unit, netcdf_data_format,            &
               netcdf_data_format_string, dots_unit, heatflux_output_unit,     &
               waterflux_output_unit, momentumflux_output_unit

    USE pegrid

    USE profil_parameter
    USE statistics
    USE stokes_drift_mod,                                                      &
        ONLY: stokes_drift_check_parameters
    USE statistics
    USE transpose_indices
    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_check_data_output, tcm_check_parameters

    IMPLICIT NONE

    CHARACTER (LEN=varnamelength)  ::  var           !< variable name
    CHARACTER (LEN=7)   ::  unit                     !< unit of variable
    CHARACTER (LEN=8)   ::  date                     !< current date string
    CHARACTER (LEN=10)  ::  time                     !< current time string
    CHARACTER (LEN=20)  ::  ensemble_string          !< string containing number of ensemble member
    CHARACTER (LEN=40)  ::  coupling_string          !< string containing type of coupling
    CHARACTER (LEN=100) ::  action                   !< flag string

    INTEGER(iwp) ::  i                               !< loop index
    INTEGER(iwp) ::  ilen                            !< string length
    INTEGER(iwp) ::  j                               !< loop index
    INTEGER(iwp) ::  k                               !< loop index
    INTEGER(iwp) ::  kk                              !< loop index
    INTEGER(iwp) ::  netcdf_data_format_save         !< initial value of netcdf_data_format
    INTEGER(iwp) ::  position                        !< index position of string

    LOGICAL     ::  found                            !< flag, true if output variable is already marked for averaging

    REAL(wp)    ::  dum                              !< dummy variable
    REAL(wp)    ::  gradient                         !< local gradient
    REAL(wp)    ::  remote = 0.0_wp                  !< MPI id of remote processor
    REAL(wp)    ::  time_to_be_simulated_from_reference_point  !< time to be simulated from reference point


!    CALL location_message( 'checking parameters', .FALSE. )

#if defined(__fullPalm)
!
!-- At first, check static and dynamic input for consistency
    CALL netcdf_data_input_check_dynamic
    CALL netcdf_data_input_check_static
#endif
    !
!-- Check for overlap combinations, which are not realized yet
    IF ( transpose_compute_overlap )  THEN
       IF ( numprocs == 1 )  STOP '+++ transpose-compute-overlap not implemented for single PE runs'
    ENDIF
!
!-- User settings for restart times requires that "restart" has been given as
!-- file activation string. Otherwise, binary output would not be saved by
!-- palmrun.
    IF (  ( restart_time /= 9999999.9_wp  .OR.  dt_restart /= 9999999.9_wp )   &
         .AND.  .NOT. write_binary )  THEN
       WRITE( message_string, * ) 'manual restart settings requires file ',    &
                                  'activation string "restart"'
       CALL message( 'check_parameters', 'PA0001', 1, 2, 0, 6, 0 )
    ENDIF


!
!-- Generate the file header which is used as a header for most of PALM's
!-- output files
    CALL DATE_AND_TIME( date, time )
    run_date = date(7:8)//'-'//date(5:6)//'-'//date(3:4)
    run_time = time(1:2)//':'//time(3:4)//':'//time(5:6)
       coupling_string = ''
       ensemble_string = ''

    WRITE ( run_description_header,                                            &
            '(A,2X,A,2X,A,A,A,I2.2,A,A,A,2X,A,A,2X,A,1X,A)' )                  &
          TRIM( version ), TRIM( revision ), 'run: ',                          &
          TRIM( run_identifier ), '.', runnr, TRIM( coupling_string ),         &
          run_date, run_time

!
!-- Check the general loop optimization method
    SELECT CASE ( TRIM( loop_optimization ) )

       CASE ( 'cache', 'vector' )
          CONTINUE

       CASE DEFAULT
          message_string = 'illegal value given for loop_optimization: "' //   &
                           TRIM( loop_optimization ) // '"'
          CALL message( 'check_parameters', 'PA0013', 1, 2, 0, 6, 0 )

    END SELECT
!
!-- Check turbulence closure setup
    CALL tcm_check_parameters
!
!-- Check flux input mode
    IF ( TRIM( flux_input_mode ) /= 'dynamic'    .AND.                         &
         TRIM( flux_input_mode ) /= 'kinematic'  .AND.                         &
         TRIM( flux_input_mode ) /= 'approximation-specific' )  THEN
       message_string = 'unknown flux input mode: flux_input_mode = "' //      &
                        TRIM( flux_input_mode ) // '"'
       CALL message( 'check_parameters', 'PA0450', 1, 2, 0, 6, 0 )
    ENDIF
!-- Set flux input mode according to approximation if applicable
    IF ( TRIM( flux_input_mode ) == 'approximation-specific' )  THEN
       IF ( TRIM( approximation ) == 'anelastic' )  THEN
          flux_input_mode = 'dynamic'
       ELSEIF ( TRIM( approximation ) == 'boussinesq' )  THEN
          flux_input_mode = 'kinematic'
       ENDIF
    ENDIF

!
!-- Check flux output mode
    IF ( TRIM( flux_output_mode ) /= 'dynamic'    .AND.                        &
         TRIM( flux_output_mode ) /= 'kinematic'  .AND.                        &
         TRIM( flux_output_mode ) /= 'approximation-specific' )  THEN
       message_string = 'unknown flux output mode: flux_output_mode = "' //    &
                        TRIM( flux_output_mode ) // '"'
       CALL message( 'check_parameters', 'PA0451', 1, 2, 0, 6, 0 )
    ENDIF
!-- Set flux output mode according to approximation if applicable
    IF ( TRIM( flux_output_mode ) == 'approximation-specific' )  THEN
       IF ( TRIM( approximation ) == 'anelastic' )  THEN
          flux_output_mode = 'dynamic'
       ELSEIF ( TRIM( approximation ) == 'boussinesq' )  THEN
          flux_output_mode = 'kinematic'
       ENDIF
    ENDIF

!
!-- set the flux output units according to flux_output_mode
    IF ( TRIM( flux_output_mode ) == 'kinematic' ) THEN
        heatflux_output_unit              = 'K m/s'
        waterflux_output_unit             = 'kg/kg m/s'
        momentumflux_output_unit          = 'm2/s2'
    ELSEIF ( TRIM( flux_output_mode ) == 'dynamic' ) THEN
        heatflux_output_unit              = 'W/m2'
        waterflux_output_unit             = 'W/m2'
        momentumflux_output_unit          = 'N/m2'
    ENDIF

!-- set time series output units for fluxes
    dots_unit(14:16) = heatflux_output_unit
    dots_unit(21)    = waterflux_output_unit
    dots_unit(19:20) = momentumflux_output_unit
    IF( momentum_advec == 'ws-scheme' .AND.                                    &
        .NOT. call_psolver_at_all_substeps  ) THEN
        message_string = 'psolver must be called at each RK3 substep when "'// &
                      TRIM(momentum_advec) // ' "is used for momentum_advec'
        CALL message( 'check_parameters', 'PA0344', 1, 2, 0, 6, 0 )
    END IF
!
!-- Advection schemes:
    IF ( momentum_advec /= 'pw-scheme'  .AND.                                  &
         momentum_advec /= 'ws-scheme'  .AND.                                  &
         momentum_advec /= 'up-scheme' )                                       &
    THEN
       message_string = 'unknown advection scheme: momentum_advec = "' //      &
                        TRIM( momentum_advec ) // '"'
       CALL message( 'check_parameters', 'PA0022', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ( momentum_advec == 'ws-scheme' .OR.  scalar_advec == 'ws-scheme' )   &
           .AND. ( timestep_scheme == 'euler' .OR.                             &
                   timestep_scheme == 'runge-kutta-2' ) )                      &
    THEN
       message_string = 'momentum_advec or scalar_advec = "'                   &
         // TRIM( momentum_advec ) // '" is not allowed with ' //              &
         'timestep_scheme = "' // TRIM( timestep_scheme ) // '"'
       CALL message( 'check_parameters', 'PA0023', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( scalar_advec /= 'pw-scheme'  .AND.  scalar_advec /= 'ws-scheme' .AND. &
         scalar_advec /= 'bc-scheme' .AND. scalar_advec /= 'up-scheme' )       &
    THEN
       message_string = 'unknown advection scheme: scalar_advec = "' //        &
                        TRIM( scalar_advec ) // '"'
       CALL message( 'check_parameters', 'PA0024', 1, 2, 0, 6, 0 )
    ENDIF

!-- Set LOGICAL switches to enhance performance
    IF ( momentum_advec == 'ws-scheme' )  ws_scheme_mom = .TRUE.
    IF ( scalar_advec   == 'ws-scheme' )  ws_scheme_sca = .TRUE.


!
!-- Timestep schemes:
    SELECT CASE ( TRIM( timestep_scheme ) )

       CASE ( 'euler' )
          intermediate_timestep_count_max = 1

       CASE ( 'runge-kutta-2' )
          intermediate_timestep_count_max = 2

       CASE ( 'runge-kutta-3' )
          intermediate_timestep_count_max = 3

       CASE DEFAULT
          message_string = 'unknown timestep scheme: timestep_scheme = "' //   &
                           TRIM( timestep_scheme ) // '"'
          CALL message( 'check_parameters', 'PA0027', 1, 2, 0, 6, 0 )

    END SELECT

    IF ( (momentum_advec /= 'pw-scheme' .AND. momentum_advec /= 'ws-scheme')   &
         .AND. timestep_scheme(1:5) == 'runge' ) THEN
       message_string = 'momentum advection scheme "' // &
                        TRIM( momentum_advec ) // '" & does not work with ' // &
                        'timestep_scheme "' // TRIM( timestep_scheme ) // '"'
       CALL message( 'check_parameters', 'PA0029', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Initializing actions must have been set by the user

    WRITE(message_string,*) 'initializing_actions = ',initializing_actions
!    CALL location_message(adjustl(trim(message_string)),.TRUE.)

    IF ( TRIM( initializing_actions ) == '' )  THEN
       message_string = 'no value specified for initializing_actions'
       CALL message( 'check_parameters', 'PA0149', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN
!
!--    No restart run: several initialising actions are possible
       action = initializing_actions
       DO  WHILE ( TRIM( action ) /= '' )
          position = INDEX( action, ' ' )
          SELECT CASE ( action(1:position-1) )

             CASE ( 'set_constant_profiles', 'set_1d-model_profiles',          &
                    'by_user', 'initialize_vortex', 			       &
                    'initialize_2D_bubble', 'initialize_3D_bubble', 'inifor' )
                action = action(position+1:)

             CASE DEFAULT
                message_string = 'initializing_action = "' //                  &
                                 TRIM( action ) // '" unknown or not allowed'
                CALL message( 'check_parameters', 'PA0030', 1, 2, 0, 6, 0 )

          END SELECT
       ENDDO
    ENDIF

    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.    &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //    &
                        ' and "set_1d-model_profiles" are not allowed ' //     &
                        'simultaneously'
       CALL message( 'check_parameters', 'PA0031', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.    &
         INDEX( initializing_actions, 'by_user' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //    &
                        ' and "by_user" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0032', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'by_user' ) /= 0  .AND.                  &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "by_user" and ' //             &
                        '"set_1d-model_profiles" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0033', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check for Stokes drift settings if Stokes forcing in ocean mode is used
    IF ( ocean .AND. stokes_force ) CALL stokes_drift_check_parameters
    !
!-- In case of no model continuation run, check initialising parameters and
!-- deduce further quantities
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN


      !
!--    Initial profiles for 1D and 3D model, respectively (u,v further below)
       pt_init = pt_surface
       sa_init = sa_surface
       !--
!--    If required, compute initial profile of the geostrophic wind
!--    (component ug)
       i = 1
       gradient = 0.0_wp
          ug_vertical_gradient_level_ind(1) = nzt+1
          ug(nzt+1) = ug_surface
          DO  k = nzt, nzb, -1
             IF ( i < 11 )  THEN
                IF ( ug_vertical_gradient_level(i) > zu(k)  .AND.              &
                     ug_vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = ug_vertical_gradient(i) / 100.0_wp
                   ug_vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   ug(k) = ug(k+1) - dzu(k+1) * gradient
                ELSE
                   ug(k)   = ug_surface - 0.5_wp * dzu(k+1) * gradient
                   ug(k+1) = ug_surface + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                ug(k) = ug(k+1)
             ENDIF
          ENDDO

!
!--    In case of no given gradients for ug, choose a zero gradient
       IF ( ug_vertical_gradient_level(1) == -9999999.9_wp )  THEN
          ug_vertical_gradient_level(1) = 0.0_wp
       ENDIF

!
!--
!--    If required, compute initial profile of the geostrophic wind
!--    (component vg)
       i = 1
       gradient = 0.0_wp
          vg_vertical_gradient_level_ind(1) = nzt+1
          vg(nzt+1) = vg_surface
          DO  k = nzt, nzb, -1
             IF ( i < 11 )  THEN
                IF ( vg_vertical_gradient_level(i) > zu(k)  .AND.              &
                     vg_vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = vg_vertical_gradient(i) / 100.0_wp
                   vg_vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   vg(k) = vg(k+1) - dzu(k+1) * gradient
                ELSE
                   vg(k)   = vg_surface - 0.5_wp * dzu(k+1) * gradient
                   vg(k+1) = vg_surface + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                vg(k) = vg(k+1)
             ENDIF
          ENDDO
!
!--    In case of no given gradients for vg, choose a zero gradient
       IF ( vg_vertical_gradient_level(1) == -9999999.9_wp )  THEN
          vg_vertical_gradient_level(1) = 0.0_wp
       ENDIF

!
!--    Let the initial wind profiles be the calculated ug/vg profiles or
!--    interpolate them from wind profile data (if given)
       IF ( u_profile(1) == 9999999.9_wp  .AND.  v_profile(1) == 9999999.9_wp )  THEN

          u_init = ug
          v_init = vg

       ELSEIF ( u_profile(1) == 0.0_wp  .AND.  v_profile(1) == 0.0_wp )  THEN

          IF ( uv_heights(1) /= 0.0_wp )  THEN
             message_string = 'uv_heights(1) must be 0.0'
             CALL message( 'check_parameters', 'PA0345', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( omega /= 0.0_wp )  THEN
             message_string = 'Coriolis force must be switched off (by setting omega=0.0)' //  &
                              ' when prescribing the forcing by u_profile and v_profile'
             CALL message( 'check_parameters', 'PA0347', 1, 2, 0, 6, 0 )
          ENDIF

          use_prescribed_profile_data = .TRUE.

          kk = 1
          u_init(0) = 0.0_wp
          v_init(0) = 0.0_wp

          DO  k = 1, nz+1

             IF ( kk < 100 )  THEN
                DO  WHILE ( uv_heights(kk+1) <= zu(k) )
                   kk = kk + 1
                   IF ( kk == 100 )  EXIT
                ENDDO
             ENDIF

             IF ( kk < 100  .AND.  uv_heights(kk+1) /= 9999999.9_wp )  THEN
                u_init(k) = u_profile(kk) + ( zu(k) - uv_heights(kk) ) /       &
                                       ( uv_heights(kk+1) - uv_heights(kk) ) * &
                                       ( u_profile(kk+1) - u_profile(kk) )
                v_init(k) = v_profile(kk) + ( zu(k) - uv_heights(kk) ) /       &
                                       ( uv_heights(kk+1) - uv_heights(kk) ) * &
                                       ( v_profile(kk+1) - v_profile(kk) )
             ELSE
                u_init(k) = u_profile(kk)
                v_init(k) = v_profile(kk)
             ENDIF

          ENDDO

       ELSE

          message_string = 'u_profile(1) and v_profile(1) must be 0.0'
          CALL message( 'check_parameters', 'PA0346', 1, 2, 0, 6, 0 )

       ENDIF

!
!--    Compute initial temperature profile using the given temperature gradients
       IF (  .NOT.  neutral )  THEN
          CALL init_vertical_profiles( pt_vertical_gradient_level_ind,          &
                                       pt_vertical_gradient_level,              &
                                       pt_vertical_gradient, pt_init,           &
                                       pt_surface, bc_pt_t_val )
       ENDIF

!
!--    If required, compute initial salinity profile using the given salinity
!--    gradients
          CALL init_vertical_profiles( sa_vertical_gradient_level_ind,          &
                                       sa_vertical_gradient_level,              &
                                       sa_vertical_gradient, sa_init,           &
                                      sa_surface, dum )

    ENDIF
      IF ( subs_vertical_gradient_level(1) /= -9999999.9_wp )  THEN
           message_string = 'Enable usage of large scale subsidence by ' //    &
                            'setting large_scale_subsidence = .T..'
          CALL message( 'check_parameters', 'PA0381', 1, 2, 0, 6, 0 )
        ENDIF

!
!-- Overwrite latitude if necessary and compute Coriolis parameter.
!-- To do - move initialization of f and fs to coriolis_mod.
    IF ( input_pids_static )  THEN
       latitude  = init_model%latitude
       longitude = init_model%longitude
    ENDIF

    f  = 2.0_wp * omega * SIN( latitude / 180.0_wp * pi )
    fs = 0.0_wp * omega * COS( latitude / 180.0_wp * pi )

!
!-- Check and set buoyancy related parameters and switches
    IF ( reference_state == 'horizontal_average' )  THEN
       CONTINUE
    ELSEIF ( reference_state == 'initial_profile' )  THEN
       use_initial_profile_as_reference = .TRUE.
    ELSEIF ( reference_state == 'single_value' )  THEN
       use_single_reference_value = .TRUE.
       IF ( pt_reference == 9999999.9_wp )  pt_reference = pt_surface
       vpt_reference = pt_reference * ( 1.0_wp + 0.61_wp * q_surface )
    ELSE
       message_string = 'illegal value for reference_state: "' //              &
                        TRIM( reference_state ) // '"'
       CALL message( 'check_parameters', 'PA0056', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Sign of buoyancy/stability terms
    atmos_ocean_sign = -1.0_wp

!
!-- Ocean version must use flux boundary conditions at the top
    IF ( ocean .AND. .NOT. use_top_fluxes )  THEN
       message_string = 'use_top_fluxes must be .TRUE. in ocean mode'
       CALL message( 'check_parameters', 'PA0042', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- In case of a given slope, compute the relevant quantities
    IF ( alpha_surface /= 0.0_wp )  THEN
       IF ( ABS( alpha_surface ) > 90.0_wp )  THEN
          WRITE( message_string, * ) 'ABS( alpha_surface = ', alpha_surface,   &
                                     ' ) must be < 90.0'
          CALL message( 'check_parameters', 'PA0043', 1, 2, 0, 6, 0 )
       ENDIF
       sloping_surface = .TRUE.
       cos_alpha_surface = COS( alpha_surface / 180.0_wp * pi )
       sin_alpha_surface = SIN( alpha_surface / 180.0_wp * pi )
    ENDIF

!
!-- Check time step and cfl_factor
    IF ( dt /= -1.0_wp )  THEN
       IF ( dt <= 0.0_wp )  THEN
          WRITE( message_string, * ) 'dt = ', dt , ' <= 0.0'
          CALL message( 'check_parameters', 'PA0044', 1, 2, 0, 6, 0 )
       ENDIF
       dt_3d = dt
       dt_fixed = .TRUE.
    ENDIF

    IF ( cfl_factor <= 0.0_wp  .OR.  cfl_factor > 1.0_wp )  THEN
       IF ( cfl_factor == -1.0_wp )  THEN
          IF ( timestep_scheme == 'runge-kutta-2' )  THEN
             cfl_factor = 0.8_wp
          ELSEIF ( timestep_scheme == 'runge-kutta-3' )  THEN
             cfl_factor = 0.9_wp
          ELSE
             cfl_factor = 0.9_wp
          ENDIF
       ELSE
          WRITE( message_string, * ) 'cfl_factor = ', cfl_factor,              &
                 ' out of range &0.0 < cfl_factor <= 1.0 is required'
          CALL message( 'check_parameters', 'PA0045', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Store simulated time at begin
    simulated_time_at_begin = simulated_time

!
!-- Store reference time for coupled runs and change the coupling flag,
!-- if ...
    IF ( simulated_time == 0.0_wp )  THEN
       IF ( coupling_start_time == 0.0_wp )  THEN
          time_since_reference_point = 0.0_wp
       ELSEIF ( time_since_reference_point < 0.0_wp )  THEN
          run_coupled = .FALSE.
       ENDIF
    ENDIF
!
!-- In case of using a prandtl-layer, calculated (or prescribed) surface
!-- fluxes have to be used in the diffusion-terms
    IF ( constant_flux_layer )  use_surface_fluxes = .TRUE.
!
!-- Check boundary conditions and set internal variables:
!-- Attention: the lateral boundary conditions have been already checked in
!-- parin
!
!-- Non-cyclic lateral boundaries require the multigrid method and Piascek-
!-- Willimas or Wicker - Skamarock advection scheme. Several schemes
!-- and tools do not work with non-cyclic boundary conditions.
    IF ( bc_lr /= 'cyclic'  .OR.  bc_ns /= 'cyclic' )  THEN
       IF ( psolver(1:9) /= 'multigrid' )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'psolver = "' // TRIM( psolver ) // '"'
          CALL message( 'check_parameters', 'PA0051', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( momentum_advec /= 'pw-scheme'  .AND.                               &
            momentum_advec /= 'ws-scheme' )  THEN

          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'momentum_advec = "' // TRIM( momentum_advec ) // '"'
          CALL message( 'check_parameters', 'PA0052', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( scalar_advec /= 'pw-scheme'  .AND.                                 &
            scalar_advec /= 'ws-scheme' )  THEN
          message_string = 'non-cyclic lateral boundaries do not allow ' //    &
                           'scalar_advec = "' // TRIM( scalar_advec ) // '"'
          CALL message( 'check_parameters', 'PA0053', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Bottom boundary condition for the turbulent Kinetic energy
    IF ( bc_e_b == 'neumann' )  THEN
       ibc_e_b = 1
    ELSEIF ( bc_e_b == '(u*)**2+neumann' )  THEN
       ibc_e_b = 2
       IF ( .NOT. constant_flux_layer )  THEN
          bc_e_b = 'neumann'
          ibc_e_b = 1
          message_string = 'boundary condition bc_e_b changed to "' //         &
                           TRIM( bc_e_b ) // '"'
          CALL message( 'check_parameters', 'PA0057', 0, 1, 0, 6, 0 )
       ENDIF
    ELSE
       message_string = 'unknown boundary condition: bc_e_b = "' //            &
                        TRIM( bc_e_b ) // '"'
       CALL message( 'check_parameters', 'PA0058', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for perturbation pressure
    IF ( bc_p_b == 'dirichlet' )  THEN
       ibc_p_b = 0
    ELSEIF ( bc_p_b == 'neumann' )  THEN
       ibc_p_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_b = "' //            &
                        TRIM( bc_p_b ) // '"'
       CALL message( 'check_parameters', 'PA0059', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( bc_p_t == 'dirichlet' )  THEN
       ibc_p_t = 0
!-- TO_DO: later set bc_p_t to neumann before, in case of nested domain
    ELSEIF ( bc_p_t == 'neumann' )  THEN
       ibc_p_t = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_t = "' //            &
                        TRIM( bc_p_t ) // '"'
       CALL message( 'check_parameters', 'PA0061', 1, 2, 0, 6, 0 )
    ENDIF

      IF ( bc_pt_b == 'dirichlet' )  THEN
          ibc_pt_b = 0
       ELSEIF ( bc_pt_b == 'neumann' )  THEN
          ibc_pt_b = 1
       ELSE
          message_string = 'unknown boundary condition: bc_pt_b = "' //        &
                           TRIM( bc_pt_b ) // '"'
          CALL message( 'check_parameters', 'PA0062', 1, 2, 0, 6, 0 )
       ENDIF

    IF ( bc_pt_t == 'dirichlet' )  THEN
       ibc_pt_t = 0
    ELSEIF ( bc_pt_t == 'neumann' )  THEN
       ibc_pt_t = 1
    ELSEIF ( bc_pt_t == 'initial_gradient' )  THEN
       ibc_pt_t = 2
    ELSE
       message_string = 'unknown boundary condition: bc_pt_t = "' //           &
                        TRIM( bc_pt_t ) // '"'
       CALL message( 'check_parameters', 'PA0063', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( ANY( wall_heatflux /= 0.0_wp )  .AND.                        &
         surface_heatflux == 9999999.9_wp )  THEN
       message_string = 'wall_heatflux additionally requires ' //     &
                        'setting of surface_heatflux'
       CALL message( 'check_parameters', 'PA0443', 1, 2, 0, 6, 0 )
    ENDIF

!
!   This IF clause needs revision, got too complex!!
    IF ( surface_heatflux == 9999999.9_wp  )  THEN
       constant_heatflux = .FALSE.
       IF ( large_scale_forcing  .OR.  land_surface  .OR.  urban_surface )  THEN
          IF ( ibc_pt_b == 0 )  THEN
             constant_heatflux = .FALSE.
          ELSEIF ( ibc_pt_b == 1 )  THEN
             constant_heatflux = .TRUE.
             surface_heatflux = 0.0_wp
          ENDIF
       ENDIF
    ELSE
       constant_heatflux = .TRUE.
    ENDIF

    IF ( top_heatflux     == 9999999.9_wp )  constant_top_heatflux = .FALSE.

    IF ( neutral )  THEN

       IF ( surface_heatflux /= 0.0_wp  .AND.                                  &
            surface_heatflux /= 9999999.9_wp )  THEN
          message_string = 'heatflux must not be set for pure neutral flow'
          CALL message( 'check_parameters', 'PA0351', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( top_heatflux /= 0.0_wp  .AND.  top_heatflux /= 9999999.9_wp )      &
       THEN
          message_string = 'heatflux must not be set for pure neutral flow'
          CALL message( 'check_parameters', 'PA0351', 1, 2, 0, 6, 0 )
       ENDIF

    ENDIF

    IF ( top_momentumflux_u /= 9999999.9_wp  .AND.                             &
         top_momentumflux_v /= 9999999.9_wp )  THEN
       constant_top_momentumflux = .TRUE.
    ELSEIF (  .NOT. ( top_momentumflux_u == 9999999.9_wp  .AND.                &
           top_momentumflux_v == 9999999.9_wp ) )  THEN
       message_string = 'both, top_momentumflux_u AND top_momentumflux_v ' //  &
                        'must be set'
       CALL message( 'check_parameters', 'PA0064', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given surface temperature implies Dirichlet boundary condition for
!-- temperature. In this case specification of a constant heat flux is
!-- forbidden.
    IF ( ibc_pt_b == 0  .AND.  constant_heatflux  .AND.                        &
         surface_heatflux /= 0.0_wp )  THEN
       message_string = 'boundary_condition: bc_pt_b = "' // TRIM( bc_pt_b ) //&
                        '& is not allowed with constant_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PA0065', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( constant_heatflux  .AND.  pt_surface_initial_change /= 0.0_wp )  THEN
       WRITE ( message_string, * )  'constant_heatflux = .TRUE. is not allo',  &
               'wed with pt_surface_initial_change (/=0) = ',                  &
               pt_surface_initial_change
       CALL message( 'check_parameters', 'PA0066', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given temperature at the top implies Dirichlet boundary condition for
!-- temperature. In this case specification of a constant heat flux is
!-- forbidden.
    IF ( ibc_pt_t == 0  .AND.  constant_top_heatflux  .AND.                    &
         top_heatflux /= 0.0_wp )  THEN
       message_string = 'boundary_condition: bc_pt_t = "' // TRIM( bc_pt_t ) //&
                        '" is not allowed with constant_top_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PA0067', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for salinity
       IF ( bc_sa_t == 'dirichlet' )  THEN
          ibc_sa_t = 0
       ELSEIF ( bc_sa_t == 'neumann' )  THEN
          ibc_sa_t = 1
       ELSE
          message_string = 'unknown boundary condition: bc_sa_t = "' //        &
                           TRIM( bc_sa_t ) // '"'
          CALL message( 'check_parameters', 'PA0068', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bc_sa_b == 'dirichlet' )  THEN
          ibc_sa_b = 0
!          CALL location_message('ib_sa_b assigned for dirichlet conditions',.TRUE.) !CB
       ELSEIF ( bc_sa_b == 'neumann' )  THEN
          ibc_sa_b = 1
!          CALL location_message('ib_sa_b assigned for neumann conditions',.TRUE.) !CB
       ELSE
          message_string = 'unknown boundary condition: bc_sa_b = "' //        &
                           TRIM( bc_sa_b ) // '"'
          CALL message( 'check_parameters', 'PA0068', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( ibc_sa_t == 1  .AND.  top_salinityflux == 9999999.9_wp )  THEN
          message_string = 'boundary condition: bc_sa_t = "' //                &
                           TRIM( bc_sa_t ) // '" requires to set ' //          &
                           'top_salinityflux'
          CALL message( 'check_parameters', 'PA0069', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( ibc_sa_b == 1  .AND.  bottom_salinityflux == 9999999.9_wp )  THEN
          message_string = 'boundary condition: bc_sa_b = "' //                &
                           TRIM( bc_sa_b ) // '" requires to set ' //          &
                           'bottom_salinityflux'
          CALL message( 'check_parameters', 'PA0069', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    A fixed salinity at the top implies Dirichlet boundary condition for
!--    salinity. In this case specification of a constant salinity flux is
!--    forbidden.
       IF ( top_salinityflux == 9999999.9_wp )  constant_top_salinityflux = .FALSE.
       IF ( ibc_sa_t == 0  .AND.  constant_top_salinityflux  .AND.             &
            top_salinityflux /= 0.0_wp )  THEN
          message_string = 'boundary condition: bc_sa_t = "' //                &
                           TRIM( bc_sa_t ) // '" is not allowed with ' //      &
                           'top_salinityflux /= 0.0'
          CALL message( 'check_parameters', 'PA0070', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( bottom_salinityflux == 9999999.9_wp )  constant_bottom_salinityflux = .FALSE.
       IF ( ibc_sa_b == 0  .AND.  constant_bottom_salinityflux  .AND.             &
            bottom_salinityflux /= 0.0_wp )  THEN
          message_string = 'boundary condition: bc_sa_b = "' //                &
                           TRIM( bc_sa_t ) // '" is not allowed with ' //      &
                           'bottom_salinityflux /= 0.0'
          CALL message( 'check_parameters', 'PA0070', 1, 2, 0, 6, 0 )
       ENDIF
!
!-- Boundary conditions for horizontal components of wind speed
    IF ( bc_uv_b == 'dirichlet' )  THEN
       ibc_uv_b = 0
    ELSEIF ( bc_uv_b == 'neumann' )  THEN
       ibc_uv_b = 1
       IF ( constant_flux_layer )  THEN
          message_string = 'boundary condition: bc_uv_b = "' //                &
               TRIM( bc_uv_b ) // '" is not allowed with constant_flux_layer'  &
               // ' = .TRUE.'
          CALL message( 'check_parameters', 'PA0075', 1, 2, 0, 6, 0 )
       ENDIF
    ELSE
       message_string = 'unknown boundary condition: bc_uv_b = "' //           &
                        TRIM( bc_uv_b ) // '"'
       CALL message( 'check_parameters', 'PA0076', 1, 2, 0, 6, 0 )
    ENDIF
      IF ( bc_uv_t == 'dirichlet' .OR. bc_uv_t == 'dirichlet_0' )  THEN
          ibc_uv_t = 0
          IF ( bc_uv_t == 'dirichlet_0' )  THEN
!
!--          Velocities for the initial u,v-profiles are set zero at the top
!--          in case of dirichlet_0 conditions
             u_init(nzt+1)    = 0.0_wp
             v_init(nzt+1)    = 0.0_wp
          ENDIF
       ELSEIF ( bc_uv_t == 'neumann' )  THEN
          ibc_uv_t = 1
       ELSEIF ( bc_uv_t == 'nested'  .OR.  bc_uv_t == 'forcing' )  THEN
          ibc_uv_t = 3
       ELSE
          message_string = 'unknown boundary condition: bc_uv_t = "' //        &
                           TRIM( bc_uv_t ) // '"'
          CALL message( 'check_parameters', 'PA0077', 1, 2, 0, 6, 0 )
       ENDIF

!
!-- Compute and check, respectively, the Rayleigh Damping parameter
    IF ( rayleigh_damping_factor == -1.0_wp )  THEN
       rayleigh_damping_factor = 0.0_wp
    ELSE
       IF ( rayleigh_damping_factor < 0.0_wp  .OR.                             &
            rayleigh_damping_factor > 1.0_wp )  THEN
          WRITE( message_string, * )  'rayleigh_damping_factor = ',            &
                              rayleigh_damping_factor, ' out of range [0.0,1.0]'
          CALL message( 'check_parameters', 'PA0078', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( rayleigh_damping_height == -1.0_wp )  THEN
       IF (  .NOT.  ocean )  THEN
          rayleigh_damping_height = 0.66666666666_wp * zu(nzt)
       ELSE
          rayleigh_damping_height = 0.66666666666_wp * zu(nzb)
       ENDIF
    ELSE
         IF ( rayleigh_damping_height > 0.0_wp  .OR.                          &
               rayleigh_damping_height < zu(nzb) )  THEN
             WRITE( message_string, * )  'rayleigh_damping_height = ',         &
                   rayleigh_damping_height, ' out of range [0.0,', zu(nzb), ']'
             CALL message( 'check_parameters', 'PA0079', 1, 2, 0, 6, 0 )
          ENDIF
    ENDIF

!
!-- Check number of chosen statistic regions
    IF ( statistic_regions < 0 )  THEN
       WRITE ( message_string, * ) 'number of statistic_regions = ',           &
                   statistic_regions+1, ' is not allowed'
       CALL message( 'check_parameters', 'PA0082', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( normalizing_region > statistic_regions  .OR.                          &
         normalizing_region < 0)  THEN
       WRITE ( message_string, * ) 'normalizing_region = ',                    &
                normalizing_region, ' must be >= 0 and <= ',statistic_regions, &
                ' (value of statistic_regions)'
       CALL message( 'check_parameters', 'PA0083', 1, 2, 0, 6, 0 )
    ENDIF

#if defined(__fullPalm)
!
!-- Set the default intervals for data output, if necessary
!-- NOTE: dt_dosp has already been set in spectra_parin
    IF ( dt_data_output /= 9999999.9_wp )  THEN
       IF ( dt_dopr           == 9999999.9_wp )  dt_dopr           = dt_data_output
       IF ( dt_dopts          == 9999999.9_wp )  dt_dopts          = dt_data_output
       IF ( dt_do3d           == 9999999.9_wp )  dt_do3d           = dt_data_output
       IF ( dt_data_output_av == 9999999.9_wp )  dt_data_output_av = dt_data_output
    ENDIF

!
!-- Set the default skip time intervals for data output, if necessary
    IF ( skip_time_dopr    == 9999999.9_wp )                                   &
                                       skip_time_dopr    = skip_time_data_output
    IF ( skip_time_do3d    == 9999999.9_wp )                                   &
                                       skip_time_do3d    = skip_time_data_output
    IF ( skip_time_data_output_av == 9999999.9_wp )                            &
                                skip_time_data_output_av = skip_time_data_output

!
!-- Check the average intervals (first for 3d-data, then for profiles)
    IF ( averaging_interval > dt_data_output_av )  THEN
       WRITE( message_string, * )  'averaging_interval = ',                    &
             averaging_interval, ' must be <= dt_data_output_av = ',           &
             dt_data_output_av
       CALL message( 'check_parameters', 'PA0085', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( averaging_interval_pr == 9999999.9_wp )  THEN
       averaging_interval_pr = averaging_interval
    ENDIF

    IF ( averaging_interval_pr > dt_dopr )  THEN
       WRITE( message_string, * )  'averaging_interval_pr = ',                 &
             averaging_interval_pr, ' must be <= dt_dopr = ', dt_dopr
       CALL message( 'check_parameters', 'PA0086', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default interval for profiles entering the temporal average
    IF ( dt_averaging_input_pr == 9999999.9_wp )  THEN
       dt_averaging_input_pr = dt_averaging_input
    ENDIF

!
!-- Set the default interval for the output of timeseries to a reasonable
!-- value (tries to minimize the number of calls of flow_statistics)
    IF ( dt_dots == 9999999.9_wp )  THEN
       IF ( averaging_interval_pr == 0.0_wp )  THEN
          dt_dots = MIN( dt_run_control, dt_dopr )
       ELSE
          dt_dots = MIN( dt_run_control, dt_averaging_input_pr )
       ENDIF
    ENDIF

!
!-- Check the sample rate for averaging (first for 3d-data, then for profiles)
    IF ( dt_averaging_input > averaging_interval )  THEN
       WRITE( message_string, * )  'dt_averaging_input = ',                    &
                dt_averaging_input, ' must be <= averaging_interval = ',       &
                averaging_interval
       CALL message( 'check_parameters', 'PA0088', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( dt_averaging_input_pr > averaging_interval_pr )  THEN
       WRITE( message_string, * )  'dt_averaging_input_pr = ',                 &
                dt_averaging_input_pr, ' must be <= averaging_interval_pr = ', &
                averaging_interval_pr
       CALL message( 'check_parameters', 'PA0089', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Determine the number of output profiles and check whether they are
!-- permissible
    DO  WHILE ( data_output_pr(dopr_n+1) /= '          ' )

       dopr_n = dopr_n + 1
       i = dopr_n

!
!--    Determine internal profile number (for hom, homs)
!--    and store height levels
       SELECT CASE ( TRIM( data_output_pr(i) ) )

          CASE ( 'u' )
             dopr_index(i) = 1
             dopr_unit(i)  = 'm/s'
             hom(:,2,1,:)  = SPREAD( zu, 2, statistic_regions+1 )
          CASE ( 'v' )
             dopr_index(i) = 2
             dopr_unit(i)  = 'm/s'
             hom(:,2,2,:)  = SPREAD( zu, 2, statistic_regions+1 )
          CASE ( 'w' )
             dopr_index(i) = 3
             dopr_unit(i)  = 'm/s'
             hom(:,2,3,:)  = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'pt' )
                dopr_index(i) = 4
                dopr_unit(i)  = 'K'
                hom(:,2,4,:)  = SPREAD( zu, 2, statistic_regions+1 )
          CASE ( 'w"u"' )
             dopr_index(i) = 6
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,6,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,6,:) = zu(1)

          CASE ( 'w*u*' )
             dopr_index(i) = 10
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,10,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"v"' )
             dopr_index(i) = 7
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,7,:) = SPREAD( zw, 2, statistic_regions+1 )
             IF ( constant_flux_layer )  hom(nzb,2,7,:) = zu(1)

          CASE ( 'w*v*' )
             dopr_index(i) = 11
             dopr_unit(i)  = TRIM ( momentumflux_output_unit )
             hom(:,2,11,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w"pt"' )
             dopr_index(i) = 8
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,8,:) = SPREAD( zw, 2, statistic_regions+1 )

          CASE ( 'w*pt*' )
             dopr_index(i) = 12
             dopr_unit(i)  = TRIM ( heatflux_output_unit )
             hom(:,2,12,:) = SPREAD( zw, 2, statistic_regions+1 )
          CASE ( 'sa' )
             IF ( .NOT. ocean )  THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 5
                dopr_unit(i)  = 'psu'
                hom(:,2,5,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF
          CASE ( 'rho_ocean' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 14
                dopr_unit(i)  = 'kg/m3'
                hom(:,2,14,:) = SPREAD( zu, 2, statistic_regions+1 )
             ENDIF
          CASE ( 'w"sa"' )
             IF (  .NOT.  ocean ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 9 
                dopr_unit(i)  = 'psu m/s'
                hom(:,2,9,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF

          CASE ( 'w*sa*' )
             IF (  .NOT. ocean  ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for ocean = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 13
                dopr_unit(i)  = 'psu m/s'
                hom(:,2,13,:) = SPREAD( zw, 2, statistic_regions+1 )
             ENDIF


          CASE DEFAULT
             IF ( unit == 'illegal' )  THEN
                IF ( data_output_pr_user(1) /= ' ' )  THEN
                   message_string = 'illegal value for data_output_pr or ' //  &
                                    'data_output_pr_user = "' //               &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0097', 1, 2, 0, 6, 0 )
                ELSE
                   message_string = 'illegal value for data_output_pr = "' //  &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0098', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

       END SELECT

    ENDDO
!
!-- Check and set steering parameters for 2d/3d data output and averaging
    i   = 1
    DO  WHILE ( data_output(i) /= ' '  .AND.  i <= 500 )
!
!--    Check for data averaging
       ilen = LEN_TRIM( data_output(i) )
       j = 0                                                 ! no data averaging
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_av' )  THEN
             j = 1                                           ! data averaging
             data_output(i) = data_output(i)(1:ilen-3)
          ENDIF
       ENDIF
!
!--    Check for cross section or volume data
       ilen = LEN_TRIM( data_output(i) )
       k = 0                                                   ! 3d data
       var = data_output(i)(1:ilen)
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_xy'  .OR.                      &
               data_output(i)(ilen-2:ilen) == '_xz'  .OR.                      &
               data_output(i)(ilen-2:ilen) == '_yz' )  THEN
             k = 1                                             ! 2d data
             var = data_output(i)(1:ilen-3)
          ENDIF

!--       Make sure that section_nn is defined
          IF ( (data_output(i)(ilen-2:ilen) == '_xy') .AND. ( ALL(section_xy == -9999) ) ) THEN
             message_string = 'to output _xy variables, the depth of at least one xy_section must be specified via the namelist parameter section_xy'
             CALL message( 'check_parameters', 'PA0561', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( (data_output(i)(ilen-2:ilen) == '_xz') .AND. ( ALL(section_xz == -9999) ) ) THEN
             message_string = 'to output _xz variables, the depth of at least one xz_section must be specified via the namelist parameter section_xz'
             CALL message( 'check_parameters', 'PA0561', 1, 2, 0, 6, 0 )
          ENDIF
          IF ( (data_output(i)(ilen-2:ilen) == '_yz') .AND. ( ALL(section_yz == -9999) ) ) THEN
             message_string = 'to output _yz variables, the depth of at least one yz_section must be specified via the namelist parameter section_yz'
             CALL message( 'check_parameters', 'PA0561', 1, 2, 0, 6, 0 )
          ENDIF

       ENDIF

!
!--    Check for allowed value and set units
       SELECT CASE ( TRIM( var ) )

          CASE ( 'e' )
            unit = 'm2/s2'

          CASE ( 'solar3d' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '^oC s^{-1}'


          CASE ( 'rho_ocean' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/m3'

          CASE ( 'alpha_T' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = '^oC^{-1}'

          CASE ( 'beta_S' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'PSU^{-1}'

          CASE ( 's' )
             IF (  .NOT.  passive_scalar )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res passive_scalar = .TRUE.'
                CALL message( 'check_parameters', 'PA0110', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'kg/m3'

          CASE ( 'sa' )
             IF (  .NOT.  ocean )  THEN
                message_string = 'output of "' // TRIM( var ) // '" requi' //  &
                                 'res ocean = .TRUE.'
                CALL message( 'check_parameters', 'PA0109', 1, 2, 0, 6, 0 )
             ENDIF
             unit = 'psu'

          CASE ( 'p', 'pt', 'u', 'v', 'w' )
             IF ( TRIM( var ) == 'p'  )  unit = 'Pa'
             IF ( TRIM( var ) == 'pt' )  unit = 'K'
             IF ( TRIM( var ) == 'u'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'v'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'w'  )  unit = 'm/s'
             CONTINUE

          CASE ( 'ghf*', 'lwp*', 'ol*', 'pra*', 'prr*', 'qsws*', 'r_a*',       &
                 'shf_sol*', 'shf*', 'ssws*', 't*', 'tsurf*', 'u*', 'z0*', 'z0h*', 'z0q*' )
             IF ( k == 0  .OR.  data_output(i)(ilen-2:ilen) /= '_xy' )  THEN
                message_string = 'illegal value for data_output: "' //         &
                                 TRIM( var ) // '" & only 2d-horizontal ' //   &
                                 'cross sections are allowed for this value'
                CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
             ENDIF
            IF ( TRIM( var ) == 'ghf*'   )  unit = 'W/m2'
             IF ( TRIM( var ) == 'lwp*'   )  unit = 'kg/m2'
             IF ( TRIM( var ) == 'ol*'    )  unit = 'm'
             IF ( TRIM( var ) == 'pra*'   )  unit = 'mm'
             IF ( TRIM( var ) == 'prr*'   )  unit = 'mm/s'
             IF ( TRIM( var ) == 'qsws*'  )  unit = 'kgm/kgs'
             IF ( TRIM( var ) == 'r_a*'   )  unit = 's/m'
             IF ( TRIM( var ) == 'shf*'   )  unit = 'K*m/s'
             IF ( TRIM( var ) == 'shf_sol*' ) unit = 'K*m/s'
             IF ( TRIM( var ) == 'ssws*'  )  unit = 'kg/m2*s'
             IF ( TRIM( var ) == 't*'     )  unit = 'K'
             IF ( TRIM( var ) == 'tsurf*' )  unit = 'K'
             IF ( TRIM( var ) == 'u*'     )  unit = 'm/s'
             IF ( TRIM( var ) == 'z0*'    )  unit = 'm'
             IF ( TRIM( var ) == 'z0h*'   )  unit = 'm'
          CASE DEFAULT

             CALL tcm_check_data_output ( var, unit, i, ilen, k )

       END SELECT
!
!--    Set the internal steering parameters appropriately
       IF ( k == 0 )  THEN
          do3d_no(j)              = do3d_no(j) + 1
          do3d(j,do3d_no(j))      = data_output(i)
          do3d_unit(j,do3d_no(j)) = unit
       ENDIF

       IF ( j == 1 )  THEN
!
!--       Check, if variable is already subject to averaging
          found = .FALSE.
          DO  k = 1, doav_n
             IF ( TRIM( doav(k) ) == TRIM( var ) )  found = .TRUE.
          ENDDO

          IF ( .NOT. found )  THEN
             doav_n = doav_n + 1
             doav(doav_n) = var
          ENDIF
       ENDIF

       i = i + 1
    ENDDO

!
!-- Averaged 2d or 3d output requires that an averaging interval has been set
    IF ( doav_n > 0  .AND.  averaging_interval == 0.0_wp )  THEN
       WRITE( message_string, * )  'output of averaged quantity "',            &
                                   TRIM( doav(1) ), '_av" requires to set a ', &
                                   'non-zero averaging interval'
       CALL message( 'check_parameters', 'PA0323', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check sectional planes and store them in one shared array
    IF ( ANY( section_xy > nz + 1 ) )  THEN
       WRITE( message_string, * )  'section_xy must be <= nz + 1 = ', nz + 1
       CALL message( 'check_parameters', 'PA0319', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( section_xz > ny + 1 ) )  THEN
       WRITE( message_string, * )  'section_xz must be <= ny + 1 = ', ny + 1
       CALL message( 'check_parameters', 'PA0320', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( ANY( section_yz > nx + 1 ) )  THEN
       WRITE( message_string, * )  'section_yz must be <= nx + 1 = ', nx + 1
       CALL message( 'check_parameters', 'PA0321', 1, 2, 0, 6, 0 )
    ENDIF
    section(:,1) = section_xy
    section(:,2) = section_xz
    section(:,3) = section_yz

!
!-- Upper plot limit for 3D arrays
    IF ( nz_do3d == -9999 )  nz_do3d = nzt + 1

!
!-- Set output format string (used in header)
    SELECT CASE ( netcdf_data_format )
       CASE ( 1 )
          netcdf_data_format_string = 'netCDF classic'
       CASE ( 2 )
          netcdf_data_format_string = 'netCDF 64bit offset'
       CASE ( 3 )
          netcdf_data_format_string = 'netCDF4/HDF5'
       CASE ( 4 )
          netcdf_data_format_string = 'netCDF4/HDF5 classic'
       CASE ( 5 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5'
       CASE ( 6 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5 classic'

    END SELECT

!
!-- Check the NetCDF data format
    IF ( netcdf_data_format > 2 )  THEN
#if defined( __netcdf4 )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 format requested but no ' //          &
                        'cpp-directive __netcdf4 given & switch '  //          &
                        'back to 64-bit offset format'
       CALL message( 'check_parameters', 'PA0171', 0, 1, 0, 6, 0 )
       netcdf_data_format = 2
#endif
    ENDIF
    IF ( netcdf_data_format > 4 )  THEN
#if defined( __netcdf4 ) && defined( __netcdf4_parallel )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 parallel output requested but no ' // &
                        'cpp-directive __netcdf4_parallel given & switch '   //&
                        'back to netCDF4 non-parallel output'
       CALL message( 'check_parameters', 'PA0099', 0, 1, 0, 6, 0 )
       netcdf_data_format = netcdf_data_format - 2
#endif
    ENDIF

!
!-- Calculate fixed number of output time levels for parallel netcdf output.
!-- The time dimension has to be defined as limited for parallel output,
!-- because otherwise the I/O performance drops significantly.
    IF ( netcdf_data_format > 4 )  THEN

!
!--    Check if any of the follwoing data output interval is 0.0s, which is
!--    not allowed for parallel output.
       CALL check_dt_do( dt_do3d,           'dt_do3d'           )
       CALL check_dt_do( dt_data_output_av, 'dt_data_output_av' )

!--    Set needed time levels (ntdim) to
!--    saved time levels + to be saved time levels.
       ntdim_3d(0) = do3d_time_count(0) + CEILING(                             &
                     ( end_time - MAX( skip_time_do3d,                         &
                                       simulated_time_at_begin )               &
                     ) / dt_do3d )
       IF ( do3d_at_begin ) ntdim_3d(0) = ntdim_3d(0) + 1

       ntdim_3d(1) = do3d_time_count(1) + CEILING(                             &
                     ( end_time - MAX( skip_time_data_output_av,               &
                                       simulated_time_at_begin )               &
                     ) / dt_data_output_av )
    ENDIF
#endif

!
!
!-- In case of non-cyclic lateral boundaries and a damping layer for the
!-- potential temperature, check the width of the damping layer
    IF ( bc_lr /= 'cyclic' ) THEN
       IF ( pt_damping_width < 0.0_wp  .OR.                                    &
            pt_damping_width > REAL( (nx+1) * dx ) )  THEN
          message_string = 'pt_damping_width out of range'
          CALL message( 'check_parameters', 'PA0124', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( bc_ns /= 'cyclic' )  THEN
       IF ( pt_damping_width < 0.0_wp  .OR.                                    &
            pt_damping_width > REAL( (ny+1) * dy ) )  THEN
          message_string = 'pt_damping_width out of range'
          CALL message( 'check_parameters', 'PA0124', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Check value range for zeta = z/L
    IF ( zeta_min >= zeta_max )  THEN
       WRITE( message_string, * )  'zeta_min = ', zeta_min, ' must be less ',  &
                                   'than zeta_max = ', zeta_max
       CALL message( 'check_parameters', 'PA0125', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check random generator
    IF ( (random_generator /= 'system-specific'      .AND.                     &
          random_generator /= 'random-parallel'   )  .AND.                     &
          random_generator /= 'numerical-recipes' )  THEN
       message_string = 'unknown random generator: random_generator = "' //    &
                        TRIM( random_generator ) // '"'
       CALL message( 'check_parameters', 'PA0135', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine upper and lower hight level indices for random perturbations
    IF ( disturbance_level_b == -9999999.9_wp )  THEN
       IF ( ocean )  THEN
          disturbance_level_b     = zu((nzt*2)/3)
          disturbance_level_ind_b = ( nzt * 2 ) / 3
       ELSE
          disturbance_level_b     = zu(nzb+3)
          disturbance_level_ind_b = nzb + 3
       ENDIF
    ELSEIF ( disturbance_level_b < zu(3) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ',                   &
                           disturbance_level_b, ' must be >= ', zu(3), '(zu(3))'
       CALL message( 'check_parameters', 'PA0126', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_b > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ',                   &
                   disturbance_level_b, ' must be <= ', zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PA0127', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_b <= zu(k) )  THEN
             disturbance_level_ind_b = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( disturbance_level_t == -9999999.9_wp )  THEN
       IF ( ocean )  THEN
          disturbance_level_t     = zu(nzt-3)
          disturbance_level_ind_t = nzt - 3
       ELSE
          disturbance_level_t     = zu(nzt/3)
          disturbance_level_ind_t = nzt / 3
       ENDIF
!    ELSEIF ( disturbance_level_t > zu(nzt-2) )  THEN
!       WRITE( message_string, * )  'disturbance_level_t = ',                   &
!                   disturbance_level_t, ' must be <= ', zu(nzt-2), '(zu(nzt-2))'
!       CALL message( 'check_parameters', 'PA0128', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_t < disturbance_level_b )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ',                   &
                   disturbance_level_t, ' must be >= disturbance_level_b = ',  &
                   disturbance_level_b
       CALL message( 'check_parameters', 'PA0129', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt
          IF ( disturbance_level_t <= zu(k+1) )  THEN
             disturbance_level_ind_t = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

!
!-- Check again whether the levels determined this way are ok.
!-- Error may occur at automatic determination and too few grid points in
!-- z-direction.
    IF ( disturbance_level_ind_t < disturbance_level_ind_b )  THEN
       WRITE( message_string, * )  'disturbance_level_ind_t = ',               &
                disturbance_level_ind_t, ' must be >= ',                       &
                'disturbance_level_ind_b = ', disturbance_level_ind_b
       CALL message( 'check_parameters', 'PA0130', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine the horizontal index range for random perturbations.
!-- In case of non-cyclic horizontal boundaries, no perturbations are imposed
!-- near the inflow and the perturbation area is further limited to ...(1)
!-- after the initial phase of the flow.

    IF ( bc_lr /= 'cyclic' )  THEN
       IF ( inflow_disturbance_begin == -1 )  THEN
          inflow_disturbance_begin = MIN( 10, nx/2 )
       ENDIF
       IF ( inflow_disturbance_begin < 0  .OR.  inflow_disturbance_begin > nx )&
       THEN
          message_string = 'inflow_disturbance_begin out of range'
          CALL message( 'check_parameters', 'PA0131', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inflow_disturbance_end == -1 )  THEN
          inflow_disturbance_end = MIN( 100, 3*nx/4 )
       ENDIF
       IF ( inflow_disturbance_end < 0  .OR.  inflow_disturbance_end > nx )    &
       THEN
          message_string = 'inflow_disturbance_end out of range'
          CALL message( 'check_parameters', 'PA0132', 1, 2, 0, 6, 0 )
       ENDIF
    ELSEIF ( bc_ns /= 'cyclic' )  THEN
       IF ( inflow_disturbance_begin == -1 )  THEN
          inflow_disturbance_begin = MIN( 10, ny/2 )
       ENDIF
       IF ( inflow_disturbance_begin < 0  .OR.  inflow_disturbance_begin > ny )&
       THEN
          message_string = 'inflow_disturbance_begin out of range'
          CALL message( 'check_parameters', 'PA0131', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( inflow_disturbance_end == -1 )  THEN
          inflow_disturbance_end = MIN( 100, 3*ny/4 )
       ENDIF
       IF ( inflow_disturbance_end < 0  .OR.  inflow_disturbance_end > ny )    &
       THEN
          message_string = 'inflow_disturbance_end out of range'
          CALL message( 'check_parameters', 'PA0132', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( random_generator == 'random-parallel' )  THEN
       dist_nxl = nxl;  dist_nxr = nxr
       dist_nys = nys;  dist_nyn = nyn
       IF ( bc_lr == 'radiation/dirichlet' )  THEN
          dist_nxr    = MIN( nx - inflow_disturbance_begin, nxr )
          dist_nxl(1) = MAX( nx - inflow_disturbance_end, nxl )
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          dist_nxl    = MAX( inflow_disturbance_begin, nxl )
          dist_nxr(1) = MIN( inflow_disturbance_end, nxr )
      ENDIF
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          dist_nyn    = MIN( ny - inflow_disturbance_begin, nyn )
          dist_nys(1) = MAX( ny - inflow_disturbance_end, nys )
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          dist_nys    = MAX( inflow_disturbance_begin, nys )
          dist_nyn(1) = MIN( inflow_disturbance_end, nyn )
       ENDIF
    ELSE
       dist_nxl = 0;  dist_nxr = nx
       dist_nys = 0;  dist_nyn = ny
       IF ( bc_lr == 'radiation/dirichlet' )  THEN
          dist_nxr    = nx - inflow_disturbance_begin
          dist_nxl(1) = nx - inflow_disturbance_end
       ELSEIF ( bc_lr == 'dirichlet/radiation' )  THEN
          dist_nxl    = inflow_disturbance_begin
          dist_nxr(1) = inflow_disturbance_end
       ENDIF
       IF ( bc_ns == 'dirichlet/radiation' )  THEN
          dist_nyn    = ny - inflow_disturbance_begin
          dist_nys(1) = ny - inflow_disturbance_end
       ELSEIF ( bc_ns == 'radiation/dirichlet' )  THEN
          dist_nys    = inflow_disturbance_begin
          dist_nyn(1) = inflow_disturbance_end
       ENDIF
    ENDIF
!
!-- Check some other 1d-model parameters
    IF ( TRIM( mixing_length_1d ) /= 'as_in_3d_model'  .AND.                   &
         TRIM( mixing_length_1d ) /= 'blackadar' )  THEN
       message_string = 'mixing_length_1d = "' // TRIM( mixing_length_1d ) //  &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0137', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( dissipation_1d ) /= 'as_in_3d_model'  .AND.                     &
         TRIM( dissipation_1d ) /= 'detering'  .AND.                           &
         TRIM( dissipation_1d ) /= 'prognostic' )  THEN
       message_string = 'dissipation_1d = "' // TRIM( dissipation_1d ) //      &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0138', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( TRIM( mixing_length_1d ) /= 'as_in_3d_model'  .AND.                   &
         TRIM( dissipation_1d ) == 'as_in_3d_model' )  THEN
       message_string = 'dissipation_1d = "' // TRIM( dissipation_1d ) //      &
                        '" requires mixing_length_1d = "as_in_3d_model"'
       CALL message( 'check_parameters', 'PA0485', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set time for the next user defined restart (time_restart is the
!-- internal parameter for steering restart events)
    IF ( restart_time /= 9999999.9_wp )  THEN
       IF ( restart_time > time_since_reference_point )  THEN
          time_restart = restart_time
       ENDIF
    ELSE
!
!--    In case of a restart run, set internal parameter to default (no restart)
!--    if the NAMELIST-parameter restart_time is omitted
       time_restart = 9999999.9_wp
    ENDIF

!
!-- Check pressure gradient conditions
    IF ( dp_external  .AND.  conserve_volume_flow )  THEN
       WRITE( message_string, * )  'Both dp_external and conserve_volume_flo', &
            'w are .TRUE. but one of them must be .FALSE.'
       CALL message( 'check_parameters', 'PA0150', 1, 2, 0, 6, 0 )
    ENDIF
    IF ( dp_external )  THEN
       IF ( dp_level_b < zu(nzb)  .OR.  dp_level_b > zu(nzt) )  THEN
          WRITE( message_string, * )  'dp_level_b = ', dp_level_b, ' is out ', &
               ' of range [zu(nzb), zu(nzt)]'
          CALL message( 'check_parameters', 'PA0151', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( .NOT. ANY( dpdxy /= 0.0_wp ) )  THEN
          WRITE( message_string, * )  'dp_external is .TRUE. but dpdxy is ze', &
               'ro, i.e. the external pressure gradient will not be applied'
          CALL message( 'check_parameters', 'PA0152', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( ANY( dpdxy /= 0.0_wp )  .AND.  .NOT.  dp_external )  THEN
       WRITE( message_string, * )  'dpdxy is nonzero but dp_external is ',     &
            '.FALSE., i.e. the external pressure gradient & will not be applied'
       CALL message( 'check_parameters', 'PA0153', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- Prevent empty time records in volume, cross-section and masked data in case
!-- of non-parallel netcdf-output in restart runs
    IF ( netcdf_data_format < 5 )  THEN
       IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
          do3d_time_count    = 0
       ENDIF
    ENDIF

!
!-- Check for valid setting of most_method
    IF ( TRIM( most_method ) /= 'circular'  .AND.                              &
         TRIM( most_method ) /= 'newton'    .AND.                              &
         TRIM( most_method ) /= 'lookup' )  THEN
       message_string = 'most_method = "' // TRIM( most_method ) //            &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0416', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Check roughness length, which has to be smaller than dz/2
    IF ( ( constant_flux_layer .OR.  &
           INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )       &
         .AND. roughness_length >= 0.5 * dz(1) )  THEN
       message_string = 'roughness_length must be smaller than dz/2'
       CALL message( 'check_parameters', 'PA0424', 1, 2, 0, 6, 0 )
    ENDIF

!    CALL location_message( 'finished', .TRUE. )
!

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check the length of data output intervals. In case of parallel NetCDF output
!> the time levels of the output files need to be fixed. Therefore setting the
!> output interval to 0.0s (usually used to output each timestep) is not
!> possible as long as a non-fixed timestep is used.
!------------------------------------------------------------------------------!

    SUBROUTINE check_dt_do( dt_do, dt_do_name )

       IMPLICIT NONE

       CHARACTER (LEN=*), INTENT (IN) :: dt_do_name !< parin variable name

       REAL(wp), INTENT (INOUT)       :: dt_do      !< data output interval

       IF ( dt_do == 0.0_wp )  THEN
          IF ( dt_fixed )  THEN
             WRITE( message_string, '(A,F9.4,A)' )  'Output at every '  //     &
                    'timestep is wanted (' // dt_do_name // ' = 0.0).&'//      &
                    'The output interval is set to the fixed timestep dt '//   &
                    '= ', dt, 's.'
             CALL message( 'check_parameters', 'PA0060', 0, 0, 0, 6, 0 )
             dt_do = dt
          ELSE
             message_string = dt_do_name // ' = 0.0 while using a ' //         &
                              'variable timestep and parallel netCDF4 ' //     &
                              'is not allowed.'
             CALL message( 'check_parameters', 'PA0081', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

    END SUBROUTINE check_dt_do




!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inititalizes the vertical profiles of scalar quantities.
!------------------------------------------------------------------------------!

    SUBROUTINE init_vertical_profiles( vertical_gradient_level_ind,            &
                                       vertical_gradient_level,                &
                                       vertical_gradient,                      &
                                       pr_init, surface_value, bc_t_val )


       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< counter
       INTEGER(iwp), DIMENSION(1:10) ::  vertical_gradient_level_ind !< vertical grid indices for gradient levels

       REAL(wp)     ::  bc_t_val      !< model top gradient
       REAL(wp)     ::  gradient      !< vertica gradient of the respective quantity
       REAL(wp)     ::  surface_value !< surface value of the respecitve quantity


       REAL(wp), DIMENSION(0:nz+1) ::  pr_init                 !< initialisation profile
       REAL(wp), DIMENSION(1:10)   ::  vertical_gradient       !< given vertical gradient
       REAL(wp), DIMENSION(1:10)   ::  vertical_gradient_level !< given vertical gradient level

       i = 1
       gradient = 0.0_wp

       IF (  .NOT.  ocean )  THEN

          vertical_gradient_level_ind(1) = 0
          DO  k = 1, nzt+1
             IF ( i < 11 )  THEN
                IF ( vertical_gradient_level(i) < zu(k)  .AND.            &
                     vertical_gradient_level(i) >= 0.0_wp )  THEN
                   gradient = vertical_gradient(i) / 100.0_wp
                   vertical_gradient_level_ind(i) = k - 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= 1 )  THEN
                   pr_init(k) = pr_init(k-1) + dzu(k) * gradient
                ELSE
                   pr_init(k) = pr_init(k-1) + dzu(k) * gradient
                ENDIF
             ELSE
                pr_init(k) = pr_init(k-1)
             ENDIF
   !
   !--       Avoid negative values
             IF ( pr_init(k) < 0.0_wp )  THEN
                pr_init(k) = 0.0_wp
             ENDIF
          ENDDO

       ELSE

          vertical_gradient_level_ind(1) = nzt+1
          DO  k = nzt, 0, -1
             IF ( i < 11 )  THEN
                IF ( vertical_gradient_level(i) >= zu(k)  .AND.            &
                     vertical_gradient_level(i) <= 0.0_wp )  THEN
                   gradient = vertical_gradient(i) / 100.0_wp
                   vertical_gradient_level_ind(i) = k + 1
                   i = i + 1
                ENDIF
             ENDIF
             IF ( gradient /= 0.0_wp )  THEN
                IF ( k /= nzt )  THEN
                   pr_init(k) = pr_init(k+1) - dzu(k+1) * gradient
                ELSE
                   pr_init(k)   = surface_value - 0.5_wp * dzu(k+1) * gradient
                   pr_init(k+1) = surface_value + 0.5_wp * dzu(k+1) * gradient
                ENDIF
             ELSE
                pr_init(k) = pr_init(k+1)
             ENDIF
   !
   !--       Avoid negative humidities
             IF ( pr_init(k) < 0.0_wp )  THEN
                pr_init(k) = 0.0_wp
             ENDIF
          ENDDO

       ENDIF

!
!--    In case of no given gradients, choose zero gradient conditions
       IF ( vertical_gradient_level(1) == -999999.9_wp )  THEN
          vertical_gradient_level(1) = 0.0_wp
       ENDIF
!
!--    Store gradient at the top boundary for possible Neumann boundary condition
       bc_t_val  = ( pr_init(nzt+1) - pr_init(nzt) ) / dzu(nzt+1)

    END SUBROUTINE init_vertical_profiles



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the bottom and top boundary conditions for humidity and scalars.
!------------------------------------------------------------------------------!

    SUBROUTINE set_bc_scalars( sq, bc_b, bc_t, ibc_b, ibc_t, err_nr_b, err_nr_t )


       IMPLICIT NONE

       CHARACTER (LEN=1)   ::  sq         !< name of scalar quantity
       CHARACTER (LEN=*)   ::  bc_b       !< bottom boundary condition
       CHARACTER (LEN=*)   ::  bc_t       !< top boundary condition
       CHARACTER (LEN=*)   ::  err_nr_b   !< error number if bottom bc is unknown
       CHARACTER (LEN=*)   ::  err_nr_t   !< error number if top bc is unknown

       INTEGER(iwp)        ::  ibc_b      !< index for bottom boundary condition
       INTEGER(iwp)        ::  ibc_t      !< index for top boundary condition

!
!--    Set Integer flags and check for possilbe errorneous settings for bottom
!--    boundary condition
       IF ( bc_b == 'dirichlet' )  THEN
          ibc_b = 0
       ELSEIF ( bc_b == 'neumann' )  THEN
          ibc_b = 1
       ELSE
          message_string = 'unknown boundary condition: bc_' // TRIM( sq ) //  &
                           '_b ="' // TRIM( bc_b ) // '"'
          CALL message( 'check_parameters', err_nr_b, 1, 2, 0, 6, 0 )
       ENDIF
!
!--    Set Integer flags and check for possilbe errorneous settings for top
!--    boundary condition
       IF ( bc_t == 'dirichlet' )  THEN
          ibc_t = 0
       ELSEIF ( bc_t == 'neumann' )  THEN
          ibc_t = 1
       ELSEIF ( bc_t == 'initial_gradient' )  THEN
          ibc_t = 2
       ELSE
          message_string = 'unknown boundary condition: bc_' // TRIM( sq ) //  &
                           '_t ="' // TRIM( bc_t ) // '"'
          CALL message( 'check_parameters', err_nr_t, 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE set_bc_scalars



!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check for consistent settings of bottom boundary conditions for humidity
!> and scalars.
!------------------------------------------------------------------------------!

    SUBROUTINE check_bc_scalars( sq, bc_b, ibc_b,                      &
                                 err_nr_1, err_nr_2,                   &
                                 constant_flux, surface_initial_change )


       IMPLICIT NONE

       CHARACTER (LEN=1)   ::  sq                       !< name of scalar quantity
       CHARACTER (LEN=*)   ::  bc_b                     !< bottom boundary condition
       CHARACTER (LEN=*)   ::  err_nr_1                 !< error number of first error
       CHARACTER (LEN=*)   ::  err_nr_2                 !< error number of second error

       INTEGER(iwp)        ::  ibc_b                    !< index of bottom boundary condition

       LOGICAL             ::  constant_flux            !< flag for constant-flux layer

       REAL(wp)            ::  surface_initial_change   !< value of initial change at the surface

!
!--    A given surface value implies Dirichlet boundary condition for
!--    the respective quantity. In this case specification of a constant flux is
!--    forbidden. However, an exception is made for large-scale forcing as well
!--    as land-surface model.
       IF ( .NOT. land_surface  .AND.  .NOT. large_scale_forcing )  THEN
          IF ( ibc_b == 0  .AND.  constant_flux )  THEN
             message_string = 'boundary condition: bc_' // TRIM( sq ) //       &
                              '_b ' // '= "' // TRIM( bc_b ) //                &
                              '" is not allowed with prescribed surface flux'
             CALL message( 'check_parameters', err_nr_1, 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF
       IF ( constant_flux  .AND.  surface_initial_change /= 0.0_wp )  THEN
          WRITE( message_string, * )  'a prescribed surface flux is not allo', &
                 'wed with ', sq, '_surface_initial_change (/=0) = ',          &
                 surface_initial_change
          CALL message( 'check_parameters', err_nr_2, 1, 2, 0, 6, 0 )
       ENDIF


    END SUBROUTINE check_bc_scalars



 END SUBROUTINE check_parameters
