MODULE palm_mod

   USE advec_ws,                                                               &
       ONLY:  ws_init_flags, ws_finalize

   USE arrays_3d

   USE control_parameters

   USE configure_3d_model

   USE cpulog,                                                                 &
        ONLY:  cpu_log, log_point, log_point_s, cpu_statistics
   USE indices

   USE fft_xy,                                                                 &
        ONLY: fft_finalize

   USE kinds

   USE pegrid

   USE ppr_1d

   USE statistics,                                                             &
       ONLY:  flow_statistics_called, hom, hom_sum, pr_palm

   IMPLICIT NONE

   ! arrays and parameters for PPR remapping
   INTEGER(iwp), PARAMETER :: ndof = 1, nvar = 4
   ! large scale forcing mode
   ! 0: none
   ! 1: use the large scale forcing term from input
   ! 2: nudging to the large scale fields
   INTEGER(iwp), PARAMETER :: LS_MODE_NONE = 0,                                &
                              LS_MODE_INPUT = 1,                               &
                              LS_MODE_NUDGING = 2

   TYPE(rmap_work) :: work
   TYPE(rmap_opts) :: opts
   TYPE(rcon_ends) :: bc_l(nvar), bc_r(nvar)

   ! local variables
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: fldSS, fldLS

   PUBLIC :: palm_init, palm_main, palm_finalize

   CONTAINS


   SUBROUTINE palm_init(dtLS, nzLS, zwLS, ptLS, saLS, uLS, vLS,                &
                        wtflux, wsflux, uwflux, vwflux, wtflux_solar,          &
                        nzSS, nxSS, nySS, dxSS, dySS,                          &
                        cfg_idealized_diurnal,                                 &
                        cfg_ideal_solar_division,                              &
                        cfg_ideal_solar_depth1,                                &
                        cfg_ideal_solar_depth2,                                &
                        cfg_spinup_time,                                       &
                        cfg_LS_mode,                                           &
                        cfg_LS_factor_tracer,                                  &
                        cfg_LS_factor_vel,                                     &
                        cfg_fCoriolis,                                         &
                        cfg_runID)

   IMPLICIT NONE

   ! flag for idealized diurnal forcing
   ! if true, wtflux_solar is treated as the maximum solar temperature flux during a day
   LOGICAL, INTENT(in) :: cfg_idealized_diurnal
   ! Run ID
   INTEGER(iwp), INTENT(in) :: cfg_runID
   ! number of vertical levels in the large scale model
   INTEGER(iwp), INTENT(in) :: nzLS
   ! number of grid cells in the small scale model (PALM)
   INTEGER(iwp), INTENT(in) :: nzSS, nxSS, nySS
   ! large scale forcing mode
   INTEGER(iwp), INTENT(in) :: cfg_LS_mode
   ! horizonal grid size in the small scale model (PALM)
   REAL(wp), INTENT(in) :: dxSS, dySS
   ! time step in the large scale model
   REAL(wp), INTENT(in) :: dtLS
   ! spinup time
   REAL(wp), INTENT(in) :: cfg_spinup_time
   ! large scale forcing factor
   REAL(wp), INTENT(in) :: cfg_LS_factor_tracer, cfg_LS_factor_vel
   ! latitude
   REAL(wp), INTENT(in) :: cfg_fCoriolis
   ! surface fluxes
   REAL(wp), INTENT(in) :: wtflux, wsflux, uwflux, vwflux, wtflux_solar
   ! solar flux absorption parameters
   REAL(wp), INTENT(in) :: cfg_ideal_solar_division, cfg_ideal_solar_depth1, cfg_ideal_solar_depth2
   ! large scale fields (T, S, U, V)
   REAL(wp), DIMENSION(nzLS), INTENT(in) :: ptLS, saLS, uLS, vLS
   ! large scale z at grid interfaces
   REAL(wp), DIMENSION(nzLS+1), INTENT(in) :: zwLS

   ! local variables
   CHARACTER (LEN=9) ::  time_to_string
   INTEGER :: k
   REAL(wp) :: bottom_depth
#ifdef __DEBUG_PALMLIB
   REAL(wp), DIMENSION(nzLS)   :: zuLS
#endif

   !-- Initialize measuring of the CPU-time remaining to the run
   CALL local_tremain_ini
   CALL cpu_log( log_point(1), 'total', 'start' )
   CALL cpu_log( log_point(2), 'initial step', 'start' )

   ! set run id
   run_id = cfg_runID

   ! initialize namelist from input
   CALL parin

   ! set large scale time step
   dt_LS = dtLS
   ! set large scale forcing mode
   mode_LS = cfg_LS_mode
   IF (mode_LS .eq. LS_MODE_NONE) THEN
      factor_LS_tracer = 0.0_wp
      factor_LS_vel = 0.0_wp
   ELSE
      factor_LS_tracer = cfg_LS_factor_tracer
      factor_LS_vel = cfg_LS_factor_vel
   ENDIF

   ! override control parameters with values from large scale model
   nz = nzSS
   nx = nxSS
   ny = nySS
   dx = dxSS
   dy = dySS
   end_time = cfg_spinup_time
   idealized_diurnal = cfg_idealized_diurnal
   averaging_interval_pr = dtLS
   dt_dopr = dtLS

   ! get vertical grid size in PALM
   bottom_depth = zwLS(nzLS+1)
   dz(1) = ABS(bottom_depth)/nz

   ! update surface fluxes
   top_momentumflux_u = uwflux
   top_momentumflux_v = vwflux
   top_heatflux = wtflux
   top_salinityflux = wsflux
   ideal_solar_heatflux = wtflux_solar
   ideal_solar_division = cfg_ideal_solar_division
   ideal_solar_efolding1 = 1.0_wp/cfg_ideal_solar_depth1
   ideal_solar_efolding2 = 1.0_wp/cfg_ideal_solar_depth2

   ! initial statistics array
   ALLOCATE( fldSS(ndof,nvar,nz) )
   ALLOCATE( fldLS(ndof,nvar,nzLS) )

   ! initialize processor topology
   CALL init_pegrid

   ! initialize grid
   CALL init_grid

   ! initialize flags for ws-scheme to degrade order of the numerics near boundaries
   CALL ws_init_flags

   ! check parameters
   CALL check_parameters

   ! update Coriolis parameters
   f = cfg_fCoriolis
   fs = 0.0_wp

   ! interpolate large scale fields on PALM grid as initial conditions
   fldLS(1,1,:) = ptLS
   fldLS(1,2,:) = saLS
   fldLS(1,3,:) = uLS
   fldLS(1,4,:) = vLS

   CALL cpu_log( log_point(9), 'interp_fields', 'start' )
   CALL interp_fields_LS2SS(nzLS, zwLS, fldLS, nz, zw(nzb:nzt), fldSS)
   CALL cpu_log( log_point(9), 'interp_fields', 'stop' )

#ifdef __DEBUG_PALMLIB
   ! z at large scale model cell centers
   DO k = 1, nzLS
      zuLS(k) = 0.5_wp * (zwLS(k) + zwLS(k+1))
   ENDDO
   write(*,*) '(palm_init) DEBUG: fldLS (initial profiles)'
   CALL print_fields(nzLS, zuLS, fldLS)
   write(*,*) '(palm_init) DEBUG: fldSS (initial profiles)'
   CALL print_fields(nz, zu(nzb+1:nzt), fldSS)
#endif

   ! set initial conditions
   pt_init(nzb+1:nzt) = fldSS(1,1,:)
   sa_init(nzb+1:nzt) = fldSS(1,2,:)
   u_init(nzb+1:nzt)  = fldSS(1,3,:)
   v_init(nzb+1:nzt)  = fldSS(1,4,:)
   pt_init(nzt+1) = pt_init(nzt)
   sa_init(nzt+1) = sa_init(nzt)
   u_init(nzt+1)  = u_init(nzt)
   v_init(nzt+1)  = v_init(nzt)
   pt_init(nzb) = pt_init(nzb+1)
   sa_init(nzb) = sa_init(nzb+1)
   u_init(nzb)  = u_init(nzb+1)
   v_init(nzb)  = v_init(nzb+1)

   ! initalize parameters
   CALL init_3d_model

   ! large scale forcing is zero at initial step
   pt_LS_forcing = 0.0_wp
   sa_LS_forcing = 0.0_wp
   u_LS_forcing  = 0.0_wp
   v_LS_forcing  = 0.0_wp

#ifdef __DEBUG_PALMLIB
   ! print header
   CALL header
#endif

   ! set start time in format hh:mm:ss
   simulated_time_chr = time_to_string( time_since_reference_point )

   ! integration of the model equations using timestep-scheme
   CALL time_integration

#ifdef __DEBUG_PALMLIB
   ! print header
   CALL header
#endif

   ! after initialization step, reset the PALM T, S, U, V fields to match LSO fields
   DO k = nzb, nzt+1
      pt(k,:,:) = pt(k,:,:) - hom(k,1,4) + pt_init(k)
      sa(k,:,:) = sa(k,:,:) - hom(k,1,5) + sa_init(k)
      u(k,:,:)  = u(k,:,:)  - hom(k,1,1) + u_init(k)
      v(k,:,:)  = v(k,:,:)  - hom(k,1,2) + v_init(k)
      hom(k,1,4) = pt_init(k)
      hom(k,1,5) = sa_init(k)
      hom(k,1,1) = u_init(k)
      hom(k,1,2) = v_init(k)
   ENDDO

   CALL cpu_log( log_point(2), 'initial step', 'stop' )

   END SUBROUTINE palm_init


   SUBROUTINE palm_main(nzLS, zwLS, ptLS, saLS, uLS, vLS,                      &
                        wtflux, wsflux, uwflux, vwflux, wtflux_solar,          &
                        ptLSF, saLSF, uLSF, vLSF,                              &
                        tend_ptLS, tend_saLS, tend_uLS, tend_vLS,              &
                        zSS, ptSS, saSS, uSS, vSS)

   IMPLICIT NONE

   ! number of vertical levels in the large scale model
   INTEGER(iwp), INTENT(in) :: nzLS
   ! surface fluxes
   REAL(wp), INTENT(in) :: wtflux, wsflux, uwflux, vwflux, wtflux_solar
   ! large scale fields (T, S, U, V)
   REAL(wp), DIMENSION(nzLS), INTENT(in) :: ptLS, saLS, uLS, vLS
   ! large scale forcing for the fields (T, S, U, V)
   REAL(wp), DIMENSION(nzLS), INTENT(in) :: ptLSF, saLSF, uLSF, vLSF
   ! large scale z at grid interfaces
   REAL(wp), DIMENSION(nzLS+1), INTENT(in) :: zwLS
   ! tendencies of large scale fields due to small scale processes
   REAL(wp), DIMENSION(nzLS), INTENT(out) :: tend_ptLS, tend_saLS, tend_uLS, tend_vLS
   ! PALM fields
   REAL(wp), DIMENSION(nzLS), INTENT(out) :: zSS, ptSS, saSS, uSS, vSS

   ! local variables
   CHARACTER (LEN=9) ::  time_to_string
   INTEGER :: k
#ifdef __DEBUG_PALMLIB
   REAL(wp), DIMENSION(nzLS)   :: zuLS
#endif

   CALL cpu_log( log_point(3), 'main step', 'start' )

   ! update end time
   end_time = end_time + dt_LS

   ! update surface fluxes
   top_momentumflux_u = uwflux
   top_momentumflux_v = vwflux
   top_heatflux = wtflux
   top_salinityflux = wsflux
   ideal_solar_heatflux = wtflux_solar

   ! large scale forcing on coarse grid
   SELECT CASE (mode_LS)
      CASE (LS_MODE_NONE)
         fldLS = 0.0_wp
      CASE (LS_MODE_INPUT)
         fldLS(1,1,:) = ptLSF
         fldLS(1,2,:) = saLSF
         fldLS(1,3,:) = uLSF
         fldLS(1,4,:) = vLSF
      CASE (LS_MODE_NUDGING)
         fldSS(1,1,:) = hom(nzb+1:nzt,1,4)
         fldSS(1,2,:) = hom(nzb+1:nzt,1,5)
         fldSS(1,3,:) = hom(nzb+1:nzt,1,1)
         fldSS(1,4,:) = hom(nzb+1:nzt,1,2)
         CALL cpu_log( log_point(9), 'interp_fields', 'start' )
         CALL interp_fields_SS2LS(nz, zw(nzb:nzt), fldSS, nzLS, zwLS, fldLS)
         CALL cpu_log( log_point(9), 'interp_fields', 'stop' )
         fldLS(1,1,:) = ptLS - fldLS(1,1,:)
         fldLS(1,2,:) = saLS - fldLS(1,2,:)
         fldLS(1,3,:) = uLS  - fldLS(1,3,:)
         fldLS(1,4,:) = vLS  - fldLS(1,4,:)
      CASE DEFAULT
         WRITE( message_string, * ) 'large scale forcing mode ', mode_LS, 'not suppported'
         CALL message( 'palm_mod', 'PA1001', 1, 2, 0, 6, 0 )
   END SELECT

   ! interpolate large scale forcing from large scale grid to PALM grid
   CALL cpu_log( log_point(9), 'interp_fields', 'start' )
   CALL interp_fields_LS2SS(nzLS, zwLS, fldLS, nz, zw(nzb:nzt), fldSS)
   CALL cpu_log( log_point(9), 'interp_fields', 'stop' )

#ifdef __DEBUG_PALMLIB
   ! z at LS cell centers
   DO k = 1, nzLS
      zuLS(k) = 0.5_wp * (zwLS(k) + zwLS(k+1))
   ENDDO
   write(*,*) '(palm_main) DEBUG: fldLS (large scale forcing)'
   CALL print_fields(nzLS, zuLS, fldLS/dt_LS)
   write(*,*) '(palm_main) DEBUG: fldSS (large scale forcing)'
   CALL print_fields(nz, zu(nzb+1:nzt), fldSS/dt_LS)
#endif

   ! large scale forcing
   pt_LS_forcing = fldSS(1,1,:)
   sa_LS_forcing = fldSS(1,2,:)
   u_LS_forcing  = fldSS(1,3,:)
   v_LS_forcing  = fldSS(1,4,:)

#ifdef __DEBUG_PALMLIB
   ! print header
   CALL header
#endif

   ! set start time in format hh:mm:ss
   simulated_time_chr = time_to_string( time_since_reference_point )

   ! integration of the model equations using timestep-scheme
   CALL time_integration

#ifdef __DEBUG_PALMLIB
   ! print header
   CALL header
#endif

   ! pass to the large scale model the flux divergence, time mean fluxes profiles are saved in hom_sum
   ! note that the average is reset every time the profile is output
   ! so the profile output frequency should be set to the time step of the large scale model
   DO k = nzb+1, nzt
      fldSS(1,1,k) = -1.0 * (hom_sum(k,15) - hom_sum(k-1,15) + hom_sum(k,14) - hom_sum(k-1,14)) * ddzw(k) + hom_sum(k,29)
      fldSS(1,2,k) = -1.0 * (hom_sum(k,17) - hom_sum(k-1,17) + hom_sum(k,16) - hom_sum(k-1,16)) * ddzw(k)
      fldSS(1,3,k) = -1.0 * (hom_sum(k,11) - hom_sum(k-1,11) + hom_sum(k,10) - hom_sum(k-1,10)) * ddzw(k)
      fldSS(1,4,k) = -1.0 * (hom_sum(k,13) - hom_sum(k-1,13) + hom_sum(k,12) - hom_sum(k-1,12)) * ddzw(k)
   ENDDO

   CALL cpu_log( log_point(9), 'interp_fields', 'start' )
   CALL interp_fields_SS2LS(nz, zw(nzb:nzt), fldSS, nzLS, zwLS, fldLS)
   CALL cpu_log( log_point(9), 'interp_fields', 'stop' )

#ifdef __DEBUG_PALMLIB
   write(*,*) '(palm_main) DEBUG: fldSS (tendencies)'
   CALL print_fields(nz, zu(nzb+1:nzt), fldSS)
   write(*,*) '(palm_main) DEBUG: fldLS (tendencies)'
   CALL print_fields(nzLS, zuLS, fldLS)
#endif

   ! profiles passed back
   tend_ptLS = fldLS(1,1,:)
   tend_saLS = fldLS(1,2,:)
   tend_uLS  = fldLS(1,3,:)
   tend_vLS  = fldLS(1,4,:)
   zSS  = zu(nzt:nzb+1:-1)
   ptSS = hom(nzt:nzb+1:-1,1,4)
   saSS = hom(nzt:nzb+1:-1,1,5)
   uSS  = hom(nzt:nzb+1:-1,1,1)
   vSS  = hom(nzt:nzb+1:-1,1,2)

   CALL cpu_log( log_point(3), 'main step', 'stop' )

   END SUBROUTINE palm_main


   SUBROUTINE palm_finalize()

   ! finalize
   CALL ws_finalize
   CALL fft_finalize

   CALL close_file( 0 )
   CALL cpu_log( log_point(1), 'total', 'stop' )
   CALL cpu_statistics

   END SUBROUTINE palm_finalize


   SUBROUTINE interp_fields_LS2SS(nzLS, zwLS, fldLS, nzSS, zwSS, fldSS)
   ! interpolation of fields from the coarse grid of the large scale model to PALM grid
   ! both fldLS and fldSS are at cell center
   ! z indices of fldSS start from the bottom
   ! z indices of fldLS start from the top

   IMPLICIT NONE

   ! input/output
   INTEGER(iwp), INTENT(in) :: nzLS, nzSS
   REAL(wp), DIMENSION(1:nzLS+1), INTENT(in)  :: zwLS
   REAL(wp), DIMENSION(1:nzSS+1), INTENT(in) :: zwSS
   REAL(wp), DIMENSION(ndof, nvar, nzLS), INTENT(in)  :: fldLS
   REAL(wp), DIMENSION(ndof, nvar, nzSS), INTENT(out) :: fldSS

   ! local variables
   INTEGER(iwp) :: il, jl
   REAL(wp), DIMENSION(1:nzSS+1) :: zwSSInv
   REAL(wp), DIMENSION(ndof, nvar, nzSS) :: fldTMP

   !-- this specifies options for the method, here is quartic interp
   opts%edge_meth = p5e_method
   opts%cell_meth = pqm_method
   opts%cell_lims = mono_limit

   bc_l(:)%bcopt = bcon_loose
   bc_r(:)%bcopt = bcon_loose

   ! inverse z in zwSS
   jl = 1
   DO il = nzSS+1, 1, -1
      zwSSInv(jl)  = zwSS(il)
      jl = jl + 1
   ENDDO

   ! initialize PPR
   CALL work%init(nzSS+1, nvar, opts)

   ! remapping
   CALL rmap1d(nzLS+1, nzSS+1, nvar, ndof,   &
               ABS(zwLS(1:nzLS+1)),ABS(zwSSInv(1:nzSS+1)), &
               fldLS, fldTMP, bc_l, bc_r, work, opts)

   ! free PPR
   CALL work%free()

   ! inverse fldTMP to get fldSS that is consistent with zwSS
   jl = 1
   DO il = nzSS, 1, -1
      fldSS(:,:,jl) = fldTMP(:,:,il)
      jl = jl + 1
   ENDDO

   END SUBROUTINE interp_fields_LS2SS


   SUBROUTINE interp_fields_SS2LS(nzSS, zwSS, fldSS, nzLS, zwLS, fldLS)
   ! interpolation of fields from PALM grid to the coarse grid of the large scale model
   ! both fldSS and fldLS are at cell center
   ! z indices of fldSS start from the bottom
   ! z indices of fldLS start from the top

   IMPLICIT NONE

   ! input/output
   INTEGER(iwp), INTENT(in) :: nzSS, nzLS
   REAL(wp), DIMENSION(1:nzSS+1), INTENT(in)  :: zwSS
   REAL(wp), DIMENSION(1:nzLS+1), INTENT(in) :: zwLS
   REAL(wp), DIMENSION(ndof, nvar, nzSS), INTENT(in)  :: fldSS
   REAL(wp), DIMENSION(ndof, nvar, nzLS), INTENT(out) :: fldLS

   ! local variables
   INTEGER(iwp) :: il, jl
   REAL(wp), DIMENSION(1:nzSS+1) :: zwSSInv
   REAL(wp), DIMENSION(ndof, nvar, nzSS) :: fldTMP

   !-- this specifies options for the method, here is quartic interp
   opts%edge_meth = p5e_method
   opts%cell_meth = pqm_method
   opts%cell_lims = mono_limit

   bc_l(:)%bcopt = bcon_loose
   bc_r(:)%bcopt = bcon_loose

   ! inverse z in zwSS
   jl = 1
   DO il = nzSS+1, 1, -1
      zwSSInv(jl)  = zwSS(il)
      jl = jl + 1
   ENDDO
   ! inverse z in fldSS
   jl = 1
   DO il = nzSS, 1, -1
      fldTMP(:,:,jl) = fldSS(:,:,il)
      jl = jl + 1
   ENDDO

   ! initialize PPR
   CALL work%init(nzSS+1, nvar, opts)

   ! remapping
   CALL rmap1d(nzSS+1, nzLS+1, nvar, ndof,   &
               ABS(zwSSInv(1:nzSS+1)),ABS(zwLS(1:nzLS+1)), &
               fldTMP, fldLS, bc_l, bc_r, work, opts)

   ! free PPR
   CALL work%free()

   END SUBROUTINE interp_fields_SS2LS


   SUBROUTINE print_fields(nz, z, fld)

   IMPLICIT NONE

   ! input/output
   INTEGER(iwp), INTENT(in) :: nz
   REAL(wp), DIMENSION(1:nz), INTENT(in) :: z
   REAL(wp), DIMENSION(ndof, nvar, nz) :: fld
   ! local variables
   INTEGER(iwp) :: k

   WRITE( *, 100 )
   DO k = 1, nz
      WRITE( *, 101 ) z(k), fld(1,1,k), fld(1,2,k), fld(1,3,k), fld(1,4,k)
   ENDDO

100 FORMAT ('       Z             T             S             U             V')
101 FORMAT (F8.2,4X,E10.4,4X,E10.4,4X,E10.4,4X,E10.4)

   END SUBROUTINE print_fields

END MODULE palm_mod
