PROGRAM test_palm_lib

! this program tests PALM in the library mode

   USE palm_mod

   USE kinds

   IMPLICIT NONE

   ! flags
   LOGICAL :: diurnal_forcing
   ! Cell ID
   INTEGER(iwp) :: iCell
   ! indices
   INTEGER(iwp) :: it
   ! error number
   INTEGER(iwp) :: ioerr
   ! MPAS number of grid cells
   INTEGER(iwp) :: nzMPAS
   ! PALM number of grid cells
   INTEGER(iwp) :: nzPALM, nxPALM, nyPALM
   ! number of time steps
   INTEGER(iwp) :: nsteps
   ! large scale forcing mode
   INTEGER(iwp) :: ls_mode
   ! PALM grid resolution
   REAL(wp) :: dxPALM, dyPALM
   ! MPAS grid resolution
   REAL(wp) :: dzMPAS
   ! MPAS time step
   REAL(wp) :: dtMPAS
   ! initialization time
   REAL(wp) :: init_time
   ! Coriolis parameter
   REAL(wp) :: fcori
   ! large scale forcing factor
   REAL(wp) :: ls_factor_tracer, ls_factor_vel
   ! surface fluxes
   REAL(wp) :: wtflux, wsflux, uwflux, vwflux
   ! solar flux
   REAL(wp) :: wtflux_solar, solar_frac, solar_depth1, solar_depth2
   ! MPAS fields (T, S, U, V)
   REAL(wp), DIMENSION(:), ALLOCATABLE :: tMPAS, sMPAS, uMPAS, vMPAS, uMPASO, vMPASO
   ! MPAS large scale forcing fields (T, S, U, V)
   REAL(wp), DIMENSION(:), ALLOCATABLE :: tLSF, sLSF, uLSF, vLSF
   ! MPAS z at edges
   REAL(wp), DIMENSION(:), ALLOCATABLE :: zEdgeMPAS, zMidMPAS
   ! Increment in MPAS fields
   REAL(wp), DIMENSION(:), ALLOCATABLE :: tend_tMPAS, tend_sMPAS, tend_uMPAS, tend_vMPAS
   ! PALM fields (T, S, U, V, z)
   REAL(wp), DIMENSION(:), ALLOCATABLE :: tPALM, sPALM, uPALM, vPALM, zPALM
   ! namelist
   NAMELIST /libpar/ nxPALM, nyPALM, nzPALM, dxPALM, dyPALM,    &
                     nzMPAS, dzMPAS, dtMPAS,                    &
                     wtflux, wsflux, uwflux, vwflux,            &
                     diurnal_forcing, wtflux_solar,             &
                     solar_frac, solar_depth1, solar_depth2,    &
                     init_time, ls_mode, ls_factor_tracer,      &
                     ls_factor_vel, fcori, nsteps, iCell

   ! default values
   nzPALM = 64
   nxPALM = 64
   nyPALM = 64
   dxPALM = 2.5_wp
   dyPALM = 2.5_wp
   nzMPAS = 40
   dzMPAS = 1.6_wp
   dtMPAS = 1800.0_wp
   wtflux = 1.19e-5_wp
   wsflux = 0.0_wp
   uwflux = 0.0_wp
   vwflux = 0.0_wp
   diurnal_forcing = .FALSE.
   wtflux_solar = 0.0_wp
   solar_frac = 0.67_wp
   solar_depth1 = 1.0_wp
   solar_depth2 = 17.0_wp
   fcori = 1.028e-4_wp
   init_time = 3600.0_wp
   ls_mode = 0
   ls_factor_tracer = 1.0_wp
   ls_factor_vel = 1.0_wp
   nsteps = 5
   iCell = 1

   ! read namelist
   OPEN( 80, FILE='PALMLIB', STATUS='OLD', FORM='FORMATTED', ACTION='READ', IOSTAT=ioerr )
   IF ( ioerr /= 0 ) THEN
      WRITE(*,* )'(test_palm_lib) ERROR: cannot open namelist PALMLIB'
      STOP
   ENDIF
   READ( 80, NML=libpar, IOSTAT=ioerr )
   IF ( ioerr /= 0 ) THEN
      WRITE(*,* )'(test_palm_mod) ERROR: cannot read namelist PALMLIB'
      STOP
   ENDIF
   CLOSE( 80 )

   ALLOCATE(tMPAS(1:nzMPAS), sMPAS(1:nzMPAS), uMPAS(1:nzMPAS), vMPAS(1:nzMPAS))
   ALLOCATE(uMPASO(1:nzMPAS), vMPASO(1:nzMPAS))
   ALLOCATE(tLSF(1:nzMPAS), sLSF(1:nzMPAS), uLSF(1:nzMPAS), vLSF(1:nzMPAS))
   ALLOCATE(tend_tMPAS(1:nzMPAS), tend_sMPAS(1:nzMPAS))
   ALLOCATE(tend_uMPAS(1:nzMPAS), tend_vMPAS(1:nzMPAS))
   ALLOCATE(tPALM(1:nzPALM), sPALM(1:nzPALM), uPALM(1:nzPALM), vPALM(1:nzPALM), zPALM(1:nzPALM))
   ALLOCATE(zEdgeMPAS(1:nzMPAS+1), zMidMPAS(1:nzMPAS))

   CALL init_z_linear(nzMPAS, dzMPAS, zEdgeMPAS, zMidMPAS)
   CALL idealized_profile_linear(nzMPAS, zMidMPAS, 0.01_wp, 293.0_wp, tMPAS)

   sMPAS = 35.0_wp
   uMPAS = 0.0_wp
   vMPAS = 0.0_wp

#ifdef __DEBUG_PALMLIB
   write(*,*) '(test_palm_mod) DEBUG: initial MPAS profiles'
   CALL print_profiles(nzMPAS, zMidMPAS, tMPAS, sMPAS, uMPAS, vMPAS)
#endif

   ! initialize arrays
   tend_tMPAS = 0.0_wp
   tend_sMPAS = 0.0_wp
   tend_uMPAS = 0.0_wp
   tend_vMPAS = 0.0_wp
   tPALM = 0.0_wp
   sPALM = 0.0_wp
   uPALM = 0.0_wp
   vPALM = 0.0_wp
   zPALM = 0.0_wp

   ! initial step
   CALL palm_init(dtMPAS, nzMPAS, zEdgeMPAS, tMPAS, sMPAS, uMPAS, vMPAS,       &
                  wtflux, wsflux, uwflux, vwflux, wtflux_solar,                &
                  nzPALM, nxPALM, nyPALM, dxPALM, dyPALM,                      &
                  diurnal_forcing, solar_frac, solar_depth1, solar_depth2,     &
                  init_time, ls_mode, ls_factor_tracer, ls_factor_vel,         &
                  fcori, iCell )

#ifdef __DEBUG_PALMLIB
   write(*,*) '(test_palm_mod) DEBUG: MPAS profiles after initial step'
   CALL print_profiles(nzMPAS, zMidMPAS, tMPAS, sMPAS, uMPAS, vMPAS)
#endif

   ! step forward
   DO it = 1, nsteps
      ! large scale forcing
      tLSF = 0.0_wp
      sLSF = 0.0_wp
      uLSF = 0.0_wp
      vLSF = 0.0_wp
      ! step forward for small scale processes
      CALL palm_main(nzMPAS, zEdgeMPAS, tMPAS, sMPAS, uMPAS, vMPAS,  &
                     wtflux, wsflux, uwflux, vwflux, wtflux_solar,   &
                     tLSF, sLSF, uLSF, vLSF,                         &
                     tend_tMPAS, tend_sMPAS, tend_uMPAS, tend_vMPAS, &
                     zPALM, tPALM, sPALM, uPALM, vPALM)

      ! update MPAS fields
      uMPASO = uMPAS
      vMPASO = vMPAS
      tMPAS = tMPAS + tend_tMPAS * dtMPAS
      sMPAS = sMPAS + tend_sMPAS * dtMPAS
      uMPAS = uMPAS + tend_uMPAS * dtMPAS - fcori * vMPASO * dtMPAS
      vMPAS = vMPAS + tend_vMPAS * dtMPAS + fcori * uMPASO * dtMPAS
      ! tMPAS = tPALM
      ! sMPAS = sPALM
      ! uMPAS = uPALM
      ! vMPAS = vPALM
   ENDDO

#ifdef __DEBUG_PALMLIB
   ! output final profiles
   write(*,*) '(test_palm_mod) DEBUG: final MPAS profiles'
   CALL print_profiles(nzMPAS, zMidMPAS, tMPAS, sMPAS, uMPAS, vMPAS)
   write(*,*) '(test_palm_mod) DEBUG: final PALM profiles'
   CALL print_profiles(nzPALM, zPALM, tPALM, sPALM, uPALM, vPALM)
#endif

   CALL palm_finalize

   DEALLOCATE(tMPAS, sMPAS, uMPAS, vMPAS, uMPASO, vMPASO)
   DEALLOCATE(tLSF, sLSF, uLSF, vLSF)
   DEALLOCATE(tend_tMPAS, tend_sMPAS)
   DEALLOCATE(tend_uMPAS, tend_vMPAS)
   DEALLOCATE(zPALM, tPALM, sPALM, uPALM, vPALM)
   DEALLOCATE(zEdgeMPAS, zMidMPAS)

END PROGRAM test_palm_lib

SUBROUTINE init_z_linear(nz, dz, zedge, zmid)

   USE kinds

   IMPLICIT NONE

   ! input/output
   INTEGER(iwp), INTENT(in) :: nz
   REAL(wp), INTENT(in) :: dz
   REAL(wp), DIMENSION(1:nz+1), INTENT(out) :: zedge
   REAL(wp), DIMENSION(1:nz), INTENT(out) :: zmid
   ! local variables
   INTEGER(iwp) :: k

   zedge(1) = 0.0_wp
   zmid(1) = -0.5_wp * dz
   DO k = 2, nz+1
      zedge(k) = zedge(k-1) - dz
   ENDDO
   zmid(2:nz) = 0.5 * ( zedge(2:nz) + zedge(3:nz+1) )

END SUBROUTINE init_z_linear

SUBROUTINE idealized_profile_linear(nz, z, dpfldz, pfl0, pfl)

   USE kinds

   IMPLICIT NONE

   ! input/output
   INTEGER(iwp), INTENT(in) :: nz
   REAL(wp), INTENT(in) :: dpfldz, pfl0
   REAL(wp), DIMENSION(1:nz), INTENT(in) :: z
   REAL(wp), DIMENSION(1:nz), INTENT(out) :: pfl
   ! local variables
   INTEGER(iwp) :: k

   pfl(1) = pfl0 + 0.5_wp * z(1) * dpfldz
   DO k = 2, nz
      pfl(k) = pfl(k-1) + ( z(k) - z(k-1) ) * dpfldz
   ENDDO

END SUBROUTINE idealized_profile_linear

SUBROUTINE print_profiles(nz, z, tpfl, spfl, upfl, vpfl)

   USE kinds

   IMPLICIT NONE

   ! input/output
   INTEGER(iwp), INTENT(in) :: nz
   REAL(wp), DIMENSION(1:nz), INTENT(in) :: z
   REAL(wp), DIMENSION(1:nz), INTENT(in) :: tpfl, spfl, upfl, vpfl
   ! local variables
   INTEGER(iwp) :: k

   WRITE( *, 100 )
   DO k = 1, nz
      WRITE( *, 101 ) z(k), tpfl(k), spfl(k), upfl(k), vpfl(k)
   ENDDO

100 FORMAT ('       Z             T             S             U             V')
101 FORMAT (F8.2,4X,F10.4,4X,F10.4,4X,F10.4,4X,F10.4)

END SUBROUTINE print_profiles
