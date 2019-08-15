!> @file palm.f90
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
! Description:
! ------------
!> Large-Eddy Simulation (LES) model for the convective boundary layer,
!> optimized for use on parallel machines (implementation realized using the
!> Message Passing Interface (MPI)). The model can also be run on vector machines
!> (less well optimized) and workstations. Versions for the different types of
!> machines are controlled via cpp-directives.
!> Model runs are only feasible using the ksh-script mrun.
!>
!> @todo create routine last_actions instead of calling lsm_last_actions etc.
!> @todo move chem_init call to init_3d_model or to check_parameters
!------------------------------------------------------------------------------!
module palm_mod

    USE advec_ws,                                                              &
            ONLY:  ws_statistics, ws_finalize

    USE arrays_3d

    USE control_parameters

    USE configure_3D_MODEL

    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point, log_point_s, cpu_statistics

    USE indices

!    USE netcdf_data_input_mod,                                                 &
!        ONLY:  netcdf_data_input_inquire_file, netcdf_data_input_init,         &
!               netcdf_data_input_surface_data, netcdf_data_input_topo

    USE fft_xy,                                                                 &
        ONLY: fft_finalize

    USE kinds

    USE make_vertical_grid

    USE ppr_1d

    USE pegrid

    USE random_generator_parallel, ONLY: deallocate_random_generator

    use surface_mod

    use tridia_solver, ONLY: tridia_deallocate

    use turbulence_closure_mod, ONLY: tcm_deallocate_arrays

    USE write_restart_data_mod,                                                &
        ONLY:  wrd_global, wrd_local

    use statistics

#if defined( __cudaProfiler )
    USE cudafor
#endif

    IMPLICIT NONE

    Real(wp),allocatable,dimension(:)   :: T_mpas2, S_mpas2, U_mpas2, V_mpas2
    Real(wp),allocatable,dimension(:)   :: Tles, Sles, Ules, Vles, zmid, zedge
    Real(wp),allocatable,dimension(:)   :: zeLES, wtLES, wsLES, wuLES, wvLES
    Real(wp),allocatable,dimension(:)   :: zeLESInv

!-- arrays and parameters for PPR remapping

    integer, parameter :: nvar = 4
    integer, parameter :: ndof = 1
    real(wp),allocatable :: fMPAS(:,:,:)
    type(rmap_work) :: work
    type(rmap_opts) :: opts
    type(rcon_ends) :: bc_l(nvar)
    type(rcon_ends) :: bc_r(nvar)

!-- Local variables
    CHARACTER(LEN=9)  ::  time_to_string  !<
    CHARACTER(LEN=10) ::  env_string      !< to store string of environment var
    INTEGER(iwp)      ::  env_stat        !< to hold status of GET_ENV
    INTEGER(iwp)      ::  myid_openmpi    !< OpenMPI local rank for CUDA aware MPI
    Real(wp) :: coeff1, coeff2

    public :: palm_init, palm_main, palm_finalize

    contains

    subroutine palm_init(nCells,nVertLevels,T_mpas,S_mpas,U_mpas,V_mpas,lt_mpas, &
                lat_mpas,maxLevels,wtflux,wtflux_solar, wsflux,uwflux, &
             vwflux,fac,dep1,dep2,dxLES, dyLES, dzLES,nxLES,nyLES,nzLES,        &
             endTime, dtDisturb, tIncrementLES,sIncrementLES,             &
             uIncrementLES,vIncrementLES,tempLES,    &
             salinityLES, uLESout, vLESout, dtLS, zLES, &
             disturbMax, disturbAmp, disturbBot, disturbTop, disturbNblocks, &
             botDepth, timeAv, constant_dz, MLD, mixed_layer_refine, restore_strength)

!
! -- Variables from MPAS
   logical,intent(in) :: mixed_layer_refine, constant_dz
   integer(iwp),dimension(nCells) :: maxLevels
   integer(iwp) :: iCell,nCells, nVertLevels, il, jl, jloc, kl, knt, nzLES, iz, disturbNblocks
   integer(iwp) :: nzMPAS, zmMPASspot, zeMPASspot, nxLES, nyLES, idx_cutoff
   Real(wp),intent(in)                             :: dtLS
   Real(wp),dimension(nVertLevels,nCells),intent(inout)   :: T_mpas, S_mpas, U_mpas, V_mpas
   Real(wp),dimension(nVertLevels,nCells),intent(inout)   :: tIncrementLES, sIncrementLES, &
                                                        uIncrementLES, vIncrementLES
   Real(wp),dimension(nzLES,nCells),intent(out)           :: tempLES, salinityLES, zLES
   Real(wp),dimension(nzLES,nCells),intent(out)           :: uLESout, vLESout
   Real(wp),dimension(nVertLevels,nCells),intent(in)      :: lt_mpas
   Real(wp),dimension(nCells) :: wtflux, wsflux, uwflux, vwflux, disturbBot
   Real(wp),dimension(nCells) :: botDepth,wtflux_solar, lat_mpas, MLD
   Real(wp) :: dxLES, dyLES, dzLES, z_fac, z_frst, z_cntr, restore_strength
   real(wp) :: z_fac1, z_fac2, z_facn, tol, test, fac, dep1, dep2
   real(wp) :: dtDisturb, endTime, thickDiff, disturbMax, disturbAmp
   real(wp) :: disturbTop, timeAv
   real(wp) :: sumValT, sumValS, sumValU, sumValV, thickVal
   real(wp) :: fLES(ndof, nvar, nzLES+1)
   real(wp) :: fMPAS(ndof, nvar, nVertLevels)
   CHARACTER(LEN=128) :: format
!-- this specifies options for the method, here is quartic interp
   opts%edge_meth = p5e_method
   opts%cell_meth = pqm_method
   opts%cell_lims = mono_limit

   bc_l(:)%bcopt = bcon_loose
   bc_r(:)%bcopt = bcon_loose

   ! call init_control_parameters

   disturbFactor = restore_strength
   dt_disturb = dtDisturb
   end_time = 3600.0_wp
   ideal_solar_division = fac
   ideal_solar_efolding1 = dep1
   ideal_solar_efolding2 = dep2
   nz = nzLES
   nx = nxLES
   ny = nyLES
   dx = dxLES
   dy = dyLES
   disturb_nblocks = disturbNblocks
   ! dt_ls is the total run time for each MPAS time step, used in time_integration()
   ! simulated_time need to be reset to 0 for each iCell loop
   dt_ls = dtLS
   dt_avg = timeAv

   disturbance_level_t = disturbTop
   disturbance_amplitude = disturbAmp
   disturbance_energy_limit = disturbMax

   allocate(zmid(nVertLevels),zedge(nVertLevels+1))
   allocate(T_mpas2(nVertLevels),S_mpas2(nVertLevels),U_mpas2(nVertLevels))
   allocate(V_mpas2(nVertLevels))

   iCell = 1

   ! vertical grid for MPAS (1:nzVertLeves)
   zmid(1) = -0.5_wp*lt_mpas(1,iCell)
   zedge(1) = 0.0_wp

   do il=2,nVertLevels
      zmid(il) = zmid(il-1) - 0.5*(lt_mpas(il-1,iCell) + lt_mpas(il,iCell))
      zedge(il) = zedge(il-1) - lt_mpas(il-1,iCell)
   enddo

   zedge(nvertLevels+1) = zedge(nVertLevels) - lt_mpas(nVertLevels,iCell)

   ! find the index of bottom
   do il=1,nVertLevels
     if(zmid(il) < minval(botDepth)) then
       zmMPASspot = il
       nzMPAS = il
       exit
     endif
   enddo

   do il=1,nVertLevels
     if(zedge(il) < minval(botDepth)) then
       zeMPASspot = il
       exit
     endif
   enddo

#if defined( __parallel )
!
!-- MPI initialisation. comm2d is preliminary set, because
!-- it will be defined in init_pegrid but is used before in cpu_log.
    CALL MPI_INIT( ierr )

    comm_palm = MPI_COMM_WORLD
    comm2d = comm_palm
!
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
!
#endif

!-- Initialize measuring of the CPU-time remaining to the run
    CALL local_tremain_ini
!
!-- Start of total CPU time measuring.
    CALL cpu_log( log_point(1), 'total', 'start' )
    CALL cpu_log( log_point(2), 'initialisation', 'start' )

    !
!-- Read control parameters from NAMELIST files and read environment-variables
    ! CALL parin
    call init_control_parameters

!-- Determine processor topology and local array indices
    CALL init_pegrid
    allocate(zu(nzb:nzt+1),zw(nzb:nzt+1))
    allocate(Tles(0:nz+1),Sles(0:nz+1),Ules(0:nz+1),Vles(0:nz+1))
    allocate(zeLES(nzb-1:nzt+1), zeLESinv(1:nz+2))
    ALLOCATE( hyp(nzb:nzt+1) )

    if (constant_dz) then
      call construct_vertical_grid_const(dzLES,nzb,nzt,zu,zw)
    else
      if (mixed_layer_refine) then
        call construct_vertical_grid_variable(zedge(zeMPASspot),zeLES,nCells,dzLES,nzLES,nzb,nzt,zw,zu,MLD)
      else
        call construct_vertical_grid_variable(zedge(zeMPASspot),zeLES,nCells,dzLES,nzLES,nzb,nzt,zw,zu)
      endif
    endif

    call work%init(nz+1,nvar,opts)
   !
!-- Generate grid parameters, initialize generic topography and further process
!-- topography information if required
    CALL init_grid
    ALLOCATE( pt_init(0:nz+1), q_init(0:nz+1), s_init(0:nz+1),        &
                       ref_state(0:nz+1), sa_init(0:nz+1), ug(0:nz+1),  &
                       u_init(0:nz+1), v_init(0:nz+1), vg(0:nz+1),       &
                       hom(0:nz+1,2,14,0:statistic_regions), &
                       hom_sum(0:nz+1,14,0:statistic_regions))!, &
    ! if( ALLOCATED(meanFields_avg)) then
    !    deallocate(meanFields_avg)
    ! endif
    allocate(meanFields_avg(0:nz+1,4))
    meanFields_avg = 0.0_wp

!-- Check control parameters and deduce further quantities
    CALL check_parameters
! interpolate mpas data to les and send to init variable
    CALL allocate_3d_arrays(nCells)
!-- Initialize all necessary variables
    CALL init_3d_model  ! need a pass through for restarts

    do iCell=1,nCells
    ! do iCell=1,1

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
#endif


    zmid(1) = -0.5_wp*lt_mpas(1,iCell)
    zedge(1) = 0.0_wp

    do il=2,maxLevels(iCell)
       zmid(il) = zmid(il-1) - 0.5*(lt_mpas(il-1,iCell) + lt_mpas(il,iCell))
       zedge(il) = zedge(il-1) - lt_mpas(il-1,iCell)
    enddo

    zedge(nvertLevels+1) = zedge(nVertLevels) - lt_mpas(maxLevels(iCell),iCell)

   do il=1,maxLevels(iCell)
     if(zmid(il) < botDepth(iCell)) then
       zmMPASspot = il-1
       nzMPAS = il-1
       exit
     endif
   enddo

   do il=1,nVertLevels
     if(zedge(il) < botDepth(iCell)) then
       zeMPASspot = il-1
       exit
     endif
   enddo

    if (constant_dz) then
      call construct_vertical_grid_const(dzLES,nzb,nzt,zu,zw)
    else
      if (mixed_layer_refine) then
        call construct_vertical_grid_variable(zedge(zeMPASspot),zeLES,nCells,dzLES,nzLES,nzb,nzt,zw,zu,MLD)
      else
        call construct_vertical_grid_variable(zedge(zeMPASspot),zeLES,nCells,dzLES,nzLES,nzb,nzt,zw,zu)
      endif
    endif


!--    In case of dirichlet bc for u and v the first u- and w-level are defined
!--    at same height.
    IF ( ibc_uv_b == 0 ) THEN
       zu(0) = zw(0)
    ENDIF

    if (iCell == 1) then
       print *, zu
       print *, zw
    endif
!
!-- Compute grid lengths.
    DO  k = 1, nzt+1
       dzu(k)  = zu(k) - zu(k-1)
       ddzu(k) = 1.0_wp / dzu(k)
       dzw(k)  = zw(k) - zw(k-1)
       ddzw(k) = 1.0_wp / dzw(k)
    ENDDO

    DO  k = 1, nzt
       dd2zu(k) = 1.0_wp / ( dzu(k) + dzu(k+1) )
    ENDDO

 !TODO add check for right / acceptable range.
    disturbance_level_b = disturbBot(iCell)
    top_momentumflux_u = uwflux(iCell)
    top_momentumflux_v = vwflux(iCell)
    top_heatflux =  wtflux(iCell)
    top_salinityflux = -wsflux(iCell)
    latitude = lat_mpas(iCell) * 180.0 / pi
    wb_solar = wtflux_solar(iCell)

    f  = 2.0_wp * omega * SIN( latitude / 180.0_wp * pi )
    fs = 0.0_wp * omega * COS( latitude / 180.0_wp * pi )

    if (iCell == 1) then
    print *, 'top momentum flux u: ', top_momentumflux_u
    print *, 'top momentum flux v: ', top_momentumflux_v
    print *, 'top heat flux: ', top_heatflux
    print *, 'top salinityflux: ', top_salinityflux
    print *, 'solar radiation: ', wb_solar
    endif

    fMPAS(1,1,:nzMPAS) = T_mpas(1:nzMPAS,iCell)
    fMPAS(1,2,:nzMPAS) = S_mpas(1:nzMPAS,iCell)
    fMPAS(1,3,:nzMPAS) = U_mpas(1:nzMPAS,iCell)
    fMPAS(1,4,:nzMPAS) = V_mpas(1:nzMPAS,iCell)

    jl=1
    do il = nzt,nzb,-1
      zeLESinv(jl) = zw(il)
      jl = jl + 1
    enddo
    zeLESinv(nzt+2) = zedge(nzMPAS+1)

    call ws_init2
    call rmap1d(nzMPAS+1,nz+2,nvar,ndof,abs(zedge(1:nzMPAS+1)),abs(zeLESinv(1:nz+2)), &
                fMPAS(:,:,:nzMPAS), fLES, bc_l, bc_r, work, opts)

    if (iCell == 1) then
    format = "(6x, F10.3, F10.3, F10.3, F10.3, F10.3)"
    print *, 'fMPAS         Z         T         S         U         V'
    do il = 1,nzMPAS
      write(*,format) 0.5*(zedge(il)+zedge(il+1)), fMPAS(1,1,il), fMPAS(1,2,il), fMPAS(1,3,il), fMPAS(1,4,il)
    enddo
    print *, ' fLES         Z         T         S         U         V'
    do il = 1,nz
      write(*,format) 0.5*(zeLESinv(il)+zeLESinv(il+1)), fLES(1,1,il), fLES(1,2,il), fLES(1,3,il), fLES(1,4,il)
    enddo
    endif


    jl = 1
    do il = nzt,nzb+1,-1
      tLSforcing(il) = fLES(1,1,jl) + 273.15
      sLSforcing(il) = fLES(1,2,jl)
      uLSforcing(il) = fLES(1,3,jl)
      vLSforcing(il) = fLES(1,4,jl)
      jl = jl + 1
    enddo

    do jl=nzt,nzb+1,-1
       pt(jl,:,:) = tLSforcing(jl)
       sa(jl,:,:) = sLSforcing(jl)
       u(jl,:,:) = uLSforcing(jl)
       v(jl,:,:) = vLSforcing(jl)
       zLES(jl,iCell) = zeLES(jl)
       tempLES(jl,iCell) = tLSforcing(jl)
       salinityLES(jl,iCell) = sLSforcing(jl)
       uLESout(jl,iCell) = uLSforcing(jl)
       vLESout(jl,iCell) = vLSforcing(jl)
       uProfileInit(jl) = uLSforcing(jl)
       vProfileInit(jl) = vLSforcing(jl)
       tProfileInit(jl) = tLSforcing(jl)
       sProfileInit(jl) = sLSforcing(jl)
    enddo

    ! no restoring for init step
    uLSforcing = 0.0_wp
    vLSforcing = 0.0_wp
    tLSforcing = 0.0_wp
    sLSforcing = 0.0_wp

    pt(nzb,:,:) = pt(nzb+1,:,:)
    sa(nzb,:,:) = sa(nzb+1,:,:)
    u(nzb,:,:) = u(nzb+1,:,:)
    v(nzb,:,:) = v(nzb+1,:,:)

    pt(nzt+1,:,:) = pt(nzt,:,:)
    sa(nzt+1,:,:) = sa(nzt,:,:)
    u(nzt+1,:,:) = u(nzt,:,:)
    v(nzt+1,:,:) = v(nzt,:,:)

    CALL init_3d_model
pt_p = pt
sa_p = sa
u_p = u
v_p = v
    CALL cpu_log( log_point(2), 'initialisation', 'stop' )
!
!-- Set start time in format hh:mm:ss
!    simulated_time_chr = time_to_string( time_since_reference_point )

! TODO: for testing in single column mode <20190724, Qing Li> !
    if (iCell == 1) then

#if defined( __cudaProfiler )
!-- Only profile time_integration
    CALL cudaProfilerStart()
#endif

!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

#if defined( __cudaProfiler )
!-- Only profile time_integration
    CALL cudaProfilerStop()
#endif
!
!-- Take final CPU-time for CPU-time analysis
    CALL cpu_log( log_point(1), 'total', 'stop' )
!    CALL cpu_statistics

   ! TODO: is there a way to merge this with the code in time_integration? <20190724, Qing Li> !
    if(average_count_meanpr /= 0) then

       meanFields_avg(:,1) = meanFields_avg(:,1) / average_count_meanpr
       meanFields_avg(:,2) = meanFields_avg(:,2) / average_count_meanpr
       meanFields_avg(:,3) = meanFields_avg(:,3) / average_count_meanpr
       meanFields_avg(:,4) = meanFields_avg(:,4) / average_count_meanpr

    endif

    endif ! iCell

    Tles = meanFields_avg(:,3)
    Sles = meanFields_avg(:,4)
    Ules = meanFields_avg(:,1)
    Vles = meanFields_avg(:,2)
 ! need to integrate over layers in mpas to get increments

 if(minval(tempLES(:,iCell)) < 100.0_wp) tempLES(:,iCell) = tempLES(:,iCell) + 273.15_wp
    tProfileInit(1:) = tempLES(:,iCell)
    sProfileInit(1:) = salinityLES(:,iCell)
    uProfileInit(1:) = uLESout(:,iCell)
    vProfileInit(1:) = vLESout(:,iCell)

    jl=1
    do il=nzt,nzb+1,-1
      fLES(1,1,jl) = (Tles(il) - tProfileInit(il)) / dtLS
      fLES(1,2,jl) = (Sles(il) - sProfileInit(il)) / dtLS
      fLES(1,3,jl) = (Ules(il) - uProfileInit(il)) / dtLS
      fLES(1,4,jl) = (Vles(il) - vProfileInit(il)) / dtLS
      jl = jl+1
    enddo
    fLES(:,:,nz+1) = fLES(:,:,nz)

    call rmap1d(nz+2,nzMPAS+1,nvar,ndof,abs(zeLESinv(1:nz+2)),abs(zedge(1:nzMPAS+1)), &
                fLES,fMPAS(:,:,:nzMPAS),bc_l,bc_r,work,opts)

    idx_cutoff = nzMPAS
    do il=1,nVertLevels
     if(zmid(il) < zw(0)) then
       idx_cutoff = il
       exit
     endif
    enddo
    fMPAS(:,:,idx_cutoff:nzMPAS) = 0.0_wp

    if (iCell == 1) then
    format = "(6x, F14.3, E14.3, E14.3, E14.3, E14.3)"
    print *, ' fLES             Z          dTdt          dSdt          dUdt          dVdt'
    do il = 1,nzLES
      write(*,format) 0.5*(zeLESinv(il)+zeLESinv(il+1)), fLES(1,1,il), fLES(1,2,il), fLES(1,3,il), fLES(1,4,il)
    enddo
    print *, 'fMPAS             Z          dTdt          dSdt          dUdt          dVdt'
    do il = 1,nzMPAS
      write(*,format) 0.5*(zedge(il)+zedge(il+1)), fMPAS(1,1,il), fMPAS(1,2,il), fMPAS(1,3,il), fMPAS(1,4,il)
    enddo
    endif

    tIncrementLES(:,iCell) = 0.0_wp
    sIncrementLES(:,iCell) = 0.0_wp
    uIncrementLES(:,iCell) = 0.0_wp
    vIncrementLES(:,iCell) = 0.0_wp
    do jl=1,nzMPAS
      tIncrementLES(jl,iCell) = fMPAS(1,1,jl)
      sIncrementLES(jl,iCell) = fMPAS(1,2,jl)
      uIncrementLES(jl,iCell) = fMPAS(1,3,jl)
      vIncrementLES(jl,iCell) = fMPAS(1,4,jl)
    enddo

    !put variables in arrays to continue
    u_restart(:,:,:,iCell) = u(:,:,:)
    v_restart(:,:,:,iCell) = v(:,:,:)
    w_restart(:,:,:,iCell) = w(:,:,:)
    pt_restart(:,:,:,iCell) = pt(:,:,:)
    sa_restart(:,:,:,iCell) = sa(:,:,:)
    e_restart(:,:,:,iCell) = e(:,:,:)
    km_restart(:,:,:,iCell) = km(:,:,:)
    kh_restart(:,:,:,iCell) = kh(:,:,:)
    u_mean_restart(:,iCell) = Ules
    v_mean_restart(:,iCell) = Vles
    t_mean_restart(:,iCell) = Tles
    s_mean_restart(:,iCell) = Sles
    call init_control_parameters
  enddo !ends icell loop

!   deallocate(zmid,zedge)
!   deallocate(T_mpas2,S_mpas2,U_mpas2)
!   deallocate(V_mpas2)
  ! do iCell=2,nCells
  !   u_restart(:,:,:,iCell) = u(:,:,:)
  !   v_restart(:,:,:,iCell) = v(:,:,:)
  !   w_restart(:,:,:,iCell) = w(:,:,:)
  !   pt_restart(:,:,:,iCell) = pt(:,:,:)
  !   sa_restart(:,:,:,iCell) = sa(:,:,:)
  !   e_restart(:,:,:,iCell) = e(:,:,:)
  !   km_restart(:,:,:,iCell) = km(:,:,:)
  !   kh_restart(:,:,:,iCell) = kh(:,:,:)

  !   u_mean_restart(:,iCell) = u_mean_restart(:,1)
  !   v_mean_restart(:,iCell) = v_mean_restart(:,1)
  !   t_mean_restart(:,iCell) = t_mean_restart(:,1)
  !   s_mean_restart(:,iCell) = s_mean_restart(:,1)
  !   tIncrementLES(:,iCell) = tIncrementLES(:,1)
  !   sIncrementLES(:,iCell) = sIncrementLES(:,1)
  !   uIncrementLES(:,iCell) = uIncrementLES(:,1)
  !   vIncrementLES(:,iCell) = vIncrementLES(:,1)
  !   tempLES(:,iCell) = tempLES(:,1)
  !   salinityLES(:,iCell) = salinityLES(:,1)
  !   uLESout(:,iCell) = uLESout(:,1)
  !   vLESout(:,iCell) = vLESout(:,1)
  !   zLES(:,iCell) = zLES(:,1)
  ! enddo

end subroutine palm_init

subroutine palm_main(nCells,nVertLevels,T_mpas,S_mpas,U_mpas,V_mpas,lt_mpas, &
             lat_mpas,maxLevels,wtflux,wtflux_solar, wsflux,uwflux, &
              vwflux,fac,dep1,dep2,dxLES,dyLES,dzLES,nxLES,nyLES,nzLES,        &
             endTime, dtDisturb, tIncrementLES,sIncrementLES,             &
             uIncrementLES,vIncrementLES,tempLES,    &
             salinityLES, uLESout, vLESout, dtLS, zLES, &
             disturbMax, disturbAmp, disturbBot, disturbTop, disturbNblocks, &
             botDepth, timeAv, constant_dz, MLD, mixed_layer_refine, restore_strength)
!
! -- Variables from MPAS
   logical,intent(in) :: mixed_layer_refine, constant_dz
   integer(iwp),dimension(nCells) :: maxLevels
   integer(iwp) :: iCell,nCells, nVertLevels, il, jl, jloc, kl, knt, nzLES, iz, disturbNblocks
   integer(iwp) :: nxLES, nyLES, nzMPAS, zmMPASspot, zeMPASspot, idx_cutoff
   Real(wp),intent(in)                             :: dtLS
   Real(wp),dimension(nVertLevels,nCells),intent(inout)   :: T_mpas, S_mpas, U_mpas, V_mpas
   Real(wp),dimension(nVertLevels,nCells),intent(inout)   :: tIncrementLES, sIncrementLES, &
                                                        uIncrementLES, vIncrementLES
   Real(wp),dimension(nzLES,nCells),intent(out)           :: tempLES, salinityLES, zLES
   Real(wp),dimension(nzLES,nCells),intent(out)           :: uLESout, vLESout
   Real(wp),dimension(nVertLevels,nCells),intent(in)      :: lt_mpas
   Real(wp),dimension(nCells) :: wtflux, wsflux, uwflux, vwflux, disturbBot
   Real(wp),dimension(nCells) :: botDepth, wtflux_solar, lat_mpas, MLD
   Real(wp) :: dxLES, dyLES, dzLES, z_fac, z_frst, z_cntr, restore_strength
   real(wp) :: z_fac1, z_fac2, z_facn, tol, test, fac, dep1, dep2
   real(wp) :: dtDisturb, endTime, thickDiff, disturbMax, disturbAmp
   real(wp) :: disturbTop, timeAv
   real(wp) :: sumValT, sumValS, sumValU, sumValV, thickVal
   real(wp) :: fLES(ndof, nvar, nzLES+1)
   real(wp) :: fMPAS(ndof, nvar, nVertLevels)
   CHARACTER(LEN=128) :: format
!-- this specifies options for the method, here is quartic interp
   opts%edge_meth = p5e_method
   opts%cell_meth = pqm_method
   opts%cell_lims = mono_limit

   bc_l(:)%bcopt = bcon_loose
   bc_r(:)%bcopt = bcon_loose

   call init_control_parameters

   create_disturbances=.false.
   dt_disturb = 0.0_wp
   end_time = endTime
   ideal_solar_division = fac
   ideal_solar_efolding1 = dep1
   ideal_solar_efolding2 = dep2
   dx = dxLES
   dy = dyLES
   nz = nzLES
   nx = nxLES
   ny = nyLES
   disturb_nblocks = disturbNblocks
   dt_ls = dtLS
   dt_avg = timeAv

   disturbFactor = restore_strength
   disturbance_level_t = disturbTop
   disturbance_amplitude = disturbAmp
   disturbance_energy_limit = disturbMax

   do iCell=1,nCells
   ! do iCell=1,1

   initializing_actions  = 'SP_run_continue'
   zmid(1) = -0.5_wp*lt_mpas(1,iCell)
   zedge(1) = 0
   do il=2,maxLevels(iCell)
      zmid(il) = zmid(il-1) - 0.5*(lt_mpas(il-1,iCell) + lt_mpas(il,iCell))
      zedge(il) = zedge(il-1) - lt_mpas(il-1,iCell)
   enddo

   zedge(nvertLevels+1) = zedge(nVertLevels) - lt_mpas(maxLevels(iCell),iCell)

   do il=1,maxLevels(iCell)
     if(zmid(il) < botDepth(iCell)) then
       zmMPASspot = il-1
       nzMPAS = il-1
       exit
     endif
   enddo

   do il=1,nVertLevels
     if(zedge(il) < botDepth(iCell)) then
       zeMPASspot = il-1
       exit
     endif
   enddo
    nzt = nzLES
    if (constant_dz) then
      call construct_vertical_grid_const(dzLES,nzb,nzt,zu,zw)
      ! call construct_vertical_grid_const(dzLES,nzLES,zeLES,zedge(zeMPASspot),nzt,nzb,zu,zw)
    else
      if (mixed_layer_refine) then
        call construct_vertical_grid_variable(zedge(zeMPASspot),zeLES,nCells,dzLES,nzLES,nzb,nzt,zw,zu,disturbBot)
      else
        call construct_vertical_grid_variable(zedge(zeMPASspot),zeLES,nCells,dzLES,nzLES,nzb,nzt,zw,zu)
      endif
    endif

!--    In case of dirichlet bc for u and v the first u- and w-level are defined
!--    at same height.
       IF ( ibc_uv_b == 0 ) THEN
          zu(0) = zw(0)
       ENDIF
!
!-- Compute grid lengths.
    DO  k = 1, nzt+1
       dzu(k)  = zu(k) - zu(k-1)
       ddzu(k) = 1.0_wp / dzu(k)
       dzw(k)  = zw(k) - zw(k-1)
       ddzw(k) = 1.0_wp / dzw(k)
    ENDDO

    DO  k = 1, nzt
       dd2zu(k) = 1.0_wp / ( dzu(k) + dzu(k+1) )
    ENDDO


!TODO add check for right / acceptable range.
    disturbance_level_b = disturbBot(iCell)
    top_momentumflux_u = uwflux(iCell)
    top_momentumflux_v = vwflux(iCell)
    top_heatflux =  wtflux(iCell)
    top_salinityflux = -wsflux(iCell)
    latitude = lat_mpas(iCell) * 180.0 / pi
    wb_solar = wtflux_solar(iCell)

    f  = 2.0_wp * omega * SIN( latitude / 180.0_wp * pi )
    fs = 0.0_wp * omega * COS( latitude / 180.0_wp * pi )

    fMPAS(1,1,:nzMPAS) = T_mpas(1:nzMPAS,iCell)
    fMPAS(1,2,:nzMPAS) = S_mpas(1:nzMPAS,iCell)
    fMPAS(1,3,:nzMPAS) = U_mpas(1:nzMPAS,iCell)
    fMPAS(1,4,:nzMPAS) = V_mpas(1:nzMPAS,iCell)

    jl=1
    do il = nzt,nzb,-1
      zeLESinv(jl) = zw(il)
      jl = jl + 1
    enddo
    zeLESinv(nzt+2) = zedge(nzMPAS+1)

    call rmap1d(nzMPAS+1,nz+2,nvar,ndof,abs(zedge(1:nzMPAS+1)),abs(zeLESinv(1:nzLES+2)), &
                fMPAS(:,:,:nzMPAS), fLES, bc_l, bc_r, work, opts)

    if (iCell == 1) then
    format = "(6x, F10.3, F10.3, F10.3, F10.3, F10.3)"
    print *, 'fMPAS         Z         T         S         U         V'
    do il = 1,nzMPAS
      write(*,format) 0.5*(zedge(il)+zedge(il+1)), fMPAS(1,1,il), fMPAS(1,2,il), fMPAS(1,3,il), fMPAS(1,4,il)
    enddo
    print *, ' fLES         Z         T         S         U         V'
    do il = 1,nzLES
      write(*,format) 0.5*(zeLESinv(il)+zeLESinv(il+1)), fLES(1,1,il), fLES(1,2,il), fLES(1,3,il), fLES(1,4,il)
    enddo
    endif

    jl = 1
    do il = nzt,nzb+1,-1
      ! tLSforcing(il) = fLES(1,1,jl) + 273.15
      ! sLSforcing(il) = fLES(1,2,jl)
      ! uLSforcing(il) = fLES(1,3,jl)
      ! vLSforcing(il) = fLES(1,4,jl)
      tLSforcing(il) = fLES(1,1,jl) + 273.15 - t_mean_restart(il,iCell)
      sLSforcing(il) = fLES(1,2,jl) - s_mean_restart(il,iCell)
      uLSforcing(il) = fLES(1,3,jl) - u_mean_restart(il,iCell)
      vLSforcing(il) = fLES(1,4,jl) - v_mean_restart(il,iCell)
      ! uLSforcing(il) = 0.0_wp
      ! vLSforcing(il) = 0.0_wp
      jl = jl + 1
    enddo

    if (iCell == 1) then
       ! format = "(5x, F10.3, F10.3)"
       ! print *, "LS forcing       max       min"
       ! write(*, format) maxval(tLSforcing), minval(tLSforcing)
       ! write(*, format) maxval(sLSforcing), minval(sLSforcing)
       ! write(*, format) maxval(uLSforcing), minval(uLSforcing)
       ! write(*, format) maxval(vLSforcing), minval(vLSforcing)
       format = "(6x, F10.3, F10.3, F10.3, F10.3, F10.3)"
       print *, '  LSF         Z         T         S         U         V'
       do il = nzt,nzb+1,-1
          write(*, format) 0.5*(zw(il)+zw(il-1)), tLSforcing(il), sLSforcing(il), uLSforcing(il), vLSforcing(il)
       enddo
    endif

    ! do jl = nzt,nzb+1,-1
! !          pt(jl,:,:) = tempLES(jl) + 273.15_wp
! !          sa(jl,:,:) = salinityLES(jl)
! !          u(jl,:,:) = uLESout(jl)
! !          v(jl,:,:) = vLESout(jl)
    !    uProfileInit(jl) = uLSforcing(jl)
    !    vProfileInit(jl) = vLSforcing(jl)
    !    tProfileInit(jl) = tLSforcing(jl)
    !    sProfileInit(jl) = sLSforcing(jl)
    ! enddo
!need to cudify this super easy collapse 3 herrel
    u(:,:,:) = u_restart(:,:,:,iCell)
    v(:,:,:) = v_restart(:,:,:,iCell)
    w(:,:,:) = w_restart(:,:,:,iCell)
    pt(:,:,:) = pt_restart(:,:,:,iCell)
    sa(:,:,:) = sa_restart(:,:,:,iCell)
    e(:,:,:) = e_restart(:,:,:,iCell)
    kh(:,:,:) = kh_restart(:,:,:,iCell)
    km(:,:,:) = km_restart(:,:,:,iCell)

    u_p = u
    v_p = v
    w_p = w
    pt_p = pt
    sa_p = sa
    e_p = e

    CALL init_3d_model
    ! CALL flow_statistics
    ! uLESout(:,iCell) = hom(:,1,1,0)
    ! vLESout(:,iCell) = hom(:,1,2,0)
    ! tempLES(:,iCell) = hom(:,1,4,0)
    ! salinityLES(:,iCell) = hom(:,1,5,0)

! TODO: for testing in single column mode <20190724, Qing Li> !
    if (iCell == 1) then

#if defined( __cudaProfiler )
!-- Only profile time_integration
    CALL cudaProfilerStart()
#endif
!
!-- Integration of the model equations using timestep-scheme
    CALL time_integration

#if defined( __cudaProfiler )
!-- Only profile time_integration
    CALL cudaProfilerStop()
#endif
!
!-- Take final CPU-time for CPU-time analysis
    CALL cpu_log( log_point(1), 'total', 'stop' )
!    CALL cpu_statistics
! flow_statistics_called = .FALSE.

! call flow_statistics
    if(average_count_meanpr /= 0) then

       meanFields_avg(:,1) = meanFields_avg(:,1) / average_count_meanpr
       meanFields_avg(:,2) = meanFields_avg(:,2) / average_count_meanpr
       meanFields_avg(:,3) = meanFields_avg(:,3) / average_count_meanpr
       meanFields_avg(:,4) = meanFields_avg(:,4) / average_count_meanpr

    endif

    endif ! iCell

    Tles = meanFields_avg(:,3)
    Sles = meanFields_avg(:,4)
    Ules = meanFields_avg(:,1)
    Vles = meanFields_avg(:,2)
 ! need to integrate over layers in mpas to get increments

 if(minval(tempLES(:,iCell)) < 100.0_wp) tempLES(:,iCell) = tempLES(:,iCell) + 273.15_wp
    ! tProfileInit(1:) = tempLES(:,iCell)
    ! sProfileInit(1:) = salinityLES(:,iCell)
    ! uProfileInit(1:) = uLESout(:,iCell)
    ! vProfileInit(1:) = vLESout(:,iCell)
    tProfileInit(1:) = t_mean_restart(1:nzt,iCell)
    sProfileInit(1:) = s_mean_restart(1:nzt,iCell)
    uProfileInit(1:) = u_mean_restart(1:nzt,iCell)
    vProfileInit(1:) = v_mean_restart(1:nzt,iCell)
    tempLES(:,iCell) = Tles(1:nzt)
    salinityLES(:,iCell) = Sles(1:nzt)
    uLESout(:,iCell) = Ules(1:nzt)
    vLESout(:,iCell) = Vles(1:nzt)

    if (iCell == 1) then
    jl=1
    do il=nzt,nzb+1,-1
      fLES(1,1,jl) = Tles(il) - 273.15_wp
      fLES(1,2,jl) = Sles(il)
      fLES(1,3,jl) = Ules(il)
      fLES(1,4,jl) = Vles(il)
      jl = jl+1
    enddo
    fLES(:,:,nz+1) = fLES(:,:,nz)
    format = "(6x, F10.3, F10.3, F10.3, F10.3, F10.3)"
    print *, 'fLES1         Z         T         S         U         V'
    sumValT = 0.0_wp
    sumValS = 0.0_wp
    sumValU = 0.0_wp
    sumValV = 0.0_wp
    do il = 1,nzLES
      write(*,format) 0.5*(zeLESinv(il)+zeLESinv(il+1)), fLES(1,1,il), fLES(1,2,il), fLES(1,3,il), fLES(1,4,il)
      sumValT = sumValT + fLES(1,1,il) * (zeLESinv(il)-zeLESinv(il+1))
      sumValS = sumValS + fLES(1,2,il) * (zeLESinv(il)-zeLESinv(il+1))
      sumValU = sumValU + fLES(1,3,il) * (zeLESinv(il)-zeLESinv(il+1))
      sumValV = sumValV + fLES(1,4,il) * (zeLESinv(il)-zeLESinv(il+1))
    enddo
    print *, '  Sum         Z         T         S         U         V'
    write(*,format) 0.0_wp, sumValT, sumValS, sumValU, sumValV
    call rmap1d(nz+2,nzMPAS+1,nvar,ndof,abs(zeLESinv(1:nzLES+2)),abs(zedge(1:nzMPAS+1)), &
                fLES,fMPAS(:,:,:nzMPAS),bc_l,bc_r,work,opts)
    sumValT = 0.0_wp
    sumValS = 0.0_wp
    sumValU = 0.0_wp
    sumValV = 0.0_wp
    print *, 'fMPAS         Z         T         S         U         V'
    do il = 1,nzMPAS
      write(*,format) 0.5*(zedge(il)+zedge(il+1)), fMPAS(1,1,il), fMPAS(1,2,il), fMPAS(1,3,il), fMPAS(1,4,il)
      sumValT = sumValT + fMPAS(1,1,il) * (zedge(il)-zedge(il+1))
      sumValS = sumValS + fMPAS(1,2,il) * (zedge(il)-zedge(il+1))
      sumValU = sumValU + fMPAS(1,3,il) * (zedge(il)-zedge(il+1))
      sumValV = sumValV + fMPAS(1,4,il) * (zedge(il)-zedge(il+1))
    enddo
    print *, '  Sum         Z         T         S         U         V'
    write(*,format) 0.0_wp, sumValT, sumValS, sumValU, sumValV

    call rmap1d(nzMPAS+1,nz+2,nvar,ndof,abs(zedge(1:nzMPAS+1)),abs(zeLESinv(1:nzLES+2)), &
                fMPAS(:,:,:nzMPAS), fLES, bc_l, bc_r, work, opts)
    print *, 'fLES2         Z         T         S         U         V'
    sumValT = 0.0_wp
    sumValS = 0.0_wp
    sumValU = 0.0_wp
    sumValV = 0.0_wp
    do il = 1,nzLES
      write(*,format) 0.5*(zeLESinv(il)+zeLESinv(il+1)), fLES(1,1,il), fLES(1,2,il), fLES(1,3,il), fLES(1,4,il)
      sumValT = sumValT + fLES(1,1,il) * (zeLESinv(il)-zeLESinv(il+1))
      sumValS = sumValS + fLES(1,2,il) * (zeLESinv(il)-zeLESinv(il+1))
      sumValU = sumValU + fLES(1,3,il) * (zeLESinv(il)-zeLESinv(il+1))
      sumValV = sumValV + fLES(1,4,il) * (zeLESinv(il)-zeLESinv(il+1))
    enddo
    print *, '  Sum         Z         T         S         U         V'
    write(*,format) 0.0_wp, sumValT, sumValS, sumValU, sumValV
    endif

    jl=1
    do il=nzt,nzb+1,-1
      fLES(1,1,jl) = (Tles(il) - tProfileInit(il)) / dtLS
      fLES(1,2,jl) = (Sles(il) - sProfileInit(il)) / dtLS
      fLES(1,3,jl) = (Ules(il) - uProfileInit(il)) / dtLS
      fLES(1,4,jl) = (Vles(il) - vProfileInit(il)) / dtLS
      jl = jl+1
    enddo
    fLES(:,:,nz+1) = fLES(:,:,nz)

    call rmap1d(nz+2,nzMPAS+1,nvar,ndof,abs(zeLESinv(1:nzLES+2)),abs(zedge(1:nzMPAS+1)), &
                fLES,fMPAS(:,:,:nzMPAS),bc_l,bc_r,work,opts)

    idx_cutoff = nzMPAS
    do il=1,nVertLevels
     if(zmid(il) < zw(0)) then
       idx_cutoff = il
       exit
     endif
    enddo
    fMPAS(:,:,idx_cutoff:nzMPAS) = 0.0_wp

    if (iCell == 1) then
    format = "(6x, F14.3, E14.3, E14.3, E14.3, E14.3)"
    print *, ' fLES             Z          dTdt          dSdt          dUdt          dVdt'
    do il = 1,nzLES
      write(*,format) 0.5*(zeLESinv(il)+zeLESinv(il+1)), fLES(1,1,il), fLES(1,2,il), fLES(1,3,il), fLES(1,4,il)
    enddo
    print *, 'fMPAS             Z          dTdt          dSdt          dUdt          dVdt'
    do il = 1,nzMPAS
      write(*,format) 0.5*(zedge(il)+zedge(il+1)), fMPAS(1,1,il), fMPAS(1,2,il), fMPAS(1,3,il), fMPAS(1,4,il)
    enddo
    endif

    tIncrementLES(:,iCell) = 0.0_wp
    sIncrementLES(:,iCell) = 0.0_wp
    uIncrementLES(:,iCell) = 0.0_wp
    vIncrementLES(:,iCell) = 0.0_wp
    do jl=1,nzMPAS
      tIncrementLES(jl,iCell) = fMPAS(1,1,jl)
      sIncrementLES(jl,iCell) = fMPAS(1,2,jl)
      uIncrementLES(jl,iCell) = fMPAS(1,3,jl)
      vIncrementLES(jl,iCell) = fMPAS(1,4,jl)
    enddo

   !put variables in arrays to continue
    u_restart(:,:,:,iCell) = u(:,:,:)
    v_restart(:,:,:,iCell) = v(:,:,:)
    w_restart(:,:,:,iCell) = w(:,:,:)
    pt_restart(:,:,:,iCell) = pt(:,:,:)
    sa_restart(:,:,:,iCell) = sa(:,:,:)
    e_restart(:,:,:,iCell) = e(:,:,:)
    km_restart(:,:,:,iCell) = km(:,:,:)
    kh_restart(:,:,:,iCell) = kh(:,:,:)
    u_mean_restart(:,iCell) = Ules
    v_mean_restart(:,iCell) = Vles
    t_mean_restart(:,iCell) = Tles
    s_mean_restart(:,iCell) = Sles
    call init_control_parameters
  enddo !ends icell loop

  ! do iCell=2,nCells
  !   u_restart(:,:,:,iCell) = u(:,:,:)
  !   v_restart(:,:,:,iCell) = v(:,:,:)
  !   w_restart(:,:,:,iCell) = w(:,:,:)
  !   pt_restart(:,:,:,iCell) = pt(:,:,:)
  !   sa_restart(:,:,:,iCell) = sa(:,:,:)
  !   e_restart(:,:,:,iCell) = e(:,:,:)
  !   km_restart(:,:,:,iCell) = km(:,:,:)
  !   kh_restart(:,:,:,iCell) = kh(:,:,:)
  !   u_mean_restart(:,iCell) = u_mean_restart(:,1)
  !   v_mean_restart(:,iCell) = v_mean_restart(:,1)
  !   t_mean_restart(:,iCell) = t_mean_restart(:,1)
  !   s_mean_restart(:,iCell) = s_mean_restart(:,1)
  !   tIncrementLES(:,iCell) = tIncrementLES(:,1)
  !   sIncrementLES(:,iCell) = sIncrementLES(:,1)
  !   uIncrementLES(:,iCell) = uIncrementLES(:,1)
  !   vIncrementLES(:,iCell) = vIncrementLES(:,1)
  !   tempLES(:,iCell) = tempLES(:,1)
  !   salinityLES(:,iCell) = salinityLES(:,1)
  !   uLESout(:,iCell) = uLESout(:,1)
  !   vLESout(:,iCell) = vLESout(:,1)
  !   zLES(:,iCell) = zLES(:,1)
  ! enddo

END subroutine palm_main

subroutine palm_finalize()
   !this dallocates lots of stuff
 DEALLOCATE( pt_init, q_init, s_init, ref_state, sa_init, ug,         &
                       u_init, v_init, vg, hom, hom_sum, meanFields_avg )

   deallocate(hor_index_bounds)

    deallocate(zu,zeLES,Tles,Sles)
    deallocate(hyp, Ules,Vles)
    deallocate(ddzu, ddzw, dd2zu, dzu, dzw, zw, ddzu_pres, nzb_s_inner,  &
               nzb_s_outer, nzb_u_inner, nzb_u_outer, nzb_v_inner,       &
               nzb_v_outer, nzb_w_inner, nzb_w_outer, nzb_diff_s_inner,  &
               nzb_diff_s_outer, wall_flags_0, advc_flags_1, advc_flags_2)

    deallocate(u_stk, v_stk, u_stk_zw, v_stk_zw )
   deallocate(zmid,zedge)
   deallocate(T_mpas2,S_mpas2,U_mpas2)
   deallocate(V_mpas2)

    call ws_finalize
    call deallocate_bc
    call deallocate_3d_variables
    call tcm_deallocate_arrays
    if (random_generator == 'random-parallel') call deallocate_random_generator
    call tridia_deallocate
    close(18)

    CALL fft_finalize
#if defined( __parallel )
    CALL MPI_FINALIZE( ierr )
#endif


end subroutine palm_finalize

subroutine init_control_parameters
    USE arrays_3d

    USE control_parameters

    USE statistics, only: flow_statistics_called
    USE kinds

    openfile = file_status(.FALSE.,.FALSE.)

    rayleigh_damping_factor = -1.0_wp
    rayleigh_damping_height = -1.0_wp
    timestep_count = 0
!        poisfft_initialized = .FALSE.
!        init_fft = .FALSE.
        psolver = 'poisfft'
        momentum_advec = 'ws-scheme'
        loop_optimization = 'vector'
        bc_e_b = 'neumann'
        bc_lr = 'cyclic'
        bc_ns = 'cyclic'
        bc_p_b = 'neumann'
        bc_p_t = 'neumann'
        bc_pt_b = 'neumann'
        bc_pt_t = 'neumann'
        bc_sa_t = 'neumann'
        bc_sa_b = 'neumann'
        bc_uv_b = 'neumann'
        bc_uv_t = 'neumann'
        ibc_uv_b = 1
        coupling_mode = 'uncoupled'
        fft_method = 'temperton-algorithm'
        topography = 'flat'
        initializing_actions = 'set_constant_profiles'
        random_generator = 'numerical-recipes'
        random_generator = 'random-parallel'
        reference_state = 'initial_profile'
        data_output = ' '
        data_output_user = ' '
        doav = ' '
        data_output_masks = ' '
        data_output_pr = ' '
        domask = ' '
        do2d = ' '
        do3d = ' '

        do3d_no(0:1) = 0
        ! meanFields_avg = 0.0_wp

        abort_mode = 1
        time_avg = 0.0_wp
        average_count_pr = 0
        average_count_meanpr = 0
        average_count_3d = 0
        current_timestep_number = 0
        coupling_topology = 0
        dist_range = 0
        doav_n = 0
        dopr_n = 0
        dopr_time_count = 0
        dopts_time_count = 0
        dots_time_count = 0
        dp_level_ind_b = 0
        dvrp_filecount = 0
        ensemble_member_nr = 0

        iran = -1234567
        length = 0
        io_group = 0
        io_blocks = 1
        masks = 0
        maximum_parallel_io_streams = -1
        mgcycles = 0
        mg_cycles = 4
        mg_switch_to_pe0_level = -1
        ngsrb = 2
        nr_timesteps_this_run = 0
        nsor = 20
        nsor_ini = 100
        normalizing_region = 0
        num_leg = 0
        num_var_fl_user = 0
        nz_do3d = -9999
        y_shift = 0
        mask_size(max_masks,3) = -1
        mask_size_l(max_masks,3) = -1
        mask_start_l(max_masks,3) = -1
        pt_vertical_gradient_level_ind(10) = -9999
        sa_vertical_gradient_level_ind(10) = -9999
 !       stokes_drift_method = -9999

 !       dz(10) = -1.0_wp
 !       dzconst = 2.5_wp
!        dt_disturb = 20.0_wp
!        dt_do3d = 9999999.9_wp
!        dt_3d = 0.01_wp

        simulated_time = 0.0_wp
        flow_statistics_called = .FALSE.
!        disturbance_created = .FALSE.
        time_disturb = 0.0_wp
        time_dopr = 0.0_wp
        time_dopr_av = 0.0_wp
        time_dots = 0.0_wp
        time_do2d_xy = 0.0_wp
        time_do2d_xz = 0.0_wp
        time_do2d_yz = 0.0_wp
        time_do3d = 0.0_wp
        time_do_av = 0.0_wp
        time_run_control = 0.0_wp

end subroutine init_control_parameters

subroutine deallocate_memory

        use pegrid

        deallocate(hor_index_bounds)

end subroutine deallocate_memory

end module palm_mod
