!> @file boundary_conds.f90
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
! 2018-10-25 cbegeman
! Treatment of dirichlet bottom boundary conditions for salinity
!
! Former revisions:
! -----------------
! $Id: boundary_conds.f90 2938 2018-03-27 15:52:42Z suehring $
! Set boundary condition for TKE and TKE dissipation rate in case of nesting
! and if parent model operates in RANS mode but child model in LES mode.
! mode
!
! 2793 2018-02-07 10:54:33Z suehring
! Removed preprocessor directive __chem
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Adjust boundary conditions for e and diss in case of TKE-e closure (TG)
! Implementation of chemistry module (FK)
!
! 2569 2017-10-20 11:54:42Z kanani
! Removed redundant code for ibc_s_b=1 and ibc_q_b=1
!
! 2365 2017-08-21 14:59:59Z kanani
! Vertical grid nesting implemented: exclude setting vertical velocity to zero
! on fine grid (SadiqHuq)
!
! 2320 2017-07-21 12:47:43Z suehring
! Remove unused control parameter large_scale_forcing from only-list
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Set boundary conditions on topography top using flag method.
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC directives removed
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1992 2016-08-12 15:14:59Z suehring
! Adjustments for top boundary condition for passive scalar
!
! 1960 2016-07-12 16:34:24Z suehring
! Treat humidity and passive scalar separately
!
! 1823 2016-04-07 08:57:52Z hoffmann
! Initial version of purely vertical nesting introduced.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme removed. microphyisics_seifert added.
!
! 1764 2016-02-28 12:45:19Z raasch
! index bug for u_p at left outflow removed
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1742 2016-01-13 09:50:06Z raasch
! bugfix for outflow Neumann boundary conditions at bottom and top
!
! 1717 2015-11-11 15:09:47Z raasch
! Bugfix: index error in outflow conditions for left boundary
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1410 2014-05-23 12:16:18Z suehring
! Bugfix: set dirichlet boundary condition for passive_scalar at model domain
! top
!
! 1399 2014-05-07 11:16:25Z heinze
! Bugfix: set inflow boundary conditions also if no humidity or passive_scalar
! is used.
!
! 1398 2014-05-07 11:15:00Z heinze
! Dirichlet-condition at the top for u and v changed to u_init and v_init also
! for large_scale_forcing
!
! 1380 2014-04-28 12:40:45Z heinze
! Adjust Dirichlet-condition at the top for pt in case of nudging
!
! 1361 2014-04-16 15:17:48Z hoffmann
! Bottom and top boundary conditions of rain water content (qr) and
! rain drop concentration (nr) changed to Dirichlet
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1257 2013-11-08 15:18:40Z raasch
! loop independent clauses added
!
! 1241 2013-10-30 11:36:58Z heinze
! Adjust ug and vg at each timestep in case of large_scale_forcing
!
! 1159 2013-05-21 11:58:22Z fricke
! Bugfix: Neumann boundary conditions for the velocity components at the
! outflow are in fact radiation boundary conditions using the maximum phase
! velocity that ensures numerical stability (CFL-condition).
! Hence, logical operator use_cmax is now used instead of bc_lr_dirneu/_neudir.
! Bugfix: In case of use_cmax at the outflow, u, v, w are replaced by
! u_p, v_p, w_p
!
! 1115 2013-03-26 18:16:16Z hoffmann
! boundary conditions of two-moment cloud scheme are restricted to Neumann-
! boundary-conditions
!
! 1113 2013-03-10 02:48:14Z raasch
! GPU-porting
! dummy argument "range" removed
! Bugfix: wrong index in loops of radiation boundary condition
!
! 1053 2012-11-13 17:11:03Z hoffmann
! boundary conditions for the two new prognostic equations (nr, qr) of the
! two-moment cloud scheme
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting
!
! 978 2012-08-09 08:28:32Z fricke
! Neumann boudnary conditions are added at the inflow boundary for the SGS-TKE.
! Outflow boundary conditions for the velocity components can be set to Neumann
! conditions or to radiation conditions with a horizontal averaged phase
! velocity.
!
! 875 2012-04-02 15:35:15Z gryschka
! Bugfix in case of dirichlet inflow bc at the right or north boundary
!
! Revision 1.1  1997/09/12 06:21:34  raasch
! Initial revision
!
!
! Description:
! ------------
!> Boundary conditions for the prognostic quantities.
!> One additional bottom boundary condition is applied for the TKE (=(u*)**2)
!> in prandtl_fluxes. The cyclic lateral boundary conditions are implicitly
!> handled in routine exchange_horiz. Pressure boundary conditions are
!> explicitly set in routines pres, poisfft, poismg and sor.
!------------------------------------------------------------------------------!
 SUBROUTINE boundary_conds


    USE arrays_3d,                                                             &
        ONLY:  e_p, pt, pt_p, sa, sa_p,                                        &
               u_init, u_p, v_init, v_p, w_p

    USE control_parameters,                                                    &
        ONLY:  dt_3d, ibc_pt_b, ibc_pt_t, ibc_sa_b, ibc_sa_t,                  &
               ibc_uv_b, ibc_uv_t, intermediate_timestep_count

    USE grid_variables,                                                        &
        ONLY:  ddx, ddy, dx, dy

    USE indices,                                                               &
        ONLY:  nx, nxl, nxlg, nxr, nxrg, ny, nyn, nyng, nys, nysg, nzb, nzt

    USE kinds

    USE pegrid

    IMPLICIT NONE

    INTEGER(iwp) ::  i  !< grid index x direction
    INTEGER(iwp) ::  j  !< grid index y direction

!
!-- Bottom boundary
    !$acc parallel present( u_p, v_p )
    IF ( ibc_uv_b == 1 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             u_p(nzb,j,i) = u_p(nzb+1,j,i)
             v_p(nzb,j,i) = v_p(nzb+1,j,i)
          ENDDO
       ENDDO
    ELSE
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
             u_p(nzb,j,i) = 0.0_wp
             v_p(nzb,j,i) = 0.0_wp
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel
!
!-- Top boundary.
    !$acc parallel present( u_p, v_p, u_init, v_init )
    IF ( ibc_uv_t == 0 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              u_p(nzt+1,j,i) = u_init(nzt+1)
              v_p(nzt+1,j,i) = v_init(nzt+1)
          ENDDO
       ENDDO
    ELSEIF ( ibc_uv_t == 1 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              u_p(nzt+1,j,i) = u_p(nzt,j,i)
              v_p(nzt+1,j,i) = v_p(nzt,j,i)
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel

!
!-- Set zero vertical velocity at top and bottom
    !$OMP PARALLEL DO PRIVATE( i, j )
    !$acc parallel present( w_p )
    !$acc loop collapse(2)
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          w_p(nzb,j,i) = 0.0_wp
          w_p(nzt:nzt+1,j,i) = 0.0_wp  !< nzt is not a prognostic level (but cf. pres)
       ENDDO
    ENDDO
    !$acc end parallel

!
!-- Temperature at bottom and top boundary.
!-- Dirichlet
    !$acc parallel present( pt_p, pt )
    IF ( ibc_pt_b == 0 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              pt_p(nzb,j,i) = pt(nzb,j,i)
          ENDDO
       ENDDO
!-- Neumann, zero-gradient
    ELSEIF ( ibc_pt_b == 1 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              pt_p(nzb,j,i) = pt_p(nzb+1,j,i)
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel

!
!-- Temperature at top boundary
    !$acc parallel present( pt_p, pt )
    IF ( ibc_pt_t == 0 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              pt_p(nzt+1,j,i) = pt(nzt+1,j,i)
          ENDDO
       ENDDO
    ELSEIF ( ibc_pt_t == 1 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              pt_p(nzt+1,j,i) = pt_p(nzt,j,i)
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel

!
!-- Boundary conditions for TKE.
!-- Generally Neumann conditions with de/dz=0 are assumed.

    !$acc parallel present( e_p )
    !$acc loop collapse(2)
    DO  i = nxlg, nxrg
       DO  j = nysg, nyng
          e_p(nzb,j,i) = e_p(nzt+1,j,i)
          e_p(nzt+1,j,i) = e_p(nzt,j,i)
       ENDDO
    ENDDO
    !$acc end parallel

!
!-- Salinity at bottom and top boundary.
!-- Dirichlet
    !$acc parallel present( sa_p, sa )
    IF ( ibc_sa_b == 0 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              sa_p(nzb,j,i) = sa(nzb,j,i)
          ENDDO
       ENDDO
!-- Neumann, zero-gradient
    ELSEIF ( ibc_sa_b == 1 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              sa_p(nzb,j,i) = sa_p(nzb+1,j,i)
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel

!
!-- Salinity at top boundary
    !$acc parallel present( sa_p, sa )
    IF ( ibc_sa_t == 0 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              sa_p(nzt+1,j,i) = sa(nzt+1,j,i)
          ENDDO
       ENDDO
    ELSEIF ( ibc_sa_t == 1 )  THEN
       !$OMP PARALLEL DO PRIVATE( i, j )
       !$acc loop collapse(2)
       DO  i = nxlg, nxrg
          DO  j = nysg, nyng
              sa_p(nzt+1,j,i) = sa_p(nzt,j,i)
          ENDDO
       ENDDO
    ENDIF
    !$acc end parallel


 END SUBROUTINE boundary_conds
