!> @file diffusion_u.f90
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
!> Diffusion term of the u-component
!> @todo additional damping (needed for non-cyclic bc) causes bad vectorization
!>       and slows down the speed on NEC about 5-10%
!------------------------------------------------------------------------------!
 MODULE diffusion_u_mod
 

    PRIVATE
    PUBLIC diffusion_u

    INTERFACE diffusion_u
       MODULE PROCEDURE diffusion_u
    END INTERFACE diffusion_u

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_u

       USE arrays_3d,                                                          &
           ONLY:  ddzu, ddzw, km, tend, u, v, w, drho_air, rho_air_zw
       
       USE control_parameters,                                                 &
           ONLY:  constant_top_momentumflux, use_top_fluxes, top_momentumflux_u
       
       USE grid_variables,                                                     &
           ONLY:  ddx, ddx2, ddy
       
       USE indices,                                                            &
           ONLY:  nxl, nxlu, nxr, nyn, nys, nzb, nzt, wall_flags_0
     
       USE kinds

!       USE surface_mod,                                                        &
!           ONLY :  surf_def_h

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  flag          !< flag to mask topography grid points
       REAL(wp)     ::  kmym          !< 
       REAL(wp)     ::  kmyp          !<
       REAL(wp)     ::  kmzm          !<
       REAL(wp)     ::  kmzp          !<
       REAL(wp)     ::  mask_bottom   !< flag to mask vertical upward-facing surface       
       REAL(wp)     ::  mask_top      !< flag to mask vertical downward-facing surface 


!-- Compute horizontal diffusion
       !$acc data present( tend ) &
       !$acc present( u, v, w ) &
       !$acc present( km ) &
       !$acc present( ddzu, ddzw, rho_air_zw, drho_air, wall_flags_0 )

       !$acc parallel
       !$acc loop
       DO  i = nxlu, nxr
          !$acc loop
          DO  j = nys, nyn
             !$acc loop
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography and wall-bounded grid points. 
!--             It is sufficient to masked only north- and south-facing surfaces, which
!--             need special treatment for the u-component. 
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmyp = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k,j+1,i)+km(k,j,i-1)+km(k,j+1,i-1) )
                kmym = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k,j-1,i)+km(k,j,i-1)+km(k,j-1,i-1) )

               tend(k,j,i) = tend(k,j,i)                                      &
                        + 2.0_wp * (                                           &
                                  km(k,j,i)   * ( u(k,j,i+1) - u(k,j,i)   )    &
                                - km(k,j,i-1) * ( u(k,j,i)   - u(k,j,i-1) )    &
                                   ) * ddx2                                    &
                        +          ( (                                         &
                            kmyp * ( u(k,j+1,i) - u(k,j,i)     ) * ddy         &
                          + kmyp * ( v(k,j+1,i) - v(k,j+1,i-1) ) * ddx         &
                                                  )                            &
                                   - (                                         &
                            kmym * ( u(k,j,i) - u(k,j-1,i) ) * ddy             &
                          + kmym * ( v(k,j,i) - v(k,j,i-1) ) * ddx             &
                                                  )                            &
                                   ) * ddy
             ENDDO
                ENDDO   
             ENDDO
!
!--          Compute vertical diffusion. In case of simulating a surface layer,
!--          respective grid diffusive fluxes are masked (flag 8) within this 
!--          loop, and added further below, else, simple gradient approach is
!--          applied. Model top is also mask if top-momentum flux is given.

       !$acc loop
       DO  i = nxlu, nxr
          !$acc loop
          DO  j = nys, nyn
             !$acc loop
             DO  k = nzb+1, nzt
!
!--             Determine flags to mask topography below and above. Flag 1 is 
!--             used to mask topography in general, and flag 8 implies 
!--             information about use_surface_fluxes. Flag 9 is used to control 
!--             momentum flux at model top.  
                mask_bottom = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k-1,j,i), 8 ) ) 
                mask_top    = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k+1,j,i), 8 ) ) *     &
                              MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k+1,j,i), 9 ) ) 
                flag        = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k,j,i), 1 ) ) 
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) )
                kmzm = 0.25_wp *                                               &
                       ( km(k,j,i)+km(k-1,j,i)+km(k,j,i-1)+km(k-1,j,i-1) )

                tend(k,j,i) = tend(k,j,i)                                      &
                        + ( kmzp * ( ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)   &
                                   + ( w(k,j,i)   - w(k,j,i-1) ) * ddx         &
                                   ) * rho_air_zw(k)   * mask_top              &
                          - kmzm * ( ( u(k,j,i)   - u(k-1,j,i)   ) * ddzu(k)   &
                                   + ( w(k-1,j,i) - w(k-1,j,i-1) ) * ddx       &
                                   ) * rho_air_zw(k-1) * mask_bottom           &
                          ) * ddzw(k) * drho_air(k) * flag
             ENDDO
                ENDDO
                ENDDO
       !$acc end parallel

       !$acc parallel
       !$acc loop collapse(2)
       DO  i = nxlu, nxr
          DO  j = nys, nyn
!
!--          Add momentum flux at model top
             IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
!                surf_s = surf_def_h(2)%start_index(j,i)
!                surf_e = surf_def_h(2)%end_index(j,i)
!                !$acc loop
!                DO  m = surf_s, surf_e

                   k   = nzt !surf_def_h(2)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                   &
                        + ( - top_momentumflux_u ) * ddzw(k) * drho_air(k)
!                ENDDO
             ENDIF

          ENDDO
       ENDDO
       !$acc end parallel
       !$acc end data

    END SUBROUTINE diffusion_u


 END MODULE diffusion_u_mod
