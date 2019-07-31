!> @file diffusion_s.f90
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
!> Diffusion term of scalar quantities (temperature and water content)
!------------------------------------------------------------------------------!
 MODULE diffusion_s_mod
 

    PRIVATE
    PUBLIC diffusion_s

    INTERFACE diffusion_s
       MODULE PROCEDURE diffusion_s
    END INTERFACE diffusion_s

 CONTAINS


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s( s, s_flux_t, s_flux_solar_t)

       USE arrays_3d,                                                          &
           ONLY:  dzw, ddzu, ddzw, kh, tend, drho_air, rho_air_zw, solar3d
       
       USE control_parameters,                                                 & 
           ONLY: use_top_fluxes, ideal_solar_division,     &
                 ideal_solar_efolding1, ideal_solar_efolding2
 
       USE grid_variables,                                                     &
           ONLY:  ddx, ddx2, ddy, ddy2
       
       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb,             &
                  nzt, wall_flags_0
       
       USE kinds

       USE surface_mod,                                                        &
           ONLY :  surf_def_h

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  zval              !< depth_variable for solar penetration
       REAL(wp) ::  flux1             !< solar flux temp variable
       REAL(wp) ::  flux2             !< solar flux temp variable
       REAL(wp) ::  flag
       REAL(wp) ::  mask_bottom       !< flag to mask vertical upward-facing surface     
       REAL(wp) ::  mask_top          !< flag to mask vertical downward-facing surface  

       REAL(wp) ::  s_flux_t           !< flux at model top

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !< 
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  s  !< 
#endif

       REAL(wp),INTENT(IN),OPTIONAL :: s_flux_solar_t  !<solar flux at sfc

!-- Compute horizontal diffusion

       !$acc data present( tend ) &
       !$acc present( s ) &
       !$acc present( solar3d ) &
       !$acc present( kh, surf_def_h ) &
       !$acc present( ddzu, ddzw, dzw, rho_air_zw, drho_air, wall_flags_0 )

       !$acc parallel
       !$acc loop
       DO  i = nxl, nxr
          !$acc loop
          DO  j = nys,nyn
             !$acc loop
             DO  k = nzb+1, nzt
                tend(k,j,i) = tend(k,j,i)                                      &
                                          + 0.5_wp * (                         &
                                     ( kh(k,j,i) + kh(k,j,i+1) )               &
                                   * ( s(k,j,i+1) - s(k,j,i)   )               &
                                   - ( kh(k,j,i) + kh(k,j,i-1) )               &
                                   * ( s(k,j,i)   - s(k,j,i-1) )               &
                                                     ) * ddx2                  &
                                          + 0.5_wp * (                         &
                                     ( kh(k,j,i) + kh(k,j+1,i) )               &
                                   * ( s(k,j+1,i) - s(k,j,i)   )               &
                                   - ( kh(k,j,i) + kh(k,j-1,i) )               &
                                   * ( s(k,j,i)   - s(k,j-1,i) )               &
                                                     ) * ddy2
             ENDDO
          ENDDO
             ENDDO

       !$acc loop
       DO i = nxl, nxr
          !$acc loop
          DO j = nys, nyn
!
!--          Compute vertical diffusion. In case that surface fluxes have been
!--          prescribed or computed at bottom and/or top, index k starts/ends at
!--          nzb+2 or nzt-1, respectively. Model top is also mask if top flux
!--          is given.
             !$acc loop
             DO  k = nzb+1, nzt
!
!--             Determine flags to mask topography below and above. Flag 0 is 
!--             used to mask topography in general, and flag 8 implies 
!--             information about use_surface_fluxes. Flag 9 is used to control 
!--             flux at model top. 
                mask_bottom = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k-1,j,i), 8 ) ) 
                mask_top    = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k+1,j,i), 8 ) ) *     &
                              MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k+1,j,i), 9 ) ) 
                flag        = MERGE( 1.0_wp, 0.0_wp,                           &
                                     BTEST( wall_flags_0(k,j,i), 0 ) )

                tend(k,j,i) = tend(k,j,i)                                      &
                                       + 0.5_wp * (                            &
                                      ( kh(k,j,i) + kh(k+1,j,i) ) *            &
                                        ( s(k+1,j,i)-s(k,j,i) )  * ddzu(k+1)  &
                                                            * rho_air_zw(k)    &
                                                            * mask_top         &
                                    - ( kh(k,j,i) + kh(k-1,j,i) ) *            &
                                      ( s(k,j,i)-s(k-1,j,i) )  * ddzu(k)    &
                                                           * rho_air_zw(k-1)  &
                                                            * mask_bottom      &
                                                  ) * ddzw(k) * drho_air(k)  !  &
                                                            !  * flag
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

       IF ( PRESENT(s_flux_solar_t )) THEN

       !$acc parallel
       !$acc loop collapse(2)
       DO i=nxl,nxr
          DO j=nys,nyn
                !LPV adding solar forcing with depth
!                IF ( PRESENT(s_flux_solar_t )) THEN
              !    m = surf_def_h(2)%start_index(j,i)
                  zval = 0.0_wp
                !$acc loop
                  DO k = nzt,nzb+1,-1
                      flux1 = (1.0_wp - ideal_solar_division)*exp(ideal_solar_efolding2*zval) + &
                                ideal_solar_division*exp(ideal_solar_efolding1*zval)
                     zval = zval - dzw(k)
                      flux2 = (1.0_wp - ideal_solar_division)*exp(ideal_solar_efolding2*zval) + &
                      ideal_solar_division*exp(ideal_solar_efolding1*zval)

                      tend(k,j,i) = tend(k,j,i) - s_flux_solar_t*(flux1 - flux2) / dzw(k)
                    
                      solar3d(k,j,i) = -s_flux_solar_t*(flux1 - flux2) / dzw(k)
                  ENDDO
          enddo
      enddo
      !$acc end parallel
      ENDIF


      if ( use_top_fluxes ) then
       !$acc parallel
       !$acc loop collapse(2)
       DO i=nxl,nxr
          DO j=nys,nyn

!--          Vertical diffusion at the last computational gridpoint along z-direction
                  k = nzt
                   tend(k,j,i) = tend(k,j,i)                                   &
                           + ( - s_flux_t ) * ddzw(k) * drho_air(k)
          ENDDO
       ENDDO
       !$acc end parallel
     endif
       !$acc end data

    END SUBROUTINE diffusion_s

 END MODULE diffusion_s_mod
