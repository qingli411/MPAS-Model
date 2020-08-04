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
! Current revisions:
! ------------------
!
!
! Former revisions:
! -----------------
! $Id: diffusion_s.f90 2759 2018-01-17 16:24:59Z suehring $
! Major bugfix, horizontal diffusion at vertical surfaces corrected.
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Adjustments to new topography and surface concept
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC version of subroutine removed
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1873 2016-04-18 14:50:06Z maronga
! Module renamed (removed _mod)
!
! 1850 2016-04-08 13:29:27Z maronga
! Module renamed
!
! 1691 2015-10-26 16:17:44Z maronga
! Formatting corrections.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1374 2014-04-25 12:55:07Z raasch
! missing variables added to ONLY list
!
! 1340 2014-03-25 19:45:13Z kanani
! REAL constants defined as wp-kind
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
! openacc loop and loop vector clauses removed
!
! 1128 2013-04-12 06:19:32Z raasch
! loop index bounds in accelerator version replaced by i_left, i_right, j_south,
! j_north
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1015 2012-09-27 09:23:24Z raasch
! accelerator version (*_acc) added
!
! 1010 2012-09-20 07:59:54Z raasch
! cpp switch __nopointer added for pointer free version
!
! 1001 2012-09-13 14:08:46Z raasch
! some arrays comunicated by module instead of parameter list
!
! Revision 1.1  2000/04/13 14:54:02  schroeter
! Initial revision
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
    SUBROUTINE diffusion_s( s, s_flux_t)

       USE arrays_3d,                                                          &
           ONLY:  dzw, ddzu, ddzw, kh, tend

       USE grid_variables,                                                     &
           ONLY:  ddx, ddx2, ddy, ddy2

       USE indices,                                                            &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt

       USE kinds

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction

       REAL(wp) ::  s_flux_t           !< flux at model top

#if defined( __nopointer )
       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !<
#else
       REAL(wp), DIMENSION(:,:,:), POINTER ::  s  !<
#endif

!-- Compute horizontal diffusion

       !$acc data present( tend ) &
       !$acc present( s, kh ) &
       !$acc present( ddzu, ddzw, dzw )

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
                                                     ) * ddy2                  &
                                       + 0.5_wp * (                            &
                                      ( kh(k,j,i) + kh(k+1,j,i) ) *            &
                                          ( s(k+1,j,i)-s(k,j,i) ) * ddzu(k+1)  &
                                    - ( kh(k,j,i) + kh(k-1,j,i) ) *            &
                                          ( s(k,j,i)-s(k-1,j,i) ) * ddzu(k)    &
                                                  ) * ddzw(k)
             ENDDO
          ENDDO
       ENDDO
       !$acc end parallel

!
!--    Add top fluxes
       !$acc parallel
       !$acc loop collapse(2)
       DO i=nxl,nxr
          DO j=nys,nyn
             tend(nzt,j,i) = tend(nzt,j,i)                                  &
                           + ( - s_flux_t ) * ddzw(nzt)
          ENDDO
       ENDDO
       !$acc end parallel
       !$acc end data

    END SUBROUTINE diffusion_s

 END MODULE diffusion_s_mod
