!> @file random_function_mod.f90
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
!> Creates the non uniform or uniform LES grid from MPAS options.
!> This routine is taken from the "numerical recipies"
!------------------------------------------------------------------------------!
 MODULE make_vertical_grid

    USE kinds

    IMPLICIT NONE

    PRIVATE

    PUBLIC construct_vertical_grid_const, construct_vertical_grid_variable,  &
           construct_vertical_grid_variable_mld

    INTEGER(iwp), PUBLIC, SAVE ::  random_iv(32)  !<
    INTEGER(iwp), PUBLIC, SAVE ::  random_iy      !<

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
    SUBROUTINE construct_vertical_grid_const(dz,nzb,nzt,zmidOUT,zedgeOUT)

       IMPLICIT NONE

       integer(iwp), intent(in) :: nzt, nzb
       real(wp), intent(in)  :: dz
       real(wp), intent(out) :: zmidOUT(nzb:nzt+1), zedgeOUT(nzb:nzt+1)
       integer(iwp) :: i

       zedgeOUT(nzt+1) = dz
       zedgeOUT(nzt) = 0.0_wp

       do i = nzt-1,nzb,-1
         zedgeOUT(i) = zedgeOUT(i+1) - dz
       enddo

       do i = nzb+1, nzt+1
         zmidOUT(i) = 0.5*( zedgeOUT(i) + zedgeOUT(i-1) )
       enddo
       zmidOUT(nzb) = zedgeOUT(nzb)

    END SUBROUTINE construct_vertical_grid_const

    SUBROUTINE construct_vertical_grid_variable(depth, dz, nzb, nzt, zmidOUT, zedgeOUT)

       IMPLICIT NONE

       integer(iwp),intent(in) :: nzb, nzt

       real(wp),intent(in) :: depth, dz

       real(wp),intent(out) :: zmidOUT(nzb:nzt+1), zedgeOUT(nzb:nzt+1)

       real(wp) :: z_cntr, z_frst, z_fac1, z_fac2, z_fac

       real(wp) :: zeINV(nzb-1:nzt+1), z_facn, dzArray(nzb:nzt+1)

       integer(iwp) :: iz, il

       ! construct a stretched stretched grid
       z_cntr = depth
       z_frst = -dz
       z_fac1 = z_cntr / z_frst
       z_fac2 = 1.0_wp / REAL(nzt,kind=wp)

       call stretch_factor(z_fac1, z_fac2, z_fac)

       zedgeOUT(nzt+1) = dz
       zedgeOUT(nzt) = 0.0_wp
       zedgeOUT(nzt-1) = -dz
       iz = 2
       do il = nzt-2,nzb,-1
         zedgeOUT(il) = zedgeOUT(nzt-1)*(z_fac**(real(iz,kind=wp)) - 1.0_wp) / (z_fac - 1.0_wp)
         iz = iz + 1
       enddo
       zedgeOUT(nzb) = z_cntr

       do il = nzb+1,nzt+1
         zmidOUT(il) = 0.5*( zedgeOUT(il) + zedgeOUT(il-1) )
       enddo
       zmidOUT(nzb) = zedgeOUT(nzb)

    END SUBROUTINE construct_vertical_grid_variable

    SUBROUTINE stretch_factor(z_fac1, z_fac2, z_fac)

       IMPLICIT NONE

       real(wp), intent(in) :: z_fac1, z_fac2
       real(wp), intent(out) :: z_fac

       integer(iwp) :: knt
       real(wp) :: tol, test, z_facn

       z_fac = 1.10_wp
       tol = 1.0E-10_wp
       test = 10.00_wp
       knt = 0
       do while (test > tol)
         knt = knt + 1
         z_facn = (z_fac1*(z_fac - 1.0_wp) + 1.0_wp)**z_fac2
         test = abs(1.0 - z_facn / z_fac)
         if(knt .gt. 500) THEN
           print *, 'cannot find stretching factor,'
           print *, 'z_fac = ',z_fac, 'z_facn = ',z_facn, 'knt = ',knt
           stop
         ENDIF
         z_fac = z_facn
       enddo

    END SUBROUTINE stretch_factor

    SUBROUTINE construct_vertical_grid_variable_mld(depth, dz, nzb, nzt, zmidOUT, zedgeOUT, MLD)

       IMPLICIT NONE

       integer(iwp),intent(in) :: nzb, nzt

       real(wp),intent(in) :: depth, dz, MLD

       real(wp),intent(out) :: zmidOUT(nzb:nzt+1), zedgeOUT(nzb:nzt+1)

       real(wp) :: z_ht, z_cntr, z_frst, z_fac1, z_fac2, z_fac, z_facbl

       real(wp) :: zeINV(nzb-1:nzt+1), z_facn, dzArray(nzb:nzt+1)

       integer(iwp) :: i, k, iz, il, nzTop, nzBot, nz1, nz2, nzb1

       nzTop = int((nzt-nzb)*2/3)
       nz1 = int(nzTop/2)
       nz2 = nzTop - nz1 - 1
       z_ht = 0.0_wp
       if(-MLD .le. depth) then
         print *, 'ERROR: MLD is bigger than bottom of LES domain, set mixed_layer_refine = .false. and rerun'
       endif
       z_ht = -MLD

       z_fac1 = -z_ht / 2.0 / dz
       z_fac2 = 1.0_wp / REAL(nz1,kind=wp)
       call stretch_factor(z_fac1, z_fac2, z_fac)
       z_facbl = z_fac

       nzb1 = nzt - nzb - nzTop + 1
       z_fac2 = 1.0_wp / REAL(nzb1,kind=wp)
       z_fac1 = -(z_ht - depth) / (-dz)
       call stretch_factor(z_fac1, z_fac2, z_fac)

       zeINV(nzb-1) = dz
       zeINV(nzb) = 0.0_wp
       zeINV(nzb+1) = -dz
       do il = nzb+2,nzb+nz1
          zeINV(il) = zeINV(nzb+1)*(z_facbl**(il) - 1) / (z_facbl - 1)
       enddo

       dzArray(nzb) = -dz
       do il = nzb+1, nz1
         dzArray(il) = zeINV(il-1) - zeINV(il)
       enddo

       zeINV(nzTop-1) = z_ht

       k = nzTop-2
       do il = nzb+1,nz2
         zeINV(k) = zeINV(k+1) + dzArray(il)
         k = k-1
       enddo

       k=2
       zeINV(nzTop) = zeINV(nzTop-1) - dz
       do il = nzTop+1,nzt
         zeINV(il) = zeINV(nzb+1)*(z_fac**(k)-1) / (z_fac - 1) + zeINV(nzTop)
         k = k+1
       enddo

       k = nzt
       do il = nzb,nzt
         zedgeOUT(k) = zeINV(il)
         k = k - 1
       enddo

       zedgeOUT(nzt+1) = dz
       zedgeOUT(nzt) = 0.0_wp
       zedgeOUT(nzt-1) = -dz

       zedgeOUT(nzb) = depth

       do il = nzb+1,nzt+1
         zmidOUT(il) = 0.5*( zedgeOUT(il) + zedgeOUT(il-1) )
       enddo
       zmidOUT(nzb) = zedgeOUT(nzb)

    END SUBROUTINE construct_vertical_grid_variable_mld

 END MODULE make_vertical_grid
