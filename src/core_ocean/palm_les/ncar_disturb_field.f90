!> @file disturb_field.f90
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
! $Id: disturb_field.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 2300 2017-06-29 13:31:14Z raasch
! NEC related code removed
! 
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Modify referenced parameter for disturb_field, instead of nzb_uv_inner, pass
! character to identify the respective grid (u- or v-grid). 
! Set perturbations within topography to zero using flags.
! 
! 2172 2017-03-08 15:55:25Z knoop
! Bugfix removed id_random_array from USE list
! 
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
! 
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable 
! 
! 1425 2014-07-05 10:57:53Z knoop
! bugfix: Parallel random number generator loop: print-statement no needed
! 
! 1400 2014-05-09 14:03:54Z knoop
! Parallel random number generator added
! 
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute 
!
! 1322 2014-03-20 16:38:49Z raasch
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
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! Revision 1.1  1998/02/04 15:40:45  raasch
! Initial revision
!
!
! Description:
! ------------
!> Imposing a random perturbation on a 3D-array.
!> On parallel computers, the random number generator is as well called for all
!> gridpoints of the total domain to ensure, regardless of the number of PEs
!> used, that the elements of the array have the same values in the same
!> order in every case. The perturbation range is steered by dist_range.
!------------------------------------------------------------------------------!
 SUBROUTINE disturb_field

    use arrays_3d,                                                             & 
        ONLY:  u, v 

    USE control_parameters,                                                    &
        ONLY:  dist_nxl, dist_nxr, dist_nyn, dist_nys, dist_range,             &
               disturbance_amplitude, disturbance_created,                     &
               disturbance_level_ind_b, disturbance_level_ind_t, iran,         &
               random_generator, topography, disturb_nblocks
                
    USE cpulog,                                                                &
        ONLY:  cpu_log, log_point
        
    USE grid_variables,                                                        &
        ONLY:  dx, dy

    USE indices,                                                               &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzb_max, &
               nzt, wall_flags_0
        
    USE kinds
    
    USE random_function_mod,                                                   &
        ONLY: random_function
        
    USE random_generator_parallel,                                             &
        ONLY:  random_number_parallel, random_seed_parallel, random_dummy,     &
               seq_random_array

    IMPLICIT NONE

    INTEGER(iwp) ::  i       !< index variable
    INTEGER(iwp) ::  j       !< index variable
    INTEGER(iwp) ::  k       !< index variable
    INTEGER(iwp) ::  nxLocal, nyLocal, iIndex, jIndex, i1, j1, randcount

    REAL(wp) ::  randomnumber(nzb+1:nzt,nysg:nyng,nxlg:nxrg),randmean,randmean2  !<
    
    REAL(wp) ::  meanTarget(nzb+1:nzt), meanCurrent(nzb+1:nzt)
    REAL(wp) ::  dist1(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !<
    REAL(wp) ::  distX(nzb:nzt+1,nysg:nyng,nxlg:nxrg)
    REAL(wp) ::  distY(nzb:nzt+1,nysg:nyng,nxlg:nxrg)

    REAL(wp) ::  vmaxx, vmag

    CALL cpu_log( log_point(20), 'disturb_field', 'start' )
!
!-- Create an additional temporary array and initialize the arrays needed
!-- to store the disturbance
    dist1 = 0.0_wp
    distX = 0.0_wp
    distY = 0.0_wp

!
!-- Create the random perturbation and store it on temporary array
       nxLocal = nxrg - nxlg + 1
       nyLocal = nyng - nysg + 1
       if (mod(nyLocal,disturb_nblocks) .ne. 0 .or. mod(nxLocal,disturb_nblocks) .ne. 0) then
               print *, 'nxLocal or nyLocal is not a multiple of disturb_nblocks, reset'
               print *, nxLocal,nyLocal,disturb_nblocks
               stop
       endif
       nxLocal = nxLocal / disturb_nblocks
       nyLocal = nyLocal / disturb_nblocks
       DO  i = 1, nxLocal
          DO  j = 1, nyLocal
!             DO  k = disturbance_level_ind_b, disturbance_level_ind_t
!                randomnumber(k,j,i) = disturbance_amplitude *                &
!                        ( random_function( iran ) - 0.5_wp )
!             enddo

              
            Do k = disturbance_level_ind_b, disturbance_level_ind_t
             randmean2 =  2.0_wp * ( random_function( iran ) - 0.5_wp )
            ! randmean2 = float (int(randmean * 1.0E10_wp + 0.5_wp)) / 1.0E10_wp

             do i1 = 1, disturb_nblocks
             do j1 = 1, disturb_nblocks
                iIndex = (i-1)*disturb_nblocks + nxlg
                jIndex = (j-1)*disturb_nblocks + nysg
          
                randomnumber(k,jIndex,iIndex) = randmean2 
             enddo   
          !      DO  k = disturbance_level_ind_b, disturbance_level_ind_t
          !         
          !         dist1(k,jIndex,iIndex) = randomnumber(k)
          !      ENDDO
             enddo
            enddo

          ENDDO
       ENDDO
       
       do k = disturbance_level_ind_b, disturbance_level_ind_t
          randmean = 0.0_wp
          randcount = 0
          DO  i = nxlg, nxrg 
              DO  j = nysg, nyng
                  randmean = randmean + randomnumber(k,j,i)
                  randcount = randcount + 1
              enddo
          enddo

          randmean = randmean/randcount

          DO  i = nxlg, nxrg 
              DO  j = nysg, nyng
              randmean2 = randomnumber(k,j,i) - randmean

              dist1(k,j,i) = randmean2 !float (int(randmean2 * 1.0E10_wp + 0.5_wp)) / 1.0E10_wp
              enddo
          enddo

          CALL exchange_horiz( dist1, nbgp )

          !calculate dpsidx (distX) dpsidy (distY)

          do i = nxlg+1, nxrg-1
             do j = nysg+1, nyng-1
                distX(k,j,i) = ( dist1(k,j,i+1) - dist1(k,j,i-1) ) / (2.0_wp*dx)
                distY(k,j,i) = ( dist1(k,j+1,i) - dist1(k,j-1,i) ) / (2.0_wp*dy)
             enddo
          enddo

          vmaxx = 0.0_wp

          CALL exchange_horiz( distX, nbgp )
          CALL exchange_horiz( distY, nbgp )

          do i = nxlg, nxrg
             do j = nysg, nyng
                vmag = sqrt(distX(k,j,i)**2.0_wp + distY(k,j,i)**2.0_wp)
                if(vmag .gt. vmaxx) vmaxx = vmag
             enddo
          enddo

          do i = nxlg, nxrg
             do j =nysg, nyng
                randmean2 = u(k,j,i) - distY(k,j,i)/vmaxx*disturbance_amplitude
                u(k,j,i) = float(int(randmean2 * 1.0E10_wp + 0.5_wp)) / 1.0E10_wp
                randmean2 = v(k,j,i) + distX(k,j,i)/vmaxx*disturbance_amplitude
                v(k,j,i) = float(int(randmean2 * 1.0E10_wp + 0.5_wp)) / 1.0E10_wp 
             enddo
          enddo
       enddo

 !
!-- Set a flag, which indicates that a random perturbation is imposed
    disturbance_created = .TRUE.


    CALL cpu_log( log_point(20), 'disturb_field', 'stop' )


 END SUBROUTINE disturb_field
