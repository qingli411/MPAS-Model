module ocn_adc_mixing_fused

  use mpas_constants
  use mpas_kind_types
  use mpas_log

  use ocn_constants
  use ocn_turbulence
  use ocn_config

  implicit none

  logical :: defineFirst, stopflag

  integer :: i1,i2

  contains

  subroutine swap_time_levels

    i1 = mod(iterCount - 1, 2) + 1
    i2 = mod(iterCount, 2) + 1

  end subroutine swap_time_levels

  subroutine dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
    integer,intent(in) :: nVertLevels, nCells
    real (kind=RKIND), dimension(2,nVertLevels,nCells), intent(in) :: activeTracers
    real (kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: alphaT, betaS
    integer :: iCell, k, ij, i

    real (kind=RKIND), dimension(nVertLevels) :: B, Bup, Bdn
    real (kind=RKIND), dimension(nVertLevels+1) :: tke, BupEdge, BdnEdge, BEdge
    real (kind=RKIND) :: sigav, tumdav, sumdav, Tup, Tdn, Sup, Sdn
    real (kind=RKIND) :: sumv, sumv0
    real (kind=RKIND), parameter :: refT = 15.0_RKIND, refS = 35.0_RKIND, minlen = 0.55_RKIND

    !NOTE: will need to convert to some form of displaced density in the mpas framework soon
    !possibly go back to a more traditional length scale formulation
    do iCell = 1, nCells
       ! level interfaces
       do k = 1, nVertLevels+1
          tke(k) = 0.5_RKIND * (u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
       end do

       ! level centers
       do k = 1, nVertLevels
          B(k) = gravity * (alphaT(k,iCell) * (activeTracers(1,k,iCell) - refT) - &
                            betaS(k,iCell)  * (activeTracers(2,k,iCell) - refS))

          sigav  = 0.5_RKIND * (areaFraction(k,iCell) + areaFraction(k+1,iCell))
          tumdav = 0.5_RKIND * (tumd(k,iCell) + tumd(k+1,iCell))
          sumdav = 0.5_RKIND * (sumd(k,iCell) + sumd(k+1,iCell))

          Tup = activeTracers(1,k,iCell) + sigav*tumdav
          Tdn = activeTracers(1,k,iCell) - (1.0_RKIND - sigav)*tumdav
          Sup = activeTracers(2,k,iCell) + sigav*sumdav
          Sdn = activeTracers(2,k,iCell) - (1.0_RKIND - sigav)*sumdav

          Bup(k) = gravity * (alphaT(k,iCell) * (Tup - refT) - betaS(k,iCell) * (Sup - refS))
          Bdn(k) = gravity * (alphaT(k,iCell) * (Tdn - refT) - betaS(k,iCell) * (Sdn - refS))

          if (k>1) then
             BupEdge(k) = 0.5_RKIND * (Bup(k-1) + Bup(k))
             BdnEdge(k) = 0.5_RKIND * (Bdn(k-1) + Bdn(k))
             BEdge(k)   = 0.5_RKIND * (B(k-1) + B(k))
          end if
       end do

       BdnEdge(nVertLevels+1) = BdnEdge(nVertLevels)
       BupEdge(nVertLevels+1) = BupEdge(nVertLevels)
       BEdge(nVertLevels+1)   = BEdge(nVertLevels)

       BdnEdge(1) = BdnEdge(2)
       BupEdge(1) = BupEdge(2)
       BEdge(1)   = BEdge(2)

       ! leve centers
       do k = 2, nVertLevels
          ! updraft length scale
          sumv  = 0.0_RKIND
          sumv0 = 0.0_RKIND
          ij=k
          lenup(k,iCell) = 0.0_RKIND
          do while(sumv0 <= tke(k) .and. ij > 1)
             ! TODO: ze -> zm?
             sumv = sumv0 + (BEdge(ij) - BupEdge(k)) * (ze(ij-1,iCell)-ze(ij,iCell))
             if (sumv > tke(k)) then
                lenup(k,iCell) = max(minlen, &
                                     lenup(k,iCell) + &
                                     abs(ze(ij-1,iCell)-ze(ij,iCell))*(tke(k)-sumv0)/(sumv-sumv0))
                exit
             else
                lenup(k,iCell) = lenup(k,iCell) + abs(ze(ij-1,iCell)-ze(ij,iCell))
             end if
             sumv0 = sumv
             ij = ij - 1
          end do

          ! downdraft length scale
          sumv  = 0.0_RKIND
          sumv0 = 0.0_RKIND
          ij=k
          lendn(k,iCell) = 0.0_RKIND
          do while(sumv0 <= tke(k) .and. ij < nVertLevels+1)
             sumv = sumv0 + (BdnEdge(k) - BEdge(ij)) * (ze(ij-1,iCell)-ze(ij,iCell))
             if (sumv > tke(k)) then
                lendn(k,iCell) = max(minlen, &
                                     lendn(k,iCell) + &
                                     abs(ze(ij-1,iCell)-ze(ij,iCell))*(tke(k)-sumv0)/(sumv-sumv0))
                exit
             else
                lendn(k,iCell) = lendn(k,iCell) + abs(ze(ij-1,iCell)-ze(ij,iCell))
             end if
             sumv0 = sumv
             ij = ij + 1
          end do

          length(k,iCell) = 2.0_RKIND * lenup(k,iCell) * lendn(k,iCell)  &
                            / (lenup(k,iCell) + lendn(k,iCell))
          ! write(*,*) k, iCell, lenup(k,iCell), lendn(k,iCell)
       enddo
       ! write(*,*) '----'
       length(1,iCell) = minlen
       lenup(1,iCell)  = minlen
       lendn(1,iCell)  = minlen
       length(nVertLevels+1,iCell) = minlen
       lenup(nVertLevels+1,iCell)  = minlen
       lendn(nVertLevels+1,iCell)  = minlen
    enddo


  end subroutine dissipation_lengths2

  subroutine compute_ADC_tends(nCells, nVertLevels, nTracers, dt, activeTracers, uvel, vvel, BVF,  &
                  uwsfc, vwsfc, wtsfc, wssfc, alphaT, betaS, fCell, boundaryLayerDepth)

      integer, intent(in) :: nCells, nVertLevels, nTracers
      real (kind=RKIND), intent(in) :: dt

      real (kind=RKIND), dimension(nTracers,nVertLevels,nCells), intent(inout) :: activeTracers
      real (kind=RKIND), dimension(nVertLevels,nCells), intent(inout) :: uvel, vvel, alphaT, betaS
      real (kind=RKIND), dimension(nCells), intent(in) :: uwsfc, vwsfc, wtsfc, wssfc, fCell
      real (kind=RKIND), dimension(nCells), intent(inout) :: boundaryLayerDepth
      real (kind=RKIND), dimension(nVertLevels,nCells), intent(inout) :: BVF
      integer :: niter, iIter, iCell, k

      real (kind=RKIND) :: dt_small

      real (kind=RKIND) :: Sw, Eav, Dav, sigav, wumdav, tumdav, sumdav
      real (kind=RKIND) :: KE, KEp1, Mcav, lenav, u2av, v2av, w2av, w3av
      real (kind=RKIND) :: w3temp, w3check, Uz, Vz, dz, vscale

      real (kind=RKIND) :: utemp, vtemp
      real (kind=RKIND) :: wb, Cval, diff, wtav, dzmid, Ksps, Sz, Tz
      real (kind=RKIND) :: lareaFraction, wstar, wbsfc
      real (kind=RKIND) :: ustarSquared, wtSumUp, wtSumDn, wsSumUp, wsSumDn

      real (kind=RKIND), dimension(nVertLevels,nCells) :: Swumd
      real (kind=RKIND), dimension(nVertLevels,nCells) :: invTau, invTauAv
      real (kind=RKIND), dimension(nVertLevels,nCells) :: w3tend
      real (kind=RKIND), dimension(nVertLevels+1,nCells) :: wttend, uttend, vttend, wstend, ustend, vstend
      real (kind=RKIND), dimension(nVertLevels+1,nCells) :: w2tend, uwtend, vwtend, u2tend, v2tend, uvtend
      real (kind=RKIND), dimension(nVertLevels+1,nCells) :: KspsUtend, KspsDtend
      real (kind=RKIND), parameter :: smallNum = 1.0e-12_RKIND

      dt_small = config_adc_timestep
      niter = dt / dt_small

      call swap_time_levels

      !on further examination build_diagnostics array can live outside the iter loop
      ! boundary conditions
      do iCell = 1, nCells
        wbsfc = gravity * (alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))
        if (wbsfc > 0.0_RKIND) then
          wstar = abs(wbsfc*boundaryLayerDepth(iCell))**(1.0_RKIND/3.0_RKIND)
        else
          wstar = 0.0_RKIND
        end if

        tumd(1,iCell) = 0.0_RKIND
        wumd(1,iCell) = 0.0_RKIND
        areaFraction(1,iCell) = 0.5_RKIND
        Mc(1,iCell) = 0.0_RKIND
        ! (71) of LR01a
        ! TODO: why different signs for temperature and salinity?
        w2t(1,iCell) = -0.3_RKIND * wstar * wtsfc(iCell)
        w2s(1,iCell) =  0.3_RKIND * wstar * wssfc(iCell)

        length(nVertLevels+1,iCell) = 1e-15_RKIND
        length(1,iCell) = ze(1,iCell) - ze(2,iCell)

        ustarSquared = sqrt(uwsfc(iCell)**2 + vwsfc(iCell)**2)
        do k = 1, 2
          u2(k,1,iCell) = 2.0_RKIND*ustarSquared + 0.3_RKIND*wstar**2
          v2(k,1,iCell) = 2.0_RKIND*ustarSquared + 0.3_RKIND*wstar**2
          uw(k,1,iCell) = -uwsfc(iCell)
          vw(k,1,iCell) = -vwsfc(iCell)
          wt(k,1,iCell) =  wtsfc(iCell)
          ws(k,1,iCell) =  wssfc(iCell)
          KE = 0.5_RKIND*(u2(i1,1,iCell) + v2(i1,1,iCell))
          eps(k,1,iCell) = 2.0_RKIND*KE**1.5_RKIND/(length(1,iCell)+smallNum)
        end do
      end do

      ! time loop
      do iIter = 1, niter

         ! update mass flux arrays
         do iCell = 1, nCells
            do k = 2, nVertLevels
               ! layer interfaces
               w3av = 0.5_RKIND*(w3(i1,k-1,iCell) + w3(i1,k,iCell))

               Sw = w3av / max(w2(i1,k,iCell), smallNum)**1.5_RKIND
               lareaFraction = 0.5_RKIND + 0.5_RKIND * Sw / sqrt(4.0_RKIND + Sw**2)
               lareaFraction = min(max(lareaFraction, 0.01_RKIND), 0.99_RKIND)

               areaFraction(k,iCell) = lareaFraction
               wumd(k,iCell) = sqrt(w2(i1,k,iCell) / (lareaFraction*(1.0_RKIND-lareaFraction)))
               Mc(k,iCell) = lareaFraction*(1.0_RKIND-lareaFraction) * wumd(k,iCell)

               tumd(k,iCell) = wt(i1,k,iCell) / (smallNum + Mc(k,iCell))
               sumd(k,iCell) = ws(i1,k,iCell) / (smallNum + Mc(k,iCell))
            end do
         end do

        do iCell = 1, nCells
          do k = 2, nVertLevels
            ! on u levels
            sigav  = 0.5_RKIND*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
            tumdav = 0.5_RKIND*(tumd(k,iCell) + tumd(k+1,iCell))
            sumdav = 0.5_RKIND*(sumd(k,iCell) + sumd(k+1,iCell))
            wumdav = 0.5_RKIND*(wumd(k,iCell) + wumd(k+1,iCell))
            w2t(k,iCell) = - sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2*tumdav
            w2s(k,iCell) = - sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2*sumdav
            !also use this loop to reset the cliptends for the step
            u2cliptend(k,iCell) = 0.0_RKIND
            v2cliptend(k,iCell) = 0.0_RKIND
            w2cliptend(k,iCell) = 0.0_RKIND
          enddo
        enddo

        ! compute dissipation length scale
        call dissipation_lengths2(nCells, nVertLevels, activeTracers, alphaT, betaS)

        ! w3 tendency
        do iCell = 1, nCells
           ! level centers
           do k = 1, nVertLevels
              Eav  = 0.5_RKIND * (Entrainment(k,iCell) + Entrainment(k+1,iCell))
              Dav  = 0.5_RKIND * (Detrainment(k,iCell) + Detrainment(k+1,iCell))

              u2av = 0.5_RKIND * (u2(i1,k,iCell) + u2(i1,k+1,iCell))
              v2av = 0.5_RKIND * (v2(i1,k,iCell) + v2(i1,k+1,iCell))
              w2av = 0.5_RKIND * (w2(i1,k,iCell) + w2(i1,k+1,iCell))

              sigav   = 0.5_RKIND * (areaFraction(k,iCell) + areaFraction(k+1,iCell))
              wumdav  = 0.5_RKIND * (wumd(k,iCell) + wumd(k+1,iCell))
              tumdav  = 0.5_RKIND * (tumd(k,iCell) + tumd(k+1,iCell))
              sumdav  = 0.5_RKIND * (sumd(k,iCell) + sumd(k+1,iCell))
              lenav   = 0.5_RKIND * (length(k,iCell) + length(k+1,iCell))

              Mcav    = sigav * (1.0_RKIND - sigav) * wumdav

              ! velocity scale
              vscale = sqrt(0.5_RKIND * (u2av+v2av+w2av))
              ! inverse of the time scale
              invTauAv(k,iCell) = vscale / lenAv

              ! layer thickness
              dz = ze(k,iCell) - ze(k+1,iCell)

              ! last term (SPS) in (40) of LR01a, also using (35) of LR01b
              Swumd(k,iCell) = - 2.0_RKIND/3.0_RKIND * ( &
                                1.0_RKIND/(1.0_RKIND-sigav) * &
                                ((1.0_RKIND-areaFraction(k,  iCell))*KspsU(i1,k,  iCell) - &
                                 (1.0_RKIND-areaFraction(k+1,iCell))*KspsU(i1,k+1,iCell)) / dz &
                                - 1.0_RKIND/sigav * &
                                (areaFraction(k,  iCell)*KspsD(i1,k,  iCell) - &
                                 areaFraction(k+1,iCell)*KspsD(i1,k+1,iCell)) / dz )

              ! first and second terms in (52) of LR01a, see also (20) of LR01b
              ! negative sign due to the different definition of sigma
              w3tend1(k,iCell) = - wumdav**3 * (Eav*(3.0_RKIND*sigav - 2.0_RKIND) + &
                                                Dav*(3.0_RKIND*sigav - 1.0_RKIND))
              ! third term in (52) of LR01a
              w3tend2(k,iCell) = - wumdav**3 * (6.0_RKIND*sigav**2 - 6.0_RKIND*sigav + 1.0_RKIND) * &
                                   ( Mc(k,iCell)-Mc(k+1,iCell) ) / dz
              ! fourth term in (52) of LR01a
              w3tend3(k,iCell) = - 1.5_RKIND * (1.0_RKIND - 2.0_RKIND*sigav) * Mcav * wumdav * &
                                   ( (1.0_RKIND - 2.0_RKIND*areaFraction(k,  iCell))*wumd(k,  iCell)**2 &
                                   - (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell))*wumd(k+1,iCell)**2 &
                                   ) / dz
              ! the buoyancy term in the sources and sinks in (52) of LR01a
              ! negative sign due to the different definition of sigma
              w3tend4(k,iCell) = - 3.0_RKIND * (1.0_RKIND - 2.0_RKIND*sigav) * Mcav * wumdav * &
                                   (1.0_RKIND - Cb_w3) * gravity * (alphaT(k,iCell)*tumdav - betaS(k,iCell)*sumdav)
              ! the SPS term in the sources and sinks in (52) of LR01a
              ! negative sign due to the different definition of sigma
              w3tend5(k,iCell) = - 3.0_RKIND * (1.0_RKIND - 2.0_RKIND*sigav) * Mcav * wumdav * Swumd(k,iCell)
              ! w3 tendency
              w3tend(k,iCell) = w3tend1(k,iCell) + w3tend2(k,iCell) + w3tend3(k,iCell) + &
                                w3tend4(k,iCell) + w3tend5(k,iCell)
           end do !nVertLevels
        enddo !nCells

        ! other third order moments
        do iCell = 1, nCells
           ! level centers
           do k = 1, nVertLevels
              ! now get all the downgradient third order moments
              Ksps = 0.5_RKIND * ( (areaFraction(k,iCell) * KspsD(i1,k,iCell) + &
                                   (1.0_RKIND-areaFraction(k,iCell)) * KspsU(i1,k,iCell)) &
                                 + (areaFraction(k+1,iCell) * KspsD(i1,k+1,iCell) + &
                                   (1.0_RKIND-areaFraction(k+1,iCell)) * KspsU(i1,k+1,iCell)))
              KE    = 0.5_RKIND * (u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
              KEp1  = 0.5_RKIND * (u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
              lenav = 0.5_RKIND * (length(k,iCell) + length(k+1,iCell))
              diff = Cmom * sqrt(0.5_RKIND*(KE + KEp1)) * lenav
              dz = ze(k,iCell) - ze(k+1,iCell)
              uw2(k,iCell) = - 1.2_RKIND * diff * (uw(i1,k,iCell) - uw(i1,k+1,iCell)) / dz
              vw2(k,iCell) = - 1.2_RKIND * diff * (vw(i1,k,iCell) - vw(i1,k+1,iCell)) / dz
              u2w(k,iCell) = - 0.6_RKIND * diff * (u2(i1,k,iCell) - u2(i1,k+1,iCell)) / dz
              v2w(k,iCell) = - 0.6_RKIND * diff * (v2(i1,k,iCell) - v2(i1,k+1,iCell)) / dz
              uvw(k,iCell) = - 0.6_RKIND * diff * (uv(i1,k,iCell) - uv(i1,k+1,iCell)) / dz

              diff = Ctherm * sqrt(0.5_RKIND*(KE + KEp1)) * lenav
              uwt(k,iCell) = - diff * (ut(i1,k,iCell) - ut(i1,k+1,iCell)) / dz
              vwt(k,iCell) = - diff * (vt(i1,k,iCell) - vt(i1,k+1,iCell)) / dz
              uws(k,iCell) = - diff * (us(i1,k,iCell) - us(i1,k+1,iCell)) / dz
              vws(k,iCell) = - diff * (vs(i1,k,iCell) - vs(i1,k+1,iCell)) / dz

          enddo !nVertLevels
        enddo !nCells

        ! all second order moment tendencies
        do iCell = 1, nCells
        ! level interfaces
           do k = 2, nVertLevels
              dzmid = (zm(k-1,iCell) - zm(k,iCell))
              dz = ze(k,iCell) - ze(k+1,iCell)
              wb = Mc(k,iCell) * gravity * (alphaT(k,iCell)*tumd(k,iCell) - betaS(k,iCell)*sumd(k,iCell))

              Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dzmid
              Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dzmid
              Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dzmid
              Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dzmid

              KE = 0.5_RKIND * (u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))

              ! velocity scale
              vscale = sqrt(KE)
              ! inverse of the time scale
              invTau(k,iCell) = vscale / length(k,iCell)

              ! See (48) of LR01a for the w2 equation
              ! dissipation by entrainment/detrainment of plumes
              ! (18) of LR01b
              w2tend1(k,iCell) = - wumd(k,iCell)**2 * (Entrainment(k,iCell) + Detrainment(k,iCell))
              ! turbulent transport
              w2tend2(k,iCell) = ( Mc(k-1,iCell) * &
                                 (1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell)) * wumd(k-1,iCell)**2 &
                                   - Mc(k+1,iCell) * &
                                 (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell)) * wumd(k+1,iCell)**2 &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell))
              ! pressure term: slow
              ! the missing w2 term is added when time-stepping the w2 equation
              w2tend3(k,iCell) = Cslow * invTau(k,iCell) * (u2(i1,k,iCell) + v2(i1,k,iCell))/3.0_RKIND
              ! buoyancy term + the buoyancy contribution in the pressure term
              w2tend4(k,iCell) = (2.0_RKIND - 4.0_RKIND/3.0_RKIND*Cb) * wb
              ! pressure term: rapid
              w2tend5(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 - alpha2) * &
                                 (uw(i1,k,iCell)*Uz + vw(i1,k,iCell)*Vz)
              ! subplume scale term
              w2tend6(k,iCell) = Mc(k,iCell) * (Swumd(k-1,iCell) + Swumd(k,iCell))
              ! w2 tendency
              w2tend(k,iCell) = w2tend1(k,iCell) + w2tend2(k,iCell) + w2tend3(k,iCell) + &
                                w2tend4(k,iCell) + w2tend5(k,iCell) + w2tend6(k,iCell)

              ! See (49) of LR01a for the wt equation
              ! sources/sinks due to entrainment/detrainment
              wttend1(k,iCell) = - (Entrainment(k,iCell) + Detrainment(k,iCell)) * &
                                   wumd(k,iCell) * tumd(k,iCell)
              ! transport term
              wttend2(k,iCell) = ( (1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell)) * &
                                   wumd(k-1,iCell) * tumd(k-1,iCell) * Mc(k-1,iCell) &
                                 - (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell)) * &
                                   wumd(k+1,iCell) * tumd(k+1,iCell) * Mc(k+1,iCell) &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell))
              ! gradient production term
              wttend3(k,iCell) = - Mc(k,iCell) * wumd(k,iCell) * Tz
              ! buoyancy term
              wttend4(k,iCell) = (1.0_RKIND - Cb_tracer) * gravity * &
                                 areaFraction(k,iCell) * (1.0_RKIND - areaFraction(k,iCell)) * &
                                 (alphaT(k,iCell)*tumd(k,iCell)**2 - &
                                  betaS(k,iCell)*tumd(k,iCell)*sumd(k,iCell))
              ! pressure term: only the rapid term here
              ! the slow term is added when time-stepping wt equation
              wttend5(k,iCell) = (alpha1_tracer - alpha2_tracer) / 2.0_RKIND * &
                                 (ut(i1,k,iCell)*Uz + vt(i1,k,iCell)*Vz)
              ! sources/sinks due to subplume terms
              wttend6(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell)) * &
                                 tumd(k,iCell) * 0.5_RKIND * (Swumd(k-1,iCell) + Swumd(k,iCell)) &
                               - Mc(k,iCell) * ( 1.0_RKIND / (1.0_RKIND - areaFraction(k,iCell)) * &
                                 ( (1.0_RKIND - areaFraction(k-1,iCell)) * wt_spsU(k-1,iCell) &
                                 - (1.0_RKIND - areaFraction(k+1,iCell)) * wt_spsU(k+1,iCell) &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell)) &
                                               - 1.0_RKIND / areaFraction(k,iCell) * &
                                 ( areaFraction(k-1,iCell) * wt_spsD(k-1,iCell) &
                                 - areaFraction(k+1,iCell) * wt_spsD(k+1,iCell) &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell)) &
                                               )
              ! wt tendency
              wttend(k,iCell) = wttend1(k,iCell) + wttend2(k,iCell) + wttend3(k,iCell) + &
                                wttend4(k,iCell) + wttend5(k,iCell) + wttend6(k,iCell)

              ! See (49) of LR01a for the wt equation
              ! sources/sinks due to entrainment/detrainment
              wstend1(k,iCell) = - (Entrainment(k,iCell) + Detrainment(k,iCell)) * &
                                   wumd(k,iCell)*sumd(k,iCell)
              ! transport term
              wstend2(k,iCell) = ( (1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell)) * &
                                   wumd(k-1,iCell) * sumd(k-1,iCell) * Mc(k-1,iCell) &
                                 - (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell)) * &
                                   wumd(k+1,iCell) * sumd(k+1,iCell) * Mc(k+1,iCell) &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell))
              ! gradient production term
              wstend2(k,iCell) = - Mc(k,iCell) * wumd(k,iCell) * Sz
              ! buoyancy term
              wstend3(k,iCell) = (1.0_RKIND - Cb_tracer) * gravity * &
                                 areaFraction(k,iCell) * (1.0_RKIND - areaFraction(k,iCell)) * &
                                 (alphaT(k,iCell)*tumd(k,iCell)*sumd(k,iCell) - &
                                  betaS(k,iCell)*sumd(k,iCell)**2)
              ! pressure term: only the rapid term here
              ! the slow term is added when time-stepping ws equation
              wstend4(k,iCell) = (alpha1_tracer - alpha2_tracer) / 2.0_RKIND * &
                                 (us(i1,k,iCell)*Uz + vs(i1,k,iCell)*Vz)
              ! sources/sinks due to subplume terms
              wstend6(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell)) * &
                                 sumd(k,iCell) * 0.5_RKIND * (Swumd(k-1,iCell) + Swumd(k,iCell)) &
                               - Mc(k,iCell) * ( 1.0_RKIND / (1.0_RKIND - areaFraction(k,iCell)) * &
                                 ( (1.0_RKIND - areaFraction(k-1,iCell)) * ws_spsU(k-1,iCell) &
                                 - (1.0_RKIND - areaFraction(k+1,iCell)) * ws_spsU(k+1,iCell) &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell)) &
                                               - 1.0_RKIND / areaFraction(k,iCell) * &
                                 ( areaFraction(k-1,iCell) * ws_spsD(k-1,iCell) &
                                 - areaFraction(k+1,iCell) * ws_spsD(k+1,iCell) &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell)) &
                                               )
              ! ws tendency
              wstend(k,iCell) = wstend1(k,iCell) + wstend2(k,iCell) + wstend3(k,iCell) + &
                                wstend4(k,iCell) + wstend5(k,iCell) + wstend6(k,iCell)

              ! uw equation
              ! transport term
              uwtend1(k,iCell) = - (uw2(k-1,iCell) - uw2(k,iCell)) / dzmid
              ! shear production + rapid pressure term: dU/dz
              ! note that the slow pressure term is in time stepping
              uwtend2(k,iCell) = 0.5_RKIND * ( (alpha0 - 4.0_RKIND/3.0_RKIND*alpha1) * KE + &
                                               (alpha1 - alpha2) * u2(i1,k,iCell) + &
                                               (alpha1 + alpha2 - 2.0_RKIND) * w2(i1,k,iCell) &
                                             ) * Uz
              ! rapid pressure term: dV/dz
              uwtend3(k,iCell) = 0.5_RKIND * (alpha1 - alpha2) * uv(i1,k,iCell) * Vz
              ! buoyancy term
              uwtend4(k,iCell) = (1.0_RKIND - Cb) * gravity * &
                                 (alphaT(k,iCell)*ut(i1,k,iCell) - betaS(k,iCell)*us(i1,k,iCell))
              ! uw tendency
              uwtend(k,iCell) = uwtend1(k,iCell) + uwtend2(k,iCell) + uwtend3(k,iCell) + &
                                uwtend4(k,iCell)

              ! vw equation
              ! transport term
              vwtend1(k,iCell) = - (vw2(k-1,iCell) - vw2(k,iCell)) / dzmid
              ! rapid pressure term: dU/dz
              ! note that the slow pressure term is in time stepping
              vwtend2(k,iCell) = 0.5_RKIND * (alpha1 - alpha2) * uv(i1,k,iCell) * Uz
              ! shear production + rapid pressure term: dV/dz
              vwtend3(k,iCell) = 0.5_RKIND * ( (alpha0 - 4.0_RKIND/3.0_RKIND*alpha1) * KE + &
                                               (alpha1 - alpha2) * v2(i1,k,iCell) + &
                                               (alpha1 + alpha2 - 2.0_RKIND) * w2(i1,k,iCell) &
                                             ) * Vz
              ! buoyancy term
              vwtend4(k,iCell) = (1.0_RKIND - Cb) * gravity * &
                                 (alphaT(k,iCell)*vt(i1,k,iCell) - betaS(k,iCell)*vs(i1,k,iCell))
              ! vw tendency
              vwtend(k,iCell) = vwtend1(k,iCell) + vwtend2(k,iCell) + vwtend3(k,iCell) + &
                                vwtend4(k,iCell)

              ! uv tendency
              uvtend(k,iCell) = &
              ! transport
                               - (uvw(k-1,iCell) - uvw(k,iCell)) / dz &
              ! shear production + rapid pressure term
              ! note that the slow pressure term is in time stepping
                               - (1.0_RKIND - 0.5_RKIND*(alpha1+alpha2)) * &
                                 (uw(i1,k,iCell)*Vz + vw(i1,k,iCell)*Uz)

              ! u2 equation
              ! transport
              u2tend1(k,iCell) = - (u2w(k-1,iCell) - u2w(k,iCell)) / dzmid
              ! shear production + rapid pressure term: dU/dz
              u2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 + alpha2 - 2.0_RKIND) * uw(i1,k,iCell) * Uz
              ! rapid pressure term: dV/dz
              u2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1 * vw(i1,k,iCell) * Vz
              ! buoyancy term
              u2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND * Cb * wb
              ! dissipation + return to isotropy
              ! note that the term proportional to u2 is in time stepping
              u2tend5(k,iCell) = - 2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
                                 Cslow * invTau(k,iCell) * (v2(i1,k,iCell) + w2(i1,k,iCell)) / 3.0_RKIND
              ! u2 tendency
              u2tend(k,iCell) = u2tend1(k,iCell) + u2tend2(k,iCell) + u2tend3(k,iCell) + &
                                u2tend4(k,iCell) + u2tend5(k,iCell)

              ! v2 equation
              ! transport
              v2tend1(k,iCell) = - (v2w(k-1,iCell) - v2w(k,iCell)) / dzmid
              ! rapid pressure term: dU/dz
              v2tend2(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1*uw(i1,k,iCell)*Uz
              ! shear production + rapid pressure term: dV/dz
              v2tend3(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 + alpha2 - 2.0_RKIND) * vw(i1,k,iCell) * Vz
              ! buoyancy term
              v2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND * Cb * wb
              ! dissipation + return to isotropy
              ! note that the term proportional to v2 is in time stepping
              v2tend5(k,iCell) = - 2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
                                 Cslow * invTau(k,iCell) * (u2(i1,k,iCell) + w2(i1,k,iCell)) / 3.0_RKIND
              ! v2 tendency
              v2tend(k,iCell) = v2tend1(k,iCell) + v2tend2(k,iCell) + v2tend3(k,iCell) + &
                                v2tend4(k,iCell) + v2tend5(k,iCell)

              ! ut tendency
              uttend(k,iCell) = &
              ! transport
                               - (uwt(k-1,iCell) - uwt(k,iCell))/dz &
              ! gradient term
                               - uw(i1,k,iCell) * Tz &
              ! rapid pressure term
                               - (1.0_RKIND - 0.5_RKIND*(alpha1_tracer+alpha2_tracer)) * wt(i1,k,iCell) * Uz

              ! vt tendency
              vttend(k,iCell) = &
              ! transport
                               - (vwt(k-1,iCell) - vwt(k,iCell)) / dz &
              ! gradient term
                               - vw(i1,k,iCell) * Tz &
              ! rapid pressure term
                               - (1.0_RKIND - 0.5_RKIND*(alpha1_tracer+alpha2_tracer)) * wt(i1,k,iCell) * Vz

              ! us tendency
              ustend(k,iCell) = &
              ! transport
                               - (uws(k-1,iCell) - uws(k,iCell)) / dz &
              ! gradient term
                               - uw(i1,k,iCell) * Sz &
              ! rapid pressure term
                               - (1.0_RKIND - 0.5_RKIND*(alpha1_tracer+alpha2_tracer)) * ws(i1,k,iCell) * Uz

              ! vs tendency
              vstend(k,iCell) = &
              ! transport
                               - (vws(k-1,iCell) - vws(k,iCell)) / dz &
              ! gradient term
                               - vw(i1,k,iCell) * Sz &
              ! rapid pressure term
                               - (1.0_RKIND - 0.5_RKIND*(alpha1_tracer+alpha2_tracer)) * ws(i1,k,iCell) * Vz

              ! output tracer statistics
              t2(i2,k,iCell) = tumd(k,iCell)**2 * areaFraction(k,iCell) * (1.0_RKIND-areaFraction(k,iCell))
              s2(i2,k,iCell) = sumd(k,iCell)**2 * areaFraction(k,iCell) * (1.0_RKIND-areaFraction(k,iCell))
              ts(i2,k,iCell) = tumd(k,iCell) * sumd(k,iCell) * areaFraction(k,iCell) * (1.0_RKIND-areaFraction(k,iCell))

           end do ! nVertLevels
        end do ! nCells

        ! subplume scale fluxes
        do iCell = 1, nCells
           ! level interfaces
           do k = 2, nVertLevels

              dzmid = (zm(k-1,iCell) - zm(k,iCell))
              Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dzmid
              Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dzmid
              Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dzmid
              Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dzmid

              ! (32) and (33) of LR01b
              if (BVF(k,iCell) <= 0.0_RKIND) then
                lenspsU(k,iCell) = dzmid
                lenspsD(k,iCell) = dzmid
              else
                lenspsU(k,iCell) = min(dzmid, 0.76_RKIND*sqrt(KspsU(i1,k,iCell)/BVF(k,iCell)))
                lenspsD(k,iCell) = min(dzmid, 0.76_RKIND*sqrt(KspsD(i1,k,iCell)/BVF(k,iCell)))
              endif

              ! (34) of LR01b
              KmU(k,iCell) = 0.1_RKIND * lenspsU(k,iCell) * sqrt(KspsU(i1,k,iCell))
              KhU(k,iCell) = (1.0_RKIND + 2.0_RKIND*lenspsU(k,iCell)/dzmid) * KmU(k,iCell)
              wt_spsU(k,iCell) = - KhU(k,iCell) * Tz
              ws_spsU(k,iCell) = - KhU(k,iCell) * Sz

              KmD(k,iCell) = 0.1_RKIND * lenspsD(k,iCell) * sqrt( KspsD(i1,k,iCell) )
              KhD(k,iCell) = (1.0_RKIND + 2.0_RKIND*lenspsD(k,iCell)/dzmid) * KmD(k,iCell)
              wt_spsD(k,iCell) = - KhD(k,iCell) * Tz
              ws_spsD(k,iCell) = - KhD(k,iCell) * Sz

              !change length scale to lenup and lendown
              ! (11) of LR01b
              Entrainment(k,iCell) = Cww_E * areaFraction(k,iCell) * (1.0_RKIND-areaFraction(k,iCell)) * &
                                     Mc(k,iCell) / lenup(k,iCell)
              Detrainment(k,iCell) = Cww_D * areaFraction(k,iCell) * (1.0_RKIND-areaFraction(k,iCell)) * &
                                     Mc(k,iCell) / lendn(k,iCell)

              ! (31) of LR01b
              if (k==2) then
                Cval = 3.96_RKIND
              else
                Cval = 0.19_RKIND+0.51_RKIND*lenspsU(k,iCell)/dzmid
              endif

              ! subplume scale TKE equation
              ! (28) of LR01b
              KspsUtend(k,iCell) = &
              ! buoyancy
                         gravity * (alphaT(k,iCell)*wt_spsU(k,iCell) - betaS(k,iCell)*ws_spsU(k,iCell)) &
              ! transport of SPS TKE and pressure term: (29) of LR01b
                       + ((KmU(k-1,iCell) + KmU(k,iCell)) * (KspsU(i1,k-1,iCell) - KspsU(i1,k,iCell)) / &
                          (ze(k-1,iCell) - ze(k,iCell)) &
                        - (KmU(k,iCell) + KmU(k+1,iCell)) * (KspsU(i1,k,iCell) - KspsU(i1,k+1,iCell)) / &
                          (ze(k,iCell) - ze(k+1,iCell)) &
                         ) / dzmid &
              ! dissipation: (30) of LR01b
                       - Cval*KspsU(i1,k,iCell)**1.5_RKIND / lenspsU(k,iCell) &
              ! large scale dissipation
                       + eps(i1,k,iCell) / (2.0_RKIND*(1.0_RKIND-areaFraction(k,iCell))) &
              ! shear production
                       + KmU(k,iCell) * (Uz**2 + Vz**2)

              if (k==2) then
                Cval = 3.96_RKIND
              else
                Cval = 0.19_RKIND+0.51_RKIND*lenspsD(k,iCell)/dzmid
              endif

              ! subplume scale TKE equation
              ! (28) of LR01b
              KspsDtend(k,iCell) = &
              ! buoyancy
                         gravity * (alphaT(k,iCell)*wt_spsD(k,iCell) - betaS(k,iCell)*ws_spsD(k,iCell)) &
              ! transport of SPS TKE and pressure term: (29) of LR01b
                       + ((KmD(k-1,iCell) + KmD(k,iCell)) * (KspsD(i1,k-1,iCell) - KspsD(i1,k,iCell)) / &
                          (ze(k-1,iCell) - ze(k,iCell)) &
                        - (KmD(k,iCell) + KmD(k+1,iCell)) * (KspsD(i1,k,iCell) - KspsD(i1,k+1,iCell)) / &
                          (ze(k,iCell) - ze(k+1,iCell)) &
                         ) / dzmid &
              ! dissipation: (30) of LR01b
                       - Cval*KspsD(i1,k,iCell)**1.5_RKIND / lenspsD(k,iCell) &
              ! large scale dissipation
                       + eps(i1,k,iCell) / (2.0_RKIND*areaFraction(k,iCell)) &
              ! shear production
                       + KmD(k,iCell) * (Uz**2 + Vz**2)

           enddo ! nVertLevels
        enddo ! nCells

        !In this step we update second moments except w3 which needs updated w2
        do iCell = 1, nCells
           do k = 2, nVertLevels
              !update second order moment tendency here
              w2(i2,k,iCell) = (w2(i1,k,iCell) + dt_small*w2tend(k,iCell)) / &
                               (1.0_RKIND + dt_small*Cslow*invTau(k,iCell)*2.0_RKIND/3.0_RKIND)
              if (w2(i2,k,iCell) < epsilon) then
                 w2cliptend(k,iCell) = epsilon-w2(i2,k,iCell)
                 w2(i2,k,iCell) = epsilon
              endif

              if (abs(w2(i2,k,iCell)) > 1.0_RKIND) then
                call mpas_log_write("ERROR: w2 out of range, w2 = $r, location = $i, $i", &
                                    MPAS_LOG_CRIT,realArgs=(/w2(i2,k,iCell)/),intArgs=(/k,iCell/))
              endif

            u2(i2,k,iCell) = (u2(i1,k,iCell) + dt_small*u2tend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow*invTau(k,iCell)*2.0_RKIND/3.0_RKIND)
            if (u2(i2,k,iCell) < 0.0_RKIND) then
              u2cliptend(k,iCell) = -u2(i2,k,iCell)
              u2(i2,k,iCell) = 0.0_RKIND
            endif

            v2(i2,k,iCell) = (v2(i1,k,iCell) + dt_small*v2tend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow*invTau(k,iCell)*2.0_RKIND/3.0_RKIND)
            if (v2(i2,k,iCell) < 0.0_RKIND) then
              v2cliptend(k,iCell) = -v2(i2,k,iCell)
              v2(i2,k,iCell) = 0.0_RKIND
            endif

            uw(i2,k,iCell) = (uw(i1,k,iCell) + dt_small*uwtend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow*invTau(k,iCell))
            vw(i2,k,iCell) = (vw(i1,k,iCell) + dt_small*vwtend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow*invTau(k,iCell))
            uv(i2,k,iCell) = (uv(i1,k,iCell) + dt_small*uvtend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow*invTau(k,iCell))
            ut(i2,k,iCell) = (ut(i1,k,iCell) + dt_small*uttend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow_tracer*invTau(k,iCell))
            vt(i2,k,iCell) = (vt(i1,k,iCell) + dt_small*vttend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow_tracer*invTau(k,iCell))
            wt(i2,k,iCell) = (wt(i1,k,iCell) + dt_small*wttend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow_tracer*invTau(k,iCell))
            us(i2,k,iCell) = (us(i1,k,iCell) + dt_small*ustend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow_tracer*invTau(k,iCell))
            vs(i2,k,iCell) = (vs(i1,k,iCell) + dt_small*vstend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow_tracer*invTau(k,iCell))
            ws(i2,k,iCell) = (ws(i1,k,iCell) + dt_small*wstend(k,iCell)) / &
                             (1.0_RKIND + dt_small*Cslow_tracer*invTau(k,iCell))
            if (abs(wt(i2,k,iCell)) > 1.0_RKIND) then
              call mpas_log_write("ERROR: wt out of range, wt = $r, location = $i, $i", &
                                  MPAS_LOG_CRIT,realArgs=(/wt(i2,k,iCell)/),intArgs=(/k,iCell/))
            endif

            if (abs(ws(i2,k,iCell)) > 1.0_RKIND) then
              call mpas_log_write("ERROR: ws out of range, ws = $r, location = $i, $i", &
                                  MPAS_LOG_CRIT,realArgs=(/ws(i2,k,iCell)/),intArgs=(/k,iCell/))
            endif

            if (abs(u2(i2,k,iCell)) > 1.0_RKIND) then
              call mpas_log_write("ERROR: u2 out of range, u2 = $r, location = $i, $i", &
                                  MPAS_LOG_CRIT,realArgs=(/u2(i2,k,iCell)/),intArgs=(/k,iCell/))
            endif

            if (abs(v2(i2,k,iCell)) > 1.0_RKIND) then
              call mpas_log_write("ERROR: v2 out of range, v2 = $r, location = $i, $i", &
                                  MPAS_LOG_CRIT,realArgs=(/v2(i2,k,iCell)/),intArgs=(/k,iCell/))
            endif

            KE = 0.5_RKIND * (u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))
            eps(i2,k,iCell) = KE**1.5_RKIND / 16.6_RKIND * length(k,iCell)
            KspsU(i2,k,iCell) = max(epsilon, KspsU(i1,k,iCell) + dt_small*KspsUtend(k,iCell))
            KspsD(i2,k,iCell) = max(epsilon, KspsD(i1,k,iCell) + dt_small*KspsDtend(k,iCell))

          enddo !nVertLevels
        enddo !nCells for second order moment tendencies

        ! update the third order w3 and mean fields
        do iCell = 1,nCells
          do k = 1,nVertLevels
            w3check = - (w2(i2,k,iCell) + w2(i2,k+1,iCell))**1.5_RKIND
            w3temp  = (w3(i1,k,iCell) + dt_small*w3tend(k,iCell)) / &
                      (1.0_RKIND + dt_small*Cslow_w3*invTauAv(k,iCell))
            w3(i2,k,iCell) = max(w3temp, w3check)

            if(abs(w3(i2,k,iCell)) > 1.0_RKIND) then
              call mpas_log_write("ERROR: w3 out of range, w3 = $r, location = $i, $i", &
                                  MPAS_LOG_CRIT,realArgs=(/w3(i2,k,iCell)/),intArgs=(/k,iCell/))
            endif
          enddo
        enddo
        iterCount = iterCount + 1
        call swap_time_levels !probably a better way to do this, with supercycling maybe
!  enddo !end iteration loop -- substepping is done.
    !now that substepping is done, apply computed fluxes to update mean fields.
    !you can collapse this loop too.
        do iCell = 1,nCells
          do k = 1,nVertLevels
            utemp = uvel(k,iCell)
            vtemp = vvel(k,iCell)
            uvel(k,iCell) = uvel(k,iCell) - dt_small*(uw(i1,k,iCell) - uw(i1,k+1,iCell)) /  &
                            (ze(k,iCell) - ze(k+1,iCell)) !+ dt_small*fCell(iCell)*vtemp

            vvel(k,iCell) = vvel(k,iCell) - dt_small*(vw(i1,k,iCell) - vw(i1,k+1,iCell)) /  &
                            (ze(k,iCell) - ze(k+1,iCell)) !- dt_small*fCell(iCell)*utemp

            wtSumUp = wt(i1,k,iCell) + areaFraction(k,iCell)*wt_spsD(k,iCell) + &
                      (1.0_RKIND - areaFraction(k,iCell))*wt_spsU(k,iCell)
            wtSumDn = wt(i1,k+1,iCell) + areaFraction(k+1,iCell)*wt_spsD(k+1,iCell) + &
                      (1.0_RKIND - areaFraction(k+1,iCell))*wt_spsU(k+1,iCell)
            wsSumUp = ws(i1,k,iCell) + areaFraction(k,iCell)*ws_spsD(k,iCell) + &
                      (1.0_RKIND - areaFraction(k,iCell))*ws_spsU(k,iCell)
            wsSumDn = ws(i1,k+1,iCell) + areaFraction(k+1,iCell)*ws_spsD(k+1,iCell) + &
                      (1.0_RKIND - areaFraction(k+1,iCell))*ws_spsU(k+1,iCell)

            activeTracers(1,k,iCell) = activeTracers(1,k,iCell) - dt_small*(wtSumUp - wtSumDn) / &
                                       (ze(k,iCell) - ze(k+1,iCell))
            activeTracers(2,k,iCell) = activeTracers(2,k,iCell) - dt_small*(wsSumUp - wsSumDn) / &
                                       (ze(k,iCell) - ze(k+1,iCell))
          enddo
        enddo
      enddo

  end subroutine compute_ADC_tends

end module ocn_adc_mixing_fused
