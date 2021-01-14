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
    real (kind=RKIND), dimension(nVertLevels+1) :: KEsps, tke, BupEdge, BdnEdge
    real (kind=RKIND) :: sigav, tumdav, sumdav, Tup, Tdn, Sup, Sdn
    real (kind=RKIND) :: s1, z1, zV, sumv, minlen
    real (kind=RKIND), parameter :: refT = 15.0_RKIND, refS = 35.0_RKIND

    !NOTE: will need to convert to some form of displaced density in the mpas framework soon
    !possibly go back to a more traditional length scale formulation
    do iCell = 1,nCells
       KEsps(:) = areaFraction(:,iCell)*KspsU(i1,:,iCell) + (1.0_RKIND - areaFraction(:,iCell))*KspsD(i1,:,iCell)
       tke(:) = 0.5_RKIND*(u2(i1,:,iCell) + v2(i1,:,iCell) + w2(i1,:,iCell)) ! + KEsps(:)

       do k = 1,nVertLevels
          B(k) = gravity * (alphaT(k,iCell) * (activeTracers(1,k,iCell) - refT) - &
                            betaS(k,iCell)  * (activeTracers(2,k,iCell) - refS))

          sigav  = 0.5_RKIND * (areaFraction(k,iCell) + areaFraction(k+1,iCell))
          tumdav = 0.5_RKIND * (tumd(k,iCell) + tumd(k+1,iCell))
          sumdav = 0.5_RKIND * (sumd(k,iCell) + sumd(k+1,iCell))

          Tup = activeTracers(1,k,iCell) + (1.0_RKIND - sigav)*tumdav
          Tdn = activeTracers(1,k,iCell) - sigav*tumdav
          Sup = activeTracers(2,k,iCell) + (1.0_RKIND - sigav)*sumdav
          Sdn = activeTracers(2,k,iCell) - sigav*sumdav

          Bup(k) = gravity * (alphaT(k,iCell) * (Tup - refT) - betaS(k,iCell) * (Sup - refS))
          Bdn(k) = gravity * (alphaT(k,iCell) * (Tdn - refT) - betaS(k,iCell) * (Sdn - refS))

          if(k>1) THEN
             BupEdge(k) = 0.5_RKIND * (Bup(k-1) + Bup(k))
             BdnEdge(k) = 0.5_RKIND * (Bdn(k-1) + Bdn(k))
          endif
       enddo

       BdnEdge(nVertLevels+1) = BdnEdge(nVertLevels)
       BupEdge(nVertLevels+1) = BupEdge(nVertLevels)

       BdnEdge(1) = BdnEdge(2)
       BupEdge(1) = BupEdge(2)

       do k=2,nVertLevels

          sumv = 0.0_RKIND
          ij=k
          lenup(k,iCell) = 0.0_RKIND
          do while(sumv <= tke(k) .and. ij < nVertLevels+1)
             sumv = sumv + (BupEdge(k) - Bup(ij))*(ze(ij-1,iCell)-ze(ij,iCell))
             lenup(k,iCell) =  lenup(k,iCell) + abs(ze(ij-1,iCell)-ze(ij,iCell))
             ij = ij + 1

             if(sumv > tke(k)) THEN
                ij = ij - 1
                s1 = sumv
                z1 = ze(ij,iCell)
                zV = ze(ij-1,iCell)

                sumv = sumv - (BupEdge(k) - Bup(ij))*(ze(ij-1,iCell)-ze(ij,iCell))
                lenup(k,iCell) = lenup(k,iCell) - abs(ze(ij,iCell)-ze(ij+1,iCell))
                lenup(k,iCell) = max(0.55_RKIND,lenup(k,iCell) + abs((z1-zV)/(s1 - sumv)*(tke(k)-sumv)))
            !    if(Bup(k-1) - Bup(k) < 0) then
            !            minlen = abs(zm(k-1,i) - zm(k,i))
            !    else
            !            minlen = 0.5_RKIND
            !    endif
            !    if(k==5) then
            !      print *, 'len pre = ',lenup(k,i)
            !    endif
!
!                lenup(k,i) = max(minlen, lenup(k,i) + tke(k)/(1.0E-15 + sumv - (BupEdge(k) -         &
!                                Bup(ij))))
!                              if(k==5) then
!                                print *, lenup(k,i),minlen, tke(k)
!                              endif
                              exit
             endif

         end do

        !find lendown
        sumv = 0.0_RKIND
        ij=k
        lendn(k,iCell) = 0.0_RKIND
        do while(sumv <= tke(k) .and. ij>1)
           sumv = sumv - (BdnEdge(k) - Bdn(ij-1))*(ze(ij-1,iCell)-ze(ij,iCell))
           lendn(k,iCell) = lendn(k,iCell) + abs(ze(ij-1,iCell)-ze(ij,iCell))
           ij = ij - 1

           if(sumv > tke(k)) THEN
              ij = ij + 1
              s1 = sumv
              z1 = ze(ij,iCell)
              zV = ze(ij-1,iCell)
              sumv = sumv + (BdnEdge(k) - Bdn(ij-1))*(ze(ij-1,iCell)-ze(ij,iCell))
              lendn(k,iCell) = lendn(k,iCell) - abs(ze(ij-1,iCell)-ze(ij,iCell))
              lendn(k,iCell) = max(0.55_RKIND,lendn(k,iCell) + ((zV-z1)/(s1 - sumv)*(tke(k) - sumv)))

           exit
           endif
        enddo

        length(k,iCell) = (2.0_RKIND*lenup(k,iCell)*lendn(k,iCell))  &
          / (lenup(k,iCell) + lendn(k,iCell))
      enddo
   enddo

   length(1,iCell) = 0.55_RKIND
   length(nVertLevels+1,iCell) = 0.55_RKIND

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
      integer :: niter, iIter,iCell,k, supercycle

      real (kind=RKIND) :: dt_small

      real (kind=RKIND), dimension(nVertLevels+1,nCells) :: w2tend
      real (kind=RKIND), dimension(nVertLevels,nCells) :: w3tend

      real (kind=RKIND) :: Sw, St, Ss, Eav, Dav, sigav, wumdav, tumdav, sumdav
      real (kind=RKIND) :: KspsUav, KspsDav, KE, Mcav, lenav,u2av,v2av,w2av
      real (kind=RKIND) :: w3temp, w3check, mval, KEsps, Uz, Vz, dz, vscale

      real (kind=RKIND) :: invLen, l, len1, len2, lenmax, integrandTop, integrandBot
      real (kind=RKIND) :: len0, len2_1, len2_2, sfcBuoy, lengthT, bvfT
      real (kind=RKIND) :: KEm1, KEp1, tauUP, tauDN, tomUP, tomDN
      real (kind=RKIND) :: tauM1, tau, tauP1, tauAV, utemp, vtemp
      real (kind=RKIND) :: B, Cval, diff, wtav, dzmid, Ksps, Sz, Tz, w4k, w4kp1, w2k, w2kp1
      real (kind=RKIND) :: lareaFraction, wstar, wbsfc, w3av, tempMoment
      real (kind=RKIND) :: sfcFrictionVelocity, wtSumUp, wtSumDn, wsSumUp, wsSumDn

      real (kind=RKIND), dimension(nVertLevels,nCells) :: Swumd
      real (kind=RKIND), dimension(nVertLevels,nCells) :: invTaupt, invTaups, invTaupv
      real (kind=RKIND), dimension(nVertLevels+1,nCells) :: wttend, wstend, uwtend, vwtend, u2tend, v2tend, &
        uttend, vttend, ustend, vstend, uvtend
      real (kind=RKIND), dimension(nVertLevels,nCells) :: epstend
      real (kind=RKIND), dimension(nVertLevels+1,nCells) :: KspsUtend, KspsDtend
      real (kind=RKIND), parameter :: smallNum = 1.0e-12_RKIND

      dt_small = config_adc_timestep
      niter = dt / dt_small

      call swap_time_levels

      !on further examination build_diagnostics array can live outside the iter loop
      ! boundary conditions
      do iCell = 1,nCells
        wbsfc = gravity * (alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))
        if (wbsfc > 0.0_RKIND) then
          wstar = abs(wbsfc*boundaryLayerDepth(iCell))**(1.0_RKIND/3.0_RKIND)
        else
          wstar = 0.0_RKIND
        endif

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

        sfcFrictionVelocity = sqrt(uwsfc(iCell)**2 + vwsfc(iCell)**2)
        do k = 1,2
          u2(k,1,iCell) = 2.0_RKIND*sfcFrictionVelocity + 0.3_RKIND*wstar**2
          v2(k,1,iCell) = 2.0_RKIND*sfcFrictionVelocity + 0.3_RKIND*wstar**2
          uw(k,1,iCell) = -uwsfc(iCell)
          vw(k,1,iCell) = -vwsfc(iCell)
          wt(k,1,iCell) =  wtsfc(iCell)
          ws(k,1,iCell) =  wssfc(iCell)
          KE = 0.5*(u2(i1,1,iCell) + v2(i1,1,iCell))
          !eps(k,1,iCell) = KE**1.5/(0.5*(ze(1,iCell) - ze(2,iCell) + 1.0E-10_RKIND))
          ! TODO
          eps(k,1,iCell) = 2.0_RKIND*KE**1.5_RKIND/(length(1,iCell)+smallNum)
        enddo
      enddo

      ! time loop
      do iIter = 1,niter
        ! We should test the use of this value versus the one above
        do iCell = 1,nCells
          ! on u levels
          sigav = 0.5_RKIND*(areaFraction(1,iCell) + areaFraction(2,iCell))
          wtav  = 0.5_RKIND*(wt(i1,1,iCell) + wt(i1,2,iCell))
          w2av  = 0.5_RKIND*(w2(i1,1,iCell) + w2(i1,2,iCell))
          ! TODO: Check the coefficient
          ! w2t(1,iCell) = (1.0_RKIND - 2.0_RKIND*sigav)*wtav*sqrt(w2av) / (EPSILON + &
          !                sigav*(1.0_RKIND-sigav))
          w2t(1,iCell) = (1.0_RKIND - 2.0_RKIND*sigav) * wtav * &
                         sqrt(w2av / (sigav*(1.0_RKIND-sigav)))
        enddo

        !Kernel 1 inlined versions of the base arrays, needed for later to make them collapsible
        do iCell = 1,nCells
          do k = 2,nVertLevels
            ! on w levels
            w3av = 0.5_RKIND*(w3(i1,k-1,iCell) + w3(i1,k,iCell))

            Sw = w3av / max(w2(i1,k,iCell), epsilon)**1.5_RKIND
            lareaFraction = 0.5_RKIND - 0.5_RKIND * Sw / sqrt(4.0_RKIND + Sw**2)
            lareaFraction = min(max(lareaFraction, 0.01_RKIND), 0.99_RKIND)

            areaFraction(k,iCell) = lareaFraction
            wumd(k,iCell) = sqrt(w2(i1,k,iCell) / (lareaFraction*(1.0_RKIND-lareaFraction)))
            Mc(k,iCell) = lareaFraction*(1.0_RKIND-lareaFraction) * wumd(k,iCell)

            tumd(k,iCell) = wt(i1,k,iCell) / (smallNum + Mc(k,iCell))
            sumd(k,iCell) = ws(i1,k,iCell) / (smallNum + Mc(k,iCell))
          enddo
        enddo

        do iCell = 1,nCells
          do k = 2,nVertLevels
            ! on u levels
            sigav  = 0.5_RKIND*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
            tumdav = 0.5_RKIND*(tumd(k,iCell) + tumd(k+1,iCell))
            sumdav = 0.5_RKIND*(sumd(k,iCell) + sumd(k+1,iCell))
            wumdav = 0.5_RKIND*(wumd(k,iCell) + wumd(k+1,iCell))
            w2t(k,iCell) = sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2*tumdav
            w2s(k,iCell) = sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2*sumdav
            !also use this loop to reset the cliptends for the step
            u2cliptend(k,iCell) = 0.0_RKIND
            v2cliptend(k,iCell) = 0.0_RKIND
            w2cliptend(k,iCell) = 0.0_RKIND
          enddo
        enddo

        !This cell loop computes w3Tend and TOMs for later tendencies.  Should be okay to collapse loops here
        do iCell = 1,nCells
          do k = 1,nVertLevels
            ! on u levels
            Eav  = 0.5_RKIND*(Entrainment(k,iCell) + Entrainment(k+1,iCell))
            Dav  = 0.5_RKIND*(Detrainment(k,iCell) + Detrainment(k+1,iCell))
            u2av = 0.5_RKIND*(u2(i1,k,iCell) + u2(i1,k+1,iCell))
            v2av = 0.5_RKIND*(v2(i1,k,iCell) + v2(i1,k+1,iCell))
            w2av = 0.5_RKIND*(w2(i1,k,iCell) + w2(i1,k+1,iCell))

            sigav   = 0.5_RKIND*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
            wumdav  = 0.5_RKIND*(wumd(k,iCell) + wumd(k+1,iCell))
            tumdav  = 0.5_RKIND*(tumd(k,iCell) + tumd(k+1,iCell))
            sumdav  = 0.5_RKIND*(sumd(k,iCell) + sumd(k+1,iCell))
            KspsUav = 0.5_RKIND*(KspsU(i1,k,iCell) + KspsU(i1,k+1,iCell))
            KspsDav = 0.5_RKIND*(KspsD(i1,k,iCell) + KspsD(i1,k+1,iCell))
            Mcav    = 0.5_RKIND*(Mc(k,iCell) + Mc(k+1,iCell))
            lenav   = 0.5_RKIND*(length(k,iCell) + length(k+1,iCell))

            KEsps = sigav*KspsUav + (1.0_RKIND - sigav)*KspsDav
            ! TODO: rename the velocity scale
            ! TODO: Should we include the subplume scale KE in the velocity scale?
            KE = sqrt((u2av+v2av+w2av) + 0.0*KEsps)

            dz = ze(k,iCell) - ze(k+1,iCell)
            !KE = sqrt(areaFraction(k,iCell)*KspsUav + (1.0 - areaFraction(k,iCell))*KspsDav)

            ! last term (SPS) in (40) of LR01a, also using (35) of LR01b
            Swumd(k,iCell) = - 2.0_RKIND/3.0_RKIND * ( &
                              1.0_RKIND/sigav * &
                              (areaFraction(k,  iCell)*KspsU(i1,k,  iCell) - &
                               areaFraction(k+1,iCell)*KspsU(i1,k+1,iCell)) / dz &
                              - 1.0_RKIND/(1.0_RKIND-sigav) * &
                              ((1.0_RKIND-areaFraction(k,  iCell))*KspsD(i1,k,  iCell) - &
                               (1.0_RKIND-areaFraction(k+1,iCell))*KspsD(i1,k+1,iCell)) / dz )

            ! first and second terms in (52) of LR01a, see also (20) of LR01b
            w3tend1(k,iCell) = wumdav**3*(Eav*(3.0_RKIND*sigav - 2.0_RKIND) + Dav*(3.0_RKIND*sigav - 1.0_RKIND))
            !LPV trying new formulation of d/dz(w4) and 3w2 d/dz(w2) -- Not convinced the LR01a equation is correct
            ! w3tend2 + w3tend3 = third + fourth terms in (52) of LR01a
            ! w4k   = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))* &
            !         (1.0_RKIND - 3.0_RKIND*areaFraction(k,iCell) + 3.0_RKIND* &
            !         areaFraction(k,iCell)**2)*wumd(k,iCell)**4
            ! w4kp1 = areaFraction(k+1,iCell)*(1.0_RKIND - areaFraction(k+1,iCell))* &
            !         (1.0_RKIND - 3.0_RKIND*areaFraction(k+1,iCell) + 3.0_RKIND* &
            !         areaFraction(k+1,iCell)**2)*wumd(k+1,iCell)**4
            ! w3tend2(k,iCell) = -(w4k - w4kp1) / dz
            !now the w3term
            ! w2k   = areaFraction(k,iCell)*(1.0_RKIND-areaFraction(k,iCell))*wumd(k,iCell)**2
            ! w2kp1 = areaFraction(k+1,iCell)*(1.0_RKIND-areaFraction(k+1,iCell))*wumd(k+1,iCell)**2
            ! w3tend3(k,iCell) = 1.5_RKIND*(w2k**2 - w2kp1**2) / dz

            ! third term in (52) of LR01a
            w3tend2(k,iCell) = - wumdav**3 * (6.0_RKIND*sigav**2 - 6.0_RKIND*sigav + 1.0_RKIND) * &
                                 ( Mc(k,iCell)-Mc(k+1,iCell) ) / dz
            ! fourth term in (52) of LR01a
            w3tend3(k,iCell) = - 1.5_RKIND * sigav * (1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav) * &
                                 wumdav**2 * &
                                 ( (1.0_RKIND - 2.0_RKIND*areaFraction(k,  iCell))*wumd(k,  iCell)**2 &
                                 - (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell))*wumd(k+1,iCell)**2 &
                                 ) / dz

            ! the pressure term in the sources and sinks in (52) of LR01a
            ! TODO: check the pressure term
            w3tend4(k,iCell) = - C_mom_w3*KE / (1E-15_RKIND+ sqrt(2.0_RKIND)*lenAv)*w3(i1,k,iCell)

            ! the buoyancy term in the sources and sinks in (52) of LR01a
            w3tend5(k,iCell) = 3.0_RKIND * (1.0_RKIND - c11) * (1.0_RKIND - 2.0_RKIND*sigav) * &
                               Mcav * wumdav * gravity * &
                               (alphaT(k,iCell)*tumdav - betaS(k,iCell)*sumdav)

            ! the SPS term in the sources and sinks in (52) of LR01a
            w3tend6(k,iCell) = 3.0_RKIND * (1.0_RKIND - 2.0_RKIND*sigav) * Mcav * wumdav * Swumd(k,iCell)

            ! w3 tendency
            w3tend(k,iCell) = w3tend1(k,iCell) + w3tend2(k,iCell) + w3tend3(k,iCell) + &
                              w3tend4(k,iCell) + w3tend5(k,iCell) + w3tend6(k,iCell)

            if (k>1 .and. k < nVertLevels) then
              ! TODO: why do we need a diffusion term?
              w3tend(k,iCell) = w3tend(k,iCell) + kappa_w3*(w3(i1,k-1,iCell) &
                 - w3(i1,k+1,iCell)) / (zm(k-1,iCell) - zm(k+1,iCell))**2
            endif

            ! now get all the downgradient TOMs
            Ksps = 0.5_RKIND*((areaFraction(k,iCell)*KspsU(i1,k,iCell) + (1.0_RKIND - &
              areaFraction(k,iCell))*KspsD(i1,k,iCell)) + (areaFraction(k+1,iCell)* &
              KspsU(i1,k+1,iCell) + (1.0_RKIND-areaFraction(k+1,iCell))* &
              KspsD(i1,k+1,iCell)))
            KE = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
            KEp1 = 0.5_RKIND*(u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
            lenav = 0.5_RKIND*(length(k,iCell) + length(k+1,iCell))
            diff = C_mom * sqrt(0.5_RKIND*(KE + KEp1)) * lenav
            dz = ze(k,iCell) - ze(k+1,iCell)
            uw2(k,iCell) = -diff*(uw(i1,k,iCell) - uw(i1,k+1,iCell)) / dz
            vw2(k,iCell) = -diff*(vw(i1,k,iCell) - vw(i1,k+1,iCell)) / dz
            u2w(k,iCell) = -diff*(u2(i1,k,iCell) - u2(i1,k+1,iCell)) / dz
            v2w(k,iCell) = -diff*(v2(i1,k,iCell) - v2(i1,k+1,iCell)) / dz
            uvw(k,iCell) = -diff*(uv(i1,k,iCell) - uv(i1,k+1,iCell)) / dz

            diff = C_therm*sqrt(0.5_RKIND*(KE + KEp1)) * lenav
            uwt(k,iCell) = -diff*(ut(i1,k,iCell) - ut(i1,k+1,iCell)) / dz
            vwt(k,iCell) = -diff*(vt(i1,k,iCell) - vt(i1,k+1,iCell)) / dz
            uws(k,iCell) = -diff*(us(i1,k,iCell) - us(i1,k+1,iCell)) / dz
            vws(k,iCell) = -diff*(vs(i1,k,iCell) - vs(i1,k+1,iCell)) / dz

          enddo !nVertLevels
        enddo !nCells

        ! Ashwath -- I'm removing this possibility from this branch since it will be near impossible to
        ! parallelize well.  It's a loop over k with an inner loop expanding underneath.
        if (config_adc_use_old_length_scale) then
		    call dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
        endif

        !next all second order moment tendencies also apply tendencies
        do iCell = 1,nCells
          do k = 2,nVertLevels
            ! on w level
            dzmid = (zm(k-1,iCell) - zm(k,iCell))
            dz = ze(k,iCell) - ze(k+1,iCell)
            B = gravity*(alphaT(k,iCell)*wt(i1,k,iCell) - betaS(k,iCell)*ws(i1,k,iCell))

            Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dzmid
            Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dzmid
            Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dzmid
            Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dzmid

            KEsps = areaFraction(k,iCell)*KspsU(i1,k,iCell) +  &
                      (1.0_RKIND-areaFraction(k,iCell))*KspsD(i1,k,iCell)
            KE = 0.5_RKIND * (u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
            ! TODO: Should we include the subplume scale KE in the velocity scale?
            vscale = sqrt(2.0_RKIND * KE)

            ! See (48) of LR01a for the w2 equation
            ! dissipation by entrainment/detrainment of plumes
            ! (18) of LR01b
            w2tend1(k,iCell) = - wumd(k,iCell)**2 * (Entrainment(k,iCell) + Detrainment(k,iCell))
            ! turbulent transport
            w2tend2(k,iCell) = - ( Mc(k-1,iCell) * &
                               (1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell)) * wumd(k-1,iCell)**2 &
                                 - Mc(k+1,iCell) * &
                               (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell)) * wumd(k+1,iCell)**2 &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell))
            ! pressure term: slow, return to isotropy
            ! the missing w2 term is added when time-stepping the w2 equation
            ! TODO: is it consisten with the u2 and v2 equations?
            w2tend3(k,iCell) = vscale / (1.0E-15_RKIND + c_pv*2.0_RKIND*sqrt(2.0_RKIND)* &
              length(k,iCell))*(u2(i1,k,iCell) + v2(i1,k,iCell))/3.0_RKIND
            ! buoyancy term + the buoyancy contribution in the pressure term
            w2tend4(k,iCell) = (2.0_RKIND - 4.0_RKIND/3.0_RKIND*C_2) * &
                               Mc(k,iCell) * gravity * &
                               (alphaT(k,iCell)*tumd(k,iCell) - betaS(k,iCell)*sumd(k,iCell))
            ! pressure term: rapid + subplume scale term
            ! TODO: use MF formulation for uw and vw as well?
            w2tend5(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 - alpha2) * &
                               (uw(i1,k,iCell)*Uz + vw(i1,k,iCell)*Vz) + &
                               Mc(k,iCell) * (Swumd(k-1,iCell) + Swumd(k,iCell))

            ! w2 tendency
            w2tend(k,iCell) = w2tend1(k,iCell) + w2tend2(k,iCell) + w2tend3(k,iCell) + &
                              w2tend4(k,iCell) + w2tend5(k,iCell)

            ! inverse of the time scale
            ! TODO: check the coefficients
            invTaupt(k,iCell) = vscale / (2.0_RKIND*sqrt(2.0_RKIND) * c_pt * length(k,iCell))
            invTaups(k,iCell) = vscale / (2.0_RKIND*sqrt(2.0_RKIND) * c_ps * length(k,iCell))
            invTaupv(k,iCell) = vscale / (2.0_RKIND*sqrt(2.0_RKIND) * c_pv * length(k,iCell))

            ! See (49) of LR01a for the wt equation
            ! sources/sinks due to entrainment/detrainment + transport term
            wttend1(k,iCell) = - (Entrainment(k,iCell) + Detrainment(k,iCell)) * &
                                 wumd(k,iCell) * tumd(k,iCell) &
                               - ( (1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell)) * &
                                   wumd(k-1,iCell) * tumd(k-1,iCell) * Mc(k-1,iCell) &
                                 - (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell)) * &
                                   wumd(k+1,iCell) * tumd(k+1,iCell) * Mc(k+1,iCell) &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell))
            ! gradient production term
            wttend2(k,iCell) = - Mc(k,iCell) * wumd(k,iCell) * Tz
            ! buoyancy term
            wttend3(k,iCell) = (1.0_RKIND - gamma1) * gravity * &
                               areaFraction(k,iCell) * (1.0_RKIND - areaFraction(k,iCell)) * &
                               (alphaT(k,iCell)*tumd(k,iCell)**2 - &
                                betaS(k,iCell)*tumd(k,iCell)*sumd(k,iCell))
            ! pressure term: rapid term
            wttend4(k,iCell) = - alpha3/4.0_RKIND * (ut(i1,k,iCell)*Uz + vt(i1,k,iCell)*Vz)
            ! pressure term:
            ! the missing wt term is added when time-stepping wt equation
            ! TODO: why do we need a diffusion term?
            wttend5(k,iCell) = kappa_FL * (wt(i1,k-1,iCell) - wt(i1,k+1,iCell)) / &
                               (ze(k-1,iCell) - ze(k+1,iCell))**2
            ! sources/sinks due to subplume terms
            wttend6(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell)) * &
                               tumd(k,iCell) * 0.5_RKIND * (Swumd(k-1,iCell) + Swumd(k,iCell)) &
                             - Mc(k,iCell) * ( 1.0_RKIND / areaFraction(k,iCell) * &
                               ( areaFraction(k-1,iCell) * wt_spsU(k-1,iCell) &
                               - areaFraction(k+1,iCell) * wt_spsU(k+1,iCell) &
                               ) / (ze(k-1,iCell) - ze(k+1,iCell)) &
                                             - 1.0_RKIND / (1.0_RKIND - areaFraction(k,iCell)) * &
                               ( (1.0_RKIND - areaFraction(k-1,iCell)) * wt_spsD(k-1,iCell) &
                               - (1.0_RKIND - areaFraction(k+1,iCell)) * wt_spsD(k+1,iCell) &
                               ) / (ze(k-1,iCell) - ze(k+1,iCell)) &
                                             )
            ! wt tendency
            wttend(k,iCell) = wttend1(k,iCell) + wttend2(k,iCell) + wttend3(k,iCell) + &
                              wttend4(k,iCell) + wttend5(k,iCell) + wttend6(k,iCell)

            ! See (49) of LR01a for the wt equation
            ! sources/sinks due to entrainment/detrainment + transport term
            wstend1(k,iCell) = - (Entrainment(k,iCell) + Detrainment(k,iCell)) * &
                                 wumd(k,iCell)*sumd(k,iCell) &
                               - ( (1.0_RKIND - 2.0_RKIND*areaFraction(k-1,iCell)) * &
                                   wumd(k-1,iCell) * sumd(k-1,iCell) * Mc(k-1,iCell) &
                                 - (1.0_RKIND - 2.0_RKIND*areaFraction(k+1,iCell)) * &
                                   wumd(k+1,iCell) * sumd(k+1,iCell) * Mc(k+1,iCell) &
                                 ) / (ze(k-1,iCell) - ze(k+1,iCell))
            ! gradient production term
            wstend2(k,iCell) = - Mc(k,iCell) * wumd(k,iCell) * Sz
            ! buoyancy term
            wstend3(k,iCell) = (1.0_RKIND - gamma1) * gravity * &
                               areaFraction(k,iCell) * (1.0_RKIND - areaFraction(k,iCell)) * &
                               (alphaT(k,iCell)*tumd(k,iCell)*sumd(k,iCell) - &
                                betaS(k,iCell)*sumd(k,iCell)**2)
            ! pressure term: rapid term
            wstend4(k,iCell) = - alpha3/4.0_RKIND * (us(i1,k,iCell)*Uz + vs(i1,k,iCell)*Vz)
            ! pressure term:
            ! the missing ws term is added when time-stepping ws equation
            ! TODO: why do we need a diffusion term?
            wsTend5(k,iCell) = kappa_FL * (ws(i1,k-1,iCell) - ws(i1,k+1,iCell)) / &
                               (ze(k-1,iCell) - ze(k+1,iCell))**2
            ! sources/sinks due to subplume terms
            wstend6(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell)) * &
                               sumd(k,iCell) * 0.5_RKIND * (Swumd(k-1,iCell) + Swumd(k,iCell)) &
                             - Mc(k,iCell) * ( 1.0_RKIND / areaFraction(k,iCell) * &
                               ( areaFraction(k-1,iCell) * ws_spsU(k-1,iCell) &
                               - areaFraction(k+1,iCell) * ws_spsU(k+1,iCell) &
                               ) / (ze(k-1,iCell) - ze(k+1,iCell)) &
                                             - 1.0_RKIND / (1.0_RKIND - areaFraction(k,iCell)) * &
                               ( (1.0_RKIND - areaFraction(k-1,iCell)) * ws_spsD(k-1,iCell) &
                               - (1.0_RKIND - areaFraction(k+1,iCell)) * ws_spsD(k+1,iCell) &
                               ) / (ze(k-1,iCell) - ze(k+1,iCell)) &
                                             )
            ! ws tendency
            wstend(k,iCell) = wstend1(k,iCell) + wstend2(k,iCell) + wstend3(k,iCell) + &
                              wstend4(k,iCell) + wstend5(k,iCell) + wstend6(k,iCell)

            ! uw equation
            ! transport term
            uwtend1(k,iCell) = - (uw2(k-1,iCell) - uw2(k,iCell)) / dzmid
            ! shear production + rapid pressure term: dU/dz
            uwtend2(k,iCell) = 0.5_RKIND * ( (alpha0 - 4.0_RKIND/3.0_RKIND*alpha1) * KE + &
                                             (alpha1 - alpha2) * u2(i1,k,iCell) + &
                                             (alpha1 + alpha2 - 2.0_RKIND) * w2(i1,k,iCell) &
                                           ) * Uz
            ! rapid pressure term: dV/dz
            uwtend3(k,iCell) = 0.5_RKIND * (alpha1 - alpha2) * uv(i1,k,iCell) * Vz
            ! buoyancy term
            ! TODO: use (1-cb) to replace beta5 here
            uwtend4(k,iCell) = beta5 * gravity * &
                               (alphaT(k,iCell)*ut(i1,k,iCell) - betaS(k,iCell)*us(i1,k,iCell))
            ! return to isotropy + diffusion
            ! TODO: move return to isotropy to time stepping?
            uwtend5(k,iCell) = - invTaupv(k,iCell) * uw(i1,k,iCell) + &
                               kappa_FL * (uw(i1,k-1,iCell) - uw(i1,k+1,iCell)) / &
                               (ze(k-1,iCell) - ze(k+1,iCell))**2
            ! uw tendency
            uwtend(k,iCell) = uwtend1(k,iCell) + uwtend2(k,iCell) + uwtend3(k,iCell) + &
                              uwtend4(k,iCell) + uwtend5(k,iCell)

            ! vw equation
            ! transport term
            vwtend1(k,iCell) = - (vw2(k-1,iCell) - vw2(k,iCell)) / dzmid
            ! rapid pressure term: dU/dz
            vwtend2(k,iCell) = 0.5_RKIND * (alpha1 - alpha2) * uv(i1,k,iCell) * Uz
            ! shear production + rapid pressure term: dV/dz
            vwtend3(k,iCell) = 0.5_RKIND * ( (alpha0 - 4.0_RKIND/3.0_RKIND*alpha1) * KE + &
                                             (alpha1 - alpha2) * v2(i1,k,iCell) + &
                                             (alpha1 + alpha2 - 2.0_RKIND) * w2(i1,k,iCell) &
                                           ) * Vz
            ! buoyancy term
            ! TODO: use (1-cb) to replace beta5 here
            vwtend4(k,iCell) = beta5 * gravity * &
                               (alphaT(k,iCell)*vt(i1,k,iCell) - betaS(k,iCell)*vs(i1,k,iCell))
            ! return to isotropy + diffusion
            ! TODO: move return to isotropy to time stepping?
            vwtend5(k,iCell) = - invTaupv(k,iCell) * vw(i1,k,iCell) + &
                                 kappa_FL * (vw(i1,k-1,iCell) - vw(i1,k+1,iCell)) / &
                                 (ze(k-1,iCell) - ze(k+1,iCell))**2
            ! vw tendency
            vwtend(k,iCell) = vwtend1(k,iCell) + vwtend2(k,iCell) + vwtend3(k,iCell) + &
                              vwtend4(k,iCell) + vwtend5(k,iCell)

            ! uv tendency
            uvtend(k,iCell) = &
            ! transport
                             - (uvw(k-1,iCell) - uvw(k,iCell)) / dz &
            ! shear production + rapid pressure term
                             - (1.0_RKIND - 0.5_RKIND*(alpha1+alpha2)) * &
                               (uw(i1,k,iCell)*Vz + vw(i1,k,iCell)*Uz) &
            ! return to isotropy
            ! TODO: move return to isotropy to time stepping?
                             - invTaupv(k,iCell) * uv(i1,k,iCell) &
            ! diffusion
                             + kappa_VAR * (uv(i1,k-1,iCell) - uv(i1,k+1,iCell)) / &
                               (ze(k-1,iCell) - ze(k+1,iCell))**2

            ! u2 equation
            ! transport
            u2tend1(k,iCell) = - (u2w(k-1,iCell) - u2w(k,iCell)) / dzmid
            ! shear production + rapid pressure term: dU/dz
            u2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 + alpha2 - 2.0_RKIND) * uw(i1,k,iCell) * Uz
            ! rapid pressure term: dV/dz
            u2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1 * vw(i1,k,iCell) * Vz
            ! buoyancy term
            ! TODO: use (1-cb) to replace beta5 here
            u2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND * (1.0_RKIND-beta5) * B
            ! return to isotropy
            ! TODO: note that the term proportional to u2 is in time stepping
            ! TODO: what is the eps term?
!            u2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
!              invTaupv(k,iCell)*(2.0_RKIND*KE/3.0_RKIND - u2(i1,k,iCell))
            u2tend5(k,iCell) = - 2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
                               invTaupv(k,iCell) * (v2(i1,k,iCell) + w2(i1,k,iCell)) / 3.0_RKIND
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
            ! TODO: use (1-cb) to replace beta5 here
            v2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND * (1.0_RKIND-beta5) * B
            ! return to isotropy
            ! TODO: note that the term proportional to v2 is in time stepping
            ! TODO: what is the eps term?
!            v2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + invTaupv(k,iCell)* &
!              (2.0_RKIND*KE/3. - v2(i1,k,iCell))
            v2tend5(k,iCell) = - 2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
                               invTaupv(k,iCell) * (u2(i1,k,iCell) + w2(i1,k,iCell)) / 3.0_RKIND
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
                             - (1.0_RKIND - alpha3) * wt(i1,k,iCell) * Uz &
            ! return to isotropy
            ! TODO: move return to isotropy to time stepping?
                             - invTaupt(k,iCell) * ut(i1,k,iCell)

            ! vt tendency
            vttend(k,iCell) = &
            ! transport
                             - (vwt(k-1,iCell) - vwt(k,iCell)) / dz &
            ! gradient term
                             - vw(i1,k,iCell) * Tz &
            ! rapid pressure term
                             - (1.0_RKIND - alpha3) * wt(i1,k,iCell) * Vz &
            ! return to isotropy
            ! TODO: move return to isotropy to time stepping?
                             - invTaupt(k,iCell) * vt(i1,k,iCell)

            ! us tendency
            ustend(k,iCell) = &
            ! transport
                             - (uws(k-1,iCell) - uws(k,iCell)) / dz &
            ! gradient term
                             - uw(i1,k,iCell) * Sz &
            ! rapid pressure term
                             - (1.0_RKIND - alpha3) * ws(i1,k,iCell) * Uz &
            ! return to isotropy
            ! TODO: move return to isotropy to time stepping?
                             - invTaups(k,iCell) * us(i1,k,iCell)

            ! vs tendency
            vstend(k,iCell) = &
            ! transport
                             - (vws(k-1,iCell) - vws(k,iCell)) / dz &
            ! gradient term
                             - vw(i1,k,iCell) * Sz &
            ! rapid pressure term
                             - (1.0_RKIND - alpha3) * ws(i1,k,iCell) * Vz &
            ! return to isotropy
            ! TODO: move return to isotropy to time stepping?
                             - invTaups(k,iCell) * vs(i1,k,iCell)

            ! TODO: not used? <
            t2(i2,k,iCell) = tumd(k,iCell)**2*areaFraction(k,iCell)* &
              (1.0_RKIND-areaFraction(k,iCell))
            s2(i2,k,iCell) = sumd(k,iCell)**2*areaFraction(k,iCell)* &
              (1.0_RKIND-areaFraction(k,iCell))
            ts(i2,k,iCell) = tumd(k,iCell)*sumd(k,iCell)* &
              areaFraction(k,iCell)*(1.0_RKIND-areaFraction(k,iCell))
            ! >

            if (.not. config_adc_use_old_length_scale) then
              !epsilon tendency is next
              KEm1  = 0.5_RKIND * (u2(i1,k-1,iCell) + v2(i1,k-1,iCell) + w2(i1,k-1,iCell))
              KE    = 0.5_RKIND * (u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
              KEp1  = 0.5_RKIND * (u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
              tauM1 = 2.0_RKIND * KEm1 / (eps(i1,k-1,iCell) + smallNum )
              tau   = 2.0_RKIND * KE / (eps(i1,k,iCell) + smallNum)
              tauP1 = 2.0_RKIND * KEp1 / (eps(i1,k+1,iCell) + smallNum)
              tauAV = 0.5_RKIND * (tauUP + tauDN)

              tomUP = 0.5577_RKIND/1.3_RKIND * (0.5_RKIND*(KEm1+KE))**2 / &
                      (0.5_RKIND*(eps(i1,k-1,iCell) + eps(i1,k,iCell)) + smallNum) * &
                      (eps(i1,k-1,iCell) - eps(i1,k,iCell)) / (ze(k-1,iCell) - ze(k,iCell))
              tomDN = 0.5577_RKIND/1.3_RKIND * (0.5_RKIND*(KE+KEp1))**2 / &
                      (0.5_RKIND*(eps(i1,k,iCell) + eps(i1,k+1,iCell)) + smallNum) * &
                      (eps(i1,k,iCell) - eps(i1,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell))

              epstend(k,iCell) = (tomUP - tomDN) / (zm(k-1,iCell) - zm(k,iCell)) - &
              !-0.8_RKIND*((tomUP / (0.5_RKIND*(tauM1 + tau) + 1.0E-10_RKIND) - tomDN / &
               ! (0.5_RKIND*(tau + tauP1) + 1.0E-10_RKIND)) / (zm(k-1,iCell) - zm(k,iCell))) - &
                                 2.88_RKIND/(tau + smallNum) * uw(i1,k,iCell) * &
                                 (uvel(k-1,iCell) - uvel(k,iCell)) / (zm(k-1,iCell) - zm(k,iCell)) - &
                                 2.88_RKIND/(tau + smallNum) * vw(i1,k,iCell) * &
                                 (vvel(k-1,iCell) - vvel(k,iCell)) / (zm(k-1,iCell) - zm(k,iCell)) + &
                                 2.88_RKIND/(tau + smallNum) * B - &
                                 3.84_RKIND * eps(i1,k,iCell)/(tau + smallNum)
            endif

            ! finally update all subplume fluxes
            ! (32) and (33) of LR01b
            if (BVF(k,iCell) <= 0.0_RKIND) then
              lenspsU(k,iCell) = dzmid
              lenspsD(k,iCell) = dzmid
            else
              lenspsU(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsU(i1,k,iCell)/BVF(k,iCell)))
              lenspsD(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsD(i1,k,iCell)/BVF(k,iCell)))
            endif

            ! (34) of LR01b
            KmU(k,iCell) = 0.1_RKIND*lenspsU(k,iCell)*sqrt( KspsU(i1,k,iCell) )
            KhU(k,iCell) = (1.0_RKIND+2.0_RKIND*lenspsU(k,iCell)/dz)*KmU(k,iCell)
            wt_spsU(k,iCell) =  -KhU(k,iCell)*Tz
            ws_spsU(k,iCell) =  -KhU(k,iCell)*Sz

            KmD(k,iCell) = 0.1_RKIND*lenspsD(k,iCell)*sqrt( KspsD(i1,k,iCell) )
            KhD(k,iCell) = (1.0_RKIND+2.0_RKIND*lenspsD(k,iCell)/dz)*KmD(k,iCell)
            wt_spsD(k,iCell) = -KhD(k,iCell)*Tz
            ws_spsD(k,iCell) = -KhD(k,iCell)*Sz

            !change length scale to lenup and lendown
            ! (11) of LR01b
            Entrainment(k,iCell) = Cww_E * areaFraction(k,iCell) * (1.0_RKIND-areaFraction(k,iCell)) * &
                                   Mc(k,iCell) / ( lendn(k,iCell) + smallNum )
            Detrainment(k,iCell) = Cww_D * areaFraction(k,iCell) * (1.0_RKIND-areaFraction(k,iCell)) * &
                                   Mc(k,iCell) / ( lenup(k,iCell) + smallNum )

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
                     - gravity * (alphaT(k,iCell)*wt_spsU(k,iCell) - betaS(k,iCell)*ws_spsU(k,iCell)) &
            ! transport of SPS TKE and pressure term: (29) of LR01b
                     + ((KmU(k-1,iCell) + KmU(k,iCell)) * (KspsU(i1,k-1,iCell) - KspsU(i1,k,iCell)) / &
                        (ze(k-1,iCell) - ze(k,iCell)) &
                      - (KmU(k,iCell) + KmU(k+1,iCell)) * (KspsU(i1,k,iCell) - KspsU(i1,k+1,iCell)) / &
                        (ze(k,iCell) - ze(k+1,iCell)) &
                       ) / dzmid &
            ! dissipation: (30) of LR01b
                     - Cval*KspsU(i1,k,iCell)**1.5_RKIND / lenspsU(k,iCell) &
            ! TODO: what is this term?
                     + eps(i1,k,iCell) / (2.0_RKIND*areaFraction(k,iCell)) &
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
                     - gravity * (alphaT(k,iCell)*wt_spsD(k,iCell) - betaS(k,iCell)*ws_spsD(k,iCell)) &
            ! transport of SPS TKE and pressure term: (29) of LR01b
                     + ((KmD(k-1,iCell) + KmD(k,iCell)) * (KspsD(i1,k-1,iCell) - KspsD(i1,k,iCell)) / &
                        (ze(k-1,iCell) - ze(k,iCell)) &
                      - (KmD(k,iCell) + KmD(k+1,iCell)) * (KspsD(i1,k,iCell) - KspsD(i1,k+1,iCell)) / &
                        (ze(k,iCell) - ze(k+1,iCell)) &
                       ) / dzmid &
            ! dissipation: (30) of LR01b
                     - Cval*KspsD(i1,k,iCell)**1.5_RKIND / lenspsD(k,iCell) &
            ! TODO: what is this term?
                     + eps(i1,k,iCell) / (2.0_RKIND*(1.0_RKIND - areaFraction(k,iCell))) &
            ! shear production
                     + KmD(k,iCell) * (Uz**2 + Vz**2)

          enddo ! nVertLevels
        enddo ! nCells

        !In this step we update second moments except w3 which needs updated w2
        do iCell = 1,nCells
          do k = 2,nVertLevels
            !update second order moment tendency here
            w2(i2,k,iCell) = (w2(i1,k,iCell) + dt_small*w2tend(k,iCell)) / &
                             (1.0_RKIND + dt_small*invTaupv(k,iCell)*2.0_RKIND/3.0_RKIND)
            if (w2(i2,k,iCell) < epsilon) then
              w2cliptend(k,iCell) = epsilon-w2(i2,k,iCell)
              w2(i2,k,iCell) = epsilon
            endif

            if (abs(w2(i2,k,iCell)) > 1.0_RKIND) then
              call mpas_log_write("ERROR: w2 out of range, w2 = $r, location = $i, $i", &
                                  MPAS_LOG_CRIT,realArgs=(/w2(i2,k,iCell)/),intArgs=(/k,iCell/))
            endif

            u2(i2,k,iCell) = (u2(i1,k,iCell) + dt_small*u2tend(k,iCell)) / &
                             (1.0_RKIND + dt_small*invTaupv(k,iCell)*2.0_RKIND/3.0_RKIND)
            if (u2(i2,k,iCell) < 0.0_RKIND) then
              u2cliptend(k,iCell) = -u2(i2,k,iCell)
              u2(i2,k,iCell) = 0.0_RKIND
            endif

            v2(i2,k,iCell) = (v2(i1,k,iCell) + dt_small*v2tend(k,iCell)) / &
                             (1.0_RKIND + dt_small*invTaupv(k,iCell)*2.0_RKIND/3.0_RKIND)
            if (v2(i2,k,iCell) < 0.0_RKIND) then
              v2cliptend(k,iCell) = -v2(i2,k,iCell)
              v2(i2,k,iCell) = 0.0_RKIND
            endif

            uw(i2,k,iCell) = uw(i1,k,iCell) + dt_small*uwtend(k,iCell)
            vw(i2,k,iCell) = vw(i1,k,iCell) + dt_small*vwtend(k,iCell)
            uv(i2,k,iCell) = uv(i1,k,iCell) + dt_small*uvtend(k,iCell)
            ut(i2,k,iCell) = ut(i1,k,iCell) + dt_small*uttend(k,iCell)
            vt(i2,k,iCell) = vt(i1,k,iCell) + dt_small*vttend(k,iCell)
            wt(i2,k,iCell) = (wt(i1,k,iCell) + dt_small*wttend(k,iCell)) / &
                             (1.0_RKIND + dt_small*invTaupt(k,iCell))
            us(i2,k,iCell) = us(i1,k,iCell) + dt_small*ustend(k,iCell)
            vs(i2,k,iCell) = vs(i1,k,iCell) + dt_small*vstend(k,iCell)
            ws(i2,k,iCell) = (ws(i1,k,iCell) + dt_small*wstend(k,iCell)) / &
                             (1.0_RKIND + dt_small*invTaups(k,iCell))
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

            if (.not. config_adc_use_old_length_scale)then
              eps(i2,k,iCell) = eps(i1,k,iCell) + dt_small*epstend(k,iCell)
              if (eps(i2,k,iCell) < 1.0E-10_RKIND) then
                eps(i2,k,iCell) = 1.0E-10_RKIND
              endif
              !Amrapalli: 16.6 here can be tuned if needed, would show up in u2 and v2 primarily
              ! also in KEspsU/D
              KE = 0.5_RKIND * (u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))
              length(k,iCell) = KE**1.5_RKIND / 16.6_RKIND * eps(i2,k,iCell)
            else
              KE = 0.5_RKIND * (u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))
              eps(:,k,iCell) = KE**1.5_RKIND / 16.6_RKIND * length(k,iCell)
            endif
            KspsU(i2,k,iCell) = max(epsilon, KspsU(i1,k,iCell) + dt_small*KspsUtend(k,iCell))
            KspsD(i2,k,iCell) = max(epsilon, KspsD(i1,k,iCell) + dt_small*KspsDtend(k,iCell))

          enddo !nVertLevels
        enddo !nCells for second order moment tendencies

        ! update the third order w3 and mean fields
        do iCell = 1,nCells
          do k = 1,nVertLevels
            w3check = (w2(i2,k,iCell) + w2(i2,k+1,iCell))**1.5_RKIND
            w3(i2,k,iCell) = min(w3(i1,k,iCell) + dt_small*w3tend(k,iCell), w3check)

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

            wtSumUp = wt(i1,k,iCell) + areaFraction(k,iCell)*wt_spsU(k,iCell) + &
                      (1.0_RKIND - areaFraction(k,iCell))*wt_spsD(k,iCell)
            wtSumDn = wt(i1,k+1,iCell) + areaFraction(k+1,iCell)*wt_spsU(k+1,iCell) + &
                      (1.0_RKIND - areaFraction(k+1,iCell))*wt_spsD(k+1,iCell)
            wsSumUp = ws(i1,k,iCell) + areaFraction(k,iCell)*ws_spsU(k,iCell) + &
                      (1.0_RKIND - areaFraction(k,iCell))*ws_spsD(k,iCell)
            wsSumDn = ws(i1,k+1,iCell) + areaFraction(k+1,iCell)*ws_spsU(k+1,iCell) + &
                      (1.0_RKIND - areaFraction(k+1,iCell))*ws_spsD(k+1,iCell)

            activeTracers(1,k,iCell) = activeTracers(1,k,iCell) - dt_small*(wtSumUp - wtSumDn) / &
                                       (ze(k,iCell) - ze(k+1,iCell))
            activeTracers(2,k,iCell) = activeTracers(2,k,iCell) - dt_small*(wsSumUp - wsSumDn) / &
                                       (ze(k,iCell) - ze(k+1,iCell))
          enddo
        enddo
      enddo

  end subroutine compute_ADC_tends

end module ocn_adc_mixing_fused
