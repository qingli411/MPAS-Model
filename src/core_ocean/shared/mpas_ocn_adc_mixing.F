module ocn_adc_mixing

  use mpas_constants
  use mpas_kind_types
  use mpas_log

  use ocn_constants
  use ocn_turbulence
  use ocn_config

  implicit none

  logical :: defineFirst, stopflag

  integer :: i1,i2

!  integer,parameter :: ntimes = 2

  contains

  subroutine swap_time_levels

    if(i1 == 1) then
      i1 = 2
    else
      i1 = 1
    endif

    if(i2 == 1) then
      i2 = 2
    else
      i2 = 1
    endif

  end subroutine swap_time_levels

  subroutine build_diagnostic_arrays(nCells,nVertLevels,BVF,wtsfc,wssfc, &
        uwsfc, vwsfc, alphaT,betaS, boundaryLayerDepth)
  !construct dTdz, dSdz, dbdz
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nCells),intent(in) :: wtsfc, wssfc, uwsfc, vwsfc
    real,dimension(nCells),intent(inout) :: boundaryLayerDepth
    real,dimension(nVertLevels,nCells),intent(in) :: BVF, alphaT, betaS
    integer :: iCell, k, idx
    real,dimension(nCells) :: wstar
    logical :: first
    real :: maximum, Q, KE

    first = .true.

    do iCell=1,nCells
      maximum = -1.0e12_RKIND
      idx =2
      first = .true.
      do k=2,nVertLevels-1

        if(BVF(k,iCell) > 1.005*maximum .and. first) then
          maximum = BVF(k,iCell)
          idx = k
        elseif(BVF(k,iCell) < maximum) then
          first = .false.
        endif
      enddo

      boundaryLayerDepth(iCell) = abs(ze(idx,iCell))
      Q = gravity*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))* &
        boundaryLayerDepth(iCell)
      if(Q > 0) then
        wstar(iCell) = abs(Q)**(1.0/3.0)
      else
        wstar(iCell) = 0.0
      endif

      u2(:,1,iCell) = 4.0*uwsfc(iCell) + 0.3*wstar(iCell)**2.0
      v2(:,1,iCell) =4.*uwsfc(iCell) + 0.3*wstar(iCell)**2.0
      uw(:,1,iCell) = -uwsfc(iCell)
      vw(:,1,iCell) = vwsfc(iCell)
      wt(:,1,iCell) = wtsfc(iCell)
      ws(:,1,iCell) = wssfc(iCell)
      KE = 0.5*(u2(i1,1,iCell) + v2(i1,1,iCell))
      eps(:,1,iCell) = KE**1.5/(0.5*(ze(1,iCell) - ze(2,iCell) + 1.0E-10_RKIND))
    enddo

  end subroutine build_diagnostic_arrays

  subroutine dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(2,nVertLevels,nCells),intent(in) :: activeTracers
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT, betaS
    integer :: i,k, ij

    real,dimension(nVertLevels) :: B, Bup, Bdo
    real,dimension(nVertLevels+1) :: tke, BupEdge, BdoEdge
    real :: sav, tudav, sudav, Tup, Tdo, Sup, Sdo
    real :: s1, z1, zV, sumv, minlen

    !NOTE: will need to convert to some form of displaced density in the mpas framework soon
    !possibly go back to a more traditional length scale formulation
    do i=1,nCells
       tke(:) = 0.5*(u2(i2,:,i) + v2(i2,:,i) + w2(i2,:,i))
       do k=1,nVertLevels
          B(k) = -9.806*(-alphaT(k,i)*(activeTracers(1,k,i) - 15.0_RKIND) + betaS(k,i)*   &
                    (activeTracers(2,k,i) - 35.0_RKIND))

          sav = 0.5_RKIND*(areaFraction(k,i) + areaFraction(k+1,i))
          tudav = 0.5_RKIND*(tumd(k,i) + tumd(k+1,i))
          sudav = 0.5*(sumd(k,i) + sumd(k+1,i))

          Tup = activeTracers(1,k,i) + (1.0_RKIND - sav)*tudav
          Tdo = activeTracers(1,k,i) - sav*tudav
          Sup = activeTracers(2,k,i) + (1.0_RKIND - sav)*sudav
          Sdo = activeTracers(2,k,i) - sav*sudav

          Bup(k) = -9.806*(-alphaT(k,i)*(Tup - 15.0) + betaS(k,i)*(Sup - 35.0))
          Bdo(k) = -9.806*(-alphaT(k,i)*(Tdo - 15.0) + betaS(k,i)*(Sdo - 35.0))

          if(k>1) THEN
             BupEdge(k) = 0.5*(Bup(k-1) + Bup(k))
             BdoEdge(k) = 0.5*(Bdo(k-1) + Bdo(k))
          endif
       enddo

       BdoEdge(nVertLevels+1) = BdoEdge(nVertLevels)
       BupEdge(nVertLevels+1) = BupEdge(nVertLevels)

       BdoEdge(1) = BdoEdge(2)
       BupEdge(1) = BupEdge(2)

       do k=2,nVertLevels

          sumv = 0
          ij=k
          lenup(k,i) = 0.0_RKIND
          do while(sumv <= tke(k) .and. ij < nVertLevels+1)
             sumv = sumv + (BupEdge(k) - Bup(ij))*(ze(ij-1,i)-ze(ij,i))
             lenup(k,i) =  lenup(k,i) + abs(ze(ij-1,i)-ze(ij,i))
             ij = ij + 1

             if(sumv > tke(k)) THEN
                ij = ij - 1
                s1 = sumv
                z1 = ze(ij,i)
                zV = ze(ij-1,i)

                sumv = sumv - (BupEdge(k) - Bup(ij))*(ze(ij-1,i)-ze(ij,i))
                lenup(k,i) = lenup(k,i) - abs(ze(ij,i)-ze(ij+1,i))
                lenup(k,i) = max(0.55,lenup(k,i) + abs((z1-zV)/(s1 - sumv)*(tke(k)-sumv)))
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
        sumv = 0
        ij=k
        lendn(k,i) = 0.0_RKIND
        do while(sumv <= tke(k) .and. ij>1)
           sumv = sumv - (BdoEdge(k) - Bdo(ij-1))*(ze(ij-1,i)-ze(ij,i))
           lendn(k,i) = lendn(k,i) + abs(ze(ij-1,i)-ze(ij,i))
           ij = ij - 1

           if(sumv > tke(k)) THEN
              ij = ij + 1
              s1 = sumv
              z1 = ze(ij,i)
              zV = ze(ij-1,i)
              sumv = sumv + (BdoEdge(k) - Bdo(ij-1))*(ze(ij-1,i)-ze(ij,i))
              lendn(k,i) = lendn(k,i) - abs(ze(ij-1,i)-ze(ij,i))
              lendn(k,i) = max(0.55,lendn(k,i) + ((zV-z1)/(s1 - sumv)*(tke(k) - sumv)))

            !    if(Bdo(k-1) - Bdo(k) < 0) then
            !            minlen = abs(zm(k-1,i) - zm(k,i))
            !    else
            !            minlen = 0.55_RKIND
            !    endif

            !   lendn(k,i) = max(minlen,lendn(k,i) + sqrt(2.0_RKIND/(1.0E-15 + abs(BdoEdge(k) -  &
            !                    Bdo(ij-1))*(tke(k) - sumv))))
              exit
           endif
        enddo

        length(k,i) = (2.0_RKIND*lenup(k,i)*lendn(k,i))  &
          / (lenup(k,i) + lendn(k,i))
      enddo
   enddo

   length(1,:) = 0.55_RKIND
   length(nVertLevels+1,:) = 0.55_RKIND

  end subroutine dissipation_lengths2

  subroutine build_length_from_dissipation_eq(nCells,nVertLevels,BVF,boundaryLayerDepth,alphaT,betaS, &
      uvel, vvel, dt)
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT, betaS, uvel, vvel
    real,dimension(nVertLevels,nCells),intent(inout) :: BVF
    real,dimension(nCells) :: boundaryLayerDepth
    real,intent(in) :: dt

    real :: invLen, l, len1, len2, lenmax, KE, integrandTop, integrandBot
    real :: len0, len2_1, len2_2, sfcBuoy, lengthT, bvfT
    real :: KEm1, KEp1, tauUP, tauDN, tomUP, tomDN
    real :: tauM1, tau, tauP1, tauAV
    real,dimension(nVertLevels) :: epstend

    integer :: k,iCell

    do iCell=1,nCells
      epstend(1) = 0.0_RKIND
      do k=2,nVertLevels
        KEm1 = 0.5_RKIND*(u2(i1,k-1,iCell) + v2(i1,k-1,iCell) + w2(i1,k-1,iCell))
        KE = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
        KEp1 = 0.5_RKIND*(u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
        tauM1 = 2.0_RKIND*KEm1 / (eps(i1,k-1,iCell) + 1.0E-15_RKIND)
        tau = 2.0_RKIND*KE / (eps(i1,k,iCell) + 1.0E-15_RKIND)
        tauP1 = 2.0_RKIND*KEp1 / (eps(i1,k+1,iCell) + 1.0E-15_RKIND)
        tomUP = u2w(k-1,iCell) + v2w(k-1,iCell) + 0.0_RKIND*w3(i1,k-1,iCell)
        tomDN = u2w(k,iCell) + v2w(k,iCell) + 0.0_RKIND*w3(i1,k,iCell)
        tauAV = 0.5_RKIND*(tauUP + tauDN)

        tomUP = 0.5577/(0.5_RKIND*(tau + tauM1)+1.0E-10_RKIND)*(eps(i1,k-1,iCell) - &
          eps(i1,k,iCell)) / (ze(k-1,iCell) - ze(k,iCell) + 1.0E-10_RKIND)
        tomDN = 0.5577/(0.5_RKIND*(tau + tauP1)+1.0E-10_RKIND)*(eps(i1,k,iCell) - &
          eps(i1,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell) + 1.0E-10_RKIND)

        epstend(k) = 2.0_RKIND*(tomUP - tomDN) / (zm(k-1,iCell) - zm(k,iCell) + 1.0E-10_RKIND) - &
        !-0.8_RKIND*((tomUP / (0.5_RKIND*(tauM1 + tau) + 1.0E-10_RKIND) - tomDN / &
         ! (0.5_RKIND*(tau + tauP1) + 1.0E-10_RKIND)) / (zm(k-1,iCell) - zm(k,iCell))) - &
          2.88/(tau + 1.0E-10_RKIND)*uw(i1,k,iCell)*(uvel(k-1,iCell) - uvel(k,iCell)) / &
          (zm(k-1,iCell) - zm(k,iCell)) - &
          2.88/(tau + 1.0E-10_RKIND)*vw(i1,k,iCell)*(vvel(k-1,iCell) - vvel(k,iCell)) / &
          (zm(k-1,iCell) - zm(k,iCell)) + &
          2.88/(tau + 1.0E-10_RKIND)*gravity*(alphaT(k,iCell)*wt(i1,k,iCell) - betaS(k,iCell)* &
          ws(i1,k,iCell)) - 3.84*eps(i1,k,iCell)/(tau + 1.0E-10_RKIND)
      enddo

      do k = 1,nVertLevels
        eps(i2,k,iCell) = eps(i1,k,iCell) + dt*epstend(k)
        if (eps(i2,k,iCell) < 1.0E-10_RKIND) then
          eps(i2,k,iCell) = 1.0E-10_RKIND
        endif
        KE = 0.5_RKIND*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))
        length(k,iCell) = KE**1.5 / eps(i2,k,iCell)
!        length(k,iCell) = max(length(k,iCell),0.55_RKIND)
!        if(-ze(k,iCell) > boundaryLayerDepth(iCell)) then
!          length(k,iCell) = min(length(k,iCell),ze(k,iCell) - ze(k+1,iCell))
!        endif
      enddo
      length(nVertLevels+1,iCell) = 1e-15_RKIND
      length(1,iCell) = length(2,iCell) !1e-1
    enddo

  end subroutine build_length_from_dissipation_eq


  subroutine build_dissipation_lengths(nCells,nVertLevels,BVF,boundaryLayerDepth,alphaT,betaS, &
      wtsfc, wssfc)
    integer,intent(in) :: nCells, nVertLevels
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT, betaS
    real,dimension(nCells),intent(in) :: wtsfc, wssfc
    real,dimension(nVertLevels,nCells),intent(inout) :: BVF
    real,dimension(nCells) :: boundaryLayerDepth

    real :: invLen, l, len1, len2, lenmax, KE, integrandTop, integrandBot
    real :: len0, len2_1, len2_2, sfcBuoy, lengthT, bvfT

    integer :: k,iCell

    do iCell=1,nCells

      sfcBuoy = gravity*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))
      integrandTop = 0.0
      integrandBot = 0.0
      k=1
      do while (-ze(k,iCell) <= boundaryLayerDepth(iCell))
        !Add subplume bits here
        KE = sqrt(0.5_RKIND*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell)))
        integrandTop = integrandTop + 0.5_RKIND*KE*abs(ze(k,iCell)**2.0_RKIND - &
          ze(k+1,iCell)**2.0_RKIND)
        integrandBot = integrandBot + KE*(ze(k,iCell) - ze(k+1,iCell))
        k = k + 1
      enddo

      do k=2,nVertLevels
        KE = 0.5_RKIND*(u2(i2,k,iCell) + v2(i2,k,iCell) + w2(i2,k,iCell))

        bvfT = max(BVF(k,iCell),0.0_RKIND)
        len0 = 0.1_RKIND*integrandTop/(integrandBot + 1.0E-15_RKIND)
        len1 = len0*(0.4_RKIND*abs(ze(k,iCell)) / (0.4_RKIND*abs(ze(k,iCell)) + len0 + 1.0E-15_RKIND))
        len2_1 = (0.4_RKIND*sfcBuoy + sqrt((0.4_RKIND*sfcBuoy)**2.0 + 4.0_RKIND*KE**2.0*bvfT)) / &
                (2.0_RKIND * sqrt(KE) * bvfT + 1.0E-15)
        len2_2 = (0.4_RKIND*sfcBuoy - sqrt((0.4_RKIND*sfcBuoy)**2.0 + 4.0_RKIND*KE**2.0*bvfT)) / &
                (2.0_RKIND * sqrt(KE) * bvfT + 1.0E-15)
        if(len2_2 < 0) len2_2 = 1.0E15_RKIND

        if(ze(k,iCell) < -boundaryLayerDepth(iCell)) THEN
          lengthT = 1.0_RKIND / len0 + 1.0_RKIND / len1 + 1.0_RKIND / len2_1 + 1.0_RKIND
        else
          lengthT = 1.0_RKIND / (sqrt(KE / ((1.0E-15_RKIND + bvfT))) + 1.0E-15_RKIND)
        endif

        length(k,iCell) = 1.0_RKIND / lengthT

        !  length(k,iCell) = config_adc_length_multiple / sqrt(1.0_RKIND / (1.0E-15_RKIND + config_adc_tau_o * sqrt(KE) *  &
        !    0.4_RKIND * abs(ze(k,iCell))) + 1.0_RKIND / (1.0E-15_RKIND + config_adc_tau_o * sqrt(KE) &
        !    * 0.1 * integrandTop / (integrandBot + 1.0E-10_RKIND)) + &
        !    max(1.0E-15_RKIND, BVF(k,iCell)) * 0.01_RKIND / (1.0E-15_RKIND + KE))

!        if(KE > EPSILON) then
!          lenmax = 0.53_RKIND*sqrt(2.0_RKIND*KE / (1.0E-15_RKIND + BVF(k,iCell)))
!        else
!          lenmax = 1.0e6_RKIND
!        endif

!        len1 = 0.4_RKIND*abs(ze(k,iCell))
!        len2 = 0.2_RKIND*integrandTop / (integrandBot + 1.0E-10_RKIND)

!        length(k,iCell) = min(1.0_RKIND/(1.0_RKIND / len1 + 1.0_RKIND / len2),lenmax)
        length(k,iCell) = max(length(k,iCell),0.55_RKIND)
        if(-ze(k,iCell) > boundaryLayerDepth(iCell)) then
          length(k,iCell) = min(length(k,iCell),ze(k,iCell) - ze(k+1,iCell))
        endif
      enddo
      length(nVertLevels+1,iCell) = 1e-15_RKIND
      length(1,iCell) = length(2,iCell) !1e-1
    enddo

  end subroutine build_dissipation_lengths

  subroutine build_areaFraction_updraft_properties(nCells,nVertLevels)
  !builds the updraft area function

  integer,intent(in) :: nCells,nVertLevels
  integer :: iCell, k
  real :: Sw, w3av, lareaFraction, wtav, wsav

  do iCell = 1,nCells
    tumd(1,iCell) = 0.0_RKIND
    wumd(1,iCell) = 0.0_RKIND
    areaFraction(1,iCell) = 0.5_RKIND
    Mc(1,iCell) = 0.0_RKIND
    do k=2,nVertLevels
      w3av = 0.5_RKIND*(w3(i2,k-1,iCell) + w3(i2,k,iCell))

      Sw = w3av / (max(w2(i2,k,iCell)**1.5_RKIND,epsilon))
      lareaFraction = 0.5_RKIND - 0.5_RKIND*Sw / sqrt(4.0_RKIND + Sw**2.0_RKIND)

      if(lareaFraction < 0.01_RKIND) lareaFraction = 0.01_RKIND
      if(lareaFraction > 0.99_RKIND) lareaFraction = 0.99_RKIND

      areaFraction(k,iCell) = lareaFraction
      wumd(k,iCell) = sqrt(w2(i2,k,iCell) / (areaFraction(k,iCell) * &
        (1.0_RKIND - areaFraction(k,iCell))))
      Mc(k,iCell) = areaFraction(k,iCell)*(1.0_RKIND - &
        areaFraction(k,iCell)) * wumd(k,iCell)
    enddo
  enddo

  end subroutine build_areaFraction_updraft_properties

  subroutine calc_scalar_updraft_properties(nCells,nVertLevels,wtsfc, wssfc, alphaT,  &
    betaS, boundaryLayerDepth, tlev)

    integer,intent(in) :: nCells, nVertLevels, tlev
    real,dimension(nCells),intent(in) :: wtsfc, wssfc, boundaryLayerDepth
    real, dimension(nVertLevels,nCells) :: alphaT, betaS

    real :: wtav, McAv, sigav, tumdav, wumdav, sumdav, wbsfc, bld, wstar
    integer :: iCell,k

    do iCell=1,nCells
      do k=2,nVertLevels

        tumd(k,iCell) = wt(tlev,k,iCell) / (1.0E-12_RKIND + Mc(k,iCell))
        sumd(k,iCell) = ws(tlev,k,iCell) / (1.0E-12_RKIND + Mc(k,iCell))
      enddo

      wbsfc = gravity*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))

      if(wbsfc > 0.0) then
        wbsfc = gravity*(alphaT(1,iCell)*wtsfc(iCell) - betaS(1,iCell)*wssfc(iCell))
        wstar = (abs(boundaryLayerDepth(iCell)*wbsfc))**(1.0_RKIND/3.0_RKIND)
        w2t(1,iCell) = -0.3_RKIND*wstar * wtsfc(iCell)
        !Below FIXME!
        w2s(1,iCell) = 0.3_RKIND*wstar * wssfc(iCell)
      else
        w2t(1,iCell) = 0.0_RKIND
        w2s(1,iCell) = 0.0_RKIND
      endif

      !try new boundary condition derived from PDF
      sigav = 0.5_RKIND*(areaFraction(1,iCell) + areaFraction(2,iCell))
      wtav = 0.5_RKIND*(wt(tlev,1,iCell) + wt(tlev,2,iCell))
      McAv = 0.5_RKIND*(w2(tlev,1,iCell) + w2(tlev,2,iCell))
      w2t(1,iCell) = (1.0_RKIND - 2.0_RKIND*sigav)*wtav*sqrt(McAv) / (EPSILON + &
        sigav*(1.0_RKIND-sigav))

      do k=2,nVertLevels
        sigav = 0.5_RKIND*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
        tumdav = 0.5_RKIND*(tumd(k,iCell) + tumd(k+1,iCell))
        sumdav = 0.5_RKIND*(sumd(k,iCell) + sumd(k+1,iCell))
        wumdav = 0.5_RKIND*(wumd(k,iCell) + wumd(k+1,iCell))
        w2t(k,iCell) = sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*tumdav
        w2s(k,iCell) = sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)*wumdav**2.0_RKIND*sumdav
      enddo

    enddo
  end subroutine calc_scalar_updraft_properties

  subroutine calc_subplume_fluxes(nCells,nVertLevels,activeTracers,uvel,vvel,BVF,   &
    alphaT,betaS,dt)
  ! builds the subplume tendency terms

  integer,intent(in) :: nCells, nVertLevels
  real,dimension(nVertLevels,nCells),intent(in) :: alphaT,betaS
  real,intent(in) :: dt
  real,dimension(nVertLevels,nCells),intent(in) :: uvel,vvel
  real,dimension(2,nVertLevels,nCells),intent(in) :: activeTracers
  real,dimension(nVertLevels,nCells),intent(in) :: BVF

  real,dimension(nVertLevels+1,nCells) :: KspsUtend, KspsDtend
  real :: Uz, Vz, Tz, Sz, B, areaFractionAv,dz,integrandTop,integrandBot, Cval
  !calculate length

  integer :: iCell, k

  do iCell = 1,nCells

    lenspsU(1,iCell) = 0.0_RKIND
    lenspsD(1,iCell) = 0.0_RKIND
    KmU(1,iCell) = 0.0_RKIND
    KhU(1,iCell) = 0.0_RKIND
    KmD(1,iCell) = 0.0_RKIND
    KhD(1,iCell) = 0.0_RKIND
    Entrainment(1,iCell) = 0.0_RKIND
    Detrainment(1,iCell) = 0.0_RKIND

    do k = 2,nVertLevels-1
      dz = (zm(k-1,iCell) - zm(k,iCell))
      !need to add length scales for Up and Down
      Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dz
      Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dz
      Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dz
      Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dz

      if(BVF(k,iCell) <= 0.0_RKIND) then
        lenspsU(k,iCell) = dz
        lenspsD(k,iCell) = dz
      else
        lenspsU(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsU(k,iCell)/BVF(k,iCell)))
        lenspsD(k,iCell) = min(dz,0.76_RKIND*sqrt(KspsD(k,iCell)/BVF(k,iCell)))
      endif

      KmU(k,iCell) = 0.1_RKIND*lenspsU(k,iCell)*sqrt( KspsU(k,iCell) )
      KhU(k,iCell) = ( 1.0_RKIND+2.0_RKIND*lenspsU(k,iCell)/dz)*KmU(k,iCell)
      wt_spsU(k,iCell) =  -KhU(k,iCell)*Tz
      ws_spsU(k,iCell) =  -KhU(k,iCell)*Sz

      KmD(k,iCell) = 0.1_RKIND*lenspsD(k,iCell)*sqrt( KspsD(k,iCell) )
      KhD(k,iCell) = ( 1.0_RKIND+2.0_RKIND*lenspsD(k,iCell)/dz)*KmD(k,iCell)
      wt_spsD(k,iCell) = -KhD(k,iCell)*Tz
      ws_spsD(k,iCell) = -KhD(k,iCell)*Sz

      Entrainment(k,iCell) = Cww_E*areaFraction(k,iCell)*(1.0_RKIND- &
        areaFraction(k,iCell))*Mc(k,iCell) / ( lenspsD(k,iCell) + EPSILON )
      Detrainment(k,iCell) = Cww_D*areaFraction(k,iCell)*(1.0_RKIND- &
        areaFraction(k,iCell))*Mc(k,iCell) / ( lenspsU(k,iCell) + EPSILON )
    enddo

    lenspsU(nVertLevels,iCell) = lenspsU(nVertLevels-1,iCell)
    lenspsD(nVertLevels,iCell) = lenspsD(nVertLevels-1,iCell)
    do k=2,nVertLevels
      epsSPS(k,iCell) = (0.5_RKIND*(u2(i2,k,iCell) + v2(i2,k,iCell) + &
        w2(i2,k,iCell)))**1.5_RKIND/length(k,iCell)

      dz = (zm(k-1,iCell) - zm(k,iCell))
      if(k==2) then
        Cval = 3.96
      else
        Cval = (0.19_RKIND+0.51_RKIND*lenspsU(k,iCell)/dz)
      endif

      KspsUtend(k,iCell) = gravity*(alphaT(k,iCell)*wt_spsU(k,iCell) - &
        betaS(k,iCell)*ws_spsU(k,ICell)) + ((KmU(k-1,iCell) +                 &
        KmU(k,iCell))*(KspsU(k-1,iCell) - KspsU(k,iCell)) /         &
        (ze(k-1,iCell) - ze(k,iCell)) - (KmU(k,iCell) +             &
        KmU(k+1,iCell)) * (KspsU(k,iCell) - KspsU(k+1,iCell)) /     &
        (ze(k,iCell) - ze(k+1,iCell))) / dz -                                 &
        Cval*KspsU(k,iCell)**1.5_RKIND/(1.0E-15 + lenspsU(k,iCell)) +                     &
        eps(i1,k,iCell) / (2.0_RKIND*areaFraction(k,iCell))

      if(k==2) then
        Cval = 3.96_RKIND
      else
        Cval = (0.19_RKIND+0.51_RKIND*lenspsD(k,iCell)/dz)
      endif

      KspsDtend(k,iCell) = gravity*(alphaT(k,iCell)*wt_spsD(k,iCell) - &
        betaS(k,iCell)*ws_spsD(k,iCell)) + ((KmD(k-1,iCell) +                 &
        KmD(k,iCell))*(KspsD(k-1,iCell) - KspsD(k,iCell)) /         &
        (ze(k-1,iCell) - ze(k,iCell)) - (KmD(k,iCell) +             &
        KmD(k+1,iCell)) * (KspsD(k,iCell) - KspsD(k+1,iCell)) /     &
        (ze(k,iCell) - ze(k+1,iCell))) / dz -                                 &
        Cval*KspsD(k,iCell)**1.5_RKIND / lenspsD(k,iCell) +                   &
        eps(i1,k,iCell) / (2.0_RKIND*(1.0_RKIND - areaFraction(k,iCell)))
    enddo
  enddo

  do iCell=1,nCells
    do k=2,nVertLevels
      KspsU(k,iCell) = KspsU(k,iCell) + dt*KspsUtend(k,iCell)
      KspsD(k,iCell) = KspsD(k,iCell) + dt*KspsDtend(k,iCell)
    enddo
  enddo

  end subroutine calc_subplume_fluxes

  subroutine diagnose_momentum_fluxes(nCells,nVertLevels,numTracers,activeTracers,uvel,vvel, &
      alphaT,betaS,dt)
! This routine diagnoses all the horizontal related momentum flux components. All assume steady state
! follows a quasi structure function approach

    integer,intent(in) :: numTracers,nCells,nVertLevels
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT,betaS
    real,dimension(nVertLevels,nCells),intent(in) :: uvel,vvel
    real,dimension(numTracers,nVertLevels,nCells),intent(in) :: activeTracers
    real,intent(in) :: dt
    real,dimension(nVertLevels,nCells) :: taupt, taups, taupv
    real,dimension(nVertLevels+1,nCells) :: wttend, wstend, uwtend, vwtend, u2tend, v2tend, &
      uttend, vttend, ustend, vstend, uvtend
    real :: B, Kps, Kpsp1, diff, lenav, Uz, Vz, Tz, Sz, sigav, sumdav
    real :: dz, dzmid, tumdav, Ksps
    integer :: iCell, k

    do iCell=1,nCells
      !compute the TOMs first.
      do k=1,nVertLevels
        Ksps = 0.5_RKIND*((areaFraction(k,iCell)*KspsU(k,iCell) + (1.0_RKIND - &
          areaFraction(k,iCell))*KspsD(k,iCell)) + (areaFraction(k+1,iCell)* &
          KspsU(k+1,iCell) + (1.0_RKIND-areaFraction(k+1,iCell))* &
          KspsD(k+1,iCell)))
        Kps = 0.5_RKIND*(u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell))
        Kpsp1 = 0.5_RKIND*(u2(i1,k+1,iCell) + v2(i1,k+1,iCell) + w2(i1,k+1,iCell))
        lenav = 0.5_RKIND*(length(k,iCell) + length(k+1,iCell))
        diff = C_mom * sqrt(0.5_RKIND*(Kps + Kpsp1)) / lenav
        dz = ze(k,iCell) - ze(k+1,iCell)
        uw2(k,iCell) = -diff*(uw(i1,k,iCell) - uw(i1,k+1,iCell)) / dz
        vw2(k,iCell) = -diff*(vw(i1,k,iCell) - vw(i1,k+1,iCell)) / dz
        u2w(k,iCell) = -diff*(u2(i1,k,iCell) - u2(i1,k+1,iCell)) / dz
        v2w(k,iCell) = -diff*(v2(i1,k,iCell) - v2(i1,k+1,iCell)) / dz
        uvw(k,iCell) = -diff*(uv(i1,k,iCell) - uv(i1,k+1,iCell)) / dz

        diff = C_therm*sqrt(0.5*(Kps + Kpsp1)) / lenav
        uwt(k,iCell) = -diff*(ut(i1,k,iCell) - ut(i1,k+1,iCell)) / dz
        vwt(k,iCell) = -diff*(vt(i1,k,iCell) - vt(i1,k+1,iCell)) / dz
        uws(k,iCell) = -diff*(us(i1,k,iCell) - us(i1,k+1,iCell)) / dz
        vws(k,iCell) = -diff*(vs(i1,k,iCell) - vs(i1,k+1,iCell)) / dz
      enddo

      do k=2,nVertLevels
        dz = ze(k,iCell) - ze(k+1,iCell)
        dzmid = zm(k-1,iCell) - zm(k,iCell)
        Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dz
        Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dz

        Kps = sqrt((u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell)))
        B = gravity*(alphaT(k,iCell)*wt(i1,k,iCell) - betaS(k,iCell)*ws(i1,k,iCell))

        taupt(k,iCell) = Kps / (sqrt(2.0)*c_pt*length(k,iCell))
        taups(k,iCell) = Kps / (sqrt(2.0)*c_ps*length(k,iCell))

        !PDF formulation will allow for fairly easy generalization of vertical tracer fluxes to
        !all tracers, write in tumd, sumd, wumd form, not sure how to do horizontal fluxes though
        wttend1(k,iCell) = -(w2t(k-1,iCell) - w2t(k,iCell)) / dzmid

        wttend2(k,iCell) = -w2(i1,k,iCell)*(activeTracers(1,k-1,iCell) &
          - activeTracers(1,k,iCell)) / dzmid
        wttend3(k,iCell) = (1.0_RKIND - gamma1)*gravity*(alphaT(k,iCell) &
          *t2(i2,k,iCell) - betaS(k,iCell)*ts(i2,k,iCell))
        wttend4(k,iCell) = -alpha3/4.0*(ut(i1,k,iCell)*Uz + &
          vt(i1,k,iCell)*Vz)
        wttend5(k,iCell) = kappa_FL*(wt(i1,k-1,iCell) -       &
              wt(i1,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0

        wttend(k,iCell) = wttend1(k,iCell) + wttend2(k,iCell) +        &
          wttend3(k,iCell) + wttend4(k,iCell) + wttend5(k,iCell)

        wstend1(k,iCell) = -(w2s(k-1,iCell) - w2s(k,iCell)) / dz
        wstend2(k,iCell) = -w2(i1,k,iCell)*(activeTracers(2,k-1,iCell) &
          - activeTracers(2,k,iCell)) / dz
        wstend3(k,iCell) = (1.0_RKIND - gamma1)*gravity*(alphaT(k,iCell)* &
          ts(i2,k,iCell) - betaS(k,iCell)*s2(i2,k,iCell))
        wstend4(k,iCell) = -alpha3/4.0*(us(i1,k,iCell)*Uz + &
          vs(i1,k,iCell)*Vz)
        wsTend5(k,iCell) = kappa_FL*(ws(i1,k-1,iCell) -     &
          ws(i1,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

        wstend(k,iCell) = wstend1(k,iCell) + wstend2(k,iCell) + &
          wstend3(k,iCell) + wstend4(k,iCell) + wstend5(k,iCell)

        taupv(k,iCell) = Kps / (c_pv*length(k,iCell))

        uwtend1(k,iCell) = -(uw2(k-1,iCell) - uw2(k,iCell)) / dz
        uwtend2(k,iCell) = 0.5*((0.8-4.0/3.0*alpha1)*0.5*Kps**2.0 +  &
          (alpha1 - alpha2)*u2(i1,k,iCell) + (alpha1 +  &
          alpha2 - 2.0_RKIND)*w2(i1,k,iCell))*Uz
        uwtend3(k,iCell) = 0.5_RKIND*(alpha1 - alpha2)*    &
          uv(i1,k,iCell)*Vz
        uwtend4(k,iCell) = beta5*gravity*(alphaT(k,iCell)*   &
          ut(i1,k,iCell) - betaS(k,iCell)*us(i1,k,iCell))
        uwtend5(k,iCell) = - 2.0_RKIND*taupv(k,iCell)*uw(i1,k,iCell) + &
          kappa_FL*(uw(i1,k-1,iCell) - uw(i1,k+1,iCell)) / &
          (ze(k-1,iCell) - ze(k+1,iCell))**2.0

        uwtend(k,iCell) = uwtend1(k,iCell) + uwtend2(k,iCell) + &
          uwtend3(k,iCell) + uwtend4(k,iCell) + uwtend5(k,iCell)

        vwtend(k,iCell) = (-(vw2(k-1,iCell) - vw2(k,iCell)) / dz +   &
          0.5_RKIND*((0.8_RKIND-4.0_RKIND/3.0_RKIND*alpha1)*0.5_RKIND*Kps**2.0_RKIND +  &
          (alpha1 - alpha2)*v2(i1,k,iCell) + (alpha1 -   &
          alpha2 - 2.0_RKIND)*w2(i1,k,iCell))*Vz + 0.5_RKIND*(alpha1 &
          - alpha2)*uv(i1,k,iCell)*Uz + beta5*gravity*       &
          (alphaT(k,iCell)*vt(i1,k,iCell) - betaS(k,iCell)*vs(i1,k,iCell))) -            &
          taupv(k,iCell)*vw(i1,k,iCell) + kappa_FL*(vw(i1,k-1,iCell) &
          - vw(i1,k+1,iCell)) / (ze(k-1,iCell) - ze(k+1,iCell))**2.0

        uvtend(k,iCell) = (-(uvw(k-1,iCell) - uvw(k,iCell)) / dz + &
          (1.0_RKIND - 0.5_RKIND*(alpha1+alpha2))*(uw(i1,k,iCell)*Vz &
          + vw(i1,k,iCell)*Uz)) - taupv(k,iCell)*uv(i1,k,iCell) +           &
          kappa_VAR*(uv(i1,k-1,iCell) - uv(i1,k+1,iCell)) /       &
          (ze(k-1,iCell) - ze(k+1,iCell))**2.0_RKIND

        u2tend1(k,iCell) = -(u2w(k-1,iCell) - u2w(k,iCell)) / dz
        u2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 + alpha2 - &
          2.0_RKIND)*uw(i1,k,iCell)*Uz
        u2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1*vw(i1,k,iCell)*Vz
        u2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND*(1.0_RKIND-beta5)*B
        u2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + &
          taupv(k,iCell)*(Kps**2.0_RKIND/3.0_RKIND - u2(i1,k,iCell))

        u2tend(k,iCell) = u2tend1(k,iCell) + u2tend2(k,iCell) + &
          u2tend3(k,iCell) + u2tend4(k,iCell) + u2tend5(k,iCell)

        v2tend1(k,iCell) = -(v2w(k-1,iCell) - v2w(k,iCell)) / dz
        v2tend2(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 +alpha2 - &
          2.0_RKIND)*vw(i1,k,iCell)*Vz
        v2tend3(k,iCell) = - 2.0_RKIND/3.0_RKIND*alpha1*uw(i1,k,iCell)*Uz
        v2tend4(k,iCell) = 2.0_RKIND/3.0_RKIND*(1-beta5)*B
        v2tend5(k,iCell) = -2.0_RKIND/3.0_RKIND*eps(i1,k,iCell) + taupv(k,iCell)* &
          (Kps**2/3. - v2(i1,k,iCell))

        v2tend(k,iCell) = v2tend1(k,iCell) + v2tend2(k,iCell) + &
          v2tend3(k,iCell) + v2tend4(k,iCell) + v2tend5(k,iCell)

        !taupt = Kps / (2.0*c_pt*len(k,iCell))
        !taups = Kps / (2.0*c_ps*len(k,iCell))

        Tz = (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / dz
        Sz = (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / dz

        uttend(k,iCell) = (-(uwt(k-1,iCell) - uwt(k,iCell))/dz -  &
          uw(i1,k,iCell)*Tz - (1.0_RKIND - alpha3)*wt(i1,k,iCell) &
          *Uz) - ut(i1,k,iCell)*taupt(k,iCell)

        vttend(k,iCell) = (-(vwt(k-1,iCell) - vwt(k,iCell))/dz -  &
          vw(i1,k,iCell)*Tz - (1.0_RKIND - alpha3)*wt(i1,k,iCell) &
          *Vz) - vt(i1,k,iCell)*taupt(k,iCell)

        ustend(k,iCell) = (-(uws(k-1,iCell) - uws(k,iCell))/dz -  &
          uw(i1,k,iCell)*Sz - (1.0_RKIND - alpha3)*ws(i1,k,iCell) &
          *Uz) - us(i1,k,iCell)*taups(k,ICell)

        vstend(k,iCell) = (-(vws(k-1,iCell) - vws(k,iCell))/dz -  &
          vw(i1,k,iCell)*Sz - (1.0_RKIND - alpha3)*ws(i1,k,iCell)* &
          Vz) - vs(i1,k,iCell)*taups(k,iCell)

        t2(i2,k,iCell) = tumd(k,iCell)**2.0_RKIND*areaFraction(k,iCell)*&
          (1.0_RKIND-areaFraction(k,iCell))
        s2(i2,k,iCell) = sumd(k,iCell)**2.0_RKIND*areaFraction(k,iCell)*&
          (1.0_RKIND-areaFraction(k,iCell))
        ts(i2,k,iCell) = tumd(k,iCell)*sumd(k,iCell)* &
          areaFraction(k,iCell)*(1.0_RKIND-areaFraction(k,iCell))

      enddo
    enddo

    do iCell=1,nCells
      u2cliptend(:,iCell) = 0.0_RKIND
      v2cliptend(:,iCell) = 0.0_RKIND
      do k=2,nVertLevels
        u2(i2,k,iCell) = u2(i1,k,iCell) + dt*u2tend(k,iCell)
        if(u2(i2,k,iCell) < 0.0_RKIND) then
          u2cliptend(k,iCell) = -u2(i2,k,iCell)
          u2(i2,k,iCell) = 0.0_RKIND
        endif

        v2(i2,k,iCell) = v2(i1,k,iCell) + dt*v2tend(k,iCell)
        if(v2(i2,k,iCell) < 0.0_RKIND) then
          v2cliptend(k,iCell) = -v2(i2,k,iCell)
          v2(i2,k,iCell) = 0.0_RKIND
        endif

        uw(i2,k,iCell) = uw(i1,k,iCell) + dt*uwtend(k,iCell)
        vw(i2,k,iCell) = vw(i1,k,iCell) + dt*vwtend(k,iCell)
        uv(i2,k,iCell) = uv(i1,k,iCell) + dt*uvtend(k,iCell)
        ut(i2,k,iCell) = ut(i1,k,iCell) + dt*uttend(k,iCell)
        wt(i2,k,iCell) = (wt(i1,k,iCell) + dt*wttend(k,iCell)) / &
          (1.0_RKIND + dt*taupt(k,iCell))
        vt(i2,k,iCell) = vt(i1,k,iCell) + dt*vttend(k,iCell)
        us(i2,k,iCell) = us(i1,k,iCell) + dt*ustend(k,iCell)
        vs(i2,k,iCell) = vs(i1,k,iCell) + dt*vstend(k,iCell)
        ws(i2,k,iCell) = (ws(i1,k,iCell) + dt*wstend(k,iCell)) / &
          (1.0_RKIND + dt*taups(k,iCell))
        if(abs(wt(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: wt out of range, wt = $r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/wt(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif

        if(abs(ws(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: ws out of range, ws = $r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/ws(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif

        if(abs(u2(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: u2 out of range, u2 = $r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/u2(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif

        if(abs(v2(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: v2 out of range, v2 = $r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/v2(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif

      enddo

    enddo

  end subroutine diagnose_momentum_fluxes

  subroutine predict_turbulent_quantities(nCells, nVertLevels, nTracers,  &
        dt, activeTracers, uvel,vvel,alphaT,betaS)
    integer,intent(in) :: nCells, nVertLevels, nTracers
    real,intent(in) :: dt
    real,dimension(nTracers,nVertLevels,nCells),intent(in) :: activeTracers
    real,dimension(nVertLevels,nCells),intent(in) :: alphaT, betaS
    real,dimension(nVertLevels,nCells) :: uvel, vvel

    real,dimension(nVertLevels+1,nCells) :: w2tend
    real,dimension(nVertLevels,nCells) :: w3tend

    real :: Sw, St, Ss, Eav, Dav, sigav, sigavp1, wumdAv, tumdAv, sumdAv, wumdAvp1, tumdAvp1, sumdAvp1
    real :: Swup, KspsUav, KspsDav, KspsUavp1, KspsDavp1, KE, Mcav, lenav,u2av,v2av,w2av
    real :: w3temp, w3check, taups, taupt, mval, KEsps, Uz, Vz, dz

    integer :: iCell, k

    do iCell = 1,nCells
      do k=1,nVertLevels
        Eav = 0.5*(Entrainment(k+1,iCell) + Entrainment(k,iCell))
        Dav = 0.5*(Detrainment(k+1,iCell) + Detrainment(k,iCell))
        u2av = 0.5*(u2(i1,k,iCell) + u2(i1,k+1,iCell))
        v2av = 0.5*(v2(i1,k,iCell) + v2(i1,k+1,iCell))
        w2av = 0.5*(w2(i1,k,iCell) + w2(i1,k+1,iCell))

        sigav = 0.5*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
        wumdav = 0.5*(wumd(k,iCell) + wumd(k+1,iCell))
        tumdav = 0.5*(tumd(k,iCell) + tumd(k+1,iCell))
        sumdav = 0.5*(sumd(k,iCell) + sumd(k+1,iCell))
        KspsUav = 0.5*(KspsU(k,iCell) + KspsU(k+1,iCell))
        KspsDav = 0.5*(KspsD(k,iCell) + KspsD(k+1,iCell))
        Mcav = 0.5*(Mc(k,iCell) + Mc(k+1,iCell))
        lenav = 0.5*(length(k,iCell) + length(k+1,iCell))
        if(k==nVertLevels) then
          sigavp1 = 0.5*(areaFraction(k,iCell))
          wumdAvp1 = 0.5*(wumd(k,iCell))
          tumdAvp1 = 0.5*(tumd(k,iCell))
          sumdAvp1 = 0.5*(sumd(k,iCell))
        else
          sigavp1 = 0.5*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
          wumdAvp1 = 0.5*(wumd(k,iCell) + wumd(k+1,iCell))
          tumdAvp1 = 0.5*(tumd(k,iCell) + tumd(k+1,iCell))
          sumdAvp1 = 0.5*(sumd(k,iCell) + sumd(k+1,iCell))
        endif

        KEsps = sigav*KspsUav+ (1.0 - sigav)*KspsDav
        KE = sqrt((u2av+v2av+w2av) + 0.0*KEsps)

        dz = ze(k,iCell) - ze(k+1,iCell)
        !KE = sqrt(areaFraction(k,iCell)*KspsUav + (1.0 - areaFraction(k,iCell))*KspsDav)
        Swup = - 2.0/3.0*(KspsU(k,iCell) - KspsU(k+1,iCell)) / dz &
           - 2.0/3.0*KspsUav*(log(areaFraction(k,iCell)) -           &
          log(areaFraction(k+1,iCell))) / dz +         &
          2.0/3.0*(KspsD(k,iCell) - KspsD(k+1,iCell)) / dz        &
         + 2.0/3.0*KspsDav*(log(1.0-areaFraction(k,iCell)) -         &
          log(1.0-areaFraction(k+1,iCell))) /  dz

        w3tend1(k,iCell) = wumdav**3.0*(Eav*(3.0*sigav - 2.0) + Dav*(3.0*sigav - 1.0))
        w3tend2(k,iCell) =  wumdav**3.0*(6.0*sigav**2.0 - 6.0*sigav + 1)*             &
              (areaFraction(k,iCell)*(1.0_RKIND - areaFraction(k,iCell))*  &
              wumd(k,iCell) - areaFraction(k+1,iCell)*(1.0_RKIND - &
              areaFraction(k+1,iCell))*wumd(k+1,iCell))/dz
        w3tend3(k,iCell) = -1.5_RKIND*sigav*(1.0_RKIND - sigav)*(1.0_RKIND - 2.0_RKIND*sigav)* &
              wumdav**2.0_RKIND*((1.0_RKIND - 2.0_RKIND*areaFraction(k,iCell))*      &
              wumd(k,iCell)**2.0_RKIND - (1.0_RKIND -                                &
              2.0_RKIND*areaFraction(k+1,iCell))*wumd(k+1,iCell)**2.0_RKIND) / dz
        w3tend4(k,iCell) = 3.0_RKIND*(1.0_RKIND - 2.0_RKIND*sigav)*Mcav*wumdav*Swup -          &
              C_mom_w3*KE / (1E-15_RKIND+ sqrt(2.0_RKIND)*lenAv)*w3(i1,k,iCell)
        w3tend5(k,iCell) = 3.0_RKIND*gravity*(alphaT(k,iCell)*w2t(k,iCell) - &
              betaS(k,iCell)*w2s(k,iCell))*0.9_RKIND

        w3tend(k,iCell) = w3tend1(k,iCell) + w3tend2(k,iCell) + w3tend3(k,iCell) + &
              w3tend4(k,iCell) + w3tend5(k,iCell)

        if(k>1 .and. k < nVertLevels) then
           w3tend(k,iCell) = w3tend(k,iCell) + kappa_w3*(w3(i1,k-1,iCell) &
              - w3(i1,k+1,iCell)) / (zm(k-1,iCell) - zm(k+1,iCell))**2.0
        endif

      enddo

      do k=2,nVertLevels
        dz = (zm(k-1,iCell) - zm(k,iCell))

        sigav = 0.5*(areaFraction(k,iCell) + areaFraction(k-1,iCell))
        wumdav = 0.5*(wumd(k,iCell) + wumd(k-1,iCell))
        tumdav = 0.5*(tumd(k,iCell) + tumd(k-1,iCell))
        sumdav = 0.5*(sumd(k,iCell) + sumd(k-1,iCell))
        KspsUav = 0.5*(KspsU(k,iCell) + KspsU(k-1,iCell))
        KspsDav = 0.5*(KspsD(k,iCell) + KspsD(k-1,iCell))
        Mcav = 0.5*(Mc(k,iCell) + Mc(k-1,iCell))

        sigavp1 = 0.5*(areaFraction(k,iCell) + areaFraction(k+1,iCell))
        KspsUavp1 = 0.5*(KspsU(k,iCell) + KspsU(k+1,iCell))
        KspsDavp1 = 0.5*(KspsD(k,iCell) + KspsD(k+1,iCell))

        Uz = (uvel(k-1,iCell) - uvel(k,iCell)) / dz
        Vz = (vvel(k-1,iCell) - vvel(k,iCell)) / dz

        KEsps = areaFraction(k,iCell)*KspsU(k,iCell) +  &
                  (1.0_RKIND-areaFraction(k,iCell))*KspsD(k,iCell)
        KE = sqrt((u2(i1,k,iCell) + v2(i1,k,iCell) + w2(i1,k,iCell)) + 0.0*KEsps)
        Swup = gravity*alphaT(k,iCell)*tumd(k,iCell) - gravity*        &
          betaS(k,iCell)*sumd(k,iCell) - 2.0_RKIND/3.0_RKIND*(1.0_RKIND/        &
          areaFraction(k,iCell)*(sigAv*KspsUav - sigavp1*KspsUavp1) / dz -     &
          1.0_RKIND/(1.0_RKIND - areaFraction(k,iCell))*((1.0_RKIND - sigav)* &
          KspsDav - (1.0 - sigavp1)*KspsDavp1) / dz)

        w2tend1(k,iCell) = -wumd(k,iCell)**2.0_RKIND*(         &
          Entrainment(k,iCell) + Detrainment(k,iCell))
        w2tend2(k,iCell) = -(w3(i1,k-1,iCell) - w3(i1,k,iCell)) / dz
        w2tend3(k,iCell) = -C_1*KE / (1.0E-15_RKIND + sqrt(2.0_RKIND)* &
          length(k,iCell))*(w2(i1,k,iCell) - KE**2.0_RKIND/3.0_RKIND)
        w2tend4(k,iCell) = 2.0_RKIND*Mc(k,iCell)*Swup + 4.0_RKIND / &
          3.0_RKIND*C_2*areaFraction(k,iCell)*(1.0_RKIND -          &
          areaFraction(k,iCell))*wumd(k,iCell)*(gravity*     &
          alphaT(k,iCell)*tumd(k,iCell) - gravity*betaS(k,iCell)*      &
          sumd(k,iCell))
        w2tend5(k,iCell) = (1.0_RKIND/3.0_RKIND*alpha1 -            &
          alpha2)*(uw(i1,k,iCell)*Uz + vw(i1,k,iCell)*Vz)

        w2tend(k,iCell) = w2tend1(k,iCell) + w2tend2(k,iCell) + &
          w2tend3(k,iCell) + w2tend4(k,iCell) + w2tend5(k,iCell)
    enddo
  end do

    do iCell=1,nCells
      do k=2,nVertLevels
        w2(i2,k,iCell) = w2(i1,k,iCell) + dt*w2tend(k,iCell)
        if(w2(i2,k,iCell) < epsilon) then
          w2cliptend(k,iCell) = epsilon-w2(i2,k,iCell)
          w2(i2,k,iCell) = epsilon
        endif

        if(abs(w2(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: w2 out of range, w2 = $r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/w2(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif
      end do

      do k=1,nVertLevels
        w3check = (w2(i2,k,iCell) + w2(i2,k+1,iCell))**1.5
        w3(i2,k,iCell) = min(w3(i1,k,iCell) + dt*w3tend(k,iCell),w3check)

        if(abs(w3(i2,k,iCell)) > 1.0_RKIND) then
          call mpas_log_write("ERROR: w3 out of range, w2 = $r, location = $i, $i", &
            MPAS_LOG_CRIT,realArgs=(/w3(i2,k,iCell)/),intArgs=(/k,iCell/))
        endif
      enddo
    enddo

  end subroutine predict_turbulent_quantities

  subroutine update_mean_fields(dt,nCells,nVertLevels,numTracers,uvel, vvel, activeTracers,fCell)

    integer,intent(in) :: nCells, nVertLevels, numTracers
    real,dimension(nVertLevels,nCells),intent(out) :: uvel,vvel
    real,dimension(2,nvertLevels,nCells),intent(out) :: activeTracers
    real,intent(in) :: dt
    real,dimension(nCells),intent(in) :: fCell

    real :: utemp, vtemp
    integer :: iCell,k

    do iCell = 1,nCells
      do k = 1,nVertLevels
        utemp = uvel(k,iCell)
        vtemp = vvel(k,iCell)
        uvel(k,iCell) = uvel(k,iCell) - dt*(uw(i2,k,iCell) - uw(i2,k+1,iCell)) /  &
                  (ze(k,iCell) - ze(k+1,iCell)) !+ dt*fCell(iCell)*vtemp

        vvel(k,iCell) = vvel(k,iCell) - dt*(vw(i2,k,iCell) - vw(i2,k+1,iCell)) /  &
                  (ze(k,iCell) - ze(k+1,iCell)) !- dt*fCell(iCell)*utemp

        activeTracers(1,k,iCell) = activeTracers(1,k,iCell) - dt*(wt(i2,k,iCell) - &
          wt(i2,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell))
     activeTracers(2,k,iCell) = activeTracers(2,k,iCell) - dt*(ws(i2,k,iCell) - &
          ws(i2,k+1,iCell)) / (ze(k,iCell) - ze(k+1,iCell))
      enddo
    enddo

  end subroutine update_mean_fields

  subroutine ADC_main_loop(nCells,nVertLevels, nTracers, dt,activeTracers, uvel, vvel, BVF,  &
      uwsfc, vwsfc, wtsfc, wssfc, alphaT, betaS, fCell, boundaryLayerDepth)

    integer,intent(in) :: nCells, nVertLevels, nTracers
    real,intent(in) :: dt

    real,dimension(nTracers,nVertLevels,nCells),intent(inout) :: activeTracers
    real,dimension(nVertLevels,nCells),intent(inout) :: uvel, vvel, alphaT, betaS
    real,dimension(nCells),intent(in) :: uwsfc, vwsfc, wtsfc, wssfc, fCell
    real,dimension(nCells),intent(inout) :: boundaryLayerDepth
    real,dimension(nVertLevels,nCells),intent(inout) :: BVF
    integer :: niter, iIter,iCell,k

    real :: dt_small

    dt_small = config_adc_timestep
    niter = dt / dt_small

    i1 = 1
    i2 = 2
!    call construct_depth_coordinate(ssh,layerThick,nCells,nVertLevels)
    do iIter=1,niter
      call build_diagnostic_arrays(nCells,nVertLevels,BVF,wtsfc,wssfc,  &
        uwsfc,vwsfc,alphaT,betaS,boundaryLayerDepth)
      call predict_turbulent_quantities(nCells, nVertLevels, nTracers, dt_small, &
        activeTracers, uvel, vvel, alphaT, betaS)
      call diagnose_momentum_fluxes(nCells,nVertLevels,nTracers,activeTracers,uvel,vvel,alphaT,betaS,dt_small)
      call build_areaFraction_updraft_properties(nCells, nVertLevels)
      call calc_scalar_updraft_properties(nCells, nVertLevels,wtsfc, wssfc, &
                                alphaT, betaS, boundaryLayerDepth, i2)
      call calc_subplume_fluxes(nCells,nVertLevels,activeTracers,uvel,vvel, BVF,alphaT,betaS,dt_small)
     ! call build_dissipation_lengths(nCells,nVertLevels,BVF, boundaryLayerDepth,alphaT,betaS,wtsfc,wssfc)
      call build_length_from_dissipation_eq(nCells,nVertLevels,BVF,boundaryLayerDepth,alphaT,betaS, &
             uvel, vvel, dt)
      !      call dissipation_lengths2(nCells,nVertLevels,activeTracers,alphaT,betaS)
      call update_mean_fields(dt_small,nCells,nVertLevels, 2, uvel,vvel,activeTracers,fCell)
      call swap_time_levels
    enddo
  end subroutine ADC_main_loop

end module ocn_adc_mixing
