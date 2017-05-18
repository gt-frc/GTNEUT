      subroutine calcxswms(teev, tiev, tnev, nem3, svion, svcx, svrec,
     .    svel, sveln)

c/    This routine does a double log interpolation to calculate various 
c/    reaction rates based on WMS data.

c/    Input parameters:
c/    ----------------
c/    teev    : electron temperature (eV)
c/    tiev    : ion temperature (eV)
c/    tnev    : neutral temperature (eV)
c/    nem3    : electron density (/m3)
      
c/    Calculated Quantities (all rates in m3/s):
c/    -----------------------------------------
c/    svion(ne,Te)   : electron impact ionization rate 
c/    svrec(ne,Te)   : recombination rate
c/    svcx(Tn,Ti)    : charge exchange rate
c/    svel(Tn,Ti)    : ion-neutral elastic scattering rate
c/    sveln(Tn)      : neutral-neutral scattering cross section

      implicit none
      include 'wmsdata.inc'

      integer itep, item, inep, inem, itip, itim, itnm, itnp, i 
      real teev, tiev, tnev, nem3, svion, svcx, svrec, svel, sveln

      real alti, altn, alte, alne, ff, gf

      alte = alog10(teev)
      alti = alog10(tiev)
      alne = alog10(nem3)
      altn = alog10(tnev)

c/    Get bounds:
c/    ----------
      if(alte.lt.tint(1)) alte = tint(1)
      if(alti.lt.tint(1)) alti = tint(1)
      if(alte.gt.tint(nti)) alte = tint(nti)
      if(alti.gt.tint(nti)) alti = tint(nti)
      if(altn.lt.tnnt(1)) altn = tnnt(1)
      if(altn.gt.tnnt(ntn)) altn = tnnt(ntn)

      if(alne.lt.znint(1)) alne = znint(1)
      if(alne.gt.znint(nne)) alne = znint(nne)

c/    Find the region (i-1, i) where the electron temperature lies:
c/    ------------------------------------------------------------
      do i = 1, nti
         itep = i
         item = i-1
         if(tint(i).gt.alte) go to 20
      enddo

   20 ff = (alte-tint(item)) / (tint(itep)-tint(item))

c/    Find the region (i-1, i) where the electron density lies:
c/    --------------------------------------------------------
      do i = 1, nne
         inep = i
         inem = i-1
         if(znint(i).gt.alne) go to 40
      enddo

   40 gf = (alne-znint(inem)) / (znint(inep)-znint(inem))


c/    Double interpolation in ne and Te for svion  and svrec:
c/    ------------------------------------------------------
      svion  = (1.-ff)*(1.-gf)*eion(inem,item)
     .       + (1.-ff)*    gf *eion(inep,item)
     .       +     ff *(1.-gf)*eion(inem,itep)
     .       +     ff *    gf *eion(inep,itep)

      svrec  = (1.-ff)*(1.-gf)*rec(inem,item)
     .       + (1.-ff)*    gf *rec(inep,item)
     .       +     ff *(1.-gf)*rec(inem,itep)
     .       +     ff *    gf *rec(inep,itep)

      svion = 10.0 ** svion
      svrec = 10.0 ** svrec

c/    Do the same thing for Tn and Ti:
c/    -------------------------------

c/    Find the region (i-1, i) where the neutral temperature lies:
c/    ------------------------------------------------------------
      do i = 1, ntn
         itnp = i
         itnm = i-1
         if(tnnt(i).gt.altn) go to 50
      enddo

   50 ff = (altn-tnnt(itnm)) / (tnnt(itnp)-tnnt(itnm))

c/    Find the region (i-1, i) where the ion temperature lies:
c/    --------------------------------------------------------
      do i = 1, nti
         itip = i
         itim = i-1
         if(tint(i).gt.alti) go to 60
      enddo

   60 gf = (alti-tint(itim)) / (tint(itip)-tint(itim))

c/    Single interpolation in Tn for n-n elastic scattering:
c/    -----------------------------------------------------
      sveln = elastn(itnm) + ff*(elastn(itnp)-elastn(itnm))

c/    Double interpolation in Tn and Ti for svel  and svcx:
c/    ----------------------------------------------------
      svel   = (1.-ff)*(1.-gf)*elast(itnm,itim)
     .       + (1.-ff)*    gf *elast(itnm,itip)
     .       +     ff *(1.-gf)*elast(itnp,itim)
     .       +     ff *    gf *elast(itnp,itip)

      svcx   = (1.-ff)*(1.-gf)*cx(itnm,itim)
     .       + (1.-ff)*    gf *cx(itnm,itip)
     .       +     ff *(1.-gf)*cx(itnp,itim)
     .       +     ff *    gf *cx(itnp,itip)

      
      sveln = 10.0 ** sveln
      svel = 10.0 ** svel
      svcx = 10.0 ** svcx

      return
      end
