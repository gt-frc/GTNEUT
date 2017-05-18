      subroutine svdegas(teev, tiev, e0ev, nem3, aion, aneut, 
     .   leh0, lchex, svione, svchx, svioni)

c/    This subroutine calculates various atomic physics cross sections
c/    using the DEGAS atomic physics data files. Prior to calling this
c/    routine, a call to DEGASREAD must be made to read the various
c/    coefficients and interpolation arrays.

c/    Written by John Mandrekas, GIT, 03/06/2001

c/    Input Variables:
c/    ---------------
c/    teev     : electron temperature (eV)
c/    tiev     : ion temperature (eV)
c/    e0ev     : neutral energy (eV)
c/    nem3     : electron density (m-3)
c/    aion     : mass number of ion species (1 for H, 2 for D, 3 for T)
c/    aneut    : mass number of neutral species
c/    leh0     : 1: e + h0 ionization rate dependent on ne
c/               2: e + h0 ionization rate independent of ne
c/    lchex    : if <=2, CX cross section depends on neutral energy
c                if = 3, CX cross section independent of neutral energy.

c/    Output variables:
c/    ----------------
c/    svione   : e + h0 -> h+ + 2e reaction rate (m^3/s)
c/    svioni   : h0 + h+ -> 2h+ + e reaction rate (m^3/s)
c/    svchx    : h+ + h0 -> h0 + h+ reaction rate (m^3/s)

      implicit none
      include 'degdat.inc'
      integer leh0, lchex
      real tiev, teev, e0ev, aion, aneut, nem3, svione, svchx, svioni

c/    Local Declarations
      integer inc, jd, ie, indxte, indxti, indxne, jindxde, 
     .        jindxp1, ie0, ke0, ke0p1
      real fuzz, necm3, zlogne, zlnte, tfact, z110den, rintd, zloge,
     .     zlogi, rinte, rinti, zti, zlnti, ze0, zlne0, zlog0, rint0,
     .     zcxrr0, zcxrr1, ziirr0, ziirr1
      double precision zsv11, zsv12, zsv21, zsv22, zsv1, zsv2, zsvion
    
      real chexrr, rionrr

      data fuzz / 1.e-10 /

      necm3 = 1.0E-6 * nem3

c/    Electron Impact Ionization rates:

c/    Find the exponent of the minimum value in the range, i.e.
c/    the density is within the range 10^jd - 10^(jd+0.5) where
c/    jd = [2*log(ne)] - 20 + 1

c/    Similarly, the minimum range for the electron temperature is 
c/    ie = 13 + (10/log(10)) * ln(Te), i.e: from 0.06 eV to 50 keV

      inc = 2*int(dhkpt(1)+fuzz) - 1
      zlogne = log10(necm3)
      jd = int(2*zlogne) - inc
      jd = max(1,jd)
      indxne = min(jd, mpdhm1)

      tfact = 10.0 / log(10.)
      zlnte = log(teev)
      ie = int(zlnte*tfact + 13.0)
      ie = max(1,ie)
      indxte = min(ie, mpehm1)

c/    Calculate interpolation coefficients for the density:

      z110den = max(zlogne, dhkpt(1))
      z110den = min(z110den, dhkpt(mpdh))
      rintd = ((z110den - dhkpt(indxne)) / 
     .        (dhkpt(indxne+1) - dhkpt(indxne)))

c/    Calculate interpolation cofficients for the electron temperature:

      zloge = zlnte
      zloge = max(ehkpt(1), zloge)
      zloge = min(zloge, ehkpt(mpeh))
      rinte = (zloge - ehkpt(indxte)) * tfact

 
c/    Check to see if rates are supposed to depend on ne:

      if (leh0.EQ.1) then
         jindxde = indxne
	 jindxp1 = indxne + 1
      else if (leh0.EQ.2) then
         jindxde = 1
	 jindxp1 = jindxde
      endif

c/    Calculate rates by performing a double logarithmic interpolation:
      
      zsv11 = wsveh(indxte, jindxde)
      zsv12 = wsveh(indxte, jindxp1)
      zsv21 = wsveh(indxte+1, jindxde)
      zsv22 = wsveh(indxte+1, jindxp1)

      zsv1 = zsv11 + rintd * (zsv12-zsv11)
      zsv2 = zsv21 + rintd * (zsv22-zsv21)

      zsvion = zsv1 + rinte * (zsv2 - zsv1)

      svione = 1.0E-06 * zsvion

c/    Charge exchange and ion impact ionization rates:
c/    ------------------------------------------------

c/    Index information for the ion temperature. Notice that the 
c/    temperature is scaled to take into account the hydrogenic
c/    isotope under consideration, assuming that the rates are
c/    derived for hydrogen (H). This version assumes a single 
c/    species background plasma.
      
c/    The temperature and energy range for these cross sections is:
c/    10^(ie-1)/10 = exp[(ie-1)*ln(10)/10], ie = 1, mpe (mpe = 48)
c/    which corresponds to the range from 1 eV to 50 keV.

      zti = tiev / aion
      zlnti = log(zti)
      ie = int(zlnti*ekptmpe1) + 1
      ie = max(1,ie)
      indxti = min(ie,mpem1)
     
c/    Interpolation factor:
      
      zlogi = zlnti
      zlogi = max(0.0, zlogi)
      zlogi = min(zlogi, ekpt(mpe))
      rinti = (zlogi - ekpt(indxti)) * tfact

c/    Do the same for the neutral energy parameters:

      ze0 = e0ev / aneut
      zlne0 = log(ze0)
      ie0 = int(zlne0*ekptmpe1) + 1
      ie0 = max(1,ie0)
      ke0 = min(ie0,mpem1)
      ke0p1 = ke0+1
      
      zlog0 = zlne0
      zlog0 = max(0.0, zlog0)
      zlog0 = min(zlog0, ekpt(mpe))
      rint0 = (zlog0 - ekpt(ke0)) * tfact

      zcxrr0 = chexrr(ke0,indxti,lchex,rinti)
      zcxrr1 = chexrr(ke0p1,indxti,lchex,rinti)

      svchx = 1.0E-6 * (zcxrr0 + rint0 * (zcxrr1-zcxrr0))

      ziirr0 = rionrr(ke0,indxti,rinti)
      ziirr1 = rionrr(ke0p1,indxti,rinti)

      svioni = 1.0E-6 * (ziirr0 + rint0 * (ziirr1-ziirr0))

      return
      end

c//////////////////////////////////////////////////////////////////////

      real function chexrr(ke0, it0, lchex, rinti)

c/    This function computes the charge exchange reaction rate
c/    for neutrals with energy index ke0 colliding with background
c/    plasma ions with energy index iti. 
c/    This version assumes CX between any two hydrogenic isotopes.

c/    Adapted from the DEGAS code by John Mandrekas, GIT, 03/09/2001
      
c/    ke0    : energy index of neutrals
c/    it0    : energy index of ions
c/    lchex  : flag determining whether svcx depends on E0
c/    rinti  : interpolation coefficient for ions

      implicit none

      include 'degdat.inc'

      integer jke0, ke0, it0, lchex, itp1
      real rinti

c/    If lchex = 3, then CX independent of neutral energy:
      
      itp1 = it0 + 1
      if (lchex.LE.2) jke0 = ke0
      if (lchex.EQ.3) jke0 = 1

      chexrr = svphcx(it0,jke0) + rinti * (svphcx(itp1,jke0) - 
     .         svphcx(it0,jke0))
      
      return
      end

c///////////////////////////////////////////////////////////////////////

      real function rionrr (ke0, it0, rinti)

c/    This function computes the ion impact ionization reaction 
c/    rate for neutrals with energy index ke0 colliding with 
c/    background plasma ions with energy index iti. 

c/    Adapted from the DEGAS code by John Mandrekas, GIT, 03/09/2001
      
c/    ke0    : energy index of neutrals
c/    it0    : energy index of ions
c/    rinti  : interpolation coefficient for ions

      implicit none

      include 'degdat.inc'

      integer jke0, ke0, it0, itp1
      real rinti
      
      itp1 = it0 + 1
      jke0 = ke0

      rionrr = svphe(it0,jke0) + rinti * (svphe(itp1,jke0) - 
     .         svphe(it0,jke0))
      
      return
      end
