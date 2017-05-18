c/    This file includes declarations for the DEGAS atomic physics data

      integer mpe, mpeh, mpdh, mpem1, mpehm1, mpdhm1
      parameter (mpe=48, mpeh=60, mpdh=15, mpehm1=mpeh-1, 
     .    mpdhm1=mpdh-1, mpem1=mpe-1)

      double precision wsveh(mpeh,mpdh), wsveh0(mpeh,mpdh),
     .    welms1(mpeh,mpdh), welms2(mpeh,mpdh), pne31(mpeh,mpdh),
     .    pne32(mpeh,mpdh), pne21(mpeh,mpdh) ,pne22(mpeh,mpdh),
     .    svphe(mpe,mpe), svphcx(mpe,mpe)
      real  dhkpt(mpdh), ehkpt(mpeh), ekpt(mpe)
      
      real erefmin, erefmax, ekptmpe1

      common /degas1/ wsveh, wsveh0, welms1 ,welms2, pne31, pne32, 
     .    pne21, pne22, svphe, svphcx, dhkpt, ehkpt, ekpt, erefmin,
     .    erefmax, ekptmpe1
