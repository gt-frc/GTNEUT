      block data xsdata

c/    This block data subprogram, initializes the cross section data
c/    compiled by W.M. Stacey and E.W. Thomas
c/    Written by John Mandrekas, 07/29/97
c/    08/07/97, jm, updated with new atomic data

c/    Description of variables:
c/    (note: units are MKS and temperatures in eV)
c/    --------------------------------------------
c/    tint          : alog10 of interpolation plasma temperatures
c/    tnnt          : alog10 of interpolation neutral temperatures
c/    znint         : alog10 of interpolation electron densities
c/    elast(Tn,Ti)  : alog10 of neutral-ion elastic scattering rate
c/    elastn(Tn)    : alog10 of neutral-neutral elastic scattering rate
c/    cx(Tn,Ti)     : alog10 of ion charge exchange rate
c/    eion(ne,Te)   : alog10 of electron impact ionization rate
c/    rec(ne,Te)    : alog10 of recombination rate

      include 'wmsdata.inc'

      data tint /
     . -1.0000E+00,  0.0000E+00,  1.0000E+00,  2.0000E+00,  3.0000E+00
     .   /

      data tnnt /
     .  0.0000E+00,  1.0000E+00,  2.0000E+00
     .   /

      data elast /
     . -1.3569E+01, -1.3337E+01, -1.3036E+01, -1.3569E+01, -1.3337E+01,
     . -1.3036E+01, -1.3337E+01, -1.3167E+01, -1.3046E+01, -1.3036E+01,
     . -1.3046E+01, -1.2796E+01, -1.3036E+01, -1.3046E+01, -1.2796E+01
     .   /

      data cx /
     . -1.4097E+01, -1.3921E+01, -1.3553E+01, -1.4097E+01, -1.3921E+01,
     . -1.3553E+01, -1.3921E+01, -1.3824E+01, -1.3538E+01, -1.3553E+01,
     . -1.3538E+01, -1.3432E+01, -1.3553E+01, -1.3538E+01, -1.3432E+01
     .   /

      data eion /
     . -2.8523E+01, -2.8523E+01, -2.8523E+01, -2.8523E+01, -2.8523E+01,
     . -1.7745E+01, -1.7745E+01, -1.7745E+01, -1.7745E+01, -1.7745E+01,
     . -1.3620E+01, -1.3620E+01, -1.3620E+01, -1.3620E+01, -1.3620E+01,
     . -1.3097E+01, -1.3097E+01, -1.3097E+01, -1.3097E+01, -1.3097E+01,
     . -1.3301E+01, -1.3301E+01, -1.3301E+01, -1.3301E+01, -1.3301E+01
     .   /

      data rec /
     . -1.7523E+01, -1.6745E+01, -1.5155E+01, -1.4222E+01, -1.3301E+01,
     . -1.8409E+01, -1.8398E+01, -1.8398E+01, -1.7886E+01, -1.7000E+01,
     . -1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01, -1.9398E+01,
     . -2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01, -2.0155E+01,
     . -2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01, -2.1000E+01
     .   /

      data elastn /
     . -1.4569E+01, -1.4167E+01, -1.3796E+01
     .   /

      data znint /
     .  1.6000E+01,  1.8000E+01,  2.0000E+01,  2.1000E+01,  2.2000E+01
     .   /

      end
