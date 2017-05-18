c/    This file contains variables for the various atomic cross
c/    sections by W.M. Stacey and E.W. Thomas
c/    Written by John Mandrekas, 07/28/97
c/    Modified to include new data, 08/07/97, jm

      integer nne, nte, nti, ntn
      parameter (nti = 5, nte = 5, nne = 5, ntn = 3)

      real elast(ntn,nti), cx(ntn,nti), eion(nne,nte), rec(nne,nte),
     .   elastn(ntn), tint(nti), tent(nte), tnnt(ntn), znint(nne)

      common /adata/ tint, tent, tnnt, znint, elastn, elast, cx, eion,
     .    rec
