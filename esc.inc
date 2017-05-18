c/////////////////////// File esc.inc //////////////////////////////

c/    This file contains variables needed in the calculation
c/    of the various escape probabilities.

      integer maxFEM, maxDiv
      parameter(maxDiv=10, maxFEM=(maxDiv+1)*maxDiv*maxSides/2+1)
      real xp(maxSides),yp(maxSides),xc,yc, subarea(maxSides),
     .    sublength(maxSides), 
     .    xa2(maxSides),ya2(maxSides),
     .    xya(maxSides), xa(maxSides), ya(maxSides),
     .    xl2(maxSides), yl2(maxSides), xyl(maxSides),
     .    xl(maxSides), yl(maxSides)
      double precision A1(maxFEM,maxFEM),B1(maxFEM,maxFEM)
	  double precision ST(maxFEM,6*maxSides),S0(maxFEM,6*maxSides)   
	  real    xt(maxFEM),yt(maxFEM)
	 
	  common /esc/A1,B1,S0,ST

	  real lam0,l01,l12,phi012
	  common /secsource/lam0,l01,l12,phi012


