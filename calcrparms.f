      subroutine calcRectParms (i, k, j)

c///////////////////////////////////////////////////////////////////////
c/                                                                     /
c/    This subroutine calculates the different parameters needed for   /
c/    the evaluation of the transmission coefficients between          /
c/    intersecting and non-intersecting sides, in straightline         /
c/    geometry.                                                        /
c/                                                                     /
c/    On input, the routine needs to know the index -i- of the region  /
c/    under consideration, and the indices -k- and -j- of the two      /
c/    sides of the region. The routine then calculates the following   /
c/    quantities:                                                      /
c/                                                                     /
c/    L_i      : The length of the "FROM" side.                        /
c/    L_j      : The length of the "TO" side.                          /
c/    theta_ij : The angle between the sides (see report).             /
c/    alpha_j  : An angle needed in the case of non-intersecting sides /
c/    L_perp   : The perpendicular distance between non-intersecting   /
c/               sides.                                                /
c/    Written by John Mandrekas, GIT, (404) 894-7730                   /
c/    Modified by Dingkang Zhang (July 18, 2002)                             /
c///////////////////////////////////////////////////////////////////////

c/    Variable declarations, etc.

      implicit none
      integer i, j, k
      include 'locGeom.inc' 
      include 'neutGlob.inc'
      include 'consts.inc'

      integer idone, it, isum, iDiff_1, lindx, n_max,
     .   indx, nTriang,k1,k2,lindx1
      real beta, omega, l1, l2, tht, thm, dw, sum_theta, wangle

      wangle(l1,l2,tht) = atan(1.0 / (l1 / (l2*sin(tht)) - 
     .   1.0 / tan(tht)))
      
c/    Check if we have intersecting or non-intersecting sides and
c/    calculate the required parameters:
      
      nph_pnts = nph

      n_max = nSides(i)
      iDiff_1 =j-k 
      if(iDiff_1.lt.0)iDiff_1=iDiff_1+nSides(i)
      L_i = lside(k,i)
      L_j = lside(j,i)
      idone = 0

c/    Case of intersecting sides:

      if (iDiff_1.eq.1) then
	 i_geom = 1
	 indx = min(k,j)
	 theta_ij = angle(k,i)
	 alpha_j=theta_ij
	 idone = 1

c/    Case of non-intersecting sides with one side between them:

      else if ((iDiff_1.eq.2).and.idone.eq.0) then
         i_geom = 2
	       k1=k+1
	       if(k1.gt.nSides(i))k1=k1-nSides(i)
	       theta_ij = angle(k,i)
	       alpha_j = angle(k,i) + angle(k1,i) - pi
	       L_perp = lside(k1,i) * sin(theta_ij)

c/    Case of non-intersecting sides with | k-j | >= 3 :

      else if ((iDiff_1.gt.2).and.idone.eq.0) then
         nTriang = iDiff_1 - 2
	 k1=k+1
	 k2=k+2
	 if(k1.gt.nSides(i))k1=k1-nSides(i)
	 if(k2.gt.nSides(i))k2=k2-nSides(i)
	 l1 = lside(k1,i)
	 l2 = lside(k2,i)
	 thm = angle(k1,i)
	 omega = 0.0
	 do it = 1, nTriang
	    dw = wangle(l1, l2, thm)
	    omega = omega + dw
	    lindx = k + it + 1
	    lindx1=lindx+1
	    if(lindx.gt.nSides(i))lindx=lindx-nSides(i)
	    if(lindx1.gt.nSides(i))lindx1=lindx1-nSides(i)
	    l1 = lside(lindx,i) * sin(thm) / sin(dw)
	    l2 = lside(lindx1,i)
	    beta = pi - dw - thm
	    thm = angle(lindx,i) - beta
         enddo
	 theta_ij = angle(k,i) - omega
	 L_perp = l1 * sin(theta_ij)
	 sum_theta = 0.0
	 do isum = 1, iDiff_1
	    k1=k+isum-1
	    if(k1.gt.nSides(i))k1=k1-nSides(i)
	    sum_theta = sum_theta + angle(k1,i)
            enddo
	    alpha_j = sum_theta - (iDiff_1 - 1) * pi
	 idone = 1
      endif

      return
      end
