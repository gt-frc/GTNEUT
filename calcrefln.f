      subroutine calcrefln

c/    This subroutine calculates the wall reflection coefficients
c/    for the case of irefl = 1. Notice that if the atomic number
c/    of a wall element is zero or negative, then the reflection 
c/    coefficient for this element is set equal to the input variable
c/    Rwall. This allows us to accomodate vacuum regions and pumping
c/    surfaces, even when we use irefl = 1

c/    Notice that in calculating the reflection coefficient, the
c/    projectile energy is taken to be the ion temperature of the 
c/    cell bounded by the wall segment.

c/    12/18/01, jm, Created

      implicit none
      include 'neutGlob.inc'

      integer iw, i_true, indx,kk,ithis
      real e0, am1, am2, z1, z2, rn, re

      if (irefl.EQ.0) then           ! No reflection model:
         do iw = 1, nWallSegm
	    refln(iw) = Rwall(iw)
	    refle(iw) = 1.0
	    fwabsorb(iw) = 1.0       ! No slow neutrals
         enddo

      else if (irefl.EQ.1) then 

         do iw = 1, nWallSegm
	    if (zwall(iw).GT.0.0) then
	       i_true = nCells + nPlasmReg + iw
	       indx = adjCell(1,i_true)
	       do kk=1,nSides(indx)
	          if(adjCell(kk,indx).eq.i_true) ithis=kk
	       enddo
	       e0 = tneut(indx,ithis)
               if(i_e0.eq.3)e0=1.5*e0
	       am1 = aneut
	       am2 = awall(iw)
	       z1 = zion
	       z2 = zwall(iw)
	       call reflect (e0, am1, am2, z1, z2, rn, re)
	       refln(iw) = rn
	       refle(iw) = re
	       Rwall(iw) = rn + (1.- rn) * (1.- fwabsorb(iw))
	    else if (zwall(iw).LE.0.0) then
	       refln(iw) = Rwall(iw)
	       refle(iw) = 1.0
	       fwabsorb(iw) = 1.0
            endif
         enddo
      endif
      return
      end
