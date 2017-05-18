      subroutine checkInput

c/    This routine checks the input for consistency, and also
c/    does some required units conversion and other chores:

      implicit none
      include 'neutGlob.inc'
      include 'comiou.inc'
      include 'consts.inc'

      integer i, ifail, ij, iw, ji, iThis, j, k, mxsds, nsds
      real oksum, sumang, anglDiff, diffl, tiny

      data tiny /1.0e-02/

c/    Make sure that the dimensions of the problem are within bounds:
c/    --------------------------------------------------------------
      if (nCells.GT.maxCell) then
	 write (6, '(1x, A30)') 'nCells > maxCell, Stopped!'
         call zstop(6,nout,2)
      endif
  
      if (nWallSegm.GT.maxWall) then
	 write (6, '(1x, A30)') 'nWallSegm > maxWall, Stopped!'
         call zstop(6,nout,2)
      endif

      if (nPlasmReg.GT.maxPlas) then
	 write (6, '(1x, A30)') 'nPlasmReg > maxPlas, Stopped!'
         call zstop(6,nout,2)
      endif

      mxsds = -1
      do i = 1, nCells
	 mxsds = max0(mxsds, nSides(i))
      enddo

      if (mxsds.GT.maxSides) then
	 write (6, '(1x, A30)') 'nSides > maxSides, Stopped!'
         call zstop(6,nout,3)
      endif


c/    Compatibility with older input files:
c/    ------------------------------------
      if (i_e0.EQ.0) i_e0 = 2
      if (v0fact.EQ.0.0) v0fact = 2.0
      if (iquad.EQ.0) iquad = 2
      if (nph.EQ.0) nph = 21

c/    Cannot have constant energy option if reflection model is on:

      if ((irefl.EQ.1).AND.(i_e0.EQ.1)) then
         write (6, '(1x,A30)') 'Cannot have i_e0 = 1 & irefl = 1'
         call zstop (6,nout,4)
      endif

c/    Check that the angles of each polygon sum up to (n-2)*pi:
     
      ifail = 0
      do i = 1, nCells
         nsds = nSides(i)
	 oksum = (nsds - 2) * 180.0
	 sumang = 0.0
	 do k = 1, nsds
	   sumang = sumang + angle(k,i)
         enddo
         anglDiff = abs(sumang-oksum)
	 if (anglDiff.gt.tiny) then
	   ifail = ifail + 1
	   write (6,1000) i
         endif
      enddo

      if (ifail.ne.0) call zstop (6,nout,1)

c/    Check for consistency of adjacent sides:
      do i = 1, nCells
        do j = 1, nSides(i)
	  ij = adjCell(j,i)
	  if (iType(ij).gt.1) then
	    ji = adjCell(1,ij)
	    if (ji.ne.i) then
	      write (6,1100) i, j
	      ifail = ifail + 1
            endif
	  else if (iType(ij).le.1) then
	    iThis = 0
	    do k = 1, nSides(ij)
	      if (adjCell(k,ij).eq.i) iThis = k 
            enddo
	    if (iThis.eq.0) then
	      write (6,1100) i, j
	      ifail = ifail + 1
            endif
	    if (iType(ij).eq.0) then
	      diffl = abs(lside(j,i) - lside(iThis,ij))
	      if (diffl.gt.tiny) then
		write (6,1200) j, i, iThis, ij
		ifail = ifail + 1
              endif
            endif
          endif
        enddo
      enddo

      if (ifail.ne.0) call zstop(6,nout,1)

c/    Convert all the angles in rads:

      do i = 1, nCells
	 do k = 1, nSides(i)
	    angle(k,i) = pi * angle(k,i) / 180.00
         enddo
      enddo

c/    Multiply lengths by scale factor:
c/   
      if (scalFact.GT.0.0) then
         do i = 1, nCells
	    do k = 1, nSides(i)
	       lside(k,i) = scalFact * lside(k,i)
            enddo
         enddo
      endif

c/    Normalization of external driving forces:
c/    ----------------------------------------
      srcNormFact = -1.0e30
      do i = 1, nCells
	srcNormFact = amax1(S_ext(i), srcNormFact)
      enddo

      do iw = 1, nWallSegm
	srcNormFact = amax1(g_ex(iw), g_ion(iw), srcNormFact)
      enddo

c/    After finding the maximum, normalize:
c/    ------------------------------------
      if (srcNormFact.gt.0.0) then
	do i = 1, nCells
	  S_ext(i) = S_ext(i) / srcNormFact
        enddo

        do iw = 1, nWallSegm
          g_ex(iw) = g_ex(iw) / srcNormFact
	  g_ion(iw) = g_ion(iw) / srcNormFact
          if(twall(iw).eq.0)twall(iw)=1.0
        enddo

      endif

c/    Set number of the angular flux expansion functions
      if(idp.eq.0)then
        nExp=1
      else if(idp.eq.1) then
        nExp=3
      endif
      if(inon.ne.1)inon=0
      mExp=nExp+inon

      if(nd0.le.0)nd0=0

      do i=1,nCells+nPlasmReg
         do j=1,nSides(i)
            tneut(i,j)=ionTemp(i)
         enddo
      enddo

      if(neitr.eq.0)neitr=1

 1000 format ('Inconsistent angles in region ', i5)
 1100 format ('Wrong neighbors! Cell ', i5, ' side ', i2)
 1200 format ('Length of side ', i2, ' of cell ', i5, '  must be equal',
     .          ' to side ', i2, ' of cell ', i5)

      return
      end
