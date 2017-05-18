      subroutine TransmCoeff(i, k, j)

c/    This subroutine calculates the first-flight transmission 
c/    coefficients, T(i)_k->j, na->nb.

c/    i    : region index
c/    k    : k-th side of region -i- (FROM)
c/    j    : j-th side of region -i- (TO)
c/ 
c/    tkji : Transmission probability FROM side k TO side j of region i 
c/           If side k is a wall segment, then tkji returns three values:
c/           tkji(1) -> fast (reflected) component      
c/           tkji(2) -> slow (wall temperature) component      
c/           tkji(3) -> Source (E_0) component, if g_ex(kw) > 0

      implicit none
      integer i, k, j, na, nb,np, np0
      include 'neutGlob.inc'
      include 'consts.inc'
      include 'locGeom.inc'

c/    Local variables:

      integer ie, kw, kcell
      real tiny, tkji(maxExp*maxExp)

      data tiny / 1.0e-06/

      idp1=idp
      inon1=inon

      kcell=adjCell(k,i)
      if(iType(kcell).eq.2)kw=kcell-(nCells+nPlasmReg)

      if(k.ne.j)call calcRectParms(i, k, j)


      if(k.eq.j.or.(abs(theta_ij-pi).le.tiny)) then      ! Set T_jji = 0.0
         do na=1,mExp
            do nb=1,mExp
               transm(k,j,na,nb,i)=0.0
            enddo
         enddo
         if(iType(kcell).eq.2)then
           do na=1,mExp
              do nb=1,mExp
                 transw_f(j,na,nb,kw)=0.0
                 transw_s(j,na,nb,kw)=0.0
                 transw_0(j,na,nb,kw)=0.0
              enddo
           enddo
         endif
         
         return
      endif

      if(iType(kcell).ge.1)then
         np=mExp
	   np0=mExp*mExp
      else
         np=mExp*mExp
      endif
 

c/    If the input variable i_e0 is equal to 1, then use the
c/    local mean free path (pre-2002 approach):

      if (i_e0.LE.1) then
         l_mfp = lmfp(i)
         call calcRect (tkji, iquad,np)
         do na=1,mExp
            do nb=1,mExp
               transm(k,j,na,nb,i)=tkji((na-1)*mExp+nb)
            enddo
         enddo
         if(iType(kcell).eq.2)then      
	   do na = 1, mExp
              do nb=1, mExp
                 transw_0(j,na,nb,kw)=tkji((na-1)*mExp+nb)
                 transw_s(j,na,nb,kw)=tkji((na-1)*mExp+nb)
                 transw_f(j,na,nb,kw)=tkji((na-1)*mExp+nb)
              enddo
	   enddo
         endif
	 return
      endif

c/    Assign local value of neutral mfp. It should correspond to the
c/    "FROM" region, to get the correct velocity.
c/
c/    If the "FROM" region is a wall segment, check to see if the 
c/    reflection model is on and if it applies to this particular 
c/    wall segment (zwall(kw) > 0), in order to take into account 
c/    fast, slow and source neutrals.


      if(iType(kcell).LT.2) then

         l_mfp = mfp(k,i)
         call calcRect (tkji, iquad, np)
	   do na=1,mExp
            do nb=1,mExp
               transm(k,j,na,nb,i)=tkji((na-1)*mExp+nb)
            enddo
         enddo

	else

	   if(((irefl.EQ.0).OR.(zwall(kw).LE.0)).and.(Rwall(kw).gt.0))then
	     l_mfp=lmfp(i)
           call calcRect (tkji, iquad, np0)
	     do na=1,mExp
              do nb=1,mExp
                 transm(k,j,na,nb,i)=tkji((na-1)*mExp+nb)
              enddo
           enddo
	   else
	     do na=1, mExp
	        do nb=1, mExp
	           transm(k,j,na,nb,i)=0
	        enddo
	     enddo
	   endif

	endif


      if (iType(kcell).EQ.2) then          ! Wall segment
	 
         if(inon.eq.0)then
	  if((irefl.EQ.0).OR.(irefl.EQ.1.AND.(zwall(kw).LE.0.0))) then
	      do na=1,mExp
               transw_f(j,1,na,kw)=transm(k,j,1,na,i)
               transw_s(j,1,na,kw)=transm(k,j,1,na,i)
            enddo
         else if ((irefl.EQ.1).AND.(zwall(kw).GT.0.0)) then ! Reflection on
	      l_mfp = mfp_wf(kw)
            call calcRect (tkji, iquad, np)
	      do na=1,mExp
               transw_f(j,1,na,kw)=tkji(na)
            enddo
	      l_mfp = mfp_ws(kw)
            call calcRect (tkji, iquad, np)
	      do na=1,mExp
               transw_s(j,1,na,kw)=tkji(na)
            enddo
         endif
         else
           if(irefl.eq.0.or.(zwall(kw).le.0)) then
             do na=1,mExp
                do nb=1,mExp
                   transw_f(j,na,nb,kw)=transm(k,j,na,nb,i)
                   transw_s(j,na,nb,kw)=transm(k,j,na,nb,i)
                enddo
              enddo
           else
             l_mfp = mfp_wf(kw)
             call calcRect (tkji, iquad, np0)
             do na=1,mExp
                do nb=1,mExp
                   transw_f(j,1,na,kw)=tkji((na-1)*mExp+nb)
                   enddo
             enddo
             l_mfp = mfp_ws(kw)
             call calcRect (tkji, iquad, np0)
             do na=1,mExp
                do nb=1,mExp
                   transw_s(j,1,na,kw)=tkji((na-1)*mExp+nb)
                enddo
             enddo
           endif
         endif

c/    If there is external flux on this wall segment, then calculate
c/    transmission coefficient for the source neutrals:
         
	 if (g_ex(kw).GT.0.0) then
	    l_mfp = mfp_w0(kw)
            call calcRect (tkji, iquad, np0)
	    do na=1,mExp
               transw_0(j,1,na,kw)=tkji(na)
            enddo
	 else
            do na=1,mExp
	       transw_0(j,1,na,kw)=0.0
            enddo
         endif	    

      endif

      return
      end
