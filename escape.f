      subroutine escape

c/    This subroutine calculates the escape probabilities of each
c/    region, using the simple rational approximation result.
c/    It also calculates the quantity lambda(i,j)

c/    04/05/99, jm : Fixed bugs in calculation of areas
c/    01/17/02, jm : Changes for first collision source
      
      implicit none
      include 'neutGlob.inc'
      include 'consts.inc'
      include 'locGeom.inc'

c/    Local declarations:
      
      integer i, j, k, kk, kw, ktype, ns, nTriang,na
      real aterm, bterm, cterm, lp, omega, sarea, sp, tht, pterm, 
     .     P_0, pterm_0s, P_0s, argmnt, l1, l2, l3, l4, lprm, 
     .     sp1, sp2, testangl
      real sauer

c/    Calculate perimeter and area of each internal cell. For cells
c/    with more than four sides, break them up into triangles. 

      do i = 1, nCells
         ns = nSides(i)
	 nTriang = ns - 2
	 lp = 0.0
	 do k = 1, ns
	    lp = lp + lside(k,i)
         enddo
         perim(i) = lp
         if (ns.eq.3) then
	    l1 = lside(1,i)
	    l2 = lside(2,i)
	    l3 = lside(3,i)
	    sp = 0.5 * (l1 + l2 + l3) 
            sarea = Sqrt(sp * (sp-l1) * (sp-l2) * (sp-l3))

         else if (ns.eq.4) then

	    l1 = lside(1,i)
	    l2 = lside(2,i)
	    l3 = lside(3,i)
	    l4 = lside(4,i)
	    tht = angle(1,i)
	    lprm = Sqrt(l1**2 + l2**2 - 2.0*l1*l2*cos(tht))
	    sp1 = 0.5 * (l1 + l2 + lprm) 
	    sp2 = 0.5 * (l3 + l4 + lprm) 
            sarea = Sqrt(sp1 * (sp1-l1) * (sp1-l2) * (sp1-lprm)) +
     .             Sqrt(sp2 * (sp2-l3) * (sp2-l4) * (sp2-lprm))

	 else if (ns.GT.4) then

            sarea = 0.0
            cterm = lside(1,i)
	    omega = 0.0
	    do j = 1, nTriang - 1
	       aterm = cterm
	       bterm = lside(j+1, i)
	       tht = angle(j, i) - omega
	       cterm = Sqrt(aterm**2 + bterm**2 - 
     .            2.0 * aterm * bterm * cos(tht))
	       sp = 0.5 * (aterm + bterm + cterm)
	       sarea = sarea + 
     .	          Sqrt(sp * (sp - aterm) * (sp - bterm) * (sp - cterm))
	       argmnt = aterm * sin(tht) / cterm
	       testangl = aterm**2-(bterm**2+cterm**2)
	       if (testangl.LE.0.0) then
	         omega = asin(argmnt)
	       else
	         omega = pi - asin(argmnt)
               endif
            enddo

	    aterm = lside(ns-1,i)
	    bterm = lside(ns,i)
	    tht = angle(ns-1,i)
	    cterm = Sqrt(aterm**2+bterm**2-2.0*aterm*bterm*cos(tht))
	    sp = 0.5 * (aterm + bterm + cterm)
	    sarea = sarea + 
     .	      Sqrt(sp * (sp - aterm) * (sp - bterm) * (sp - cterm))

         endif
         
	 area(i) = sarea

c/    Calculate now the escape probability of region -i- :
	 
	 pterm = 4.0 * area(i) / (lmfp(i) * perim(i))

	 if (iescp.EQ.0) then
	    P_0 = 1.0 / (1.0 + pterm) ! Wigner's form
	 else 
	    P_0 = sauer(pterm)        ! Sauer-like approximation
         endif

	 pEscp0(i) = P_0
	 pEscp(i) = P_0 / (1.0 - A_cx(i) * (1.0 - P_0))

c/    Calculate terms needed for first collision source:

c/    Check first for volumetric sources in region -i- :
         
	 if (S_ext(i).GT.0.0) then

            if (ifrstcol.EQ.1) then
	       pterm_0s = 4.0 * area(i) / (lmfp0(i) * perim(i))

	       if (iescp.EQ.0) then
	          P_0s = 1.0 / (1.0 + pterm_0s)
	       else 
	          P_0s = sauer(pterm_0s)
               endif

	       pEscp0s(i) = P_0s

            else

	       pEscp0s(i) = pEscp0(i)

            endif

	 endif   

c/    Now, go over the sides of each cell. For now,  the escape 
c/    probability of first-collided neutrals is set equal to pEscp0(i),
c/    but the option remains for a treatment that will take into account
c/    first collision effects. 

	 do kk = 1, ns
	    k = adjCell(kk,i)
	    ktype = iType(k)
	    if (ktype.LE.1) then         ! internal cell or plasma 
	       pEscp0k(kk,i) = pEscp0(i)
	    else if (ktype.EQ.2) then      !  Wall
	       kw = k - (nCells + nPlasmReg)
	       pEscpw0(kw) = pEscp(i)
	       pEscpwf(kw) = pEscp(i)
	       pEscpws(kw) = pEscp(i)
	    endif
         enddo

c/    Calculate now the quantity lambda(i,j):
	 
	 do j = 1, nSides(i)
	   lambda(i,j) = lside(j,i) / perim(i)
         enddo

         if(nd0.eq.0)then
         do j=1, nSides(i)
            pEscpk(j,i)=pEscp(i)
            do k=1, nSides(i)
               lambdak(j,i,k,1)=lambda(i,k)
               do na=2, mExp
                  lambdak(j,i,k,na)=lambda(i,k)
                  if(na.lt.(2+inon))then
                    lambdak(j,i,k,na)=lambda(i,k)
     .                *colterm(j,na,i)/colterm(j,1,i)
                  endif
               enddo
            enddo
            if(iType(adjCell(j,i)).eq.2)then
              kw=adjcell(j,i)-(nCells+nPlasmReg)
              do k=1, nSides(i)
                 lambdawf(kw,k,1)=lambda(i,k)
                 lambdaws(kw,k,1)=lambda(i,k)
                 lambdaw0(kw,k,1)=lambda(i,k)
               do na=2, mExp
                  lambdawf(kw,k,na)=lambda(i,k)
                  lambdaws(kw,k,na)=lambda(i,k)
                  lambdaw0(kw,k,na)=lambda(i,k)
                  if(na.lt.(2+inon))then
                    lambdawf(kw,k,na)=lambdawf(kw,k,na)
     .                *stransw_f(kw,na)/stransw_f(kw,1)
                    lambdaws(kw,k,na)=lambdaws(kw,k,na)
     .                *stransw_s(kw,na)/stransw_s(kw,1)
                    lambdaw0(kw,k,na)=lambdaw0(kw,k,na)
     .                *stransw_0(kw,na)/stransw_0(kw,1)
                  endif
               enddo
              enddo
            endif
         enddo
         endif

      enddo

      return
      end

c//////////////////////////////////////////////////////////////////////

      real function sauer(x)

c/    This function calculates the Escape Probability from a region
c/    using a Sauer-like approximation, with an exponent equal to 2.09.
c/    This value of the exponent was determined from numerical simulations
c/    and comparison with Monte Carlo codes. 
c/    Reference: Roberto Rubilar Ph.D. thesis, June 2000, page 47.

      implicit none
      real x, n, term

      data n /2.0931773/

      term = (1.0 + x / n)
       
      sauer = (1.0/x) * (1.0 - 1.0 / term**n)

      return
      end
