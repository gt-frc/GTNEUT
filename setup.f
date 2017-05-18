      subroutine setup

c/    Sparse Matrix version (using the UMFPACK library)

c/    This subroutine sets up the elements of the solution matrix:
c/
c/    12/14/01: Changed definition of npos(i,j) so that second index
c/              now refers to the j-th side of cell i, to reduce
c/              memory requirements.
c/
c/    01/14/02: Eliminated core plasma and wall fluxes from
c/              solution vector, to conserve memory and to
c/              make implementation of energy-dependent wall
c/              reflection easier. Wall and core plasma fluxes
c/              are now constructed from internal cell fluxes.
c/
c/    06/06/03: Sparse matrix version
c/    06/06/03: Correct treatment of ion fluxes

      implicit none
      include 'neutGlob.inc'
      include 'consts.inc'

c/    Local variables:
 
      integer i, iCount, i_row, i_col, i_elmnt, j, jj, jtype, jw, 
     .        j_rel, k, kk, ktype, kw, k_rel, kw_rel, ll, nfrst,mm,
     .        na, nb, ll0, i_sparsec(maxElmns),iflag,nab,nb0,na0

      real cterm, uterm, cterm_p, uterm_p, cterm_w0, uterm_w0, uterm_wf,
     .     cterm_wf, uterm_ws, cterm_ws 

      data nfrst /0/

c/    INITIALIZATIONS, etc:
c/    --------------------
c/    Count equations and non-zero matrix elements. 
c/    Initialize the position holder and the vectors storing the 
c/    sparse matrix triplet information:   

      if (nfrst.ge.0) then

c/    Count equations and unknowns:

         nEqs = 0
         do i = 1, nCells
            nEqs=nEqs+nSides(i)*mExp
         enddo

c/    Find out the number of non-zero matrix elements. Recall that the
c/    equation at each interface couples the flux from this interface to
c/    the adjacent cell with the fluxes from all adjacent interfaces to
c/    the cell. If the adjacent cell is a wall segment or a core plasma
c/    region, then the flux from that element to the cell is expressed
c/    in terms of the flux from the cell to the element, so it should not
c/    be counted.


c/    POSITION INDEX:
 
c/    Establish a one-to-one correspondence between the 
c/    flux variable g(i->j) and the solution vector x(n), so that:
c/    g(i->j) == x( npos(i,j) ). Notice, that unlike the notation
c/    used in the theory, here the second index for both the flux
c/    and the position array refer to the j-th side of cell i:
      
         iCount = 0
         do i = 1, nCells
            do j = 1, nSides(i)
               do na=1,mExp
	          iCount = iCount + 1
	          npos(i,j,na) = iCount
               enddo
            enddo
         enddo
         
c/    Zero out sparse matrix arrays (array triplet):
         
	 do i = 1, maxElmns
	    a_sparse(i) = 0.0
	    i_sparse(i) = 0
            i_sparsec(i)=0
	    i_sparse(maxElmns+i) = 0
         enddo

c/    Zero out the RHS array:
	 
	 do i = 1, nEqs
	    b(i) = 0.0
         enddo

         nfrst = 1

      endif               ! END of IF loop for nfrst

c/    MAIN PART OF SETUP ROUTINE: 
c/    --------------------------
c/
c/    Place the non-zero elements of the sparse matrix M(nEqs,nEqs) into
c/    the holding vector a_sparse(nElmnts) using the storage scheme used
c/    by the UMFPACK solver. The index array i_sparse holds the row and
c/    column position of each a_sparse element, i.e., i_sparse(i) is the
c/    row index of a_sparse(i) and i_sparse(nElmnts+i) is the column index
c/    of a_sparse(i). 
c/    
      i_elmnt = 0

      do i = 1, nCells
	 do jj = 1, nSides(i)
            do nb=1,mExp
	       i_row = npos(i,jj,nb)
	       j = adjCell(jj,i)
	       jtype = iType(j)

c/    We first consider the coefficient for the flux exiting cell
c/    i from side jj. Due to our adopted index assignment scheme, 
c/    this is the diagonal element of row i_row, i.e. M(i_row,i_row).
c/
c/    If the cell j that is adjacent to side jj of cell i is an internal
c/    cell, then this coefficient is equal to 1. If the adjacent cell j 
c/    is a core plasma region or a wall segment, then we must take into 
c/    account the contribution from the elimination of fluxes FROM these 
c/    elements in terms of the flux from i onto these elements:

	       if(jtype.EQ.0.or.(nb.gt.(1+inon)))then 
c		         regular plasma cell
c                  or the higher partial current moments 
	       
	          i_elmnt = i_elmnt + 1
	          a_sparse(i_elmnt) = 1.0
	          i_sparse(i_elmnt) = i_row
	          i_sparsec(i_elmnt) = i_row
            
	       else if (jtype.EQ.1) then        ! core plasma cell
               
	         j_rel = j - nCells

c/    Only collided term from a side to itself:
                 if(nb.eq.1)then
                   iflag=1
                 else
                   iflag=-1
                 endif
                 nb0=nab(nb,nb)

	         cterm_p = colterm(jj,1,i)*A_cxk(jj,i)* 
     .	            *pEscpk(jj,i)*lambdak(jj,i,jj,nb0)*iflag  

	         i_elmnt = i_elmnt + 1

	         a_sparse(i_elmnt) = 1.0 - albedo(j_rel) * cterm_p

	         i_sparse(i_elmnt) = i_row
	         i_sparsec(i_elmnt) = i_row

	       else if (jtype.EQ.2) then        ! wall segment
               
	         jw = j - (nCells +  nPlasmReg)
                 nb0=nab(nb,nb)
                 if(nb.eq.1)then
                   iflag=1
                 else
                   iflag=-1
                 endif

c/    No uncollided fluxes from this element to itself. 
c/    Calculate the collided contributions:

                 if(nb.eq.1)then
	           cterm_w0 = stransw_0(jw,1)*
     .		   A_cxw0(jw)*pEscpw0(jw)*lambdaw0(jw,jj,nb0)*iflag
                 endif

	         cterm_wf = stransw_f(jw,1)*refln(jw) * 
     .		  A_cxwf(jw)*pEscpwf(jw)*lambdawf(jw,jj,nb0)*iflag

	         cterm_ws =  stransw_s(jw,1) * (1.- refln(jw)) * 
     .		  (1.- fwabsorb(jw)) * A_cxws(jw) * 
     .            pEscpws(jw)*lambdaws(jw,jj,nb0)*iflag

	         i_elmnt = i_elmnt + 1

	         a_sparse(i_elmnt) = 1.0 - cterm_ws - cterm_wf

	         i_sparse(i_elmnt) = i_row
	         i_sparsec(i_elmnt) = i_row


c/    In case of wall segment, we should add any wall contributions to 
c/    the RHS vector b(i_row) from external and ion fluxes:

                 if(nb.eq.1)then
	           b(i_row)=b(i_row)+g_ex(jw)*cterm_w0 + 
     .	              g_ion(jw) * (cterm_wf + cterm_ws)
                 endif

	      endif        ! END of IF loop for diagonal element

c/    Now, we proceed with fluxes from all sides, taking care
c/    not to double-count the core plasma and wall contributions from
c/    side jj that we assigned to the diagonal element already.
	    
	      do kk = 1, nSides(i)
	         k = adjCell(kk,i)
	         ktype = iType(k)

               
	         if (ktype.EQ.0) then          ! regular plasma cell

c/    Find corresponding column for this flux:


		    do ll = 1, nSides(k)
		       if (adjCell(ll,k).EQ.i) then
                          ll0=ll
                       endif
                    enddo

                    do na=1,mExp
                       i_col=npos(k,ll0,na)

	            uterm = transm(kk,jj,na,nb,i)
                    nb0=nab(na,nb)

                    if(nb.le.(1+inon))then
                      if(na.le.(1+inon))then
                        na0=1
                      else
                        na0=na
                      endif

		      cterm=colterm(kk,na0,i)*A_cxk(kk,i)* 
     .		        pEscpk(kk,i)*lambdak(kk,i,jj,nb0)
                    else
                       cterm=0.0
                    endif             ! end  if nb==1

                  if(abs(uterm+cterm).gt.0)then
	            i_elmnt = i_elmnt + 1

	            a_sparse(i_elmnt) = - (uterm + cterm)

	            i_sparse(i_elmnt) = i_row
	            i_sparsec(i_elmnt) = i_col
	            endif

                    enddo                      !End of the na loop

	         else if (ktype.EQ.1) then          ! core plasma cell

c/    Make sure we don't double-count the contribution from "jj":

                   do na=1,1+inon
	          
		    if (kk.NE.jj.or.(na.ne.nb))then
		     
		     i_col = npos(i,kk,na)
		     k_rel = k - nCells
		     uterm_p =  transm(kk,jj,na,nb,i)
                     nb0=nab(na,nb)

                     if(na.eq.1)then
                       iflag=1
                     else
                       iflag=-1
                     endif


	            cterm_p =colterm(kk,1,i)*A_cxk(kk,i)* 
     .	             pEscpk(kk,i)*lambdak(kk,i,jj,nb0)
          

	             if(abs(uterm+cterm).gt.0)then

	             i_elmnt = i_elmnt + 1

	             a_sparse(i_elmnt)=-albedo(k_rel)*(uterm_p+cterm_p)
     .                     *iflag

	             i_sparse(i_elmnt) = i_row
	             i_sparsec(i_elmnt) = i_col
	             endif

		   endif                          ! kk == jj if loop
                  enddo                           ! end of do loop over na

	       else if (ktype.EQ.2) then          ! Wall segment

		  kw = k - (nCells + nPlasmReg)

	      if((irefl.EQ.0).OR.(irefl.EQ.1.AND.(zwall(kw).LE.0.0)))then
c/	                                  ! mirror or vacauum boundary
              do na=1, mExp
	           if(kk.ne.jj.or.(na.ne.nb))then
	           i_col=npos(i,kk,na)
                   nb0=nab(na,nb)
	           uterm_wf=refln(kw)*transm(kk,jj,na,nb,i)
	           if(na.eq.1)then
	             uterm_w0=transw_0(jj,1,nb,kw)
	           else
	             uterm_w0=0
	           endif
                   if(na.le.(1+inon))then
                     na0=1
                   else
                     na0=na
                   endif
                 if(nb.le.(1+inon))then

		    cterm_wf =stransw_f(kw,na0)*refln(kw)
     .                *A_cxwf(kw)*pEscpwf(kw)*lambdawf(kw,jj,nb0)
		    cterm_w0 = stransw_0(kw,na) *
     .		        A_cxw0(kw)*pEscpw0(kw)*lambdaw0(kw,jj,nb0)
                  else
                     cterm_wf=0.0
	             cterm_w0=0.0
                  endif             ! end  if nb==1
	            if(na.eq.2)then   ! odd moments
	              uterm_wf=-uterm_wf
	              cterm_wf=-cterm_wf
	            endif
                    if(inon.eq.1.and.(na.eq.3))then
                      uterm_wf=-uterm_wf
                      cterm_wf=-cterm_wf
                    endif

	            if(abs(uterm_wf+cterm_wf).gt.0)then
	              i_elmnt = i_elmnt + 1

	              a_sparse(i_elmnt) = - (uterm_wf + cterm_wf)

	              i_sparse(i_elmnt) = i_row
	              i_sparsec(i_elmnt) = i_col

	            endif                !end of uterm+cterm.ne.0
	            if(na.eq.1)then
	            b(i_row) = b(i_row) + 
     .		        g_ex(kw) * (uterm_w0 + cterm_w0) +
     .		        g_ion(kw) * (uterm_wf  + 
     .                  cterm_wf )
	            endif
                  


	           endif                !end of kk.ne.jj.and.na.ne.nb

	        enddo                   !end of na loop


	      else            !  a martial reflection wall segment


c/    Make sure we don't double-count the contribution from "jj":
	         
                  do na=1,1+inon
                     nb0=nab(na,nb) 
		     if(kk.NE.jj.or.(na.ne.nb))then
		     
		     i_col = npos(i,kk,na)
                     if(na.eq.2)then
                       iflag=-1
                     else
                       iflag=1
                     endif

c/    Uncollided contributions first:

		     uterm_w0=transw_0(jj,na,nb,kw)
		     uterm_wf=refln(kw)*transw_f(jj,na,nb,kw)
		     uterm_ws=(1.-refln(kw))*(1.-fwabsorb(kw))*
     .		        transw_s(jj,na,nb,kw)

c/    Collided contributions:

                     if(nb.le.(1+inon))then
                       if(na.eq.1)then
		         cterm_w0 =  stransw_0(kw,na) *
     .		          A_cxw0(kw)*pEscpw0(kw)*lambdaw0(kw,jj,nb0)
                       else
                         cterm_w0=0.0
                       endif

		       cterm_wf =stransw_f(kw,1)*refln(kw) * 
     .		          A_cxwf(kw)*pEscpwf(kw)*lambdawf(i,jj,nb0)

		       cterm_ws=stransw_s(kw,1)*(1.-refln(kw))*
     .		          (1.- fwabsorb(kw))*A_cxws(kw) * 
     .                    pEscpws(kw)*lambdaws(kw,jj,nb0)
                      else
                        cterm_w0=0
                        cterm_wf=0
                        cterm_ws=0
                      endif                        !end if nb==1

	           if(abs(uterm_wf+uterm_ws+cterm_wf+cterm_ws).gt.0)then

	             i_elmnt = i_elmnt + 1

	             a_sparse(i_elmnt) = - (uterm_wf + uterm_ws +
     .       		        cterm_wf + cterm_ws)

	             i_sparse(i_elmnt) = i_row
	             i_sparsec(i_elmnt) = i_col
	             endif

c/    Assign contributions to RHS vector b:
                     if(na.eq.1)
     .               b(i_row) = b(i_row) + 
     .		        g_ex(kw) * (uterm_w0 + cterm_w0) +
     .		        g_ion(kw) * (uterm_wf + uterm_ws + 
     .                  cterm_wf + cterm_ws)

		  endif                      ! kk == jj if loop
                  enddo                   !End of do loop na

	      endif                           ! irefl==1

               endif                              ! End of "ktype" loop

            enddo                                 ! End of "kk" do loop

c/    Add any contributions to RHS from volumetric sources:

            if(nb.eq.1)then

	      b(i_row) = b(i_row) + 
     .	                 S_ext(i) * (pEscp0s(i) * lambda(i,jj) +
     .                   (1.- Pescp0s(i)) * A_cx0(i) * pEscp(i) * 
     .                   lambda(i,jj))
             endif

          enddo  ! End of Do loop over nb=1,mExp

         enddo   ! End of DO loop over "jj" sides of cell -i-

      enddo      ! End of DO loop over cells -i-

      nElmnts=i_elmnt
      do i=1,nElmnts
         i_sparse(nElmnts+i)=i_sparsec(i)
      enddo

      return
      end

      integer function nab(na,nb)
      implicit none
      include 'neutGlob.inc'
      integer na, nb 
      if(inon.eq.0)then
        if(nb.eq.1)then
          nab=1
        else
          nab=5
        endif
      else if(nd0.eq.0)then
        if(nb.eq.1)then
          nab=na
        else
          nab=5
        endif
      else if(nd0.gt.0)then
        if(na.le.(1+inon))then
          if(nb.le.(1+inon))then
            nab=(na-1)*2+nb
          else
            nab=5
          endif
        else
          if(nb.eq.1)then
             nab=1
          else
             nab=5
          endif
        endif
      endif
      return
      end
        
