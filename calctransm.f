      subroutine calcTransm

c/    This subroutine sets up the calculation for the first-flight
c/    transmission coefficients for each internal cell and wall 
c/    segment. Unlike the theory, where T_kji is defined as the 
c/    transmission coefficient FROM cell -k- to cell -j- through 
c/    cell -i-, in the code the transmission coefficients are defined
c/    in terms of the sides of cell -i-, i.e. T_kji is the transmission
c/    coefficient FROM side -k- TO side -j- of cell -i-. This conserves
c/    memory and makes the code easier to understand. 
c/
c/    For internal and core plasma cells, this routine calculates:
c/    
c/    transm(k,j,i) : Transmission coefficient from side -k- to side -j-
c/                    of cell -i-
c/
c/    If the -FROM- side -k- is a wall segment, then three transmission 
c/    coefficients are computed:

c/    transw_f(j,kw)  : transmission coeff. for fast neutrals
c/    transw_s(j,kw)  : transmission coeff. for slow neutrals
c/    transw_0(j,kw)  : transmission coeff. for source neutrals
c/
c/    Notice that if sides k and j are identical, transm is set to 0.
c/
c/    In addition, the sum of the transmission coefficients from a side
c/    to all other sides is also computed, since it appears frequently
c/    in the balance equations:

c/    stransm(k,i)  = Sum_over_all_j (transm(k,j,i))
c/    stransw_f(kw) = Sum_over_all_j (transw_f(j,kw))
c/    stransw_s(kw) = Sum_over_all_j (transw_s(j,kw))
c/    stransw_0(kw) = Sum_over_all_j (transw_0(j,kw))

      implicit none
      include 'neutGlob.inc'

c/    Local variables:
c/
      integer i, k, kk, kw, ll, jj, na
      real   sumtw_f, sumtw_s, sumtw_0, sumt


      do i=1, nCells
         do kk = 1, nSides(i)
            do jj = 1, nSides(i)
               call TransmCoeff(i, kk, jj)
            enddo
         enddo
      enddo

c/    Calculate the summations over all sides:

      do i=1, nCells
        do kk=1,nSides(i)
           do na=1, mExp
              sumt=0.0
              do jj=1,nSides(i)
                 sumt=sumt+transm(kk,jj,na,1,i)
              enddo
              if(na.eq.1)then
                if(sumt.gt.1)sumt=1.0
                sumt=1.0-sumt
              else
                sumt=-sumt
              endif
              colterm(kk,na,i)=sumt
           enddo
        enddo
      enddo
      
      do i = 1, nCells
	 do kk = 1, nSides(i)
	    k = adjCell(kk,i)
            if (iType(k).EQ.2) then
	       kw = k - (nCells + nPlasmReg)
               do na=1,mExp
	       sumtw_f = 0.0
	       sumtw_s = 0.0
	       sumtw_0 = 0.0
	       do ll = 1, nSides(i)
		  sumtw_f = sumtw_f + transw_f(ll,na,1,kw)
		  sumtw_s = sumtw_s + transw_s(ll,na,1,kw)
		  sumtw_0 = sumtw_0 + transw_0(ll,na,1,kw)
               enddo
               if(na.eq.1)then
	         if (sumtw_f.GT.1.0) sumtw_f = 1.0
	         if (sumtw_s.GT.1.0) sumtw_s = 1.0
	         if (sumtw_0.GT.1.0) sumtw_0 = 1.0
                 sumtw_f=1-sumtw_f
                 sumtw_s=1-sumtw_s
                 sumtw_0=1-sumtw_0
               else
                 sumtw_f=-sumtw_f
                 sumtw_s=-sumtw_s
                 sumtw_0=-sumtw_0
               endif
	       stransw_f(kw,na) = sumtw_f
	       stransw_s(kw,na) = sumtw_s
	       stransw_0(kw,na) = sumtw_0
               enddo
            endif
	 enddo                                     ! DO loop over sides
      enddo                                        ! DO loop over cells

      return
      end
