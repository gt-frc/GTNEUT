      subroutine solvers

c/    This subroutine solves the system of TEP equations using the
c/    sparse matrix unsymmetric-pattern multifrontal package UMFPACK
c/    version 2.2.1. This package has been developed by Tim Davis at
c/    the University of Florida and is equivalent to the MA38 routine
c/    from the Harwell Subroutine Library (HSL). More information about
c/    the package and the method can be found at:
c/
c/    http://www.cise.ufl.edu/research/sparse/umfpack/
c/
c/    The UMFPACK package is supposed to be used for educational, research
c/    and benchmarking purposes by non-profit organizations and the U.S.
c/    government. Please, read the README file included in the UMFAPCK 2.2.1
c/    distribution for restrictions on its use.
c/
c/    This solver subroutine uses the single-precision version of UMFPACK
c/    but extension to double precision is trivial.
c/
c/    06/10/03, Created by John Mandrekas for GTNEUT
c/    04/20/17, Modified to use the MA38 libraries due to the unavailability of the necessary version of UMFPACK
c/              MA38 subroutines are drop in replacements. No other changes to the code were made. - Max Hill

      implicit none
      include 'neutGlob.inc'
      include 'comiou.inc'

c/    UMFPACK declarations:
c/    --------------------
      integer lvalue, lindex
c/    ---------------------------------------------
c/    original value parameter (lvalue = 5000000, lindex = 5000000) modified by zwf 1/09
c/    ---------------------------------------------
      parameter (lvalue = 9500000, lindex = 9500000)

c/    Variables for UMS21I and following routines:
c/    -------------------------------------------
      integer keep(20), icntl(20)
      real cntl(10)

c/    Variables for UMS2FA and following routines:
c/    --------------------------------------------
      integer job_fa, index_umf(lindex), info_umf(40)
      real value_umf(lvalue), rinfo(40)
      logical transa

c/    Variables for the UMS2SO routine:
c/    --------------------------------
      integer job_so
      real x_vec(maxEqs), w_umf(4*maxEqs)
      logical transc

c/    Local variables:
c/    ---------------
      integer i, jj, j, jpos, kk, k, kw, ll, lk, index, ipl, i_true,
     .   iw, iw_true, ktype, k_indx, na, iflag,mExp0,nab,nb0,na0,irep,
     .   in
      real denom, uflux, cflux_n, cterm, g_net, cplfactr,
     .     inzrate_n, gterm_k, src_ion_n, src_ion,
     .     neutDens0, tflux, slnon


c/    Assign values to the value and index arrays:

      do i = 1, nElmnts
         index_umf(i) = i_sparse(i)
         index_umf(nElmnts+i) = i_sparse(nElmnts+i)
         value_umf(i) = a_sparse(i)
      enddo

c/    First step is to call the initialization routine UMS21I. This
c/    routine initializes user-controllable parameters to default values,
c/    and non-user controllable parameters:

c/    4/20/2017 - Replacing UMS21I with the equivelant MA38I - Maxwell Hill

      call MA38I (keep, cntl, icntl)

c/    Notice that there are options for long diagnostic output. This is
c/    controlled by the value of ICNTL(3). The default (accepted here)
c/    is icntl(3) = 2, which prints error messages and terse diagnostics
c/    in the file with unit equal to icntl(2).

      icntl(2) = nsdbug
      icntl(3) = 2

c/    We next have to call the general analysis and factorization
c/    routine UMS2FA, which finds a sparsity-preserving and numerically
c/    acceptable pivot order and computes the LU factors, PAQ = LU.

c/    4/20/2017 - Replacing UMS2FA with the equivelant MA38A - Maxwell Hill

      job_fa = 1          ! preserve input matrix
      transa = .false.    ! Factorize A, not A'

      call MA38A (nEqs, nElmnts, job_fa, transa, lvalue, lindex,
     .   value_umf, index_umf, keep, cntl, icntl, info_umf, rinfo)

c/    We now check info(1) to make sure that it is not negative,
c/    indicating that no errors were found. If it's not equal to
c/    zero, we stop the code:

      if (info_umf(1).LT.0) then
         write (6,1000) info_umf(1)
 1000    format (1x, 'Error in the UMS2FA routine: info(1) = ', i2)
         write (nsdbug, 1050) info_umf(1)
 1050 format (1x, 'Value of info(1) from UMS2FA routine = ',
     .   i2)
         return
      endif

c/    Finally, we are ready to call the routine that actually solves the
c/    A*x = b sparse linear system of equations:

c/    Notice that iterative refinement is only performed if job_so = 0,
c/    job_fa = 1, and icntl(8) > 0.

c/    The following 2 flags ensure that we solve A*x = b:

c/    4/20/2017 - Replacing UMS2SO with the equivalent MA38C - Maxwell Hill

      transc = .false.
      job_so = 0
      icntl(8) = isparsitr    ! Steps of iterative refinement.

      call MA38C (nEqs, job_so, transc, lvalue, lindex, value_umf,
     .     index_umf, keep, b, x_vec, w_umf, cntl, icntl, info_umf,
     .     rinfo)

c/    Check the value of info_umf(1) flag. Write it into the nsdbug file.
c/    Terminate routine if negative:

      write (nsdbug, 1100) info_umf(1)
 1100 format (1x, 'Value of info(1) from UMS2SO routine = ',
     .   i2)

      if (info_umf(1).LT.0) then
         write (6,1200) info_umf(1)
 1200    format (1x, 'Error in the UMS2SO routine: info(1) = ', i2)
         return
      endif

c/    We now accept the solution, and calculate various parameters from
c/    the solution vector x_vec:

c/    Calculate internal fluxes from solution vector:

      do i = 1, nCells
        do jj = 1, nSides(i)
           do na=1,mExp
              index = npos(i,jj,na)
              gflux(i,jj,na) = x_vec(index)
           enddo
        enddo
      enddo

c/    Now, calculate fluxes from plasma and wall segments, using albedo
c/    and reflection conditions:

c/    Fluxes from Core plasma (albedo BC):

      if (nPlasmReg.NE.0) then
         do ipl = 1, nPlasmReg
            i_true = nCells + ipl
            do jj = 1, nSides(i_true)
               j = adjCell(jj,i_true)
               do kk = 1, nSides(j)
                  if (adjCell(kk,j).EQ.i_true) then
                     do na=1,mExp
                        if(na.eq.(1+inon))then
                          gflux(i_true,jj,na)=albedo(ipl)*gflux(j,kk,na)
                          if(na.eq.2)then
                            gflux(i_true,jj,na)=-gflux(i_true,jj,na)
                          endif
                        else
                          gflux(i_true,jj,na)=0
                        endif
                     enddo
                  endif
               enddo
            enddo
         enddo
      endif

c/    Fluxes from Wall Segments:

      do iw = 1, nWallSegm
         iw_true = nCells + nPlasmReg + iw
         j = adjCell(1, iw_true)
         do ll = 1, nSides(j)
            if (adjCell(ll,j).EQ.iw_true) then
c            write(*,*)g_ex(iw),gflux(j,ll),refln(iw),fwabsorb(iw)
               gflux (iw_true,1,1) = g_ex(iw) +
     .        (gflux(j,ll,1) + g_ion(iw)) *
     .        (refln(iw) + (1.- refln(iw)) * (1.- fwabsorb(iw)))
             endif
         enddo
      enddo

c/    Now we calculate the uncollided and collided fluxes (including
c/    the contribution of the first collision source) plus the ionization
c/    rate in each cell, from which we can determine the neutral density:
c/    Calculate the sum of temperature for

      do i = 1, nCells
         do jj = 1, nSides(i)
         if(uflux.lt.0.0)then
	     print*, 'uflux = ',uflux
            print*, 'i     = ',i
            print*, 'jj    = ',jj
         endif
            uflux = 0.0                     ! Uncollided
            do kk = 1, nSides(i)
               k = adjCell(kk,i)
               ktype = iType(k)
               if (ktype.LT.2) then         ! Not a wall
                  do lk = 1, nSides(k)
                     if (adjCell(lk,k).EQ.i) then
                        do na=1,mExp
                           uflux=uflux+gflux(k,lk,na)*
     .                           transm(kk,jj,na,1,i)
                        enddo
                     endif
                  enddo
               else if (ktype.EQ.2) then     ! Wall segment
                  kw = k - (nCells + nPlasmReg)
                  uflux = uflux+g_ex(kw)*transw_0(jj,1,1,kw) +
     .                    (gflux(i,kk,1) + g_ion(kw)) *
     .                    (refln(kw) * transw_f(jj,1,1,kw) +
     .                    (1.- refln(kw)) * (1.- fwabsorb(kw)) *
     .                       transw_s(jj,1,1,kw))
               endif
            enddo

	if(uflux.eq.0)then
                   print*, 'i = ',i,' j = ',jj
	            print*, 'uflux      = ',uflux
                   print*, 'ktype =',ktype
                   print*, 'transm(kk,jj,na,1,i) =',transm(kk,jj,na,1,i)
                   print*, 'gflux(,k,lk,na) = ',gflux(k,lk,na)

c/		     pause
	endif

            gflux_u(i,jj) = uflux           ! Total uncollided flux

c/    Collided fluxes and ionization rates:

            cflux_n = 0.0                   ! Rest of collided
            inzrate_n = 0.0

            do kk = 1, nSides(i)
               k = adjCell(kk,i)
               ktype = iType(k)
               if (ktype.LT.2)then
c			          ! Not a wall
                 do na=1,mExp
                    nb0=nab(na,1)
                    if(na.le.(1+inon))then
                      cterm=colterm(kk,1,i)
                    else
                      cterm =colterm(kk,na,i)
                    endif
                    do lk = 1, nSides(k)
                       if (adjCell(lk,k).EQ.i) then
                          k_indx = lk
                       endif
                    enddo

                    gterm_k = gflux(k, k_indx,na) * cterm


                    cflux_n = cflux_n + gterm_k *  pEscpk(kk,i)*
     .                 A_cxk(kk,i)*lambdak(kk,i,jj,nb0)

                    if(na.gt.1.and.na.le.(1+inon))then
                      slnon=0.0
                      do in=1,nSides(i)
                         slnon=slnon+lambdak(kk,i,in,nb0)
                      enddo
                      inzrate_n=inzrate_n-gterm_k*slnon*A_cxk(kk,i)
                    else
                      inzrate_n = inzrate_n + gterm_k *
     .                        (1.- pEscpk(kk,i)*A_cxk(kk,i))
                    endif
                 enddo


               else if (ktype.EQ.2) then         ! Wall segment

                  kw = k - (nCells + nPlasmReg)
                  nb0=nab(1,1)


                  cflux_n = cflux_n + g_ex(kw)*stransw_0(kw,1)*
     .               A_cxw0(kw)* pEscpw0(kw)*lambdaw0(kw,jj,nb0)+
     .               (gflux(i,kk,1)+g_ion(kw))*(stransw_f(kw,1)*
     .               refln(kw)*A_cxwf(kw)*pEscpwf(kw)*
     .               lambdawf(kw,jj,nb0)+stransw_s(kw,1)*(1-refln(kw))*
     .               (1.- fwabsorb(kw))*A_cxws(kw)*pEscpws(kw)
     .               *lambdaws(kw,jj,nb0))

                  inzrate_n = inzrate_n + (g_ex(kw) *
     .               stransw_0(kw,1)*(1-A_cxw0(kw)*pEscpw0(kw))
     .               +(gflux(i,kk,1)+g_ion(kw)) *
     .               (refln(kw)*stransw_f(kw,1) *
     .               (1.- A_cxwf(kw)*pEscpwf(kw))+(1.- refln(kw)) *
     .               (1.- fwabsorb(kw))*stransw_s(kw,1)*
     .               (1.-A_cxws(kw)*pEscpws(kw))))



	           if(((irefl.eq.0.or.(zwall(kw).le.0)).and.
     .         		   Rwall(kw).gt.0).or.inon.eq.1)then
                         if(irefl.eq.1.and.(zwall(kw).gt.0))then
                           mExp0=1+inon
                         else
                           mExp0=mExp
                         endif
		         do na=2, mExp0
                          nb0=nab(na,1)
			  if(na.eq.2)then
	                    iflag=-1
                            na0=na
                            if(inon.eq.1)na0=1
		          else if(inon.eq.1.and.(na.eq.3))then
                            iflag=-1
                            na0=na
                          else
			    iflag=1
                            na0=na
			  endif

                      cterm =colterm(kk,na0,i)
                      do lk = 1, nSides(k)
                         if (adjCell(lk,k).EQ.i) then
                            k_indx = lk
                         endif
                      enddo

                      gterm_k = gflux(i,kk,na)*cterm
     .                         *iflag*refln(kw)


                      cflux_n = cflux_n+gterm_k*
     .                   A_cxwf(kw)*pEscpwf(kw)*lambdawf(kw,jj,nb0)

                    if(na.gt.1.and.na.le.(1+inon))then
                      slnon=0.0
                      do in=1,nSides(i)
                         slnon=slnon+lambdawf(kw,in,nb0)
                      enddo
                      inzrate_n=inzrate_n-gterm_k*slnon*A_cxwf(kw)
                    else
                      inzrate_n = inzrate_n + gterm_k *
     .                        (1.- A_cxwf(kw)*pEscpwf(kw))
                    endif



	             enddo     !end of do loop over na

	           endif       !end of IF loop over irefl

               endif         ! end of IF loop over ktype

            enddo            ! end of DO loop over kk side summation

            gflux_cn(i,jj) = cflux_n

            gflux_c(i,jj) =  cflux_n

         enddo               ! end of DO loop over jj side summation

         inzRate(i) = inzrate_n

c/    External source contributions, if any:

         if (S_ext(i).GT.0.0) then

            src_ion_n = (1.- A_cx0(i)*pEscp(i))
            src_ion = S_ext(i) * (1. - pEscp0s(i)) * src_ion_n

            inzRate(i) = inzRate(i) + src_ion

c/    Must also include contribution to collided fluxes:

            do jj = 1, nSides(i)

               gflux_cn(i,jj)=gflux_cn(i,jj)+S_ext(i)*(pEscp0s(i)
     .            +A_cx0(i)*(1.- pEscp0s(i))*pEscp(i))*lambda(i,jj)

               gflux_c(i,jj) = gflux_cn(i,jj)

            enddo

         endif






c/    Calculate neutral density from ionization rate, assuming that
c/    svion_tot does not depend strongly on the neutral energy:

         neutDens0 = srcNormFact * inzRate(i) /(area(i) * ionDens(i)
     .      * svion_tot(i))
	   neutDens(i)=neutDens0

      enddo                  ! end of DO loop over internal cells


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c/    Now we calculate the average temperature of neutrals from region i
c/    side jj

      do irep=1,100

      do i = 1, nCells
         do jj = 1, nSides(i)

            tflux = 0.0
c/		                                Uncollided
            do kk = 1, nSides(i)
               k = adjCell(kk,i)
               ktype = iType(k)
               if (ktype.LT.2) then         ! Not a wall
                  do lk = 1, nSides(k)
                     if (adjCell(lk,k).EQ.i) then
                        do na=1,mExp
                           tflux=tflux+gflux(k,lk,na)*
     .                    transm(kk,jj,na,1,i)/sqrt(tneut(k,lk))
                        enddo
                     endif
                  enddo
               else if (ktype.EQ.2) then     ! Wall segment
                  kw = k - (nCells + nPlasmReg)
                  tflux=tflux+g_ex(kw)*transw_0(jj,1,1,kw)/sqrt(eneut)+
     .                (gflux(i,kk,1)+g_ion(kw))*refln(kw)
     .                /sqrt(tneut(i,kk)/refln(kw)
     .                 *refle(kw)) * transw_f(jj,1,1,kw) +
     .                 (gflux(i,kk,1) + g_ion(kw))/sqrt(twall(kw))
     .                 *(1.- refln(kw)) * (1.- fwabsorb(kw)) *
     .                     transw_s(jj,1,1,kw)
               endif
            enddo


c/    Collided fluxes and ionization rates:


            do kk = 1, nSides(i)
               k = adjCell(kk,i)
               ktype = iType(k)
               if (ktype.LT.2)then
c			          ! Not a wall
                 do na=1,mExp
                    nb0=nab(na,1)
                    if(na.le.(1+inon))then
                      cterm=colterm(kk,1,i)
                    else
                      cterm =colterm(kk,na,i)
                    endif
                    do lk = 1, nSides(k)
                       if (adjCell(lk,k).EQ.i) then
                          k_indx = lk
                       endif
                    enddo

                    gterm_k = gflux(k, k_indx,na) * cterm


                    tflux=tflux+gterm_k*pEscpk(kk,i)/sqrt(ionTemp(i))*
     .                 A_cxk(kk,i)*lambdak(kk,i,jj,nb0)

                 enddo


               else if (ktype.EQ.2) then         ! Wall segment

                  kw = k - (nCells + nPlasmReg)
                  nb0=nab(1,1)


                  tflux = tflux + (g_ex(kw)*stransw_0(kw,1)*
     .               A_cxw0(kw)* pEscpw0(kw)*lambdaw0(kw,jj,nb0)+
     .               (gflux(i,kk,1)+g_ion(kw))*(stransw_f(kw,1)*
     .               refln(kw)*A_cxwf(kw)*pEscpwf(kw)*
     .               lambdawf(kw,jj,nb0)+stransw_s(kw,1)*(1-refln(kw))*
     .               (1.- fwabsorb(kw))*A_cxws(kw)*pEscpws(kw)
     .               *lambdaws(kw,jj,nb0)))/sqrt(ionTemp(i))


	           if(((irefl.eq.0.or.(zwall(kw).le.0)).and.
     .         		   Rwall(kw).gt.0).or.inon.eq.1)then
                         if(irefl.eq.1.and.(zwall(kw).gt.0))then
                           mExp0=1+inon
                         else
                           mExp0=mExp
                         endif
		         do na=2, mExp0
                          nb0=nab(na,1)
			  if(na.eq.2)then
	                    iflag=-1
                            na0=na
                            if(inon.eq.1)na0=1
		          else if(inon.eq.1.and.(na.eq.3))then
                            iflag=-1
                            na0=na
                          else
			    iflag=1
                            na0=na
			  endif

                      cterm =colterm(kk,na0,i)
                      do lk = 1, nSides(k)
                         if (adjCell(lk,k).EQ.i) then
                            k_indx = lk
                         endif
                      enddo

                      gterm_k = gflux(i,kk,na)*cterm
     .                         *iflag*refln(kw)


                      tflux = tflux+gterm_k/sqrt(ionTemp(i))
     .                   *A_cxwf(kw)*pEscpwf(kw)*lambdawf(kw,jj,nb0)



	             enddo     !end of do loop over na

	           endif       !end of IF loop over irefl

               endif         ! end of IF loop over ktype

            enddo            ! end of DO loop over kk side summation
            if (S_ext(i).GT.0.0) then

                tflux=tflux+S_ext(i)*(pEscp0s(i)*sqrt(eneut)+
     .           1.0/sqrt(ionTemp(i))*
     .           A_cx0(i)*(1.- pEscp0s(i))*pEscp(i))*lambda(i,jj)


            endif

            tneut(i,jj)=((gflux_u(i,jj)+gflux_c(i,jj))/tflux)**2



         enddo               ! end of DO loop over jj side summation


      enddo                  ! end of DO loop over internal cells
      enddo                  ! endo of DO loop over irep

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c/    Calculate neutral density using alternative approach for check:

      do i = 1, nCells
        g_net = 0.0
        do kk = 1, nSides(i)
          k = adjCell(kk,i)
          do ll = 1, nSides(k)
             if (adjCell(ll,k).EQ.i) jpos = ll
          enddo
          g_net = g_net + gflux(i,kk,1) - gflux(k,jpos,1)
        enddo
       neutDens0 = srcNormFact * (S_ext(i) - g_net) /
     .   (area(i) * ionDens(i) * svion_tot(i))
c       neutDens(i)=neutDens0
      enddo

      return
      end
