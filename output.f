      subroutine output(iflag)

c/    This routine writes the output for the neutrals code.
c/    Written by John Mandrekas, GIT, 09/26/93

c/    iflag  : If equal to 0, write debugging information for 
c/             the i_inp = 1 case (rectangular regions)
c/           : If equal to 1, normal output.

      implicit none
      include 'neutGlob.inc'
      include 'comiou.inc'
      include 'consts.inc'

c/    Define local variables:
      integer i, iw, j, jj, k, kk, kw, l, ll, iflag, na,nb,icn,nduedge
      real inzDens, delta_x, delta_y, averden,iRatezwf




c/    Write debugging information for rectangular regions (automatic
c/    input generation) if iflag = 0:
      
      if (iflag.EQ.0) then
         delta_x = lside(1,1)
         delta_y = lside(2,1)
         write (ndbug, 500) nCells, nWallSegm, delta_x, delta_y
         write(ndbug, '(1x)')
         do i = 1, nCells                                      
            do j = 1, 4
               write (ndbug, 550) j, i, adjCell(j,i)
            enddo                                              
         enddo                                                 
         do i = nCells+1, nCells + nWallSegm                               
            write (ndbug, 600) i, adjCell(1,i)                    
         enddo                                                 
         write(ndbug, '(1x)')
         return
      endif                       ! End of IF loop over rectangular input
                                                      
c/    If the first element of the array prntOrdr is negative, then assign
c/    natural order:

      if(prntOrdr(1).lt.0) then
        do i = 1, nCells
          prntOrdr(i) = i
        enddo
      endif

c/    Write debugging output to file "neut.dbg" (unit ndbug) if idbug=1:
      
      if (idbug.eq.1) then
         do i = 1, nCells
           write (ndbug, '(1x, i5)') i
           write (ndbug, 1000) svion_e(i), svion_i(i), svion_tot(i),
     .        sv_cx(i), lmfp(i), A_cx(i), pEscp(i), area(i), perim(i)
           do ll = 1, nSides(i)
             l = adjCell(ll,i)
c             write (ndbug, 1050) i, l, lambda(i,ll)
             do kk=1,nSides(i)
                write(ndbug,1055)ll,i,kk,lambdak(ll,i,kk,1)
             enddo
           enddo
           if (i_e0.NE.1) then
              do ll = 1, nSides(i)
                l = adjCell(ll,i)
                if (iType(l).NE.2) then
                   write (ndbug, 1060) l, i, mfp(ll,i)
                   write (ndbug, 1070) l, i, A_cxk(ll,i)
                else if (iType(l).EQ.2) then
                   kw = l - (nCells+nPlasmReg)
                   write (ndbug, 1061) l, i, mfp_wf(kw)
                   write (ndbug, 1071) l, i, A_cxwf(kw)
                   write (ndbug, 1062) l, i, mfp_ws(kw)
                   write (ndbug, 1072) l, i, A_cxws(kw)
                   if (g_ex(kw).GT.0.0) then
                      write (ndbug,1063) i, i, mfp_w0(kw)
                      write (ndbug,1073) l, i, A_cxw0(kw)
                   endif
                endif
              enddo
           endif
           do kk = 1, nSides(i)
              k = adjCell(kk,i)
              if (iType(k).NE.2) then
                 do jj = 1, nSides(i)
                   j = adjCell(jj,i)
                   if (j.NE.k) then
                     do na=1,mExp
                        do nb=1,mExp
                           write(ndbug,1100)k,j,i,transm(kk,jj,na,nb,i)
                        enddo
                     enddo
                   endif
                 enddo
              else if (iType(k).EQ.2) then
                 kw = k - (nCells+nPlasmReg)
                 do jj = 1, nSides(i)
                    j = adjCell(jj,i)
                    if (j.NE.k) then
                      do na=1,mExp
                         write(ndbug,1110)k,j,i,transw_f(jj,1,na,kw)
                      enddo
                      do na=1,mExp
                         write(ndbug,1111)k,j,i,transw_s(jj,1,na,kw)
                      enddo
                      if (g_ex(kw).GT.0.0) then
                        do na=1,mExp
                           write(ndbug,1112)k,j,i,transw_0(jj,1,na,kw)
                        enddo
                      endif
                    endif
                 enddo
              endif
           enddo
           write (ndbug, '(1x)')
         enddo                       ! End of DO loop over nCells

         do i = 1, nPlasmReg
            j = nCells + i
            write (ndbug, 1200) j, albedo(i)
         enddo

c/    If reflection model is on, print out particle and energy 
c/    reflection coefficients:

         if (irefl.NE.0) then
            write (ndbug, '(1x)')
            write (ndbug, 1250)
            do iw = 1, nWallSegm
               write (ndbug, 1255) iw, refln(iw), refle(iw), Rwall(iw)
            enddo
         endif

c/    Details of the Solution Matrix:

         write (ndbug, '(1x)')
         write (ndbug, '(1x, A27)') 'SOLUTION MATRIX INFORMATION'
         write (ndbug, '(1x)')
         write (ndbug, '(1x, A15)') 'Position Index:'
         write (ndbug, '(1x)')
         do i = 1, nCells 
            do jj = 1, nSides(i)
              j = adjCell(jj, i)
              do na=1,mExp
                 write(ndbug,1300)i,j,npos(i,jj,na)
              enddo
            enddo
         enddo

c/      Write the non-zero elements of the sparse solution matrix:

         write (ndbug, '(1x)')
         write (ndbug, '(1x, A31)') 'Elements of coefficient matrix:'
         write (ndbug, '(1x)')
         
         write (ndbug, 1400)
         do i = 1, nElmnts
            write (ndbug, 1450) i, a_sparse(i), i_sparse(i), 
     .        i_sparse(nElmnts+i)
         enddo
         write (ndbug, '(1x)')
         do i = 1, nEqs
           write (ndbug, 1500) i, b(i)
         enddo

      endif   ! END of IF loop for DEBUG output

c/    Write regular output to file neut.out (nout):
c/###################DING'S VERSION MISSING srcNormFact##############
c/    NEUTRAL DENSITIES AND IONIZATION RATES:
      write (nout, 2000)
      write (nout, '(1x)')
      do i = 1, nCells
        j = prntOrdr(i)
c/      recompute inzDens
c/      inzDens = srcNormFact * inzRate(j) / area(j)
c/      use method that works
        inzDens =neutDens(j)*svion_tot(j)*ionDens(j)
c/      used definition from Mandrekas Paper
        iRatezwf = srcNormFact * inzRate(j)
        write (nout, 2100) j, i, neutDens(j), inzDens, iRatezwf
        write (tomat,3200) j,neutDens(j)
        write (tomat,3211) j,inzDens
        write (tomat,3212) j,srcNormFact*inzRate(j)
        write (tomat,3213) j,elecTemp(j)
        write (tomat,3214) j,ionTemp(j)
        write (tomat,3215) j,elecDens(j)
        write (tomat,3216) j,ionDens(j) 
        write (tomat, 3219) i, area(i)
        write (tomat, 3221) i, svion_tot(i)
        write (tomat, 3220) i, neutDens(i)*(area(i) * ionDens(i) 
     .      * svion_tot(i))  
        write (tomat, 3235) i, neutDens(i)*(area(i) * ionDens(i) 
     .      * sv_cx(i))  
        write (tomat, 3230) i, lmfp(i)                           
      enddo

c/    Write to terminal:
      do i = 1, nCells
        j = prntOrdr(i)
c/      recompute inzDens
c/      inzDens = srcNormFact * inzRate(j) / area(j)
c/      use method that works
        inzDens =neutDens(j)*svion_tot(j)*ionDens(j)
c/      used definition from Mandrekas Paper
        iRatezwf = srcNormFact * inzRate(j)
        write (6, 2100) j, i, neutDens(j), inzDens, iRatezwf
      enddo
c/###############################################################
c/##################COPIED FROM MANDREKAS VERSION###############
c/zwf    NEUTRAL DENSITIES AND IONIZATION RATES:
c/zwf      write (nout, 2000)
c/zwf      write (nout, '(1x)')
c/zwf      do i = 1, nCells
c/zwf	j = prntOrdr(i)
c/zwf	inzDens = srcNormFact * inzRate(j) / area(j)
c/zwf	write (nout, 2100) j, i, neutDens(j), inzDens, 
c/zwf     .    	srcNormFact * inzRate(j)
c/zwf	write (tomat,3200) j,neutDens(j)
c/zwf       write (tomat,3210) j,srcNormFact * inzRate(j)
c/zwf       write (tomat,3211) j,inzDens 
c/zwf      enddo
c/zwf
c/zwf    Write to terminal:
c/zwf      do i = 1, nCells
c/zwf	j = prntOrdr(i)
c/zwf	inzDens = srcNormFact * inzRate(j) / area(j)
c/zwf	write (6, 2100) j, i, neutDens(j), inzDens, 
c/zwf     .	   srcNormFact * inzRate(j)
c/zwf      enddo
c/#############################################################
      write (nout, '(1x)')
c/    FLUXES (TOTAL, COLLIDED & UNCOLLIDED)
      write (nout, 2200)
      write (nout, '(1x)')
      do i = 1, nTotal 
        do jj = 1, nSides(i)
          j = adjCell(jj, i)
          write (nout, 2300) i, j,srcNormFact*gflux(i,jj,1), 
     .       srcNormFact*gflux_c(i,jj), srcNormFact * gflux_u(i,jj)
          write (tomat,3218) i, j,srcNormFact*gflux(i,jj,1)
        enddo
      enddo


c/    Write power balance parameters:
      write (nout, '(1x)')
      write (nout, 3000)
      write (nout, 3100) fluxn_tot, fluxi_tot, volsrc_tot, 
     .       coreloss_tot, wloss_tot, ionloss_tot, tot_part_in,
     .       tot_part_out, partbalnerr
      write(tomat, 3217) nxleg1,nxleg2,nxcore1,nxcore2,
     .nycore1,nysol1,nxxpt,nxmod,Timeslice,Shotnumber

c/ THIS SECTION WAS PREVENTING IT FROM COMPILING AND 
c/ I DON'T CARE ABOUT UEDGE ANYWAY. COMMENTING OUT.
c/ -MAX HILL 4/20/2017
 
c/      open (nduedge, file = 'loadnd', status = 'unknown')
c/      NY = nycore1 + nysol1 + 1
c/      NX = nxcore1 + nxcore2 + 4*nxxpt + nxleg1 + nxleg2 

c/c/    WRITE ndense to file that loads ndense into UEDGE
c/         write(nduedge,3222) NX,NY
c/         write(nduedge,3226) NX,NY
c/         write(nduedge,3227) NX,NY
c/         write(nduedge,3228) NX,NY
c/         write(nduedge,3232) NX,NY
c/         write(nduedge,3233) NX,NY

c/      do j = 1,NY
c/         do i = 1, NX
c/          icn = i + NX * (j-1)
c/          write(nduedge,3223) i,j,neutDens(icn)
c/	   write(nduedge,3224) i,j,neutDens(icn)*(area(icn) * ionDens(icn) 
c/     .      * svion_tot(icn))
c/	   write(nduedge,3225) i,j,neutDens(icn)*svion_tot(icn)*ionDens(icn)
c/	   write(nduedge,3229) i,j,A_cx(icn)
c/	   write(nduedge,3231) i,j,lmfp(icn)
c/          write(nduedge,3234) i,j,neutDens(i)*(area(i) * ionDens(i) 
c/     .      * sv_cx(i))  


c/         enddo
c/      enddo
c/
c/ END OF MAX EDITS 4/20/2017

 500  format (1x, 'nCells    = ', i4/
     .        1x, 'nWallSegm = ', i4/
     .        1x, 'delta_x   = ', e12.5, ' m',/,
     .        1x, 'delta_y   = ', e12.5, ' m')

 550  format (1x, 'adj(',i1,',',i5,') = ', i5)              
 600  format (1x, 'adj(1,',i5,') = ', i5)                   

 1000 format(1x, '<sv>_e    :', e12.5, ' m^3/s'/
     .       1x, '<sv>_i    :', e12.5, ' m^3/s'/
     .       1x, '<sv>_itot :', e12.5, ' m^3/s'/
     .       1x, '<sv>_cx   :', e12.5, ' m^3/s'/
     .       1x, 'l_mfp     :', e12.5, ' m'/
     .       1x, 'A_cx      :', f8.5, /
     .       1x, 'pEsc      :', e12.5/
     .       1x, 'Area      :', e12.5, ' m^2'/
     .       1x, 'Perimeter :', e12.5, ' m')
 1050 format(1x, 'lambda(',i5, ',', i5,') = ', f8.5)
 1055 format(1x, 'lambdak(',i1,',',i5,',',i1,') =',f8.5)
 1060 format(1x, 'mfp(',i5, ',', i5,')    = ', f8.5, ' m')
 1061 format(1x, 'mfp_wf(',i5, ',', i5,') = ', f8.5, ' m')
 1062 format(1x, 'mfp_ws(',i5, ',', i5,') = ', f8.5, ' m')
 1063 format(1x, 'mfp_w0(',i5, ',', i5,') = ', f8.5, ' m')
 1070 format(1x, 'A_cxk(',i5, ',', i5,')    = ', f8.5, ' m')
 1071 format(1x, 'A_cxwf(',i5, ',', i5,') = ', f8.5, ' m')
 1072 format(1x, 'A_cxws(',i5, ',', i5,') = ', f8.5, ' m')
 1073 format(1x, 'A_cxw0(',i5, ',', i5,') = ', f8.5, ' m')

 1100 format(1x, 'From ', i5, ' to ', i5, ' through ', i5, ' : ', e12.5)
 1110 format(1x, 'From ', i5, ' to ', i5, ' through ', i5, ' : ', e12.5,
     .            ' (fast)')
 1111 format(1x, 'From ', i5, ' to ', i5, ' through ', i5, ' : ', e12.5,
     .            ' (slow)')
 1112 format(1x, 'From ', i5, ' to ', i5, ' through ', i5, ' : ', e12.5,
     .            ' source)')
 1200 format(1x, 'Plasma albedo in region ', i5, ' : ', e12.5)
 1250 format (1x, 'Wall Reflection Coefficients',/
     .        2x, 'iw', 7x, 'RN', 11x, 'RE', 10x, 'Rwall')
 1255 format (1x, i5, 3x, 3(f8.5, 5x))
 1300 format(1x, 'npos(', i5, ',', i5, ') = ', i5)
 1400 format (4x, 'i', 4x, 'a_sparse', 7x, 'i_row', 3x, 'i_col')
 1450 format(1x, i4, 2x, e12.5, 4x, i4, 4x, i4)
 1500 format(1x, 'b(', i5, ') = ', e12.5)

 2000 format(1x, 'i', 3x, 'Region', 3x, 'Neutral Density (#/m3)', 2x, 
     .  'Ionization Density', 2x, 'Ionization Rate (#/s)')
 2100 format(1x, i5, 4x, i5, 8x, e12.5, 12x, e12.5, 12x, e12.5)
 2105 format(1x,e12.5)
 2200 format(17x, 'Total', 7x, 'Collided', 6x, 'Uncollided')
 2300 format(1x, i5, '-->', i5, ':', 3(2x, e12.5))
 3000 format (1x, 'Global Particle Balance')
 3100 format (1x, 'Total external neutral flux = ', e12.5, ' #/s'/
     .        1x, 'Total recycling ion flux    = ', e12.5, ' #/s'/
     .        1x, 'Total volumetric source     = ', e12.5, ' #/s'/
     .        1x, 'Lost to core plasma         = ', e12.5, ' #/s'/
     .        1x, 'Wall absorption or pumping  = ', e12.5, ' #/s'/
     .        1x, 'Lost by ionization          = ', e12.5, ' #/s'/
     .        1x, 'Total particles in          = ', e12.5, ' #/s'/
     .        1x, 'Total particles out         = ', e12.5, ' #/s'/
     .        1x, 'Particle balance error      = ', f8.3, ' %')
 3200 format('Cell.ndense(',i5,') = ',  4x, e12.5,';')     
 3210 format('Cell.gflux(',i5,') = ',  4x, e12.5,';')  
 3211 format('Cell.iDenz(',i5,') = ',  4x, e12.5,';')  
 3212 format('Cell.iRate(',i5,') = ',  4x, e12.5,';')    
 3213 format('Cell.elecTemp(',i5,') = ',  4x, e12.5,';')   
 3214 format('Cell.ionTemp(',i5,') = ',  4x, e12.5,';')        
 3215 format('Cell.elecDens(',i5,') = ',  4x, e12.5,';')        
 3216 format('Cell.ionDens(',i5,') = ',  4x, e12.5,';')       
 3217 format (1x, 'nxleg1    = ', i5/
     .        1x, 'nxleg2    = ', i5/
     .        1x, 'nxcore1   = ', i5/
     .        1x, 'nxcore2   = ', i5/
     .        1x, 'nycore1   = ', i5/
     .        1x, 'nysol1    = ', i5/
     .        1x, 'nxxpt     = ', i5/
     .        1x, 'nxmod     = ', i5/
     .        1x, 'Timeslice = ', i5/
     .        1x, 'Shotnumber= ', i5)
 3218 format('Cell.gflux(',i5,',',i5,') = ', e12.5,';')  
 3219 format('Cell.area(',i5,') = ',  4x, e12.5,';')  
 3220 format('Cell.recomputed_irate(',i5,') = ',  4x, e12.5,';') 
 3221 format('Cell.svion_tot(',i5,') = ',  4x, e12.5,';')   
 3222 format('real gtndense(',i5,',', i5,')')
 3223 format('gtndense(',i5,',', i5,')=',e12.5)
 3224 format('gtirate(',i5,',', i5,')=',e12.5)
 3225 format('gtizdense(',i5,',', i5,')=',e12.5)
 3226 format('real gtirate(',i5,',', i5,')')
 3227 format('real gtizdense(',i5,',', i5,')')
 3228 format('real acx(',i5,',', i5,')')
 3229 format('acx(',i5,',', i5,')=',e12.5)
 3230 format('Cell.mfp(',i5,') = ',  4x, e12.5,';') 
 3231 format('gtmfp(',i5,',', i5,')=',e12.5)
 3232 format('real gtmfp(',i5,',', i5,')')
 3233 format('real crate(',i5,',', i5,')')
 3234 format('crate(',i5,',', i5,')=',e12.5)
 3235 format('Cell.crate(',i5,') = ',  4x, e12.5,';')    
      return
      end
