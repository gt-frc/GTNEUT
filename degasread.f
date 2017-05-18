      subroutine degasread

c/    This subroutine reads the DEGAS atomic physics data files

      implicit none

      include 'degdat.inc'
      include 'comiou.inc'
      character*80 zdummy
      integer jd, jt, je
      real tlogfact

c/    Create index arrays for interpolation:

      dhkpt(1) = 10.0
      do jd = 2, mpdh
	 dhkpt(jd) = dhkpt(jd-1) + 0.5
      enddo
      
      tlogfact = log(10.0) / 10.0
      do jt = 1, mpeh
	 ehkpt(jt) = real(jt-13)*tlogfact
      enddo

c/    Same for the ion cross sections:
      do je = 1, mpe
	 ekpt(je) = real(je-1)*tlogfact
      enddo
      erefmin = exp(ekpt(1))
      erefmax = exp(ekpt(mpe))
      ekptmpe1 = real(mpe-1) / ekpt(mpe)

c/    Read the contents of the file ehr1.dat (electron impact ionization)

      read(ioeh, 9013) zdummy

      do jd=1,mpdh
         read(ioeh,9013) zdummy
         read(ioeh,9012) (wsveh(jt,jd),jt=1,mpeh)
      enddo

      read(ioeh,9013) zdummy

      do jd=1,mpdh
         read(ioeh,9013) zdummy
         read(ioeh,9012) (wsveh0(jt,jd),jt=1,mpeh)
      enddo

      read(ioeh,9013) zdummy

      do jd=1,mpdh
         read(ioeh,9013) zdummy
         read(ioeh,9012) (welms1(jt,jd),jt=1,mpeh)
      enddo

      read(ioeh,9013) zdummy

      do jd=1,mpdh
         read(ioeh,9013) zdummy
         read(ioeh,9012) (welms2(jt,jd),jt=1,mpeh)
      enddo

      read(ioeh,9013) zdummy

      do jd=1,mpdh
         read(ioeh,9013) zdummy
         read(ioeh,9012) (pne31(jt,jd),jt=1,mpeh)
      enddo

      read(ioeh,9013) zdummy

      do jd=1,mpdh
         read(ioeh,9013) zdummy
         read(ioeh,9012) (pne32(jt,jd),jt=1,mpeh)
      enddo

      read(ioeh,9013) zdummy

      do jd=1,mpdh
         read(ioeh,9013) zdummy
         read(ioeh,9012) (pne21(jt,jd),jt=1,mpeh)
      enddo

      read(ioeh,9013) zdummy

      do jd=1,mpdh
         read(ioeh,9013) zdummy
         read(ioeh,9012) (pne22(jt,jd),jt=1,mpeh)
      enddo

c/    Read the contents of the file cxionh.dat (CX and ion impact ioniz.)
c/    (The first index corresponds to the temperature of the background
c/    ions and the second to the energy of the neutrals)

      read(iocx,9015) ((svphe(jt,je),jt=1,mpe),je=1,mpe)
      read(iocx,9015) ((svphcx(jt,je),jt=1,mpe),je=1,mpe)



9012  format(10(6(1x,e12.5)/))
9013  format(a80)
9015  format(/,8(6(1x,e12.5)/))

      return
      end
