      subroutine t_ij (x,numfun,funval)

c/    This function returns the integral of the transmission coefficient
c/    over the angle phi for a given value of ksi_i = x (see the paper in 
c/    Nuclear Fusion 34 (1994) 1385)

c/    John Mandrekas, GIT, 12/11/97
c/    04/07/1999, jm: Introduced Bickley function for correct evaluation
c/                    of the transmission coefficient

      implicit none

      include 'locGeom.inc'

      integer maxph, i, numfun, j, nExp,mExp
      parameter (maxph = 51, nExp=3)
      real funval(numfun)
      real x, tvec(maxph,numfun),sina,cosa,phi_min,phi_max,d_phi,phi, 
     .   ksi_j, el_of_phi, pi, tanth, sinth,costh, reslt(numfun)
      real bickley3, bickley2, bickley1
      real ki1,ki2,ki3,ki4,ki5
      real sphi0,cphi0,sphi1,cphi1,x0, phi0,phi1
      real x1, x2

      data pi /3.141592654/

      save pi

c/    Calculate maximum and minimum angles for integration:
c/    ----------------------------------------------------
      
      if (nph_pnts.GT.maxph) nph_pnts = maxph
      mExp=2*idp1+1+inon1


      tanth = tan (theta_ij)
      sinth = sin (theta_ij)
      costh = cos (theta_ij)
      x1=x/L_i-0.5

      if (i_geom.EQ.1) then
	 phi_min = atan (tanth / (1.0 - x * tanth / (L_j * sinth)))
	 phi_max = pi

      else if (i_geom.EQ.2) then

         sina = sin (alpha_j)
         cosa = cos (alpha_j)
         phi_min = atan ((L_perp + L_j * sina) /
     .      (L_perp / tanth  + L_j * cosa - x))
         phi_max = atan (L_perp / (L_perp / tanth - x))

      endif

      if (phi_min.LT.0.0) phi_min = pi + phi_min
      if (phi_max.LT.0.0) phi_max = pi + phi_max

      if(phi_min.gt.phi_max) phi_min=0.0

      d_phi = (phi_max - phi_min) / float (nph_pnts - 1)

      do i = 1, nph_pnts
	 phi = phi_min + (i-1) * d_phi
	 if (i_geom.EQ.1) then
	    el_of_phi = x * sinth / sin (phi - theta_ij)
	 else
           if(phi.gt.1.0e-10)then
	      ksi_j =((L_perp /tan(theta_ij) - x) * tan (phi)-L_perp) /
     .         (sina - cosa * tan(phi))
	      el_of_phi = (L_perp + ksi_j*sina)/ sin (phi)
           else
              el_of_phi=L_perp*costh/sinth-x-L_perp*cosa/sina
           endif
	 endif

         x0=el_of_phi/l_mfp
         phi0=phi-pi/2
         cphi0=cos(phi0)
         ki3=bickley3(x0)
         x2=sin(phi_max-phi)*el_of_phi/sin(phi_max-alpha_j)/L_j-0.5


          

         tvec(i,1)=cphi0*ki3




         if(numfun.gt.1)then
           if(inon1.eq.0)then
             phi1=phi0-alpha_j
             sphi0=sin(phi0)
             cphi1=cos(phi1)
             sphi1=sin(phi1)
             ki1=bickley1(x0)
             ki2=bickley2(x0)
             ki4=(2.0*ki2+x0*(ki1-ki3))/3.0
             tvec(i,2)=sphi1*cphi0*ki4
             tvec(i,3)=cphi0*cphi1*ki4
             if(numfun.gt.mExp)then
               ki5=(3.0*ki3+x0*(ki2-ki4))/4.0
               tvec(i,4)=sphi0*cphi0*ki4
               tvec(i,5)=sphi1*sphi0*cphi0*ki5
               tvec(i,6)=cphi1*sphi0*cphi0*ki5
               tvec(i,7)=cphi0*cphi0*ki4
               tvec(i,8)=sphi1*cphi0*cphi0*ki5
               tvec(i,9)=cphi1*cphi0*cphi0*ki5
             endif
           else    !if inon==1
             x2=sin(phi_max-phi)*el_of_phi/sin(phi_max-alpha_j)/L_j-0.5
             if(idp1.eq.0)then
               tvec(i,2)=cphi0*ki3*x2
               tvec(i,3)=cphi0*ki3*x1
               tvec(i,4)=cphi0*ki3*x1*x2
             else
               phi1=phi0-alpha_j
               sphi0=sin(phi0)
               cphi1=cos(phi1)
               sphi1=sin(phi1)
               ki1=bickley1(x0)
               ki2=bickley2(x0)
               ki4=(2.0*ki2+x0*(ki1-ki3))/3.0
               ki5=(3.0*ki3+x0*(ki2-ki4))/4.0
               tvec(i,2)=cphi0*ki3*x2
               tvec(i,3)=sphi1*cphi0*ki4
               tvec(i,4)=cphi0*cphi1*ki4
               tvec(i,5)=cphi0*ki3*x1
               tvec(i,6)=cphi0*ki3*x1*x2
               tvec(i,7)=sphi1*cphi0*ki4*x1
               tvec(i,8)=cphi0*cphi1*ki4*x1
               tvec(i,9)=sphi0*cphi0*ki4
               tvec(i,10)=sphi0*cphi0*ki4*x2
               tvec(i,11)=sphi1*sphi0*cphi0*ki5
               tvec(i,12)=cphi1*sphi0*cphi0*ki5
               tvec(i,13)=cphi0*cphi0*ki4
               tvec(i,14)=cphi0*cphi0*ki4*x2
               tvec(i,15)=sphi1*cphi0*cphi0*ki5
               tvec(i,16)=cphi1*cphi0*cphi0*ki5

             endif  !End of if idp1
               


           endif    !End of if inon 
         endif      !End of if numfun>1
          

      enddo

c/    Perform integration wrt the angle phi:
c/    -------------------------------------
      call simpson (tvec, d_phi, nph_pnts, reslt,numfun)
      
      do j=1,numfun
         funval(j) = reslt(j)
      enddo

      return
      end
