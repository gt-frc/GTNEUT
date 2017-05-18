      subroutine reflect (e0, am1, am2, z1, z2, rn, re)

c/    This subroutine calculates the particle and energy reflection
c/    coefficients for normal incidence, using fits depending on the
c/    projectile and target properties.
c/    
c/    12/18/01, Created by John Mandrekas and Dingkang Zhang

c/    References:
c/    ----------
c/    W. Eckstein, J. Nucl. Materials, 248 (1997) 1
c/    E.W. Thomas, R.K. Janev and J. Smith, Nucl. Instrum. & Methods,
c/       B69 (1992) 427.
c/
c/    Parameters:
c/    ----------
c/    e0    : projectile energy in keV
c/    am1   : projectile mass (amu)
c/    am2   : target mass
c/    z1    : projectile atomic number
c/    z2    : target atomic number
c/    rn    : particle reflection coefficient, RN
c/    re    : energy reflection coefficient, RE
c/
c/    Notice, that the average energy of reflected (backscattered) neutrals
c/    is e0 * RE / RN, where RE and RN are calculated at E = e0


      implicit none
      integer i, indx
      real e0, am1, am2, z1, z2, rn, re
      real e, mu, zfactr, epsln
      real an1(5), an2(5), an3(5), an4(5), an5(5), an6(5)
      real ae1(5), ae2(5), ae3(5), ae4(5), ae5(5), ae6(5)
      real r_n, r_e

      data e /2.71828/
      data an1 /0.02129, 0.3680, 0.5173, 0.6192, 0.8250/
      data ae1 /0.001445, 0.2058, 0.4222, 0.4484, 0.6831/
      data an2 /16.39, 2.985, 2.549, 20.01, 21.41/
      data ae2 /404.7, 3.848, 3.092, 27.16, 27.16/
      data an3 /26.39, 7.122, 5.325, 8.922, 8.606/
      data ae3 /73.73, 19.07, 13.17, 15.66, 15.66/
      data an4 /0.9131, 0.5802, 0.5719, 0.6669, 0.6425/
      data ae4 /0.6519, 0.4872, 0.5393, 0.6598, 0.6598/
      data an5 /6.249, 4.211, 1.094, 1.864, 1.907/
      data ae5 /34.66, 15.13, 4.464, 7.967, 7.967/
      data an6 /2.550, 1.597, 1.933, 1.899, 1.927/
      data ae6 /1.971, 1.638, 1.877, 1.822, 1.822/

      save e


      r_n(i) = an1(i) * alog (an2(i)*epsln + e) / 
     .   (1.+ an3(i)*epsln**an4(i) + an5(i)*epsln**an6(i))

      r_e(i) = ae1(i) * alog (ae2(i)*epsln + e) / 
     .   (1.+ ae3(i)*epsln**ae4(i) + ae5(i)*epsln**ae6(i))

c/    Calculate reduced energy:
      
      mu = am2 / am1

      zfactr = 1.0 / (z1 * z2 * sqrt(z1**0.67 + z2**0.67))
      epsln = 32.55 * mu * zfactr * e0 / (1.+mu)

      indx = 0
      if (mu.EQ.1) then
	 indx = 1
      else if (mu.EQ.3)  then
	 indx = 2
      else if ((mu.GE.6.0).AND.(mu.LE.7.0)) then
	 indx = 3
      else if ((mu.GE.12.0).AND.(mu.LE.15.0)) then
	 indx = 4
      else if (mu.GE.20.0) then
	 indx = 5
      endif

      if (indx.NE.0) then
	 rn = r_n(indx)
	 re = r_e(indx)
      endif

      return
      end
