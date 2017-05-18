      subroutine simpson(f, h, n, s, m)

c///////////////////////////////////////////////////////////////////////
c/                                                                     /
c/    This subroutine calculates the integral of the function -f-      /
c/    which is defined at -n- discrete points with equal step -h-      /
c/    using  Simpson's rule. If the number of points -n- is even,      /
c/    the trapezoidal rule is used for the last two points.            /
c/    Note: use this routine for cases, where the function is given    /
c/    in discrete points. If an analytic form of the function is       /
c/    available, it's better to use an adaptive scheme.                /
c/    Created by John Mandrekas, GIT, 1/14/92                          /
c/    Cosmetic changes, 06/21/95                                       /
c/                                                                     /
c///////////////////////////////////////////////////////////////////////

      implicit none
      integer i, n, m, j,max
      parameter (max=51)
      real s(m), h, f(max,m)

      do j=1,m
         s(j) = 0.0
         if(n.LE.1) return

         do i = 3, n, 2
            s(j)=s(j)+f(i-2,j)+ 4.0 * f(i-1,j)+f(i,j)
         enddo
         s(j) = s(j) / 3.0

         if (mod(n,2).EQ.0) then
            s(j)=s(j)+ 0.5*(f(n-1,j)+f(n,j))
         endif

         s(j) = s(j) * h
      enddo

      return
      end
