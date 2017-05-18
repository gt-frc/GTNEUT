      subroutine calcRect (tkji,iquad1,numfun)

c/    This subroutine calculates the transmission coefficient
c/    used in the TEP theory for rectangular regions.

c/    01/13/98, jm: new version using direct integration over the 
c/                  angle -phi-
c/    06/05/98, jm: options for # of points in quadrature integrations
      
      implicit none
      include 'locGeom.inc'
      include 'neutGlob.inc'

      integer iquad1,numfun,i,j,k
      real tkji(numfun), ss(numfun), pi,coef(maxExp), temp,temp1

      external t_ij

      data pi /3.1415926/
      save pi

      data coef /1,3.46410162,2,4.2426406871/
      save coef


c/    Call the appropriate quadrature integration routine:
c/    ZWF added new quadrature set n = 5,10
      if (iquad1.EQ.55) then
         call qgauss6 (t_ij, 0.0, L_i, ss,numfun)
      else if (iquad1.EQ.0) then
         call qgauss10 (t_ij, 0.0, L_i, ss,numfun)
      else if (iquad1.EQ.1) then
         call qgauss20 (t_ij, 0.0, L_i, ss,numfun)
      else if (iquad1.EQ.2) then
         call qgauss40 (t_ij, 0.0, L_i, ss,numfun)      
      else if (iquad1.EQ.3) then
         call qgauss60 (t_ij, 0.0, L_i, ss,numfun)
      else if (iquad1.EQ.4) then
         call qgauss80 (t_ij, 0.0, L_i, ss,numfun)
      else if (iquad1.EQ.5) then
         call qgauss100 (t_ij, 0.0, L_i, ss,numfun)
      endif

      if(inon.eq.0)then
        if(numfun.eq.1)then
           tkji(1)=2.0*coef(1)*ss(1)/(pi*L_i)
           return
        endif
        if(numfun.eq.mExp)then
          do i=1,mExp
             if(i.eq.1)then
               temp=coef(i)
             else
               temp=coef(i+1)
             endif
             tkji(i)=2.0*temp*ss(i)/(pi*L_i)
          enddo
          tkji(3)=tkji(3)-2.8284271247*tkji(1)
          return
        endif
        if(numfun.eq.mExp*mExp)then
          do i=1,mExp
             if(i.eq.1)then
               temp=coef(i)
             else
               temp=coef(i+1)
             endif
             do j=1,mExp
                 k=(i-1)*mExp+j
                 if(j.eq.1)then
                   temp1=coef(1)
                 else
                   temp1=coef(j+1)
                 endif
                 tkji(k)=2.0*temp*temp1*
     .             ss(k)/(pi*L_i)
              enddo
           enddo
           tkji(3)=tkji(3)-2.8284271247*tkji(1)
           tkji(6)=tkji(6)-2.8284271247*tkji(4)
           tkji(9)=tkji(9)-2.8284271247*tkji(7)
     .                    -2.8284271247*tkji(3)
           tkji(7)=tkji(7)-2.8284271247*tkji(1)
           tkji(8)=tkji(8)-2.8284271247*tkji(2)
        endif
      else
        do i=1, mExp
           do j=1, mExp
              k=(i-1)*mExp+j
              tkji(k)=2.0*coef(i)*coef(j)*ss(k)/(pi*L_i)
           enddo
        enddo
        if(idp.eq.1)then
           tkji(4)=tkji(4)-2.8284271247*tkji(1)
           tkji(8)=tkji(8)-2.8284271247*tkji(5)
           tkji(12)=tkji(12)-2.8284271247*tkji(9)
           tkji(16)=tkji(16)-2.8284271247*tkji(13)
     .                    -2.8284271247*tkji(4)
           tkji(13)=tkji(13)-2.8284271247*tkji(1)
           tkji(14)=tkji(14)-2.8284271247*tkji(2)
           tkji(15)=tkji(15)-2.8284271247*tkji(3)


        endif  !End of if idp==1
      endif    !End of if inon


      return
      end
