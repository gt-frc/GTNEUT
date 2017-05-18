	subroutine FEM(i,nd)
c/    This subroutine use Finite Element Method (FEM) to solve diffuision
c/    equation
        implicit none
	include 'neutGlob.inc'
	include 'esc.inc'

	external calssource

c/      local variables
	integer i, nd, j, k, j0, k1, l, m, l1, l2, lt
	integer n,n1,l0, ld, i0, np, n0, del
	integer info, lda, ldb, nrhs, ipiv(maxFEM), neq
        integer ktype, kw, ns
	real xtr(3,maxSides), ytr(3,maxSides), H1(maxSides,6),
     .    H2(maxSides,6), H3(maxSides,3), a(3,maxSides),
     .    b0(3,maxSides), c(3,maxSides),sva, svtr, ang0, pi,
     .    xp1, yp1, xp2, yp2, current(maxSides, 4), totcur,
     .    xb1, xb2, temp, temps, temp0, current0(maxSides,4),
     .    currents(maxSides,4),totcur0,totcurs,lam_old
        double precision SE(maxFEM,6*maxSides)

	real l02, phimax, ssource(2)

	svtr=1.0/lmfp(i)
	sva=svtr*(1-A_cx(i))
        ns=nSides(i)

	pi=3.14159254
c/    Assign the coordinates for each sub-region (triangle)
	xp(1)=lside(1,i)
	yp(1)=0
        xc=xp(1)
	yc=yp(1)
	sublength(1)=lside(1,i)
	ang0=pi
	do j=2, ns
	   xp(j)=xp(j-1)+lside(j-1,i)*cos(ang0)
	   yp(j)=yp(j-1)+lside(j-1,i)*sin(ang0)
	   xc=xc+xp(j)
	   yc=yc+yp(j)
	   sublength(j)=lside(j,i)
	   j0=j-1
	   if(j0.lt.1) j0=j0+ns
	   ang0=ang0-(pi-angle(j0,i))
	enddo
	xc=xc/ns
	yc=yc/ns
	do j=1, ns
	   j0=j+1
	   if(j0.gt.ns)j0=j0-ns
	   xtr(1,j)=xp(j)
	   ytr(1,j)=yp(j)
	   xtr(2,j)=xp(j0)
	   ytr(2,j)=yp(j0)
	   xtr(3,j)=xc
	   ytr(3,j)=yc

c/    calculate the area of each triangle
	  subarea(j)=0
	  do k=1,3
	     k1=k+1
	     if(k1.gt.3)k1=k1-3
	     subarea(j)=subarea(j)+xtr(k,j)*ytr(k1,j)
     .                -xtr(k1,j)*ytr(k,j)
	  enddo
	  subarea(j)=-subarea(j)/2
	enddo

	do m=1, ns
	   do l=1,3
	      l1=l+1
	      l2=l+2
	      if(l1.gt.3)l1=l1-3
	      if(l2.gt.3)l2=l2-3
	      a(l,m)=-(xtr(l1,m)*ytr(l2,m)-xtr(l2,m)*ytr(l1,m))
     .		  /(2*subarea(m))
	      b0(l,m)=-(ytr(l1,m)-ytr(l2,m))/(2*subarea(m))
	      c(l,m)=-(xtr(l2,m)-xtr(l1,m))/(2*subarea(m))
	   enddo
	enddo

c/    calcuate the average qualities
      do m=1, ns
	   xa2(m)=(xtr(1,m)**2+xtr(2,m)**2+xtr(3,m)**2+xtr(1,m)*xtr(2,m)
     .         +xtr(1,m)*xtr(3,m)+xtr(2,m)*xtr(3,m))/6
	   ya2(m)=(ytr(1,m)**2+ytr(2,m)**2+ytr(3,m)**2+ytr(1,m)*ytr(2,m)
     .         +ytr(1,m)*ytr(3,m)+ytr(2,m)*ytr(3,m))/6
	   xya(m)=(2*xtr(1,m)*ytr(1,m)+2*xtr(2,m)*ytr(2,m)+
     .	    2*xtr(3,m)*ytr(3,m)+xtr(1,m)*ytr(2,m)+xtr(2,m)*ytr(1,m)
     .        +xtr(1,m)*ytr(3,m)+xtr(3,m)*ytr(1,m)
     .        +xtr(2,m)*ytr(3,m)+xtr(3,m)*ytr(2,m))/12
	   xa(m)=(xtr(1,m)+xtr(2,m)+xtr(3,m))/3
	   ya(m)=(ytr(1,m)+ytr(2,m)+ytr(3,m))/3
         xl2(m)=(xtr(1,m)**2+xtr(2,m)**2+xtr(1,m)*xtr(2,m))/3
         yl2(m)=(ytr(1,m)**2+ytr(2,m)**2+ytr(1,m)*ytr(2,m))/3
	   xyl(m)=(2*xtr(1,m)*ytr(1,m)+2*xtr(2,m)*ytr(2,m)
     .       +xtr(1,m)*ytr(2,m)+xtr(2,m)*ytr(1,m))/6
	   xl(m)=(xtr(1,m)+xtr(2,m))/2
	   yl(m)=(ytr(1,m)+ytr(2,m))/2
	enddo


c/    Calculate the elements of Diffusion matrix
c/    First, calculate the matrix for the big triangles
      do m=1, ns
	   do l=1,6
	      if(l.lt.4)then
	        l1=l
	        l2=l
	      else
		    if(l.lt.6) then
	          l1=1
	          l2=l-2
	        else
	          l1=2
	          l2=3
	        endif
	      endif
	      H1(m,l)=subarea(m)*(b0(l1,m)*b0(l2,m)+c(l1,m)*c(l2,m))
     .              /(3*svtr)
		  H2(m,l)=subarea(m)*(a(l1,m)*a(l2,m)+xa2(m)*b0(l1,m)*b0(l2,m)
     .       +ya2(m)*c(l1,m)*c(l2,m)+(b0(l1,m)*c(l2,m)+b0(l2,m)*c(l1,m))
     .        *xya(m)+(a(l1,m)*b0(l2,m)+a(l2,m)*b0(l1,m))*xa(m)
     .        +(a(l1,m)*c(l2,m)+a(l2,m)*c(l1,m))*ya(m))/nd**2
	   enddo
	   do l=1,3
	      if(l.lt.3)then
	        l1=l
	        l2=l
	      else
	        l1=1
	        l2=2
	      endif
	      H3(m,l)=0.5*sublength(m)*(a(l1,m)*a(l2,m)+b0(l1,m)*b0(l2,m)
     .        *xl2(m)+c(l1,m)*c(l2,m)*yl2(m)+(a(l1,m)*b0(l2,m)
     .        +a(l2,m)*b0(l1,m))*xl(m)+(a(l1,m)*c(l2,m)+a(l2,m)*c(l1,m))
     .         *yl(m)+(b0(l1,m)*c(l2,m)+b0(l2,m)*c(l1,m))*xyl(m))/nd
	   enddo
      enddo

c/    Calculate the total points
	lt=1+ns*nd*(nd+1)/2
	
c/    Initialize matrix
	do l=1,lt
	   do m=1,lt
	      B1(l,m)=0
		  A1(l,m)=0
	   enddo
	enddo
	
	do n=1,ns
	   l=n+1
	   l0=1
	   n1=n+1
	   if(n1.gt.ns) n1=n1-ns
	   l1=n1+1
	   B1(l,l)=B1(l,l)+H2(n,1)
	   B1(l1,l1)=B1(l1,l1)+H2(n,2)
	   A1(l,l)=A1(l,l)+sva*H2(n,1)+H1(n,1)
	   A1(l1,l1)=A1(l1,l1)+sva*H2(n,2)+H1(n,2)

	   B1(l,l0)=B1(l,l0)+H2(n,5)
	   B1(l1,l0)=B1(l1,l0)+H2(n,6)
	   B1(l0,l)=B1(l0,l)+H2(n,5)
	   B1(l0,l1)=B1(l0,l1)+H2(n,6)
	   A1(l,l0)=A1(l,l0)+sva*H2(n,5)+H1(n,5)
	   A1(l1,l0)=A1(l1,l0)+sva*H2(n,6)+H1(n,6)
	   A1(l0,l)=A1(l0,l)+sva*H2(n,5)+H1(n,5)
	   A1(l0,l1)=A1(l0,l1)+sva*H2(n,6)+H1(n,6)

	   B1(l0,l0)=B1(l0,l0)+H2(n,3)
	   A1(l0,l0)=A1(l0,l0)+H1(n,3)+sva*H2(n,3)
	   
	   B1(l,l1)=B1(l,l1)+H2(n,4)
	   B1(l1,l)=B1(l1,l)+H2(n,4)
	   A1(l,l1)=A1(l,l1)+sva*H2(n,4)+H1(n,4)  			    
 	   A1(l1,l)=A1(l1,l)+sva*H2(n,4)+H1(n,4)  
	enddo
	
	ld=4
	do i0=2,nd
	   do n=1, ns
	      do m=1,i0
		     l0=(n-1)*i0+m
			 l=l0+ns*(i0-1)*i0/2+1
			 l1=l+1
			 l2=l-ld
			 if(l0.eq.i0*ns)then
			   l1=l1-i0*ns
			   l2=l2-(i0-1)*ns
			 endif
			 B1(l,l)=B1(l,l)+H2(n,1)
			 B1(l1,l1)=B1(l1,l1)+H2(n,2)
			 B1(l2,l2)=B1(l2,l2)+H2(n,3)
			 A1(l,l)=A1(l,l)+sva*H2(n,2)+H1(n,1)
			 A1(l1,l1)=A1(l1,l1)+sva*H2(n,2)+H1(n,2)    
			 A1(l2,l2)=A1(l2,l2)+sva*H2(n,3)+H1(n,3)
			 
			 B1(l,l1)=B1(l,l1)+H2(n,4)    
			 B1(l1,l)=B1(l1,l)+H2(n,4)    
			 B1(l,l2)=B1(l,l2)+H2(n,5)    
			 B1(l2,l)=B1(l2,l)+H2(n,5)    
			 B1(l1,l2)=B1(l1,l2)+H2(n,6)    
			 B1(l2,l1)=B1(l2,l1)+H2(n,6)  
			 A1(l,l1)=A1(l,l1)+sva*H2(n,4)+H1(n,4)  
			 A1(l1,l)=A1(l1,l)+sva*H2(n,4)+H1(n,4)  
			 A1(l,l2)=A1(l,l2)+sva*H2(n,5)+H1(n,5)  
			 A1(l2,l)=A1(l2,l)+sva*H2(n,5)+H1(n,5)  
			 A1(l1,l2)=A1(l1,l2)+sva*H2(n,6)+H1(n,6)  
			 A1(l2,l1)=A1(l2,l1)+sva*H2(n,6)+H1(n,6)
		  enddo  

	      do m=1,i0-1
	         l0=(n-1)*i0+m+1
	         l2=l0+ns*(i0-1)*i0/2+1
	         l=l2-ld
	         l1=l-1
	         if(l0.eq.i0*ns)l=l-(i0-1)*ns
	         B1(l,l)=B1(l,l)+H2(n,1)
	         B1(l1,l1)=B1(l1,l1)+H2(n,2)
	         B1(l2,l2)=B1(l2,l2)+H2(n,3)
	         A1(l,l)=A1(l,l)+sva*H2(n,1)+H1(n,1)
	         A1(l1,l1)=A1(l1,l1)+sva*H2(n,2)+H1(n,2)
	         A1(l2,l2)=A1(l2,l2)+sva*H2(n,3)+H1(n,3)


	         B1(l,l1)=B1(l,l1)+H2(n,4)
	         B1(l1,l)=B1(l1,l)+H2(n,4)
	         B1(l,l2)=B1(l,l2)+H2(n,5)
	         B1(l2,l)=B1(l2,l)+H2(n,5)
	         B1(l1,l2)=B1(l1,l2)+H2(n,6)
	         B1(l2,l1)=B1(l2,l1)+H2(n,6)
	         A1(l,l1)=A1(l,l1)+sva*H2(n,4)+H1(n,4)
	         A1(l1,l)=A1(l1,l)+sva*H2(n,4)+H1(n,4)
	         A1(l,l2)=A1(l,l2)+sva*H2(n,5)+H1(n,5)
	         A1(l2,l)=A1(l2,l)+sva*H2(n,5)+H1(n,5)
	         A1(l1,l2)=A1(l1,l2)+sva*H2(n,6)+H1(n,6)
	         A1(l2,l1)=A1(l2,l1)+sva*H2(n,6)+H1(n,6)
	      enddo
	      ld=ld+1

	   enddo
	enddo	
	
C/	set up the boundary
	do n=1, ns
	   if(n.eq.ns)then
	     n1=1
	   else
	     n1=n+1
	   endif
	   do i0=1,nd
	      l0=(n-1)*nd+i0
		  l=l0+ns*(nd-1)*nd/2+1
		  l1=l+1
		  if(l0.eq.nd*ns) l1=l1-nd*ns
	      A1(l,l)=A1(l,l)+H3(n,1)
		  A1(l1,l1)=A1(l1,l1)+H3(n,2)
		  
		  A1(l,l1)=A1(l,l1)+H3(n,3)
		  A1(l1,l)=A1(l1,l)+H3(n,3)
	   enddo
	enddo
	

c/    Calculate the first collision source density
c/    First, calculate the coordinates for each point
	xt(1)=xc
	yt(1)=yc
	np=1
	do n=1,nd
	   do l=1,ns
	      xp1=xc+(xtr(1,l)-xc)*n/nd
	      yp1=yc+(ytr(1,l)-yc)*n/nd
	      xp2=xc+(xtr(2,l)-xc)*n/nd
	      yp2=yc+(ytr(2,l)-yc)*n/nd
	      do j=1,n
	         np=np+1
	         xt(np)=xp1+(j-1)*(xp2-xp1)/n
                 yt(np)=yp1+(j-1)*(yp2-yp1)/n
	      enddo
	   enddo
	enddo

	n0=ns*(nd-1)*nd/2+1

        do l=1,lt
           do m=1,6*ns
              s0(l,m)=0.0
           enddo
        enddo

	do m=1, ns
           ktype=iType(adjCell(m,i))
	   if(kType.lt.2)then
	     lam0=mfp(m,i)
	   else
             kw=adjCell(m,i)-(nCells+nPlasmReg)
	     lam0=mfp_wf(kw)
	   endif


	   do l=1,lt
	      del=l-n0-(m-1)*nd
	      if(m.eq.ns.and.(l.eq.(n0+1)))del=del+ns*nd
	      if((del.gt.0).and.(del.le.(nd+1)))then
	        s0(l,2*m-1)=1
	        s0(l,2*m)=(0.5-(del-1.0)/nd)*3.46410162
                if(ktype.eq.2)then
                  if(irefl.eq.1.and.(zwall(kw).gt.0))then
                    s0(l,2*m-1+2*ns)=1
                    s0(l,2*m+2*ns)=s0(l,2*m)
                  endif
                  if(g_ex(kw).gt.0)then
                    s0(l,2*m-1+4*ns)=1
                    s0(l,2*m+4*ns)=s0(l,2*m)
                  endif
                endif
	      else
	        l01=sqrt((xtr(1,m)-xt(l))**2+(ytr(1,m)-yt(l))**2)
	        l12=sublength(m)
	        l02=sqrt((xtr(2,m)-xt(l))**2+(ytr(2,m)-yt(l))**2)
	        phi012=acos((l01**2+l12**2-l02**2)/(2*l01*l12))
	        phimax=acos((l01**2+l02**2-l12**2)/(2*l01*l02))
                call qgauss20(calssource,0,phimax,ssource,2)
	        S0(l,2*m-1)=ssource(1)/pi
	        S0(l,2*m)=ssource(2)/pi*3.46410162
                if(ktype.eq.2)then
                  if(irefl.eq.1.and.(zwall(kw).gt.0))then
                    lam_old=lam0
                    lam0=mfp_ws(kw)
                    call qgauss20(calssource,0,phimax,ssource,2)
                    S0(l,2*m-1+2*ns)=ssource(1)/pi
                    S0(l,2*m+2*ns)=ssource(2)/pi*3.46410162
                    lam0=lam_old
                  endif
                  if(g_ex(kw).gt.0)then
                    lam_old=lam0
                    lam0=mfp_w0(kw)
                    call qgauss20(calssource,0,phimax,ssource,2)                
                    S0(l,2*m-1+4*ns)=ssource(1)/pi
                    S0(l,2*m+4*ns)=ssource(2)/pi*3.46410162
                    lam0=lam_old
                  endif
                endif
	      endif
	   enddo
	enddo

	do l1=1, lt
	   do m=1, 6*ns
	      ST(l1,m)=0.0
	      do l2=1,lt
	         ST(l1,m)=ST(l1,m)+B1(l1,l2)*S0(l2,m)
	      enddo
	      SE(l1,m)=ST(l1,m)
	   enddo
	enddo

	lda=maxFEM
	ldb=maxFEM
	neq=lt
	nrhs=6*ns

	call DGESV(neq, nrhs, A1, lda, ipiv, SE, ldb, info)

        if(info.ne.0) then
          write (*, *)'inf0=', info
        endif

	n0=ns*nd*(nd-1)/2+1
	do m=1, ns
           ktype=iType(adjCell(m,i))
           if(ktype.eq.2)kw=adjCell(m,i)-(nCells+nPlasmReg)
	   totcur=0
           totcurs=0
           totcur0=0
	   do n=1, ns
	      do j=1,2
	         current(n,2*j-1)=0
	         current(n,2*j)=0
                 currents(n,2*j-1)=0
                 currents(n,2*j)=0
                 current0(n,2*j-1)=0
                 current0(n,2*j)=0
	         do i0=1, nd
	            l0=(n-1)*nd+i0
	            l=l0+n0
	            l1=l+1
	            if(l0.eq.nd*ns)l1=l1-nd*ns
	            j0=j+(m-1)*2
	            current(n,2*j-1)=current(n,2*j-1)+SE(l,j0)+SE(l1,j0)
		    xb1=(i0-1.0)/nd-0.5
		    xb2=xb1+1.0/nd
                    temp=2*(xb1*SE(l,j0)+xb2*SE(l1,j0))
     .                 +xb1*SE(l1,j0)+xb2*SE(l,j0)
		    current(n,2*j)=current(n,2*j)+temp
                    if(ktype.eq.2)then
                      if(irefl.eq.1.and.zwall(kw).gt.0)then
                        currents(n,2*j-1)=currents(n,2*j-1)
     .                    +SE(l,j0+2*ns)+SE(l1,j0+2*ns)
                        temps=2*(xb1*SE(l,j0+2*ns)+xb2*SE(l1,j0+2*ns))
     .                    +xb1*SE(l1,j0+2*ns)+xb2*SE(l,j0+2*ns)
                        currents(n,2*j)=currents(n,2*j)+temps
                      endif
                      if(g_ex(kw).gt.0)then
                        current0(n,2*j-1)=current0(n,2*j-1)
     .                    +SE(l,j0+4*ns)+SE(l1,j0+4*ns)
                        temp0=2*(xb1*SE(l,j0+4*ns)+xb2*SE(l1,j0+4*ns))
     .                    +xb1*SE(l1,j0+4*ns)+xb2*SE(l,j0+4*ns)
                        current0(n,2*j)=current0(n,2*j)+temp0
                      endif
                    endif
	         enddo
		 current(n,2*j-1)=current(n,2*j-1)*sublength(n)/(2.0*nd)
		 current(n,2*j)=current(n,2*j)*sublength(n)
     .		        /(1.73205081*nd)
                 if(ktype.eq.2)then
                   if(irefl.eq.1.and.zwall(kw).gt.0)then
                     currents(n,2*j-1)=currents(n,2*j-1)*
     .                  sublength(n)/(2.0*nd)
                     currents(n,2*j)=currents(n,2*j)*sublength(n)
     .                  /(1.73205081*nd)
                   endif
                   if(g_ex(kw).gt.0)then
                     current0(n,2*j-1)=current0(n,2*j-1)*
     .                  sublength(n)/(2.0*nd)
                     current0(n,2*j)=current0(n,2*j)*sublength(n)
     .                  /(1.73205081*nd)
                   endif
                 endif
	      enddo

	      totcur=totcur+current(n,1)
              if(ktype.eq.2)then
                if(irefl.eq.1.and.zwall(kw).gt.0)then
                  totcurs=totcurs+currents(n,1)
                endif
                if(g_ex(kw).gt.0)then
                  totcur0=totcur0+current0(n,1)
                endif
              endif
	   enddo

	   do n=1, ns
	      do j=1, 4
	         lambdak(m,i,n,j)=current(n,j)/totcur
                 if(ktype.eq.2)then
                   lambdawf(kw,n,j)=lambdak(m,i,n,j)
                   if(irefl.eq.1.and.zwall(kw).gt.0)then
                     lambdaws(kw,n,j)=currents(n,j)/totcurs
                   endif
                   if(g_ex(kw).gt.0)then
                     lambdaw0(kw,n,j)=current0(n,j)/totcur0
                   endif
                 endif
	      enddo
	   enddo

	   pEscpk(m,i)=0
           if(ktype.eq.2)then
             pEscpws(kw)=0
             pEscpw0(kw)=0
           endif
	   do l=1,lt
	      pEscpk(m,i)=pEscpk(m,i)+ST(l,(m-1)*2+1)
	   enddo
	   pEscpk(m,i)=totcur/(2*pEscpk(m,i))
           if(ktype.eq.2)then
             pEscpwf(kw)=pEscpk(m,i)
             if(irefl.eq.1.and.zwall(kw).gt.0)then
               do l=1,lt
                  pEscpws(kw)=pEscpws(kw)+ST(l,(m-1)*2+1+2*ns)
               enddo
               pEscpws(kw)=totcurs/(2*pEscpws(kw))
             endif
             if(g_ex(kw).gt.0)then
               do l=1,lt
                  pEscpw0(kw)=pEscpw0(kw)+ST(l,(m-1)*2+1+4*ns)
               enddo
               pEscpw0(kw)=totcur0/(2*pEscpw0(kw))
             endif
           endif
	enddo

	if(i.eq.-5)then
	write(*,*)'xp,yp='
	write(*,*)(xp(j),j=1,ns)
	write(*,*)(yp(j),j=1,ns)

	write(*,*)'length,area='
	write(*,*)(sublength(m),m=1,ns)
	write(*,*)(subarea(m),m=1,ns)

	write(*,*)'H1='
	do m=1,ns
	   write(*,*)(H1(m,j),j=1,6)
	enddo
	write(*,*)'H2='
	do m=1,ns
	   write(*,*)(H2(m,j),j=1,6)
	enddo
	write(*,*)'H3='
	do m=1,ns
	   write(*,*)(H3(m,j),j=1,3)
	enddo


	write(*,*)'A1='
	do m=1,lt
	   write(*,*)(A1(m,l),l=1,lt)
	enddo
	write(*,*)'B1='
	do m=1,lt
	   write(*,*)(B1(m,l),l=1,lt)
	enddo
	write(*,*)'S0=',S0
	do m=1,2*ns
	write(*,*)(S0(l,m),l=1,lt)
	enddo
	write(*,*)'ST='
	do m=1,2*ns
	write(*,*)(ST(l,m),l=1,lt)
	enddo
	write(*,*)'SE='
	do l=1,lt
	write(*,*)(SE(l,m),m=1,2*ns)
	enddo


	endif

	if(i.eq.-5)then
	write(*,*)'current='
	do m=1,ns
	write(*,*)(current(m,j),j=1,4)
	enddo
	do m=1,ns
	write(*,*)'m=',m
	write(*,*)'pEscpk=',pEscpk(m,i)
	do l=1,ns
	write(*,*)'from m=',m,' to l=',l, ' lambdak(m,i,l,j)='
	write(*,*)(lambdak(m,i,l,j),j=1,4)
	enddo
	enddo
        write(*,*)'pEscpwf=',pEscpwf(kw)
        write(*,*)'pEscpws=',pEscpws(kw)
        write(*,*)'pEscpw0=',pEscpw0(kw)

        do m=1,ns
           write(*,*)'m=',m,' lambdawf=',(lambdawf(kw,m,n),n=1,4)
        enddo
        do m=1,ns
           write(*,*)'m=',m,' lambdaws=',(lambdaws(kw,m,n),n=1,4)
        enddo
        do m=1,ns
           write(*,*)'m=',m,' lambdaw0=',(lambdaw0(kw,m,n),n=1,4)
        enddo

	endif

1000    format (1x, 'ERROR IN DGESV (LAPACK)! Value of INFO = ', i2)

	end	 
	
	
	subroutine calssource(x,numfun,funval)

	implicit none
		
	include 'neutGlob.inc'
	include 'esc.inc'	    
	integer numfun
	real x, funval(numfun)
	real bickley2, l, x1

	l=sin(phi012)/sin(x+phi012)*l01/lam0
	x1=0.5-sin(x)/sin(x+phi012)*l01/l12
	funval(1)=bickley2(l)
	if(numfun.gt.1)funval(2)=x1*funval(1)

	return
	end

	



