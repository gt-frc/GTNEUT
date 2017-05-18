      subroutine rectinp

c/    This subroutine sets up the required input parameters for a 
c/    rectangular region bounded by a wall and with uniform background
c/    plasma parameters. The description of the input variables affecting 
c/    the input generation are repeated here for convenience.

c/    05/17/2000, John Mandrekas, created.
c/   

c/    Lx       : Length of x-side of rectangle (m)
c/    Ly       : Length of y-side of rectangle (m)
c/    NX       : Number of cells in the x-direction
c/    NY       : Number of cells in the y-direction
c/    ne_fixed : Background plasma electron density (/m^3) if no
c/               gradients are present. If no gradients are present and
c/               if ne_fixed <= 0, then the n_e for each cell is read from
c/               the input file.
c/    igradneh : If equal to 1, linear variation in the electron density
c/               in the horizontal direction (X) determined by the values
c/               of ne_lft and ne_rgt (see below)
c/    igradnev : If equal to 1, linear variation in the electron density
c/               in the vertical direction (Y) determined by the values
c/               of ne_btm and ne_top (see below)
c/    ne_lft   : Electron density at the left boundary (/m3)
c/    ne_rgt   : Electron density at the right boundary (/m3)
c/    ne_top   : Electron density at the top boundary (/m3)
c/    ne_btm   : Electron density at the bottom boundary (/m3)
c/               The same holds for the other background plasma quantities
c/               which are assigned depending on the values of the input
c/               parameters:
c/               ni_fixed, igradnih, igradniv, ni_lft, ni_rgt,
c/               ni_top, ni_btm for the ion density,
c/               te_fixed, igradteh, igradtev, te_lft, te_rgt, te_top,
c/               te_btm for the electron temperature and,
c/               ti_fixed, igradtih, igradtiv, ti_lft, ti_rgt, ti_top,
c/               ti_btm for the ion temperature.
c/               
c/    S_0      : Background volumetric source (/s)

c/    The following input variables designate and assign various 
c/    parameters to the top, bottom, left and right walls. It is 
c/    still possible to use the original input form if more "fine-
c/    tuning" is desired:
c/
c/    r_lft, r_rgt, r_top, r_btm : Reflection coefficients
c/    g_lft, g_rgt, g_top, g_btm : External incoming currents (#/s)
c/    flx_lft, flx_rgt, flx_top, flx_btm : ion fluxes (#/s)

      implicit none
      include 'neutGlob.inc'
      include 'comiou.inc'

c/    Local Variable Declarations:

      integer i, j, ij, k, iCell, iWall, iWall2, iw
      real delta_x, delta_y, ne_at_i, ne_at_j, ni_at_i, ni_at_j, 
     .     te_at_i, te_at_j, ti_at_i, ti_at_j, delta_ne, delta_ni,
     .     delta_te, delta_ti

c/    Calculate the number of internal cells, wall segments and plasma
c/    regions:
      
      nCells = NX * NY
      nWallSegm = 2* (NX + NY)
      nPlasmReg = 0

c/    Check Dimensions:
c/    ----------------
      if (nCells.GT.maxCell) then
	 write (6, '(1x, A30)') 'nCells > maxCell, Stopped!'
         call zstop(6,nout,2)
      endif
      if (nWallSegm.GT.maxWall) then
	 write (6, '(1x, A30)') 'nWallSegm > maxWall, Stopped!'
         call zstop(6,nout,2)
      endif

      delta_x = Lx / float(NX)
      delta_y = Ly / float(NY)

c/    Assign side lengths, angles, cell types etc.
c/    (Notice that side 1 is assumed to be the bottom side, and then
c/    we follow the  clock-wise direction.)

c/    Internal cells first:

      do iCell = 1, nCells
	 iType(iCell) = 0
	 nSides(iCell) = 4
	 do k = 1, 4
	    angle(k,iCell) = 90.0
         enddo
	 lside(1,iCell) = delta_x
	 lside(3,iCell) = delta_x
	 lside(2,iCell) = delta_y
	 lside(4,iCell) = delta_y
	 S_ext(iCell) = S_0
      enddo

c/    Assign background plasma parameters. For each of the four
c/    background plasma quantities (ne, ni, Te, Ti) we can have:
c/
c/    - Fixed values : ne_fixed, etc. positive and igradneh and
c/      igradnev equal to zero
c/    - Linear variation in the horizontal (igradneh = 1) or vertical
c/      (igradnev = 1) directions with input determined boundary 
c/      values (ne_lft, ne_rgt for igradneh = 1 or ne_btm, ne_top
c/      for igradnev = 1)
c/    - Manual input of the plasma parameters in each cell (igradnev = 
c/      igradneh = 0 and ne_fixed <=0)

      
c/    Electron density:

      if ((igradneh+igradnev).EQ.0) then    ! No gradients
	 if (ne_fixed.GT.0.0) then          ! Fixed background
	    do iCell = 1, nCells
	       elecDens(iCell) = ne_fixed
            enddo
         endif
      else if (igradneh.NE.0) then          ! Gradient in horizontal dir.
         delta_ne = (ne_rgt - ne_lft) / float(NX-1)
         if(iexp.eq.1)delta_ne=(ne_rgt/ ne_lft)**(1.0/float(NX-1))	 
	 do i = 1, NX
	    ne_at_i = ne_lft + (i-1) * delta_ne
            if(iexp.eq.1)ne_at_i=ne_lft*delta_ne**(i-1)
	    do j = 1, NY
	       ij = NY * (i-1) + j
	       elecDens(ij) = ne_at_i
            enddo
         enddo
      else if (igradnev.NE.0) then           ! Gradient in vertical dir.
	 delta_ne = (ne_top - ne_btm) / float(NY-1)
         if(iexp.eq.1)delta_ne=(ne_top/ ne_btm)**(1.0/float(NY-1))      
	 do j = 1, NY
            ne_at_j = ne_btm + (j-1) * delta_ne
            if(iexp.eq.1)ne_at_j=ne_btm*delta_ne**(j-1)
	    do i = 1, NX
	       ij = NY * (i-1) + j
	       elecDens(ij) = ne_at_j
            enddo
         enddo
      endif

c/    Ion density:

      if ((igradnih+igradniv).EQ.0) then    ! No gradients
	 if (ni_fixed.GT.0.0) then          ! Fixed background
	    do iCell = 1, nCells
	       ionDens(iCell) = ni_fixed
            enddo
         endif
      else if (igradnih.NE.0) then          ! Gradient in horizontal dir.
         delta_ni = (ni_rgt - ni_lft) / float(NX-1)	
         if(iexp.eq.1)delta_ni=(ni_rgt/ni_lft)**(1.0/float(NX-1)) 
	 do i = 1, NX
	    ni_at_i = ni_lft + (i-1) * delta_ni
            if(iexp.eq.1)ni_at_i = ni_lft* delta_ni**(i-1)
	    do j = 1, NY
	       ij = NY * (i-1) + j
	       ionDens(ij) = ni_at_i
            enddo
         enddo
      else if (igradniv.NE.0) then           ! Gradient in vertical dir.
	 delta_ni = (ni_top - ni_btm) / float(NY-1)
         if(iexp.eq.1)delta_ni=(ni_top/ni_btm)**(1.0/float(NY-1))
	 do j = 1, NY
            ni_at_j = ni_btm + (j-1) * delta_ni
            if(iexp.eq.1)ni_at_j=ni_btm *delta_ni**(j-1)
	    do i = 1, NX
	       ij = NY * (i-1) + j
	       ionDens(ij) = ni_at_j
            enddo
         enddo
      endif

c/    Electron Temperature:

      if ((igradteh+igradtev).EQ.0) then    ! No gradients
	 if (te_fixed.GT.0.0) then          ! Fixed background
	    do iCell = 1, nCells
	       elecTemp(iCell) = te_fixed
            enddo
         endif
      else if (igradteh.NE.0) then          ! Gradient in horizontal dir.
         delta_te = (te_rgt - te_lft) / float(NX-1)
         if(iexp.eq.1)delta_te=(te_rgt/te_lft)**(1.0/float(NX-1))	 
	 do i = 1, NX
	    te_at_i = te_lft + (i-1) * delta_te
            if(iexp.eq.1)te_at_i=te_lft*delta_te**(i-1)
	    do j = 1, NY
	       ij = NY * (i-1) + j
	       elecTemp(ij) = te_at_i
            enddo
         enddo
      else if (igradtiv.NE.0) then           ! Gradient in vertical dir.
	 delta_te = (te_top - te_btm) / float(NY-1)
         if(iexp.eq.1)delta_te=(te_top/te_btm)**(1.0/float(NY-1))
	 do j = 1, NY
            te_at_j = te_btm + (j-1) * delta_te
            if(iexp.eq.1)te_at_j=te_btm*delta_te**(j-1)
	    do i = 1, NX
	       ij = NY * (i-1) + j
	       elecTemp(ij) = te_at_j
            enddo
         enddo
      endif

c/    Ion Temperature:

      if ((igradtih+igradtiv).EQ.0) then    ! No gradients
	 if (ti_fixed.GT.0.0) then          ! Fixed background
	    do iCell = 1, nCells
	       ionTemp(iCell) = ti_fixed
            enddo
         endif
      else if (igradtih.NE.0) then          ! Gradient in horizontal dir.
         delta_ti = (ti_rgt - ti_lft) / float(NX-1)
         if(iexp.eq.1)delta_ti=(ti_rgt/ti_lft)**(1.0/float(NX-1))	 
	 do i = 1, NX
	    ti_at_i = ti_lft + (i-1) * delta_ti
            if(iexp.eq.1)ti_at_i=ti_lft*delta_ti**(i-1)
	    do j = 1, NY
	       ij = NY * (i-1) + j
	       ionTemp(ij) = ti_at_i
            enddo
         enddo
      else if (igradtiv.NE.0) then           ! Gradient in vertical dir.
	 delta_ti = (ti_top - ti_btm) / float(NY-1)
         if(iexp.eq.1)delta_ti=(ti_top/ti_btm)**(1.0/float(NY-1))
	 do j = 1, NY
            ti_at_j = ti_btm + (j-1) * delta_ti
            if(iexp.eq.1)ti_at_j=ti_btm*delta_ti**(j-1)
	    do i = 1, NX
	       ij = NY * (i-1) + j
	       ionTemp(ij) = ti_at_j
            enddo
         enddo
      endif

      if ((ialphaDen).EQ.1) then    ! Expoential
	  do i = 1, NX
          ni_at_i = 1+alpha*exp(-10*((i-0.5) / float(NX)-0.5)**2)	 
	    do j = 1, NY
	       ij = NY * (i-1) + j
	       ionDens(ij) = ni_at_i*ni_fixed
	       elecDens(ij) = ni_at_i*ne_fixed
            enddo
         enddo
      endif

      if ((ialphaAll).EQ.1) then    ! Expoential
	  do i = 1, NX
             ni_at_i =(((i-0.5)/NX-0.5)*Lx)**2
	    do j = 1, NY
             ni_at_j =(((j-0.5)/NY-0.5)*Ly)**2
             ni_at_j= exp(alpha*(ni_at_i+ni_at_j))   	 
	       ij = NY * (i-1) + j
	       ionDens(ij) = ni_at_j*ni_fixed
	       elecDens(ij) = ni_at_j*ne_fixed
	       ionTemp(ij) = ni_at_j*ti_fixed
	       elecTemp(ij) = ni_at_j*te_fixed
            enddo
         enddo
      endif


      do iWall = nCells + 1, nCells + nWallSegm
	 iType(iWall) = 2
	 nSides(iWall) = 1
      enddo

c/    Now, assign neighbors. Start first with the cells that are bounded
c/    by the first wall. Notice, that the convention that we have adopted
c/    here is that cell 1 is the bottom leftmost cell. Then, we move in
c/    columns from bottom to top and left to right. The first wall segment
c/    is the left side of the first cell, and then we move clockwise,
c/    i.e, left boundary (LB) to Top, to Right Boundary (RB) to Bottom.

c/    Left Boundary:
c/    -------------
      iWall = nCells
      do j = 1, NY
         iWall = iWall + 1
         iCell = j
         adjCell(1,iWall) = iCell
         adjCell(2,iCell) = iWall
         adjCell(4,iCell) = iCell + NY
         if (j.EQ.1) then
            adjCell(1,iCell) = nCells + nWallSegm
            adjCell(1,nCells+nWallSegm) = iCell
         else
            adjCell(1,iCell) = iCell - 1
         endif
         if (j.EQ.NY) then
            adjCell(3,iCell) = nCells + NY + 1
            adjCell(1,nCells+NY+1) = iCell
         else
            adjCell(3,iCell) = iCell + 1
         endif
      enddo

c/    Top Boundary (except cells at the edges):
c/    ----------------------------------------
      iWall = iWall + 1
      do i = 2, NX - 1
         iWall = iWall + 1
         iCell = i * NY
         adjCell(1,iWall) = iCell
         adjCell(3,iCell) = iWall
         adjCell(2,iCell) = (i-1) * NY
         adjCell(4,iCell) = (i+1) * NY
         adjCell(1,iCell) = iCell - 1
      enddo

c/    Right Boundary (scan from top to bottom):
c/    ----------------------------------------
      iWall = iWall + 1
      do j = 1, NY
         iWall = iWall + 1
         iCell = nCells - j + 1
         adjCell(1,iWall) = iCell
         adjCell(4,iCell) = iWall
         if(NX.gt.1)adjCell(2,iCell) = iCell - NY
         if (j.EQ.1) then
            iWall2 = nCells + NY + NX
            adjCell(3,iCell) = iWall2
            adjCell(1,iWall2) = iCell
         else
            adjCell(3,iCell) = iCell + 1
         endif
         if (j.EQ.NY) then
            iWall2 = nCells + NY + NX + NY + 1
            adjCell(1,iCell) = iWall2
            adjCell(1,iWall2) = iCell
         else
            adjCell(1,iCell) = iCell - 1
         endif
      enddo

c/    Bottom boundary (scan from right to left, except edge cells):
c/    ------------------------------------------------------------
      iWall = iWall + 1
      do i = 2, NX - 1
         iWall = iWall + 1
         iCell = (NX-1)*NY + 1 - (i-1) * NY
         adjCell(1,iWall) = iCell
         adjCell(1,iCell) = iWall
         adjCell(2,iCell) = iCell - NY
         if(NY.gt.1)adjcell(3,iCell) = iCell + 1
         adjCell(4,iCell) = iCell + NY
      enddo

c/    Now we do the internal cells:
c/    Notice, that iCell(i,j) = (i-1)*NY + j, in general)
c/    ---------------------------------------------------
      do i = 2, NX - 1
         do j = 2, NY - 1
            iCell = (i-1) * NY + j
            adjCell(1,iCell) = iCell - 1
            adjCell(3,iCell) = iCell + 1
            adjCell(2,iCell) = iCell - NY
            adjCell(4,iCell) = iCell + NY
         enddo
      enddo

c/    We now assign various coefficients and fluxes on the walls. 
c/    We also allow for original input assignement for these
c/    quantities:

c/    Left wall:
c/    ----------
      do iw = 1, NY
	 Rwall(iw) = Rwall(iw) + r_lft
	 g_ex(iw) = g_ex(iw) + g_lft
	 g_ion(iw) = g_ion(iw) + flx_lft
      enddo

c/    Top boundary:
c/    ------------
      do iw = NY + 1, NY + NX
	 Rwall(iw) = Rwall(iw) + r_top
	 g_ex(iw) = g_ex(iw) + g_top
	 g_ion(iw) = g_ion(iw) + flx_top
      enddo

c/    Right Wall:
c/    ----------
      do iw = NY+NX+1, NX + 2*NY
	 Rwall(iw) = Rwall(iw) + r_rgt
	 g_ex(iw) = g_ex(iw) + g_rgt
	 g_ion(iw) = g_ion(iw) + flx_rgt
      enddo

c/    Bottom Wall:
c/    -----------
      do iw = NX + 2*NY + 1, 2*(NX+NY)
	 Rwall(iw) = Rwall(iw) + r_btm
	 g_ex(iw) = g_ex(iw) + g_btm
	 g_ion(iw) = g_ion(iw) + flx_btm
      enddo

      return
      end
