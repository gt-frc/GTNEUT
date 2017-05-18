c/                       G T N E U T

c///////////////////////////////////////////////////////////////////////
c/                                                                     /
c/    This is a 2-D neutral transport code based on the TEP            /
c/    (transmission/escape probability) method.                        /
c/                                                                     /
c/    Reference: W.M. Stacey, J. Mandrekas, "A Transmission/Escape     /
c/    Probabilities Model for Neutral Particle Transport in the Outer  /
c/    Regions of a Diverted Tokamak," Nuclear Fusion, 34 (1994) 1385.  /
c/                                                                     /
c/    04/01/99, jm: New 1999 version with only straight-line geometry  /
c/    and corrected expressions for the transmission coefficients.     /
c/                                                                     /
c/    07/28/99, jm: Several changes prompted by the interface with     /
c/    the UEDGE code:                                                  /
c/                                                                     /
c/     - Removed iTypes 3 and 4. Now, all material interfaces are of   /
c/       iType 2.                                                      /
c/                                                                     /
c/     - Changed the index of several wall related quantities (g_ex,   /
c/       g_ion, Rwall) from absolute to relative, i.e. index           /
c/       1 is now the first wall segment. **NOTICE**: This WILL break  /
c/       old input files! The relation between old absolute indices    /
c/       and new relative ones is: irel = iabs - (nCells + nPlasmReg)  /
c/                                                                     /
c/     - Changed the 3D array transm(j,k,i) to represent the           /
c/       transmission coefficient from side j to side k of cell i,     /
c/       in order to drastically reduce the storage requirements.      /
c/                                                                     /
c/    05/18/00, jm, introduced option for automatic input generation   /
c/    for rectangular regions (i_inpt = 1)                             /
c/                                                                     /
c/    11/06/00, jm, introduced scale factor to multiply lenghts for    /
c/    manually prepapred input files                                   /
c/                                                                     /
c/    03/13/01, jm, introduced option to use DEGAS atomic rates        /
c/                                                                     /
c/    10/08/01, jm, replaced cross section evaluation with Janev's     /
c/    fits.                                                            /
c/                                                                     /
c/    12/14/01, jm, changes in several arrays (npos, gflux, etc.) so   /
c/    that now the second index refers to the side of the cell rather  /
c/    than the neighboring cell, to reduce memory requirements.        /
c/                                                                     /
c/    12/18/01, jm, new reflection model (irefl = 1)                   /
c/                                                                     /
c/    01/14/02, jm, eliminated core plasma and wall fluxes from        /
c/                  solution vector, to conserve memory and to         /
c/                  make implementation of energy-dependent wall       /
c/                  reflection easier. Wall and core plasma fluxes     /
c/                  are now constructed from internal cell fluxes      /
c/                                                                     /
c/    06/09/03, jm, moved details of transmission coefficients calc.   /
c/                  into a separate routine (calctransm)               /
c/                                                                     /
c/    06/10/03, jm, Sparse matrix version using the UMFPACK 2.2.1      /
c/                  library                                            /
c/                                                                     /
c/    02/03/09,zwf, added new gauss quadrature set n=10
c/    02/25/09,zwf, eneut_v no longer used when using the ANE
c/                  instead, the local ionTemp is used. This is better /
c/                  for modeling a recombination source                /
c///////////////////////////////////////////////////////////////////////
c/
c/    INPUT VARIABLES:
c/    ---------------
c/    i_inp    : Flag, determining the input options.
c/               = 0, Original input format
c/               = 1, Automatic input generation for rectangular regions.
c/
c/    Input variables needed when i_inp = 0:
c/    -------------------------------------
c/    Notice that in numbering the various regions, the internal cells
c/    come first, then any plasma regions and finally the wall segments!
c/
c/    nCells       : Number of (internal) regions (or cells)
c/    nPlasmReg    : Number of plasma regions
c/    nWallSegm    : Number of wall segments
c/    iType(i)     : Type of region -i- . Options are:
c/                   0 = regular (internal) cell
c/                   1 = plasma
c/                   2 = wall
c/    nSides(i)    : Number of sides for cell -i-
c/    lside(k,i)   : Length of side -k- for rect. cell -i- [m]
c/    angle(k,i)   : Angle between k and k+1 sides for rect. cell -i- [deg]
c/    adjCell(k,i) : Index of adjacent cell to side -k- of cell -i-
c/    scalFact     : Scale factor to multiply lengths (equal to 1 by default)
c/
c/    The following plasma parameters must be defined for the internal
c/    regions (cells) AND for the plasma regions, in order to calculate
c/    the albedo coefficients of the plasma regions:
c/
c/    elecTemp(i)  : Electron temperature in cell -i- [keV]
c/    ionTemp(i)   : Ion temperature in cell -i- [keV]
c/    elecDens(i)  : Electron density in cell -i- [m^-3]
c/    ionDens(i)   : Ion density in cell -i- [m^-3]
c/
c/    S_ext(i)     : External volumetric neutral source in cell i, [#/s]
c/
c/    The index of the following input variables starts at the first wall
c/    segment:
c/
c/    g_ex(iw)     : External surface flux at wall segment -iw-  [#/s]
c/                   (neutrals are supposed to be launched with energy
c/                   eneut)
c/    Rwall(iw)    : Reflection coefficient of wall segment -iw-
c/    g_ion(iw)    : Ion flux on wall segment -iw- (usually provided
c/                   by a plasma code), [#/s]
c/    fwabsorb(iw) : Wall absorption coefficient (if irefl = 1)
c/    awall(iw)    : Atomic mass of wall material (irefl = 1)
c/    zwall(iw)    : Atomic number of wall material (irefl = 1)
c/    twall(iw)    : Wall temperature (keV) of wall (irefl = 1)
c/
c/    Input variables needed when i_inp = 1 (rectangular regions):
c/    -----------------------------------------------------------
c/    Lx           : Length of x-side of rectangle (m)
c/    Ly           : Length of y-side of rectangle (m)
c/    NX           : Number of cells in the x-direction
c/    NY           : Number of cells in the y-direction
c/    ne_fixed     : Background plasma electron density (/m^3) if no
c/                   gradients are present. If no gradients are present
c/                   and if ne_fixed <= 0, then the n_e for each cell is
c/                   read from the input file.
c/    igradneh     : If equal to 1, linear variation in the electron
c/                   density in the horizontal direction (X) determined
c/                   by the values of ne_lft and ne_rgt (see below)
c/    igradnev     : If equal to 1, linear variation in the electron
c/                   density in the vertical direction (Y) determined
c/                   by the values of ne_btm and ne_top (see below).
c/    ne_lft       : Electron density at the left boundary (/m3)
c/    ne_rgt       : Electron density at the right boundary (/m3)
c/    ne_top       : Electron density at the top boundary (/m3)
c/    ne_btm       : Electron density at the bottom boundary (/m3)
c/                   The same holds for the other background plasma
c/                   quantities which are assigned depending on the
c/                   values of the input parameters:
c/                   ni_fixed, igradnih, igradniv, ni_lft, ni_rgt,
c/                   ni_top, ni_btm for the ion density,
c/                   te_fixed, igradteh, igradtev, te_lft, te_rgt,
c/                   te_top, te_btm for the electron temperature and,
c/                   ti_fixed, igradtih, igradtiv, ti_lft, ti_rgt, ti_top,
c/                   ti_btm for the ion temperature.
c/
c/    te_fixed     : Background plasma electron temperature (keV)
c/    ti_fixed     : Background plasma ion temperature (keV)
c/    ne_fixed     : Background plasma electron density (/m3)
c/    ni_fixed     : Background plasma ion density (/m3)
c/    S_0          : Background volumetric source (/s)

c/    The following input variables designate and assign various
c/    parameters to the top, bottom, left and right walls. It is
c/    still possible to use the original input form if more "fine-
c/    tuning" is desired:
c/
c/    r_lft, r_rgt, r_top, r_btm : Reflection coefficients
c/    g_lft, g_rgt, g_top, g_btm : External incoming currents (#/s)
c/    flx_lft, flx_rgt, flx_top, flx_btm : ion fluxes (#/s)
c/
c/    Common input for all cases:
c/    --------------------------
c/    zion         : Atomic number of ion species
c/    aion         : Mass number of ion species
c/    aneut        : Mass number of neutral (1,2 or 3)
c/    i_e0         : Flag for the calculation of the neutral velocity:
c/                   if equal to 1, constant E0 = eneut
c/                   if equal to 2, use local ion temperature (E0 = T_i)
c/    icosn        : If equal to 1, assume a cosine distribution function at
c/                   the interfaces. If equal to 0, assume a uniform
c/                   distribution function
c/    v0fact       : v0 = Sqrt(v0fact*E0 / aneut)
c/    eneut        : Neutral energy for g_ex neutrals or for constant energy
c/                   option [keV]
c/    eneut_v      : Energy of volumetric neutrals [keV]
c/                  (If the source is due to recombination set ifrstcol = 0.)
c/    prntOrdr     : Integer array affecting the order of printout
c/    iatdat       : Flag determining the atomic reaction rates to use:
c/                   0 : Original (Janev) formulation
c/                   1 : DEGAS reaction rates
c/                   2 : Thomas/Stacey rates
c/    ifjsv        : If > 0, use Freeman and Jones fits for electron impact
c/                   ionization reaction rates (for iatdat = 0 case)
c/    irefl        : Flag for reflection model. If equal to 0, then
c/                   wall reflection is controlled by the input array
c/                   Rwall. If irefl = 1, then the reflection coefficient
c/                   is calculated using fits that depend on the projectile
c/                   and target properties. Rwall can still be used (e.g.,
c/                   to model vacuum regions) for wall segments with
c/                   negative zwall. See the routine calcrefln for more
c/                   details.
c/    ifrstcol     : If equal to 1, treat first collision and volumetric
c/                   neutrals separately
c/
c/    The following two variables are only relevant if iatdat = 1:
c/    leh0         : 1: e + h0 ionization rate dependent on ne
c/                   2: e + h0 ionization rate independent of ne
c/    lchex        : if <=2, CX cross section depends on neutral energy
c/                 if = 3, CX cross section independent of neutral energy.
c/
c/    iquad        : index determining the number of integration points
c/                   used in some of the transmission coefficient
c/                   integrals. Current choices are:
c/                   iquad = 1-5 for 20 - 100 points.
c/                   new quadrature set found from:
c/                   http://www.holoborodko.com/pavel/?page_id=679
c/                   Quadrature Library is under LGPL license
c/                   zwf added new quadrature set for n=10
c/                   iquad = 0 for 10 points
c/                   iquad = 55 for 6 points (probably not useful)
c/
c/    nph          : Number of points for angular integration
c/    iescp        : Flag determining how we calculate the escape
c/                   probability, P0.
c/                   = 0 : Original Wigner form
c/                   = 1 : Sauer-like approximation with optimized exponent
c/    isparsitr    : Number of steps for iterative improvement of Sparse
c/                   system solution (=0: no iterative improvement)
c/
c/    CONFUSION ALERT:
c/    ---------------
c/    Notice that all the quantities in GTNEUT (volume sources, surface
c/    sources, ionization rates, neutrals) are TOTAL quantities in each
c/    region per unit toroidal length. To get densities, we have to
c/    divide each of these quantities by the area of each cell.

c/    Notice convention on numbering sides of cells:
c/    ----------------------------------------------
c/    Straightline polygons: 1 is the bottom side (if more than one,
c/    take the one facing left), the rest clockwise. The plasma region
c/    can consist of sub-regions or can be a single multi-sided region.

c/    EXIT CODES:
c/    ----------
c/    0  : Normal exit
c/    1  : Inconsistent input parameters
c/    2  : Need to increase dimensions in neutGlob.inc

      implicit none
      include 'neutGlob.inc'
      include 'comiou.inc'
      include 'consts.inc'

	integer i,j,converge

c/    Input namelists and IO:
c/    ----------------------
      namelist /inp/ i_inp, nCells, nWallSegm, nPlasmReg, iType,
     .   nSides, lside, angle, adjCell, scalFact, idbug,
     .   elecTemp, ionTemp, elecDens, ionDens, S_ext, Rwall, g_ion,
     .   i_e0, eneut, eneut_v, zion, aion, aneut, g_ex, prntOrdr,
     .   iatdat, lchex, leh0, v0fact, icosn, iquad, nph,
     .   iescp, irefl, ifjsv, twall, fwabsorb, zwall, awall, ifrstcol,
     .   isparsitr,idp, inon, nd0, neitr , Shotnumber, Timeslice,nxleg1,
     .   nxcore1,nxcore2,nycore1,nysol1,nxxpt,nxmod,nxleg2


      namelist /inp1/ NX, NY, Lx, Ly, te_fixed, ti_fixed, ne_fixed,
     .   ni_fixed, S_0, r_lft, r_rgt, r_top, r_btm, g_lft, g_rgt, g_top,
     .   g_btm, flx_lft, flx_rgt, flx_top, flx_btm, igradneh, igradnev,
     .   igradnih, igradniv, igradteh, igradtev, igradtih, igradtiv,
     .   ne_lft, ne_rgt, ne_top, ne_btm, ni_lft, ni_rgt, ni_top, ni_btm,
     .   te_lft, te_rgt, te_top, te_btm, ti_lft, ti_rgt, ti_top, ti_btm,
     .   iexp, ialphaAll, ialphaDen, alpha

      open (nin, file = 'toneut', status = 'old')

      read (nin,inp)

      if (i_inp.EQ.1) read (nin, inp1)

      open (nout, file = 'neut.out', status = 'unknown')
      open (tomat, file = 'nmat.m', status = 'unknown')

      if (idbug.eq.1) then
        open (ndbug, file = 'neut.dbg', status = 'unknown')
      endif

c/    Open file for error and diagnostic messages from the UMFPACK
c/    library:

      open (nsdbug, file = 'umferr.dat', status = 'unknown')

c/    Read the DEGAS atomic data files if iatdat = 1:

      if (iatdat.EQ.1) then
	 open (ioeh, file = 'ehr1.dat', status = 'old')
	 open (iocx, file = 'cxionh.dat', status = 'old')
         call degasread
	 close (ioeh)
	 close (iocx)
      endif

      if (i_inp.EQ.1) then
	 call rectinp
	 if (idbug.EQ.1) call output(0)
      endif

c/    START CALCULATIONS:
c/    ------------------
c/    Total number of regions:

      nTotal = nCells + nWallSegm + nPlasmReg

c/    Check if input data are consistent:

      call checkInput
      print *, 'neitr = ', neitr

c/    START CONVERGENCE LOOP
      do j=1, neitr
c/          START LOOP OVER SPECIES
c/            do spec=1,2
c/              Call the wall reflection model:

                call calcrefln

c/              Calculate the neutral mean free path in each region:

                call calcmfp
c/                if (spec.EQ.1) then
c/                    call calcmfpad
c/                elseif (spec.EQ.2) then
c/                    call calcmfpmd
c/                endif

c/              Calculate first-flight transmission coefficients:

                call calctransm

c/              Calculate the Escape probabilities:

                call escape

c/              FEM
                if(nd0.gt.0)then
                    do i=1, nCells
                        call fem(i,nd0)
                    enddo
                endif

c/              Set up the solution matrix:

                call setup

c/              Solve the linear system of Equations:

                call solvers
                print *, 'value of j is', j
                do i=1,nCells
                    print *,'i, neutdens = ',i,neutdens(i)
                enddo

                print *,''
c/          end species loop
c/            enddo
c/    end convergence loop
      enddo

c/    Calculate global particle balance:
      call pbalance
c/    END OF LOOP OVER SPECIES

c/    CONVERGENCE CHECK
c/    end do
c/    END OF CONVERGENCE LOOP

c/    Write output:
      call output(1)



      call zstop (6, nout, 0)


      end
