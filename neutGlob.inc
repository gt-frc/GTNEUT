c/////////////////// NEUTGLOB.INC //////////////////////////////////////
c/    The big COMMON BLOCK with the declarations of the global
c/    variables of the code.

      integer maxCell, maxWall, maxPlas, maxCPl, maxSides, maxTot,
     .   maxEqs, maxElmns, maxExp

      parameter (maxCell = 60000,
     .           maxWall = 15000,
     .           maxPlas = 15000,
     .           maxSides = 5,
     .           maxExp = 4,
     .           maxTot = maxCell + maxPlas + maxWall,
     .           maxCPl = maxCell + maxPlas,
     .           maxEqs = maxCell * maxSides*maxExp,
     .           maxElmns = maxEqs * (1+maxSides*maxExp))

      integer iType(maxTot), nSides(maxTot), adjCell(maxSides, maxTot)

      real lside(maxSides, maxCell), angle(maxSides,maxCell),scalFact,
     .    tneut(maxCPl,maxSides)

      common /neut01/ iType, nSides, lside, angle, adjCell, scalFact,
     .    tneut


      integer i_inp, nCells, nWallSegm, nPlasmReg, nEqs, nTotal, iquad,
     .   nph, icosn, iescp, ifrstcol,idp,nExp,inon,mExp, nd0, neitr,
     .   Shotnumber, Timeslice,nxleg1,
     .   nxcore1,nxcore2,nycore1,nysol1,nxxpt,nxmod,nxleg2,spec

      real transm(maxSides, maxSides, maxExp, maxExp, maxCPl),
     .  transw_f(maxSides,maxExp,maxExp,maxWall),
     .  transw_s(maxSides,maxExp,maxExp,maxWall),
     .  transw_0(maxSides,maxExp,maxExp,maxWall),
     .  stransw_f(maxWall,maxExp),
     .  stransw_s(maxWall,maxExp), stransw_0(maxWall,maxExp),
     .  colterm(maxSides,maxExp,maxCPl)

      common /neut02/ i_inp, nCells, nWallSegm, nPlasmReg, nTotal, nEqs,
     .   iquad, nph, icosn, iescp, ifrstcol,transm,transw_f,idp,nExp,
     .   transw_s, transw_0, stransw_f, stransw_s, stransw_0,inon,mExp,
     .   nd0,colterm, neitr,
     .   Shotnumber, Timeslice,nxleg1,
     .   nxcore1,nxcore2,nycore1,nysol1,nxxpt,nxmod,nxleg2,spec

c/    Sources, wall & divertor properties, etc:

      integer irefl
      real S_ext(maxCell), Rwall(maxWall), g_ion(maxWall),
     .  g_ex(maxWall), rcycle(maxWall), twall(maxWall), awall(maxWall),
     .  zwall(maxWall), fwabsorb(maxWall), refln(maxWall),
     .  refle(maxWall)
      common /neut03/ irefl, S_ext, Rwall, g_ion, g_ex, rcycle, twall,
     .   awall, zwall, fwabsorb, refln, refle

c/    Variables related to printing:
      integer prntOrdr(maxTot), idbug
      common /comprn/ prntOrdr, idbug

c/    Some calculated geometric parameters:

      real area(maxCell), perim(maxCell)
      common /neut04/ area, perim

c/    Plasma parameters:

      integer i_e0
      real elecTemp(maxCPl), ionTemp(maxCPl), elecDens(maxCPl),
     .   ionDens(maxCPl), eneut, eneut_v, zion, aion, aneut, v0fact

      common /plasma/ i_e0, elecTemp, ionTemp, elecDens, ionDens, eneut,
     .   eneut_v, zion, aion, aneut, v0fact

c/    Atomic physics parameters (cross sections, etc)

      integer iatdat, ifjsv, lchex, leh0
      real svion_e(maxCPl), svion_i(maxCPl), svion_tot(maxCPl),
     .   sv_cx(maxCPl), lmfp(maxCPl), A_cx(maxCPl), albedo(maxPlas),
     .   lmfp0(maxCPl), mfp(maxSides, maxCPl), A_cxk(maxSides, maxCPl),
     .   mfp_w0(maxWall), A_cx0(maxCPl), mfp_wf(maxWall),
     .   mfp_ws(maxWall), A_cxw0(maxWall), A_cxwf(maxWall),
     .   A_cxws(maxWall)

      common /xsect/ svion_e, svion_i, svion_tot, sv_cx, lmfp, A_cx,
     .   albedo, lmfp0, A_cx0, mfp, A_cxk, mfp_w0, mfp_wf, mfp_ws,
     .   A_cxw0, A_cxwf, A_cxws, iatdat, lchex, leh0, ifjsv

c/    Escape - probability quantities:

      real pEscp(maxCell), pEscp0(maxCell), pEscp0s(maxCell),
     .     pEscpw0(maxWall), pEscpwf(maxWall), pEscpws(maxWall),
     .     lambda(maxTot,maxSides), pEscp0k(maxSides,maxCell),
     .     lambdak(maxSides,maxCell,maxSides,5),
     .     pEscpk(maxSides,maxCell), lambdawf(maxWall,maxSides,5),
     .     lambdaws(maxWall,maxSides,5),
     .     lambdaw0(maxWall,maxSides,5)

      common /escap/ pEscp, pEscp0, pEscp0s, pEscpw0, pEscpwf, pEscpws,
     .    lambda, pEscp0k, lambdak, pEscpk, lambdaw0,lambdawf,lambdaws

c/    Variables related to the matrix inversion:

      integer npos(maxCell, maxSides,maxExp)
      integer nElmnts, isparsitr, i_sparse(2*maxElmns)
      real b(maxEqs), a_sparse(maxElmns)
      common /smatrx/ nElmnts, isparsitr, npos, b, i_sparse, a_sparse

c/    Fluxes, ionization rates, etc:

      real srcNormFact, gflux(maxTot, maxSides,maxExp),
     .     gflux_u(maxTot,maxSides), gflux_c(maxTot, maxSides),
     .     gflux_c1(maxTot,maxSides), gflux_cn(maxTot, maxSides),
     .     inzRate(maxCell), neutDens(maxCell)

      common /fluxes/ srcNormFact, gflux, gflux_c, gflux_c1, gflux_cn,
     .   gflux_u, inzRate, neutDens

c/    Global particle balance parameters:

      real fluxn_tot, fluxi_tot, volsrc_tot, coreloss_tot, wloss_tot,
     .   ionloss_tot, tot_part_in, tot_part_out, partbalnerr

      common /pbaln/ fluxn_tot, fluxi_tot, volsrc_tot, coreloss_tot,
     .       wloss_tot, ionloss_tot, tot_part_in, tot_part_out,
     .       partbalnerr

c/    The following COMMON block contains the input variables needed
c/    for the input generation for rectangular regions.

      integer NX, NY, igradteh, igradtev, igradtih, igradtiv, igradneh,
     .   igradnev, igradnih, igradniv, iexp, ialphaDen,ialphaAll
      real Lx, Ly, te_fixed, ti_fixed, ne_fixed, ni_fixed, S_0, r_lft,
     .   r_rgt, r_top, r_btm, g_lft, g_rgt, g_top, g_btm, flx_lft,
     .   flx_rgt, flx_top, flx_btm, rcy_lft, rcy_rgt, rcy_top, rcy_btm,
     .   te_lft, te_rgt, te_top, te_btm, ti_lft, ti_rgt, ti_top, ti_btm,
     .   ne_lft, ne_rgt, ne_top, ne_btm, ni_lft, ni_rgt, ni_top, ni_btm,
     .   alpha


      common /rectangl/ NX, NY, igradteh, igradtev, igradtih, igradtiv,
     .   igradneh, igradnev, igradnih, igradniv, Lx, Ly, te_fixed,
     .   ti_fixed, ne_fixed, ni_fixed, S_0, r_lft, r_rgt, r_top, r_btm,
     .   g_lft, g_rgt, g_top, g_btm, flx_lft, flx_rgt, flx_top, flx_btm,
     .   rcy_lft, rcy_rgt, rcy_top, rcy_btm, te_lft, te_rgt, te_top,
     .   te_btm, ti_lft, ti_rgt, ti_top, ti_btm, ne_lft, ne_rgt, ne_top,
     .   ne_btm, ni_lft, ni_rgt, ni_top, ni_btm, iexp,ialphaDen,
     .   ialphaAll, alpha

