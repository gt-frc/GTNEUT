c/////////////////////// File locGeom.inc //////////////////////////////

c/    This file contains geometric variables needed in the calculation
c/    of the various transmission coefficients.

      integer i_geom, nph_pnts, lcosn, inon1, idp1
      real L_i, L_j, L_perp, l_mfp, theta_ij, alpha_j, R1, R2,
     .   dTheta, L_rad
      common /rcgeom/i_geom, nph_pnts, L_i, L_j, L_perp, l_mfp, 
     .   lcosn, theta_ij, alpha_j, R1, R2, dTheta, L_rad,
     .   idp1, inon1
