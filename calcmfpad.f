      subroutine calcmfpad

c/    The purpose of this routine is to calculate the neutral mean free
c/    path and the charge exchange fraction in each cell, as well as the
c/    albedo for the plasma regions. Two different mean free paths are
c/    calculated: A local one (lmfp(i)) which is used for the
c/    calculation of escape probabilities and for the constant energy
c/    case, and another one which is calculated using the velocity of
c/    the neutrals from the side from which they entered (mfp(kk,i)).
c/
c/    The routine includes options for using various cross section
c/    libraries, and for whether to take into account first-collision
c/    effects.


c/    Record of changes:
c/    -----------------
c/    08/14/97, jm: Option for E.W. Thomas cross sections
c/    06/05/98, jm: Options for the calculation of neutral velocity
c/    03/13/01, jm: Option for using DEGAS atomic rates (iatdat = 1)
c/    10/08/01, jm: Replaced original cross section routines with
c/                  the latest fits from Janev's book.
c/    01/07/02, jm: Significant changes to calculate the mean free
c/                  path taking into account the correct neutral
c/                  velocity.
c/    06/03/03, jm: Rearranged parts of the routine for easier
c/                  understanding. Removed old -imfptr- option.
c/    06/03/03, jm: The energy of volumetric neutrals is now
c/                  eneut_v. The old input variable -eneut- is
c/                  reserved for the external flux neutrals, as
c/                  well as for the constant energy case (i_e0 = 1).

      implicit none
      include 'neutGlob.inc'
      include 'locGeom.inc'
      include 'consts.inc'

      integer i, j, k, kk, kw, jj
      real v0, dsvione, dsvioni, dsvcx, invmfp, afactr, teev, tiev,
     .     nem3, tnev, svion, svcx, svel, sveln, svrec, E_0, svione,
     .     svioni, svcxi, svefj, d_ratio, v0k, v0f, v0s, svion_k,
     .     sv_cx_k, svion_tot_k, dsvcxk, dsvionik, svcxk, svreck,
     .     svelk, svelnk, invmfpk, e_slow, e_fast, tnfev, tnsev, tn0ev,
     .     svion_w, sv_cx_w0, sv_cx_wf, sv_cx_ws, svion_tot_w0,
     .     svion_tot_ws, svion_tot_wf, dsvcxw0, dsvcxwf, dsvcxws,
     .     dsvioniw0, dsvioniwf, dsvioniws, svcxw0, svcxwf, svcxws,
     .     svrecw0, svrecwf, svrecws, svelw0, svelwf, svelws, svelnw0,
     .     svelnwf, svelnws, svion_w0, svion_wf, svion_ws, invmfpwf,
     .     invmfpws, invmfpw0, E_s0, tnev0, svion_i0, sv_cx0,
     .     dsvione0, dsvcx0, dsvioni0, svion0, svrec0,
     .     svcx0, svel0, sveln0, vs0, invmfp0, svion_tot0, albdfit_ad


c/    Main DO loop over internal cells and plasma regions:

      do i = 1, nCells + nPlasmReg

c/    We first calculate the neutral mean-free-path and the number
c/    of secondary neutrals per collision (charge exchange fraction)
c/    using the local parameters in each cell, including the neutral
c/    energy. If the i_e0 flag is equal to 1 (which means that we
c/    do not have the wall reflection model) then the neutral energy
c/    is equal to the input variable eneut. Otherwise, eneut is set
c/    equal to the local ion temperature. If the reflection model is
c/    not on, then this is the only calculation that we need since we
c/    do not have different energy groups.

         if (i_e0.EQ.1) then
	    E_0 = eneut
	 else if (i_e0.EQ.2) then
	    E_0 = ionTemp(i)
         else if (i_e0.GE.3) then
            E_0 = 1.5*ionTemp(i)
         endif

	 teev = 1000.0 * elecTemp(i)
	 tiev = 1000.0 * ionTemp(i)
	 tnev = 1000.0 * E_0
	 nem3 = elecDens(i)
	 d_ratio = nem3 / ionDens(i)

c/    Calculate rates using appropriate model:

	 if (iatdat.EQ.0) then                    ! Janev

	   if (ifjsv.EQ.0) then
	      svion_e(i) = svione(teev)
	   else
	      svion_e(i) = svefj(teev)            ! Freeman-Jones
           endif

           svion_i(i) = svioni(aion, tiev, aneut, tnev)
           sv_cx(i) = svcxi(aion, tiev, aneut, tnev)
	   svion_tot(i) = svion_i(i) + d_ratio * svion_e(i)

	 else if (iatdat.EQ.1) then                ! DEGAS

	   call svdegas (teev, tiev, tnev, nem3, aion, aneut,
     .	      leh0, lchex, dsvione, dsvcx, dsvioni)
	   svion_e(i) = dsvione
	   svion_i(i) = dsvioni
	   sv_cx(i) = dsvcx
	   svion_tot(i) = svion_i(i) + d_ratio * svion_e(i)

	 else if (iatdat.EQ.2) then               ! Thomas / Stacey

	   call calcxswms(teev, tiev, tnev, nem3, svion, svcx, svrec,
     .       svel, sveln)
	   svion_e(i) = svion
	   sv_cx(i) = svcx + svel
	   svion_tot(i) = d_ratio * svion_e(i)

	 endif

         v0 = sqrt(v0fact * E_0 * xj7kv / (protMass * aneut))
	 invmfp = ionDens(i) * (svion_tot(i) + sv_cx(i)) / v0
	 lmfp(i) = 1.0 / invmfp
	 A_cx(i) = sv_cx(i) / (sv_cx(i) + svion_tot(i))

c/    If we are running a case with constant neutral energy (usually
c/    for diagnostic purposes), then we set all other mean-free-paths
c/    and charge exchange fractions equal to the local ones:

         if (i_e0.EQ.1) then

            lmfp0(i) = lmfp(i)   ! for volumetric neutrals
            A_cx0(i) = A_cx(i)

            do kk = 1, nSides(i)
	       k = adjCell(kk,i)
	       if (iType(k).NE.2) then
	          mfp(kk,i) = lmfp(i)
	          A_cxk(kk,i) = A_cx(i)
	       else if (iType(k).EQ.2) then
		  kw = k - (nCells + nPlasmReg)
		  mfp_w0(kw) = lmfp(i)
		  mfp_wf(kw) = lmfp(i)
		  mfp_ws(kw) = lmfp(i)
		  A_cxw0(kw) = A_cx(i)
		  A_cxwf(kw) = A_cx(i)
		  A_cxws(kw) = A_cx(i)
	       endif

            enddo    ! end of do loop over sides of -i-

c/    For all other cases, we calculate the mean-free-paths and
c/    charge exchange fractions taking into account the proper
c/    neutral energy:

         else

            do kk = 1, nSides(i)
	       k = adjCell(kk,i)
	       do jj=1,nSides(k)
	          if(adjCell(jj,k).eq.i)j=jj
	       enddo
c/    Adjacent cell is NOT a wall segment:

	       if (iType(k).NE.2) then
	         tnev = 1000.0 * tneut(k,j)
                 if(i_e0.eq.3)tnev=1.5*tnev
                 v0k = sqrt(v0fact * tnev * xj7ev / (protMass * aneut))

	         if (iatdat.EQ.0) then               ! Janev

c/    Since ion-impact ionization cross section is very small, use
c/    the already calculated value:

                   svion_k = svion_i(i)
                   sv_cx_k = svcxi(aion, tiev, aneut, tnev)
	           svion_tot_k = svion_k + d_ratio * svion_e(i)

	         else if (iatdat.EQ.1) then          ! DEGAS

	           call svdegas (teev, tiev, tnev, nem3, aion, aneut,
     .	              leh0, lchex, dsvione, dsvcxk, dsvionik)
	           svion_k = dsvionik
	           sv_cx_k = dsvcxk
	           svion_tot_k = svion_k + d_ratio * svion_e(i)

	         else if (iatdat.EQ.2) then          ! Thomas / Stacey

	           call calcxswms(teev, tiev, tnev, nem3, svion,
     .		      svcxk, svreck, svelk, svelnk)
	           sv_cx_k = svcxk + svelk
	           svion_tot_k = d_ratio * svion_e(i)

	         endif

	         invmfpk = ionDens(i) * (svion_tot_k + sv_cx_k) / v0k
	         mfp(kk,i) = 1.0 / invmfpk

c/    For the calculation of the charge exchange fraction, determine
c/    if we should use first collision effects:

	         if (ifrstcol.EQ.0) then
	            A_cxk(kk,i) = A_cx(i)
	         else
	            A_cxk(kk,i) = sv_cx_k / (sv_cx_k + svion_tot_k)
                 endif

c/    Adjacent cell is a wall segment:

	       else if (iType(k).EQ.2)   then

		  kw = k - (nCells + nPlasmReg)
		  if ((irefl.EQ.0).OR.((irefl.EQ.1).AND.
     .   		  (zwall(kw).LE.0.0))) then
		     E_0 = eneut
		     e_slow = ionTemp(i)
		     e_fast = ionTemp(i)
                  else if ((irefl.EQ.1).AND.(zwall(kw).GT.0.0)) then
		     E_0 = eneut
		     e_slow = twall(kw)
		     e_fast = refle(kw) * tneut(i,kk) / refln(kw)

                  endif
                  if(i_e0.eq.3)then
                    e_slow=1.5*e_slow
                    e_fast=1.5*e_fast
                  endif

		  tn0ev = 1000.0 * E_0
	          tnfev = 1000.0 * e_fast
	          tnsev = 1000.0 * e_slow

                  v0  = sqrt(v0fact * tn0ev * xj7ev / (protMass*aneut))
                  v0f = sqrt(v0fact * tnfev * xj7ev / (protMass*aneut))
                  v0s = sqrt(v0fact * tnsev * xj7ev / (protMass*aneut))

	          if (iatdat.EQ.0) then               ! Janev

                    svion_w = svion_i(i)
                    sv_cx_w0 = svcxi(aion, tiev, aneut, tn0ev)  ! Source
                    sv_cx_wf = svcxi(aion, tiev, aneut, tnfev)  ! Fast
                    sv_cx_ws = svcxi(aion, tiev, aneut, tnsev)  ! Slow

	            svion_tot_wf = svion_w + d_ratio * svion_e(i)
	            svion_tot_ws = svion_tot_wf
	            svion_tot_w0 = svion_tot_wf

	          else if (iatdat.EQ.1) then          ! DEGAS

c/    Source neutrals:
	            call svdegas (teev, tiev, tn0ev, nem3, aion, aneut,
     .	              leh0, lchex, dsvione, dsvcxw0, dsvioniw0)
c/    Fast neutrals:
	            call svdegas (teev, tiev, tnfev, nem3, aion, aneut,
     .	              leh0, lchex, dsvione, dsvcxwf, dsvioniwf)

c/    Slow neutrals:
	            call svdegas (teev, tiev, tnsev, nem3, aion, aneut,
     .	              leh0, lchex, dsvione, dsvcxws, dsvioniws)

	            svion_wf = dsvioniwf
	            svion_ws = dsvioniws
	            svion_w0 = dsvioniw0
	            sv_cx_wf = dsvcxwf
	            sv_cx_ws = dsvcxws
	            sv_cx_w0 = dsvcxw0

	            svion_tot_wf = svion_wf + d_ratio * svion_e(i)
	            svion_tot_ws = svion_ws + d_ratio * svion_e(i)
	            svion_tot_w0 = svion_w0 + d_ratio * svion_e(i)

	          else if (iatdat.EQ.2) then          ! Thomas / Stacey

c/    Source neutrals:
	            call calcxswms(teev, tiev, tn0ev, nem3, svion,
     .		        svcxw0, svrecw0, svelw0, svelnw0)

c/    Fast neutrals:
	            call calcxswms(teev, tiev, tnfev, nem3, svion,
     .		        svcxwf, svrecwf, svelwf, svelnwf)

c/    Slow neutrals:
	            call calcxswms(teev, tiev, tnsev, nem3, svion,
     .		        svcxws, svrecws, svelws, svelnws)

	            sv_cx_w0 = svcxw0 + svelw0
	            sv_cx_wf = svcxwf + svelwf
	            sv_cx_ws = svcxws + svelws
	            svion_tot_ws = d_ratio * svion_e(i)
	            svion_tot_w0 = svion_tot_ws
	            svion_tot_wf = svion_tot_ws

	          endif

	          invmfpw0 = ionDens(i)*(svion_tot_w0 + sv_cx_w0)/v0
	          invmfpwf = ionDens(i)*(svion_tot_wf + sv_cx_wf)/v0f
	          invmfpws = ionDens(i)*(svion_tot_ws + sv_cx_ws)/v0s

	          mfp_w0(kw) = 1.0 / invmfpw0
	          mfp_wf(kw) = 1.0 / invmfpwf
	          mfp_ws(kw) = 1.0 / invmfpws

c/    For the charge exchange fractions, we check again whether
c/    we should separate out first collisions:

                  if (ifrstcol.EQ.0) then
                    A_cxw0(kw) = A_cx(i)
		    A_cxwf(kw) = A_cx(i)
		    A_cxws(kw) = A_cx(i)
                  else
	            A_cxw0(kw) = sv_cx_w0 / (sv_cx_w0 + svion_tot_w0)
	            A_cxws(kw) = sv_cx_ws / (sv_cx_ws + svion_tot_ws)
	            A_cxwf(kw) = sv_cx_wf / (sv_cx_wf + svion_tot_wf)
	          endif

	       endif           ! end of IF loop for interface type

            enddo              ! end of DO loop over sides of cell -i-

c/    Case of volumetric sources:
c/    If S_ext(i) > 0, and if we want first collision effects,
c/    we must calculate the mean free path and charge exchange
c/    fraction for the first-collision neutrals assuming that their
c/    energy is equal to eneut_v

            if ((S_ext(i).GT.0.0).AND.(ifrstcol.EQ.1)) then

c/	       E_s0 = ionTemp(i)
c/#############################################
	       E_s0 = eneut_v
		 if(i_e0.eq.3)E_s0=1.5*ionTemp(i)
c/Volumetric Source now uses local ionTemp instead of eneut_v
c/This is because S_ex a recombination source.
c/E_s0=1.5*ionTemp(i)
c/E_s0=1.5*eneut_v
c/#############################################
c/               if(i_e0.eq.3)E_s0=1.5*ionTemp(i)
	       tnev0 = 1000.0 * E_s0

	       if (iatdat.EQ.0) then

                  svion_i0 = svioni(aion, tiev, aneut, tnev0)
                  sv_cx0 = svcxi(aion, tiev, aneut, tnev0)
	          svion_tot0 = svion_i(i) + d_ratio * svion_e(i)

	       else if (iatdat.EQ.1) then

	          call svdegas (teev, tiev, tnev0, nem3, aion, aneut,
     .  	       leh0, lchex, dsvione0, dsvcx0, dsvioni0)
	          svion_i0 = dsvioni0
	          sv_cx0 = dsvcx0
	          svion_tot0 = svion_i0 + d_ratio * svion_e(i)

	       else if (iatdat.EQ.2) then

	          call calcxswms(teev, tiev, tnev0, nem3, svion0, svcx0,
     .	             svrec0, svel0, sveln0)
	          sv_cx0 = svcx0 + svel0
	          svion_tot0 = d_ratio * svion_e(i)

	       endif

               vs0 = sqrt(v0fact * E_s0 * xj7kv / (protMass * aneut))

	       invmfp0 = ionDens(i) * (svion_tot0 + sv_cx0) / vs0
	       lmfp0(i) = 1.0 / invmfp0
	       A_cx0(i) = sv_cx0 / (sv_cx0 + svion_tot0)

            endif      ! End of IF loop on volumetric source

         endif                 ! end of IF loop for i_e0

      enddo   ! end of DO loop over cells



c/    Calculate albedo for plasma edge:

      if (nPlasmReg.GT.0) then
         do i = 1, nPlasmReg
	    j = nCells + i
            albedo(i) = albdfit_ad(A_cx(j))
         enddo
      endif

      return
      end

      real function albdfit_ad(x)

c/    This function calculates a fit to the albedo coefficient based
c/    on data from Monte Carlo simulations performed by Dingkang Zhang
c/    at Georgia Tech. The results agree also with the semi-analytical
c/    derivation of the albedo coefficient by P. Rafalski, in
c/    Nuclear Sci. & Eng., 19 (1964) 378.
c/
c/    Prepared by John Mandrekas, GIT, January 2004.

c/    y=(a+cx+ex^2+gx^3)/(1+bx+dx^2+fx^3)

c/    x       : charge exchange fraction, 0 <= x <= 1.0
c/    albdfit_ad : the albedo coefficient

      implicit none
      real a, b, c, d, e, f, g, x, y

      data a/5.9720174E-4/, b/-2.46848679/, c/0.2045041/,
     .     d/1.9744939/, e/-0.3818644/, f/-0.505836/, g/0.1769341/

      y = (a + x*(c + x*(e + x*g))) / (1.0 + x*(b + x*(d + f*x)))

      albdfit_ad = y

      return
      end

