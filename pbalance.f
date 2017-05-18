      subroutine pbalance

c/    This subroutine performs a global particle balance, as
c/    an additional test of the correctness of the solution.

c/    John Mandrekas, June 2003

      implicit none
      include 'neutGlob.inc'

c/    Local variables:
      
      integer i, i_true, j, jj, kk, kpl, kw, kw_true, ll
      real sum_gext, sum_gion, sum_ion, sum_vol, sum_core, sum_wall

c/    First, calculate the gains:
c/    ---------------------------
      
c/    Ion and neutral fluxes from wall segments:

      sum_gext = 0.0
      sum_gion = 0.0
      do kw = 1, nWallSegm
	 sum_gext = sum_gext + g_ex(kw)
	 sum_gion = sum_gion + g_ion(kw)
      enddo

c/    Volumetric source contributions, if any:    
      
      sum_vol = 0.0
      do i = 1, nCells
	 sum_vol = sum_vol + S_ext(i)
      enddo

c/    Now, calculate the losses:
c/    --------------------------

c/    Ionization loss first:
      
      sum_ion = 0.0
      do i = 1, nCells
	 sum_ion = sum_ion + inzRate(i)
      enddo

c/    Losses through the core plasma-interface:
      
      sum_core = 0.0
      if (nPlasmReg.GT.0) then
         do kpl = 1, nPlasmReg
	    i_true = nCells + kpl
	    do jj = 1, nSides(i_true)
	       j = adjCell(jj,i_true)
	       do kk = 1, nSides(j)
		  if (adjCell(kk,j).EQ.i_true) then
		     sum_core = sum_core +gflux(j,kk,1) * 
     .                  (1.0 - albedo(kpl))
                  endif
               enddo
            enddo
         enddo
      endif

c/    Wall losses due to wall absorption or pumping:
      
      sum_wall = 0.0
      do kw = 1, nWallSegm
	 kw_true = nCells + nPlasmReg + kw
	 j = adjCell(1,kw_true)
	 do ll = 1, nSides(j)
	    if (adjCell(ll,j).EQ.kw_true) then
               sum_wall = sum_wall + fwabsorb(kw) * (1.0 - refln(kw)) *
     .            (gflux(j,ll,1) + g_ion(kw))
            endif
         enddo
      enddo

c/    Calculate global particle balance parameters:
      
      fluxn_tot = srcNormFact *sum_gext        ! external surface fluxes (gas puffing)
      fluxi_tot = srcNormFact *sum_gion        ! ion fluxes (recycling)
      volsrc_tot = srcNormFact *sum_vol        ! volumetric sources
      coreloss_tot = srcNormFact *sum_core     ! lost to the core plasma region
      wloss_tot = srcNormFact *sum_wall        ! lost by wall absorption or pumping
      ionloss_tot = srcNormFact *sum_ion       ! lost by ionization

c/    Total particles entering volume:

      tot_part_in = (fluxn_tot + fluxi_tot + volsrc_tot)

c/    Total particles lost in volume:

      tot_part_out = (coreloss_tot + wloss_tot + ionloss_tot)

c/    Relative percentage error:

      partbalnerr = 100.0 * (tot_part_in - tot_part_out) / tot_part_in


      return
      end
