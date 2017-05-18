      real function svione (te)

c/    This function evaluates the electron impact ionization 
c/    reactivity (<sigma*v>) as a function of the electron temperature
c/    using the logarithmic polynomial approximation of Janev, et al.
c/
c/    Reference:
c/    ---------
c/    R.K. Janev et al. 'Elementary Processes in Hydrogen-Helium Plasmas,
c/    Springer-Verlag, 1987 
c/
c/    e + H(1s) --> e + H+ + e (Reaction 2.1.5)
c/    Notice that only six significant digits have been kept, per
c/    Janev's suggestion (page 256)
c/
c/    Written by John Mandrekas, GIT, 10/05/2001 for the GTNEUT code
c/
c/    Parameters:
c/    ----------
c/    te      : electron temperature (eV)
c/    svione  : electron impact ionization rate (m^3/s)
c/
c/    Comments:
c/    --------
c/    The electron temperature range is: 0.1 eV <= te < 20 keV

      implicit none
      integer nfit,i 
      parameter (nfit = 9)
      real bn(nfit)
      real te, t, alnte, sum
     

      data bn
     . / -3.271396786375e+01,  1.353655609057e+01, -5.739328757388e+00,
     .	  1.563154982022e+00, -2.877056004391e-01,  3.482559773737e-02,
     .	 -2.631976175590e-03,  1.119543953861e-04, -2.039149852002e-06/

c/    Make sure electron temperature is within range of validity:
      
      t = te
      if (te.LT.0.1) t = 0.1
      if (te.GT.2.0e4) t = 2.0e4

      alnte = alog(t)

      sum = 0.0
      do i = 1, nfit
           sum = sum + bn(i) * (alnte**(i-1))
      enddo

      svione = 1.0e-6 * exp(sum)

      return
      end 

c///////////////////////////////////////////////////////////////////////
c/    Ion Impact Ionization:

      real function svioni (aion, ti, aneut, e0)

c/    This function evaluates the ion impact ionization reactivity
c/    (<sigma*v>) as a function of the ion temperature and the neutral
c/    energy using the double logarithmic polynomial approximation of 
c/    Janev, et al.
c/
c/    Reference:
c/    ---------
c/    R.K. Janev et al. 'Elementary Processes in Hydrogen-Helium Plasmas,
c/    Springer-Verlag, 1987 
c/
c/    p + H(1s) --> p + H+ + e (Reaction 3.1.6)
c/
c/    Written by John Mandrekas, GIT, 10/05/2001 for the GTNEUT code
c/
c/    Parameters:
c/    ----------
c/    aion  : ion mass (1 for H, 2 for D, etc.)
c/    ti    : ion temperature (eV)
c/    aneut : neutral mass (1 for H, 2 for D, etc.)
c/    e0    : ion temperature (eV)
c/
c/    Output:
c/    ------
c/    svioni : ion impact ionization reaction rate (m^3/s)
c/
c/    Comments:
c/    --------
c/    The ion temperature  and neutral energy ranges are: 
c/    0.1 eV <= Ti, E0 <= 20 keV
c/    Energies and temperatures are scaled by mass ratios


      implicit none
      integer i, j, nfit
      parameter(nfit = 9)
      real an(nfit,nfit)
      real ti, e0, aion, aneut, dblsum, alne0, alnti, t, e

      data ((an(i,j), j = 1, nfit), i = 1, nfit)
     . /-1.617454916209e+02, 1.021458246570e+02, -5.712267930902e+01,
     .   2.140540272484e+01,-4.767517412803e+00, 6.293295208376e-01,
     .  -4.858173640838e-02, 2.031177914273e-03,-3.557982934756e-05,
     .   1.767238902030e+01,-7.102574692619e+01, 4.246688953154e+01,
     .  -1.128638171243e+01, 1.661679851896e+00,-1.476754423056e-01,
     .   8.175790218529e-03,-2.732531913524e-04, 4.398387454014e-06,
     .  -4.334843983767e+01, 3.855259623260e+01,-1.316883030631e+01,
     .   2.145592145856e+00,-1.467281287038e-01,-2.915256218527e-03,
     .   1.092542891192e-03,-6.205102802216e-05, 1.158798945435e-06,
     .   2.464254915383e+01,-1.283426276878e+01, 2.369698902002e+00,
     .  -1.506665823159e-01,-8.144926683660e-03, 2.231505500086e-03,
     .  -2.210941355372e-04, 1.310924337643e-05,-3.431837053957e-07,
     .  -5.439093405254e+00, 2.357085001656e+00,-2.961732508220e-01,
     .  -9.917174972226e-04, 1.935894665907e-03,-1.679264493005e-05,
     .   5.532386419162e-08,-1.121430499351e-06, 5.960280736984e-08,
     .   5.959975304236e-01,-2.391382925527e-01, 2.789277301925e-02,
     .   8.562387824450e-05,-1.340759667335e-04,-5.927455645560e-06,
     .   5.820264508685e-07, 7.694068657107e-08,-4.972708712807e-09,
     .  -3.361958123977e-02, 1.289667246580e-02,-1.858739201548e-03,
     .   9.235982885753e-05, 9.875232214392e-06,-1.680823118052e-06,
     .   3.019916624608e-08, 6.889325889968e-09,-3.171970185702e-10,
     .   8.706597041685e-04,-3.140899683782e-04, 7.343984485463e-05,
     .  -8.601564864429e-06,-6.467790579320e-07, 1.734797315767e-07,
     .   2.523651535182e-09,-1.719633613108e-09, 7.332933714195e-11,
     .  -6.359765062372e-06, 1.742836004704e-06,-1.235536456998e-06,
     .   2.257852760280e-07, 1.608335682237e-08,-3.855914336143e-09,
     .  -3.556222618473e-10, 7.627265694554e-11,-2.960493966948e-12/


c/    Scale with isotope mass (in case of D or T):
      
      t = ti / aion
      e = e0 / aneut

c/    Make sure we are within limits of validity:

      if (t.LT.0.1) t = 0.1
      if (e.LT.0.1) e = 0.1
      if (t.GT.2.0e4) t = 2.0e4
      if (e.GT.2.0e4) e = 2.0e4

      alnti = alog(t)
      alne0 = alog(e)
      
      dblsum = 0.0
      do i = 1, nfit
         do j = 1, nfit
          dblsum = dblsum + an(i,j)*(alne0**(i-1))*(alnti**(j-1))
         enddo
      enddo

      svioni = 1.0e-6 * exp(dblsum)

      return
      end 

c///////////////////////////////////////////////////////////////////////
c/    Charge Exchange:

      real function svcxi (aion, ti, aneut, e0)

c/    This function evaluates the ion charge exchange  reactivity
c/    (<sigma*v>) as a function of the ion temperature and the neutral
c/    energy using the double logarithmic polynomial approximation of 
c/    Janev, et al.
c/
c/    Reference:
c/    ---------
c/    R.K. Janev et al. 'Elementary Processes in Hydrogen-Helium Plasmas,
c/    Springer-Verlag, 1987 
c/
c/    p + H(1s) --> H(1s) + p  (Reaction 3.1.8)
c/
c/    Written by John Mandrekas, GIT, 10/05/2001 for the GTNEUT code
c/
c/    Parameters:
c/    ----------
c/    aion  : ion mass (1 for H, 2 for D, etc.)
c/    ti    : ion temperature (eV)
c/    aneut : neutral mass (1 for H, 2 for D, etc.)
c/    e0    : ion temperature (eV)
c/
c/    Output:
c/    ------
c/    svcxi : ion impact ionization reaction rate (m^3/s)
c/
c/    Comments:
c/    --------
c/    The ion temperature  and neutral energy ranges are: 
c/    0.1 eV <= Ti, E0 <= 20 keV
c/    Energies and temperatures are scaled by mass ratios


      implicit none
      integer i, j, nfit
      parameter(nfit = 9)
      real an(nfit,nfit)
      real ti, e0, aion, aneut, dblsum, alne0, alnti, t, e

      data ((an(i,j), j = 1,nfit), i = 1, nfit)
     . /-1.829079581680e+01, 2.169137615703e-01, 4.307131243894e-02,
     . -5.754895093075e-04,-1.552077120204e-03,-1.876800283030e-04,
     .  1.125490270962e-04,-1.238982763007e-05, 4.163596197181e-07,
     .  1.640252721210e-01,-1.106722014459e-01, 8.948693624917e-03,
     .  6.062141761233e-03,-1.210431587568e-03,-4.052878751584e-05,
     .  2.875900435895e-05,-2.616998139678e-06, 7.558092849125e-08,
     .  3.364564509137e-02,-1.382158680424e-03,-1.209480567154e-02,
     .  1.075907881928e-03, 8.297212635856e-04,-1.907025662962e-04,
     .  1.338839628570e-05,-1.171762874107e-07,-1.328404104165e-08,
     .  9.530225559189e-03, 7.348786286628e-03,-3.675019470470e-04,
     . -8.119301728339e-04, 1.361661816974e-04, 1.141663041636e-05,
     . -4.340802793033e-06, 3.517971869029e-07,-9.170850253981e-09,
     . -8.519413589968e-04,-6.343059502294e-04, 1.039643390686e-03,
     .  8.911036876068e-06,-1.008928628425e-04, 1.775681984457e-05,
     . -7.003521917385e-07,-4.928692832866e-08, 3.208853883734e-09,
     . -1.247583860943e-03,-1.919569450380e-04,-1.553840717902e-04,
     .  3.175388949811e-05, 1.080693990468e-05,-3.149286923815e-06,
     .  2.318308730487e-07, 1.756388998863e-10,-3.952740758950e-10,
     .  3.014307545716e-04, 4.075019351738e-05, 2.670827249272e-06,
     . -4.515123641755e-06, 5.106059413591e-07, 3.105491554749e-08,
     . -6.030983538280e-09,-1.446756795654e-10, 2.739558475782e-11,
     . -2.499323170044e-05,-2.850044983009e-06, 7.695300597935e-07,
     .  2.187439283954e-07,-1.299275586093e-07, 2.274394089017e-08,
     . -1.755944926274e-09, 7.143183138281e-11,-1.693040208927e-12,
     .  6.932627237765e-07, 6.966822400446e-08,-3.783302281524e-08,
     . -2.911233951880e-09, 5.117133050290e-09,-1.130988250912e-09,
     .  1.005189187279e-10,-3.989884105603e-12, 6.388219930167e-14/

c/    Scale with isotope mass (in case of D or T):
      
      t = ti / aion
      e = e0 / aneut

c/    Make sure we are within limits of validity:

      if (t.LT.0.1) t = 0.1
      if (e.LT.0.1) e = 0.1
      if (t.GT.2.0e4) t = 2.0e4
      if (e.GT.2.0e4) e = 2.0e4

      alnti = alog(t)
      alne0 = alog(e)
      
      dblsum = 0.0
      do i = 1, nfit
         do j = 1, nfit
          dblsum = dblsum + an(i,j)*(alne0**(i-1))*(alnti**(j-1))
         enddo
      enddo

      svcxi = 1.0e-6 * exp(dblsum)

      return
      end 

c/    Freeman and Jones <sigma*v>_e for compatibility with older results:

      real function svefj(te)

c/    This function calculates the electron impact ionization rate
c/    of neutral hydrogenic atoms using the Freeman and Jones fitting
c/    formula for compatibility with older calculations. It is better
c/    to use the more recent Janev fits for the various atomic rates.

c/    Reference:
c/    ---------
c/    R.L. Freeman, E.M. Jones, Culham Laboratory Report, clm-r137, 1974.

c/    Variables:
c/    ---------
c/    te  : Electron temperature (eV)
c/    svefj : <sigma*v>_e in m^3/s

      implicit none
      real te
      real tkev, tl, a(7)

      data a/-3.140212e+1, -3.024379e-1, -5.616546e-2,
     .       +7.902886e-3, +1.246713e-3, +2.217222e-4, -9.486967e-5/

      save a

      tkev = 1.0e-3 * te

c/    Fit from Freeman and Jones is valid for 1 eV < Te < 100 keV:
c/    Local value of temperature is forced to be within these limits
      
      if (tkev.LT.1.0e-3) tkev = 1.0e-3
      if (tkev.GT.100.0) tkev = 100.0

      tl = alog(tkev)

      svefj = exp(a(1)+tl*(a(2)+tl*(a(3)+tl*(a(4)+tl*(a(5)+tl*(a(6)
     .		+ tl*a(7)))))))
      
      return
      end
