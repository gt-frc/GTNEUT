c///////////////////// GAUSSIAN QUADRATURE ROUTINES ////////////////////
c/                                                                     /
c/    The following subroutines perform Gauss-Legendre integration for /
c/    the function -func- from -a- to -b- for different numbers of     /
c/    integration points. Since the integration points are symmetric,  /
c/    the arrays for x and w have a dimension of (number of points/2)  /
c/    The name of each routine, reflects the total number of interior  /
c/    integration points. For example, qgauss20 means that we have     /
c/    2X10 = 20 interior integration points.                           /
c/                                                                     /
c/    Each routine returns as ss the integral of the function -func-   /
c/    between -a- and -b-                                              /
c/                                                                     /
c/    Written by John Mandrekas, GIT, 08/12/93                         /
c/                                                                     /
c/    Record of Changes:                                               /
c/    -----------------                                                /
c/    06/08/98, jm: Changed x_i and w_i to double precision            / 
c/                                                                     /
c/    References:                                                      /
c/    ----------                                                       /
c/    1) Numerical Recipes, Chapter 4.5                                /
c/    2) Abramowitz and Stegun, Handbook of Mathematical Functions.    /
c/                                                                     /
c/    Modified by Zach Friis Feb 03, 2009                              /
c/    Added new Quadrature Set for n = 10,n=5                          /
c/    Experience oscilattions with higher order quadrature sets        /
c/    on certain problems.                                             / 
c////////////////////////// QGAUS6 ////////////////////////////////////

      subroutine qgauss6 (func, a, b, ss,numfun)
     
      implicit none
      external func

      integer numfun
      real a, b, ss(numfun)
      integer i, npoints, j
      real dx, xm, xr, funv1(numfun),funv2(numfun)
      
      parameter (npoints = 3)
      double precision x_i(npoints), w_i(npoints)
      save x_i, w_i
      
      data x_i/0.2386191860831969086305017, 0.6612093864662645136613996,
     .         0.9324695142031520278123016/

      data w_i/0.4679139345726910473898703, 0.3607615730481386075698335,
     .         0.1713244923791703450402961/

      xm = 0.5 * (b+a)
      xr = 0.5 * (b-a)
      do j=1,numfun
         ss(j)=0.0
      enddo
      do i = 1, npoints
         dx = xr * x_i(i)
         call func(xm+dx,numfun,funv1)
         call func(xm-dx,numfun,funv2)
         do j=1,numfun
            ss(j) = ss(j)+w_i(i)*(funv1(j)+funv2(j))
         enddo
      enddo
      do j=1,numfun
         ss(j) = xr * ss(j)
      enddo
      return
      end      
c////////////////////////// QGAUS10 ////////////////////////////////////

      subroutine qgauss10 (func, a, b, ss,numfun)
     
      implicit none
      external func

      integer numfun
      real a, b, ss(numfun)
      integer i, npoints, j
      real dx, xm, xr, funv1(numfun),funv2(numfun)
      
      parameter (npoints = 5)
      double precision x_i(npoints), w_i(npoints)
      save x_i, w_i
      
      data x_i/0.1488743389816312108848260, 0.4333953941292471907992659,
     .         0.6794095682990244062343274, 0.8650633666889845107320967,
     .         0.9739065285171717200779640/

      data w_i/0.2955242247147528701738930, 0.2692667193099963550912269,
     .         0.2190863625159820439955349, 0.1494513491505805931457763,
     .         0.0666713443086881375935688/

      xm = 0.5 * (b+a)
      xr = 0.5 * (b-a)
      do j=1,numfun
         ss(j)=0.0
      enddo
      do i = 1, npoints
         dx = xr * x_i(i)
         call func(xm+dx,numfun,funv1)
         call func(xm-dx,numfun,funv2)
         do j=1,numfun
            ss(j) = ss(j)+w_i(i)*(funv1(j)+funv2(j))
         enddo
      enddo
      do j=1,numfun
         ss(j) = xr * ss(j)
      enddo
      return
      end      

c////////////////////////// QGAUS20 ////////////////////////////////////

      subroutine qgauss20 (func, a, b, ss,numfun)
     
      implicit none
      external func

      integer numfun
      real a, b, ss(numfun)
      integer i, npoints, j
      real dx, xm, xr, funv1(numfun),funv2(numfun)

      parameter (npoints = 10)
      double precision x_i(npoints), w_i(npoints)
      save x_i, w_i

      data x_i/0.0765265211, 0.2277858511, 0.3737060887, 0.5108670019,
     .         0.6360536807, 0.7463319064, 0.8391169718, 0.9122344282,
     .         0.9639719272, 0.9931285991/

      data w_i/0.1527533871, 0.1491729864, 0.1420961093, 0.1316886384,
     .         0.1181945319, 0.1019301198, 0.0832767415, 0.0626720483,
     .         0.0406014298, 0.0176140071/

      xm = 0.5 * (b+a)
      xr = 0.5 * (b-a)
      do j=1,numfun
         ss(j)=0.0
      enddo
      do i = 1, npoints
         dx = xr * x_i(i)
         call func(xm+dx,numfun,funv1)
         call func(xm-dx,numfun,funv2)
         do j=1,numfun
            ss(j) = ss(j)+w_i(i)*(funv1(j)+funv2(j))
         enddo
      enddo
      do j=1,numfun
         ss(j) = xr * ss(j)
      enddo

      return
      end                                                              

c////////////////////////// QGAUS40 ////////////////////////////////////

      subroutine qgauss40 (func, a, b, ss,numfun)
     
      implicit none
      external func

      integer numfun
      real a, b, ss(numfun)
      integer i, npoints, j
      real dx, xm, xr, funv1(numfun),funv2(numfun)

      parameter (npoints = 20)
      double precision x_i(npoints), w_i(npoints)
      save x_i, w_i

      data x_i/3.8772417506051D-02, 0.11608407067526, 0.19269758070137,
     .         0.26815218500725, 0.34199409082576, 0.41377920437160,
     .         0.48307580168618, 0.54946712509513, 0.61255388966798,
     .         0.67195668461418, 0.72731825518993, 0.77830565142652,
     .         0.82461223083331, 0.86595950321226, 0.90209880696887,
     .         0.93281280827868, 0.95791681921379, 0.97725994998377,
     .         0.99072623869946, 0.99823770971056 /

      data w_i/7.7505947978425D-02, 7.7039818164248D-02, 
     .   7.6110361900626D-02, 7.4723169057968D-02, 7.2886582395804D-02,
     .   7.0611647391287D-02, 6.7912045815234D-02, 6.4804013456601D-02,
     .   6.1306242492929D-02, 5.7439769099392D-02, 5.3227846983937D-02,
     .   4.8695807635072D-02, 4.3870908185673D-02, 3.8782167974472D-02,
     .   3.3460195282546D-02, 2.7937006980016D-02, 2.2245849194167D-02,
     .   1.6421058381908D-02, 1.0498284531153D-02, 4.5212770985331D-03/

      xm = 0.5 * (b+a)
      xr = 0.5 * (b-a)
      do j=1,numfun
         ss(j)=0.0
      enddo
      do i = 1, npoints
         dx = xr * x_i(i)
         call func(xm+dx,numfun,funv1)
         call func(xm-dx,numfun,funv2)
         do j=1,numfun
            ss(j) = ss(j)+w_i(i)*(funv1(j)+funv2(j))
         enddo
      enddo
      do j=1,numfun
         ss(j) = xr * ss(j)
      enddo

      return
      end                                                              

c////////////////////////// QGAUS60 ////////////////////////////////////

      subroutine qgauss60 (func, a, b, ss,numfun)
     
      implicit none
      external func

      integer numfun
      real a, b, ss(numfun)
      integer i, npoints, j
      real dx, xm, xr, funv1(numfun),funv2(numfun)

      parameter (npoints = 30)
      double precision x_i(npoints), w_i(npoints)
      save x_i, w_i

      data x_i/2.5959772301248D-02, 7.7809333949537D-02, 
     .   0.12944913539694, 0.18073996487343, 0.23154355137603,
     .   0.28172293742326, 0.33114284826845, 0.37967005657680,
     .   0.42717374158308, 0.47352584176171, 0.51860140005857,
     .   0.56227890075394, 0.60444059704851, 0.64497282848948,
     .   0.68376632738136, 0.72071651335573, 0.75572377530659,
     .   0.78869373993226, 0.81953752616215, 0.84817198478593,
     .   0.87451992264690, 0.89851031081005, 0.92007847617763,
     .   0.93916627611642, 0.95572225584000, 0.96970178876505,
     .   0.98106720175260, 0.98978789522222, 0.99584052511884,
     .   0.99921012322744/

      data w_i/ 5.1907877631221D-02, 5.1767943174910D-02,
     .   5.1488451500981D-02, 5.1070156069856D-02, 5.0514184532509D-02,
     .   4.9822035690550D-02, 4.8995575455757D-02, 4.8037031819971D-02,
     .   4.6948988848912D-02, 4.5734379716114D-02, 4.4396478795787D-02,
     .   4.2938892835936D-02, 4.1365551235585D-02, 3.9680695452381D-02,
     .   3.7888867569243D-02, 3.5994898051084D-02, 3.4003892724946D-02,
     .   3.1921219019296D-02, 2.9752491500789D-02, 2.7503556749925D-02,
     .   2.5180477621521D-02, 2.2789516943998D-02, 2.0337120729457D-02,
     .   1.7829901014207D-02, 1.5274618596784D-02, 1.2678166476812D-02,
     .   1.0047557182267D-02, 7.3899311633457D-03, 4.7127299269533D-03,
     .   2.0268119688734D-03 /

      xm = 0.5 * (b+a)
      xr = 0.5 * (b-a)
      do j=1,numfun
         ss(j)=0.0
      enddo
      do i = 1, npoints
         dx = xr * x_i(i)
         call func(xm+dx,numfun,funv1)
         call func(xm-dx,numfun,funv2)
         do j=1,numfun
            ss(j) = ss(j)+w_i(i)*(funv1(j)+funv2(j))
         enddo
      enddo
      do j=1,numfun
         ss(j) = xr * ss(j)
      enddo

      return
      end                                                              

c////////////////////////// QGAUS80 ////////////////////////////////////

      subroutine qgauss80 (func, a, b, ss,numfun)
     
      implicit none
      external func

      integer numfun
      real a, b, ss(numfun)
      integer i, npoints, j
      real dx, xm, xr, funv1(numfun),funv2(numfun)

      parameter (npoints = 40)
      double precision x_i(npoints), w_i(npoints)
      save x_i, w_i

      data x_i/1.9511383256794D-02, 5.8504437152421D-02,
     .   9.7408398441585D-02, 0.13616402280914, 0.17471229183265,
     .   0.21299450285767, 0.25095235839227, 0.28852805488451,
     .   0.32566437074770, 0.36230475349949, 0.39839340588197,
     .   0.43387537083176, 0.46869661517054, 0.50280411188879,
     .   0.53614592089713, 0.56867126812271, 0.60033062282975,
     .   0.63107577304687, 0.66085989898612, 0.68963764434203,
     .   0.71736518536210, 0.74400029758360, 0.76950242013504,
     .   0.79383271750461, 0.81695413868146, 0.83883147358026,
     .   0.85943140666311, 0.87872256767821, 0.89667557943877,
     .   0.91326310257176, 0.92845987717245, 0.94224276130987,
     .   0.95459076634363, 0.96548508904380, 0.97490914058573,
     .   0.98284857273863, 0.98929130249976, 0.99422754096569,
     .   0.99764986439824, 0.99955382265163/

      data w_i/ 3.9017813656307D-02, 3.8958395962770D-02,
     .   3.8839651059052D-02, 3.8661759774076D-02, 3.8424993006959D-02,
     .   3.8129711314478D-02, 3.7776364362001D-02, 3.7365490238731D-02,
     .   3.6897714638276D-02, 3.6373749905836D-02, 3.5794393953416D-02,
     .   3.5160529044748D-02, 3.4473120451754D-02, 3.3733214984612D-02,
     .   3.2941939397645D-02, 3.2100498673488D-02, 3.1210174188115D-02,
     .   3.0272321759558D-02, 2.9288369583268D-02, 2.8259816057277D-02,
     .   2.7188227500486D-02, 2.6075235767565D-02, 2.4922535764115D-02,
     .   2.3731882865930D-02, 2.2505090246332D-02, 2.1244026115782D-02,
     .   1.9950610878142D-02, 1.8626814208299D-02, 1.7274652056269D-02,
     .   1.5896183583726D-02, 1.4493508040509D-02, 1.3068761592401D-02,
     .   1.1624114120798D-02, 1.0161766041103D-02, 8.6839452692603D-03,
     .   7.1929047681150D-03, 5.6909224513905D-03, 4.1803131246950D-03,
     .   2.6635335895127D-03, 1.1449500031869D-03 /

      xm = 0.5 * (b+a)
      xr = 0.5 * (b-a)
      do j=1,numfun
         ss(j)=0.0
      enddo
      do i = 1, npoints
         dx = xr * x_i(i)
         call func(xm+dx,numfun,funv1)
         call func(xm-dx,numfun,funv2)
         do j=1,numfun
            ss(j) = ss(j)+w_i(i)*(funv1(j)+funv2(j))
         enddo
      enddo
      do j=1,numfun
         ss(j) = xr * ss(j)
      enddo

      return
      end                                                              

c////////////////////////// QGAUS100 ///////////////////////////////////

      subroutine qgauss100 (func, a, b, ss,numfun)
     
      implicit none
      external func

      integer numfun
      real a, b, ss(numfun)
      integer i, npoints, j
      real dx, xm, xr, funv1(numfun),funv2(numfun)

      parameter (npoints = 50)
      double precision x_i(npoints), w_i(npoints)
      save x_i, w_i

      data x_i/ 1.5628984421543D-02, 4.6871682421592D-02,
     .  7.8068582813437D-02, 0.10918920358006, 0.14020313723611,
     .  0.17108008053860, 0.20178986409574, 0.23230248184497,
     .  0.26258812037150, 0.29261718803847, 0.32236034390053,
     .  0.35178852637242, 0.38087298162463, 0.40958529167830,
     .  0.43789740217203, 0.46578164977336, 0.49321078920819,
     .  0.52015801988176, 0.54659701206509, 0.57250193262138,
     .  0.59784747024718, 0.62260886020371, 0.64676190851413,
     .  0.67028301560314, 0.69314919935580, 0.71533811757306,
     .  0.73682808980202, 0.75759811851971, 0.77762790964950,
     .  0.79689789239031, 0.81538923833918, 0.83308387988840,
     .  0.84996452787959, 0.86601468849716, 0.88121867938502,
     .  0.89556164497073, 0.90902957098253, 0.92160929814533,
     .  0.93328853504308, 0.94405587013626, 0.95390078292549,
     .  0.96281365425582, 0.97078577576371, 0.97780935848692,
     .  0.98387754070606, 0.98898439524299, 0.99312493703744,
     .  0.99629513473313, 0.99849195063960, 0.99971372677344/

      data w_i/ 3.1255423453863D-02, 3.1224884254849D-02,
     . 3.1163835696210D-02, 3.1072337427567D-02, 3.0950478850491D-02,
     . 3.0798379031153D-02, 3.0616186583980D-02, 3.0404079526455D-02,
     . 3.0162265105169D-02, 2.9890979593333D-02, 2.9590488059913D-02,
     . 2.9261084110638D-02, 2.8903089601125D-02, 2.8516854322395D-02,
     . 2.8102755659101D-02, 2.7661198220792D-02, 2.7192613446577D-02,
     . 2.6697459183571D-02, 2.6176219239546D-02, 2.5629402910208D-02,
     . 2.5057544481580D-02, 2.4461202707957D-02, 2.3840960265968D-02,
     . 2.3197423185254D-02, 2.2531220256336D-02, 2.1843002416247D-02,
     . 2.1133442112528D-02, 2.0403232646209D-02, 1.9653087494435D-02,
     . 1.8883739613375D-02, 1.8095940722128D-02, 1.7290460568324D-02,
     . 1.6468086176145D-02, 1.5629621077546D-02, 1.4775884527441D-02,
     . 1.3907710703719D-02, 1.3025947892972D-02, 1.2131457662979D-02,
     . 1.1225114023186D-02, 1.0307802574869D-02, 9.3804196536944D-03,
     . 8.4438714696690D-03, 7.4990732554647D-03, 6.5469484508452D-03,
     . 5.5884280038651D-03, 4.6244500634204D-03, 3.6559612013180D-03,
     . 2.6839253715535D-03, 1.7093926535180D-03, 7.3463449050559D-04/

      xm = 0.5 * (b+a)
      xr = 0.5 * (b-a)
      do j=1,numfun
         ss(j)=0.0
      enddo
      do i = 1, npoints
         dx = xr * x_i(i)
         call func(xm+dx,numfun,funv1)
         call func(xm-dx,numfun,funv2)
         do j=1,numfun
            ss(j) = ss(j)+w_i(i)*(funv1(j)+funv2(j))
         enddo
      enddo
      do j=1,numfun
         ss(j) = xr * ss(j)
      enddo
     

      return
      end                                                              
