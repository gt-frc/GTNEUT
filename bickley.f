c/    This file contains functions for the evaluation of the Bickley
c/    functions of order 3 and 4, that are needed in the calculation
c/    of the transmission coefficients.

      real function bickley3(x)

c/    This is an approximation for the Ki3 Bickley-Naylor
c/    function. 
c/    Ref: F.G. Lether, J. Quant. Spectrosc. Radiat. Transfer,
c/         43 (1990) 187.
c/    Created by John Mandrekas, on 4/7/99 for the GTNEUT code.

      implicit none
      real x, numer, denom, Ki3

      numer = 1.88571 + x*(1.56855 + 0.202569*x)
      denom = 2.40135 + x*(1.43711 + 0.161532*x)

      Ki3 = exp(-x) * numer / (Sqrt(x+1.0) * denom)

      bickley3 = Ki3

      return
      end


      real function bickley2(x)
      implicit none
      real x, numer,denom
      numer=1.00026+x*(6.18878+2.76828*x)
      denom=1+x*(6.25525+2.21050*x)
      bickley2=numer/denom*exp(-x)/sqrt(1+x)
      return
      end

      real function bickley1(x)
      implicit none
      real x,numer,denom
      numer=0.0048385+x*(0.6646409+x*(5.451198+3.974590*x))
      denom=0.003080874+x*(0.4370382+x*(4.530696
     .     +x*(7.855812+3.173373*x)))
      bickley1=exp(-x)*sqrt(1+x)*numer/denom
      return
      end
