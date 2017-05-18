      block data ndata

c/    This is a block data file that initializes several variables that
c/    are in COMMON.

      implicit none
      include 'comiou.inc'
      include 'consts.inc'

      data nin /10/, nout /12/, ioeh /14/, iocx /15/, ndbug /20/,
     .   nsdbug /22/

      data xj7kv /1.6022e-16/, xj7ev /1.6022e-19/, 
     .     elecMass /9.1094e-31/, protMass  /1.6726e-27/, 
     .     pi /3.141592654/

      end
