      subroutine zstop(ntty, nout, n)

c///////////////////////////////////////////////////////////////////////
c/                                                                     /      
c/    This subroutine terminates the run, and writes the reason for    /
c/    termination as determined by the flag n, to the terminal (ntty)  /
c/    and the main output file (nout).                                 /
c/    Written by John Mandrekas, GIT, 12-12-93, for the 0-D code.      /
c/                                                                     /
c/////////////////////////////////////////////////////////////////////// 
      
      implicit none
      
      integer ntty, nout, n
      
      write (ntty, 5000) n
      write (nout, 5000) n

 5000 format (1x, 'Exit code = ', i2)
 
      stop
      end