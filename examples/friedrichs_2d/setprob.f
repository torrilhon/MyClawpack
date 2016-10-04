      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname
      common /cparam/ tau, Re
c
c     # Set the relaxation parameter
c     # These values are passed in a common block
c

      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
      read(7,*) tau
      read(7,*) Re

      return
      end
