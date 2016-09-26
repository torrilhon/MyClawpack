      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname
      common /cparam/ tau,p1,p2, fsource
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
      read(7,*) tau
      read(7,*) p1
      read(7,*) p2
      read(7,*) fsource

      return
      end
