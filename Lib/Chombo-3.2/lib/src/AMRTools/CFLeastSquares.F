
#include "CONSTANTS.H"
#include "REAL.H"

      real_t function xpower1d0avg(lower, upper, p)
!     
!     xpower1d0avg(lower, upper p)
!     returns the average of x^p - C on the 1D interval (lower:upper),
!     where C = 0 unless p is even > 0, in which case C = (1/2)^p / (p+1).
!
!     The purpose of including C is to make the average over all
!     fine cells zero with the exception of p = 0.
!     If p = 0 then this is (upper - lower) / (upper - lower) / (p + 1) = 1.
      implicit none
!     --- arguments ---
      real_t lower, upper
      integer p
!     --- locals ---
      real_t avg

!     Integral(x^p) = (upper^(p+1) - lower^(p+1)) / (p + 1)
!     <x^p> = Integral(x^p) / (upper - lower)
!           = ((upper^(p+1) - lower^(p+1)) / (upper - lower)) / (p+1)
      avg = (upper**(p+1) - lower**(p+1)) / (upper - lower)
      if ((p .gt. 0) .and. (mod(p, 2) .eq. 0)) then
         avg = avg - (half**p)
      endif
      avg = avg / ((p + 1) * one)
      xpower1d0avg = avg

      return
      end


      real_t function xpower1dcoarseind0avg(off, p)
!     
!     xpower1dcoarseind0avg(offset, p)
!     returns the average of x^p - C over the cell offset by "off"
!     where C = 0 unless p is even > 0, in which case C = (1/2)^p / (p+1).
!
!     The purpose of including C is to make the average over all
!     fine cells zero with the exception of p = 0.
!     If p = 0 then this works out to be 1.
      implicit none
!     --- arguments ---
      integer off, p
!     --- locals ---
      real_t lower, upper, avg
      real_t xpower1d0avg

      lower = (off*one) - half
      upper = (off*one) + half
      avg = xpower1d0avg(lower, upper, p)
      xpower1dcoarseind0avg = avg

      return
      end


      subroutine xpowercoarseind0avg(
     &     nnbrs, npwrs, nbrs, pwrs, xpower)
!
!     xpowercoarseind0avg(nnbrs, npwrs, nbrs, pwrs, xpower)
!     returns xpower(1:nnbrs, 1:npwrs)
!     where xpower(inbr, p) = <x^p - C> over neighbor indexed by inbr
      implicit none
!     --- arguments ---
      integer nnbrs, npwrs
      integer nbrs(nnbrs, CH_SPACEDIM), pwrs(npwrs, CH_SPACEDIM)
      real_t xpower(nnbrs, npwrs)
!     --- locals ---
      integer inbr, ipwr, nbr, pwr, idim
      real_t xpower1dcoarseind0avg

      do inbr = 1, nnbrs
         do ipwr = 1, npwrs
            xpower(inbr, ipwr) = one
         end do
      end do

      do idim = 1, CH_SPACEDIM
         do inbr = 1, nnbrs
            nbr = nbrs(inbr, idim)
            do ipwr = 1, npwrs
               pwr = pwrs(ipwr, idim)
               xpower(inbr, ipwr) = xpower(inbr, ipwr) *
     &              xpower1dcoarseind0avg(nbr, pwr)
            end do
         end do
      end do

      return
      end


      real_t function xpower1dfineind0avg(off, p, nref)
!     
!     xpower1dfineind0avg(off, p, nref)
!     returns average of x^p over interval [off, off+1]/nref - 1/2 .
!
!     off in range 0:nref-1, corresponding to being inside coarse box 0.
!
      implicit none
!     --- arguments ---
      integer off, p, nref
!     --- locals ---
      real_t lower, upper, avg
      real_t xpower1d0avg

!     lower coordinate of coarse box 0 is -1/2
      lower = ( off     *one) / (nref*one) - half
      upper = ((off + 1)*one) / (nref*one) - half
      avg = xpower1d0avg(lower, upper, p)
      xpower1dfineind0avg = avg

      return
      end


      subroutine xpowerfineind0avg(
     &     nfine, npwrs, nref, pwrs, fine, xpower)
!
!     Returns xpower(ifine, ipwr) = average of x^(pwrs(ipwr)) over cell ifine.
!
      implicit none
!     --- arguments ---
      integer nfine, npwrs, nref
      integer pwrs(npwrs, CH_SPACEDIM)
      integer fine(nfine, CH_SPACEDIM)
      real_t xpower(nfine, npwrs)
!     --- locals ---
      integer ifine, ipwr, pwr, idim, off
      real_t xpower1dfineind0avg

      do ifine = 1, nfine
         do ipwr = 1, npwrs
            xpower(ifine, ipwr) = one
         end do
      end do

      do idim = 1, CH_SPACEDIM
         do ifine = 1, nfine
            off = fine(ifine, idim)
            do ipwr = 1, npwrs
               pwr = pwrs(ipwr, idim)
               xpower(ifine, ipwr) = xpower(ifine, ipwr) *
     &              xpower1dfineind0avg(off, pwr, nref)
            end do
         end do
      end do
            
      return
      end


      integer function findnbr0(nnbrs, nbrs)
!
!     Returns inbr in nbrs(1:nnbrs, 1:D) such that nnbrs(inbr, 1:D) = 0.
!     This is the index of the central coarse cell.
!
      implicit none
!     --- arguments ---
      integer nnbrs
      integer nbrs(nnbrs, CH_SPACEDIM)
!     --- locals ---
      integer inbr, idim
      logical found

      do inbr = 1, nnbrs
         found = .true.
         do idim = 1, CH_SPACEDIM
            found = found .and. (nbrs(inbr, idim) .eq. 0)
         end do
         if (found) then
            findnbr0 = inbr
         endif
      end do
      
      return
      end


      subroutine coarsefineleastsquares(
     &     nref, nnbrs, npwrs, nfine,
     &     nbrs, pwrs, fine, mat)
!
!     We have calculated nbrs in FourthOrderInterpStencil::define, so
!     we can then send nbrs to this subroutine.
!
!     We OMIT the zero power from pwrs.
!     Note that we need fine(nfine, CH_SPACEDIM).
!
!     Output is mat(nfine, nnbrs).  We should reorder this
!     to an FArrayBox on base fine box, with nnbrs components.
!
      implicit none
!     --- arguments ---
      integer nref, nnbrs, npwrs, nfine
      integer nbrs(nnbrs, CH_SPACEDIM)
      integer pwrs(npwrs, CH_SPACEDIM)
      integer fine(nfine, CH_SPACEDIM)
      real_t mat(nfine, nnbrs)
!     --- locals ---
      real_t xpowc(nnbrs, npwrs)
      real_t ata(npwrs, npwrs)
      real_t work(npwrs)
      real_t coeffs(npwrs, nnbrs)
      real_t sumrow(npwrs)
      real_t xpowf(nfine, npwrs)
      integer ipiv(npwrs), info, nbr0
      integer inbr, ipwr, ifine
      character*1 prec
      integer findnbr0

#ifdef CH_USE_DOUBLE
      prec = 'D'
#else
      prec = 'S'
#endif

!     Set xpowc(inbr, ipwr) to
!     average of x^(pwrs(ipwr, 1:D)) on cell nbrs(inbr, 1:D).
      call xpowercoarseind0avg(nnbrs, npwrs, nbrs, pwrs, xpowc)

!     Set ata(1:npwrs, 1:npwrs) =
!     xpowc(1:nnbrs, 1:npwrs)^T * xpowc(1:nnbrs, 1:npwrs) .
!     
!     That is:
!     ata(ipwr1, ipwr2) = 
!     sum_inbr ( xpowc(inbr, ipwr1) * xpowc(inbr, ipwr2) )
#ifdef CH_USE_DOUBLE
      call DGEMM('T', 'N', npwrs, npwrs, nnbrs, one, xpowc,
     &     nnbrs, xpowc, nnbrs, zero, ata, npwrs)
#else
      call SGEMM('T', 'N', npwrs, npwrs, nnbrs, one, xpowc,
     &     nnbrs, xpowc, nnbrs, zero, ata, npwrs)
#endif

!     Replace ata(1:npwrs, 1:npwrs) by its LU decomposition.
!     We do this because LAPACK requires it to compute the inverse.
#ifdef CH_USE_DOUBLE
      call DGETRF(npwrs, npwrs, ata, npwrs, ipiv, info)
#else
      call SGETRF(npwrs, npwrs, ata, npwrs, ipiv, info)
#endif
      if (info .lt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GETRF returned ', info
         print *, 'Illegal value in argument number ', -info
         call MAYDAY_ERROR()
      elseif (info .gt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GETRF returned ', info
         print *, 'Singular matrix, 0 at diagonal element ', info
         call MAYDAY_ERROR()
      endif

!     Replace ata(1:npwrs, 1:npwrs) by its inverse.
#ifdef CH_USE_DOUBLE
      call DGETRI(npwrs, ata, npwrs, ipiv, work, npwrs, info)
#else
      call SGETRI(npwrs, ata, npwrs, ipiv, work, npwrs, info)
#endif
      if (info .lt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GETRI returned ', info
         print *, 'Illegal value in argument number ', -info
         call MAYDAY_ERROR()
      elseif (info .gt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GETRI returned ', info
         print *, 'Singular matrix, 0 at diagonal element ', info
         call MAYDAY_ERROR()
      endif

!     Set coeffs(1:npwrs, 1:nnbrs) =
!     ata(1:npwrs, 1:npwrs) * xpowc(1:nnbrs, 1:npwrs)^T .
!
!     That is:
!     coeffs(ipwr, inbr) =
!     sum_ipwr2 ( ata(ipwr, ipwr2) * xpowc(inbr, ipwr2) )
#ifdef CH_USE_DOUBLE
      call DGEMM('N', 'T', npwrs, nnbrs, npwrs, one, ata,
     &     npwrs, xpowc, nnbrs,0, coeffs, npwrs)
#else
      call SGEMM('N', 'T', npwrs, nnbrs, npwrs, one, ata,
     &     npwrs, xpowc, nnbrs,0, coeffs, npwrs)
#endif

!     Get sumrow(1:npwrs), sums of the rows of coeffs:
!     sumrow(ipwr) = sum( coeffs(ipwr, 1:nnbrs) )
!                  = sum_inbr ( coeffs(ipwr, inbr) )
      do ipwr = 1, npwrs
         sumrow(ipwr) = 0
         do inbr = 1, nnbrs
            sumrow(ipwr) = sumrow(ipwr) + coeffs(ipwr, inbr)
         end do
      end do

!     nbr0 is index in nbrs(1:nnbrs, :) of the central coarse cell.
!     This is stored as nbrs(nbr0, 1:D) = 0
      nbr0 = findnbr0(nnbrs, nbrs)

!     Reduce coeffs(1:npwrs, nbr0) by sumrow(1:npwrs):
!     coeffs(1:npwrs, nbr0) -= sumrow(1:npwrs)
      do ipwr = 1, npwrs
         coeffs(ipwr, nbr0) = coeffs(ipwr, nbr0) - sumrow(ipwr)
      end do
      
!     Set xpowf(ifine, ipwr) to
!     average of x^(pwrs(ipwr, 1:D)) on cell ifine.
      call xpowerfineind0avg(nfine, npwrs, nref, pwrs, fine, xpowf)

!     Set mat(1:nfine, 1:nnbrs) =
!     xpowf(1:nfine, 1:npwrs) * coeffs(1:npwrs, 1:nnbrs).
!
!     That is:
!     Set mat(ifine, inbr) =
!     sum_ipwr ( xpowf(ifine, ipwr) * coeffs(ipwr, inbr) )
#ifdef CH_USE_DOUBLE
      call DGEMM('N', 'N', nfine, nnbrs, npwrs, one, xpowf,
     &     nfine, coeffs, npwrs, 0, mat, nfine)
#else
      call SGEMM('N', 'N', nfine, nnbrs, npwrs, one, xpowf,
     &     nfine, coeffs, npwrs, 0, mat, nfine)
#endif

!     Increment mat(1:nfine, nbr0) by 1:
!     mat(1:nfine, nbr0) += 1.
      do ifine = 1, nfine
         mat(ifine, nbr0) = mat(ifine, nbr0) + one
      end do
      
      return
      end
