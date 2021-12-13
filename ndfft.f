      program nd_fft
c
c This is the f90 scalarized code
c Modified to include Open MP Optimizations -10/22/01 
c

      implicit none

      integer, parameter :: numprocs = 1
      integer, parameter :: base = 2
!      include 'parameter.h'
      integer, parameter :: num1=16
      integer, parameter :: log1=4
      integer, parameter :: num2=16
      integer, parameter :: log2=4
      integer, parameter :: num3=16
      integer, parameter :: log3=4
      integer, parameter :: dim = 3
**************
      integer ivec(0:num1*num2*num3 -1)
      integer opl, opr
**************
      complex zvec(0: num1*num2*num3 -1)
      integer logs(0:dim-1),svec(0:dim-1)
      integer tvec(0:dim-1),iotav(0:dim -1)
      integer rlogs(0:dim-1),rsvec(0:dim-1),rpscan(0:dim-1)
      integer i,ind,k,l,m
      integer sp,total

! initializations 
!
      logs = (/log1,log2,log3/)
      svec = (/num1,num2,num3/)
      total = product(svec)
!
**
*  Initialize zvec. This is the vector we do the fft on.
**
      do i =0,total-1 
         zvec(i) = 2*i
      end do

*********
*      MoA-Psi design for 3-d fft. The design is
*        dimension independent, i.e. it's an n-d fft.
*      LRM  1/04.
*********
       iotav =(/((i),i=0,dim-1)/)
*** Reverse shape, logs, etc. since its column major in Fortran
       rsvec = rev(svec)
       rlogs = rev(logs)
       rpscan= rev(pscan(svec)) 
*******do the fft over all dimensions
       do i = 0,dim -1
*******
* This is the permutation vector used in the transpose
* Notice how i is used
******
       tvec = cat(take(i,iotav),
     &           rev(take(-(dim-i),iotav)))
******
* These are the coefficients used in the outerproduct to
* generate all indices need.
******
       opl = rsvec(tvec(0))
       ivec(total-opl:total-1)=rpscan(tvec(0))
     &           *(/((i),i=0,rsvec(tvec(0))-1)/) 
       do k = 0,dim-2 
          ind = 0
         do l = 0,(product(take((k+1),rsvec(tvec))))-1
           do m = 0,rsvec(tvec(k+1))-1
             opl = product(take((k+2),rsvec(tvec)))
             opr = product(take((k+1),rsvec(tvec)))
             ivec((total-opl)+ind)= ivec((total-opr)+l)
     &           +m*rpscan(tvec((k+1)))
             ind=ind+1
           end do
         end do
      end do
***********************
* Now all indices are available. Do the FFT. 
***********************
      sp= rsvec(tvec(dim-1))        
      do k = 0,total-1,sp
        zvec(ivec(k:k+sp-1))=fft(zvec(ivec(k:k+sp-1)),
     &     rsvec(tvec(dim-1)),rlogs(tvec(dim-1))  )
      end do
**********************
      end do 
**********************
!      print*,'TRANSFORMED ARRAY'
!      do k = 0,total-1
!          write(6,"(f12.2,f12.2)") zvec(k)
!      end do
       
*********
      contains
      function pscan(svec)
      integer, intent(in) :: svec(0:)
      integer pscan(0:size(svec)-1)
********
      pscan(0) = 1
      do i =1,size(svec)-1
      pscan(i)=pscan(i-1)*svec(i-1)
      end do      
      end function
      function cat(ivecl,ivecr)
      integer, intent(in) ::  ivecl(0:),ivecr(0:)
      integer cat(0:size(ivecl)+size(ivecr)-1)
      integer i
      do i =0,(size(ivecl)+size(ivecr))-1
         if (i < size(ivecl)) then
            cat(i) = ivecl(i)
         else
            cat(i) = ivecr(i-size(ivecl) )
         end if

      end do 
      end function
      function rot(sigma,xvec)
      integer, intent(in) :: sigma
      integer, intent(in) :: xvec(0:)
      integer i
      integer rot(0:size(xvec)-1)
      do i = 0,size(xvec)-1
         rot(i) = xvec(mod(i+sigma,size(xvec)))
      end do
      end function
      function rev(ivec)
      integer, intent(in) :: ivec(0:)
      integer rev(0:size(ivec)-1)  
      integer i
         do i = 0,size(ivec)-1
         rev(i) = ivec(size(ivec)-(i+1)) 
         end do
      end function 
      function take(sigma,ivec)
      integer, intent(in) :: sigma
      integer, intent(in) :: ivec(0:)
      integer take(0:abs(sigma)-1)  
      if (sigma > 0) then
         take(0:sigma-1) = ivec(0:sigma -1)
      else
         take(0:abs(sigma)-1) = ivec(size(ivec)+sigma: 
     &         size(ivec)+sigma+ abs(sigma)-1)
      end if
      end function
!
      function fft(z_p,n,log2n)
      integer n
      integer log2n
      complex z_p(0:(n-1))
      complex fft(0:size(z_p)-1)
! i_ and j_ are used as indicies or loop variables
      integer i_, j_, breakpoint, np, msize, p_
      integer g, L, Lnew
      integer ind,ind2
      real  pi
      complex i, c, d
! ind <-> index vector
! This is an in-place algorithm
      complex  ww_p(0:n/2-1)

! This is the initial condition 
      g=0
      L=INT(2**g)
      np = numprocs
      msize = n/np
! We need pi and the imaginary number i
      i = (0.0,1.0)
      pi = 3.14159
!
! We see this below in the assignment to "c".
      do i_=0,n-1
      fft(i_)=z_p(bit_rev(i_,log2n))
      end do
!
! This is the major loop which performs the fft. There
! are log2n loops and the algorithm is "in-place".
!
! There are now two copies of this loop to run in parallel 10/22/01
      
! One will run  1,breakpoint-1 and the other breakpoint,log2n

!!$OMP PARALLEL PRIVATE(ww_p, z_p, Lnew, c, d, i_)
!!$OMP DO
      do p_ = 0, np-1
        z_p(msize*p_:msize*(p_+1)-1:1)=fft(msize*p_:msize*(p_+1)-1:1)
        do g = 1,breakpoint-1
          Lnew = INT(2**g) 
!
! Note the limits on j_. The do loop that follows has the same 
! bounds as the j_ loop below it.
! 
          do j_=0,Lnew/2-1 
            ww_p(j_) = EXP((2*pi*i*j_)/Lnew)
          end do
          do i_= (p_*msize), (p_+1)*msize-1,Lnew
            do j_=0,Lnew/2-1
              c = ww_p(j_)*z_p(i_+j_+Lnew/2)               
              d=z_p(i_+j_)
              z_p((/ i_+j_ /)) = (/ d+c /)
              z_p((/ i_+j_+Lnew/2 /)) = (/ d-c /)
             end do
           end do
c       L=Lnew
       end do
      fft(msize*p_:msize*(p_+1)-1:1)=z_p(msize*p_:msize*(p_+1)-1:1)
      end do
!!$OMP END DO
!!$OMP END PARALLEL 
!



c The second loop 10/22/01 -

!!$OMP PARALLEL PRIVATE(ww_p, z_p, Lnew, c, d, i_)
!!$OMP DO
      do p_ = 0, np-1
        z_p(p_:n-1:np)=fft(p_:n-1:np)
        do g = breakpoint,log2n
          Lnew = INT(2**g) 
!!
!! Note the limits on j_. The do loop that follows has the same 
!! bounds as the j_ loop below it.
!! 
          do j_=0,Lnew/2-1 
            ww_p(j_) = EXP((2*pi*i*j_)/Lnew)
          end do
          do i_=0,n-1,Lnew
            do j_=p_,Lnew/2-1,np
              c = ww_p(j_)*z_p(i_+j_+Lnew/2)               
              d = z_p(i_+j_)
              z_p((/ i_+j_ /)) = (/ d+c /)
              z_p((/ i_+j_+Lnew/2 /)) = (/ d-c /)
            end do
          end do
c         L=Lnew
       end do
      fft(p_:n-1:np)=z_p(p_:n-1:np)
      end do
!!$OMP END DO
!!$OMP END PARALLEL


       end function fft

!
      function bit_rev(i_,number_bits)
! computes bit reversal of a non-negative integer
      integer i_           ! value whose bits are to be reversed
      integer bit_rev
      integer number_bits  ! number of bits in i_
! constraints: number_bits > 1; 0 <= i_ < 2**number_bits
! way of computing bit reversal based on Van Loan
! (Algorithm 1.5.1 on page 39)
      integer j_,m_,q_,s_
      j_ = 0
      m_ = i_
      do q_ = 0,number_bits-1
            s_ = m_/2
            j_ = 2*j_+m_-2*s_
            m_ = s_
      end do
      bit_rev = j_
      end function bit_rev
      end program nd_fft
