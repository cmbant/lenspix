!MPI Healpix routines, including generation of spin-s maps, 
!mapping gradients of scalars and the exact and approx weak lensed CMB 
!Antony Lewis Sept 2005, Based on Healpix 1.2 or 2.x

!Requires linking to Healpix libraries:
!See http://www.eso.org/science/healpix/

!Sign conventions follow Healpix/CMBFAST/CAMB
!Most easily used using HealpixObj.f90 wrapper routines
!Compile with -DMPIPIX to use MPI

!Performance could be improved by changing the theta-CPU sharing
!However scaling is quite good up to 50 or so processors for high res transforms
!Temporary arrays use more memory than non-MPI routines
!For compatibility with Healpix input/output alm arrays are not packed (2*waste of mem)

!Jan 2005: improved/fixed polarization lens rotation factors. Minor fixes.
!Sept 2005: fixed bug in map2polalm

module MPIstuff
implicit none
double precision starttime
#ifdef MPIPIX
    include "mpif.h"
    integer ::  DebugMsgs = 1
    integer MPIstatus(MPI_STATUS_SIZE), ierr
    integer SP_MPI,CSP_MPI 
#endif

end module MPIstuff

module spinalm_tools
  use utilities, only: die_alloc
  use healpix_types
  use healpix_fft, only : real_fft
  IMPLICIT none

  Type HealpixInfo
     integer :: nside, lmax, nalms_max, Lastlmax
     logical pol
     REAL(KIND=DP), DIMENSION(:,:), Pointer :: w8ring_TQU
     INTEGER(I8B), DIMENSION(:), pointer :: istart_south, istart_north  
     COMPLEX(DPC),DIMENSION(:), pointer :: trig
     REAL(DP), DIMENSION(:), Pointer :: recfac, Lambda_slm
     integer MpiID, MPISize, MpiStat, last_nph
     integer(I4B), dimension(:), pointer :: ith_start, ith_end
     integer, dimension(:), pointer :: North_Start, North_Size, South_Start, South_Size
  end type HealpixInfo

  integer, parameter :: EB_sign = -1
    !definition: for pol {}_2 a_{lm} = EB_sign*(E + iB)_{lm}
    !EB_sign = -1 corresponds to Healpix and CAMB/CMBFAST conventions

  logical :: mmax_approx = .true.
   !when true, uses that fact that don't need high m near the poles because the legendre 
   !functions are tiny for m >> l sin(theta)

  ! keep everything private unless stated otherwise
  private
  ! define large and small numbers used to renormalise the recursion on the Legendre Polynomials
  real(KIND=DP), private, PARAMETER :: FL_LARGE = 1.0e30_dp
  real(KIND=DP), private, PARAMETER :: FL_SMALL = 1.0e-30_dp
  real(KIND=DP), private :: OVFLOW, UNFLOW, ScaleFactors(-10:10)  
 
  ! make (front end) routines public
  public :: spinalm2map, alm2GradientMap, map2spinalm,scalalm2map, mmax_approx, HealpixInfo, &
            HealpixInit,HealpixFree, map2scalalm, a_ix, scalalm2LensedMap, &
            alm2Lensedmap, map2polalm, polalm2map, alm2LensedQuadContrib, EB_sign, &
            alm2LensedmapInterp   
contains

 function GeteTime()
      use MPIStuff
      double precision GeteTime
#ifndef MPIPIX
      real etime,tarray(2)
      external etime
      GeteTime = etime(tarray)
#else
      GeteTime = MPI_WTime()
#endif
 end function GeteTime

  function a_ix(lmax, l, m) result(index)
    integer, intent(in) :: lmax, l, m
    integer :: index
    index = (m*(2*lmax-m+1))/2 + l + 1
  end function a_ix

   subroutine HealpixInit(H, nside, lmax, HasPol, w8dir, method)
    USE fitstools, ONLY : getsize_fits, input_map
    use MPIStuff
    use healpix_types
    Type (HealpixInfo) :: H
    Integer, optional, intent(in) :: method
    logical, intent(in), optional :: HasPol
    character(LEN=*), optional, intent(in) :: w8dir
    integer, intent(in) :: nside, lmax
    real(dp) dw8, phi0, logOVFLOW
    integer npol, n_rings
    character(LEN=120) :: sstr, filename
    REAL(SP), DIMENSION(:,:), allocatable :: w8
    integer ith, nph, i, m, status, ierror, delta, st
    ! Changed for new division between threads
    Real(sp) :: mean_pix !Mean number of pixels in each section of northern hemisphere
    Integer :: pixels, row
    Integer :: division = 1 
      !Determines whether to give each section
      ! equal numbers of rows (1), or equal numbers of pixels (2)
      ! (2) is much faster for 'exact' lensing

    CHARACTER(LEN=*), PARAMETER :: code = 'HealpixInit'
 

#ifndef MPIPIX
        call HealpixFree(H)
!If MPI must call healpixFree manually
#endif    
       nullify(H%recfac,H%Lambda_slm)

       call HealpixInitTrig(H,nside,lmax)

       npol = 1
       if (present(HasPol)) then
        npol = 3
       end if
       H%pol = HasPol
  

    allocate(H%w8ring_TQU(1:2*nside,1:max(1,npol)))
 
    if (present(w8dir)) then
         allocate(w8(1:2*nside,1:max(1,npol)))
 
         write (sstr,"(I5.5)") nside
         filename= trim(w8dir)//"weight_ring_n"//trim(sstr)//".fits"

         n_rings = 2 * nside
     
         if (getsize_fits(filename) /= n_rings) then
            write (*,*) 'HealpixInit:wrong file'//trim(filename)
            stop
         endif
     
         if (HasPol) then
            call input_map(filename, w8, n_rings, 3, fmissval=0.0_sp)
         else
            call input_map(filename, w8, n_rings, 1, fmissval=0.0_sp)
         endif

         H%w8ring_TQU =  1 + w8
         deallocate(w8)
     else

       H%w8ring_TQU=1
  
     endif

!Get factors for making well behaved Ylm
    OVFLOW=exp(log(FL_LARGE))
    UNFLOW=exp(log(FL_SMALL))
    logOVFLOW=log(FL_LARGE)

    do i=-10,10
     ScaleFactors(i) = exp(i*logOVFLOW)
    end do

!Mpi properties
    H%MpiID = 0; H%MpiSize = 1
    H%MpiStat = 0


#ifdef MPIPIX
        if (SP==KIND(1.d0)) then
         SP_MPI = MPI_DOUBLE_PRECISION
         CSP_MPI = MPI_DOUBLE_COMPLEX
        else if (SP == KIND(1.)) then
         SP_MPI = MPI_REAL
         CSP_MPI=  MPI_COMPLEX
        else
         stop 'Unknown SP KIND for MPI'
        end if
               
        call mpi_comm_rank(mpi_comm_world,H%MpiID,ierror)
        if (ierror/=MPI_SUCCESS) stop 'HealpixInit: MPI rank'

        call mpi_comm_size(mpi_comm_world,H%MpiSize,ierror)
     !   call mpi_buffer_attach(H%MPIBuffer,MPIBUfSize, ierror)
#endif

!Sectioning of the sphere between threads
!This is not very optimal
!Should give more rows to each node when section is near the pole

! New version SJS 15/12/2004 for equal pixels per thread
! mean_pix is the average number of pixels given to each thread

#ifdef MPIPIX
        If (present(method)) division = method
#endif
        pixels = 0
        mean_pix = nside*(6*nside+2)/H%MpiSize

        delta = (2*nside)/H%MpiSize
        st = 1 + mod(2*nside,H%MpiSize)

       allocate(H%ith_start(0:H%MpiSIze-1), H%ith_end(0:H%MpiSIze-1), H%North_Start(0:H%MpiSIze-1), &
         H%North_Size(0:H%MpiSIze-1), H%South_Start(0:H%MpiSIze-1), H%South_Size(0:H%MpiSIze-1))
       H%ith_start = 1
       H%ith_end = 2*nside


        do i= 0, H%MpiSize -1 
    
           If (division == 2) Then

              ! New method - divide into ~equal numbers of pixels
              ! Should be significantly faster if using 'exact' lensing method
              ! very marginaly faster if using interpolation method 
              ! ideally need a third method which gives less to the poles
              ! if doing interpolation
              if (i == 0) then
                 H%ith_start(i) = 1
              else
                 H%ith_start(i) = H%ith_end(i-1) + 1
              end if

              row = H%ith_start(i)-1
              do while (pixels .LT. (i+1.0)*mean_pix)
                 row = row + 1
                 nph = 4*nside
                 if (row .LT. nside) nph = 4*row
                 pixels = pixels + nph
              end do
              H%ith_end(i) = row
              If (i == (H%MpiSize-1)) H%ith_end(i) = 2*nside

           Else If (division == 1) Then
              !divide into equal numbers of rows

            if (i == 0) then
             !Do more, but poles are faster anyway
              H%ith_start(i) = 1
              H%ith_end(i) = st + delta-1 
            else
              H%ith_start(i) = st + i*delta
              H%ith_end(i) =  H%ith_start(i) +  delta -1
            end if
            
           Else
              Stop 'HealpixInit : Unknown method'
           End If


            if (H%ith_end(i)< nside) then  
                  nph = 4*H%ith_end(i)
               else                  
                  nph = 4*nside
            endif
       
            H%North_start(i) = H%istart_north(H%ith_start(i)-1)
            H%North_Size(i) = H%istart_north(H%ith_end(i)-1) + nph &
                              -H%North_start(i)
        
            if (H%ith_start(i) < nside) then  
                  nph = 4*H%ith_start(i)
               else                   
                  nph = 4*nside
            endif
            if (H%ith_end(i) == nside*2) then
             H%South_start(i) = H%istart_south(H%ith_end(i)-1)
            else
             H%South_start(i) = H%istart_south(H%ith_end(i))
            end if
            H%South_Size(i) = H%istart_south(H%ith_start(i)) + nph &
                              - H%South_start(i)


        end do
#ifdef MPIPIX
        if (H%MpiID>0) call MessageLoop(H)
#endif
 
  end subroutine HealpixInit

   subroutine HealpixInitTrig(H, nside, lmax)
    use MPIStuff
    use healpix_types
    Type (HealpixInfo) :: H
    integer, intent(in) :: lmax, nside
    integer ith, i, m, status, nph
    real(dp) Phi0
    CHARACTER(LEN=*), PARAMETER :: code = 'HealpixTrig'


       nullify(H%trig)
       H%last_nph = -1
       H%lmax = lmax
       H%nalms_max = ((lmax+1)*(lmax+2))/2
       H%nside = nside
       H%Lastlmax = 0

        ALLOCATE(H%istart_north(0:2*nside),stat = status)
        if (status /= 0) call die_alloc(code,'istart_north')

        ALLOCATE(H%istart_south(0:2*nside),stat = status)
        if (status /= 0) call die_alloc(code,'istart_south')

        H%istart_north(0)=0
        H%istart_south(0)=12*int(nside,I8B)**2
        do ith=1,2*nside
           if (ith.lt.nside) then  ! polar cap (north)
              nph = 4*ith
           else                   ! tropical band (north) + equator
              nph = 4*nside
           endif
           H%istart_north(ith)=H%istart_north(ith-1)+nph
           H%istart_south(ith)=H%istart_south(ith-1)-nph
        enddo

   end subroutine HealpixInitTrig

   subroutine HealpixInfo_GetTrig(H, nph)
     Type (HealpixInfo) :: H
     integer, intent(in) :: nph    
     integer status, m
     real(dp) phi0    

     if (H%last_nph /= nph) then
 
       deallocate(H%trig,stat = status)
       ALLOCATE(H%trig(0:max(2*H%nside,H%lmax)),stat = status) 
 
       H%trig=1
       phi0=PI/DBLE((nph/4)*4)
       do m=0,max(2*H%nside,H%lmax)
          H%trig(m)= CMPLX( DCOS(m*phi0), DSIN(m*phi0), kind=DP)
       enddo
       H%last_nph = nph
      end if

   end  subroutine HealpixInfo_GetTrig

  
   function ScaleFactor(i)
    integer, intent(in) :: i
    real(dp) :: ScaleFactor
    
     if (i>-10) then
        ScaleFactor = ScaleFactors(i)
     else
        ScaleFactor = 0
     end if

   end function ScaleFactor

   subroutine HealpixFree(H)
    Type (HealpixInfo) :: H
    integer status
#ifdef MPIPIX
    if (H%MpiID == 0) call SendMessages(H, 'EXIT')
#endif    
    deallocate(H%w8ring_TQU, stat = status)
    deallocate(H%istart_north, H%istart_south, stat = status)
    deallocate(H%trig, stat = status)
    deallocate(H%recfac, stat = status)
    deallocate(H%ith_start, H%ith_end,H%North_Start, H%North_Size, & 
           H%South_Start, H%South_Size, stat = status)
    nullify(H%w8ring_TQU)
   end subroutine HealpixFree


   subroutine HealpixInitRecfac(H,nlmax)
     Type (HealpixInfo) :: H
     INTEGER(I4B), intent(in):: nlmax
     integer(I8B) :: m, l, fm2, fl2
     integer status, a_ix
     integer l2, m2

     if (H%MpiID > 0 .and. associated(H%recfac) .and. nlmax == H%Lastlmax) return   
     call HealpixFreeRecfac(H)  
     H%Lastlmax = nlmax
     deallocate(H%recfac,stat= status)
     ALLOCATE(H%recfac(((nlmax+1)*(nlmax+2))/2),stat = status)    
     if (status /= 0) call die_alloc('HealpixInitRecfac','recfac')

     a_ix = 0
     do m = 0, nlmax
      m2 = m**2
      do l = m, nlmax
        a_ix = a_ix + 1
        l2 = (l+1)**2
        H%recfac(a_ix) = SQRT( real(4 * l2 - 1,dp) / (l2-m2) )
      end do
     end do 
 
   end subroutine HealpixInitRecfac

   
          
  subroutine HealpixFreeRecfac(H)
      Type (HealpixInfo) :: H
         integer status

         if (H%MpiID > 0) return   !cache it as have loads of memory

         deallocate(H%recfac,stat= status)


  end subroutine HealpixFreeRecfac

  
  subroutine spinring_synthesis(H,nlmax,datain,nph,dataout,kphi0,mmax_ring)
 !Don't fully follow the signs here, but the answer is correct
 !Note no point using FFTW etc as FFT is a negligible fraction of computation cost
    !=======================================================================
    !     RING_SYNTHESIS
    !       called by alm2map
    !       calls     real_fft
    !
    !     dataout(j) = sum_m datain(m) * exp(i*m*phi(j)) 
    !     with phi(j) = j*2pi/nph + kphi0*pi/nph and kphi0 =0 or 1
    !
    !     as the set of frequencies {m} is larger than nph, 
    !     we wrap frequencies within {0..nph-1}
    !     ie  m = k*nph + m' with m' in {0..nph-1}
    !     then
    !     noting bw(m') = exp(i*m'*phi0) 
    !                   * sum_k (datain(k*nph+m') exp(i*k*pi*kphi0))
    !        with bw(nph-m') = CONJ(bw(m')) (if datain(-m) = CONJ(datain(m)))
    !     dataout(j) = sum_m' [ bw(m') exp (i*j*m'*2pi/nph) ]
    !                = Fourier Transform of bw
    !        is real
    !
    !         NB nph is not necessarily a power of 2
    !
    !=======================================================================

    Type (HealpixInfo) :: H

    INTEGER(I4B) :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax
    INTEGER(I4B), INTENT(IN) :: mmax_ring
    INTEGER(I4B), INTENT(IN) :: nph, kphi0

    COMPLEX(DPC), DIMENSION(0:nlmax), INTENT(IN) :: datain
    REAL(SP),     DIMENSION(0:nph-1), INTENT(OUT)     :: dataout
    REAL(DP),     DIMENSION(0:nph-1)     :: data
    REAL(DP) :: phi0
    INTEGER(I4B) :: i,iw,ksign,m,k,kshift
    COMPLEX(DPC), DIMENSION(0:nph-1) :: bw
    COMPLEX(DPC) :: dat
    INTEGER(I4B) :: status

   
    !=======================================================================

    call HealpixInfo_GetTrig(H, nph)

    nsmax = H%nside
    ksign = + 1

    kshift = (-1)**kphi0  ! either 1 or -1
    bw(0:nph-1) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)

    !     all frequencies [-m,m] are wrapped in [0,nph-1]
    bw(0)=datain(0)
    do m  = 1, mmax_ring                        ! in -nlmax, nlmax
       iw = MODULO(m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + datain(m)*(kshift**k)  ! complex number
       iw = MODULO(-m, nph)  ! between 0 and nph-1  = m', F90 intrisic
       k  = (-m - iw) / nph                ! number of 'turns'
       bw(iw) = bw(iw) + CONJG(datain(m))*(kshift**k)  ! complex number
    enddo
    !     kshift**k = 1       for even turn numbers
    !               = 1 or -1 for odd  turn numbers : results from the shift in space

    !     applies the shift in position <-> phase factor in Fourier space
    data(0)=REAL(bw(0))

    !Data is in packed storage
    do iw = 1, nph/2  -1
       m = ksign*(iw)
       if(kphi0==1) then
          dat =bw(iw) * H%trig(m)
       else
          dat =bw(iw)
       endif
       data(iw*2-1 ) = REAL(dat)
       data(iw*2) = AIMAG(dat)

    enddo
!    nph is always even for Healpix
    iw=nph/2
   m = ksign*(iw)
   if(kphi0==1) then
       dat =bw(iw) * H%trig(m)
    else
      dat =bw(iw)
    endif
    data(iw*2-1) = REAL(dat)
 
    call real_fft (data, backward=.true.)
    !     ^^^^^^^^^^^^
    dataout=REAL(data(0:nph-1))

    RETURN
  END subroutine spinring_synthesis
  
  subroutine alm2GradientMap(H, inlmax, alm, map_QU)
 !Get the map of the gradient of alm (know pure E, so quicker than general routine)
 !internally use EB_sign=1 convention, though result is independent
    use alm_tools
    use MPIstuff
    Type (HealpixInfo) :: H
    INTEGER(I4B), INTENT(IN) :: inlmax 
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:)  :: alm
    COMPLEX(SPC), INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map_QU
    COMPLEX(SPC), DIMENSION(:), pointer :: map2
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0, nlmax

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: factor_1, factor_2
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U

    CHARACTER(LEN=*), PARAMETER :: code = 'ALM2GRADIENTMAP'
    COMPLEX(DPC) ::  b_north_Q(0:H%lmax), b_north_U(0:H%lmax)
    COMPLEX(DPC) ::  b_south_Q(0:H%lmax), b_south_U(0:H%lmax)
    INTEGER(I4B) :: status,par_lm, a_ix

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(:),   ALLOCATABLE :: ringR, ringI
    integer mmax_ring, nalms
    double precision Initime
    !=======================================================================

     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiID==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    call MPI_BCAST(nlmax,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms),stat = status )
     if (status /= 0) call die_alloc(code,'alm2')
     if (H%MpiID==0) call Alm2PackAlm(alm,alm2,nlmax)
    
#ifdef MPIPIX
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code//' Got alm ',H%MpiID, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1), stat = status)
     if (status /= 0) call die_alloc(code,'map2')    
#else
     map2 => map_QU
#endif

    ALLOCATE(lam_fact(nalms),stat = status)    
 
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(ringR(0:4*nsmax-1),ringI(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    call HealpixInitRecfac(H,nlmax)
    call GetLamfact(lam_fact, nlmax)
   
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )
    !     --------------------------------------------

    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)   ! 0 <= cos theta < 1

       !        cos(theta) in the pixelisation scheme
       if (ith < nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2 !cos theta
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       
       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       b_north_Q(0:nlmax) = 0
       b_north_U(0:nlmax) = 0
       b_south_Q(0:nlmax) = 0
       b_south_U(0:nlmax) = 0
       
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1
       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = -1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)

          lam_lm = corfac*lam_mm/OVFLOW     !  actual lambda_mm      
          
          a_ix = a_ix + 1
          if (m >=1) then
!normal_l cancels with gradient, sign from gradient 
              lambda_x = - lam_lm * fm / sth 
              lambda_w = -lambda_x * cth

              b_n_Q =  lambda_w * alm2(a_ix)
              b_s_Q =  par_lm * b_n_Q

              b_n_U = (0,-1)* lambda_x * alm2(a_ix)
              b_s_U = -par_lm * b_n_U

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac*lam_mm/OVFLOW ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl

             a_ix = a_ix + 1

             a_w = -1 / sth
             lambda_x = a_w * fm  * lam_lm
             lambda_w = a_w * (lam_fact(a_ix)*lam_lm1m - fl*cth*lam_lm) 

             factor_1 =  lambda_w * alm2(a_ix)
             b_n_Q = b_n_Q +          factor_1 
             b_s_Q = b_s_Q + par_lm * factor_1 ! X has a diff. parity

             factor_2 =   (0,1) *  lambda_x * alm2(a_ix) 
             b_n_U = b_n_U - factor_2
             b_s_U = b_s_U + par_lm * factor_2

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif

          enddo

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U

       enddo

       call spinring_synthesis(H,nlmax, b_north_Q, nph, ringR, kphi0,mmax_ring)
       call spinring_synthesis(H,nlmax, b_north_U, nph, ringI, kphi0,mmax_ring)
       map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) = cmplx(RingR(0:nph-1),RingI(0:nph-1))

       if (ith  <  2*nsmax) then
          call spinring_synthesis(H,nlmax, b_south_Q, nph, ringR, kphi0,mmax_ring)
          call spinring_synthesis(H,nlmax, b_south_U, nph, ringI, kphi0,mmax_ring)
          map2(H%istart_south(ith):H%istart_south(ith)+nph-1) = cmplx(RingR(0:nph-1),RingI(0:nph-1))
       endif

    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(ringR,ringI)
    deallocate(alm2)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiID
    StartTime = Getetime()
    call MPI_GATHERV(map2(H%North_Start(H%MpiID)),H%North_Size(H%MpiID),CSP_MPI, &
       map_QU,H%North_Size,H%North_Start,CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2(H%South_Start(H%MpiID)),H%South_Size(H%MpiID),CSP_MPI, &
       map_QU,H%South_Size,H%South_Start,CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiID, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2) 
#endif

  end subroutine alm2GradientMap

  subroutine GetLamfact(lam_fact, nlmax)
  real(dp) lam_fact(*), fm2
  integer, intent(in) :: nlmax
  integer a_ix,l,m

    a_ix = 0
    do m = 0, nlmax
      fm2 = real(m,dp) **2
      a_ix = a_ix + 1
      do l = m+1, nlmax
        a_ix = a_ix + 1
        lam_fact(a_ix) = SQRT( (2 * l + 1) / real(2*l - 1,dp) * (l**2-fm2))
      end do
   end do 

  end subroutine GetLamFact


  subroutine spinalm2map(H,inlmax, alm_EB, map_QU, inspin)
    use alm_tools
    use MPIstuff
    Type (HealpixInfo) :: H

    INTEGER(I4B), INTENT(IN) :: inlmax, inspin
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_EB
    COMPLEX(SPC), INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map_QU
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: EB
    COMPLEX(SPC), DIMENSION(:), pointer :: map2
    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U

    CHARACTER(LEN=*), PARAMETER :: code = 'SPINALM2MAP'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_north_Q, b_north_U
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_south_Q, b_south_U
    INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax, spin

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(:),   ALLOCATABLE :: ringR, ringI
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    integer inf(2), a_ix, nalms
    double precision Initime
    !=======================================================================

    !     --- allocates space for arrays ---


     nsmax = H%nside
     nlmax = inlmax
     spin = inspin

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiID==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
     inf(1) = nlmax; inf(2) = spin
     call MPI_BCAST(inf,size(inf),MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
     nlmax = inf(1); spin = inf(2)
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(EB(2,nalms))
     if (H%MpiID==0) call EB2PackEB(alm_EB,EB,nlmax)
    
#ifdef MPIPIX
     call MPI_BCAST(EB,SIze(EB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiID, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1), stat = status)    
     if (status /= 0) call die_alloc(code,'map2')
#else
     map2 => map_QU 
#endif

    if (spin<1 .or. spin>3) stop 'Only spin 1 to 3 supported'

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')

    ALLOCATE(b_north_Q(0:nlmax),&
         &   b_north_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_north')

    ALLOCATE(b_south_Q(0:nlmax),&
         &   b_south_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_south')

    ALLOCATE(ringR(0:4*nsmax-1),ringI(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')


    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamFact(lam_fact, nlmax)
   if (spin ==2 ) lam_fact = lam_fact * 2 !HealPix polarization def
   

    normal_l = 0.0_dp
    do l = spin, nlmax
       fl = DBLE(l)
       if (spin==1) then 
        normal_l(l) = EB_sign*SQRT( 1 / ((fl+1.0_dp)*fl) ) 
       else if (spin==2) then
        normal_l(l) = EB_sign*SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
       else if (spin==3) then
        normal_l(l) = EB_sign*SQRT( 1 / ((fl+3.0_dp)*(fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)*(fl-2.0_dp)) ) 
       end if 
    enddo

    !     ----- set the whole map to zero ------
!    map_QU = 0.0
    !     --------------------------------------
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !     --------------------------------------------

    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)      ! 0 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
       if (ith < nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2 !cos theta
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = (-1)**spin  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac*lam_mm/OVFLOW     !  actual lambda_mm      
          a_ix = a_ix + 1
          !l=m special case
          if (m >=spin) then

              if (spin==1) then
                !swapped
                lambda_x = normal_l(m) * lam_lm * fm / sth 
                lambda_w = -lambda_x * cth
              else if (spin==2) then
                lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
                lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              else if (spin==3) then
                lambda_x = normal_l(m) * lam_lm / sth * fm*(fm-1)*(fm-2) 
                lambda_w = lambda_x * cth * ( 1 - 4*one_on_s2)
                lambda_x = lambda_x * (4*one_on_s2 - 3)             
              end if 
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

              ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
              b_n_Q =  lambda_w * EB(1,a_ix) + zi_lam_x * EB(2,a_ix)
              b_s_Q =  par_lm*(lambda_w * EB(1,a_ix) - zi_lam_x * EB(2,a_ix))

              ! - alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
              b_n_U = lambda_w * EB(2,a_ix) - zi_lam_x * EB(1,a_ix)
              b_s_U = par_lm*(lambda_w * EB(2,a_ix) + zi_lam_x * EB(1,a_ix))

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac*lam_mm/OVFLOW ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl
             a_ix = a_ix + 1
             if (l>=spin .and. corfac /= 0) then

                 if (spin==1) then
                    a_w = normal_l(l) / sth
                    lambda_x = a_w * fm  * lam_lm
                    lambda_w = a_w * (lam_fact(a_ix)*lam_lm1m - fl*cth*lam_lm) 
                 else if (spin==2) then
                     a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                     b_w =  c_on_s2 * lam_fact(a_ix)
                     a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                     lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                     lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 else if (spin==3) then
                     a_w = normal_l(l) /sth
                     b_w = (l-1)*(l-2)
                     lambda_x = a_w*fm*( (one_on_s2*(4*fm2-(12*l-8)) - 3 * b_w)*lam_lm + &
                          12*lam_fact(a_ix) * c_on_s2  * lam_lm1m)  
                     lambda_w = a_w*( (l*b_w - one_on_s2*(8*l+fm2*(4*l-12))) * cth * lam_lm  &
                        - lam_fact(a_ix)*( fl2 + fl +6 - (8+4*fm2)*one_on_s2) * lam_lm1m)

                 end if
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
                 factor_1 =  lambda_w * EB(1,a_ix)
                 factor_2 =  zi_lam_x * EB(2,a_ix) ! X is imaginary
                 b_n_Q = b_n_Q +           factor_1 + factor_2
                 b_s_Q = b_s_Q + par_lm * (factor_1 - factor_2)! X has a diff. parity

                 !- alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
                 factor_1 =   lambda_w * EB(2,a_ix) 
                 factor_2 =   zi_lam_x * EB(1,a_ix) ! X is imaginary
                 b_n_U = b_n_U +           factor_1 - factor_2
                 b_s_U = b_s_U + par_lm * (factor_1 + factor_2)! X has a diff. parity
             end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif

          enddo

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U

       enddo

       call spinring_synthesis(H,nlmax, b_north_Q, nph, ringR, kphi0,mmax_ring)
       call spinring_synthesis(H,nlmax, b_north_U, nph, ringI, kphi0,mmax_ring)
       map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) = cmplx(RingR(0:nph-1),RingI(0:nph-1))
  
       if (ith  <  2*nsmax) then
          call spinring_synthesis(H,nlmax, b_south_Q, nph, ringR, kphi0,mmax_ring)
          call spinring_synthesis(H,nlmax, b_south_U, nph, ringI, kphi0,mmax_ring)
          map2(H%istart_south(ith):H%istart_south(ith)+nph-1) = cmplx(RingR(0:nph-1),RingI(0:nph-1))

       endif

    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call HealpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l)
    DEALLOCATE(b_north_Q,b_north_U)
    DEALLOCATE(b_south_Q,b_south_U)
    DEALLOCATE(ringR,ringI)
    deallocate(EB)

#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiID
    StartTime = Getetime()
    call MPI_GATHERV(map2(H%North_Start(H%MpiID)),H%North_Size(H%MpiID),CSP_MPI, &
       map_QU,H%North_Size,H%North_Start,CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2(H%South_Start(H%MpiID)),H%South_Size(H%MpiID),CSP_MPI, &
       map_QU,H%South_Size,H%South_Start,CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
!    call MPI_REDUCE(map2,map,size(map),SP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiID, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2) 
#endif

  end subroutine spinalm2map

  subroutine map2spinalm(H, inlmax, map_QU,alm_EB, spinin, cos_theta_cut)
    use MPIStuff
    Type (HealpixInfo) :: H

    INTEGER(I4B)  :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    COMPLEX(SPC), INTENT(OUT),  DIMENSION(:,:,:) :: alm_EB
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: map_QU
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: EB, EB2
    COMPLEX(SPC), DIMENSION(:), pointer :: map2

    integer, intent(in) :: spinin
    REAL(DP),     INTENT(IN) :: cos_theta_cut

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax, spin       ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: omega_pix
    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x

    CHARACTER(LEN=*), PARAMETER :: code = 'MAP2SPINALM'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_nQ, phas_nU
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_sQ, phas_sU
    INTEGER(I4B) mmax_ring, status, par_lm, a_ix
    integer buf(2), nalms

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    LOGICAL   :: keep_it
    double precision Initime
    !=======================================================================

    nsmax = H%nside
    nlmax = inlmax
    spin = spinin
     
#ifdef MPIPIX
     if (cos_theta_cut/=-1) stop 'cos_theta_cut /= -1'
     if (H%MpiID==0) then 
      if(DebugMsgs>0) print *,code //': Sending to farm '
      call SendMessages(H,code)
      map2 => map_QU
    else
       allocate(map2(0:12*H%nside**2-1),stat = status) 
       if (status /= 0) call die_alloc(code,'map2')   
    end if

    StartTime = getetime()    
    iniTime = StartTime
    buf(1)=nlmax
    buf(2)=spin
    call MPI_BCAST(buf,2,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
    nlmax = buf(1)
    spin = buf(2)
    call MPI_SCATTERV(map_QU,H%North_Size, H%North_Start, &
       CSP_MPI, map2(H%North_Start(H%MpiID)),H%North_Size(H%MpiID),CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_SCATTERV(map_QU,H%South_Size, H%South_Start, &
       CSP_MPI, map2(H%South_Start(H%MpiID)),H%South_Size(H%MpiID),CSP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if(DebugMsgs>1) print *,code //' Scattered ',H%MpiID, GeteTime() - StartTime
#else
    map2 => map_QU
#endif

    nalms = ((nlmax+1)*(nlmax+2))/2   

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')

    ALLOCATE(phas_nQ(0:nlmax),&
         &   phas_nU(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_n')

    ALLOCATE(phas_sQ(0:nlmax),&
         &   phas_sU(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_s')

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    if (spin<1 .or. spin>3) stop 'Only spin 1 to 3 supported'

    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamfact(lam_fact, nlmax)
   if (spin==2) lam_fact = lam_fact*2

   allocate(EB(2,nalms),stat = status)
   if (status /= 0) call die_alloc(code,'EB')
   EB = 0    
       
   omega_pix = EB_sign * pi / (3 * nsmax * real(nsmax,dp))

    normal_l = 0.0_dp
    do l = spin, nlmax
       fl = DBLE(l)
       if (spin==1) then 
        normal_l(l) = omega_pix * SQRT( 1 / ((fl+1.0_dp)*fl) ) 
       else if (spin==2) then
        normal_l(l) = omega_pix * SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
       else if (spin==3) then
        normal_l(l) = omega_pix * SQRT( 1 / ((fl+3.0_dp)*(fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)*(fl-2.0_dp)) ) 
       end if 
    enddo

    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    
    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)

       phas_nQ=0; phas_sQ=0;phas_nU=0;phas_sU=0

       if (ith  <=  nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          cth = DBLE(2*nsmax-ith) * dth2
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2

       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       keep_it = (ABS(cth) > cos_theta_cut) ! part of the sky out of the symmetric cut
       if (keep_it) then
          ring(0:nph-1) = real(map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1)) 
          call spinring_analysis(H,nlmax, ring, nph, phas_nQ, kphi0, mmax_ring)
          ring(0:nph-1) = aimag(map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1)) 
          call spinring_analysis(H,nlmax, ring, nph, phas_nU, kphi0, mmax_ring)
       endif

       if (ith  <  2*nsmax .and. keep_it) then
          ring(0:nph-1) = real(map2(H%istart_south(ith):H%istart_south(ith)+nph-1)) 
          call spinring_analysis(H,nlmax, ring, nph, phas_sQ, kphi0, mmax_ring)
          ring(0:nph-1) = aimag(map2(H%istart_south(ith):H%istart_south(ith)+nph-1)) 
          call spinring_analysis(H,nlmax, ring, nph, phas_sU, kphi0, mmax_ring)
       endif

       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

       if (keep_it) then ! avoid un-necessary calculations (EH, 09-2001)
          lam_mm = sq4pi_inv
          scalem=1
          a_ix = 0
          do m = 0, mmax_ring
             fm  = DBLE(m)
             f2m = 2.0_dp * fm
             fm2 = fm * fm
             fm_on_s2 = fm * one_on_s2

             !           ---------- l = m ----------
             par_lm = (-1)**spin   ! = (-1)^(l+m+s)
             if (m  >=  1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm) < UNFLOW) then
                lam_mm=lam_mm*OVFLOW
                scalem=scalem-1
             endif

             a_ix = a_ix+1

             corfac = ScaleFactor(scalem)
             lam_lm = corfac*lam_mm/OVFLOW
        
            if (m >=spin) then
            
              if (spin==1) then
                lambda_x = normal_l(m) * lam_lm * fm / sth 
                lambda_w = -lambda_x * cth
              else if (spin==2) then
                lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
                lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              else if (spin==3) then
                lambda_x = normal_l(m) * lam_lm / sth * fm*(fm-1)*(fm-2) 
                lambda_w = lambda_x * cth * ( 1 - 4*one_on_s2)
                lambda_x = lambda_x * (4*one_on_s2 - 3)             
              end if 
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)
              
                 EB(1,a_ix) = EB(1,a_ix) &
                  &                 + lambda_w * (phas_nQ(m) + par_lm*phas_sQ(m)) &
                  &                 + zi_lam_x * (phas_nU(m) - par_lm* phas_sU(m))

                 EB(2,a_ix) = EB(2,a_ix) &
                  &                 + lambda_w * (phas_nU(m) + par_lm*phas_sU(m)) &
                  &                 - zi_lam_x * (phas_nQ(m) - par_lm*phas_sQ(m))
            
            end if

             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = H%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec
             do l = m+1, nlmax
                par_lm = - par_lm  ! = (-1)^(l+m)
                lam_lm1m=lam_lm ! actual lambda_l-1,m (useful for polarisation)
                lam_lm   = lam_2*corfac*lam_mm/OVFLOW ! actual lambda_lm (OVFLOW factors removed)
                fl  = DBLE(l)
                fl2 = fl * fl

                a_ix = a_ix + 1
             if (l>=spin .and. corfac /= 0) then
                 !Corfac=0 guarantees lam(l-1) is also v close to zero
                  
                 if (spin==1) then
                    a_w = normal_l(l) / sth
                    lambda_x = a_w * fm  * lam_lm
                    lambda_w = a_w * (lam_fact(a_ix)*lam_lm1m - fl*cth*lam_lm) 
                 else if (spin==2) then
                     a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                     b_w =  c_on_s2 * lam_fact(a_ix)
                     a_x =  2*(l-1)* cth * lam_lm
                     lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                     lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 else if (spin==3) then
                     a_w = normal_l(l) /sth
                     b_w = (l-1)*(l-2)
                     lambda_x = a_w*fm*( (one_on_s2*(4*fm2-(12*l-8)) - 3 * b_w)*lam_lm + &
                          12*lam_fact(a_ix) * c_on_s2  * lam_lm1m)  
                     lambda_w = a_w*( (l*b_w - one_on_s2*(8*l+fm2*(4*l-12))) * cth * lam_lm  &
                        - lam_fact(a_ix)*( fl2 + fl +6 - (8+4*fm2)*one_on_s2) * lam_lm1m)

                 end if

                zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 EB(1,a_ix) = EB(1,a_ix) &
                     &          + lambda_w * (phas_nQ(m) + par_lm*phas_sQ(m)) &
                     &          + zi_lam_x * (phas_nU(m) - par_lm*phas_sU(m))
                 EB(2,a_ix) = EB(2,a_ix) &
                     &         +  lambda_w * (phas_nU(m) + par_lm*phas_sU(m)) &
                     &         - zi_lam_x  * (phas_nQ(m) - par_lm*phas_sQ(m))

              end if ! l allowed by spin or zero

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_2)  >  OVFLOW) then
                   lam_0=lam_0/OVFLOW
                   lam_1=lam_1/OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = ScaleFactor(scalem+scalel)
                elseif (abs(lam_2)  <  UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = ScaleFactor(scalem+scalel)
                endif

             enddo ! loop on l
          enddo ! loop on m
       endif ! test on cut sky
    enddo ! loop on theta

    !     --------------------
    !     free memory and exit
    !     --------------------
    call HealpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l)
    DEALLOCATE(phas_nQ,phas_nU)
    DEALLOCATE(phas_sQ,phas_sU)
    DEALLOCATE(ring)
#ifdef MPIPIX
    if (H%MpiID>0) deallocate(map2)
    allocate(EB2(2,nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'EB2')
    StartTime = Getetime()
    call MPI_REDUCE(EB,EB2,size(EB),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    if (DebugMsgs>1) print *,code//' done reduce ', H%MpiID, GeteTime() -StartTime
    deallocate(EB)
    if (H%MpiID == 0) call PackEB2EB(EB2,alm_EB, nlmax)
    deallocate(EB2)
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code //' Time: ',GeteTime() - IniTime
#else
    call PackEB2EB(EB,alm_EB, nlmax)
    deallocate(EB)
#endif

  END subroutine map2spinalm

  subroutine spinring_analysis(H, nlmax, datain,nph,dataout,kphi0,mmax_ring)
    !=======================================================================
    !     ring_analysis
    !       called by map2alm
    !       calls     real_fft
    !
    !     integrates (data * phi-dependence-of-Ylm) over phi
    !     --> function of m can be computed by FFT
    !     with  0<= m <= npoints/2 (: Nyquist)
    !     because the data is real the negative m are the conjugate of the 
    !     positive ones
    !=======================================================================
    Type (HealpixInfo) :: H

    INTEGER(I4B) :: nsmax
    INTEGER(I4B), INTENT(IN) :: nlmax
    INTEGER(I4B), INTENT(IN) :: mmax_ring
    INTEGER(I4B), INTENT(IN) :: nph, kphi0

    REAL(DP),     DIMENSION(0:nph-1), INTENT(IN)  :: datain
    COMPLEX(DPC), DIMENSION(0:nlmax), INTENT(OUT) :: dataout
    INTEGER(I4B) :: i,m,im_max,ksign
    REAL(DP) :: phi0
    REAL(DP), DIMENSION(0:nph-1) :: data
    INTEGER(I4B) :: status

    !-----------------------------------------------------------------------

    call HealpixInfo_GetTrig(H,nph)

    nsmax = H%nside
 
    ksign = - 1
    data=0.
    data(0:nph-1)=datain(0:nph-1)

    call real_fft(data, backward=.false.)

    im_max = MIN(nph/2,mmax_ring)
    dataout(0)=CMPLX(data(0),0.0_dp,kind=DP)

    do i = 1, im_max*2-3, 2
       dataout((i+1)/2) = CMPLX( data(i), data(i+1),kind= DP) 
    enddo

    if(im_max==nph/2) then
       dataout(im_max)= CMPLX( data(nph-1),0.0_dp,kind=DP)
    else
       dataout(im_max)= CMPLX( data(2*im_max-1),data(2*im_max),kind=DP)
    endif

    if(im_max==mmax_ring) goto 1000

    do i =  im_max+1,min(nph,mmax_ring)
       dataout(i) = conjg(dataout(2*im_max-i) )
    end do

    if(min(nph,mmax_ring)==mmax_ring) goto 1000

    do i =  2*im_max+1,mmax_ring
       dataout(i) = dataout(mod(i,2*im_max)) 
    end do

1000 continue

    if(kphi0==1)then
       do i =  0,mmax_ring
          m = ksign*i
          dataout(i)=dataout(i)* CONJG(H%trig(-m))
       enddo
    end if


  END subroutine spinring_analysis

    !=======================================================================
  subroutine scalalm2map(H, inlmax, alm, map)
    !=======================================================================
    !     computes a map form its alm for the HEALPIX pixelisation
    !      for the Temperature field
    !     map(theta,phi) = sum_l_m a_lm Y_lm(theta,phi)
    !                    = sum_m {e^(i*m*phi) sum_l a_lm*lambda_lm(theta)}
    !
    !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
    !
    !     * the recurrence of Ylm is the standard one (cf Num Rec)
    !     * the sum over m is done by FFT
    !
    !              -------------------------------
    !          precomputes the Lambda_lm recurrence factor 
    !      and is therefore ~50% faster than previous versions
    !              -------------------------------
    !
    !=======================================================================
    use MPIstuff
    Type (HealpixInfo) :: H

    INTEGER(I4B) :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm
    REAL(SP),     INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map

    REAL(SP),     DIMENSION(:), pointer :: map2
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax          ! alm related
    INTEGER(I4B) :: nph, kphi0                         ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2
    REAL(DP) :: f2m, fm2, fl2, corfac
    COMPLEX(DPC) :: b_n, b_s, factor

    CHARACTER(LEN=*), PARAMETER :: code = 'SCALALM2MAP'
    COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
    INTEGER(I4B) :: mmax_ring, status, par_lm
    integer nalms, a_ix

    REAL(SP), DIMENSION(0:4*H%nside-1) :: ring
    double precision Initime
    !=======================================================================
  
     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiID==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    call MPI_BCAST(nlmax,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms))
     if (H%MpiID==0) call Alm2PackAlm(alm,alm2,nlmax)
    
#ifdef MPIPIX 
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 

     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiID, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1), stat = status)    
     if (status /= 0) call die_alloc(code,'map2')
#else
     map2 => map 
#endif

    call HealpixInitRecfac(H,nlmax)
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    do ith =H%ith_start(H%MpiID), H%ith_end(H%MpiID)  
       !        cos(theta) in the pixelisation scheme

       if (ith.lt.nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)

       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       lam_mm = sq4pi_inv ! lambda_00
       scalem=1
       a_ix = 0
       do m = 0, mmax_ring
          f2m = 2.0_dp * m

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m)
          if (m >= 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)
  
          lam_lm = lam_mm*corfac/OVFLOW
          a_ix = a_ix + 1
          b_n = lam_lm * alm2(a_ix)
          b_s = b_n

          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp 
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)

             lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors 
             a_ix = a_ix + 1
             factor = lam_lm * alm2(a_ix)
             b_n = b_n +          factor
             b_s = b_s + par_lm * factor

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2) > OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2) .lt. UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif
          enddo

          b_north(m) = b_n
          b_south(m) = b_s

       enddo

       call spinring_synthesis(H,nlmax,b_north,nph,ring,kphi0,mmax_ring)   ! north hemisph. + equator
       map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) = ring(0:nph-1)
       
       if (ith < 2*nsmax) then
          call spinring_synthesis(H,nlmax,b_south,nph,ring,kphi0,mmax_ring) ! south hemisph. w/o equat
          map2(H%istart_south(ith):H%istart_south(ith)+nph-1) = ring(0:nph-1)
       endif

    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixFreeRecFac(H)
    deallocate(alm2)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code //' Gather ',H%MpiID
    
    StartTime = Getetime()
    call MPI_GATHERV(map2(H%North_Start(H%MpiID)),H%North_Size(H%MpiID),SP_MPI, &
       map,H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2(H%South_Start(H%MpiID)),H%South_Size(H%MpiID),SP_MPI, &
       map,H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
!    call MPI_REDUCE(map2,map,size(map),SP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    if (DebugMsgs>1) print *,code //' Done Gather ',H%MpiID, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code //' Time :', GeteTime()-IniTime
    deallocate(map2)

#endif
  end subroutine scalalm2map


  subroutine map2scalalm(H,inlmax, map, alm, cos_theta_cut)
    !=======================================================================
    !     computes the a(l,m) from a map for the HEALPIX pixelisation
    !      for the Temperature field
    !     a(l,m) = int T(theta,phi) Y_lm(theta,phi)^* dtheta*dphi
    !            = int dtheta lambda_lm(theta) 
    !                  * int dphi T(theta,phi) e^(-i*m*phi)
    !
    !     where Y_lm(theta,phi) = lambda(theta) * e^(i*m*phi)
    !
    !     * the recurrence of Ylm is the standard one (cf Num Rec)
    !     * the integral over phi is done by FFT
    !
    !     cos_theta_cut (>0) is the cosine of the 
    !     symmetric cut applied to the sky
    !     if it is <= 0 no cut is applied
    !
    !     NB : these a(l,m) have to be multiplied by the pixel solid angle
    !      to give the correct coefficients
    !
    !             -------------------------------
    !         precomputes the Lambda_lm recurrence factor
    !      and is therefore ~50% faster than previous versions
    !     the multiplication by omega_pix is done in the routine
    !             -------------------------------
    !
    !=======================================================================
    use MPIStuff
    Type (HealpixInfo) :: H
    INTEGER(I4B)  :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    REAL(SP),     INTENT(IN),  DIMENSION(0:12*H%nside**2-1), target :: map
    COMPLEX(SPC), INTENT(OUT), DIMENSION(:,:,:) :: alm
    COMPLEX(SPC),   DIMENSION(:), allocatable :: alm2,alm3
    REAL(SP),     DIMENSION(:), pointer :: map2
    
    REAL(DP),     INTENT(IN) :: cos_theta_cut

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax, a_ix       ! alm related
    INTEGER(I4B) :: nph, kphi0   ! map related

    REAL(DP) :: omega_pix
    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2
    REAL(DP) :: f2m, fm2, fl2, corfac

    CHARACTER(LEN=*), PARAMETER :: code = 'MAP2SCALALM'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_n
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_s
    INTEGER(I4B) :: mmax_ring, status, par_lm, nalms

    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    LOGICAL   :: keep_it
    double precision Initime
    !=======================================================================


    nsmax = H%nside
    nlmax = inlmax
     
#ifdef MPIPIX
     if (H%MpiID==0) then 
      if(DebugMsgs>1) print *,code //': Sending to farm '
      IniTime = GeteTime()
      call SendMessages(H,code)
      map2 => map
    else
       allocate(map2(0:12*H%nside**2-1))    
    end if

    StartTime = getetime()    
    call MPI_BCAST(nlmax,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 

    call MPI_SCATTERV(map,H%North_Size, H%North_Start, &
       SP_MPI, map2(H%North_Start(H%MpiID)),H%North_Size(H%MpiID),SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_SCATTERV(map,H%South_Size, H%South_Start, &
       SP_MPI, map2(H%South_Start(H%MpiID)),H%South_Size(H%MpiID),SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if(DebugMsgs>1) print *,code //' Scattered ',H%MpiID, GeteTime() - StartTime
!   call MPI_BCAST(map,SIze(Map),SP_MPI, 0, MPI_COMM_WORLD, ierr) 
#else
    map2 => map
#endif


    !     --- allocates space for arrays ---

    ALLOCATE(phas_n(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_n')

    ALLOCATE(phas_s(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_s')

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

    call HealpixInitRecfac(H,nlmax)

    nalms = ((nlmax+1)*(nlmax+2))/2   
    allocate(alm2(nalms), stat=status)
    if (status /= 0) call die_alloc(code,'alm2')

    alm2 = 0

    !     -------------------------------------------

    omega_pix = pi / (3.0_dp * nsmax * nsmax)

    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)

       phas_n(0:nlmax) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)   ! North    m >= 0
       phas_s(0:nlmax) = CMPLX(0.0_dp, 0.0_dp, KIND=DP)   ! South    m >= 0

       if (ith .le. nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          cth = DBLE(2*nsmax-ith) * dth2
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif


       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       keep_it = (ABS(cth).gt.cos_theta_cut) ! part of the sky out of the symmetric cut


       if (keep_it) then
          ring(0:nph-1) = map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1) * H%w8ring_TQU(ith,1)
          call spinring_analysis(H,nlmax, ring, nph, phas_n, kphi0, mmax_ring)
  
       if (ith .lt. 2*nsmax) then
          ring(0:nph-1) = map2(H%istart_south(ith):H%istart_south(ith)+nph-1) * H%w8ring_TQU(ith,1)
          call spinring_analysis(H,nlmax, ring, nph, phas_s, kphi0,mmax_ring)
       endif

       endif


       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

       if (keep_it) then ! avoid un-necessary calculations (EH, 09-2001)

          lam_mm = sq4pi_inv * omega_pix ! lambda_00 * norm
          scalem=1 
          a_ix = 0
          do m = 0, mmax_ring
             f2m = 2.0_dp * m
         
             !           ---------- l = m ----------
             par_lm = 1  ! = (-1)^(l+m)
             if (m .ge. 1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm).lt.UNFLOW) then
                lam_mm=lam_mm*OVFLOW
                scalem=scalem-1
             endif
            corfac = ScaleFactor(scalem)
  
             lam_lm = corfac*lam_mm/OVFLOW
             a_ix = a_ix + 1
             alm2(a_ix) = alm2(a_ix) + lam_lm * (phas_n(m) + phas_s(m))
             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = H%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec
             do l = m+1, nlmax
                par_lm = - par_lm  ! = (-1)^(l+m)
                
   
                lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors
                a_ix = a_ix+1
                alm2(a_ix) = alm2(a_ix) + lam_lm * (phas_n(m) + par_lm*phas_s(m))

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_2) .gt. OVFLOW) then
                   lam_0=lam_0/OVFLOW
                   lam_1=lam_1/OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = ScaleFactor(scalem+scalel)
                elseif (abs(lam_2) .lt. UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = ScaleFactor(scalem+scalel)
                endif
             enddo ! loop on l

          enddo ! loop on m
       endif ! test on cut sky
    enddo ! loop on theta

    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(phas_n)
    DEALLOCATE(phas_s)
    DEALLOCATE(ring)
    call HealpixFreeRecFac(H)
    if (H%MpiID>0) deallocate(map2)
#ifdef MPIPIX
    allocate(alm3(nalms), stat = status)
    if (status /= 0) call die_alloc(code,'alm3')

    StartTime = Getetime()
    call MPI_REDUCE(alm2,alm3,size(alm2),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    if (DebugMsgs>1) print *,code //': done reduce ', H%MpiID, GeteTime() -StartTime
    deallocate(alm2)
    if (H%MpiID == 0) then
      call PackAlm2Alm(alm3, alm, nlmax)
      print *,code // ' Time: ', GeteTime() - IniTime
    end if
    deallocate(alm3)
#else
    call PackAlm2Alm(alm2, alm, nlmax)
    deallocate(alm2)
#endif

  END subroutine map2scalalm


    !=======================================================================
  subroutine scalalm2LensedMap(H, inlmax, alm, grad_phi_map, map)
   !Get lensed map without using series expansion.
   !No FFT so does not scale with ln, so slow
    use MPIstuff
    Type (HealpixInfo) :: H

    INTEGER(I4B) :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm
    REAL(SP),     INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
    REAL(SP),     DIMENSION(:), pointer :: map2

    COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel, nlmax          ! alm related
    INTEGER(I4B) :: nph, kphi0                         ! map related

    REAL(DP) :: cth0,sth0,cth, sth, dth1, dth2, dst1
    real(DP) :: phi
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2
    REAL(DP) :: f2m, fm2, fl2, corfac
    REAL(DP) :: grad_len, sinc_grad_len
    COMPLEX(DPC) :: b_n, factor

    CHARACTER(LEN=*), PARAMETER :: code = 'SCALALM2LENSEDMAP'
    INTEGER(I4B) :: mmax_ring, status, par_lm
    integer nalms, a_ix, ring_ix
    double precision Initime
    !=======================================================================
  
     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiID==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    call MPI_BCAST(nlmax,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 

#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms))
     if (H%MpiID==0) then
       call Alm2PackAlm(alm,alm2,nlmax)
       grad_phi => grad_phi_map
     else
       allocate(grad_phi(0:12*H%nside**2-1))    
     end if

    
#ifdef MPIPIX 
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 

     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiID, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1), stat = status)    
     if (status /= 0) call die_alloc(code,'map2')
#else
     map2 => map 
#endif

    map2 = 0

    call HealpixInitRecfac(H,nlmax)
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    do ith =H%ith_start(H%MpiID), H%ith_end(H%MpiID)  
       !        cos(theta) in the pixelisation scheme


       if (ith.lt.nsmax) then  ! polar cap (north)
          cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth0 = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
       endif
    

       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth0)))
       else
        mmax_ring = nlmax
       end if

       do ring_ix = H%istart_north(ith-1),H%istart_north(ith-1)+nph-1
       
!       cth = cos(theta + real(grad_phi(ring_ix)))
!       sth = sin(theta + real(grad_phi(ring_ix))) 
!       phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph + aimag(grad_phi(ring_ix))/sth0
        phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph

        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
            sinc_grad_len = sin(grad_len)/grad_len
            cth =  cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
            sth = sqrt((1._dp-cth)*(1._dp+cth))
            if (sth > 1e-10_dp) then
             phi = phi + asin(aimag(grad_phi(ring_ix))*sinc_grad_len/ sth ) 
            end if
        else
         cth=cth0
         sth=sth0
        end if

       if (kphi0==1) phi=phi + pi/nph

       lam_mm = sq4pi_inv ! lambda_00
       scalem=1
       a_ix = 0
       do m = 0, mmax_ring
          f2m = 2.0_dp * m

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m)
          if (m >= 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)
  
          lam_lm = lam_mm*corfac/OVFLOW
          a_ix = a_ix + 1
          b_n = lam_lm * alm2(a_ix)
          
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp 
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)

             lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors 
             a_ix = a_ix + 1
             factor = lam_lm * alm2(a_ix)
             b_n = b_n + factor

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2) > OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2) .lt. UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif
          enddo

          if (m==0) then
           map2(ring_ix) =  real(b_n)  
          else
           map2(ring_ix) =  map2(ring_ix) + 2*real(b_n*cmplx(cos(m*phi),sin(m*phi)))  
          end if
       enddo

      enddo !ring_ix (phi)

      if (ith < 2*nsmax) then
       do ring_ix = H%istart_south(ith),H%istart_south(ith)+nph-1
       
       !Think of spherical triangle with points 0, (theta0,phi0), (theta, phi)
        phi = (ring_ix-H%istart_south(ith))*2*pi/nph
        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
         sinc_grad_len = sin(grad_len)/grad_len
         cth =  -cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
         sth = sqrt((1._dp-cth)*(1._dp+cth))
         if (sth > 1e-10_dp) then
          phi = phi +asin(aimag(grad_phi(ring_ix))*sinc_grad_len/ sth ) 
         end if

       else
        cth=-cth0
        sth=sth0
       end if
       !cth = cos(pi - theta + real(grad_phi(ring_ix)))
       !sth = sin(pi - theta + real(grad_phi(ring_ix))) 

       if (kphi0==1) phi=phi + pi/nph

       lam_mm = sq4pi_inv ! lambda_00
       scalem=1
       a_ix = 0
       do m = 0, mmax_ring
          f2m = 2.0_dp * m

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m)
          if (m >= 1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif
          if (abs(lam_mm).lt.UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)
  
          lam_lm = lam_mm*corfac/OVFLOW
          a_ix = a_ix + 1
          b_n = lam_lm * alm2(a_ix)
          
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp 
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)

             lam_lm = lam_2*corfac*lam_mm/OVFLOW ! Remove OVFLOW-factors 
             a_ix = a_ix + 1
             factor = lam_lm * alm2(a_ix)
             b_n = b_n + factor

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2) > OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2) .lt. UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif
          enddo

          if (m==0) then
           map2(ring_ix) =  real(b_n)  
          else
           map2(ring_ix) =  map2(ring_ix) + 2*real(b_n*cmplx(cos(m*phi),sin(m*phi)))  
          end if
       enddo
       
      enddo !ring_ix (phi)
     end if !not middle theta

    enddo    ! loop on cos(theta)

    call healpixFreeRecFac(H)
    deallocate(alm2)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code //' Gather ',H%MpiID
    StartTime = Getetime()
    call MPI_GATHERV(map2(H%North_Start(H%MpiID)),H%North_Size(H%MpiID),SP_MPI, &
       map,H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2(H%South_Start(H%MpiID)),H%South_Size(H%MpiID),SP_MPI, &
       map,H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if (DebugMsgs>1) print *,code //' Done Gather ',H%MpiID, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code //' Time :', GeteTime()-IniTime
    if (H%MpiID/=0) deallocate(grad_phi) 
    deallocate(map2)

#endif
  end subroutine scalalm2LensedMap


  subroutine alm2Lensedmap(H,inlmax, alm_TEB, grad_phi_map,map_TQU)
    use alm_tools
    use MPIstuff
    Type (HealpixInfo) :: H

    INTEGER(I4B), INTENT(IN) :: inlmax
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_TEB
    REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1,3), target :: map_TQU
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: TEB

    REAL(SP),     DIMENSION(:,:), pointer :: map2
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map

    COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi

    INTEGER(I4B) :: i, l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related
    REAL(DP) :: cth0,sth0

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n, b_Q, b_U, exp_m_phi, gammfac
    REAL(DP) :: Re, Im, grad_len, sinc_grad_len, phi, gamma, DeltaPhi, alpha
    INTEGER(I4B) :: ring_ix  

    CHARACTER(LEN=*), PARAMETER :: code = 'ALM2LENSEDMAP'
    INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax
    integer, parameter :: spin=2

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    integer inf(2), a_ix, nalms
    double precision Initime
    real tmp1,tmp2,tmp3
    !=======================================================================

    !     --- allocates space for arrays ---


     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiID==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
     call MPI_BCAST(nlmax,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
#endif

     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(TEB(3,nalms))
     if (H%MpiID==0) then
        call TEB2PackTEB(alm_TEB,TEB,nlmax)
       grad_phi => grad_phi_map
     else
        allocate(grad_phi(0:12*H%nside**2-1))    
     end if

#ifdef MPIPIX
     call MPI_BCAST(TEB,SIze(TEB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiID, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1,3), stat = status)    
     if (status /= 0) call die_alloc(code,'map2')
#else
     map2 => map_TQU 
#endif


    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')


   call HealpixInitRecfac(H,nlmax)
   call GetLamFact(lam_fact, nlmax)
   lam_fact = lam_fact * 2 !HealPix polarization def

    normal_l = 0.0_dp
    do l = 2, nlmax
       fl = DBLE(l)
       normal_l(l) = EB_sign * SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo


    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)      ! 0 <= cos theta < 1

       if (ith < nsmax) then  ! polar cap (north)
          cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth0 = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
       endif
   
       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth0)))
       else
        mmax_ring = nlmax
       end if
   
       do ring_ix = H%istart_north(ith-1),H%istart_north(ith-1)+nph-1

        phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph

        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
            sinc_grad_len = sin(grad_len)/grad_len
            cth =  cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
            sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
            DeltaPhi = asin(aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )
            phi = phi  +  DeltaPhi 
        else
         DeltaPhi = 0
         cth=cth0
         sth=sth0
        endif
       
        one_on_s2 = 1.0_dp / sth**2 
        c_on_s2 = cth * one_on_s2

        if (kphi0==1) phi=phi + pi/nph

       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac*lam_mm/OVFLOW     !  actual lambda_mm      
          a_ix = a_ix + 1

           b_n = lam_lm * TEB(1,a_ix)
    
           !l=m special case
           if (m >=2) then

              lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
              lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)
              b_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
              b_U = lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
  
          else
             b_Q=0
             b_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac*lam_mm/OVFLOW ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl
             a_ix = a_ix + 1

             b_n = b_n + lam_lm * TEB(1,a_ix)

             if (l>=2 .and. corfac /= 0) then

                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 factor_1 =  lambda_w * TEB(2,a_ix)
                 factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
                 b_Q = b_Q +           factor_1 + factor_2
  
                 factor_1 =   lambda_w * TEB(3,a_ix) 
                 factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
                 b_U = b_U +           factor_1 - factor_2
              end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif

          enddo

          if (m==0) then
           map2(ring_ix,1) =  real(b_n)
           map2(ring_ix,2) =  real(b_Q)
           map2(ring_ix,3) =  real(b_U)  
          else
           exp_m_phi = cmplx(cos(m*phi),sin(m*phi)) 
           map2(ring_ix,1) =  map2(ring_ix,1) + 2*real(b_n*exp_m_phi) 
           map2(ring_ix,2) =  map2(ring_ix,2) + 2*real(b_Q*exp_m_phi)
           map2(ring_ix,3) =  map2(ring_ix,3) + 2*real(b_U*exp_m_phi)  
          end if

       enddo

     !Put in factor for change of co-ordinate axes between deflected and original point
          if (grad_len > 1e-20_dp) then
             Re = real(grad_phi(ring_ix))
             Im = aimag(grad_phi(ring_ix))
             gamma = grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
             if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
             gamma = Im/gamma
             !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
             gammfac = cmplx(map2(ring_ix,2),map2(ring_ix,3))* &
               cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
                    2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))             
             map2(ring_ix,2) = real(gammfac)
             map2(ring_ix,3) = aimag(gammfac)
         end if

      enddo !ring_ix (phi)

      if (ith < 2*nsmax) then
       do ring_ix = H%istart_south(ith),H%istart_south(ith)+nph-1
    
        phi = (ring_ix-H%istart_south(ith))*2*pi/nph
        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
            sinc_grad_len = sin(grad_len)/grad_len
            cth =  -cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
            sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
            DeltaPhi = asin(aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )
            phi = phi + DeltaPhi 
        else
         DeltaPhi = 0
         cth=-cth0
         sth=sth0
        endif
        if (kphi0==1) phi=phi + pi/nph
        one_on_s2 = 1.0_dp / sth**2 
        c_on_s2 = cth * one_on_s2


        lam_mm = sq4pi_inv ! lamda_00
        scalem=1

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac*lam_mm/OVFLOW     !  actual lambda_mm      
          a_ix = a_ix + 1
          b_n = lam_lm * TEB(1,a_ix)
         
          !l=m special case
          if (m >=2) then

              lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
              lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

              b_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
              b_U =  lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
  
          else
             b_Q=0
             b_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac*lam_mm/OVFLOW ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl
             a_ix = a_ix + 1

             b_n = b_n + lam_lm * TEB(1,a_ix)

             if (l>=2 .and. corfac /= 0) then

                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 factor_1 =  lambda_w * TEB(2,a_ix)
                 factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
                 b_Q = b_Q +           factor_1 + factor_2
  
                 factor_1 =   lambda_w * TEB(3,a_ix) 
                 factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
                 b_U = b_U +           factor_1 - factor_2
               end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif

          enddo

          if (m==0) then
           map2(ring_ix,1) =  real(b_n)
           map2(ring_ix,2) =  real(b_Q)
           map2(ring_ix,3) =  real(b_U)  
          else
           exp_m_phi = cmplx(cos(m*phi),sin(m*phi)) 
           map2(ring_ix,1) =  map2(ring_ix,1) + 2*real(b_n*exp_m_phi) 
           map2(ring_ix,2) =  map2(ring_ix,2) + 2*real(b_Q*exp_m_phi)
           map2(ring_ix,3) =  map2(ring_ix,3) + 2*real(b_U*exp_m_phi)  
          end if
       end do

     !Put in factor for change of co-ordinate axes between deflected and original point
          if (grad_len > 1e-20_dp) then
             Re = real(grad_phi(ring_ix))
             Im = aimag(grad_phi(ring_ix))
             gamma = grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
             if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
             gamma = Im/gamma
             !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
             gammfac = cmplx(map2(ring_ix,2),map2(ring_ix,3))* &
               cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
                    2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))             
             map2(ring_ix,2) = real(gammfac)
             map2(ring_ix,3) = aimag(gammfac)
         end if

      end do
      end if  
       
    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call HealpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l)
    deallocate(TEB)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiID
    StartTime = Getetime()
    do i=1,3
     call MPI_GATHERV(map2(H%North_Start(H%MpiID),i),H%North_Size(H%MpiID),SP_MPI, &
       map_TQU(:,i),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
     call MPI_GATHERV(map2(H%South_Start(H%MpiID),i),H%South_Size(H%MpiID),SP_MPI, &
       map_TQU(:,i),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    end do
    if (H%MpiID/=0) deallocate(grad_phi)
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiID, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2) 
#endif

  end subroutine alm2Lensedmap


 !=======================================================================
  subroutine ang2pix_ring8(nside, costheta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map 
    !     resolution parameter nside
    !=======================================================================
    !AL: Uses full I8B integer pixel range, takes in cos(theta) rather than theta
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I8B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN) ::  costheta, phi

    INTEGER(KIND=I8B) ::  nl4, jp, jm
    REAL(KIND=DP) ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER(KIND=I8B) ::  ir, ip, kshift

    z = costheta
    za = ABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)

    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(.5_dp+tt)
       temp2 = nside*.75_dp*z
       jp = int(temp1-temp2) ! index of  ascending edge line 
       jm = int(temp1+temp2) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 1 - modulo(ir,2_I8B) ! kshift=1 if ir even, 0 otherwise

       nl4 = 4*nside
       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*int(nside-1,I8B) + nl4*(ir-1) + ip 

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       tmp = nside * SQRT( 3.0_dp*(1.0_dp - za) )

       jp = INT(tp          * tmp ) ! increasing edge line index
       jm = INT((1.0_dp - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir )     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._dp) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 12*int(nside,I8B)**2 - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine ang2pix_ring8


  subroutine alm2LensedmapInterp(H,inlmax, alm_TEB, grad_phi_map, map_TQU, nside_factor)
  !This routine is designed to be used over a cluster of say 30+ nodes, at least
  !for high resolution maps. Without a cluster memory use will probably be > 2GB.
  !It currently just re-maps pixels, without any interpolation or accounting for pixel shape etc
  !nside_fac=4 at nside=1024 is probably sufficient before Planck, nside_fac=8 for Planck at 0.5%.
    use MPIstuff
    Type (HealpixInfo) :: H, H_res

    INTEGER(I4B), INTENT(IN) :: inlmax
    INTEGER(I4B), INTENT(IN), optional :: nside_factor

    integer nsmax
    integer  :: nside_fac = 8
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_TEB
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
    REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1,3), target :: map_TQU
    REAL(SP), DIMENSION(:,:), pointer :: map2N, map2S
    REAL(SP), DIMENSION(:,:), pointer :: high_resN,high_resS
    COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: TEB

    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U
    COMPLEX(DPC) :: b_n, b_s

    CHARACTER(LEN=*), PARAMETER :: code = 'ALM2LENSEDMAPINTERP'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_north_Q, b_north_U
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_south_Q, b_south_U
    COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
    INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(:), allocatable ::  ring
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    integer(I4B) i, a_ix, nalms, border_inc, high_th_start,high_th_end, high_nside
    integer(I8B) ipix, ring_ix, resN_start, resS_start
    COMPLEX(DPC) :: gammfac
    REAL(DP) :: grad_len, sinc_grad_len, phi, gamma, DeltaPhi, alpha
    REAL(DP) :: Re, Im, cth0, sth0, theta
    integer :: buf(3)
    double precision Initime
    !=======================================================================

    !     --- allocates space for arrays ---

     nsmax = H%nside
     nlmax = inlmax
     if (present(nside_factor)) nside_fac = nside_factor

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiID==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    buf(1) = nlmax
    buf(2) = nside_fac
#endif

    if (H%MpiID ==0) then
     high_nside = nsmax*nside_fac
 !border_inc=number of high-res pixels we need to go outside zero-lensing border
 !1.18 is approx 3Pi/8 which is the large-n_side vertical overdensity of pixels near the 
 !equator relative to the average. Could speed by putting in position dependent
 !factor. Corrected AL: 30 Sept 2004
     border_inc = int(maxval(abs(real(grad_phi_map)))/ PI *4*high_nside * 1.18) + 1  
    end if

#ifdef MPIPIX
    buf(3) = border_inc
    call MPI_BCAST(buf,3,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
    nlmax = buf(1)
    nside_fac = buf(2)
    border_inc = buf(3)
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(TEB(3,nalms))
     if (H%MpiID==0) then
        call TEB2PackTEB(alm_TEB,TEB,nlmax)
     end if
     
#ifdef MPIPIX
     call MPI_BCAST(TEB,SIze(TEB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) then
      print *,code //': Got alm ',H%MpiID, GeteTime() - StartTime
      StartTime = geteTime()
     end if
#endif

    high_nside = H%nside*nside_fac
    call HealpixInitTrig(H_res,high_nside,nlmax)
    H_res%MpiID = 1
 
    high_th_start = max((H%ith_start(H%MpiID)-1) * nside_fac - border_inc, 1)
    high_th_end  =  min(H%ith_end(H%MpiID)  * nside_fac + border_inc, 2*high_nside)
    if (high_th_end  < high_nside) then  ! polar cap (north)
              nph = 4*high_th_end
           else                   
              nph = 4*high_nside
    endif
    resN_start = H_res%istart_north(high_th_start-1) 
    allocate(high_resN(0:H_res%istart_north(high_th_end-1)+nph-1 - resN_start,3))
    if (high_th_start  < high_nside) then  ! polar cap (north)
              nph = 4*high_th_start
           else                   
              nph = 4*high_nside
    endif
    resS_start = H_res%istart_south(min(high_th_end,2*high_nside-1)) 
    allocate(high_resS(0:H_res%istart_south(high_th_start)+nph-1 - resS_start,3))    

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')

    ALLOCATE(b_north_Q(0:nlmax),&
         &   b_north_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_north')

    ALLOCATE(b_south_Q(0:nlmax),&
         &   b_south_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_south')

    ALLOCATE(ring(0:4*H%nside*nside_fac-1), stat = status)
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamFact(lam_fact, nlmax)
   lam_fact = lam_fact * 2 !HealPix polarization def

    normal_l = 0.0_dp
    do l = 2, nlmax
       fl = DBLE(l)
       normal_l(l) = EB_sign*SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo
 
    dth1 = 1.0_dp / (3.0_dp*DBLE(high_nside)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(high_nside))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(high_nside) )

    do ith = high_th_start, high_th_end      ! 0 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
    
       if (ith < high_nside) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*high_nside-ith) * dth2 !cos theta
          nph = 4*high_nside
          kphi0 = MOD(ith+1-high_nside,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac*lam_mm/OVFLOW     !  actual lambda_mm      

          a_ix = a_ix + 1

          b_n = lam_lm * TEB(1,a_ix)
          b_s = b_n

          !l=m special case
          if (m >=2) then
              lambda_w = - 2.0_dp *(normal_l(m) * lam_lm * (fm - fm2) ) * ( one_on_s2 - 0.5_dp )
              lambda_x =  ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

              b_n_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
              b_s_Q =  par_lm*(lambda_w * TEB(2,a_ix) - zi_lam_x * TEB(3,a_ix))

              b_n_U = lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
              b_s_U = par_lm*(lambda_w * TEB(3,a_ix) + zi_lam_x * TEB(2,a_ix))

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac*lam_mm/OVFLOW ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl
             a_ix = a_ix + 1

             factor = lam_lm * TEB(1,a_ix)
             b_n = b_n +          factor
             b_s = b_s + par_lm * factor

             if (l>=2 .and. corfac /= 0) then

                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
                 factor_1 =  lambda_w * TEB(2,a_ix)
                 factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
                 b_n_Q = b_n_Q +           factor_1 + factor_2
                 b_s_Q = b_s_Q + par_lm * (factor_1 - factor_2)! X has a diff. parity

                 !- alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
                 factor_1 =   lambda_w * TEB(3,a_ix) 
                 factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
                 b_n_U = b_n_U +           factor_1 - factor_2
                 b_s_U = b_s_U + par_lm * (factor_1 + factor_2)! X has a diff. parity
             end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif

          enddo

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U
          b_north(m) = b_n
          b_south(m) = b_s

       enddo

       call spinring_synthesis(H_res,nlmax,b_north,nph,ring,kphi0,mmax_ring)   
       high_resN(H_res%istart_north(ith-1)-resN_start:H_res%istart_north(ith-1)+nph-1-resN_start,1) = ring(0:nph-1)
       call spinring_synthesis(H_res,nlmax, b_north_Q, nph, ring, kphi0,mmax_ring)
       high_resN(H_res%istart_north(ith-1)-resN_start:H_res%istart_north(ith-1)+nph-1-resN_start,2) = ring(0:nph-1)
       call spinring_synthesis(H_res,nlmax, b_north_U, nph, ring, kphi0,mmax_ring)
       high_resN(H_res%istart_north(ith-1)-resN_start:H_res%istart_north(ith-1)+nph-1-resN_start,3) = ring(0:nph-1)
  
       if (ith  <  2*high_nside) then
          call spinring_synthesis(H_res,nlmax, b_south, nph, ring, kphi0,mmax_ring)
          high_resS(H_res%istart_south(ith)-resS_start:H_res%istart_south(ith)+nph-1-resS_start,1) = ring(0:nph-1)
          call spinring_synthesis(H_res,nlmax, b_south_Q, nph, ring, kphi0,mmax_ring)
          high_resS(H_res%istart_south(ith)-resS_start:H_res%istart_south(ith)+nph-1-resS_start,2) = ring(0:nph-1)
          call spinring_synthesis(H_res,nlmax, b_south_U, nph, ring, kphi0,mmax_ring)
          high_resS(H_res%istart_south(ith)-resS_start:H_res%istart_south(ith)+nph-1-resS_start,3) = ring(0:nph-1)
       endif

    enddo    ! loop on cos(theta)


    !     --------------------
    !     free memory and exit
    !     --------------------
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l, ring)
    DEALLOCATE(b_north_Q,b_north_U)
    DEALLOCATE(b_south_Q,b_south_U)
    deallocate(TEB)
    deallocate(H_res%trig)

    deallocate(H%recfac,stat= status)
    nullify(H%recfac)

     if (H%MpiID==0) then
       grad_phi => grad_phi_map
     else
        allocate(grad_phi(0:12*H%nside**2-1))    
        grad_phi = 0
     end if
#ifdef MPIPIX
     call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) then
      print *,code //': Got grad_phi ',H%MpiID, GeteTime() - StartTime
     end if
#endif


#ifdef MPIPIX
     allocate(map2N(H%North_Start(H%MpiID):H%North_Start(H%MpiID)+H%North_Size(H%MpiID)-1,3), stat = status)
     if (status /= 0) call die_alloc(code,'map2N')
     allocate(map2S(H%South_Start(H%MpiID):H%South_Start(H%MpiID)+H%South_Size(H%MpiID)-1,3),stat = status)
     if (status /= 0) call die_alloc(code,'map2S')
#else
     map2N => map_TQU
     map2S => map_TQU
#endif


    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)      ! 0 <= cos theta < 1

       if (ith < nsmax) then  ! polar cap (north)
          cth0 = 1.0_dp  - DBLE(ith)**2 * dth1
          nph = 4*ith
          kphi0 = 1
          sth0 = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth0 = DBLE(2*nsmax-ith) * dth2
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth0 = DSQRT((1.0_dp-cth0)*(1.0_dp+cth0)) ! sin(theta)
       endif

       do ring_ix = H%istart_north(ith-1),H%istart_north(ith-1)+nph-1

        phi = (ring_ix-H%istart_north(ith-1))*2*pi/nph

        grad_len = abs(grad_phi(ring_ix))
        if (grad_len>0) then
            sinc_grad_len = sin(grad_len)/grad_len
            cth =  cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
            sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
            DeltaPhi = asin(aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )
            phi = phi  +  DeltaPhi 
        else
         DeltaPhi = 0
         cth=cth0
         sth=sth0
        endif

        if (kphi0==1) phi=phi + pi/nph

        call ang2pix_ring8(high_nside, cth, phi, ipix)
        if (ipix >= H_res%istart_south(2*high_nside-1)) then
          map2N(ring_ix,:) = high_resS(ipix-resS_start,:)  
        else
          map2N(ring_ix,:) = high_resN(ipix-resN_start,:)
        end if         

     !Put in factor for change of co-ordinate axes between deflected and original point
          if (grad_len > 1e-20_dp) then
             Re = real(grad_phi(ring_ix))
             Im = aimag(grad_phi(ring_ix))
             gamma = grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
             if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
             gamma = Im/gamma
             !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
             gammfac = cmplx(map2N(ring_ix,2),map2N(ring_ix,3))* &
               cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
                    2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))             
             map2N(ring_ix,2) = real(gammfac)
             map2N(ring_ix,3) = aimag(gammfac)
         end if

       end do !ring_ix
      
      if (ith < 2*nsmax) then
       do ring_ix = H%istart_south(ith),H%istart_south(ith)+nph-1
            phi = (ring_ix-H%istart_south(ith))*2*pi/nph
            grad_len = abs(grad_phi(ring_ix))
            if (grad_len>0) then
                sinc_grad_len = sin(grad_len)/grad_len
                cth =  -cos(grad_len) * cth0 - sinc_grad_len *sth0*real(grad_phi(ring_ix)) 
                sth = max(1e-10_dp,sqrt((1._dp-cth)*(1._dp+cth)))
                DeltaPhi = asin(aimag(grad_phi(ring_ix))*sinc_grad_len/ sth )
                phi = phi + DeltaPhi 
            else
             DeltaPhi = 0
             cth=-cth0
             sth=sth0
            endif
           if (kphi0==1) phi=phi + pi/nph
           call ang2pix_ring8(high_nside, cth, phi, ipix)
           if (ipix >= H_res%istart_south(2*high_nside-1)) then
              map2S(ring_ix,:) = high_resS(ipix-resS_start,:)  
           else
              map2S(ring_ix,:) = high_resN(ipix-resN_Start,:)
           end if         

          
          if (grad_len > 1e-20_dp) then
             Re = real(grad_phi(ring_ix))
             Im = aimag(grad_phi(ring_ix))
             gamma = grad_len*sin(grad_len)*cth0/sth0 + Re*cos(grad_len)
             if (abs(gamma) < 1e-20_dp) gamma = 1e-20_dp  
             gamma = Im/gamma
             !use identity for cos(2(tan^{-1} A - tan^{-1} B)) and similarly sin
             gammfac = cmplx(map2S(ring_ix,2),map2S(ring_ix,3))* &
               cmplx( 2*((Re + Im*gamma)/grad_len)**2/(1+gamma**2) -1 ,  &
                    2*(Re+gamma*Im)*(Im - gamma*Re)/grad_len**2/(1+gamma**2))             
             map2S(ring_ix,2) = real(gammfac)
             map2S(ring_ix,3) = aimag(gammfac)
         end if

       end do
      end if      

    end do !ith

    deallocate(high_resS,high_resN)
    call HealpixFree(H_res)

#ifdef MPIPIX
    if (H%MpiID/=0) then
        deallocate(grad_phi)    
    end if
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiID, GeteTime()-StartTime
    StartTime = Getetime()
    do i=1,3
    call MPI_GATHERV(map2N(H%North_Start(H%MpiID),i),H%North_Size(H%MpiID),SP_MPI, &
       map_TQU(:,i),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2S(H%South_Start(H%MpiID),i),H%South_Size(H%MpiID),SP_MPI, &
       map_TQU(:,i),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    end do
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiID, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2N,map2S) 
#endif

  end subroutine alm2LensedmapInterp



  subroutine alm2LensedQuadContrib(H, inlmax, alm, grad_phi_map, map_T)
 !Get lensed map using fourth order series expansion
 !*** Note this does not give a map with an accurate lensed C_l at l>~ 1200***
    use alm_tools
    use MPIstuff
    Type (HealpixInfo) :: H
    INTEGER(I4B), INTENT(IN) :: inlmax 
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:)  :: alm
    REAL(SPC), INTENT(OUT), DIMENSION(0:12*H%nside**2-1), target :: map_T
    COMPLEX(SPC), INTENT(IN), DIMENSION(0:12*H%nside**2-1), target :: grad_phi_map
    REAL(SPC), DIMENSION(:), pointer :: map2
    COMPLEX(SPC),  DIMENSION(:), pointer :: grad_phi
    COMPLEX(SPC), DIMENSION(:), allocatable :: alm2

    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0, nlmax

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: a_w, b_w, a_x, fac
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n, b_s, b_n4,b_s4     
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U
    COMPLEX(DPC) :: b_n_Q3, b_s_Q3, b_n_U3, b_s_U3
    COMPLEX(DPC) :: b_n_Q4, b_s_Q4, b_n_U4, b_s_U4
    COMPLEX(DPC) :: b_n_Q33, b_s_Q33, b_n_U33, b_s_U33
    COMPLEX(DPC) :: b_n_Q44, b_s_Q44, b_n_U44, b_s_U44
    COMPLEX(DPC) :: b_n_Q2, b_s_Q2, b_n_U2, b_s_U2

    CHARACTER(LEN=*), PARAMETER :: code = 'ALM2LENSEDQUADCONTRIB'
    COMPLEX(DPC) ::  b(0:H%lmax,2), b4(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q(0:H%lmax,2), b_U(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q2(0:H%lmax,2), b_U2(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q3(0:H%lmax,2), b_U3(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q33(0:H%lmax,2), b_U33(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q44(0:H%lmax,2), b_U44(0:H%lmax,2)
    COMPLEX(DPC) ::  b_Q4(0:H%lmax,2), b_U4(0:H%lmax,2)

    REAL(DP) :: LastX(4),LastW(4),W(4),X(4)
    
    INTEGER(I4B) :: status,par_lm, a_ix

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(:),   ALLOCATABLE :: ringR, ringI
    integer NS, mmax_ring, nalms, ix
    double precision Initime
    !=======================================================================

     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiID==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
    call MPI_BCAST(nlmax,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(alm2(nalms),stat = status )
     if (status /= 0) call die_alloc(code,'alm2')
     if (H%MpiID==0) then
       call Alm2PackAlm(alm,alm2,nlmax)
       grad_phi => grad_phi_map
     else
       allocate(grad_phi(0:12*H%nside**2-1)) 
     end if
#ifdef MPIPIX
     call MPI_BCAST(alm2,SIze(alm2),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     call MPI_BCAST(grad_phi,SIze(grad_phi),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code//' Got alm ',H%MpiID, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1), stat = status)
     if (status /= 0) call die_alloc(code,'map2')    
#else
     map2 => map_T
#endif

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(ringR(0:4*nsmax-1),ringI(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    call HealpixInitRecfac(H,nlmax)
    call GetLamfact(lam_fact, nlmax)
!Don't put in 2 factor for spin 2 here, but below   
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )
    !     --------------------------------------------

    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)   ! 0 <= cos theta < 1

       !        cos(theta) in the pixelisation scheme
       if (ith < nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2 !cos theta
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       
       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       lam_mm = sq4pi_inv ! lamda_00
       scalem=1
       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2
          !           ---------- l = m ----------
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)

          lam_lm = corfac*lam_mm/OVFLOW     !  actual lambda_mm      
          
          a_ix = a_ix + 1

          par_lm = 1  

!Del^2 T term

          b_n = -(fm2+fm)*lam_lm * alm2(a_ix)
          b_s = b_n

!4th order
          fac = (2-3*(fm2+fm))
          b_n4 = factor * b_n
          b_s4 = b_n4

          par_lm = -1
!Grad T term
          if (m >=1) then
!normal_l cancels with gradient, sign from gradient 
              X(1) =  lam_lm * fm / sth 
              W(1) =  -X(1) * cth

              b_n_Q =  -W(1) * alm2(a_ix)
              b_s_Q =  par_lm * b_n_Q

              b_n_U = cmplx(0._dp,X(1)) * alm2(a_ix)
              b_s_U = -par_lm * b_n_U

              b_n_Q3 =  -fac * W(1) * alm2(a_ix)
              b_s_Q3 =  par_lm * b_n_Q

              b_n_U3 = cmplx(0._dp,X(1)*fac) * alm2(a_ix)
              b_s_U3 = -par_lm * b_n_U

          else
             W(1) = 0
             X(1) = 0

             b_n_Q=0; b_s_Q=0; b_n_U=0; b_s_U=0
             b_n_Q3=0; b_s_Q3=0; b_n_U3=0; b_s_U3=0
          end if

!eth eth T term
          par_lm = 1
          if (m >=2) then
              W(2) = -( lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
              X(2) =  ( lam_lm * (fm - fm2) ) *   2.0_dp * c_on_s2
              b_n_Q2 =  W(2)* alm2(a_ix) 
              b_s_Q2 =  par_lm* b_n_Q2

              b_n_U2 =  cmplx(0.0_dp, -X(2)) * alm2(a_ix)
              b_s_U2 =  -par_lm*b_n_U2  
 
 !4th order spin 2 term
              b_n_Q4 =  (8 - 4*(fm+fm2)) * W(2)* alm2(a_ix) 
              b_s_Q4 =  par_lm* b_n_Q4

              b_n_U4 =  cmplx(0.0_dp, -(8 - 4*(fm+fm2))*X(2)) * alm2(a_ix)
              b_s_U4 =  -par_lm*b_n_U4  

          else
             W(2)=0;X(2)=0
             b_n_Q2=0;b_s_Q2=0;b_n_U2=0;b_s_U2=0
             b_n_Q4=0;b_s_Q4=0;b_n_U4=0;b_s_U4=0
          end if

!Irreducible spin 3 term
          par_lm = -1
          if (m >=3) then
              X(3) =  lam_lm / sth * fm*(fm-1)*(fm-2) 
              W(3) = X(3) * cth * ( 1 - 4*one_on_s2)
              X(3) = X(3) * (4*one_on_s2 - 3)             

              b_n_Q33 =  W(3) * alm2(a_ix) 
              b_s_Q33 =  par_lm* b_n_Q2

              b_n_U33 =  cmplx(0.0_dp, -X(3)) * alm2(a_ix)
              b_s_U33 =  -par_lm*b_n_U2  

          else
             W(3)=0;X(3)=0
             b_n_Q33=0;b_s_Q33=0;b_n_U33=0;b_s_U33=0
          end if

          par_lm = 1
  
           if (m>=4 ) then
   
             W(4) = ((l-3)*(X(3) - cth*W(3))) /(sth)
             X(4) = ((l-3)*(W(3) - cth*X(3))) /(sth)
             factor_1 =  W(4) * alm2(a_ix)
             b_n_Q44 =  factor_1
             b_s_Q44 = par_lm * factor_1 

             factor_2 =  cmplx(0.0_dp, -X(4)) * alm2(a_ix) 
             b_n_U44 =  factor_2
             b_s_U44 =  -par_lm * factor_2
           else 
             b_n_Q44=0;b_s_Q44=0;b_n_U44=0;b_s_U44=0
           end if


!Keep par_lm correct for spin zero
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m
             LastX = X
             LastW = W 
             lam_lm = lam_2 * corfac*lam_mm/OVFLOW ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl

             a_ix = a_ix + 1

!Del^2 T
             factor = -(fl2+fl)*lam_lm * alm2(a_ix)
             b_n = b_n +          factor
             b_s = b_s + par_lm * factor

!4th order
             fac = (2-3*(fl2+fl))
             factor = factor * fac 
             b_n4 = b_n4 + factor
             b_s4 = b_s4 + par_lm*factor

!grad T
             par_lm = -par_lm
             a_w = 1 / sth
             X(1) = a_w * fm  * lam_lm
             W(1) = a_w * (lam_fact(a_ix)*lam_lm1m - fl*cth*lam_lm) 

             factor_1 =  -W(1) * alm2(a_ix)
             b_n_Q = b_n_Q +          factor_1 
             b_s_Q = b_s_Q + par_lm * factor_1 ! X has a diff. parity

             factor_2 =   cmplx(0._dp,-X(1)) * alm2(a_ix) 
             b_n_U = b_n_U - factor_2
             b_s_U = b_s_U + par_lm * factor_2

!grad [(3 grad^2 + 2)T]
             factor_1 = - fac * W(1) * alm2(a_ix)
             b_n_Q3 = b_n_Q3 +          factor_1 
             b_s_Q3 = b_s_Q3 + par_lm * factor_1 ! X has a diff. parity

             factor_2 =   cmplx(0._dp,-fac * X(1)) * alm2(a_ix) 
             b_n_U3 = b_n_U3 - factor_2
             b_s_U3 = b_s_U3 + par_lm * factor_2


!eth eth T
             par_lm = -par_lm
             if (l>=2 .and. corfac /= 0) then

                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
    !put in 2 factor here
                 b_w =  2* c_on_s2 * lam_fact(a_ix)
                 a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 W(2) =  ( a_w * lam_lm + b_w * lam_lm1m ) 
    !and here
                 X(2) =  fm_on_s2 * ( 2* lam_fact(a_ix) * lam_lm1m - a_x)

                 factor_1 =  W(2)* alm2(a_ix)
                 b_n_Q2 = b_n_Q2 +  factor_1
                 b_s_Q2 = b_s_Q2 + par_lm * factor_1 

                 factor_2 =  cmplx(0.0_dp, X(2)) * alm2(a_ix) 
                 b_n_U2 = b_n_U2  - factor_2
                 b_s_U2 = b_s_U2 + par_lm * factor_2
    
    !4th order spin 2 term

                 factor_1 = (8-4*(fl2+fl))* W(2)* alm2(a_ix)
                 b_n_Q4 = b_n_Q4 +  factor_1
                 b_s_Q4 = b_s_Q4 + par_lm * factor_1 

                 factor_2 =  cmplx(0.0_dp, (8-4*(fl2+fl))* X(2)) * alm2(a_ix) 
                 b_n_U4 = b_n_U4  - factor_2
                 b_s_U4 = b_s_U4 + par_lm * factor_2

              else
               X(2)=0;W(2)=0;
              end if
!Irreducible spin 3 term
             par_lm = -par_lm
             if (l>=3 .and. corfac /= 0) then
   
           W(3) = ((l+2)*lam_fact(a_ix)*LastW(2) + (l-2)*(m*X(2) - cth*l*W(2))) /(l*sth)
           X(3) = ((l+2)*lam_fact(a_ix)*LastX(2) + (l-2)*(m*W(2) - cth*l*X(2))) /(l*sth)
!              a_w = 1._dp /sth
!                 b_w = (l-1)*(l-2)
!                 lambda_x = a_w*fm*( (one_on_s2*(4*fm2-(12*l-8)) - 3 * b_w)*lam_lm + &
!                         12*lam_fact(a_ix) * c_on_s2  * lam_lm1m)  
 !                lambda_w = a_w*( (l*b_w - one_on_s2*(8*l+fm2*(4*l-12))) * cth * lam_lm  &
 !                       - lam_fact(a_ix)*( fl2 + fl +6 - (8+4*fm2)*one_on_s2) * lam_lm1m)
                 factor_1 =  W(3) * alm2(a_ix)
                 b_n_Q33 = b_n_Q33 +  factor_1
                 b_s_Q33 = b_s_Q33 + par_lm * factor_1 

                 factor_2 =  cmplx(0.0_dp, X(3)) * alm2(a_ix) 
                 b_n_U33 = b_n_U33  - factor_2
                 b_s_U33 = b_s_U33 + par_lm * factor_2
             else
              W(3)=0;X(3)=0
             end if
             par_lm = -par_lm

!Irreducible spin 4 term
           if (l>=4 .and. corfac /= 0) then
   
             W(4) = ((l+3)*lam_fact(a_ix)*LastW(3) + (l-3)*(m*X(3) - cth*l*W(3))) /(l*sth)
             X(4) = ((l+3)*lam_fact(a_ix)*LastX(3) + (l-3)*(m*W(3) - cth*l*X(3))) /(l*sth)
             factor_1 =  W(4) * alm2(a_ix)
             b_n_Q44 = b_n_Q44 +  factor_1
             b_s_Q44 = b_s_Q44 + par_lm * factor_1 

             factor_2 =  cmplx(0.0_dp, X(4)) * alm2(a_ix) 
             b_n_U44 = b_n_U44  - factor_2
             b_s_U44 = b_s_U44 + par_lm * factor_2
     
           end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif

          enddo

          b(m,1) = b_n; b(m,2) = b_s
          b4(m,1) = b_n4; b4(m,2) = b_s4

          b_Q(m,1) = b_n_Q; b_Q(m,2) = b_s_Q; b_U(m,1) = b_n_U; b_U(m,2) = b_s_U
          b_Q2(m,1) = b_n_Q2; b_Q2(m,2) = b_s_Q2; b_U2(m,1) = b_n_U2; b_U2(m,2) = b_s_U2
          b_Q3(m,1) = b_n_Q3; b_Q3(m,2) = b_s_Q3; b_U3(m,1) = b_n_U3; b_U3(m,2) = b_s_U3
          b_Q4(m,1) = b_n_Q4; b_Q4(m,2) = b_s_Q4; b_U4(m,1) = b_n_U4; b_U4(m,2) = b_s_U4
          b_Q33(m,1) = b_n_Q33; b_Q33(m,2) = b_s_Q33; b_U33(m,1) = b_n_U33; b_U33(m,2) = b_s_U33
          b_Q44(m,1) = b_n_Q44; b_Q44(m,2) = b_s_Q44; b_U44(m,1) = b_n_U44; b_U44(m,2) = b_s_U44

       enddo

       ix = H%istart_north(ith-1)
       do NS = 1,2  
!Fourth order
   !spin zero
           call spinring_synthesis(H,nlmax, b4(:,NS), nph, ringR, kphi0,mmax_ring)
           map2(ix:ix+nph-1) =  &
           (real(grad_phi(ix:ix+nph-1))**2+aimag(grad_phi(ix:ix+nph-1))**2)**2*RingR(0:nph-1)
            !factor 1/(4!) goes in below
   !irreducible spin 4
           call spinring_synthesis(H,nlmax, b_Q44(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U44(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = map2(ix:ix+nph-1) + real( &
             grad_phi(ix:ix+nph-1)**4*cmplx(RingR(0:nph-1),-RingI(0:nph-1)) )
     
   !spin 2 cross
           call spinring_synthesis(H,nlmax, b_Q4(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U4(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = (map2(ix:ix+nph-1) + &
             (real(grad_phi(ix:ix+nph-1))**2+aimag(grad_phi(ix:ix+nph-1))**2) * &
             real( grad_phi(ix:ix+nph-1)**2*cmplx(RingR(0:nph-1),-RingI(0:nph-1)) ))/8


!Cubic reducible term
           call spinring_synthesis(H,nlmax, b_Q3(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U3(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = map2(ix:ix+nph-1) + ((real(grad_phi(ix:ix+nph-1))*RingR(0:nph-1) + &
                              aimag(grad_phi(ix:ix+nph-1))*RingI(0:nph-1)) * &
                       (real(grad_phi(ix:ix+nph-1))**2+aimag(grad_phi(ix:ix+nph-1))**2))
!Cubic irreducible term
           call spinring_synthesis(H,nlmax, b_Q33(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U33(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = (map2(ix:ix+nph-1) + real( &
             -grad_phi(ix:ix+nph-1)**3*cmplx(RingR(0:nph-1),-RingI(0:nph-1)) )) / (3*2)
                   !factor of 1/4 goes in with quadratic addition


    !Quadratic
      ! Re(eth eth T ethb phi ethb phi)
           call spinring_synthesis(H,nlmax, b_Q2(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U2(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = map2(ix:ix+nph-1) + &
                   (real(grad_phi(ix:ix+nph-1))+ aimag(grad_phi(ix:ix+nph-1)))* &
                   (real(grad_phi(ix:ix+nph-1))- aimag(grad_phi(ix:ix+nph-1))) * RingR(0:nph-1) &
              + 2*RingI(0:nph-1) * real(grad_phi(ix:ix+nph-1))*aimag(grad_phi(ix:ix+nph-1)) 
      !grad^2 T |Grad phi|^2
           call spinring_synthesis(H,nlmax, b(:,NS), nph, ringR, kphi0,mmax_ring)
           map2(ix:ix+nph-1) =0.25_dp * ( map2(ix:ix+nph-1) + &
               (real(grad_phi(ix:ix+nph-1))**2+aimag(grad_phi(ix:ix+nph-1))**2)*RingR(0:nph-1))

    !Linear
     !grad T dot grad phi
           call spinring_synthesis(H,nlmax, b_Q(:,NS), nph, ringR, kphi0,mmax_ring)
           call spinring_synthesis(H,nlmax, b_U(:,NS), nph, ringI, kphi0,mmax_ring)
           map2(ix:ix+nph-1) = map2(ix:ix+nph-1) + real(grad_phi(ix:ix+nph-1))*RingR(0:nph-1) + &
                              aimag(grad_phi(ix:ix+nph-1))*RingI(0:nph-1)

           if (ith  >=  2*nsmax) exit
           ix = H%istart_south(ith)
      end do

    enddo    ! loop on cos(theta)

    !     --------------------
    !     free memory and exit
    !     --------------------
    call healpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(ringR,ringI)
    deallocate(alm2)
#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiID
    StartTime = Getetime()
    call MPI_GATHERV(map2(H%North_Start(H%MpiID)),H%North_Size(H%MpiID),SP_MPI, &
       map_T,H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2(H%South_Start(H%MpiID)),H%South_Size(H%MpiID),SP_MPI, &
       map_T,H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    if (H%MpiID/=0) deallocate(grad_phi)
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiID, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2) 
#endif

  end subroutine alm2LensedQuadContrib


  subroutine map2polalm(H, inlmax, map_TQU,alm_TEB, cos_theta_cut)
    use MPIStuff
    Type (HealpixInfo) :: H

    INTEGER(I4B)  :: nsmax
    INTEGER(I4B), INTENT(IN) :: inlmax
    COMPLEX(SPC), INTENT(OUT),  DIMENSION(:,:,:) :: alm_TEB
    REAL(SP), INTENT(IN), DIMENSION(0:12*H%nside**2-1,3), target :: map_TQU
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: TEB, TEB2
    REAL(SP), DIMENSION(:,:), pointer :: map2

    REAL(DP),     INTENT(IN) :: cos_theta_cut

    Integer :: nlmax !base integer type for MPI compat
    INTEGER(I4B) :: l, m, ith, scalem, scalel        ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: omega_pix
    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x

    CHARACTER(LEN=*), PARAMETER :: code = 'MAP2POLALM'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_nQ, phas_nU
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_sQ, phas_sU
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_n
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: phas_s

    INTEGER(I4B) mmax_ring, status, par_lm, a_ix
    integer i, nalms

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: ring
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    LOGICAL   :: keep_it
    double precision Initime
    !=======================================================================

    nsmax = H%nside
    nlmax = inlmax
     
#ifdef MPIPIX
     if (cos_theta_cut/=-1) stop 'cos_theta_cut /= -1'
     if (H%MpiID==0) then 
      if(DebugMsgs>0) print *,code //': Sending to farm '
      call SendMessages(H,code)
      map2 => map_TQU
    else
       allocate(map2(0:12*H%nside**2-1,3),stat = status) 
       if (status /= 0) call die_alloc(code,'map2')   
    end if

    StartTime = getetime()    
    iniTime = StartTime
    call MPI_BCAST(nlmax,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
    do i=1,3
    call MPI_SCATTERV(map_TQU(:,i),H%North_Size, H%North_Start, &
       SP_MPI, map2(H%North_Start(H%MpiID),i),H%North_Size(H%MpiID),SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_SCATTERV(map_TQU(:,i),H%South_Size, H%South_Start, &
       SP_MPI, map2(H%South_Start(H%MpiID),i),H%South_Size(H%MpiID),SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    end do
    if(DebugMsgs>1) print *,code //' Scattered ',H%MpiID, GeteTime() - StartTime
#else
    map2 => map_TQU
#endif

    nalms = ((nlmax+1)*(nlmax+2))/2   

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')


    ALLOCATE(phas_n(0:nlmax), phas_nQ(0:nlmax),&
         &   phas_nU(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_n')

    ALLOCATE(phas_s(0:nlmax), phas_sQ(0:nlmax),&
         &   phas_sU(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'phas_s')

    ALLOCATE(ring(0:4*nsmax-1),stat = status) 
    if (status /= 0) call die_alloc(code,'ring')

    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamfact(lam_fact, nlmax)
   lam_fact = lam_fact*2

   allocate(TEB(3,nalms),stat = status)
   if (status /= 0) call die_alloc(code,'TEB')
   TEB = 0    
       
   omega_pix = pi / (3 * nsmax * real(nsmax,dp))

    normal_l = 0.0_dp
    do l = 2, nlmax
       fl = DBLE(l)
        normal_l(l) = EB_sign * SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo

    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !-----------------------------------------------------------------------
    !           computes the integral in phi : phas_m(theta)
    !           for each parallele from north to south pole
    !-----------------------------------------------------------------------
    
    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)

       phas_nQ=0; phas_sQ=0;phas_nU=0;phas_sU=0; phas_n=0; phas_s=0

       if (ith  <=  nsmax-1) then      ! north polar cap
          nph = 4*ith
          kphi0 = 1 
          cth = 1.0_dp  - DBLE(ith)**2 * dth1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                            ! tropical band + equat.
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          cth = DBLE(2*nsmax-ith) * dth2
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2

       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       keep_it = (ABS(cth) > cos_theta_cut) ! part of the sky out of the symmetric cut
       if (keep_it) then
          ring(0:nph-1) = map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,1) * H%w8ring_TQU(ith,1)
          call spinring_analysis(H,nlmax, ring, nph, phas_n, kphi0, mmax_ring)
          ring(0:nph-1) = map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,2) * H%w8ring_TQU(ith,2)
          call spinring_analysis(H,nlmax, ring, nph, phas_nQ, kphi0, mmax_ring)
          ring(0:nph-1) = map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,3) * H%w8ring_TQU(ith,3) 
          call spinring_analysis(H,nlmax, ring, nph, phas_nU, kphi0, mmax_ring)
       endif

       if (ith  <  2*nsmax .and. keep_it) then
          ring(0:nph-1) = map2(H%istart_south(ith):H%istart_south(ith)+nph-1,1) * H%w8ring_TQU(ith,1)
          call spinring_analysis(H,nlmax, ring, nph, phas_s, kphi0, mmax_ring)
          ring(0:nph-1) = map2(H%istart_south(ith):H%istart_south(ith)+nph-1,2) * H%w8ring_TQU(ith,2)
          call spinring_analysis(H,nlmax, ring, nph, phas_sQ, kphi0, mmax_ring)
          ring(0:nph-1) = map2(H%istart_south(ith):H%istart_south(ith)+nph-1,3) * H%w8ring_TQU(ith,3)
          call spinring_analysis(H,nlmax, ring, nph, phas_sU, kphi0, mmax_ring)
       endif

       !-----------------------------------------------------------------------
       !              computes the a_lm by integrating over theta
       !                  lambda_lm(theta) * phas_m(theta)
       !                         for each m and l
       !-----------------------------------------------------------------------

       if (keep_it) then ! avoid un-necessary calculations (EH, 09-2001)
          lam_mm = sq4pi_inv * omega_pix 
          scalem=1
          a_ix = 0
          do m = 0, mmax_ring
             fm  = DBLE(m)
             f2m = 2.0_dp * fm
             fm2 = fm * fm
             fm_on_s2 = fm * one_on_s2

             !           ---------- l = m ----------
             par_lm = 1   ! = (-1)^(l+m+s)
             if (m  >=  1) then ! lambda_0_0 for m>0
                lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
             endif

             if (abs(lam_mm) < UNFLOW) then
                lam_mm=lam_mm*OVFLOW
                scalem=scalem-1
             endif

             a_ix = a_ix+1

             corfac = ScaleFactor(scalem)
             lam_lm = corfac*lam_mm/OVFLOW

             TEB(1,a_ix) = TEB(1,a_ix) + lam_lm * (phas_n(m) + phas_s(m))
             if (m >=2) then
            
              lambda_w = - ( normal_l(m) * lam_lm * (fm - fm2) ) * ( 2.0_dp * one_on_s2 - 1.0_dp )  
              lambda_x = ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)
              
                 TEB(2,a_ix) = TEB(2,a_ix) &
                  &                 + lambda_w * (phas_nQ(m) + par_lm*phas_sQ(m)) &
                  &                 + zi_lam_x * (phas_nU(m) - par_lm* phas_sU(m))

                 TEB(3,a_ix) = TEB(3,a_ix) &
                  &                 + lambda_w * (phas_nU(m) + par_lm*phas_sU(m)) &
                  &                 - zi_lam_x * (phas_nQ(m) - par_lm*phas_sQ(m))
            
             end if

             !           ---------- l > m ----------
             lam_0 = 0.0_dp
             lam_1 = 1.0_dp
             scalel=0
             a_rec = H%recfac(a_ix)
             lam_2 = cth * lam_1 * a_rec
             do l = m+1, nlmax
                par_lm = - par_lm  ! = (-1)^(l+m)
                lam_lm1m=lam_lm ! actual lambda_l-1,m (useful for polarisation)
                lam_lm   = lam_2*corfac*lam_mm/OVFLOW ! actual lambda_lm (OVFLOW factors removed)
                fl  = DBLE(l)
                fl2 = fl * fl

                a_ix = a_ix + 1

                TEB(1,a_ix) = TEB(1,a_ix) + lam_lm * (phas_n(m) + par_lm*phas_s(m))

             if (l>=2 .and. corfac /= 0) then
                 !Corfac=0 guarantees lam(l-1) is also v close to zero
                  
                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  2*(l-1)* cth * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
  
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 TEB(2,a_ix) = TEB(2,a_ix) &
                     &          + lambda_w * (phas_nQ(m) + par_lm*phas_sQ(m)) &
                     &          + zi_lam_x * (phas_nU(m) - par_lm*phas_sU(m))
                 TEB(3,a_ix) = TEB(3,a_ix) &
                     &         +  lambda_w * (phas_nU(m) + par_lm*phas_sU(m)) &
                     &         - zi_lam_x  * (phas_nQ(m) - par_lm*phas_sQ(m))

              end if ! l allowed by spin or zero

                lam_0 = lam_1 / a_rec
                lam_1 = lam_2
                a_rec = H%recfac(a_ix)
                lam_2 = (cth * lam_1 - lam_0) * a_rec

                if (abs(lam_2)  >  OVFLOW) then
                   lam_0=lam_0/OVFLOW
                   lam_1=lam_1/OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec
                   scalel=scalel+1
                   corfac = ScaleFactor(scalem+scalel)
                elseif (abs(lam_2)  <  UNFLOW) then
                   lam_0=lam_0*OVFLOW
                   lam_1=lam_1*OVFLOW
                   lam_2 = (cth * lam_1 - lam_0) * a_rec 
                   scalel=scalel-1
                   corfac = ScaleFactor(scalem+scalel)
                endif

             enddo ! loop on l
          enddo ! loop on m
       endif ! test on cut sky
    enddo ! loop on theta

    !     --------------------
    !     free memory and exit
    !     --------------------
    call HealpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l)
    DEALLOCATE(phas_nQ,phas_nU)
    DEALLOCATE(phas_sQ,phas_sU)
    DEALLOCATE(phas_n)
    DEALLOCATE(phas_s)
    DEALLOCATE(ring)
#ifdef MPIPIX
    if (H%MpiID>0) deallocate(map2)
    allocate(TEB2(3,nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'TEB2')
    StartTime = Getetime()
    call MPI_REDUCE(TEB,TEB2,size(TEB),CSP_MPI,MPI_SUM,0,MPI_COMM_WORLD,l) 
    if (DebugMsgs>1) print *,code//' done reduce ', H%MpiID, GeteTime() -StartTime
    deallocate(TEB)
    if (H%MpiID == 0) call PackTEB2TEB(TEB2,alm_TEB, nlmax)
    deallocate(TEB2)
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code //' Time: ',GeteTime() - IniTime
#else
    call PackTEB2TEB(TEB,alm_TEB, nlmax)
    deallocate(TEB)
#endif

  END subroutine map2polalm


  subroutine polalm2map(H,inlmax, alm_TEB, map_TQU)
    use MPIstuff
    Type (HealpixInfo) :: H

    INTEGER(I4B), INTENT(IN) :: inlmax
    integer nsmax
    COMPLEX(SPC), INTENT(IN),  DIMENSION(:,:,:) :: alm_TEB
    REAL(SP), INTENT(OUT), DIMENSION(0:12*H%nside**2-1,3), target :: map_TQU
    COMPLEX(SPC), DIMENSION(:,:), allocatable :: TEB
    REAL(SP), DIMENSION(:,:), pointer :: map2
    INTEGER(I4B) :: l, m, ith, scalem, scalel          ! alm related
    INTEGER(I4B) :: nph, kphi0 ! map related

    REAL(DP) :: cth, sth, dth1, dth2, dst1
    REAL(DP) :: a_rec, lam_mm, lam_lm, lam_lm1m, lam_0, lam_1, lam_2
    REAL(DP) :: fm, f2m, fm2, fl, fl2, corfac
    REAL(DP) :: c_on_s2, fm_on_s2, one_on_s2
    REAL(DP) :: lambda_w, lambda_x, a_w, b_w, a_x
    COMPLEX(DPC) :: zi_lam_x
    COMPLEX(DPC) :: factor, factor_1, factor_2
    COMPLEX(DPC) :: b_n_Q, b_s_Q, b_n_U, b_s_U
    COMPLEX(DPC) :: b_n, b_s

    CHARACTER(LEN=*), PARAMETER :: code = 'POLALM2MAP'
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_north_Q, b_north_U
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE ::  b_south_Q, b_south_U
    COMPLEX(DPC), DIMENSION(0:H%lmax) :: b_north,b_south
    INTEGER(I4B) :: mmax_ring,status,par_lm, nlmax

    REAL(DP), DIMENSION(:), ALLOCATABLE :: lam_fact
    REAL(SP), DIMENSION(0:4*H%nside-1) :: ring
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: normal_l
    integer i, a_ix, nalms
    double precision Initime
    !=======================================================================

    !     --- allocates space for arrays ---


     nsmax = H%nside
     nlmax = inlmax

#ifdef MPIPIX
    StartTime = Getetime()
    iniTime = StartTime
    if (H%MpiID==0) then 
     print *,code //': Sending to farm ' 
     call SendMessages(H,code)
    end if
     call MPI_BCAST(nlmax,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
#endif
     nalms = ((nlmax+1)*(nlmax+2))/2   
     allocate(TEB(3,nalms))
     if (H%MpiID==0) call TEB2PackTEB(alm_TEB,TEB,nlmax)
    
#ifdef MPIPIX
     call MPI_BCAST(TEB,SIze(TEB),CSP_MPI, 0, MPI_COMM_WORLD, ierr) 
     if(DebugMsgs>1) print *,code //': Got alm ',H%MpiID, GeteTime() - StartTime
     allocate(map2(0:12*nsmax**2-1,3), stat = status)    
     if (status /= 0) call die_alloc(code,'map2')
#else
     map2 => map_TQU 
#endif

    ALLOCATE(lam_fact(nalms),stat = status)    
    if (status /= 0) call die_alloc(code,'lam_fact')

    ALLOCATE(normal_l(0:nlmax),stat = status)    
    if (status /= 0) call die_alloc(code,'normal_l')

    ALLOCATE(b_north_Q(0:nlmax),&
         &   b_north_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_north')

    ALLOCATE(b_south_Q(0:nlmax),&
         &   b_south_U(0:nlmax),stat = status) 
    if (status /= 0) call die_alloc(code,'b_south')

    !     ------------ initiate arrays ----------------

   call HealpixInitRecfac(H,nlmax)
   call GetLamFact(lam_fact, nlmax)
   lam_fact = lam_fact * 2 !HealPix polarization def
   

    normal_l = 0.0_dp
    do l = 2, nlmax
       fl = DBLE(l)
       normal_l(l) = EB_sign*SQRT( 1/ ((fl+2.0_dp)*(fl+1.0_dp)*fl*(fl-1.0_dp)) ) 
    enddo

 
    dth1 = 1.0_dp / (3.0_dp*DBLE(nsmax)**2)
    dth2 = 2.0_dp / (3.0_dp*DBLE(nsmax))
    dst1 = 1.0_dp / (SQRT(6.0_dp) * DBLE(nsmax) )

    !     --------------------------------------------

    do ith = H%ith_start(H%MpiID), H%ith_end(H%MpiID)      ! 0 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
       if (ith < nsmax) then  ! polar cap (north)
          cth = 1.0_dp  - DBLE(ith)**2 * dth1  !cos theta
          nph = 4*ith
          kphi0 = 1
          sth = SIN( 2.0_dp * ASIN( ith * dst1 ) ) ! sin(theta)
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nsmax-ith) * dth2 !cos theta
          nph = 4*nsmax
          kphi0 = MOD(ith+1-nsmax,2)
          sth = DSQRT((1.0_dp-cth)*(1.0_dp+cth)) ! sin(theta)
       endif
       one_on_s2 = 1.0_dp / sth**2 ! 1/sin^2
       c_on_s2 = cth * one_on_s2
       !        -----------------------------------------------------
       !        for each theta, and each m, computes
       !        b(m,theta) = sum_over_l>m (lambda_l_m(theta) * a_l_m) 
       !        ------------------------------------------------------
       !        lambda_mm tends to go down when m increases (risk of underflow)
       !        lambda_lm tends to go up   when l increases (risk of overflow)
       lam_mm = sq4pi_inv ! lamda_00
       scalem=1

       if (mmax_approx) then
        mmax_ring = min(nlmax,max(40,nint(1.25*nlmax*sth)))
       else
        mmax_ring = nlmax
       end if

       a_ix = 0
       do m = 0, mmax_ring
          fm  = DBLE(m)
          f2m = 2.0_dp * fm
          fm2 = fm * fm
          fm_on_s2 = fm * one_on_s2

          !           ---------- l = m ----------
          par_lm = 1  ! = (-1)^(l+m+s)
          if (m  >=  1) then ! lambda_0_0 for m>0
             lam_mm = -lam_mm*sth*dsqrt((f2m+1.0_dp)/f2m)
          endif

          if (abs(lam_mm) < UNFLOW) then
             lam_mm=lam_mm*OVFLOW
             scalem=scalem-1
          endif
          corfac = ScaleFactor(scalem)
  
          ! alm_T * Ylm : Temperature
          lam_lm = corfac*lam_mm/OVFLOW     !  actual lambda_mm      

          a_ix = a_ix + 1

          b_n = lam_lm * TEB(1,a_ix)
          b_s = b_n

          !l=m special case
          if (m >=2) then
              lambda_w = - 2.0_dp *(normal_l(m) * lam_lm * (fm - fm2) ) * ( one_on_s2 - 0.5_dp )
              lambda_x =  ( normal_l(m) * lam_lm * (fm - fm2) ) *   2.0_dp *   c_on_s2
              
              zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

              b_n_Q =  lambda_w * TEB(2,a_ix) + zi_lam_x * TEB(3,a_ix)
              b_s_Q =  par_lm*(lambda_w * TEB(2,a_ix) - zi_lam_x * TEB(3,a_ix))

              b_n_U = lambda_w * TEB(3,a_ix) - zi_lam_x * TEB(2,a_ix)
              b_s_U = par_lm*(lambda_w * TEB(3,a_ix) + zi_lam_x * TEB(2,a_ix))

          else
             b_n_Q=0
             b_s_Q=0
             b_n_U=0
             b_s_U=0
          end if
          !           ---------- l > m ----------
          lam_0 = 0.0_dp
          lam_1 = 1.0_dp
          scalel=0
          a_rec = H%recfac(a_ix)
          lam_2 = cth * lam_1 * a_rec
          do l = m+1, nlmax
             par_lm = - par_lm  ! = (-1)^(l+m+s)
             lam_lm1m=lam_lm  ! actual lambda_l-1,m 
             lam_lm = lam_2 * corfac*lam_mm/OVFLOW ! actual lambda_lm, OVFLOW factors removed
             fl  = DBLE(l)
             fl2 = fl * fl
             a_ix = a_ix + 1

             factor = lam_lm * TEB(1,a_ix)
             b_n = b_n +          factor
             b_s = b_s + par_lm * factor

             if (l>=2 .and. corfac /= 0) then

                 a_w =  2* (fm2 - fl) * one_on_s2 - (fl2 - fl)
                 b_w =  c_on_s2 * lam_fact(a_ix)
                 a_x =  2.0_dp * cth * (fl-1.0_dp) * lam_lm
                 lambda_w =  normal_l(l) * ( a_w * lam_lm + b_w * lam_lm1m ) 
                 lambda_x =  normal_l(l) * fm_on_s2 * ( lam_fact(a_ix) * lam_lm1m - a_x)
                 zi_lam_x = CMPLX(0.0_dp, lambda_x, KIND=DP)

                 ! alm_G * Ylm_W - alm_C * Ylm_X : Polarisation Q
                 factor_1 =  lambda_w * TEB(2,a_ix)
                 factor_2 =  zi_lam_x * TEB(3,a_ix) ! X is imaginary
                 b_n_Q = b_n_Q +           factor_1 + factor_2
                 b_s_Q = b_s_Q + par_lm * (factor_1 - factor_2)! X has a diff. parity

                 !- alm_G * Ylm_X - alm_C * Ylm_W : Polarisation U
                 factor_1 =   lambda_w * TEB(3,a_ix) 
                 factor_2 =   zi_lam_x * TEB(2,a_ix) ! X is imaginary
                 b_n_U = b_n_U +           factor_1 - factor_2
                 b_s_U = b_s_U + par_lm * (factor_1 + factor_2)! X has a diff. parity
             end if

             lam_0 = lam_1 / a_rec
             lam_1 = lam_2
             a_rec = H%recfac(a_ix)
             lam_2 = (cth * lam_1 - lam_0) * a_rec

             if (abs(lam_2)  >  OVFLOW) then
                lam_0=lam_0/OVFLOW
                lam_1=lam_1/OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec
                scalel=scalel+1
                corfac = ScaleFactor(scalem+scalel)
             elseif (abs(lam_2)  <  UNFLOW) then
                lam_0=lam_0*OVFLOW
                lam_1=lam_1*OVFLOW
                lam_2 = (cth * lam_1 - lam_0) * a_rec 
                scalel=scalel-1
                corfac = ScaleFactor(scalem+scalel)
             endif

          enddo

          b_north_Q(m) = b_n_Q 
          b_south_Q(m) = b_s_Q 
          b_north_U(m) = b_n_U
          b_south_U(m) = b_s_U
          b_north(m) = b_n
          b_south(m) = b_s

       enddo

       call spinring_synthesis(H,nlmax,b_north,nph,ring,kphi0,mmax_ring)   
       map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,1) = ring(0:nph-1)
       call spinring_synthesis(H,nlmax, b_north_Q, nph, ring, kphi0,mmax_ring)
       map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,2) = ring(0:nph-1)
       call spinring_synthesis(H,nlmax, b_north_U, nph, ring, kphi0,mmax_ring)
       map2(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,3) = ring(0:nph-1)
  
       if (ith  <  2*nsmax) then
          call spinring_synthesis(H,nlmax, b_south, nph, ring, kphi0,mmax_ring)
          map2(H%istart_south(ith):H%istart_south(ith)+nph-1,1) = ring(0:nph-1)
          call spinring_synthesis(H,nlmax, b_south_Q, nph, ring, kphi0,mmax_ring)
          map2(H%istart_south(ith):H%istart_south(ith)+nph-1,2) = ring(0:nph-1)
          call spinring_synthesis(H,nlmax, b_south_U, nph, ring, kphi0,mmax_ring)
          map2(H%istart_south(ith):H%istart_south(ith)+nph-1,3) = ring(0:nph-1)
       endif

    enddo    ! loop on cos(theta)


    !     --------------------
    !     free memory and exit
    !     --------------------
    call HealpixFreeRecfac(H)
    DEALLOCATE(lam_fact)
    DEALLOCATE(normal_l)
    DEALLOCATE(b_north_Q,b_north_U)
    DEALLOCATE(b_south_Q,b_south_U)

    deallocate(TEB)

#ifdef MPIPIX
    if(DebugMsgs>1) print *,code//' Gather ',H%MpiID
    StartTime = Getetime()
    do i=1,3
    call MPI_GATHERV(map2(H%North_Start(H%MpiID),i),H%North_Size(H%MpiID),SP_MPI, &
       map_TQU(:,i),H%North_Size,H%North_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    call MPI_GATHERV(map2(H%South_Start(H%MpiID),i),H%South_Size(H%MpiID),SP_MPI, &
       map_TQU(:,i),H%South_Size,H%South_Start,SP_MPI, 0 ,MPI_COMM_WORLD, ierr)
    end do
    if(DebugMsgs>1) print *,code //' Done Gather ',H%MpiID, Getetime()-StartTime
    if (DebugMsgs>0 .and. H%MpiID==0) print *,code // ' Time: ', GeteTime() - iniTime
    deallocate(map2) 
#endif

  end subroutine polalm2map



 subroutine PackAlm2Alm(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(Out), DIMENSION(1:1,0:nlmax,0:nlmax) :: almout
   COMPLEX(SPC), INTENT(IN), DIMENSION(((nlmax+1)*(nlmax+2))/2) :: almin

   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(1,l,m) = almin(a_ix)
     end do
    end do  

 end  subroutine PackAlm2Alm

 subroutine Alm2PackAlm(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(in), DIMENSION(1:1,0:nlmax,0:nlmax) :: almin
   COMPLEX(SPC), INTENT(out), DIMENSION(((nlmax+1)*(nlmax+2))/2) :: almout
   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(a_ix) = almin(1,l,m) 
     end do
    end do  

 end  subroutine Alm2PackAlm


 subroutine PackEB2EB(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(Out), DIMENSION(1:2,0:nlmax,0:nlmax) :: almout
   COMPLEX(SPC), INTENT(IN), DIMENSION(1:2,((nlmax+1)*(nlmax+2))/2) :: almin

   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(:,l,m) = almin(:,a_ix)
     end do
    end do  

 end  subroutine PackEB2EB

 subroutine EB2PackEB(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(in), DIMENSION(1:2,0:nlmax,0:nlmax) :: almin
   COMPLEX(SPC), INTENT(out), DIMENSION(1:2,((nlmax+1)*(nlmax+2))/2) :: almout
   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(:,a_ix) = almin(:,l,m) 
     end do
    end do  

 end  subroutine EB2PackEB


 subroutine TEB2PackTEB(almin,TEBout,nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(in), DIMENSION(1:3,0:nlmax,0:nlmax) :: almin
   COMPLEX(SPC), INTENT(out), DIMENSION(1:3,((nlmax+1)*(nlmax+2))/2) :: TEBout
  
   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      TEBout(:,a_ix) = almin(:,l,m) 
     end do
    end do  

 end  subroutine TEB2PackTEB


subroutine PackTEB2TEB(almin,almout, nlmax)
   integer, intent (in) :: nlmax
   COMPLEX(SPC), INTENT(Out), DIMENSION(1:3,0:nlmax,0:nlmax) :: almout
   COMPLEX(SPC), INTENT(IN), DIMENSION(1:3,((nlmax+1)*(nlmax+2))/2) :: almin

   integer a_ix, m, l

    a_ix = 0
    do m = 0, nlmax
     do l=m, nlmax
      a_ix = a_ix + 1
      almout(:,l,m) = almin(:,a_ix)
     end do
    end do  

 end  subroutine PackTEB2TEB


#ifdef MPIPIX

  subroutine MessageLoop(H)
   use MPIStuff
    Type (HealpixInfo) :: H
    character (LEN=64) :: Msg
    REAL(SP),   DIMENSION(1) :: dummymap
    REAL(SP),   DIMENSION(1,3) :: dummymapTQU
    COMPLEX(SPC),   DIMENSION(1) :: dummymapC
    COMPLEX(SPC),   DIMENSION(1,1,1)  :: dummyalm
    integer i

     do
      Msg = ''
      call MPI_BCAST(Msg,64,MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)  !***

      if (DebugMsgs>1) print *,'Got message ', H%MpiID, ' '//trim(Msg)

      if (Msg=='EXIT') return
      if (Msg == 'SCALALM2MAP') then
         call scalalm2map(H,i,dummyalm, dummymap)
      else if (Msg == 'ALM2GRADIENTMAP') then
         call alm2GradientMap(H, i , dummyalm, dummymapC)
      else if (Msg == 'SPINALM2MAP') then
        call SpinAlm2Map(H,i,dummyalm, dummymapC, 1)
      else if (Msg == 'MAP2SCALALM') then
        call map2scalalm(H, i, dummymap, dummyalm, -1.d0)   
      else if (Msg=='MAP2SPINALM') then
        call map2spinalm(H, i, dummymapC, dummyalm, 0, -1.d0)   
      else if (Msg=='SCALALM2LENSEDMAP') then
         call scalalm2LensedMap(H,i,dummyalm, dummymapC, dummymap)
      else if (Msg=='ALM2LENSEDMAP') then
         call alm2LensedMap(H,i,dummyalm, dummymapC, dummymapTQU)
      else if (Msg=='ALM2LENSEDMAPINTERP') then
         call alm2LensedmapInterp(H,i,dummyalm, dummymapC, dummymapTQU)
      else if (Msg=='MAP2POLALM') then
         call Map2PolAlm(H,i,dummymapTQU, dummyalm, -1.d0)
      else if (Msg=='POLALM2MAP') then
         call polalm2map(H,i,dummyalm, dummymapTQU)
       else if (Msg=='ALM2LENSEDQUADCONTRIB') then
          call alm2LensedQuadContrib(H, i, dummyalm, dummymapC, dummymap)
      end if
     end do 
   
  end subroutine MessageLoop
 
  subroutine SendMessages(H, MsgIn)
   use MPIStuff
   Type (HealpixInfo) :: H
   CHARACTER(LEN=*), intent(in) :: MsgIn
   CHARACTER(LEN=64) :: Msg
    
    Msg = MsgIn
    if (DebugMsgs>1) print *,'Send messages '//trim(Msg)
    call MPI_BCAST(Msg,64,MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr) !***
   ! same time as BCAST in MessaegLoop

  end  subroutine SendMessages

#endif
end module spinalm_tools
