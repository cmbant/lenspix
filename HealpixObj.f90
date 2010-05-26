module HealpixObj
 use healpix_types, ONLY: SP,DP,I4B,SPC
 USE head_fits, ONLY : add_card, get_card
 USE pix_tools, ONLY :  npix2nside, nside2npix, query_disc
 USE fitstools, ONLY : getsize_fits, input_map, read_par, read_dbintab, write_asctab, &
   dump_alms,write_bintab
 USE spinalm_tools
 implicit none

 REAL(SP) :: fmissval = -1.6375e-30

 integer, parameter :: ord_ring = 1, ord_nest = 2

 integer, parameter ::  C_T = 1, C_E = 2, C_B = 3, C_C = 4

 logical :: RandInited = .false.

 integer, parameter :: nospinmap = 0

 type HealpixMap

  integer(I4B) npix, nmaps, ordering, nside, type
  REAL(SP), DIMENSION(:,:),   POINTER :: TQU
  COMPLEX(SP), dimension(:), pointer :: SpinField
!  REAL(SP), DIMENSION(:), pointer :: SpinQ, SpinU
  REAL(SP), dimension(:), pointer :: Phi
  integer :: spin
  logical :: HasPhi
 end type HealpixMap


 type HealpixAlm

   integer(I4B) lmax, npol, spin
   logical HasPhi
   COMPLEX(KIND=SPC), DIMENSION(:,:,:), pointer :: TEB    !T,E,B index, l,m 
   COMPLEX(KIND=SPC), DIMENSION(:,:,:), pointer :: SpinEB  ! E,B index, l, m
   COMPLEX(KIND=SPC), DIMENSION(:,:,:), pointer :: Phi !Has to be 3D array for Healpix defs

 end type HealpixAlm

 type HealpixPower
   !Raw Cls in milliK^2 units
   !read from text files in l(l+1)C_l/(2pi microK^2)
   !PhiCl is lensing potential phi-phi and phi-T. 
   !Phi is read in above units, but stored here as dimensionless
   integer(I4B) lmax
   logical pol, lens
   REAL(SP), DIMENSION(:,:),   POINTER :: Cl
   REAL(SP), DIMENSION(:,:),  POINTER :: PhiCl  

 end  type HealpixPower

 type HaarComponents
   
   integer order
   logical pol
   type(HealpixMap) degraded
   type(HealpixMap), dimension(:), pointer :: details
  
 end  type HaarComponents

contains

   subroutine HealpixPower_Init(P, almax, pol, dolens)
    Type(HealpixPower) P
    integer, intent(in) :: almax
    logical, intent(in) :: pol
    logical, intent(in), optional :: dolens
      
     call HealpixPower_Free(P)
     P%lmax = almax
     P%pol = pol
     P%lens = .false.
     if (present(dolens)) then
       P%lens = dolens
     end if
     if (pol) then
       allocate(P%Cl(0:almax,4))  
     else
       allocate(P%Cl(0:almax,1))  
     end if
     if (P%lens) then
      allocate(P%PhiCl(0:almax,2))      
     end if

   end subroutine HealpixPower_Init

   subroutine HealpixPower_Free(P)
    Type(HealpixPower) P
    integer status

    deallocate(P%Cl, stat = status)
    deallocate(P%PhiCl, stat = status)
    nullify(P%Cl, P%PhiCl)

   end subroutine HealpixPower_Free


   subroutine HealpixPower_ReadFromTextFile(P, f, lmax, pol, dolens)
     use AMLutils
     Type(HealpixPower) P
     character(LEN=*), intent(IN) :: f
     integer, intent(in) :: lmax
     logical,intent(in), optional :: dolens, pol
     integer i,l, ColNum
     real(DP) T,E,B,TE, scal, phi, phiT
     logical tensors, dolens2, pol2
     integer, parameter :: io_unit= 1
     character(LEN=200) :: InLine
     real(sp) test(8)

     open(unit=io_unit,file=f,form='formatted',status='old')

     call HealpixPower_Free(P)
     dolens2=.false.
     if (present(dolens)) then
        dolens2 = dolens
     end if
     pol2 = .false.
     if (present(pol)) then
        pol2 = pol
     end if
     call HealpixPower_Init(P, lmax, pol = pol2, dolens = dolens2)
 
    !See how many columns - includes magnetic polarization if four columns
      read(io_unit,'(a)') InLine
      do l=8,4,-1
         Colnum=l 
         read(InLine,*, end=110) test(1:l)
         exit
110      cycle

      end do 
      
      tensors = colnum==5 .or. colnum==7
      if (P%lens .and. colnum /= 6 .and. colnum/=7) stop 'can''t indentify text phi cl columns'
       
      rewind io_unit

      B=0
      P%Cl = 0
  
      if (dolens2) P%PhiCl = 0
        
      do i=0,lmax
       if (tensors) then
         if (P%lens) then
          read (io_unit,*,end=118) l, T, E, B , TE, phi, phiT
         else
          read (io_unit,*,end=118) l, T, E, B , TE
         end if
       else
         if (P%lens) then
          read (io_unit,*,end=118) l, T, E, TE, phi, phiT
         else
          read (io_unit,*,end=118) l, T, E, TE
         end if
       end if
       if (l<= lmax .and. l>=2) then
         scal=twopi/(l*(l+1))/1e6
         P%Cl(l,1) = T*scal
         if (pol2) then
         P%Cl(l,2) = E*scal
         P%Cl(l,3) = B*scal
         P%Cl(l,4) = TE*scal
         end if
         if (P%Lens) then
             P%PhiCl(l,1) = phi/real(l,dp)**4/1e12/2.726**2
             P%PhiCl(l,2) = phiT/real(l,dp)**3/1e9/2.726
         end if
       end if
       enddo
     118 close(io_unit)

  end subroutine HealpixPower_ReadFromTextFile


   subroutine HealpixPower_Write(P,fname) 
     use AMLutils
     Type(HealpixPower) P
     character(Len=*), intent(in) :: fname
     integer l

     call CreateTxtFile(fname,1)
     do l=0,P%lmax
       if (P%pol) then
         write (1,'(1I7,4E15.5)') l,l*(l+1)*P%Cl(l,:)/twopi *1e6
       else
         write (1,'(1I7,1E15.5)') l,l*(l+1)*P%Cl(l,1)/twopi *1e6
       end if
     end do

     close(1)

   end  subroutine HealpixPower_Write

  subroutine HealpixPower_Smooth(P,fwhm, sgn)
   Type(Healpixpower) :: P
   integer l, sn
   integer, intent(in), optional :: sgn
   real(dp) xlc,sigma2,fwhm

   if (present(sgn)) then
    sn = sgn
   else
    sn = -1
   end if
   
   xlc= 180*sqrt(8.*log(2.))/3.14159
   sigma2 = (fwhm/xlc)**2

   do l=2,P%lmax
     P%Cl(l,:) =  P%Cl(l,:)*exp(sn*l*(l+1)*sigma2)
   end do

  end subroutine HealpixPower_Smooth

   subroutine HealpixAlm2Power(A,P)
    Type(HealpixPower) P
    Type(HealpixAlm) :: A
    integer l,i,ix

    call HealpixPower_Init(P,A%lmax,A%npol==3)
    
    P%Cl(0:1,:) = 0
    do l=0, P%lmax
     if (l<2) then
     ix= 1
    else
     ix = A%npol
     end if
     do i = 1, ix
      P%Cl(l,i) = ( REAL(A%TEB(i,l,0))**2 &
            + 2.*SUM(A%TEB(i,l,1:l)*CONJG(A%TEB(i,l,1:l)) )) / (2.*l + 1.)
     end do
     if (ix==3) then
        P%Cl(l,4) = ( REAL(A%TEB(1,l,0))*REAL(A%TEB(2,l,0)) &
            + 2.*SUM(real(A%TEB(1,l,1:l)*CONJG(A%TEB(2,l,1:l))) ) &
           ) / (2.*l + 1.)
      end if
    end do 

   end subroutine HealpixAlm2Power


   subroutine HealpixAlm_Init(A,almax,npol,spinmap, HasPhi)
     Type(HealpixAlm) :: A
     integer, intent(in) :: almax
     integer, intent(in), optional :: npol, spinmap
     logical, intent(in), optional :: HasPhi
     integer status

    if (almax > 3000) stop 'HealpixAlm_Init: Not sure can handle that large l_max!'

     call HealpixAlm_Free(A)

     A%lmax = almax
     if (present(npol)) then
      A%npol = npol
     else
      A%npol = 1
     end if

     if (A%npol /= 0) then 
      ALLOCATE(A%TEB(1:A%npol, 0:almax, 0:almax),stat = status)
      if (status /= 0) stop 'No Mem: HealpixAlm_Init'
      A%TEB=0
     end if

     if (present(spinmap)) then
         if (spinmap /= nospinmap) then
          if (spinmap<1 .or. spinmap > 3) stop 'Spin must be 0<spin<4'
          ALLOCATE(A%SpinEB(2, 0:almax, 0:almax),stat = status)
          if (status /= 0) stop 'No Mem: HealpixAlm_Init'
         end if
         A%spin= spinmap     
     else
      A%spin = nospinmap      
     end if

     if (present(HasPhi)) then
       A%HasPhi = HasPhi
     else
       A%HasPhi = .false.   
     end if

     if (A%HasPhi) then
          ALLOCATE(A%Phi(1:1,0:almax, 0:almax),stat = status)
          A%Phi = 0
          if (status /= 0) stop 'No Mem: HealpixAlm_Init'
     end if

   end subroutine HealpixAlm_Init


  subroutine HealpixAlm_Assign(AOut, Ain)
   Type(HealpixAlm) :: AOut, Ain
   integer status

   call HealpixAlm_Free(AOut)
   Aout = Ain 
   nullify(AOut%TEB, AOut%SpinEB, AOut%Phi)
   if (Ain%npol>0) then
      ALLOCATE(AOut%TEB(1:AOut%npol, 0:AOut%lmax, 0:AOut%lmax),stat = status)
      if (status /= 0) stop 'No Mem: HealpixAlm_Assign'
     AOut%TEB= Ain%TEB
   end if
   if (AIn%spin /= nospinmap) then
        ALLOCATE(AOut%SpinEB(2, 0:AOut%lmax, 0:AOut%lmax),stat = status)
        if (status /= 0) stop 'No Mem: HealpixAlm_Assign'
        AOut%SpinEB = Ain%SpinEB
   end if
   if (AIn%HasPhi) then
       ALLOCATE(AOut%Phi(1:1,0:AOut%lmax, 0:AOut%lmax),stat = status)
      if (status /= 0) stop 'No Mem: HealpixMap_Assign'
      AOut%Phi = Ain%Phi   
   end if 
  
  end subroutine HealpixAlm_Assign


  subroutine HealpixAlm_Free(A)
   Type(HealpixAlm) :: A
   integer status

     deallocate(A%TEB,stat=status)
     deallocate(A%SpinEB,stat=status)
     nullify(A%TEB)
     nullify(A%SpinEB)
     deallocate(A%Phi,stat=status)
     nullify(A%Phi)

  end subroutine HealpixAlm_Free

  subroutine HealpixAlm_PhiOnly(A)
   Type(HealpixAlm) :: A
   integer status

     deallocate(A%TEB,stat=status)
     deallocate(A%SpinEB,stat=status)
     nullify(A%TEB)
     nullify(A%SpinEB)
     A%spin = nospinmap      
     A%npol = 0 
  end subroutine 


  subroutine HealpixAlm_GradientOf(A, B, field, updown)
   type(HealpixAlm) :: A, B
   integer l
   character(LEN=*), intent(in) :: field
   character(LEN=*), intent(in), optional :: updown
   logical Div  
   integer spin
   
   if (field(1:1) /= 'S') then
    spin = 1
   else
    if (.not. present(updown))  stop 'HealpixAlm_GradientOf: Must say which derivative'
    Div = updown(1:1) == 'D'

    if (Div) then
  !    spin = A%spin-1
      spin = 1
    else
!     spin = A%spin+1
      spin = 3
    end if
   end if


   call HealpixAlm_Init(B,A%lmax,0,spinmap= spin)

   B%SpinEB = 0
    do l=B%spin, A%lmax
      if (field(1:1)=='P') then
       B%SpinEB(1,l,0:l) = - EB_sign*sqrt(real(l*(l+1),dp))*A%Phi(1,l,0:l)
      else if (field(1:1)=='T') then
       B%SpinEB(1,l,0:l) = - EB_sign*sqrt(real(l*(l+1),dp))*A%TEB(1,l,0:l)
      else if (field(1:1)=='S') then
         if (Div) then
        !Divergence
          B%SpinEB(:,l,0:l) =  sqrt(real((l+2)*(l-1),dp))*A%TEB(2:3,l,0:l)
          else
          !STF of outer product
         B%SpinEB(:,l,0:l) =  -sqrt(real((l+3)*(l-2),dp))*A%TEB(2:3,l,0:l)
         end if
      else
       stop 'HealpixAlm_GradientOf: Unknown field'
      end if 
    end do

  end subroutine HealpixAlm_GradientOf

  subroutine HealpixAlm_PolToSpin2(A)
   Type(HealpixAlm) :: A
   integer status

   deallocate(A%SpinEB, stat=status)
   ALLOCATE(A%SpinEB(2, 0:A%lmax, 0:A%lmax),stat = status)
   A%SpinEB(1:2,:,:) = A%TEB(2:3,:,:)
   A%spin = 2
  end   subroutine HealpixAlm_PolToSpin2

  subroutine HealpixAlm_Spin2ToPol(A)
   Type(HealpixAlm) :: A
   integer :: status

   if (A%npol==0) then
      A%npol = 3
      ALLOCATE(A%TEB(1:A%npol, 0:A%lmax, 0:A%lmax),stat = status)
      if (status /= 0) stop 'No Mem: HealpixAlm_Spin2ToPol'
      A%TEB=0
   end if
   A%TEB(2:3,:,:) = A%SpinEB(1:2,:,:)

  end subroutine HealpixAlm_Spin2ToPol



  subroutine HealpixAlm_Smooth(A,fwhm, sgn)
   Type(HealpixAlm) :: A
   integer l, sn
   integer, intent(in), optional :: sgn
   real(dp) xlc,sigma2,fwhm
   
   if (present(sgn)) then
     sn = sgn
   else
     sn = -1
   end if
   xlc= 180*sqrt(8.*log(2.))/3.14159
   sigma2 = (fwhm/xlc)**2

   do l=2,A%lmax
     A%TEB(:,l,:) =  A%TEB(:,l,:)*exp(sn*l*(l+1)*sigma2/2)
   end do

  end subroutine HealpixAlm_Smooth

  subroutine HealpixMap_GetAzimCut(M, npix,rad, theta,phi)
   !1 inside disc radius rad centred at theta, phi (radians)
    use pix_tools
    Type(HealpixMap) :: M
    real(dp), intent(in) :: rad, theta, phi
    integer, intent(in) :: npix
    real(dp) vec(3)
    
    call ang2vec(theta,phi,vec)
    call HealpixMap_GetAzimCutVec(M, npix, rad, vec)
    
  end subroutine HealpixMap_GetAzimCut

  subroutine HealpixMap_GetAzimCutVec(M, npix,rad, Vec)
   !inside disc radius rad centred at vec
    Type(HealpixMap) :: M
    real(dp), intent(in) :: rad
    integer, intent(in) :: npix
    integer, dimension(:), allocatable :: listpix
    real(dp), intent(in) :: vec(3)
    integer nlist,i
  
    call HealpixMap_Init(M,npix,1)
    allocate(listpix(0:npix-1))
    call query_disc(M%nside,vec, rad, listpix,nlist)
    M%TQU = 0
    M%TQU(listpix(0:nlist-1),1) = 1
    deallocate(listpix)     

  end subroutine HealpixMap_GetAzimCutVec


  function HealpixMap_Vec2pix(M, vec) result(pix)
    use pix_tools
    Type(HealpixMap) :: M 
    integer :: pix
    real(dp), intent(in) :: vec(3)
    
    if (M%Ordering == ord_ring) then
     call vec2pix_ring(M%nside,vec,pix) 
    else
     call vec2pix_nest(M%nside,vec,pix) 
    end if

  end function HealpixMap_Vec2pix


  subroutine HealpixMap_Pix2Vec(M, pix,vec)
    use pix_tools
    Type(HealpixMap) :: M
 
    integer, intent(in) :: pix
    real(dp), intent(out) :: vec(3)
    
    if (M%Ordering == ord_ring) then
     call pix2vec_ring(M%nside,pix, vec) 
    else
     call pix2vec_nest(M%nside,pix, vec) 
    end if

  end subroutine HealpixMap_Pix2Vec

  subroutine HealpixMap_Pix2Vertex(M, pix,vertex)
    use pix_tools
    Type(HealpixMap) :: M
 
    integer, intent(in) :: pix
    real(dp), intent(out) :: vertex(3,4)
    real(dp) :: vec(3)
    
    if (M%Ordering == ord_ring) then
     call pix2vec_ring(M%nside,pix, vec, vertex) 
    else
     call pix2vec_nest(M%nside,pix, vec, vertex) 
    end if

  end subroutine HealpixMap_Pix2Vertex


  function HealpixMap_Ang2Pix(M, theta, phi) result(pix)
    use pix_tools
    Type(HealpixMap) :: M
    real(dp), intent(in):: theta, phi
    integer pix
    
    if (M%Ordering == ord_ring) then
     call Ang2Pix_ring(M%nside,theta,phi,pix) 
    else
     call Ang2Pix_nest(M%nside,theta,phi,pix) 
    end if

  end function HealpixMap_Ang2Pix

  subroutine Healpix_GetRotation(R, theta, phi, chi)
!Rotate phi about z axis, then rotation by theta about new y axis, then chi about new z axis
    real(dp), intent(in) :: theta, phi, chi 
    real(dp), intent(inout) :: R(3,3)
  
    stop 'Healpix_GetRotation: You''ll have to check this routine'

    R(1,1) = cos(phi)*cos(theta)*cos(chi) - sin(phi)*sin(chi)
    R(1,2) = sin(phi)*cos(theta)*cos(chi) + cos(phi)*sin(chi)
    R(1,3) = -sin(theta)*cos(chi)
    R(2,1) = -cos(phi)*cos(theta)*sin(chi) - sin(phi)*cos(chi)
    R(2,2) = -sin(phi)*cos(theta)*sin(chi) + cos(phi)*cos(chi)
    R(2,3) = sin(phi)*sin(chi)
    R(3,1) = cos(phi)*sin(theta)
    R(3,2) = sin(phi)*sin(theta)
    R(3,3) = cos(theta)

  end subroutine Healpix_GetRotation

 
  subroutine HealpixMap_Rotate(M, MR, theta, phi, chi)
!This is very crude for temp
   use pix_tools
    Type(HealpixMap) :: M, MR
    real(dp), intent(in) :: theta, phi, chi 
    real(dp) t,p
    integer i, ix
    real(dp) vec(3), R(3,3)

    stop 'Don''t use this'
    call Healpix_GetRotation(R, theta, phi, chi)
    call HealpixMap_Init(MR, M%npix, M%nmaps)
    call HealpixMap_ForceRing(M)
    do i=0, MR%npix -1
      call pix2vec_ring(MR%nside, i, vec)
      vec = matmul(R,vec)
      call vec2pix_ring(MR%nside, vec, ix)
      MR%TQU(i,:) = M%TQU(ix,:)
    end do

  end subroutine HealpixMap_Rotate

  subroutine HealpixAlm_Sim(A, P, seed, HasPhi, dopol)
   use random
   use alm_tools
   use ran_tools
   Type(HealpixAlm) :: A
   Type(HealpixPower) :: P
   integer, intent(in), optional :: seed
   logical, intent(in), optional :: HasPhi,dopol
   integer l,m
   logical wantphi
   integer wantpol
   real(sp) xamp, corr, tamp, Bamp

   if (present(seed)) then
      call InitRandom(seed)
   else
     if (.not. RandInited) call InitRandom
     RandInited = .true.
   end if

   wantpol = 1

   if (present(dopol)) then
     if (dopol) wantpol = 3
   end if

   if (present(HasPhi)) then
     wantphi= HasPhi
     if (wantphi .and. .not. associated(P%PhiCl)) stop 'HealpixAlm_Sim: PhiCl not present'
   else
     wantphi = .false.
   end if

   call HealpixAlm_Init(A,P%lmax, wantpol, HasPhi=wantphi)

   A%TEB=0
   do l=2, P%lmax
      A%TEB(1,l,0) =Gaussian1()* sqrt(P%Cl(l,1))
      tamp = sqrt(P%Cl(l,1)/2)
      do m = 1, l
       A%TEB(1,l,m) =cmplx(Gaussian1(),Gaussian1())*tamp
      end do 
   end do

    if (wantpol >= 3) then  
      !polarization, E corrolated to T

      do l = 2, P%lmax 
        if (p%cl(l,C_t) == 0) then
          tamp = 1.0  !Prevent divide by zero - TE should also be zero
        else
          tamp = p%cl(l,C_T)
        end if
        corr = p%cl(l,C_C)/tamp
        xamp = sqrt(P%cl(l,C_E) - corr*P%cl(l,C_C))
        Bamp = sqrt(P%cl(l,C_B))
        A%TEB(2,l,0) = Corr*A%TEB(1,l,0) + Gaussian1()*xamp
        A%TEB(3,l,0)=  Bamp*Gaussian1()
        xamp = xamp / sqrt(2.0)
        Bamp = Bamp /sqrt(2.0)
        do m =1, l
          A%TEB(2,l,m) = corr*A%TEB(1,l,m) + cmplx(Gaussian1(),Gaussian1())*xamp
          A%TEB(3,l,m) = Bamp*cmplx(Gaussian1(),Gaussian1())      
        end do

      end do

    end if 
   if (wantphi) then
     !Make phi with correct correlation to T
     A%Phi=0
     do l=2, P%lmax
      if (P%Cl(l,1)==0) then
        tamp = 1.0
      else
        tamp=P%Cl(l,1)
      end if
      corr = P%PhiCl(l,2)/tamp
      xamp = sqrt(max(0._sp,P%PhiCl(l,1) - corr*P%PhiCl(l,2)))
      A%Phi(1,l,0) = corr*A%TEB(1,l,0) + Gaussian1()*xamp
      xamp=  xamp/sqrt(2.0)
      do m = 1, l
        A%Phi(1,l,m) = corr*A%TEB(1,l,m) + cmplx(Gaussian1(),Gaussian1())*xamp
      end do
     end do 

   end if

  end subroutine HealpixAlm_Sim



  subroutine HealpixAlm_SimPhi(A, P, seed)
   use random
   use alm_tools
   use ran_tools
   Type(HealpixAlm) :: A
   Type(HealpixPower) :: P
   integer, intent(in), optional :: seed
   
   integer l,m

   if (present(seed)) then
      call InitRandom(seed)
   else
     if (.not. RandInited) call InitRandom
     RandInited = .true.
   end if
   call HealpixAlm_Init(A,P%lmax, 0,HasPhi = .true.)
   if (.not. P%lens) stop 'must have phi power spectrum'

   A%Phi(:,0:1,:)=0 !So what about the dipole??
   do l=2, P%lmax
      A%Phi(1,l,0) =Gaussian1()* sqrt(P%PhiCl(l,1))
      do m = 1, l
       A%Phi(1,l,m) =cmplx(Gaussian1(),Gaussian1())* sqrt(P%PhiCl(l,1)/2)
      end do 
   end do

  end subroutine HealpixAlm_SimPhi


  subroutine HealpixMap_Init(M, npix, nmaps, nested, spinmap, HasPhi)

    Type(HealpixMap) :: M
    integer, intent(in) :: npix 
    integer, intent(in), optional :: nmaps, spinmap
    logical, intent(in), optional :: nested, HasPhi
    integer status

     M%npix = npix

     call HealpixMap_Free(M)
        
     if (present(nmaps)) then
       M%nmaps = nmaps
     else
       M%nmaps = 1
     end if

     if (present(HasPhi)) then
       M%HasPhi = HasPhi
     else
       M%HasPhi = .false.
     end if

     if (M%nmaps /= 0) call HealpixMap_AllocateTQU(M,M%nmaps) 
 
     M%ordering = ord_ring
     if (present(nested)) then
        if (nested) M%ordering = ord_nest
     end if
  
     M%nside = npix2nside(npix)
     if (present(spinmap)) then
      M%spin = spinmap
     else
      M%spin = nospinmap
     end if
     if (M%spin /= nospinmap) then
       if (M%spin<1 .or. M%spin > 3) stop 'Spin must be 0<spin<4'
       ALLOCATE(M%SpinField(0:M%npix-1),stat = status)
       M%spin = spinmap
     end if
    
    if (M%HasPhi) call HealpixMap_AllocatePhi(M)

  end subroutine HealpixMap_Init


  subroutine HealpixMap_AllocatePhi(M)
    Type(HealpixMap) :: M
    integer status
 
      ALLOCATE(M%Phi(0:M%npix-1),stat = status)
      if (status /=0) stop 'HealpixMap_AllocatePhi: allocate'
      M%HasPhi = .true.

  end  subroutine HealpixMap_AllocatePhi

  subroutine HealpixMap_AllocateTQU(M, nmaps)
    Type(HealpixMap) :: M
    integer, intent(in) :: nmaps
    integer status
 
      M%nmaps = nmaps
      deallocate(M%TQU, stat=status)
      ALLOCATE(M%TQU(0:M%npix-1,M%nmaps),stat = status)
      if (status /= 0) stop 'HealpixMap_AllocateTQU: allocate'
      M%TQU= 0
 
  end  subroutine HealpixMap_AllocateTQU


  subroutine HealpixMap_DeAllocateTQU(M)
    Type(HealpixMap) :: M
    integer status
 
    deallocate(M%TQU, stat=status)
    nullify(M%TQU)
    M%nmaps = 0

  end subroutine HealpixMap_DeAllocateTQU

  subroutine HealpixMap_PhiOnly(M)
    Type(HealpixMap) :: M
    integer status
       
   deallocate(M%TQU, stat=status)
   nullify(M%TQU)
   deallocate(M%SpinField, stat=status)
   nullify(M%SpinField)
   M%spin = nospinmap
   M%nmaps = 0   

  end subroutine HealpixMap_PhiOnly
  

  subroutine HealpixMap_Assign(MOut, Min)
   Type(HealpixMap) :: MOut, Min
   integer status

   call HealpixMap_Free(MOut)
   Mout = Min 
   nullify(MOut%TQU)
   if (Min%nmaps>0) then
     allocate(MOut%TQU(0:Min%npix-1,Min%nmaps),stat = status)
     if (status /= 0) stop 'No Mem: HealpixMap_Assign'
     MOut%TQU = Min%TQU
   end if
   nullify(MOut%SpinField,MOut%phi)
   if (MIn%spin /= nospinmap) then
     allocate(MOut%SpinField(0:Min%npix-1),stat = status)
     if (status /= 0) stop 'No Mem: HealpixMap_Assign'
     MOut%SpinField = Min%SpinField
     !MOut%SpinQ => MOut%SpinField(:,1)
     !MOut%SpinU => MOut%SpinField(:,2) 
   end if
   if (MIn%HasPhi) then
     allocate(MOut%Phi(0:Min%npix-1),stat = status)
     if (status /= 0) stop 'No Mem: HealpixMap_Assign'
     MOut%Phi = Min%Phi   
   end if 

  end subroutine HealpixMap_Assign


  subroutine HealpixMap_Read(OutMAP,fname)
   CHARACTER(LEN=80), DIMENSION(1:120) :: header_in
   character(LEN=*), intent(in) :: fname
   integer status, j
   
   Type(HealpixMap) OutMap 

   call HealpixMap_Free(OutMap)
    
   OutMap%npix = getsize_fits(fname, nmaps=OutMap%nmaps, ordering=OutMap%ordering, nside=OutMap%nside,&
        type=OutMap%type)
   if ((OutMap%ordering /=  ord_ring).and.(OutMap%ordering /= ord_nest)) then
     PRINT*,'The ordering scheme of the map must be RING or NESTED.'
     PRINT*,'No ordering specification is given in the FITS-header!'
     stop 1
    endif
  if (OutMap%nside /= npix2nside(OutMap%npix)) then ! accept only full sky map
     print*,'FITS header keyword NSIDE = ',OutMap%nside,' does not correspond'
     print*,'to the size of the map!'
     stop 1
   endif

   call HealpixMap_AllocateTQU(OutMap,OutMap%nmaps) 
   OutMap%HasPhi = .false.
   OutMap%spin = nospinmap
 
   call input_map(fname, OutMAP%TQU, OutMap%npix, OutMap%nmaps, &
       &   fmissval=fmissval, header= header_in)

!!To do, boring...
  ! do j=1,nmaps
  !   call get_card(header_in, trim(numcat('TTYPE',j)), ttype(j))
  !   call get_card(header_in, trim(numcat('TUNIT',j)), tunit(j))
  ! enddo

  end subroutine HealpixMap_Read

  subroutine HealpixMap_Write(M, fname)
    Type(HealpixMap), intent(in) :: M
    character(LEN=*), intent(in) :: fname
    CHARACTER(LEN=80), DIMENSION(1:120) :: header
    integer j, nlheader

    header = ' '
    call add_card(header,'COMMENT','-----------------------------------------------')
    call add_card(header,'COMMENT','     Sky Map Pixelisation Specific Keywords    ')
    call add_card(header,'COMMENT','-----------------------------------------------')
    call add_card(header,'PIXTYPE','HEALPIX','HEALPIX Pixelisation')
    if (M%ordering == ord_ring) then
     call add_card(header,'ORDERING','RING',  'Pixel ordering scheme, either RING or NESTED')
    else
     call add_card(header,'ORDERING','NESTED',  'Pixel ordering scheme, either RING or NESTED')
    endif
    call add_card(header,'NSIDE'   ,M%nside,   'Resolution parameter for HEALPIX')
    call add_card(header,'FIRSTPIX',0,'First pixel # (0 based)')
    call add_card(header,'LASTPIX',M%npix-1,'Last pixel # (0 based)')
    call add_card(header) ! blank line
    call add_card(header,'CREATOR','HEALPixObj',        'Software creating the FITS file')
    call add_card(header,'INDXSCHM','IMPLICIT',' Indexing : IMPLICIT or EXPLICIT')
    call add_card(header,'GRAIN', 0, ' Grain of pixel indexing') ! full sky
    if (M%nmaps == 3) then
     call add_card(header,'POLAR',.true.," Polarisation included (True/False)")
    else
     call add_card(header,'POLAR',.false.," Polarisation included (True/False)")
    endif
    call add_card(header) ! blank line
    call add_card(header,"TTYPE1", "TEMPERATURE","Temperature map")
    call add_card(header,"TUNIT1", "mK", "map unit")
    call add_card(header)
    if (M%nmaps == 3) then
       call add_card(header,"TTYPE2", "Q-POLARISATION","Q Polarisation map")
       call add_card(header,"TUNIT2", "mK", "map unit")
       call add_card(header)
       call add_card(header,"TTYPE3", "U-POLARISATION","U Polarisation map")
       call add_card(header,"TUNIT3", "mK", "map unit")
       call add_card(header)
    endif
    call add_card(header,"COMMENT","*************************************")

   nlheader = SIZE(header)
    call write_bintab(M%TQU, M%npix, M%nmaps, header, nlheader, fname)

  end  subroutine HealpixMap_Write


  subroutine HealpixMap_Free(M)
   Type(HealpixMap) :: M
   integer status  

   call HealpixMap_DeAllocateTQU(M)

   deallocate(M%SpinField, M%Phi, stat=status)
   nullify(M%SpinField,M%Phi)

  end subroutine HealpixMap_Free

  subroutine HealpixMap_ForceRing(M)
    USE pix_tools, ONLY : convert_nest2ring
   Type(HealpixMap) :: M
   integer i

   if (M%ordering /= ord_ring) then
    do i=1, M%nmaps
      call convert_nest2ring (M%nside, M%TQU(:,i))
     end do
     if (M%spin/= nospinmap) then
       stop 'ring not done for pol'
!       call convert_nest2ring (M%nside, M%SpinField(:,1))
!       call convert_nest2ring (M%nside, M%SpinField(:,2))
      end if
      if (M%HasPhi) then
       call convert_nest2ring (M%nside, M%Phi)
      end if
   end if
   M%ordering = ord_ring

  end subroutine HealpixMap_ForceRing

  
  subroutine HealpixMap_ForceNest(M)
    USE pix_tools, ONLY : convert_ring2nest
   Type(HealpixMap) :: M
   integer i

   if (M%ordering /= ord_nest) then
     do i=1, M%nmaps
      call convert_ring2nest (M%nside, M%TQU(:,i))
     end do
     if (M%spin/= nospinmap) then
      stop 'nest not done for pol'
       !call convert_ring2nest (M%nside, M%SpinField(:,1))
       !call convert_ring2nest (M%nside, M%SpinField(:,2))
      end if

      if (M%HasPhi) then
       call convert_ring2nest (M%nside, M%Phi)
      end if
     end if
   M%ordering = ord_nest

  end subroutine HealpixMap_ForceNest


  subroutine HealpixMapMulCut(InMap,CutMap,OutMap, map_ix, missval)
    Type(HealpixMap), intent(in) :: InMap,CutMap
    Type(HealpixMap), intent(out) :: OutMap
    integer, intent(in), optional :: map_ix
    real(sp), intent(in), optional :: missval
    integer i, j, ix

   if (present(map_ix)) then
     ix = map_ix
   else
     ix = 2
   end if
   call HealpixMap_ForceRing(InMap)
   call HealpixMap_ForceRing(CutMap)
   if (InMap%npix /= CutMap%npix) stop 'Map size mismatch'
   call HealpixMap_Init(OutMap,InMap%npix,InMap%nmaps)
   outMap%ordering  = ord_ring
   if (CutMap%nmaps < ix ) stop 'not enough maps'
   if (present(missval)) then
 
      do j=0, InMap%npix -1
        if (CutMap%TQU(j,ix) == 0) then
          OutMap%TQU(j,:) = missval
        else       
          OutMap%TQU(j,:) = InMap%TQU(j,:)
        end if
     end do
 
   else

    do i=1,InMap%nmaps
     OutMap%TQU(:,i) = InMap%TQU(:,i) * CutMap%TQU(:,ix)
    end do

   end if
  end subroutine HealpixMapMulCut

   subroutine HealpixMap_MulCutFile(InMap,CutFile,OutMap)
    Type(HealpixMap), intent(in) :: InMap
    Type(HealpixMap), intent(out) :: OutMap
    character(len=*), intent(in) :: CutFile
    Type(HealpixMap) CutMap

    call HealpixMap_Read(CutMap,Cutfile)
    call HealpixMapMulCut(InMap,CutMap,OutMap)
    call HealpixMap_Free(CutMap)

   end subroutine HealpixMap_MulCutFile


   subroutine HealpixMap2alm(H, M,A, almax,theta_cut_deg,map_ix, dopol)

     Type (HealpixInfo) :: H
     Type(HealpixMap), intent(in) :: M
     integer, intent(in) :: almax
     real(dp), intent(in), optional :: theta_cut_deg
      integer, intent(in), optional :: map_ix
      logical, intent(in), optional :: dopol
     Type(HealpixAlm) :: A
     integer npol, ix
     real(dp) cos_theta_cut

     call HealpixMap_ForceRing(M)

     npol = 1
     if (present(dopol)) then
      if (dopol) npol =3
     else
      if (M%nmaps ==3) npol = 3
     end if

     if (M%nmaps ==0) npol = 0

     if (present(theta_cut_deg)) then
       cos_theta_cut =  SIN(theta_cut_deg/180.d0*3.14159265358)
       if (theta_cut_deg < 0) cos_theta_cut = -1
      else
       cos_theta_cut = -1
     end if

     call HealpixAlm_Init(A,almax, npol, M%spin, M%HasPhi)
 
     if (npol==1 .and. M%nmaps>0) then
      ix = 1
      if (present(map_ix)) ix = map_ix
      call map2scalalm(H, almax, M%TQU(:,ix), A%TEB,cos_theta_cut)
     else if (npol ==3 .and. M%nmaps>0) then
      if (present(map_ix)) stop ' cannot have polarization and multiple map indices'
      call map2polalm(H, almax, M%TQU, A%TEB,cos_theta_cut)
     end if  

     if (M%HasPhi) then
       call map2scalalm(H, almax, M%Phi, A%Phi,cos_theta_cut)
     end if

     if (M%spin /= nospinmap) then
       !Haven't computed ring weights for odd spins
        call map2spinalm(H,almax, M%SpinField,  A%SpinEB, M%spin,cos_theta_cut)
     end if  
 
   end subroutine HealpixMap2alm


  subroutine HealpixAlm2GradientMap(H, A, M, npix, What)
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: M
     Type(HealpixAlm), intent(in) :: A
     Type(HealpixAlm) :: AT
     integer, intent(in) :: npix
     character(LEN=*), intent(in) :: What
   
    call HealpixMap_Init(M,npix,nmaps = 0, spinmap = 1)
    if (What(1:1) == 'P') then
     if (.not. A%HasPhi) stop 'HealpixAlm2GradientMap: No phi field'
     call alm2GradientMap(H, A%lmax, A%Phi,M%SpinField)
    else if (What(1:1) == 'T') then
     call HealpixAlm_Init(AT, A%lmax,npol = 0, HasPhi = .true.)
     AT%Phi = A%TEB(1:1,:,:)
     call alm2GradientMap(H,A%lmax, AT%Phi,M%SpinField)
     call HealpixAlm_Free(AT)
    else
     stop 'HealpixAlm2GradientMap: unknown field'
    end if
  end subroutine HealpixAlm2GradientMap


  subroutine  HealpixExactLensedMap(H,A, M, npix)
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: GradPhi, M
     integer, intent(in) :: npix
     Type(HealpixAlm), intent(in) :: A
      
     call HealpixAlm2GradientMap(H,A,GradPhi,npix,'PHI')
     call HealpixExactLensedMap_GradPhi(H,A, GradPhi,M)
     call HealpixMap_Free(GradPhi)

  end subroutine  HealpixExactLensedMap

  subroutine  HealpixExactLensedMap_GradPhi(H,A, GradPhi,M)
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: GradPhi, M
     Type(HealpixAlm), intent(in) :: A
      
     call HealpixMap_Init(M,GradPhi%npix,nmaps = A%npol)
     if (GradPhi%spin /=1) stop 'HealpixExactLensedMap: GradPhi must be spin 1 field'
     if (A%npol ==1) then
      call scalalm2LensedMap(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU)
     else
      call alm2LensedMap(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU)
     end if
  end subroutine  HealpixExactLensedMap_GradPhi

  subroutine  HealpixInterpLensedMap(H,A, M, npix, factor)
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: GradPhi, M
     integer, intent(in) :: npix
     Type(HealpixAlm), intent(in) :: A
     integer, intent(in), optional :: factor
     integer fact

     fact = 8 !sufficient for 0.5% accuracy to l=2000
     if (present(factor)) fact = factor
     call HealpixAlm2GradientMap(H,A,GradPhi,npix,'PHI')
     call HealpixInterpLensedMap_GradPhi(H,A, GradPhi,M, fact)
     call HealpixMap_Free(GradPhi)

  end subroutine  HealpixInterpLensedMap

  subroutine  HealpixInterpLensedMap_GradPhi(H,A, GradPhi,M, factor)
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: GradPhi, M
     Type(HealpixAlm), intent(in) :: A
     integer, intent(in) :: factor
      
     call HealpixMap_Init(M,GradPhi%npix,nmaps = A%npol)
     if (GradPhi%spin /=1) stop 'HealpixExactLensedMap: GradPhi must be spin 1 field'
     if (A%npol ==1) then
      stop 'HealpixInterpLensedMap_GradPhi: currently only with pol'
     else
      call alm2LensedmapInterp(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU, factor)
     end if
  end subroutine  HealpixInterpLensedMap_GradPhi


  subroutine  HealpixQuadLensedMap(H,A, M, npix)
   ! ** Note does not give sky with accurate lensed C_l at l>~1200 **
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: GradPhi, M
     integer, intent(in) :: npix
     Type(HealpixAlm), intent(in) :: A
      
     call HealpixAlm2GradientMap(H,A,GradPhi,npix,'PHI')
     call HealpixQuadLensedMap_GradPhi(H,A, GradPhi,M)
     call HealpixMap_Free(GradPhi)

  end subroutine  HealpixQuadLensedMap


  subroutine  HealpixQuadLensedMap_GradPhi(H,A, GradPhi,M)
   ! ** Note does not give sky with accurate lensed C_l  at l>~1200 **
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: GradPhi, M
     Type(HealpixAlm), intent(in) :: A
      
     call HealpixMap_Init(M,GradPhi%npix,nmaps = A%npol)
     if (GradPhi%spin /=1) stop 'HealpixExactLensedMap: GradPhi must be spin 1 field'
     if (A%npol ==1) then
      call alm2LensedQuadContrib(H, A%lmax, A%TEB, GradPhi%SpinField, M%TQU(:,1))
     else
     stop 'not done yet'
     end if
  end subroutine  HealpixQuadLensedMap_GradPhi
   

  subroutine HealpixAlm2Map(H, A, M, npix, DoPhi, DoT)
     Type (HealpixInfo) :: H
     Type(HealpixMap) :: M
     Type(HealpixAlm), intent(in) :: A
     integer, intent(in) :: npix
     logical, intent(in), optional :: DoPhi, DoT
     logical Phi
     integer npol

     if (present(DoPhi)) then
      Phi = DoPhi .and. A%HasPhi
     else
      Phi = A%HasPhi
      end if

     if (present(DoT)) then
       if (DoT) then 
           npol = A%npol 
       else  
           npol = 0 
       end if
     else
       npol = A%npol
      end if
       
    call HealpixMap_Init(M,npix,npol, spinmap = A%spin, HasPhi = Phi)

    if (M%nmaps > 0) then
     if (npol>1) then
        call polalm2map(H,A%lmax, A%TEB,M%TQU)
      else
        call scalalm2map(H,A%lmax, A%TEB,M%TQU(:,1))
      end if    
    end if
    if (M%spin /= nospinmap) then
        call spinalm2map(H,A%lmax, A%SpinEB,M%SpinField, A%spin)
    end if
    if (Phi) then
        call scalalm2map(H,A%lmax, A%Phi,M%Phi)
    end if

  end subroutine HealpixAlm2Map

  subroutine HealpixMap_PolToSpin2Field(M)
      Type(HealpixMap) :: M
      integer status

       deallocate(M%SpinField,stat =status)
       ALLOCATE(M%SpinField(0:M%npix-1),stat = status)
       M%spin = 2
       M%SpinField = cmplx(M%TQU(:,2),M%TQU(:,3))

  end  subroutine HealpixMap_PolToSpin2Field

  subroutine HealpixMap_Spin2FieldTopol(M, delspin)
      Type(HealpixMap) :: M
      logical, intent(in), optional :: delspin
      integer status

       if (M%nmaps==0) then
         call HealpixMap_AllocateTQU(M,3) 
       end if
       M%TQU(:,2) = real(M%SpinField)
       M%TQU(:,3) = aimag(M%SpinField)
       if (present(delspin)) then
       if (delspin) then
          deallocate(M%SpinField)
           M%spin = nospinmap
        end if
       end if

  end  subroutine HealpixMap_Spin2FieldToPol


  subroutine HealpixMap_udgrade(M, Mout, nside_out, pessimistic)
     use udgrade_nr
     Type(HealpixMap) :: M, Mout
     integer, intent(in) :: nside_out
     logical, intent(in), optional :: pessimistic
     integer i
     logical pess
     

    if (present(pessimistic)) then 
      pess = pessimistic
    else
      pess = .false.
    end if

    call HealpixMap_ForceNest(M)
    call HealpixMap_Init(Mout,nside2npix(nside_out), M%nmaps, nested=.true.)
     do i=1, M%nmaps
       call sub_udgrade_nest(M%TQU(:,i),M%nside,MOut%TQU(:,i),nside_out, fmissval, pess)
    end do

  end subroutine HealpixMap_udgrade

  subroutine HealpixMap_HaarTransform(M,Mdegrade,Mdetail)
   implicit none
   Type(HealpixMap) :: M, Mdegrade, Mdetail
   integer i, j
   real :: bas1(4) = (/ -1, -1,  1,  1 /) /4.
  ! real :: bas2(4) = (/ -1,  1,  0,  0 /) /2.
  ! real :: bas3(4) = (/  0,  0, -1,  1 /) /2.
   real :: bas2(4) = (/ -1,  1,  1,  -1 /)  /4.
   real :: bas3(4) = (/ -1,  1, -1,   1 /) /4.
   if (M%npix <=12) stop 'Map too coarse to Haar transform'
   call HealpixMap_ForceNest(M)
   call HealpixMap_Init(Mdetail,M%npix/4, M%nmaps*3, nested = .true.)
   call HealpixMap_udgrade(M,Mdegrade,M%nside/2, pessimistic = .true.)
   
   do j=1, M%nmaps
    do i=0,Mdetail%npix -1
      if (any(M%TQU(i*4:i*4+3,j) == fmissval)) then
        Mdetail%TQU(i,(j-1)*4+1:(j-1)*4+3) = fmissval
      else
       Mdetail%TQU(i,(j-1)*4+1) = sum(M%TQU(i*4:i*4+3,j)*bas1)
       Mdetail%TQU(i,(j-1)*4+2) = sum(M%TQU(i*4:i*4+3,j)*bas2)
       Mdetail%TQU(i,(j-1)*4+3) = sum(M%TQU(i*4:i*4+3,j)*bas3)
      end if
    end do
   end do

  end subroutine HealpixMap_HaarTransform


  subroutine HealpixMap_HaarReconstruct(Mdegrade, Mdetail, M)
     Type(HealpixMap) :: M, Mdegrade, Mdetail
   integer i, j
   real :: bas1(4) = (/ -1, -1,  1,  1 /) 
!   real :: bas2(4) = (/ -1,  1,  0,  0 /) 
!   real :: bas3(4) = (/  0,  0, -1,  1 /) 
   real :: bas2(4) = (/ -1,  1,  1,  -1 /) 
   real :: bas3(4) = (/ -1,  1, -1,  1 /) 

   if (Mdegrade%npix /= Mdetail%npix) stop 'map size mismatch'
   if (mod(Mdetail%nmaps,3)/=0) stop 'detail not sets of 3 maps' 
   call HealpixMap_ForceNest(Mdegrade)
   call HealpixMap_ForceNest(Mdetail)
   call HealpixMap_Init(M,Mdegrade%npix*4, Mdegrade%nmaps/3, nested = .true.)
   call HealpixMap_udgrade(Mdegrade,M,Mdegrade%nside*2)
   
   do j=1, M%nmaps
    do i=0,Mdetail%npix -1
        if (any(M%TQU(i*4:i*4+3,j) == fmissval)) then
          M%TQU(i*4:i*4+3,j) = fmissval
        else
         M%TQU(i*4:i*4+3,j) = M%TQU(i*4:i*4+3,j) + Mdetail%TQU(i,(j-1)*4+1)*bas1 &
          + Mdetail%TQU(i,(j-1)*4+2)*bas2 + Mdetail%TQU(i,(j-1)*4+3)*bas3
        end if
    end do
   end do

  end subroutine HealpixMap_HaarReconstruct

  subroutine HaarComponents_Free(C) 
    Type(HaarComponents) :: C
    integer i, status

    if (associated(C%details)) then
     do i=1, C%order
      call HealpixMap_Free(C%details(i))
     end do
     deallocate(C%details)
     nullify(C%details)
    end if
    call HealpixMap_Free(C%degraded)
 
  end subroutine HaarComponents_Free

  subroutine HealpixMap_HaarComponents(M,C, order)
     Type(HealpixMap) :: M, D
     Type(HaarComponents) :: C
     integer, intent(in) :: order
     integer i
  
     C%order = order
     allocate(C%details(order))
     call HealpixMap_Assign(D,M)
     do i=1, order
       call HealpixMap_HaarTransform(D,C%degraded,C%details(i))     
       if (i/= order) call HealpixMap_Assign(D, C%degraded)
     end do
  end subroutine HealpixMap_HaarComponents

  
   subroutine HealpixMap_FromHaarComponents(C,M)
     Type(HealpixMap) :: M, D
     Type(HaarComponents) :: C
     integer i
  
     call HealpixMap_Assign(D,C%Degraded)
     do i=C%order,1,-1
       call HealpixMap_HaarReconstruct(D, C%details(i), M)
       if (i/=1) call HealpixMap_Assign(D, M)
     end do

  end subroutine HealpixMap_FromHaarComponents


 subroutine HealpixMap_HarrPowerMap(C, M, nside_power, nside_map)
!Map the power in the three details making up each pixel is nside_power averaged at resolution 
!nside_map

   Type(HaarComponents) :: C
   Type(HealpixMap) :: M
   integer, intent(in) :: nside_power, nside_map
   integer i, fac, j, pixcount
  

   if (C%degraded%nmaps /= 1) stop 'Only does power spectra for single map'
   if (nside_power < nside_map) stop 'Can only make map of power at lower resolution'

   call HealpixMap_Init(M, nside2npix(nside_map), nested = .true.)
   fac = nside2npix(nside_power) / M%npix
   j = 1
   do while (C%details(j)%nside > nside_power)
    j=j+1
    if (j> C%order) stop 'Haar coeffs not available'
   end do
   if (C%details(j)%nside /= nside_power) stop 'Haar coeffs not available'
   do i=0, M%npix-1
     pixcount =  count( C%details(j)%TQU(i*fac:(i+1)*fac - 1,:) /= fmissval) 
     if (pixcount /= 0) then
      M%TQU(i,1) = sum(C%details(j)%TQU(i*fac:(i+1)*fac - 1,:)**2, &
              mask =  (C%details(j)%TQU(i*fac:(i+1)*fac - 1,:) /= fmissval) ) / pixcount
     else
      M%TQU(i,1) = fmissval
     end if         
   end do      

  
 end subroutine HealpixMap_HarrPowerMap

  subroutine HealpixMap_HarrPowerSpec(C, P)
   Type(HaarComponents) :: C
   real(dp) P(C%order)
   integer i, pix

   if (C%degraded%nmaps /= 1) stop 'Only does power spectra for single map'
   do i=1, C%order
    P(i) = sum(C%details(i)%TQU(:,:)**2, mask = C%details(i)%TQU(:,:) /= fmissval) / &
           count(C%details(i)%TQU(:,:)/=fmissval)
   end do      
  
 end subroutine HealpixMap_HarrPowerSpec

end module HealpixObj








