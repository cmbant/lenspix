!Module of useful routines and definitions

 module Lists
  !Currently implements lists of strings and lists of arrays of reals
  implicit none
  type real_pointer
    real, dimension(:), pointer :: p
  end type real_pointer

  type double_pointer
    double precision, dimension(:), pointer :: p
  end type double_pointer

  type String_pointer
    character, dimension(:), pointer :: p
  end type String_pointer


  Type TList_RealArr
    integer Count
    integer Delta
    integer Capacity
    type(Real_Pointer), dimension(:), pointer :: Items
  end Type TList_RealArr

  Type TStringList
    integer Count
    integer Delta
    integer Capacity
    type(String_Pointer), dimension(:), pointer :: Items
  end Type TStringList

 contains
 
   subroutine TList_RealArr_Init(L)
    Type (TList_RealArr) :: L
    
     L%Count = 0
     L%Capacity = 0
     L%Delta = 64
     nullify(L%items)

   end subroutine TList_RealArr_Init

   subroutine TList_RealArr_Clear(L)
    Type (TList_RealArr) :: L
    integer i, status
     
     do i=L%Count,1,-1 
       deallocate (L%Items(i)%P, stat = status)
     end do
    deallocate (L%Items, stat = status)
    L%Count = 0
    L%Capacity = 0

   end subroutine TList_RealArr_Clear

    
   subroutine TList_RealArr_Add(L, P)
    Type (TList_RealArr) :: L
    real, intent(in) :: P(:)
    integer s
  
    if (L%Count == L%Capacity) call TList_RealArr_SetCapacity(L, L%Capacity + L%Delta)
    s = size(P)
    L%Count = L%Count + 1
    allocate(L%Items(L%Count)%P(s))
    L%Items(L%Count)%P = P

   end subroutine TList_RealArr_Add

   subroutine TList_RealArr_SetCapacity(L, C)
    Type (TList_RealArr) :: L
    integer C
    type(Real_Pointer), dimension(:), pointer :: TmpItems
    
    if (L%Count > 0) then
      if (C < L%Count) stop 'TList_RealArr_SetCapacity: smaller than Count'
      allocate(TmpItems(L%Count))
      TmpItems = L%Items(1:L%Count)
      deallocate(L%Items)
      allocate(L%Items(C))
      L%Items(1:L%Count) = TmpItems
      deallocate(TmpItems)
    else
     allocate(L%Items(C))
    end if  
    L%Capacity = C
   end subroutine TList_RealArr_SetCapacity

   subroutine TList_RealArr_Delete(L, i)
    Type (TList_RealArr) :: L
    integer, intent(in) :: i
    integer status
     
     deallocate(L%items(i)%P, stat = status)
     if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
     L%Count = L%Count -1
     
   end subroutine TList_RealArr_Delete



   subroutine TList_RealArr_Thin(L, i)
    Type (TList_RealArr) :: L
    integer, intent(in) :: i
    integer newCount
    type(Real_Pointer), dimension(:), pointer :: TmpItems
    
    if (L%Count > 1) then
      newCount = (L%Count-1)/i+1
      allocate(TmpItems(newCount))
      TmpItems = L%Items(1:i:L%Count)
      deallocate(L%Items)
      L%Capacity = newCount
      allocate(L%Items(L%Capacity))
      L%Items(1:i:L%Count) = TmpItems
      L%Count = newCount
      deallocate(TmpItems)
    end if    
   end subroutine TList_RealArr_Thin

   subroutine TList_RealArr_ConfidVal(L, ix, limfrac, ix1, ix2, Lower, Upper)
   !Taking the ix'th entry in each array to be a sample, value for which
   !limfrac of the items between ix1 and ix2 (inc) are above or below
   !e.g. if limfrac = 0.05 get two tail 90% confidence limits
     Type (TList_RealArr) :: L
     integer, intent(IN) :: ix
     real, intent(IN) :: limfrac
     real, intent(OUT), optional :: Lower, Upper
     integer, intent(IN), optional :: ix1,ix2
     integer b,t,samps
     real pos, d
     type(Real_Pointer), dimension(:), pointer :: SortItems
    
     b=1
     t=L%Count
     if (present(ix1)) b = ix1
     if (present(ix2)) t = ix2
     samps = t - b + 1
  
     allocate(SortItems(samps))
     SortItems = L%Items(b:t)
     call QuickSortArr_Real(SortItems, 1, samps, ix)
     if (present(Lower)) then
       pos = (samps-1)*limfrac + 1 
       b = max(int(pos),1)
       Lower = SortItems(b)%P(ix)
      if (b < samps .and. pos>b) then
       d = pos - b
       Lower = Lower*(1 - d) + d * SortItems(b+1)%P(ix) 
      end if
     end if
     if (present(Upper)) then
      pos = (samps-1)*(1.-limfrac) + 1
      b = max(int(pos),1)
      Upper = SortItems(b)%P(ix)
      if (b < samps .and. pos>b) then
       d = pos - b
       Upper = Upper*(1 - d) + d * SortItems(b+1)%P(ix) 
      end if
     end if
   
     deallocate(SortItems)

    end subroutine TList_RealArr_ConfidVal



   subroutine TStringList_Init(L)
    Type (TStringList) :: L
    
     L%Count = 0
     L%Capacity = 0
     L%Delta = 128

   end subroutine TStringList_Init

   subroutine TStringList_Clear(L)
    Type (TStringList) :: L
    integer i, status
     
     do i=L%Count,1,-1 
       deallocate (L%Items(i)%P, stat = status)
     end do
    deallocate (L%Items, stat = status)
    L%Count = 0
    L%Capacity = 0

   end subroutine TStringList_Clear

    
   subroutine TStringList_Add(L, P)
    Type (TStringList) :: L
    character(LEN=*), intent(in) :: P
    integer s
  
    if (L%Count == L%Capacity) call TStringList_SetCapacity(L, L%Capacity + L%Delta)
    s = len_trim(P)
    L%Count = L%Count + 1
    allocate(L%Items(L%Count)%P(s))
    L%Items(L%Count)%P = P(1:s)

   end subroutine TStringList_Add

   subroutine TStringList_SetCapacity(L, C)
    Type (TStringList) :: L
    integer C
    type(String_Pointer), dimension(:), pointer :: TmpItems
    
    if (L%Count > 0) then
      if (C < L%Count) stop 'TStringList_SetCapacity: smaller than Count'
      allocate(TmpItems(L%Count))
      TmpItems = L%Items(1:L%Count)
      deallocate(L%Items)
      allocate(L%Items(C))
      L%Items(1:L%Count) = TmpItems
      deallocate(TmpItems)
    else
     allocate(L%Items(C))
    end if  
    L%Capacity = C
  
   end subroutine TStringList_SetCapacity

   subroutine TStringList_Delete(L, i)
    Type (TStringList) :: L
    integer, intent(in) :: i
    integer status
     
     deallocate(L%items(i)%P, stat = status)
     if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
     L%Count = L%Count -1
     
   end subroutine TStringList_Delete

      recursive subroutine QuickSortArr_Real(Arr, Lin, R, index)
      !Sorts an array of pointers to arrays of reals by the value of the index'th entry
      integer, intent(in) :: Lin, R, index
      type(real_pointer), dimension(*) :: Arr
      integer I, J, L
      real P
      type(real_pointer) :: temp
  
      L = Lin
      do

      I = L
      J = R
      P = Arr((L + R)/2)%p(index)
   
      do
      do while (Arr(I)%p(index) <  P) 
         I = I + 1
      end do
    
      do while (Arr(J)%p(index) > P) 
         J = J - 1
      end do

      if (I <= J) then
     
       Temp%p => Arr(I)%p
       Arr(I)%p => Arr(J)%p
       Arr(J)%p => Temp%p
       I = I + 1
       J = J - 1
      end if
      if (I > J) exit
      
      end do
    if (L < J) call QuickSortArr_Real(Arr, L, J, index);
    L = I
    if (I >= R) exit
    end do

    end subroutine QuickSortArr_Real



    recursive subroutine QuickSortArr(Arr, Lin, R, index)
      !Sorts an array of pointers to arrays of reals by the value of the index'th entry
      integer, intent(in) :: Lin, R, index
      type(double_pointer), dimension(*) :: Arr
      integer I, J, L
      double precision P
      type(double_pointer) :: temp
  
      L = Lin
      do

      I = L
      J = R
      P = Arr((L + R)/2)%p(index)
   
      do
      do while (Arr(I)%p(index) <  P) 
         I = I + 1
      end do
    
      do while (Arr(J)%p(index) > P) 
         J = J - 1
      end do

      if (I <= J) then
     
       Temp%p => Arr(I)%p
       Arr(I)%p => Arr(J)%p
       Arr(J)%p => Temp%p
       I = I + 1
       J = J - 1
      end if
      if (I > J) exit
      
      end do
    if (L < J) call QuickSortArr(Arr, L, J, index);
    L = I
    if (I >= R) exit
    end do

    end subroutine QuickSortArr


 end module Lists

  module AMLutils
  use Lists
!DEC$ IF .false.        
#ifdef DECONLY
!DEC$ ENDIF    
   !Comment out if linking to LAPACK/MKL separetly 
   !CXML only has LAPACK 2.0
    include 'CXML_INCLUDE.F90'
!DEC$ IF .false.
#endif
!DEC$ ENDIF
  

!DEC$ IF .false.        
#ifdef NAGF95
        use F90_UNIX
#endif
!DEC$ ENDIF     

!DEC$ IF .false.
#ifndef NAGF95
!DEC$ ENDIF
        integer iargc
        external iargc
!DEC$ IF .false.
#endif
!DEC$ ENDIF

#ifdef MPI
    include "mpif.h"
#endif


  integer :: Feedback = 1


  double precision, parameter :: pi=3.14159265358979323846264338328d0, &
      twopi=2*pi, fourpi=4*pi
  double precision, parameter :: root2 = 1.41421356237309504880168872421d0, sqrt2 = root2
  double precision, parameter :: log2 = 0.693147180559945309417232121458d0

  real, parameter :: pi_r = 3.141592653, twopi_r = 2*pi_r, fourpi_r = twopi_r*2

  integer, parameter :: tmp_file_unit = 50
  logical :: flush_write = .true.
    !True means no data lost on crashes, but may make it slower

  contains


  function GetParamCount()
 
    GetParamCount = iargc() 

  end function GetParamCount

  function GetParam(i)

   character(LEN=120) GetParam
   integer, intent(in) :: i
 
   if (iargc() < i) then
     GetParam = ''
   else
    call getarg(i,GetParam)
   end if
  end function GetParam

  function concat(S1,S2,S3,S4,S5,S6)
   character(LEN=*), intent(in) :: S1, S2
   character(LEN=*), intent(in) , optional :: S3, S4, S5, S6
   character(LEN = 1000) concat

   concat = trim(S1) // S2
   if (present(S3)) then
    concat = trim(concat) // S3
     if (present(S4)) then
       concat = trim(concat) // S4
       if (present(S5)) then
         concat = trim(concat) // S5
           if (present(S6)) then
             concat = trim(concat) // S6
           end if
       end if
     end if
   end if

  end function concat


  subroutine WriteS(S)
   character(LEN=*), intent(in) :: S

    write (*,*) trim(S)
 
  end subroutine WriteS

  function numcat(S, num)
   character(LEN=*) S
   character(LEN=120) numcat, numstr
   integer num

   write (numstr, *) num
   numcat = trim(S) // trim(adjustl(numstr))
   !OK, so can probably do with with a format statement too... 
  end function numcat


  function IntToStr(I)
   integer , intent(in) :: I
   character(LEN=30) IntToStr

   write (IntToStr,*) i
   IntToStr = adjustl(IntToStr)

  end function IntToStr

  
  function RealToStr(R, figs)
   real, intent(in) :: R
   integer, intent(in), optional :: figs
   character(LEN=30) RealToStr

    if (abs(R)>=0.001) then
     write (RealToStr,'(f12.6)') R

   RealToStr = adjustl(RealToStr)
   if (present(figs)) then
    RealToStr = RealToStr(1:figs)
   else
    RealToStr = RealToStr(1:6)  
   end if

    else
     if (present(figs)) then
      write (RealToStr,'(E'//trim(numcat('(',figs))//'.2)') R
     else
      write (RealToStr,'(G9.2)') R
     end if
     RealToStr = adjustl(RealToStr)
    end if
        

  end function RealToStr
  
  function IndexOf(aval,arr, n)
     integer, intent(in) :: n, arr(n), aval
     integer i

     do i=1,n
        if (arr(i)==aval) then
          IndexOf= i
          return
        end if
     end do
    IndexOf = 0

  end function IndexOf

  function MaxIndex(arr, n)
     integer, intent(in) :: n
     real, intent(in) :: arr(n)
     integer locs(1:1), MaxIndex

     locs = maxloc(arr(1:n))
     MaxIndex = locs(1)

  end function MaxIndex

 
  function MinIndex(arr, n)
     integer, intent(in) :: n
     real, intent(in) :: arr(n)
     integer locs(1:1), MinIndex

     locs = minloc(arr(1:n))
     MinIndex = locs(1)

  end function MinIndex


  function ExtractFilePath(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=120) ExtractFilePath
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
       if (aname(i:i)=='/') then
          ExtractFilePath = aname(1:i)
          return
       end if
    end do
    ExtractFilePath = ''

  end function ExtractFilePath

 function ExtractFileName(aname)
    character(LEN=*), intent(IN) :: aname
    character(LEN=120) ExtractFileName
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
       if (aname(i:i)=='/') then
          ExtractFileName = aname(i+1:len)
          return
       end if
    end do
    ExtractFileName = aname

  end function ExtractFileName

  subroutine DeleteFile(aname)
    character(LEN=*), intent(IN) :: aname

     open(unit = tmp_file_unit, file = aname, err = 2)
     close(unit = tmp_file_unit, status = 'DELETE')
 2   return
    
  end subroutine DeleteFile

  subroutine FlushFile(aunit)
    integer, intent(IN) :: aunit

!DEC$ IF .false.        
#ifdef IBMXL
     call flush_(aunit)
#else
!DEC$ ENDIF     
     call flush(aunit)
!DEC$ IF .false.
#endif
!DEC$ ENDIF
    
    
  end subroutine FlushFile


  function FileExists(aname)
    character(LEN=*), intent(IN) :: aname
    logical FileExists
        
        inquire(file=aname, exist = FileExists)

  end function FileExists

 subroutine OpenFile(aname, aunit,mode)
   character(LEN=*), intent(IN) :: aname,mode
   integer, intent(in) :: aunit


   open(unit=aunit,file=aname,form=mode,status='old', err=500)
   return

500 write(*,*) 'File not found: '//trim(aname)
    stop

 end subroutine OpenFile
 
 
 subroutine OpenTxtFile(aname, aunit)
   character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: aunit

   call OpenFile(aname,aunit,'formatted')
 
 end subroutine OpenTxtFile

 subroutine CreateTxtFile(aname, aunit)
    character(LEN=*), intent(IN) :: aname
   integer, intent(in) :: aunit

   call CreateFile(aname,aunit,'formatted')

 end subroutine CreateTxtFile

 
 subroutine CreateFile(aname, aunit,mode)
    character(LEN=*), intent(IN) :: aname,mode
   integer, intent(in) :: aunit

     open(unit=aunit,file=aname,form=mode,status='replace', err=500)

   return

500 write(*,*) 'Error creating file '//trim(aname)
    stop

 end subroutine CreateFile


 function LastFileLine(aname)
   character(LEN=*), intent(IN) :: aname
   character(LEN = 5000) LastFileLine, InLine

 
   InLine = ''
   call OpenTxtFile(aname,tmp_file_unit)
   do
   read(tmp_file_unit,'(a)', end = 200) InLine
 
   end do
 
200  close(tmp_file_unit)

   LastFileLine = InLine
 
 end function LastFileLine

 

      subroutine spline_real(x,y,n,y2)

      integer, intent(in) :: n
      real, intent(in) :: x(n),y(n)
      real, intent(out) :: y2(n)
      integer i,k
      real p,qn,sig,un
      real, dimension(:), allocatable :: u

       
      allocate(u(1:n))
  
        y2(1)=0
        u(1)=0
        
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.0 
   
        y2(i)=(sig-1.0)/p
      
         u(i)=(6.0*((y(i+1)-y(i))/(x(i+ &
         1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
         u(i-1))/p
      end do
        qn=0.0
        un=0.0

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      deallocate(u)
  
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      end subroutine spline_real


      subroutine spline_double(x,y,n,y2)

      integer, intent(in) :: n
      double precision, intent(in) :: x(n),y(n)
      double precision, intent(out) :: y2(n)
      integer i,k
      double precision p,qn,sig,un
      double precision, dimension(:), allocatable :: u

       
      allocate(u(1:n))
  
        y2(1)=0
        u(1)=0
        
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2 
   
        y2(i)=(sig-1)/p
      
         u(i)=(6*((y(i+1)-y(i))/(x(i+ &
         1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
         u(i-1))/p
      end do
      qn=0
      un=0

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      deallocate(u)
  
!  (C) Copr. 1986-92 Numerical Recipes Software, adapted.
      end subroutine spline_double


      function LogGamma(x)
        real LogGamma
        real, intent(in) :: x
        integer i, j
        real r

        i = nint(x*2)
        if (abs(i-x*2) > 1e-4) stop 'LogGamma function for half integral only'
        if (mod(i,2) == 0) then
           r=0
           do j = 2, i/2-1
              r = r + log(real(j))
           end do
           LogGamma = r
        else
           r = log(pi)/2
           do j = 1, i-2 , 2
             r = r+ log(j/2.0)
           end do
           LogGamma = r
        end if

      end function LogGamma


    subroutine GetThreeJs(thrcof,l2in,l3in,m2in,m3in)
      !Recursive evaluation of 3j symbols. Does minimal error checking on input parameters.
      implicit none
      integer, parameter :: dl = KIND(1.d0)
      integer, intent(in) :: l2in,l3in, m2in,m3in
      real(dl), dimension(*) :: thrcof
      
      INTEGER, PARAMETER :: i8 = 8
 
      integer(i8) :: l2,l3,m2,m3
      integer(i8) :: l1, m1, l1min,l1max, lmatch, nfin, a1, a2
      
      real(dl) :: newfac, oldfac, sumfor, c1,c2,c1old, dv, denom, x, sum1, sumuni
      real(dl) :: x1,x2,x3, y,y1,y2,y3,sum2,sumbac, ratio,cnorm, sign1, thresh
      integer i,ier, index, nlim, sign2
      integer nfinp1,nfinp2,nfinp3, lstep, nstep2,n
      real(dl), parameter :: zero = 0._dl, one = 1._dl
      real(dl), parameter ::  tiny = 1.0d-10, srtiny=1.0d-5, huge = 1.d10, srhuge = 1.d5
  
    ! routine to generate set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3)

    ! by recursion from l1min = max(abs(l2-l3),abs(m1)) 
    !                to l1max = l2+l3
    ! the resulting 3j-coeffs are stored as thrcof(l1-l1min+1)

    ! to achieve the numerical stability, the recursion will proceed
    ! simultaneously forwards and backwards, starting from l1min and l1max
    ! respectively.
    !
    ! lmatch is the l1-value at which forward and backward recursion are matched.
    !
    ! ndim is the length of the array thrcof
    !
    ! ier = -1 for all 3j vanish(l2-abs(m2)<0, l3-abs(m3)<0 or not integer)
    ! ier = -2 if possible 3j's exceed ndim
    ! ier >= 0 otherwise

      l2=l2in
      l3=l3in
      m2=m2in
      m3=m3in
      newfac = 0
      lmatch = 0
      m1 = -(m2+m3)

    ! check relative magnitude of l and m values
      ier = 0
 
      if (l2 < abs(m2) .or. l3 < m3) then
      ier = -1
      stop 'error ier = -1'
      return
      end if

    ! limits for l1
      l1min = max(abs(l2-l3),abs(m1))
      l1max = l2+l3

      if (l1min >= l1max) then
       if (l1min/=l1max) then
       ier = -1
        stop 'error ier = -1'
       return
       end if

    ! reached if l1 can take only one value, i.e.l1min=l1max
      thrcof(1) = (-1)**abs(l2+m2-l3+m3)/sqrt(real(l1min+l2+l3+1,dl))
      return

      end if

      nfin = l1max-l1min+1
 
    ! starting forward recursion from l1min taking nstep1 steps
      l1 = l1min
      thrcof(1) = srtiny
      sum1 = (2*l1 + 1)*tiny

      lstep = 1

30    lstep = lstep+1
      l1 = l1+1

      oldfac = newfac
      a1 = (l1+l2+l3+1)*(l1-l2+l3)*(l1+l2-l3)
      a2 = (l1+m1)*(l1-m1)*(-l1+l2+l3+1)
      newfac = sqrt(a2*real(a1,dl))
      if (l1 == 1) then
         !IF L1 = 1  (L1-1) HAS TO BE FACTORED OUT OF DV, HENCE
         c1 = -(2*l1-1)*l1*(m3-m2)/newfac
      else

       dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 + l1*(l1-1)*(m3-m2)
       denom = (l1-1)*newfac

       if (lstep > 2) c1old = abs(c1)
       c1 = -(2*l1-1)*dv/denom

      end if

      if (lstep<= 2) then

    ! if l1=l1min+1 the third term in the recursion eqn vanishes, hence
       x = srtiny*c1
       thrcof(2) = x
       sum1 = sum1+tiny*(2*l1+1)*c1*c1
       if(lstep==nfin) then
          sumuni=sum1
          go to 230
       end if
       goto 30

      end if

      c2 = -l1*oldfac/denom

    ! recursion to the next 3j-coeff x  
      x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
      thrcof(lstep) = x
      sumfor = sum1
      sum1 = sum1 + (2*l1+1)*x*x
      if (lstep/=nfin) then

    ! see if last unnormalised 3j-coeff exceeds srhuge
      if (abs(x) >= srhuge) then
     
         ! REACHED IF LAST 3J-COEFFICIENT LARGER THAN SRHUGE
         ! SO THAT THE RECURSION SERIES THRCOF(1), ... , THRCOF(LSTEP)
         ! HAS TO BE RESCALED TO PREVENT OVERFLOW
     
         ier = ier+1
         do i = 1, lstep
            if (abs(thrcof(i)) < srtiny) thrcof(i)= zero
            thrcof(i) = thrcof(i)/srhuge
         end do

         sum1 = sum1/huge
         sumfor = sumfor/huge
         x = x/srhuge

      end if

    ! as long as abs(c1) is decreasing, the recursion proceeds towards increasing
    ! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
    ! detected, the recursion direction is reversed.

     if (c1old > abs(c1)) goto 30

     end if !lstep/=nfin

    ! keep three 3j-coeffs around lmatch for comparison with backward recursion

      lmatch = l1-1
      x1 = x
      x2 = thrcof(lstep-1)
      x3 = thrcof(lstep-2)
      nstep2 = nfin-lstep+3

    ! --------------------------------------------------------------------------
    !
    ! starting backward recursion from l1max taking nstep2 stpes, so that
    ! forward and backward recursion overlap at 3 points 
    ! l1 = lmatch-1, lmatch, lmatch+1

      nfinp1 = nfin+1
      nfinp2 = nfin+2
      nfinp3 = nfin+3
      l1 = l1max
      thrcof(nfin) = srtiny
      sum2 = tiny*(2*l1+1)
 
      l1 = l1+2
      lstep=1

      do
      lstep = lstep + 1
      l1= l1-1

      oldfac = newfac
      a1 = (l1+l2+l3)*(l1-l2+l3-1)*(l1+l2-l3-1)
      a2 = (l1+m1-1)*(l1-m1-1)*(-l1+l2+l3+2)
      newfac = sqrt(a1*real(a2,dl))

      dv = -l2*(l2+1)*m1 + l3*(l3+1)*m1 +l1*(l1-1)*(m3-m2)

      denom = l1*newfac
      c1 = -(2*l1-1)*dv/denom
      if (lstep <= 2) then

         ! if l2=l2max+1, the third term in the recursion vanishes
     
         y = srtiny*c1
         thrcof(nfin-1) = y
         sumbac = sum2
         sum2 = sum2 + tiny*(2*l1-3)*c1*c1

         cycle

      end if

      c2 = -(l1-1)*oldfac/denom

    ! recursion to the next 3j-coeff y
      y = c1*thrcof(nfinp2-lstep)+c2*thrcof(nfinp3-lstep)

      if (lstep==nstep2) exit
  
      thrcof(nfinp1-lstep) = y
      sumbac = sum2
      sum2 = sum2+(2*l1-3)*y*y

    ! see if last unnormalised 3j-coeff exceeds srhuge
      if (abs(y) >= srhuge) then
     
         ! reached if 3j-coeff larger than srhuge so that the recursion series
         ! thrcof(nfin),..., thrcof(nfin-lstep+1) has to be rescaled to prevent overflow
     
         ier=ier+1
         do i = 1, lstep
            index=nfin-i+1
            if (abs(thrcof(index)) < srtiny) thrcof(index)=zero
            thrcof(index) = thrcof(index)/srhuge
         end do

         sum2=sum2/huge
         sumbac=sumbac/huge

      end if

      end do

    ! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
    ! corresponding backward recursion vals y1, y2, y3

      y3 = y
      y2 = thrcof(nfinp2-lstep)
      y1 = thrcof(nfinp3-lstep)

    ! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal error

      ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
      nlim = nfin-nstep2+1

      if (abs(ratio) >= 1) then

       thrcof(1:nlim) = ratio*thrcof(1:nlim) 
       sumuni = ratio*ratio*sumfor + sumbac

      else

      nlim = nlim+1
      ratio = 1/ratio
      do n = nlim, nfin
         thrcof(n) = ratio*thrcof(n)
      end do
      sumuni = sumfor + ratio*ratio*sumbac

      end if
    ! normalise 3j-coeffs

230  cnorm = 1/sqrt(sumuni)

    ! sign convention for last 3j-coeff determines overall phase

      sign1 = sign(one,thrcof(nfin))
      sign2 = (-1)**(abs(l2+m2-l3+m3))
      if (sign1*sign2 <= 0) then
        cnorm = -cnorm
      end if
      if (abs(cnorm) >= one) then
         thrcof(1:nfin) = cnorm*thrcof(1:nfin)
         return
      end if

      thresh = tiny/abs(cnorm)

      do n = 1, nfin
         if (abs(thrcof(n)) < thresh) thrcof(n) = zero
         thrcof(n) = cnorm*thrcof(n)
      end do
      return 

    end subroutine GetThreeJs



  end module AMLutils
 
  

module Random
 double precision  gset
 integer iset
 integer :: rand_inst = 0 

contains
   
  subroutine initRandom(i)
  use AMLUtils
  implicit none
  external rmarin
  integer, optional, intent(IN) :: i
  integer kl,ij
  character(len=10) :: fred
  real :: klr
      if (present(i)) then
       iset=0
       kl = 9373
       ij = i
      else
       call system_clock(count=ij)
       ij = mod(ij + rand_inst*100, 31328)
       call date_and_time(time=fred)
       read (fred,'(e10.3)') klr
       kl = mod(int(klr*1000), 30081)       
      end if

      if (Feedback > 0 ) write(*,'(" Random seeds:",1I6,",",1I6," rand_inst:",1I4)")') ij,kl,rand_inst
      call rmarin(ij,kl)
  end subroutine initRandom

  subroutine RandIndices(indices, nmax, n)
    integer, intent(in) :: nmax, n
    integer indices(n),i, ix
    integer tmp(nmax)
    double precision ranmar
    external ranmar

    if (n> nmax) stop 'Error in RandIndices, n > nmax'
    do i=1, nmax
       tmp(i)=i
    end do
    do i=1, n
       ix = int(ranmar()*(nmax +1 -i)) + 1
       indices(i) = tmp(ix)
       tmp(ix) = tmp(nmax+1-i)
    end do

  end subroutine RandIndices

  double precision function GAUSSIAN1()
    implicit none
    double precision R, V1, V2, FAC
    double precision ranmar
    external ranmar
 

        if (ISET==0) then
        R=2
        do while (R >= 1.d0)
        V1=2.d0*ranmar()-1.d0
        V2=2.d0*ranmar()-1.d0
        R=V1**2+V2**2
        end do
        FAC=sqrt(-2.d0*log(R)/R)
        GSET=V1*FAC
        GAUSSIAN1=V2*FAC
        ISET=1
      else
        GAUSSIAN1=GSET
        ISET=0
      endif
   
      end function GAUSSIAN1




end module Random


! This random number generator originally appeared in ''Toward a Universal 
! Random Number Generator'' by George Marsaglia and Arif Zaman. 
! Florida State University Report: FSU-SCRI-87-50 (1987)
! 
! It was later modified by F. James and published in ''A Review of Pseudo-
! random Number Generators'' 
! 
! THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
!    (However, a newly discovered technique can yield 
!        a period of 10^600. But that is still in the development stage.)
!
! It passes ALL of the tests for random number generators and has a period 
!   of 2^144, is completely portable (gives bit identical results on all 
!   machines with at least 24-bit mantissas in the floating point 
!   representation). 
! 
! The algorithm is a combination of a Fibonacci sequence (with lags of 97
!   and 33, and operation "subtraction plus one, modulo one") and an 
!   "arithmetic sequence" (using subtraction).
!
! On a Vax 11/780, this random number generator can produce a number in 
!    13 microseconds. 
!======================================================================== 
!
!      PROGRAM TstRAN
!     INTEGER IJ, KL, I
! Thee are the seeds needed to produce the test case results
!      IJ = 1802
!      KL = 9373
!
!
! Do the initialization
!      call rmarin(ij,kl)
!
! Generate 20000 random numbers
!      do 10 I = 1, 20000
!         x = RANMAR()
!10    continue
!
! If the random number generator is working properly, the next six random
!    numbers should be:
!          6533892.0  14220222.0  7275067.0
!    6172232.0  8354498.0   10633180.0
!           
!           
!        
!      write(6,20) (4096.0*4096.0*RANMAR(), I=1,6)
!20    format (3f12.1)
!      end
!
      subroutine RMARIN(IJ,KL)
! This is the initialization routine for the random number generator RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
!The random number sequences created by these two seeds are of sufficient 
! length to complete an entire calculation with. For example, if sveral 
! different groups are working on different parts of the same calculation,
! each group could be assigned its own IJ seed. This would leave each group
! with 30000 choices for the second seed. That is to say, this random 
! number generator can create 900 million different subsequences -- with 
! each subsequence having a length of approximately 10^30.
!
! Use IJ = 1802 & KL = 9373 to test the random number generator. The
! subroutine RANMAR should be used to generate 20000 random numbers.
! Then display the next six random numbers generated multiplied by 4096*4096
! If the random number generator is working properly, the random numbers
!    should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0
      double precision U(97), C, CD, CM
      integer I97, J97
    
!      INTEGER IRM(103)
      
      common /RASET1/ U, C, CD, CM, I97, J97
      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or. &
         KL .lt. 0  .or.  KL .gt. 30081 ) then
          print '(A)', ' The first random number seed must have a value  between 0 and 31328'
          print '(A)',' The second seed must have a value between 0 and   30081'
            stop
      endif
      I = mod(IJ/177, 177) + 2
      J = mod(IJ    , 177) + 2
      K = mod(KL/169, 178) + 1
      L = mod(KL,     169) 
      do 2 II = 1, 97
         S = 0.0
         T = 0.5
         do 3 JJ = 1, 24
            M = mod(mod(I*J, 179)*K, 179)
            I = J
            J = K
            K = M
            L = mod(53*L+1, 169)
            if (mod(L*M, 64) .ge. 32) then
               S = S + T
            endif
            T = 0.5 * T
3        continue
         U(II) = S
2     continue
      C = 362436.0 / 16777216.0
      CD = 7654321.0 / 16777216.0
      CM = 16777213.0 /16777216.0
      I97 = 97
      J97 = 33
    
      return
      end

      double precision function RANMAR()
! This is the random number generator proposed by George Marsaglia in 
! Florida State University Report: FSU-SCRI-87-50
! It was slightly modified by F. James to produce an array of pseudorandom
! numbers.
      double precision U(97), C, CD, CM
      integer I97, J97
    
      common /RASET1/ U, C, CD, CM, I97, J97
!      INTEGER IVEC
         UNI = U(I97) - U(J97)
         if( UNI .lt. 0.0 ) UNI = UNI + 1.0
         U(I97) = UNI
         I97 = I97 - 1
         if(I97 .eq. 0) I97 = 97
         J97 = J97 - 1
         if(J97 .eq. 0) J97 = 97
         C = C - CD
         if( C .lt. 0.0 ) C = C + CM
         UNI = UNI - C
         if( UNI .lt. 0.0 ) UNI = UNI + 1.0 ! bug?
         RANMAR = UNI
      return
      end



