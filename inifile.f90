!Module to read in name/value pairs from a file, with each line of the form line 'name = value'
!Should correctly interpret FITS headers
!Antony Lewis (http://cosmologist.info/). Released to the public domain.
!This version October 2002.

module IniFile
 implicit none
 public
  character (LEN=80), dimension(:), allocatable :: Keys, Vals
  integer :: NumKeys = 0, ArrSize = 0
  logical :: Ini_fail_on_not_found = .false., SlashComments = .false. 
 
contains
  
  subroutine realloc(NewSize)
     integer, intent(IN) :: NewSize
     character (LEN=80), dimension(:), allocatable :: TMP

     if (NewSize==0) then
        if (allocated(Keys)) deallocate(Keys)
        if (allocated(Vals)) deallocate(Vals)
     else
      if (allocated(Keys)) then
      allocate(TMP(1:NewSize))
      TMP=''
      TMP(1:min(NewSize,size(Keys, DIM=1)))=Keys(:)
      deallocate(Keys)
      allocate(Keys(1:NewSize))
      Keys=TMP
      deallocate(TMP)
      else
         allocate(Keys(1:NewSize))
         Keys=''
      end if
      if (allocated(Vals)) then
      allocate(TMP(1:NewSize))
      TMP=''
      TMP(1:min(NewSize,size(Vals, DIM=1)))=Vals(:)
      deallocate(Vals)
      allocate(Vals(1:NewSize))
      Vals=TMP
      deallocate(TMP)
      else
         allocate(Vals(1:NewSize))
         Vals=''
      end if
      
     end if

     ArrSize=NewSize

  end subroutine realloc

  subroutine Ini_AddLine(AInLine)
    character (LEN=*), intent(IN) :: AInLine
    integer EqPos, slashpos, lastpos
    character (LEN=120) :: S, InLine

      InLine=trim(adjustl(AInLine))
      EqPos = scan(InLine,'=')
      if (EqPos/=0 .and. InLine(1:1)/='#' .and. InLine(1:7) /= 'COMMENT' ) then
   
         NumKeys=NumKeys+1      
  
         if (NumKeys > ArrSize) then
            ArrSize=ArrSize+100
            call realloc(ArrSize)
         end if

         Keys(NumKeys) = trim(InLine(1:EqPos-1))
         S = adjustl(InLine(EqPos+1:)) 
           if (SlashComments) then
           slashpos=scan(S,'/')
           if (slashpos /= 0) then
              S  = S(1:slashpos-1)
           end if
         end if
         lastpos=len_trim(S)
         if (S(1:1)=='''' .and. S(lastpos:lastpos)=='''') then
           S = S(2:lastpos-1)
         end if
         Vals(NumKeys)=trim(S)
   
      end if

  end subroutine Ini_AddLine

  subroutine Ini_Open(filename, unit_id,  error, slash_comments)
     character (LEN=*), intent(IN) :: filename
     integer, intent(IN) :: unit_id
     logical, intent(OUT) :: error
     logical, optional, intent(IN) :: slash_comments
     character (LEN=120) :: InLine
    
    if (present(slash_comments)) then
     SlashComments = slash_comments
    else
     SlashComments = .false.
    end if
     NumKeys=0

    
    open(unit=unit_id,file=filename,form='formatted',status='old', err=500)
   
    ArrSize=50
 
    call realloc(ArrSize)
    do 
      read (unit_id,'(a)',end=400) InLine
      if (InLine == 'END') exit;
      if (InLine /= '') call Ini_AddLine(InLine) 
    end do

400 close(unit_id)
    error=.false.
    return

500 error=.true.

  end subroutine Ini_Open

  subroutine Ini_Open_Fromlines(Lines, NumLines, slash_comments)
    integer, intent(IN) :: NumLines
    character (LEN=*), dimension(NumLines), intent(IN) :: Lines
    logical, intent(IN) :: slash_comments
    integer i

    SlashComments = slash_comments
    call realloc(NumLines)
    NumKeys=0
    do i=1,NumLines
       call Ini_AddLine(Lines(i))
    end do
  

  end  subroutine Ini_Open_Fromlines

  subroutine Ini_Close
    
    call realloc(0)

  end  subroutine Ini_Close
  
  function Ini_Read_String(Key, NotFoundFail)
   character (LEN=80) Ini_Read_String
   character (LEN=*), intent(IN) :: Key
   logical, optional, intent(IN) :: NotFoundFail
   integer i

   do i=1, NumKeys
      if (Key == Keys(i)) then
         Ini_Read_String = Vals(i)
         return
      end if
   end do
   Ini_Read_String=''
   if (Ini_fail_on_not_found) then
      write(*,*) 'key not found : '//Key
      stop
   end if
   if (present(NotFoundFail)) then
      if (NotFoundFail) then
         write(*,*) 'key not found : '//Key
         stop
      end if
   end if

  end function Ini_Read_String


  function Ini_Read_Int(Key, Default)
   integer Ini_Read_Int
   integer, optional, intent(IN) :: Default
   character  (LEN=*), intent(IN) :: Key
  
   character(LEN=80) :: S
   
   S = Ini_Read_String(Key,.not. present(Default))
   if (S == '') then
      if (.not. present(Default)) then
        write(*,*) 'no value for key: '//Key
        stop
      end if
      Ini_Read_Int = Default
   else
   read (S,*, err = 10) Ini_Read_Int
   end if

  return

10 write (*,*) 'error reading integer for key: '//Key
   stop

  end function Ini_Read_Int

  function Ini_Read_Double(Key, Default)
   double precision Ini_Read_Double 
   double precision, optional, intent(IN) :: Default
   character (LEN=*), intent(IN) :: Key
   character(LEN=80) :: S
   
   S = Ini_Read_String(Key,.not. present(Default))
   if (S == '') then
      if (.not. present(Default)) then
        write(*,*) 'no value for key: '//Key
        stop
      end if
      Ini_Read_Double = Default
   else
   read (S,*, err=10) Ini_Read_Double
   end if

  return

10 write (*,*) 'error reading double for key: '//Key
   stop

  end function Ini_Read_Double

    function Ini_Read_Real(Key, Default)
    real Ini_Read_Real 
    real, optional, intent(IN) :: Default
    character (LEN=*), intent(IN) :: Key
    character(LEN=80) :: S
   
   S = Ini_Read_String(Key,.not. present(Default))
   if (S == '') then
      if (.not. present(Default)) then
        write(*,*) 'no value for key: '//Key
        stop
      end if
      Ini_Read_Real = Default
   else
   read (S,*, err=10) Ini_Read_Real
   end if

  return

10 write (*,*) 'error reading double for key: '//Key
   stop

  end function Ini_Read_Real

  function Ini_Read_Logical(Key, Default)
   logical Ini_Read_Logical
   logical, optional, intent(IN) :: Default
   character  (LEN=*), intent(IN) :: Key
  
   character(LEN=80) :: S
   
   S = Ini_Read_String(Key,.not. present(Default))
   if (S == '') then
      if (.not. present(Default)) then
        write(*,*) 'no value for key: '//Key
        stop
      end if
      Ini_Read_Logical = Default
   else
   read (S,*, err = 10) Ini_Read_Logical
   end if

  return

10 write (*,*) 'error reading logical for key: '//Key
   stop
  end function Ini_Read_Logical

end module IniFile

