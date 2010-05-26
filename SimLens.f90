!Simple program demonstrating how to generate a simulated lensed map
!AL, Feb 2004
program SimLensCMB
 use HealpixObj
 use HealpixVis
 use Random
 use spinalm_tools
 use IniFile
 use AMLUtils
 implicit none
 Type(HealpixInfo)  :: H
 Type(HealpixMap)   :: M, GradPhi
 Type(HealpixPower) :: P
 Type(HealpixAlm)   :: A
 
 integer            :: nside, lmax, npix
 character(LEN=80)  :: w8name = '../healpix-1.2/data/'
 character(LEN=80)  :: cls_file, cls_lensed_file
 integer, parameter :: lens_interp =1, lens_exact = 2
 integer :: lens_method = lens_interp
 integer :: i, interp_factor
 logical :: err
#ifdef MPIPIX
 call mpi_init(i)
#endif

 Ini_Fail_On_Not_Found = .true.
 call Ini_Open(GetParam(1), 3,err)
 if (err) then
#ifdef MPIPIX
    call mpi_finalize(i)
#endif
   stop 'No ini'
 end if
 w8name = Ini_Read_String('w8dir')
 nside  = Ini_Read_Int('nside')
 npix = 12*nside**2

 lmax   = Ini_Read_Int('lmax')  
 cls_file = Ini_Read_String('cls_file')
 cls_lensed_file = Ini_Read_String('cls_lensed_file')
 lens_method = Ini_Read_Int('lens_method')

 Ini_Fail_On_Not_Found = .false.
 if (lens_method == lens_interp) interp_factor = Ini_Read_Int('interp_factor',8)

 call Ini_Close

 if (nside > 1024) then
  call HealpixInit(H,nside, lmax,.true.) 
  !the ring weight files are currently only supplied for nside <= 1024: set to unity
 else
  call HealpixInit(H,nside, lmax,.true., w8dir=w8name) 
 end if

 if (H%MpiID ==0) then !if we are main thread
  !All but main thread stay in HealpixInit

  call HealpixPower_ReadFromTextFile(P,cls_file,lmax,pol=.true.,dolens = .true.)
  !Reads in unlensed C_l text files as produced by CAMB (or CMBFAST if you aren't doing lensing)

  call HealpixAlm_Sim(A, P, HasPhi=.true., dopol = .true.)

  call HealpixAlm2GradientMap(H,A, GradPhi,npix,'PHI')

  if (lens_method == lens_exact) then
   call HealpixExactLensedMap_GradPhi(H,A,GradPhi,M)
  else if (lens_method == lens_interp) then
   call HealpixInterpLensedMap_GradPhi(H,A,GradPhi, M, interp_factor)
  else
   stop 'unknown lens_method'
  end if

  call HealpixMap2Alm(H,M, A, lmax, dopol = .true.)
  !This automatically frees previous content of A, and returns new one

  call HealpixAlm2Power(A,P)
  call HealpixAlm_Free(A)
  !Note usually no need to free objects unless memory is short

  call HealpixPower_Write(P,cls_lensed_file)

  !Save map to .fits file
  !call HealpixMap_Write(M, '!lensed_map.fits')
 
  end if

#ifdef MPIPIX
    call HealpixFree(H)
    call mpi_finalize(i)
#endif


end program SimLensCMB
