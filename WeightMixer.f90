!Code for simple simultation of hybrid Pseudo-Cl estiamtors
!Assumes for simicity that polarization and temperature noise
!and mask are the same (noise proportional with ENoiseFac)
module WeightMixing
 use PseudoCl
 use HealpixObj
 use HealpixVis
 use Random
 use spinalm_tools
 use IniFile
 use AMLUtils
 use PseudoCl
 use MatrixUtils
 implicit none
 character(LEN=*), parameter :: nulltest='nulltest'
 character(LEN=*), parameter :: pixtest='pixtest'

 character(LEN=256) :: data_dir = 'data/'
 character(LEN=256) :: w8name = '../Healpix_2.00/data/'
 integer, parameter :: max_detectors = 4 
 integer:: nside, healpix_res, lmin,lmax, npix
 logical :: want_pol, est_noise
 real(dp) :: apodize_mask_fwhm
 logical :: sim_noise = .true., sim_signal = .true.
 logical :: fullsky_test = .false.
 logical :: uniform_noise = .false.
 logical :: cross_spectra = .false.
 logical :: no_cache = .false.
 logical :: do_lowl_separation = .true. !project lowl map out of high-l one
 real(dp) :: ENoiseFac 
 real(dp) :: white_NL, white_NL_P !N_l when testing with white noise
 integer :: l_exact = 20
 integer :: l_exact_margin = 20 !50
 real(dp)  :: WellSupported = 0.99_dp
 real(dp) :: HighlNoiseBoost = 1.0_dp 
 real(dp) :: fake_noise = 0._dp
 character(LEN=3) :: map_unit = 'muK'
 real(dp) :: map_scale = 1._dp !mutiply map to get muK
 logical :: get_mean_likes = .false.

 real(dp) :: noise_inv_fwhm = 0.d0
 real :: point_source_A = 0.d0
 real :: point_source_A_frac_error = 0.d0

 real(sp) :: noise
 character(LEN=256) :: file_stem
 character(LEN=256) :: beam_filename_format, noise_filename_format, combined_filename_format, year_filename_format
 character(LEN=256) :: fits_mask !0 and 1 foreground mask
 character(LEN=256) :: processing_mask !optional mask says where hit counts are zero and must be zeroed
 character(LEN=256) :: inv_noise_map !Map to smooth to get 'inverse noise' weight map if noise_inv_fwhm /=0

!Testing
 character(LEN=256) :: check_cls_file1, check_cls_file2

 integer vec_size
 integer nchannels, nweights, nyears, InYears
 Type(TChannel), target, dimension(:), allocatable :: Channels
 Type(HealPixMap), dimension(:), target, allocatable :: WeightMaps
 Type (TCouplingMatrix), dimension(:), pointer :: Coupler, XiMatrices
 Type(TCovMatSet), allocatable :: HybridMix(:)
 Type(HealpixMap) :: SmoothedNoise 
contains


 subroutine AnalyseMap(H, Maps, HybridP)
  Type(HealpixInfo) :: H
  Type(HealpixMap) :: Maps(:)
  Type(HealpixPower) :: HybridP, CHat
  Type(HealpixCrossPowers) :: CovPowers, MaskPowers, PCls
 
  Type(HealPixMap):: WMaps(nweights*nchannels)
  integer chan1, chan2, channel,nmaps
  integer Pix, ix, x,y
  integer i,j
  real(dp) StTime
      
    nmaps = nweights*nchannels
  
    ix = 0
    do channel = 1, nchannels
    do i=1, nweights
          ix =ix + 1
          call HealpixMap_Nullify(WMaps(ix))
          call HealpixMapMulCut(Maps(Channel),WeightMaps(i),WMaps(ix), 1)
            !   call HealpixVis_Map2ppmfile(WeightMaps(i), 'outfiles/weightmap.ppm')
    end do
    end do
                
    call HealpixMapSet2CrossPowers(H, WMaps, PCls, nmaps, lmax)
    call HealpixMapArray_Free(WMaps)
                
    StTime = GeteTime()  
    call HealpixPower_Init(HybridP,lmax, want_pol)
    call HealpixPower_Nullify(CHat)
            i=0
            Pix = 0
            do chan1 = 1, nchannels
            do x=1, nweights
             i=i+1
             j=0
            do chan2 = 1, nchannels
             do y=1, nweights
                j=j+1
                if (j>i) cycle
                Pix = Pix + 1
                ix = sym_ix(nweights, x,y)
                print *, chan1, x,chan2,y,ix
                
                call PseudoCl_GetCHat(Coupler(ix),PCls%Ps(i,j), Chat)
                call TBeam_PowerSmooth2(Channels(chan1)%Beam,Channels(chan2)%Beam,CHat,+1)
                HybridP%Cl(lmin:lmax,C_T)= HybridP%Cl(lmin:lmax,C_T) + &
                    MatMul(HybridMix(1)%Cov(Pix)%C, CHat%Cl(lmin:lmax,C_T))
                if (vec_size >=3) then
                 HybridP%Cl(lmin:lmax,C_C)= HybridP%Cl(lmin:lmax,C_C) + &
                    MatMul(HybridMix(2)%Cov(Pix)%C, CHat%Cl(lmin:lmax,C_C))
                 HybridP%Cl(lmin:lmax,C_E)= HybridP%Cl(lmin:lmax,C_E) + &
                    MatMul(HybridMix(3)%Cov(Pix)%C, CHat%Cl(lmin:lmax,C_E))
                 if (vec_size >=4) then
                 HybridP%Cl(lmin:lmax,C_B)= HybridP%Cl(lmin:lmax,C_B) + &
                     MatMul(HybridMix(4)%Cov(Pix)%C, CHat%Cl(lmin:lmax,C_B))
                 end if
                end if
            end do
            end do 
            end do
            end do
                
            call HealpixCrossPowers_Free(PCls)
   end subroutine AnalyseMap


 function FormatFilename(FString, Channel, Detector, Year) result (formatted)
   character(LEN=*), intent(in) :: FString
   character(LEN=*), intent(in), optional ::  Channel
   integer, intent(in), optional :: Detector, Year
   character(LEN=1024) formatted
   
   formatted = FString
   
   call StringReplace('%RES%',IntToStr(healpix_res), formatted)
   if (present(Channel)) call StringReplace('%CHANNEL%',Channel,formatted)
   if (present(Detector))  call StringReplace('%DA%',IntToStr(Detector),formatted)
   if (present(Year))  call StringReplace('%YEAR%',IntToStr(Year),formatted)
   
 end function FormatFilename

 function CacheName(FName, nostem) result (cache)
  character(LEN=*) :: FName
  character(LEN=1024) :: cache
  logical nostem
  
  cache = data_dir
  if (.not. nostem) cache = trim(cache)//trim(ExtractFileName(file_stem))
  cache = trim(cache)//trim(ExtractFileName(Fname))
  if (ExtractFileExt(cache)/= '.fits') cache = trim(cache)//'.fits'
  
 end function CacheName

 subroutine ReadBeams(C, beamfile, lmax)
 !Just weight uniformly for now for effective combined map
  Type(TChannel) :: C
  character(LEN=*), intent(in) :: beamfile
  character(LEN=256) :: file
  integer, intent(in) :: lmax
  integer i
  Type(TBeam) B
  
  allocate(C%DetectorBeams(C%Count))
  
  do i=1, C%Count
   file = FormatFilename(beamfile,C%Name, i)
   C%DetectorBeams(i)%beam_transfer = .true.
   print *, 'reading beam file: '//trim(file)
   call TBeam_ReadFile(C%DetectorBeams(i),file,lmax)
   if (i==1) then
    call TBeam_ReadFile(C%Beam,file,lmax)
   else
    C%Beam%Beam= C%Beam%Beam + C%DetectorBeams(i)%Beam
   end if
  end do
  C%Beam%Beam = C%Beam%Beam/C%Count
 
 end  subroutine ReadBeams
 
 subroutine LoadPixelWindow(pixlw,nside)
  use alm_tools
  integer nside
  real(dp):: pixlw(0:,1:)
  character(LEN=5) :: sstr

   write(sstr,'(I4.4)') nside
   call pixel_window(pixlw, windowfile=trim(w8name)//"pixel_window_n"//trim(sstr)//".fits")
 
 end subroutine LoadPixelWindow
 
 subroutine PixSmoothBeams(lmax)
  real(dp), allocatable :: pixlw(:,:)
  integer npw, i, lmax, channel
  
   npw = 4*nside + 1 
   if (lmax > npw-1) call MpiStop('lmax too large for pixel window')
   allocate(pixlw(0:npw-1,1:3))
   call LoadPixelWindow(pixlw, nside)
   do channel=1, nchannels
    do i=1, Channels(channel)%Count
     Channels(channel)%DetectorBeams(i)%Beam= Channels(channel)%DetectorBeams(i)%Beam * pixlw(0:lmax,1)
    end do
    Channels(channel)%Beam%Beam= Channels(channel)%Beam%Beam * pixlw(0:lmax,1)
   end do
   deallocate(pixlw)

 end subroutine PixSmoothBeams

 function r_nu(fGhz)
 real(dp):: r_nu, fGhz
 real(dp) :: x
 real(dp), parameter :: const =  1.76054182e-11  !in seconds
 !For WMAP5 point source
 !x = h nu/ k T_CMB
  
  x = const*(fGhz*1d9)
  r_nu = (exp(x)-1)**2/x**2/exp(x)
  
 end function r_nu

 subroutine AzimuthalMap(H,fname)
  Type(HealpixInfo) :: H
  Type(HealpixMap) :: Obs, Hits
  character(LEN=*) :: fname
  Type(HealpixAlm) :: Alm
  integer ith , nph,l
  
       call HealpixMap_Read(Obs,fname)
       call HealpixMap_Init(Hits, Obs%npix,nested=Obs%ordering==ord_nest,pol = .false.)
!Note using index 2
       Hits%TQU(:,1) =  1786.0**2/Obs%TQU(:,2)
       call HealpixMap_ForceRing(Hits)
       
       if (.false.) then
       do ith = H%ith_start(H%MpiId), H%ith_end(H%MpiId)   
      
       if (ith < H%nside) then  ! polar cap (north)
          nph = 4*ith
       else                   ! tropical band (north) + equator
          nph = 4*H%nside
       endif
       
         Hits%TQU(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,1) = &
            sum(dble(Hits%TQU(H%istart_north(ith-1):H%istart_north(ith-1)+nph-1,1)))/nph
        if (ith  <  2*H%nside) then
           Hits%TQU(H%istart_south(ith):H%istart_south(ith)+nph-1,1) = &
               sum(Hits%TQU(H%istart_south(ith):H%istart_south(ith)+nph-1,1))/nph
        end if
        
       end do     
       end if
       
      call HealpixMap2Alm(H, Hits, Alm, H%lmax)   
      call HealpixAlm2Map(H, Alm, Hits, 12*H%nside**2)   
      print *,'min/max', minval(Hits%TQU),maxval(Hits%TQU)
      call HealpixVis_Map2ppmfile(Hits, 'outfiles/azimmap.ppm')
           
      do l=0, H%lmax
        write(*,*) l, Alm%TEB(1,l,0)
      end do  
      
      call MpiStop()

end subroutine AzimuthalMap

 subroutine ReadYearMaps(fname, M, DA)
  character(LEN=*), intent(in):: fname
  character(LEN=256) :: aname
  Type(HealpixMap) :: M(nyears)
  integer year, DA
 
   do year = 1, nyears 
    aname = FormatFilename(fname, '', DA, year)
    call HealpixMap_Nullify(M(year))
    call HealPixMap_read(M(year), aname)
   end do
 
 end subroutine ReadYearMaps

 subroutine GetWeightedDataMaps(WMaps)
 Type(HealpixMap) :: WMaps(:)
 integer detector, ix, weight, year, channel
 character(LEN=256) :: aname
 Type(HealpixMap) :: M

 print *,'Reading maps to get weighted map array'
  ix = 0
  do channel = 1, nchannels  
   print *, 'channel '//trim(Channels(channel)%Name)

   do Detector=1, Channels(channel)%Count
     print *, 'Detector', Detector

       do year = 1, nyears 
       
                aname = FormatFilename(year_filename_format, Channels(Channel)%Name, Detector, year)
                call HealpixMap_Nullify(M)
                call HealPixMap_read(M, aname)
                call HealpixMap_SetToIndexOnly(M,1)
                M%TQU = M%TQU*map_scale/mK
                
                do weight = 1, nweights
                 ix = ix + 1
                 call HealpixMap_Nullify(WMaps(ix))
                 call HealpixMapMulCut(M,WeightMaps(weight),WMaps(ix), 1)
                end do
                call HealpixMap_Free(M)
                
       end do !year
   end do !detector
  end do !channel

 end subroutine GetWeightedDataMaps

 subroutine GetWeightedDataNulltestMaps(WMaps)
 !Get all linear combinations of years that should have zero signal
 !nyears maps, make from nyears+1 maps (i.e. actual number of years is nyears+1)
 Type(HealpixMap) :: WMaps(:), M, YearMaps(nyears+1)
 integer detector, i, ix, weight, year, channel
 character(LEN=256) :: aname

 print *,'Making signal-null maps and weighted array'
  ix = 0
  do channel = 1, nchannels  
   print *, 'channel '//trim(Channels(channel)%Name)

   do Detector=1, Channels(channel)%Count
     print *, 'Detector', Detector

       do year = 1, nyears+1 
       
                aname = FormatFilename(year_filename_format, Channels(Channel)%Name, Detector, year)
                call HealpixMap_Nullify(yearMaps(year))
                call HealPixMap_read(yearMaps(year), aname)
                call HealpixMap_SetToIndexOnly(yearMaps(year),1)
                yearMaps(year)%TQU = yearMaps(year)%TQU*map_scale/mK
 
       end do         

       do year = 1, nyears 
                call HealpixMap_Assign(M, yearMaps(1))
                M%TQU = M%TQU*real(nyears-year+1)/(nyears-year+2)
                do i=year+1,nyears+1
                 M%TQU = M%TQU - yearMaps(i)%TQU/(nyears-year+2)
                end do

                do weight = 1, nweights
                 ix = ix + 1
                 call HealpixMap_Nullify(WMaps(ix))
                 call HealpixMapMulCut(M,WeightMaps(weight),WMaps(ix), 1)
                end do

                call HealpixMap_Free(M)
                
       end do !year
   end do !detector
  end do !channel

 end subroutine GetWeightedDataNulltestMaps

 
  subroutine SimulateChannelsUnlensed(H, MapArray, P, nulltest)           
   Type(HealpixInfo) :: H
   Type(HealpixMap) :: MapArray(:), SignalMap, YearMaps(nyears)
   Type(HealpixPower) :: P, WhiteP
   integer channel, detector, year, nullmaps, ix, i
   logical, intent(in) :: nulltest
   Type(HealpixAlm) :: A, PtSrcA, SmoothA
   
   call HealpixAlm_Sim(A, P, HasPhi=.false., dopol = want_pol)
   if (point_source_A /= 0) then 
     !Crude idealization where point source are just Gaussian white random field.
    call HealpixPower_Assign(WhiteP,P)
    WhiteP%CL=1
    call HealpixAlm_Sim(PtSrcA, WhiteP, HasPhi=.false., dopol = .false.)
    call HealpixPower_Free(WhiteP)
   end if
   print *,'doing each map'   
   ix = 0
   do channel=1, nchannels
    do detector = 1, Channels(channel)%Count 
     
     if (sim_signal) then
        call HealpixMap_Nullify(SignalMap)
        call HealpixAlm_Assign(SmoothA, A)
        if (point_source_A /= 0) SmoothA%TEB(1,:,:) = SmoothA%TEB(1,:,:) +  &
                                    sqrt(Channels(channel)%PtSrcA) * PtSrcA%TEB(1,:,:) 
        call HealpixAlm_Smooth_beam(SmoothA,Channels(channel)%DetectorBeams(Detector)%Beam)
        call HealpixAlm2Map(H, SmoothA, SignalMap, nside2npix(H%nside))
        call HealpixAlm_Free(SmoothA)
     end if

     if (nulltest) then
      nullmaps = nyears-1
      do year = 1, nyears
       call HealpixMap_Nullify(YearMaps(Year))
       if (sim_signal) then
        call HealpixMap_Assign(YearMaps(Year), SignalMap)
       else
        call HealpixMap_Init(YearMaps(Year), npix, pol = want_pol) 
       end if
       if (sim_noise) call HealpixMap_AddUncorrelatedNoise(YearMaps(Year), &
                           Channels(channel)%DetectorYearNoiseMaps(detector,year))
      end do
       do year = 1, nullMaps 
                ix = ix+1
                call HealpixMap_Nullify(MapArray(ix))
                call HealpixMap_Assign(MapArray(ix), yearMaps(1))
                MapArray(ix)%TQU = MapArray(ix)%TQU*real(nyears-year)/(nyears-year+1)
                do i=year+1,nyears
                 MapArray(ix)%TQU = MapArray(ix)%TQU - yearMaps(i)%TQU/(nyears-year+1)
                end do
                
       end do !year
       call HealpixMapArray_Free(YearMaps)
      
     else

      do year = 1, nyears
      ix = ix+1
       call HealpixMap_Nullify(MapArray(ix))
       if (sim_signal) then
        call HealpixMap_Assign(MapArray(ix), SignalMap)
       else
        call HealpixMap_Init(MapArray(ix), npix, pol = want_pol) 
       end if
       if (sim_noise) call HealpixMap_AddUncorrelatedNoise(MapArray(ix), &
                           Channels(channel)%DetectorYearNoiseMaps(detector,year))
      end do

     end if

     if (sim_signal) call HealpixMap_Free(SignalMap)
     
    end do
   end do 
   
   call HealpixAlm_Free(A)
   if (point_source_A /= 0) call HealpixAlm_Free(PtSrcA)
   
  end  subroutine SimulateChannelsUnlensed


  subroutine SimCombinedChannelsUnlensed(H, MapArray, P)           
   Type(HealpixInfo) :: H
   Type(HealpixMap) :: MapArray(:)
   Type(HealpixPower) :: P, WhiteP
   integer channel
   Type(HealpixAlm) :: A, PtSrcA, SmoothA
   
   call HealpixAlm_Sim(A, P, HasPhi=.false., dopol = want_pol)
   if (point_source_A /= 0) then 
    call HealpixPower_Assign(WhiteP,P)
    WhiteP%CL=1
    call HealpixAlm_Sim(PtSrcA, WhiteP, HasPhi=.false., dopol = .false.)
    call HealpixPower_Free(WhiteP)
   end if
   
   print *,'doing each map'   
   do channel=1, nchannels
       call HealpixMap_Nullify(MapArray(Channel))
       if (.not. sim_signal) then
        call HealpixMap_Init(MapArray(channel), npix, pol = want_pol) 
       else
        call HealpixAlm_Assign(SmoothA, A)
        if (point_source_A /= 0) SmoothA%TEB = SmoothA%TEB +  sqrt(Channels(Channel)%PtSrcA) * PtSrcA%TEB 
        call HealpixAlm_Smooth_beam(SmoothA,Channels(channel)%Beam%Beam)
        call HealpixAlm2Map(H, SmoothA, MapArray(Channel), nside2npix(H%nside))
        call HealpixAlm_Free(SmoothA)
       end if
       if (sim_noise) call HealpixMap_AddUncorrelatedNoise(MapArray(channel), Channels(channel)%NoiseMap)
   end do 
   
   call HealpixAlm_Free(A)
   if (point_source_A /= 0) call HealpixAlm_Free(PtSrcA)
   
  end  subroutine SimCombinedChannelsUnlensed
   

   
 subroutine GetWeightedMaps(WMaps, MapArray)
 Type(HealpixMap) :: WMaps(:), MapArray(:)
 integer detector, nmaps, weight, year, channel
 integer ix

  nmaps = 0
  ix = 0
  do channel = 1, nchannels  
   do Detector=1, Channels(channel)%Count
       do year = 1, nyears 
             ix =ix + 1
             do weight = 1, nweights
                 nmaps = nmaps + 1
                 call HealpixMap_Nullify(WMaps(nmaps))
                 call HealpixMapMulCut(MapArray(ix),WeightMaps(weight),WMaps(nmaps), 1)
             end do
       end do !year
   end do !detector
  end do !channel

 end subroutine GetWeightedMaps
             

 
 function TotYearWeightMaps(ny)
  integer channel, TotYearWeightMaps
  integer, intent(in), optional :: ny   
  TotYearWeightMaps =0
  do channel=1, nchannels
   TotYearWeightMaps = TotYearWeightMaps + Channels(channel)%Count
  end do
  if (present(ny)) then
   TotYearWeightMaps= TotYearWeightMaps * nweights* ny
  else
   TotYearWeightMaps= TotYearWeightMaps * nweights* nyears
  end if
 end  function TotYearWeightMaps
 
 subroutine CrossPCls(H, Wmaps, HybridP, dowrite)
  Type(HealpixInfo) :: H
  logical dowrite
  integer nmaps, Detector, detector2, weight, weight2
  integer channel,channel2,year, year2
  Type(HealpixPower)::HybridP
  Type(HealpixMap) :: M
  Type(HealpixMap) WMaps(:)
  Type(HealpixCrossPowers) :: PCls
  Type(HealpixPower) :: CHat, TotCl(nchannels*nweights*(nchannels*nweights+1)/2)
  integer i, ix,ix2
  real(dp), allocatable :: TotCount(:,:)
  real(dp) :: clWeight(0:lmax)
  integer njoint, jointix, jointix2, index, Pix, l
  
  print *,'getting  CrossPCls'
  
  nmaps = TotYearWeightMaps()
  
  njoint = nchannels*nweights
  allocate(TotCount(0:lmax,nchannels*nweights*(nchannels*nweights+1)/2))
  TotCount = 0
  do i = 1, njoint*(njoint+1)/2   
   call HealpixPower_Nullify(TotCl(i))
   call HealpixPower_Init(TotCl(i),lmax, want_pol)
  end do
  
  call HealpixPower_Nullify(Chat)

  print *, 'Getting cross powers'

  call HealpixMapSet2CrossPowers(H, WMaps, PCls, nmaps, lmax)
  do i=1, nmaps
     call HealpixMap_Free(WMaps(i))
  end do

  print *, 'Getting Chat'
  ix =0
  do channel = 1, nchannels  
    do Detector=1, Channels(channel)%Count
       do year = 1, nyears 
          do weight = 1, nweights
            ix = ix + 1
            ix2 =0
            do channel2 = 1, nchannels  
                do Detector2=1, Channels(channel2)%Count
                do year2 = 1, nyears 
                    do weight2 = 1, nweights
                    ix2 = ix2 + 1
                    if (ix2>ix) cycle
                    if (channel==channel2 .and. year==year2 .and. detector==detector2) cycle 
                    
                    call PseudoCl_GetCHat(Coupler(sym_ix(nweights,weight,weight2)),PCls%Ps(ix,ix2), Chat)
                    call TBeam_PowerSmooth2(Channels(channel)%DetectorBeams(Detector),&
                        Channels(channel2)%DetectorBeams(Detector2),CHat,+1)

 !Put all together as though from just channel and weight; don't attempt optical detector weighting
                       jointix = (channel-1)*nweights + weight
                       jointix2 = (channel2-1)*nweights + weight2
                       index =    sym_ix(njoint,jointix,jointix2)
                       !Just use inverse noise weight as though C_l were full sky uni-weighted
                       ClWeight(0:lmax) = (Channels(channel)%DetectorBeams(Detector)%Beam * &
                                          Channels(channel2)%DetectorBeams(Detector2)%Beam   &
                                        / ( Channels(channel)%sig0(Detector) * Channels(channel2)%sig0(Detector2)))**2
                       do l=2,lmax
                        TotCl(index)%Cl(l,:) = TotCl(index)%Cl(l,:) + Chat%Cl(l,:)*ClWeight(l)
                       end do
                       TotCount(:,index) = TotCount(:,index) + ClWeight
                    end do
                end do
                end do
            end do

               
          end do
       end do
     end do
   end do

  call HealpixCrossPowers_Free(PCls)

 call HealpixPower_Nullify(HybridP)
 call HealpixPower_Init(HybridP,lmax, want_pol)
       
 print *,'Combining to hybrid'
       
 jointix=0          
 Pix = 0
 do channel=1, nchannels
  do weight = 1, nweights
   jointix=jointix+1
   jointix2=0
   do channel2=1, nchannels
    do weight2 = 1, nweights
           jointix2=jointix2 + 1
           if (jointix2 > jointix) cycle
           Pix = Pix + 1      
           do l=2,lmax
            TotCl(Pix)%Cl(l,:) = TotCl(Pix)%Cl(l,:)/TotCount(l,Pix)
           end do
           if (dowrite) call HealpixPower_Write(TotCl(Pix),concat(trim(file_stem)//'_vec',vec_size,'_c',channel,'_w',weight, &
                             '_cc',channel2,'_ww',weight2,'crossPcl.dat'))
 
           HybridP%Cl(lmin:lmax,C_T)= HybridP%Cl(lmin:lmax,C_T) + &
                    MatMul(HybridMix(1)%Cov(Pix)%C, TotCl(Pix)%Cl(lmin:lmax,C_T))
           if (vec_size >=3) then
                 HybridP%Cl(lmin:lmax,C_C)= HybridP%Cl(lmin:lmax,C_C) + &
                    MatMul(HybridMix(2)%Cov(Pix)%C, TotCl(Pix)%Cl(lmin:lmax,C_C))
                 HybridP%Cl(lmin:lmax,C_E)= HybridP%Cl(lmin:lmax,C_E) + &
                    MatMul(HybridMix(3)%Cov(Pix)%C, TotCl(Pix)%Cl(lmin:lmax,C_E))
                 if (vec_size >=4) then
                 HybridP%Cl(lmin:lmax,C_B)= HybridP%Cl(lmin:lmax,C_B) + &
                     MatMul(HybridMix(4)%Cov(Pix)%C, TotCl(Pix)%Cl(lmin:lmax,C_B))
                 end if
          end if       
    end do
   end do      
  end do
 end do  
 
  deallocate(TotCount)
  do i = 1, njoint*(njoint+1)/2   
   call HealpixPower_Free(TotCl(i))
  end do
 
 end subroutine CrossPCls

 function PtScrcA(chan1,chan2)
  real PtScrcA, alpha
  integer, intent(in) :: chan1, chan2
   alpha = 0 !-0.09
   PtScrcA = point_source_A*r_nu(Channels(chan1)%Ghz)*r_nu(Channels(chan2)%Ghz)*&
          (Channels(chan1)%Ghz*Channels(chan2)%Ghz/40.7**2)**(alpha-2)
 end function PtScrcA

 subroutine CrossPClsChannel(H, TotCl, C,  weight)
  Type(HealpixInfo) :: H
  character(LEN=256) :: fname
  Type(TChannel) :: C   
  integer Detector, weight
  integer year, year2
  Type(HealpixMap) :: M(nyears)
  Type(HealpixMap) WMaps(nyears*C%Count)
  Type(HealpixCrossPowers) :: PCls
  Type(HealpixPower) :: CHat, TotCl
  integer nmaps,nest,ix, ix2, Detector2

  print *,'getting channel cross PCls '//trim(C%Name)
  
  fname = FormatFilename(year_filename_format, C%Name)

  call HealpixPower_Init(TotCl,lmax, want_pol)
  call HealpixPower_Nullify(Chat)
 
  nmaps = nyears * C%Count
  nest = 0
  ix= 0  
  do Detector=1, C%Count

    print *, 'Detector', Detector
    call ReadYearMaps(fname, M, Detector)
    do year=1, nyears
          ix = ix+1
          call HealpixMap_Nullify(WMaps(ix))
          call HealpixMap_SetToIndexOnly(M(year),1)
          M(year)%TQU = M(year)%TQU*map_scale/mK
          call HealpixMapMulCut(M(year),WeightMaps(weight),WMaps(ix), 1)
    end do
    call HealpixMapArray_Free(M)

  end do  
    
   print *, 'Getting cross powers'

   call HealpixMapSet2CrossPowers(H, WMaps, PCls, nmaps, lmax)
   call HealpixMapArray_Free(WMaps)
   print *, 'Getting Chat'

  ix =0
  do Detector=1, C%Count
   do year = 1, nyears
    ix = ix+1
    ix2=0
    do Detector2=1, C%Count
     do year2 = 1, nyears
       ix2 = ix2+1
       if (ix2 >= ix) cycle
       call PseudoCl_GetCHat(Coupler(sym_ix(nweights,weight,weight)),PCls%Ps(ix,ix2), Chat)
       call TBeam_PowerSmooth2(C%DetectorBeams(Detector),C%DetectorBeams(Detector2),CHat,+1)
       nest = nest + 1
       TotCl%Cl = TotCl%Cl + Chat%Cl
    end do
   end do 
   end do
   end do
   call HealpixCrossPowers_Free(PCls)
   TotCl%Cl = TotCl%Cl / nest
  
 end subroutine CrossPClsChannel


 subroutine EstAndSetNoise(C)
  Type(TChannel) :: C
  character(LEN=256) :: fname
  integer DA
  integer year, year2 , i
  integer (KIND=8) :: missing,pixcount
  Type(HealpixMap) :: Mask
  Type(HealpixMap) :: M(nyears)
  real(dp) :: sig0(nyears*(nyears+1)/2)
  real sigvar, hit1,hit2
  integer ix
 
  character(LEN=256) :: cachefile
  
   cachefile = concat(trim(data_dir)//'noise_sig0_nside_',nside,'_'//trim(C%Name)) 
   if (FileExists(CacheFile)) then
   
    print *,'Reading cached noise'
    i=new_file_unit()
    call OpenTxtFile(cachefile, i)
     do DA = 1, C%Count
      read (i,*) C%sig0(DA)
      print *, 'DA', DA, 'sig0 = ', C%sig0(DA)
     end do
    call CloseFile(i)
   
   else

   call HealpixMap_Nullify(Mask)
   call HealpixMap_Read(Mask, fits_mask) 

   fname = FormatFilename(year_filename_format, C%Name)
   print *,'Getting noise: '//trim(C%Name)
   do DA = 1, C%Count
   
   call ReadYearMaps(fname, M, DA)
  
   pixcount = 0
   missing=0
   sig0=0
   ix = 0
   do year = 1, nyears 
   do year2 = year+1,nyears
   ix = ix + 1
   do i=0, M(1)%npix-1
    if (Mask%TQU(i,2)>0) then
     hit1=M(year)%TQU(i,2) 
     hit2=M(year2)%TQU(i,2)
     if (hit1>0 .and. hit2>0) then
      pixcount = pixcount+1
      sig0(ix) = sig0(ix) + (M(year)%TQU(i,1)-M(year2)%TQU(i,1))**2/( 1/hit1+ 1/hit2 ) 
     else 
      missing = missing + 1
     end if
    end if
   end do
 
   end do
   end do
   
   sigvar = sum(sig0)/pixcount
   C%sig0(DA) = sqrt(sigvar)

   print *, 'DA', DA, 'sig0 = ', C%sig0(DA)
!   print *, 'ignored frac', real(missing)/(pixcount+missing)
   call HealpixMapArray_Free(M)

   end do
   
   call HealpixMap_Free(Mask)

    i=new_file_unit()
    call CreateTxtFile(cachefile, i)
     do DA = 1, C%Count
      write (i,*) C%sig0(DA)
     end do
    call CloseFile(i)
   
   end if
   
 end  subroutine EstAndSetNoise


 subroutine CombineYears(C)
  Type(TChannel) :: C
  character(LEN=256) :: aname, outname
  integer year, Detector
  Type(HealpixMap) :: M, MTot
  Type(TBeam) :: Beam
  !Assume second column is hit count

  do Detector = 1,C%Count

  outname = FormatFilename(combined_filename_format, C%Name, Detector)
  if (.not. FileExists(outname)) then  
  print *,'Combining years '//trim(C%Name), Detector
  do year = 1, nyears
   aname = FormatFilename(year_filename_format, C%Name, Detector, Year)

   print *,'reading '//trim(aname)
   if (year==1) then
    call HealPixMap_read(Mtot, aname)
     MTot%TQU(:,1) = MTot%TQU(:,1)*Mtot%TQU(:,2) 
   else
    call HealPixMap_read(M, aname)
    MTot%TQU(:,1) = MTot%TQU(:,1) + M%TQU(:,1)*M%TQU(:,2) 
    MTot%TQU(:,2) = MTot%TQU(:,2) + M%TQU(:,2) 
   end if  
  end do
  call HealpixMap_Free(M)
  where (MTot%TQU(:,2)>0)
   MTot%TQU(:,1) = MTot%TQU(:,1)/MTot%TQU(:,2)
  end where
  
  print *,'writing '//trim(outname)
  call HealpixMap_Write(Mtot, outname) 
  call HealpixMap_Free(MTot)
  
  end if
  
  end do

 end  subroutine CombineYears

 subroutine ProcessDatamap(H, C, M)
  Type(HealpixInfo) :: H
  Type(TChannel), target :: C
  Type(HealpixMap) ::  M
  Type(HealpixMap) :: AMap
  character(LEN=1024) :: cache_name, map_fname, fname   
  integer i
      
     call HealpixMap_Nullify(M)
     if (combined_filename_format=='') call MpiStop('No Data file')
     map_fname = FormatFilename(combined_filename_format, C%Name)
     cache_name = CacheName(map_fname, .true.)
     call StringReplace('%DA%','allDA_signal',cache_name)
         
     if (FileExists(cache_name)) then
                print *,'reading cached data map: '// trim(C%Name)
                call HealpixMap_Read(M,cache_name)
     else
              call healpixMap_Init(M, npix, pol = want_pol)
            
              call HealpixMap_ForceRing(M)
              do i=1, C%Count   
                 fname = FormatFilename(combined_filename_format, C%Name,i)
                 call HealpixMap_Read(AMap, fname)
                 call HealpixMap_ForceRing(AMap)
                 M%TQU(:,1) = M%TQU(:,1) + AMap%TQU(:,1)
              end do
               M%TQU(:,1) = M%TQU(:,1)/C%Count * map_scale/mK
               call HealpixMap_Free(AMap)
               print *,'writing '//trim(cache_name)
               call HealpixMap_Write(M,cache_name)
         end if

 end subroutine ProcessDatamap

 subroutine GetSmoothedNoise(H)
  Type(HealpixInfo) :: H
  Type(HealpixMap) :: Hits
  character(LEN=256):: cache_name
  
  
  call HealpixMap_Nullify(SmoothedNoise)

 if (inv_noise_map/='') then

     if (.not. fullsky_test) then
     
     if (noise_inv_fwhm/=0.d0) then

      cache_name = CacheName('smoothed_noise',.false.)
       
      if (.not. no_cache .and. FileExists(cache_name)) then
        print *,'Reading smoothed noise'
        call HealpixMap_Read(SmoothedNoise,cache_name)
      else
       print *,'Getting smoothed noise map'
       call HealpixMap_Read(Hits,inv_noise_map)
       call HealpixMap_SetToIndexOnly(Hits, 2)
       Hits%TQU = 1/Hits%TQU
       call HealpixMap_udgrade(Hits, SmoothedNoise, nside, pessimistic=.false.)
       call HealpixMap_Smooth(H, SmoothedNoise, SmoothedNoise, 2*lmax, noise_inv_fwhm)
       call DeleteFile(cache_name)
       call HealpixMap_Write(SmoothedNoise,cache_name)
   
      end if
       print *,'min/max smoothed noise = ', minval(SmoothedNoise%TQU(:,1))/maxval(SmoothedNoise%TQU(:,1))    
 
     else
      print *, 'using non-smoothed noise'
      call HealpixMap_Read(SmoothedNoise,inv_noise_map)
      call healpixMap_SetToIndexOnly(SmoothedNoise,2)
      call HealpixMap_ForceRing(SmoothedNoise)
      where (SmoothedNoise%TQU >0) 
      SmoothedNoise%TQU = 1/SmoothedNoise%TQU
      end where
     end if        

     end if

    else
     call HealpixMap_Assign(SmoothedNoise,Channels(1)%NoiseMap)
    end if

 end subroutine GetSmoothedNoise

 subroutine ProcessNoiseMaps(H, C)
  Type(HealpixInfo) :: H
  Type(TChannel), target :: C
  Type(HealpixMap), pointer ::  NoiseMap
  Type(HealpixMap) :: AMap
  character(LEN=1024) :: cache_name, noise_fname, fname   
  integer minnoise
  integer i, year
      
        NoiseMap =>  C%NoiseMap

        call HealpixMap_Nullify(NoiseMap)
        call HealpixMap_Nullify(AMap)
        
        if (noise_filename_format/='') then
        
         noise_fname = FormatFilename(noise_filename_format, C%Name)
         cache_name = CacheName(noise_fname, .true.)
         call StringReplace('%DA%','allDA',cache_name)
         
         if (FileExists(cache_name) .and. .not. cross_spectra) then
                print *,'reading cached noise map: '// trim(C%Name)
                call HealpixMap_Read(NoiseMap,cache_name)
         else
              if (cross_spectra) then
                 allocate(C%DetectorYearNoiseMaps(C%Count,nyears))
              end if
              do i=1, C%Count   
                 fname = FormatFilename(noise_filename_format, C%Name, i) 
                 if (i==1) then
                   call GetNoiseMap(H, NoiseMap, fname, C%sig0(i))
                  else
                   call GetNoiseMap(H, AMap, fname, C%sig0(i))
                    NoiseMap%TQU =  NoiseMap%TQU + AMap%TQU
                  end if

                  if (cross_spectra) then
                   do year = 1, nyears
                    call HealpixMap_Nullify(C%DetectorYearNoiseMaps(i,year))
                    fname = FormatFilename(year_filename_format, C%Name, i, year)
                    call GetNoiseMap(H, C%DetectorYearNoiseMaps(i,year), fname, C%sig0(i))
                   end do
                  end if
              end do
!               NoiseMap%TQU = 1._dp/NoiseMap%TQU
               NoiseMap%TQU = NoiseMap%TQU/C%Count**2
               call HealpixMap_Free(AMap)
               call HealpixMap_ForceRing(NoiseMap)
               
               if (.not. cross_spectra) then
                print *,'writing '//trim(cache_name)
                call HealpixMap_Write(NoiseMap,cache_name)
               end if 
         end if
         call HealpixMap_ForceRing(NoiseMap)
        
         print *,trim(C%Name)//' min/max noise = ', &
            minval(NoiseMap%TQU(:,1)),maxval(NoiseMap%TQU(:,1))    
 
        else

          call HealpixMap_Init(NoiseMap, npix,pol = .false.)
               NoiseMap%TQU(:,1) = noise*NoiseMap%npix/(HO_fourpi)
    
        end if 
      
       if (uniform_noise) then
      !Isotropic white noise, no cut
          print *,'Doing uniform white noise'
          NoiseMap%TQU(:,1) = 1/(sum(1/dble(NoiseMap%TQU(:,1)), mask=NoiseMap%TQU(:,1)>0)/count(NoiseMap%TQU(:,1)>0))  
          white_NL = NoiseMap%TQU(1,1)*HO_fourpi/NoiseMap%npix
       else
          white_NL =  1/(sum(1/dble(NoiseMap%TQU(:,1)), mask=NoiseMap%TQU(:,1)>0)/count(NoiseMap%TQU(:,1)>0)) &
                    *HO_fourpi/NoiseMap%npix
       end if

  
      if (want_pol) then
            !same noise for testing
           white_NL_P = ENoiseFac*white_NL
           print *,'White_NL_P = ',white_NL

           call HealpixMap_AddPol(NoiseMap)             
          NoiseMap%TQU(:,2) = NoiseMap%TQU(:,1)*ENoiseFac 
          NoiseMap%TQU(:,3) = NoiseMap%TQU(:,1)*ENoiseFac
      end if


 end subroutine ProcessNoiseMaps


 subroutine GetNoiseMapPlanck(H, NoiseMap, noise_map)
  character(LEN=*), intent(in) :: noise_map
  Type(HealpixInfo) :: H
  Type(HealpixMap) :: NoiseMap 
  Type(HealpixMap) :: Hits
  real(dp) :: meanhit, booknoise
  real(sp) :: minnoise
      print *,'Converting noise map planck, probably broken now'
      
      if (Channels(1)%beam%beam_transfer) call MpiStop('not done beams for Planck')
      
            call HealpixMap_Nullify(NoiseMap)
            call HealpixMap_Read(NoiseMap, noise_map)
            print *,'counts nside =', NoiseMap%nside

            call HealpixMap_Nullify(Hits)
            call HealpixMap_Init(Hits, NoiseMap%npix,nested=NoiseMap%ordering==ord_nest,pol = .false.)
                
!            if (want_pol) then
 !               Hits%TQU(:,2) = NoiseMap%TQU(:,7)/mK**2 
  !              Hits%TQU(:,3) = NoiseMap%TQU(:,9)/mK**2
   !         end if

            Hits%TQU(:,1) =  NoiseMap%TQU(:,4)/mK**2
                
            print *,'rescaling noise to science case'    
            Hits%TQU= Hits%TQU/16   !!!Scale to bluebook
                
            meanhit = sum(dble(Hits%TQU(:,1)))/Hits%npix
            print *,'Mean pix noise =', meanhit
            booknoise =  mK*sqrt(meanhit * ( HO_fourpi/Hits%npix)/ (7.1/60./180.*HO_pi)**2 )/2.726
            print *,'Mean sigma per 7.1fhwm^2 muK/K:',booknoise
                
            call HealpixMap_Free(NoiseMap)
            !  call HealpixVis_Map2ppmfile(Hits, 'outfiles/noisemap.ppm')
                
            call HealpixMap_udgrade(Hits, NoiseMap, nside, pessimistic=.false.)
            NoiseMap%TQU = NoiseMap%TQU* (real(npix)/hits%npix)      
                
            if (nside /= Hits%nside) then
                print *,'smoothing High-res map'
                minnoise = minval(NoiseMap%TQU(:,1))
!               call HealpixMap_Smooth(H, NoiseMap, NoiseMap, lmax, apodize_mask_fwhm)
 
 !Be careful smoothing Q/U noise, don't want to miss monopole and dipole
                call HealpixMap_Smooth(H, NoiseMap, NoiseMap, 2*lmax, apodize_mask_fwhm)
                where (NoiseMap%TQU < minnoise) 
                 NoiseMap%TQU = minnoise
                end where
                
            end if

            call HealpixMap_Free(Hits)
            call HealpixMap_ForceRing(NoiseMap)      

    end subroutine GetNoiseMapPlanck

 subroutine GetNoiseMap(H, NoiseMap, noise_map, noise_sig0)
  character(LEN=*), intent(in) :: noise_map
  Type(HealpixInfo) :: H
  Type(HealpixMap) :: NoiseMap 
  Type(HealpixMap) :: Hits
  real(dp) :: meanhit, booknoise
  real(sp) :: minnoise
  real(dp), intent(in) :: noise_sig0
  integer i
      print *,'converting noise'  // trim(noise_map)
      if (want_pol) call MpiStop('WMAP only TT at the mo')
     
            call HealpixMap_Nullify(NoiseMap)
            call HealpixMap_Read(NoiseMap, noise_map)
            print *,'counts nside =', NoiseMap%nside

            call HealpixMap_Nullify(Hits)
            call HealpixMap_Init(Hits, NoiseMap%npix,nested=NoiseMap%ordering==ord_nest,pol = .false.)
            print *,'Number of zero-hit pixels', count( NoiseMap%TQU(:,2)==0)    
              
                   
            where (NoiseMap%TQU(:,2)>0)       
             Hits%TQU(:,1) = noise_sig0**2*map_scale**2/mK**2/ NoiseMap%TQU(:,2) 
            end where
            
            meanhit = sum(dble(Hits%TQU(:,1)), mask=Hits%TQU(:,1)>0)/count(Hits%TQU(:,1)>0)
            print *,'Mean unmasked pix noise =', meanhit
            
            booknoise =  mK*sqrt(meanhit * ( HO_fourpi/Hits%npix)/ (7.1/60./180.*HO_pi)**2 )/2.726
            print *,'Mean sigma per 7.1fhwm^2 muK/K:',booknoise
                
            call HealpixMap_Free(NoiseMap)
            !  call HealpixVis_Map2ppmfile(Hits, 'outfiles/noisemap.ppm')
                
            call HealpixMap_udgrade(Hits, NoiseMap, nside, pessimistic=.false.)
            NoiseMap%TQU = NoiseMap%TQU* (real(npix)/hits%npix)      
                
            if (nside /= Hits%nside) then
              call MpiStop('noise nside mismatch')                
            end if

            call HealpixMap_Free(Hits)
            call HealpixMap_ForceRing(NoiseMap)      

    end subroutine GetNoiseMap


   subroutine GetLensedFourpointCov(Cov, PUnlensed, PLensed, lmax, lmax_lensed)
   !See astro-ph/0105117; code not used, was just checking
    type(HealpixPower) :: PUnlensed, PLensed
    integer, intent (in) :: lmax, lmax_lensed
    integer l1, l2, lplus, lminus, L, sgn
    Type(TCovMat) :: Cov
    real(dp) threej(0:lmax*2),llp1,llp12, Lfac    
    real(dp) phiphi(lmax_lensed), tmp    
    
    do L=1, lmax_lensed
      phiphi(L) = (2*L+1)*PUnLensed%PhiCl(L,1)/(2*HO_fourpi)
    end do
    
    allocate(Cov%C(2:lmax,2:lmax))

    do l1 = 2, lmax
     llp1 = l1*(l1+1)
     do l2 = 2,l1
      llp12 = l2*(l2+1)
 
       lplus =  min(l1+l2,lmax_lensed)
       lminus = abs(l1-l2)
       call GetThreeJs(threej(lminus:),l1,l2,0,0)          
       tmp=0
       sgn = (-1)**(l1+l2+max(1,lminus)) 
       do L=max(1,lminus), lplus
        lfac = L*(L+1)
        tmp = tmp + sgn*phiphi(L)*( PLensed%Cl(l2,C_T)* (Lfac + llp12 - llp1) + &
           sgn* PLensed%Cl(l1,C_T)*(Lfac + llp1 - llp12))**2*threej(L)**2
        sgn = -sgn
       end do
      Cov%C(l1,l2) = tmp
      Cov%C(l2,l1) = tmp
     end do
   end do  
   
   end subroutine GetLensedFourpointCov          


    subroutine ppm_masked_map(WM,WeightMap, fname, range) 
     use HealpixVis
     Type(HealpixMap) :: WM, WeightMap
     character(LEN=*) :: fname 
     type(HealpixPPM):: ppm, ppmmask
     real, intent(in), optional :: range
     integer i,j
     HealpixVis_force_range = .true.
     if (present(range)) then
      HealpixVis_abs_max = range
     else
      HealpixVis_abs_max = 500.
     end if
     call HealpixVis_Map2ppm(WM,ppm, n=4,plot = splot_none,symmetric = .true.)
     HealpixVis_force_range = .false.
     call HealpixVis_Map2ppm(WeightMap,ppmmask, n=1,plot = splot_none,symmetric = .false.)    
     do i=1,800
      do j=1,400
       if (ppmmask%rgb(3,i,j)/=0 .and.  ppmmask%rgb(1,i,j)==0) then
         ppm%rgb(:,i,j)=0
       end if  
      end do
     end do 
    call HealpixVis_ppm_write(ppm,  fname)
    call HealpixVis_ppm_Free(ppm)
    HealpixVis_force_range = .false.
    end subroutine ppm_masked_map

   subroutine TestExactLike(H,WeightMap, NoiseMap, M, fid_cl_file, sim_cl_file)
    use CutSkyAsymm
    Type(HealpixInfo) :: H
    Type(HealpixMap) :: ProjMap, WeightMap, NoiseMap, NWMap,M, WM
    character(LEN=*), intent(in) :: fid_cl_file,sim_cl_file
    Type (VecArray) :: MapModes(3), SNModes(3)
    Type (HealpixAlm) :: MaskA, NoiseA, A, MapA, MapAProj
    Type(HealpixPower) :: PFid, P, ProjCl
    Type(HealpixPower):: CheckCls1, CheckCls2
    integer l,lmax, nmodes
    Type(ProjMat) :: Proj, TheoryProj, DataProj, SNProj
    Type(ProjMatPol) :: ProjPol, TheoryProjPol, SNProjPol
    integer  npol
    real(dp), allocatable :: modes(:)
    Type(TCovMat) :: Cov, BigNoiseCov, BigCov, FiducialChol, NoiseCov,ModeNoise,JustNoiseCov,&
      NoiseCovM,HighlNoiseCov,MixNoise
    Type(TComplexCovMat) :: PolCov, JustPolNoiseCov, PolNoiseCov, PolHighlNoiseCov, PolModeNoise
    real(dp), allocatable :: diag(:), BigModes(:), ModeVec(:)
    Type(AsymmCouplings) PseudoNoiseCov, TestW
    real(sp), allocatable :: AlmVec(:), AlmVec2(:) 
    real(dp) chisq,amp, like, term
    real(dp) StTime
    real(dp) highlScaleB,highlScaleT,highlScaleE, highlScaleC
    real(dp), allocatable :: Linv(:,:), tmp(:,:), tmpSN(:,:), TmpDataProj(:,:)
    Complex(dp), allocatable :: PolLinv(:,:), Poltmp(:,:), PoltmpSN(:,:)
    complex(dp) AMode
    integer i,j
    character(LEN=150) :: exact_file
    logical :: noise_cut = .false.
    
    
    if (Channels(1)%beam%beam_transfer) call MpiStop('not done beams for exact')
    
    asymm_pol = want_pol
    
    lmax = l_exact+l_exact_margin
    npol = 1
    if (want_pol) npol=3

    print *,'Doing exact with lmax = ',lmax    
    exact_file = concat(trim(data_dir)//'Proj_lexact',l_exact,'_margin',l_exact_margin)
    if (want_pol) exact_file=concat(exact_file,'_pol')
    if (uniform_noise) exact_file=concat(exact_file,'_uninoise') 
    
    exact_file = concat(exact_file,'.dat')    
    
    call HealpixMap2Alm(H,M,MapA,lmax)    
    call HealpixAlm2Power(MapA,P)

    call ppm_masked_map(M,WeightMap, 'outfiles/weighted_map.ppm') 
  
    
    if (.not. no_cache .and. FileExists(exact_file)) then
        print *, 'Reading cached projection matrix and modes' 
        call OpenFile(exact_file,1,'unformatted')
     
        read (1) TheoryProj%nl, TheoryProj%nr
        allocate(SNModes(1)%V(TheoryProj%nr))
        read(1) SNModes(1)%V
        allocate(TheoryProj%M(TheoryProj%nl,TheoryProj%nr))
        nullify(TheoryProj%RootDiag)
        do i=1, TheoryProj%nr
         read(1) TheoryProj%M(:,i)
        end do

        Read(1) highlScaleT
        
        if (want_pol) then
            read(1) highlScaleE, highlScaleC, highlScaleB
            read (1) TheoryProjPol%nl, TheoryProjPol%nr
            allocate(SNModes(2)%V(TheoryProjPol%nr))
            allocate(SNModes(3)%V(TheoryProjPol%nr))
            read(1) SNModes(2)%V,SNModes(3)%V
            nullify(TheoryProjPol%RootDiag)
            allocate(TheoryProjPol%EProj(TheoryProjPol%nl,TheoryProjPol%nr))
            allocate(TheoryProjPol%BProj(TheoryProjPol%nl,TheoryProjPol%nr))
            do i=1, TheoryProjPol%nr
             read(1) TheoryProjPol%EProj(:,i)
             read(1) TheoryProjPol%BProj(:,i)
            end do
        end if
        
       nmodes = TheoryProj%nr
       if (want_pol) then
        nmodes= nmodes + TheoryProjPol%nr*2
       end if
       
       print *,'nmodes = ',nmodes
       
        allocate(bigNoiseCov%C(nmodes,nmodes))
        allocate(HighlNoiseCov%C(nmodes,nmodes))
        do i=1,nmodes
          read(1) BigNoiseCov%C(1:i,i)
          read(1) HighlNoiseCov%C(1:i,i)          
        end do
        do i=1,nmodes
         do j=i+1, nmodes
           BigNoiseCov%C(j,i) = BigNoiseCov%C(i,j) 
           HighlNoiseCov%C(j,i) = HighlNoiseCov%C(i,j) 
         end do
        end do 
          
        read (1) i
        if (i/=252353) stop 'Bad file' 
   
        close(1)
        print *,'Read file'
   
    else

        print *, 'frac weight < 1e-5', count(WeightMap%TQU(:,1)<1e-5)/real(WeightMap%npix)
        call HealpixMap_Assign(WM, M)
        WM%TQU(:,1)=M%TQU(:,1)*WeightMap%TQU(:,1)
        if (want_pol) then
         WM%TQU(:,2)=M%TQU(:,2)*WeightMap%TQU(:,1)
         WM%TQU(:,3)=M%TQU(:,3)*WeightMap%TQU(:,1)
        end if
     
        
        call HealpixMap2Alm(H,WM,MapA,lmax)    
        call healpixMap_Free(WM)
    
        call HealpixMap2Alm(H,WeightMap, MaskA, lmax*2)

        call HealpixMap_Assign(NWMap,WeightMap)
  
        NWMap%TQU(:,1) = WeightMap%TQU(:,1)**2*NoiseMap%TQU(:,1)
        NoiseScale = (HO_fourpi/dble(NoiseMap%npix)) 
        call HealpixMap2Alm(H,NWMap, NoiseA, lmax*2)
        call HealpixMap_Free(NWMap)

!$ call OMP_SET_NUM_THREADS(3)

        
        print *,'CutSkyAsymm_ModeMatrixFromMap'
        StTime = GeteTime()  
        call CutSkyAsymm_ModeMatrixFromMap(MaskA, Proj, ProjPol, lmax, WellSupported, MapModes, MapA)
        print *,'CutSkyAsymm_ModeMatrixFromMap time',  GeteTime()   - StTime
   
        call HealpixPower_ReadFromTextFile(PFid,fid_cl_file,lmax,pol=want_pol,dolens = .false.)
        call TBeam_PowerSmooth(Channels(1)%beam,PFid,-1)

    !Get Noise    
        StTime = GeteTime()  
        print *,'Get noise covariance'
        call CutSkyAsymm_GetNoiseCovariance(Proj,ProjPol, JustNoiseCov, JustPolNoiseCov,NoiseA, lmax, PseudoNoiseCov)
        if (.not. do_lowl_separation) deallocate(PseudoNoiseCov%WASymm)

        print *,'CutSkyAsymm_GetNoiseCovariance time',  GeteTime()   - StTime
     
    !add noise from high l
        print *,'get high l covariance'
         StTime = GeteTime()
         call CutSkyAsymm_GetCovariance(Proj, HighlNoiseCov, PFid, l_exact+1, lmax, .true.)
         if (fake_noise/=0) then
          print *,'adding fake noise to diagonal noise var'
          do i=1,Proj%nr
            JustNoiseCov%C(i,i) = JustNoiseCov%C(i,i) + fake_noise/NoiseScale
          end do         
         end if
         allocate(noiseCov%C(Proj%nr,Proj%nr))
         NoiseCov%C = NoiseScale*JustNoiseCov%C + HighlNoiseBoost*HighlNoiseCov%C
         print *,'CutSkyAsymm_GetCovariance time',  GeteTime()   - StTime
     
        deallocate(HighlNoiseCov%C)

        StTime = GeteTime() 
        call CutSkyAsymm_GetCovariance(Proj, Cov, PFid, 2,l_exact,.false.)
        print *,'CutSkyAsymm_GetCovariance time',  GeteTime()   - StTime

       print *,'Doing temperature S+N'       
      !S+N
       StTime = GeteTime() 
       
        allocate(Linv(Proj%nr,Proj%nr))
   
        Linv = NoiseCov%C + Cov%C
        deallocate(NoiseCov%C)
        
        ! S + N = L L^T
        ! Get Linv = [L^{-1}]^T, so Linv is upper triangular
        call Matrix_CholeskyRootInverse(Linv, transpose=.true.) 
        print *,'Cholesky root time',  GeteTime()   - StTime

        allocate(tmpSN(Proj%nr,Proj%nr))
        StTime = GeteTime() 
        call Matrix_RotateSymm(Cov%C, LInv, Proj%nr,  tmpSN, triangular = .true.)
        deallocate(Cov%C)
        print *,'Rotate time',  GeteTime()   - StTime
        StTime = GeteTime() 
        !Get modes of Signal/(signal+ noise)
        call CutSkyAsymm_GetSupportedModes(tmpSN, TheoryProj, 1e-3_dp) 
        print *,'CutSkyAsymm_GetSupportedModes',  GeteTime()   - StTime
      
        deallocate(TheoryProj%RootDiag)
        nullify(TheoryProj%RootDiag)

        StTime = GeteTime() 
        deallocate(tmpSN)
        SNProj%nl = Proj%nr
        SNProj%nr=TheoryProj%nr
        allocate(SNProj%M(SNProj%nl,SNProj%nr))
        !SNProj(polix)%M = Proj^T L^{-1}
        call Matrix_Mult(LInv,TheoryProj%M,SNProj%M) 
        deallocate(LInv)
 
! Get projected noise

        StTime = GeteTime() 
        allocate(ModeNoise%C(TheoryProj%nr,TheoryProj%nr))
!        !tmpSN = <nn^T> = [ Proj^T L^{-1}] NoiseCov [ Proj^T L^{-1}]^T
        call Matrix_RotateSymm(JustNoiseCov%C, SNProj%M, TheoryProj%nr,  ModeNoise%C)
        deallocate(JustNoiseCov%C)
        print *,'noise cov',  GeteTime()   - StTime
     
!Data
       allocate(SNModes(1)%V(SNProj%nr))
       do i = 1, SNProj%nr
        SNModes(1)%V(i) = dot_product(SNProj%M(:,i),MapModes(1)%V)
       end do
       deallocate(MapModes(1)%V)
   
!Get theory coupling matrix
        StTime = GeteTime() 
        allocate(tmpSN(TheoryProj%nl,TheoryProj%nr))
        do i = 1, Proj%nr
           tmpSN(i,:) = Proj%RootDiag(i)*SNProj%M(i,:)
        end do

        if (do_lowl_separation) then
            allocate(TmpDataProj(TheoryProj%nl,TheoryProj%nr))
            do i=1, Proj%nr
            TmpDataProj(i,:) = SNProj%M(i,:)/Proj%RootDiag(i)
            end do
        end if
        deallocate(SNProj%M)

        print *,'getting big theory proj'
        TheoryProj%nl = Proj%nl
        allocate(TheoryProj%M(TheoryProj%nl,TheoryProj%nr))

        call Matrix_Mult(Proj%M, tmpSN,TheoryProj%M)

        deallocate(tmpSN)
  
        if (do_lowl_separation) then
            print *,'getting big data proj'
            DataProj%nl = Proj%nl
            DataProj%nr = TheoryProj%nr
            allocate(DataProj%M(DataProj%nl,DataProj%nr))
            call Matrix_Mult(Proj%M, tmpDataProj,DataProj%M)
            deallocate(tmpDataProj)
        end if
     
      deallocate(Proj%M)
       
 !!!!Doing polarization      
      if (want_pol) then      
       print *,'Doing Polarization'        

                !for projection into S+N eigenmodes just approx pol cov as with no mixing from EE
                StTime = GeteTime()  
                !Get real part
                call CutSkyAsymm_UnmixedPolCov(ProjPol, PolHighlNoiseCov, PFid, l_exact+1, lmax)
                allocate(PolnoiseCov%C(ProjPol%nr,ProjPol%nr))
                PolNoiseCov%C = (ENoiseFac*NoiseScale)*JustPolNoiseCov%C + HighlNoiseBoost*PolHighlNoiseCov%C
                print *,'CutSkyAsymm_UnmixedPolCov time',  GeteTime()   - StTime
                deallocate(PolHighlNoiseCov%C)

                StTime = GeteTime() 
                call CutSkyAsymm_UnmixedPolCov(ProjPol, PolCov, PFid, 2,l_exact)
                print *,'CutSkyAsymm_UnmixedPolCov time',  GeteTime()   - StTime

                !S+N
                StTime = GeteTime() 
                
                allocate(PolLinv(ProjPol%nr,ProjPol%nr))

                PolLinv = PolNoiseCov%C + PolCov%C
                deallocate(PolNoiseCov%C)
                
                ! S + N = L L^\dag
                ! Get Linv = [L^{-1}]^\dag, so Linv is upper triangular
                call Matrix_CCholeskyRootInverse(PolLinv, dagger=.true.) 
                print *,'CCholesky root time',  GeteTime()   - StTime

                
                allocate(PoltmpSN(ProjPol%nr,ProjPol%nr))
                StTime = GeteTime() 
                call Matrix_CRotateSymm(PolCov%C, PolLInv, ProjPol%nr,  PoltmpSN, triangular = .true.)
                deallocate(PolCov%C)
                print *,'Rotate time',  GeteTime()   - StTime
                StTime = GeteTime() 
                !Get modes of Signal/(signal+ noise)
                call CutSkyAsymm_GetSuppModesPol(PoltmpSN, TheoryProjPol,1e-4_dp) 
                print *,'CutSkyAsymm_GetSuppModesPol',  GeteTime()   - StTime
                
                deallocate(PoltmpSN)
                deallocate(TheoryProjPol%RootDiag)
                nullify(TheoryProjPol%RootDiag)

                StTime = GeteTime() 

                SNProjPol%nl = ProjPol%nr
                SNProjPol%nr= TheoryProjPol%nr
                allocate(SNProjPol%WComp(SNProjpol%nl,SNProjPol%nr))
                !SNProj(polix)%M = Proj^\dag L^{-1}
                call Matrix_CMult(PolLInv,TheoryProjPol%WComp,SNProjPol%WComp) 
                deallocate(PolLInv)
                deallocate(TheoryProjPol%WComp)
                
            ! Get projected noise

                StTime = GeteTime() 
                print *,'Getting polmodenoise'
                allocate(PolModeNoise%C(TheoryProjPol%nr,TheoryProjPol%nr))
            !        !tmpSN = <nn^T> = [ Proj^T L^{-1}] NoiseCov [ Proj^T L^{-1}]^T
                call Matrix_CRotateSymm(JustPolNoiseCov%C, SNProjPol%WComp, TheoryProjPol%nr,  PolModeNoise%C)
                deallocate(JustPolNoiseCov%C)
                print *,'pol noise cov',  GeteTime()   - StTime
                
            !Data
                allocate(SNModes(2)%V(SNProjPol%nr))
                allocate(SNModes(3)%V(SNProjPol%nr))
              
                do i = 1, SNProjPol%nr
                 AMode = dot_product(cmplx(MapModes(2)%V,MapModes(3)%V),SNProjPol%WComp(:,i))
                  !does sum(cjg(A)*B)
                 SNModes(2)%V(i) = real(AMode)
                 SNModes(3)%V(i) = -aimag(AMode)
                end do
                deallocate(MapModes(2)%V)
                deallocate(MapModes(3)%V)

            !Get theory coupling matrix
                StTime = GeteTime() 
                allocate(PoltmpSN(TheoryProjPol%nl,TheoryProjPol%nr))
                do i = 1, ProjPol%nr
                    PoltmpSN(i,:) = ProjPol%RootDiag(i)*SNProjPol%WComp(i,:)
                end do
                deallocate(SNProjPol%WComp)

                print *,'getting big theory proj'
                TheoryProjPol%nl = ProjPol%nl
                allocate(TheoryProjPol%WComp(TheoryProjPol%nl,TheoryProjPol%nr))
                call Matrix_CMult(ProjPol%WComp, PoltmpSN,TheoryProjPol%WComp)
                deallocate(ProjPol%WComp)
                nullify(ProjPol%WComp)
                deallocate(PoltmpSN)

                allocate(TheoryProjPol%EProj(TheoryProjPol%nl,TheoryProjPol%nr))
                allocate(TheoryProjPol%BProj(TheoryProjPol%nl,TheoryProjPol%nr))
                TheoryProjPol%EProj = real(TheoryProjPol%WComp)  !'E'^T = E^T EProj + B^T BProj
                TheoryProjPol%BProj = aimag(TheoryProjPol%WComp)
                deallocate(TheoryProjPol%WComp)
                nullify(TheoryProjPol%WComp)
                   
      end if

!!!Big vectors
       nmodes = TheoryProj%nr
       if (want_pol) then
        nmodes= nmodes + TheoryProjPol%nr*2
       end if

       StTime = GeteTime() 
       print *,'Get full covariance'
       
       ModeNoise%C = NoiseScale*ModeNoise%C
       PolModeNoise%C = (ENoiseFac*NoiseScale)*PolModeNoise%C
       allocate(bigNoiseCov%C(nmodes,nmodes))
       BigNoiseCov%C=0
       BigNoiseCov%C(1:TheoryProj%nr,1:TheoryProj%nr) = ModeNoise%C
       deallocate(ModeNoise%C)
       if (want_pol) then
        BigNoiseCov%C(TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr,TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr) = & 
           + real(PolModeNoise%C)
       
        BigNoiseCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr,&
              TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr) = real(PolModeNoise%C)

        !BE           
        BigNoiseCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr,&
              TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr) = aimag(PolModeNoise%C)
       
       !EB    
        BigNoiseCov%C(TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr, &
                TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr) = - aimag(PolModeNoise%C)

        deallocate(PolModeNoise%C)
       end if

       print *,'get full high l covariance'         
       call CutSkyAsymm_GetFullCovariance(TheoryProj,TheoryProjPol, HighlNoiseCov, PFid, l_exact+1,lmax)
       
       print *,'Get full high l covariance time',  GeteTime()   - StTime
       
       print *,'writing files'
        
        call CreateFile(exact_file,1,'unformatted')
        write(1) TheoryProj%nl, TheoryProj%nr
        write(1) SNModes(1)%V 
        do i=1, TheoryProj%nr
         write(1) TheoryProj%M(1:TheoryProj%nl,i)
        end do

       highlScaleT = PFid%Cl(l_exact+1,C_T)
       write(1) highlScaleT
        
       if (want_pol) then
            highlScaleB = PFid%Cl(l_exact+1,C_B)
            highlScaleE = PFid%Cl(l_exact+1,C_E)
            highlScaleC = PFid%Cl(l_exact+1,C_C)
 
            write(1) highlScaleE, highlScaleC, highlScaleB

            write (1) TheoryProjPol%nl, TheoryProjPol%nr
            write(1) SNModes(2)%V,SNModes(3)%V
 !           write(1) DiagNoise(2)%V 
            do i=1, TheoryProjPol%nr 
             write(1) TheoryProjPol%EProj(1:TheoryProjPol%nl,i)
             write(1) TheoryProjPol%BProj(1:TheoryProjPol%nl,i)
            end do
        end if


        do i=1,nmodes
          write(1) BigNoiseCov%C(1:i,i)
          write(1) HighlNoiseCov%C(1:i,i)
        end do
       i=252353
       write(1) i
      close(1)
    end if    


   StTime = GeteTime() 

  if (do_lowl_separation) then
   if (want_pol) stop 'not done pol high l projection'
   !Get high-l map with rubbish projected

     print *,'getting lowl map'
   !Theory part
     print *,'theory part of cov'
     allocate(AlmVec(TheoryProj%nl))
     AlmVec =  matmul(TheoryProj%M,SNModes(1)%V)   
     call HealpixVectorMultPower(AlmVec,PFid, lmax)

     call CutSkyAsymm_GetCoupling(TestW, MaskA, lmax, plusonly = .false.)
     AlmVec = matmul(TestW%WAsymm, AlmVec)

   !Noise part
     print *,'noise part of cov'
     allocate(AlmVec2(DataProj%nl))
     AlmVec2 =  matmul(DataProj%M,SNModes(1)%V)   
     AlmVec2 = NoiseScale*matmul(PseudoNoiseCov%WAsymm, AlmVec2)
     AlmVec = AlmVec + AlmVec2
     deallocate(ALmVec2)
     print *, DataProj%nl, TheoryProj%nl, (lmax+1)**2
     call HealpixVector2Alm(AlmVec, MapAProj, lmax)
     call HealpixAlm2Map(H, MapAProj, ProjMap, M%npix)
     call HealpixVis_Map2ppmfile(ProjMap, concat('outfiles/',l_exact,'_lmax',lmax,'lowl_map.ppm'),symmetric=.true.)
     call ppm_masked_map(ProjMap,WeightMap, concat('outfiles/',l_exact,'_lmax',lmax,'lowl_map_masked.ppm'), 200.) 

     ProjMap%TQU = M%TQU - ProjMap%TQU
     call ppm_masked_map(ProjMap,WeightMap, concat('outfiles/',l_exact,'_lmax',lmax,'highl_map.ppm')) 

     MapAProj%TEB(:,0:l_exact,:) = 0
     call HealpixAlm2Map(H, MapAProj, ProjMap, M%npix)
     call HealpixVis_Map2ppmfile(ProjMap, concat('outfiles/',l_exact,'_lmax',lmax,'lowl_map_highl.ppm'),symmetric=.true.)
     call ppm_masked_map(ProjMap,WeightMap, concat('outfiles/',l_exact,'_lmax',lmax,'lowl_map_highl_masked.ppm'),200.) 
   
     call HealpixVector2Alm(AlmVec, MapAProj, lmax)
     call HealpixMap2Alm(H,M,MapA,lmax)    
     MapAProj%TEB = MapA%TEB - MapAProj%TEB
     MapAProj%TEB(:,l_exact+1:lmax,:) = 0
     call HealpixAlm2Map(H, MapAProj, ProjMap, M%npix)
     call ppm_masked_map(ProjMap,WeightMap, concat('outfiles/',l_exact,'_lmax',lmax,'lowl_missing_masked.ppm'),200.) 
      
     
     deallocate(AlmVec)
     call HealpixAlm2Power(MapAProj,ProjCl)
     call HealpixPower_Write(ProjCl,'outfiles/lowl_cl.dat')
     call HealpixAlm_Free(MapAProj)
     call HealpixAlm2Power(MapA,ProjCl)
     call HealpixPower_Write(ProjCl,'outfiles/pcl.dat')

     deallocate(PseudoNoiseCov%WASymm)
   end if

   !Chop to just range of interest
  TheoryProj%nl =(l_exact+1)**2
  if (want_pol) TheoryProjPol%nl = (l_exact+1)**2-4
 
   nmodes = TheoryProj%nr
   if (want_pol) then
    nmodes = nmodes + TheoryProjPol%nr*2 
   end if

   print *,'testing likelihood'
   if (check_cls_file1 /='') then
    call HealpixPower_ReadFromTextFile(checkCls1,check_cls_file1,lmax,pol=want_pol,dolens = .false.)
    call HealpixPower_ReadFromTextFile(checkCls2,check_cls_file2,lmax,pol=want_pol,dolens = .false.)
    call HealpixPower_ReadFromTextFile(PFid,check_cls_file1,lmax,pol=want_pol,dolens = .false.) !dummy allocate
   end if

   allocate(BigModes(nmodes))

   if (get_mean_likes) then
        print *,'Getting mean log likelihoods'
        call HealpixPower_ReadFromTextFile(PFid,sim_cl_file,lmax,pol=want_pol,dolens = .false.)
       !!!
    !    PFid%Cl(2:lmax,C_B)=0
        !For full sky result
        P%Cl = PFid%Cl
        P%Cl(2:lmax,C_T) = P%Cl(2:lmax,C_T) + white_NL
        if (want_pol) then
         P%Cl(2:lmax,C_E) = P%Cl(2:lmax,C_E) + white_NL_P
         P%Cl(2:lmax,C_B) = P%Cl(2:lmax,C_B) + white_NL_P
        end if
        !For modes
        StTime = GeteTime() 
        call CutSkyAsymm_GetFullCovariance(TheoryProj,TheoryProjPol, FiducialChol, PFid, 2,l_exact)
        print *, 'full covariance time', GeteTime()  - StTime
        do j=1,nmodes
         FiducialChol%C(:,j) = FiducialChol%C(:,j)+BigNoiseCov%C(:,j) + HighlNoiseCov%C(:,j)
        end do
        call Matrix_Cholesky(FiducialChol%C)  !Note upper triangular is not zeroed
        print *,'got mean matrix'
   end if
   
   do i=0, 20

    if (check_cls_file1 /='') then
     !Check interpolation between two models
      PFid%Cl = CheckCls2%Cl*(i/20.) + (1-i/20.)*CheckCls1%Cl
    else
     call HealpixPower_ReadFromTextFile(PFid,sim_cl_file,lmax,pol=want_pol,dolens = .false.)
      ! amp = 1+(i-8)/real(l_exact)/2
       amp = i/10.
     
      if (amp<0) cycle
!      PFid%Cl(2:l_exact,C_B) =  PFid%Cl(2:l_exact,C_B)*amp
     
   !   amp = 1+(i-8)/real(l_exact)/2.
      PFid%Cl(2:lmax,C_B) =  PFid%Cl(2:lmax,C_B)*amp
     end if

    call CutSkyAsymm_GetFullCovariance(TheoryProj,TheoryProjPol, BigCov, PFid, 2,l_exact)

        do j=1,nmodes
         BigCov%C(:,j) = BigCov%C(:,j)+BigNoiseCov%C(:,j) 
        end do

   !Scale high l
       BigCov%C(1:TheoryProj%nr,1:TheoryProj%nr) = &
        BigCov%C(1:TheoryProj%nr,1:TheoryProj%nr)  &
          + PFid%Cl(l_exact+1,C_T)/highlScaleT*HighlNoiseCov%C(1:TheoryProj%nr,1:TheoryProj%nr)
     
       if (want_pol) then
        !TE
        BigCov%C(1:TheoryProj%nr,TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr) = & 
          BigCov%C(1:TheoryProj%nr,TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr)  & 
          + sqrt(PFid%Cl(l_exact+1,C_T)/highlScaleT*PFid%Cl(l_exact+1,C_E)/highlScaleE) * &
          HighlNoiseCov%C(1:TheoryProj%nr,TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr)
       
        BigCov%C(TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr,1:TheoryProj%nr) = &
         transpose(BigCov%C(1:TheoryProj%nr,TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr)) 
 
       
       !EE
        BigCov%C(TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr,TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr) = & 
         BigCov%C(TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr,TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr)  & 
          + PFid%Cl(l_exact+1,C_E)/highlScaleE* &
          HighlNoiseCov%C(TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr,TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr)
  
       !BB
        BigCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr,&
                       TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr) = &
              BigCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr, &
                       TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr)  &
             + PFid%Cl(l_exact+1,C_B)/highlScaleB* &
               HighlNoiseCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr, &
              TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr)              

       
        BigCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr,&
              TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr) = &
         BigCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr,&
              TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr)  & 
             + PFid%Cl(l_exact+1,C_E)/highlScaleE* &
               HighlNoiseCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr,&
              TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr)
       !EB    
        BigCov%C(TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr, &
                 TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr) = &
           transpose(BigCov%C(TheoryProj%nr+TheoryProjPol%nr+1:TheoryProj%nr+2*TheoryProjPol%nr,&
              TheoryProj%nr+1:TheoryProj%nr+TheoryProjPol%nr))

        end if
        
   
     if (get_mean_likes) then
         like= 0
        call  Matrix_CholeskyRootInverse(BigCov%C)
        do j=1, nmodes
         like = like - log(BigCov%C(j,j))
       end do
       call Matrix_MultTri(BigCov%C,FiducialChol%C,'Right')
       like = like + sum(BigCov%C**2)/2 
     else  
        do j=1, TheoryProj%nr
        BigModes(j) = SNModes(1)%V(j)
        end do
        if (want_pol) then
        do j=1, TheoryProjPol%nr
        BigModes(TheoryProj%nr+j) = SNModes(2)%V(j)
        BigModes(TheoryProj%nr+TheoryProjPol%nr+j) = SNModes(3)%V(j)
        end do
        end if  
!         print *,'testsum', sum(BigCov%C**2)
        like = Matrix_GaussianLogLike(BigCov%C,BigModes)
     end if
 
    deallocate(BigCov%C)
    chisq = 0
    do l=2, l_exact
     PFid%Cl(l,C_T) = PFid%Cl(l,C_T) + white_NL
      if (want_pol) then
       PFid%Cl(l,C_E) = PFid%Cl(l,C_E) + white_NL_P
       PFid%Cl(l,C_B) = PFid%Cl(l,C_B) + white_NL_P
       
       term = PFid%Cl(l,C_T)*PFid%Cl(l,C_E) - PFid%Cl(l,C_C)**2
       
       chisq = chisq + (2*l+1)* ( &
       (PFid%Cl(l,C_T)*P%Cl(l,C_E) + PFid%Cl(l,C_E)*P%Cl(l,C_T) - 2 *PFid%Cl(l,C_C)*P%Cl(l,C_C))/term &
        + log( term/ (P%Cl(l,C_T)*P%Cl(l,C_E) - P%Cl(l,C_C)**2)) -2)
       
        chisq = chisq + (2*l+1)* (P%Cl(l,C_B)/PFid%Cl(l,C_B) &
                +log(PFid%Cl(l,C_B)/P%Cl(l,C_B)) - 1)
            
      else
       chisq = chisq + (2*l+1)* (P%Cl(l,C_T)/PFid%Cl(l,C_T) + log(PFid%Cl(l,C_T)))
      end if
    end do
    print *, i, chisq/2, like  

end do
 
    print *, 'Time for likes', GeteTime()-StTime
    call MpiStop()
 
   end subroutine TestExactLike


   subroutine TestPix(H, Mask, cls_file)
    !Effect of masked pixel window averaging seems OK modelled by healpix window function at 0.2% for nside=1024 l<800
    !But should also check for inverse noise mask
    !Comment out line freeing mask before calling
    use alm_tools
    real(dp), allocatable :: pixlw(:,:)
    Type(Healpixinfo) :: H
    character(LEN=*) :: cls_file
    Type(HealpixPower) :: Av,P, PHat
    integer npw, simno
    Type (HealpixAlm) :: A
    Type (HealpixMap) :: GradPhi,M, MLow, Mask
    integer nsims, highnside
    logical :: lensed = .true.
    
    print *,'Testing pixel window on cut sky, lmax = ', lmax
    
    highnside=  4*nside
    call HealpixPower_Init(Av,lmax, .false.)
    nsims = 30
    do simno = 1, nsims
    print *,'sim ', simno, 'making high res...'
    call HealpixInit(H, highnside, lmax,.true., w8dir=w8name,method=division_balanced) 
    if (H%MpiID ==0) then !if we are main thread
     
     npw = 4*highnside + 1 
     if (lmax > npw-1) call MpiStop('lmax too large for pixel window')
     allocate(pixlw(0:npw-1,1:3))
     call LoadPixelWindow(pixlw, highnside)
     call HealpixPower_ReadFromTextFile(P,cls_file,lmax,pol=.true.,dolens =lensed)
     P%Cl(:,C_T) = P%Cl(:,C_T)* pixlw(2:lmax,1)**2
     deallocate(pixlw)    
     call HealpixAlm_Sim(A, P, HasPhi=lensed, dopol = want_pol)
     if (lensed) then
               call HealpixMap_Nullify(GradPhi)
               call HealpixAlm2GradientMap(H,A, GradPhi, nside2npix(highnside) ,'PHI')
               call HealpixInterpLensedMap_GradPhi(H,A,GradPhi, M, 1.0, interp_cyl)
               call HealpixMap_Free(GradPhi)
               !call HealpixMap2Alm(H, M, A, lmax)
               !call HealpixAlm2Power(A, PHat)
               !call HealpixPower_Write(PHat, 'outfiles/pix_simulated.dat')
     else
       call HealpixAlm2Power(A, PHat)
       call HealpixPower_Write(PHat, 'outfiles/pix_simulated.dat')
       call HealpixAlm2Map(H, A, M, nside2npix(highnside))
     end if
     call HealpixAlm_Free(A)
     print *,'degrade'
     call HealpixMap_udgrade(M, MLow, nside, pessimistic=.false.)
    end if
    
    call HealpixFree(H)
    
 
    call HealpixInit(H, nside, lmax,.true., w8dir=w8name,method=division_balanced) 
    if (H%MpiID ==0) then !if we are main thread

     print *,'getting degraded power'
     call HealpixMap_ForceRing(Mask)
     MLow%TQU(:,1) = MLow%TQU(:,1) * MAsk%TQU(:,1)
     print *,'alm'
     call HealpixMap2Alm(H, MLow, A, H%lmax)
     call HealpixAlm2Power(A, PHat)
     PHat%Cl(2:lmax,1)= matmul(Coupler(1)%InvT(2:lmax,2:lmax),PHat%Cl(2:lmax,1))
    
     call HealpixPower_Write(PHat, 'outfiles/pix_simulated_low.dat')

     npw = 4*H%nside + 1 
     allocate(pixlw(0:npw-1,1:3))
     call LoadPixelWindow(pixlw, H%nside)
     PHat%Cl(0:lmax,1) = Phat%Cl(0:lmax,1)/pixlw(0:lmax,1)**2
     Av%Cl(:,1) = Av%Cl(:,1) + PHat%Cl(:,1)
     call HealpixPower_Write(PHat, 'outfiles/pix_simulated_low_pixbeam.dat')
     deallocate(pixlw)
    end if   

    call HealpixFree(H)
    
    end do !sims
     
    if (H%MpiID ==0) then
     print *,'writing av, lmax - ', lmax
     Av%Cl = Av%Cl/nsims
     call HealpixPower_Write(Av, 'outfiles/pix_simulated_low_pixbeam_av.dat')
     call MpiStop()
    end if 

   end subroutine TestPix


   subroutine PowersToRows
    integer f, i
    Type (HealpixPower) :: P
    character(LEN=256) :: fname
    real(dm) :: cl(0:1100,2)
    f = new_file_unit()
    call CreateTxtFile('outfiles/auto20_samples.dat',f)
    do i=1,5000
     fname = concat('c:\tmp\data\auto20_lmax1100_nside1024_chanWV_cls_hybrid_sim',i,'.dat')
     call Matrix_Read(fname, cl)
     call Matrix_WriteFileRow(f, real(Cl(2:1100,2),dm),1100-2+1 )
    end do 
    call CloseFile(f)
      
   end subroutine PowersToRows
   
end module WeightMixing


program WeightMixer
 use PseudoCl
 use HealpixObj
 use HealpixVis
 use Random
 use spinalm_tools
 use IniFile
 use AMLUtils
 use PseudoCl
 use MatrixUtils
 use WeightMixing
 
 implicit none
 Type(HealpixMap), allocatable ::  WMaps(:)
 Type(HealpixInfo)  :: H
 Type(HealpixMap)   :: M, CutM, Mask, GradPhi
 Type(HealpixMap), allocatable, target   :: MapArray(:)
 Type(HealpixPower) :: PBinnedDiag, PUnlensed, PFid, P, PSim, &
            CAvg, CSkew, CAvgTemp,CVar,Chat, HybridNoise, HybridP, DataCl
 Type(HealpixPower) :: SigNoiseCov, SignalCov, NoiseCov, NEff, fSkyEff
 Type(HealpixPower), dimension(:), allocatable ::CovP, MaskP
 Type(HealpixAlm)   :: A, A2
 Type (TCouplingMatrix), target :: dummycoupler(1)

 Type(HealpixMap), pointer :: Map1,Map2, amap, NoiseMap 
 Type(HealpixCrossPowers) :: CovPowers, MaskPowers 
 
 Type(HealpixMap) :: BinaryMask
 integer rand_seed
 character(LEN=32) :: weights(5)
 character(LEN=256)  :: input_data_dir,  NuMStr,cache_name,mask_fname
 character(LEN=256)  :: l_stem, cls_file, cls_unlensed_sim_file, cls_unlensed_sim_file_tensor, &
        fid_cls_file, out_file_base, out_file_root, sim_map_file, analysis_root, anastem
 character(LEN=256)  :: beamfile,sim_stem, covstem
 integer :: pol_vec_size, i, and_seed
 logical :: err, debug_files
 
 real(dp) ::  meanN, A_PtSc
 Type(TCovMatPolArray), target  :: CovArr(3)
 Type(TCovMat), pointer :: ACov
 Type(TCovMat) HyCov(3)
 Type(TCovMat) :: LensedCov,ACovHybrid
 real(DP) noise_fac, chisq, mask_lowl_fwhm
 real(dp), allocatable :: tmp(:,:), tmp2(:,:),tmpcov(:,:), simcov(:,:)
 real(dp), allocatable :: bigVec(:)
 real(SP) fac
 integer l,ix,j
 integer ix1,ix2, chan1,chan2,x,y
 
 Type (HealpixPower) :: PointSourceP
 logical :: get_covariance, get_sim, get_analysis, get_hybrid = .true.
 logical :: wmap_like_files = .false.
 logical :: cache_couple = .false.
 integer MpiID, MpiSize, nsim, sim_no
 integer ncovpower, ncrossweights, ncl_tot
 integer polx,poly,cov_x, cov_y
 integer lmin_hybrid_calc
 integer nl, lmin_hybrid, offset
 logical :: lensed_sim
 logical :: do_exact_like = .false.
 logical :: beam_transfer
 Type(TCovMatArray) :: CovInv
 double precision StTime
 integer channel, file_unit
 integer justsignal
 logical :: get_signal_covariance = .false.
 logical :: get_noise_covariance = .false.
 character(LEN=16) :: action = ''
 integer cl_sample_unit
! character(LEN=2000) :: inline
! real count1
! integer count
! 
! call OpenTxtFile('C:\Work\cosmomc\WMAP5\newtau_1.txt',1)
! call CreateTxtFile('z:\chain.txt',2)
! do i=1, 2000
!  read (1,'(a)') inline
!  read(inline, *) count1
!  count = nint(count1)
!  do j=1, count
!   write (2,'(a)') trim(inline)
!  end do
! end do 
! stop


#ifdef MPIPIX
 call mpi_init(i)
 
#endif

  call GetMpiStat(MpiID, MpiSize)

 Ini_Fail_On_Not_Found = .true.
 call Ini_Open(GetParam(1), 3,err)
 if (err) then
#ifdef MPIPIX
    call mpi_finalize(i)
#endif
   stop 'No ini'
 end if

 nside  = Ini_Read_Int('nside')
 npix = 12*nside**2
 if (nside==1024) then 
  healpix_res = 10
 else if (nside==512) then
  healpix_res = 9
 else if (nside == 2048) then
  healpix_res = 11
 else
  call MpiStop('unknown nside res')
 end if   

 nyears = Ini_read_int('nyears')
 InYears = nyears
 nweights = Ini_read_int('nweights')
 do i=1, nweights
   weights(i) = Ini_Read_String_Array('weights',i)
 end do
 input_data_dir = Ini_read_String('input_data_dir')
 data_dir = Ini_read_String('data_dir')

 fullsky_test = Ini_read_Logical('fullsky_test')
 uniform_noise = Ini_Read_Logical('uniform_noise')
 
 lmax   = Ini_Read_Int('lmax')  
 lmin =  Ini_Read_int('lmin')
 nl = lmax-lmin+1
 lmin_hybrid = Ini_Read_Int('lmin_hybrid')
 lmin_hybrid_calc = max(lmin,lmin_hybrid-200)

 map_unit = Ini_read_String('map_unit')
 if (map_unit == 'mK') then
  map_scale = 1d3
 else if (map_unit == 'muK') then
  map_scale = 1.d0
 else
  call MpiStop('unknown map unit')  
 end if
 fid_cls_file = Ini_read_String('fid_cls_file')

 out_file_base = Ini_Read_String('out_file_root')
 out_file_root = 'outfiles/' //trim(out_file_base) 
 analysis_root = Ini_read_String('analysis_root',.false.)
 
 want_pol = Ini_Read_Logical('want_pol')
 rand_seed = Ini_Read_Int('rand_seed')
 
 !Read in arcminute, internally in degrees
 apodize_mask_fwhm = Ini_read_Double('apodize_mask_fwhm')/60
 noise_inv_fwhm =  Ini_read_Double('noise_inv_fwhm')/60
 est_noise = Ini_read_logical('est_noise')
 pol_vec_size= Ini_read_int('pol_vec_size')

 point_source_A = Ini_read_Real('point_source_A')
 point_source_A_frac_error = Ini_read_Real('point_source_A_frac_error')
 
 noise = Ini_read_real('noise')/mK**2
 fake_noise = Ini_read_real('fake_noise')/mK**2
 ENoiseFac = Ini_read_real('noise_E_fac')

 beam_transfer = Ini_read_Logical ('beam_transfer')

 year_filename_format = trim(input_data_dir)//Ini_Read_String('year_filename_format')

 !combined_filename_format is computed maps for each detector added over years
 combined_filename_format = trim(data_dir)//trim(ExtractFileName(year_filename_format))
 call StringReplace('%YEAR%',trim(IntToStr(nyears))//'years',combined_filename_format)
 noise_filename_format = combined_filename_format
 
 beam_filename_format = trim(input_data_dir)//Ini_Read_String('beam_filename_format')
 
 nchannels = Ini_read_Int('nchannels')
 allocate(Channels(nchannels))
 do i = 1, nchannels
  Channels(i)%Name = Ini_Read_String_Array('channel_name',i)
  Channels(i)%Count = Ini_read_int_array('channel_count',i)
  Channels(i)%Ghz = Ini_read_Double_Array('channel_Ghz',i)
  allocate(Channels(i)%sig0(Channels(i)%Count))
  NumStr = Ini_read_String_Array('noise_sig0',i)
  read(NumStr,*) Channels(i)%sig0
  Channels(i)%Beam%beam_transfer = beam_transfer
  if (.not. beam_transfer) then
    Channels(i)%Beam%fwhm = Ini_read_real('fwhm')/60  !Read in arcminute, internally in degrees
  end if
   
 end do

 Ini_Fail_On_Not_Found = .false.
 
 wmap_like_files = Ini_read_logical('wmap_like_files', .false.)
 
 get_covariance = Ini_read_logical('get_covariance')
 
 if (get_covariance) then
  get_signal_covariance = Ini_read_logical('get_signal_covariance')
  get_noise_covariance = Ini_read_logical('get_noise_covariance')
 end if 
 
 get_sim =  Ini_read_logical('get_sim')
 get_analysis = Ini_read_Logical('get_analysis')
 
 if (get_sim) then
   sim_noise = Ini_read_Logical('sim_noise', .true.)
   sim_signal = Ini_read_Logical('sim_signal',.true.)
   nsim = Ini_Read_Int('nsim');
   sim_stem = Ini_read_String('sim_root')
   lensed_sim = Ini_Read_Logical('lensed_sim')
   if (.not. sim_noise .and. .not. sim_signal .or. nsim<1) call mpiStop('Nothing to simulate')
 end if
 
 cross_spectra = Ini_Read_Logical('cross_spectra')
 
 w8name = Ini_Read_String('w8dir')
 
 if (uniform_noise) then
  nweights = 1
  lmin_hybrid_calc = lmin 
 end if  
 Ini_Fail_On_Not_Found= .false.
 
 action = Ini_read_String('action')
 
 debug_files = Ini_read_Logical('debug_files',.false.)

 check_cls_file1 = Ini_read_String('check_cls_file1')
 check_cls_file2 = Ini_read_String('check_cls_file2')
 
 fits_mask = trim(input_data_dir)//Ini_read_String('unapodized_mask')
 processing_mask = trim(input_data_dir)//Ini_read_String('processing_mask')
 inv_noise_map = trim(input_data_dir)//Ini_Read_String('inv_noise_map')
 
 cls_file = Ini_Read_String('cls_file')
 cls_unlensed_sim_file = Ini_read_string('cls_unlensed_sim_file')
 cls_unlensed_sim_file_tensor =  Ini_read_string('cls_unlensed_sim_file_tensor')

 do_exact_like = Ini_Read_Logical('do_exact_like')
 
 if (do_exact_like) then
  WellSupported  = Ini_read_real('min_support')
  print *,'min_support = ', WellSupported
  l_exact = Ini_read_Int('l_exact')
  l_exact_margin = Ini_read_Int('l_exact_margin')
 end if
 
 call Ini_Close

 l_stem = trim(IntToStr(lmax))
 if (want_pol) l_stem=trim(l_stem)//'_pol'
 
 file_stem =  concat(trim(out_file_root)//'_lmax',lmax,'_nside',nside,'_chan')
 do channel=1, nchannels
  file_stem = trim(file_stem)//trim(Channels(channel)%name)
 end do
 if (sim_stem /='') then
  sim_stem =  trim(out_file_root)//'_'//trim(sim_stem)//'_lmax'//trim(IntToStr(lmax)) &
                  //'_nside'//trim(IntTOStr(nside))
 else
   sim_stem = file_stem
 end if
 if (want_pol) file_stem=trim(file_stem)//'_pol'
  
 call SetIdlePriority();

 if (w8name=='') then
  call MpiStop('no w8dir specified')
 end if


 fits_mask = FormatFilename(fits_mask)
 processing_mask = FormatFilename(processing_mask)
 

!############### Do Stuff #####################

 ncrossweights = nweights*(nweights+1)/2
 ncl_tot = nweights*nchannels*(nweights*nchannels+1)/2
 if (want_pol) then
    vec_size=pol_vec_size
    else
    vec_size=1
 end if

 !Get azimuthally averaged noise map
! call HealpixInit(H,nside, lmax,.true., w8dir=w8name,method=division_balanced) 
! if (H%MpiID ==0) then !if we are main thread
!   call AzimuthalMap(H, inv_noise_map)
! end if
! call MpiStop('')

 call HealpixInit(H,nside, 2*lmax,.true., w8dir=w8name,method=division_balanced) 

 if (H%MpiID ==0) then !if we are main thread
 
        print *,'Using nside = ', nside
        print *,'Point source A = ', point_source_A, ' +- ',point_source_A_frac_error*100,'%'

        do channel = 1, nchannels
          if (Channels(channel)%Beam%beam_transfer) then
            beamfile = FormatFilename(beam_filename_format,Channels(channel)%Name) 
            call ReadBeams(Channels(channel), beamfile, lmax)
          end if
          Channels(channel)%PtSrcA =  PtScrcA(channel,channel)
          print *,'Channel '//trim(Channels(Channel)%Name)//' point source A = ', Channels(channel)%PtSrcA 
        end do
        call PixSmoothBeams(lmax)

        if (apodize_mask_fwhm/=0.d0) then
        
        mask_fname = trim(data_dir)//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))//'_fwhm' &
          //trim(RealToStr(real(apodize_mask_fwhm),4))//'_mask.fits'
        file_stem = trim(file_stem)//'_apo'
        l_stem = trim(l_stem)//'_apo'
        
        if (.not. FileExists(mask_fname)) then  
          
            print *,'Reading BinaryMask'
            call HealpixMap_Nullify(BinaryMask)
            call HealpixMap_Read(BinaryMask, fits_mask) 
            !Note cut map is the second map, first map is empty

            if (BinaryMask%nside /= nside) print *, 'Upgrading '
            call HealpixMap_udgrade(BinaryMask, Mask, nside, pessimistic=.false.)
            call HealpixMap_Free(BinaryMask)
            call HealpixMap_ForceRing(Mask)

            print *,'generating apodised high-res mask'  
            call HealpixMap2Alm(H,Mask, A, (lmax*3)/2, map_ix=2)
            call HealpixAlm_Smooth(A,apodize_mask_fwhm)
            call HealpixAlm2Map(H,A,Mask,npix)

            print *,'minmax mask = ', minval(Mask%TQU), maxval(Mask%TQU)
!Note to get variance right we need w^2 to have not too much high frequency.
            Mask%TQU  = max(0.,Mask%TQU*1.2-0.1)  
            Mask%TQU  = min(1.,Mask%TQU)  
 
            call HealpixMap_Smooth(H, Mask, Mask, 2*lmax, apodize_mask_fwhm/2)
            print *,'minmax mask = ', minval(Mask%TQU), maxval(Mask%TQU)
            Mask%TQU  = max(0.,Mask%TQU)
            
            if (processing_mask /='') then
             !Make sure exactly zero in zero-hit pixels
             call HealpixMap_Nullify(BinaryMask)
             call HealpixMap_Read(BinaryMask, processing_mask) 
             call HealpixMap_ForceRing(BinaryMask)
             if (BinaryMask%npix /= Mask%npix) Call MpiStop('processing mask not right res')             
             where (BinaryMask%TQU(:,2)==0) 
              Mask%TQU(:,1)=0
             end where
             call HealpixMap_Free(BinaryMask)
            end if  
            
            call HealpixVis_Map2ppmfile(Mask, 'outfiles/mask.ppm')

            call HealpixMap_Write(MAsk,mask_fname)
          
            print *,'minmax mask = ', minval(Mask%TQU), maxval(Mask%TQU)
        else    
           
           print *,'Reading cached mask'
           call HealpixMap_Nullify(Mask)
           call HealpixMap_Read(Mask,mask_fname)
            print *,'minmax mask = ', minval(Mask%TQU), maxval(Mask%TQU)
       
        end if
           
        else
         !Don't apodize mask
            print *,'Reading BinaryMask'
             call HealpixMap_Read(BinaryMask, fits_mask) 
            !Note cut map is the second map, first map is empty

            print *, 'Upgrading '
            call HealpixMap_udgrade(BinaryMask, Mask, nside, pessimistic=.true.)
            call HealpixMap_Free(BinaryMask)
            call HealpixMap_SetToIndexOnly(Mask,2)
            call HealpixMap_ForceRing(Mask)

      end if
    

    do channel= 1, nchannels
        call CombineYears(Channels(channel))
        if (est_noise) call EstAndSetNoise(Channels(channel))
        call ProcessNoiseMaps(H, Channels(Channel))
    end do 
    call GetSmoothedNoise(H)

    if (fullsky_test) then
           print *,'Doing full sky (no cut)'
           Mask%TQU=1
    end if 
      
    if (fake_noise/=0) then
         print *,'Adding fake noise', fake_noise
         white_NL = white_NL + fake_noise
         !Only add to temp
    end if

       
   allocate(WeightMaps(nweights))
   do i=1,nweights
     call HealpixMap_Nullify(WeightMaps(i))
     call HealpixMap_Assign(WeightMaps(i),Mask)
     
     if (weights(i)=='uniform') then

     else if (weights(i)=='invnoise') then

      meanN = sum(SmoothedNoise%TQU(:,1), mask=SmoothedNoise%TQU(:,1) > 0)/count(SmoothedNoise%TQU(:,1) > 0)
      print *,'mean smoothed noise = ', meanN 
      print *, 'sum mask', sum(WeightMaps(i)%TQU(:,1)) !this stops ifort 10.1 crashing on next line ???
      where (SmoothedNoise%TQU(:,1) > 0)
       WeightMaps(i)%TQU(:,1) = WeightMaps(i)%TQU(:,1)*meanN/( SmoothedNoise%TQU(:,1) )
      end where
  
     else if (weights(i)=='mixedinvnoise') then

       meanN = sum(SmoothedNoise%TQU(:,1), mask=SmoothedNoise%TQU(:,1) > 0)/count(SmoothedNoise%TQU(:,1) > 0)
       WeightMaps(i)%TQU(:,1) = WeightMaps(i)%TQU(:,1)*meanN/(SmoothedNoise%TQU(:,1)+meanN)
     
     else
       call MpiStop('unknown weight: '//trim(weights(i)))
     end if
   
   end do
  
   call HealpixMap_Free(SmoothedNoise)    

   if (action/= pixtest .and. .not. do_exact_like) call HealpixMap_Free(Mask)    
 
   print *,'Getting new weights power'
   call HealpixMapSet2CrossPowers(H, WeightMaps, MaskPowers, nweights, lmax*2)
    
   if (get_covariance) then 
       print *,'getting weight powers for the covariance'
       call PseudoCl_WeightsToCovPowers(H, WeightMaps, Channels, CovPowers, nweights, lmax*2)
   end if

end if  !MpiID=0

if (do_exact_like) then
 !Harmonic low-l exact likelihood (under devel)
 
 print *,'Doing exact with l_exact = ', l_exact, 'margin =', l_exact_margin
 sim_map_file = concat(trim(data_dir)//'sim_map',lmax,'_nside',nside)
 if (want_pol) sim_map_file = concat(sim_map_file,'_pol')
 if (uniform_noise) sim_map_file=concat(sim_map_file,'_uninoise') 
 sim_map_file = concat(sim_map_file,'.fits')
 call InitRandom(rand_seed)
 RandInited = .true.

 if (get_sim .or. .not. FileExists(sim_map_file)) then
  print *,'Simulating map'
  call DeleteFile(sim_map_file)
  call HealpixPower_ReadFromTextFile(PUnlensed,cls_file,lmax,pol=.true.,dolens = .false.)
  call TBeam_PowerSmooth(Channels(1)%Beam,PUnlensed,-1)
  call HealpixAlm_Sim(A, PUnlensed, HasPhi=.false., dopol = want_pol)
  call HealpixAlm2Map(H,A,M, npix)
  if (fake_noise/=0) then
      NoiseMap%TQU(:,1) =  NoiseMap%TQU(:,1) + fake_noise*NoiseMap%npix/(HO_fourpi)
  end if     
  call HealpixMap_AddUncorrelatedNoise(M, NoiseMap)
  if (fake_noise/=0) then
      NoiseMap%TQU(:,1) =  NoiseMap%TQU(:,1) - fake_noise*NoiseMap%npix/(HO_fourpi)
  end if     
  call HealpixMap_Write(M, sim_map_file)     
  print *,'done map simulation'
 else
  call HealpixMap_Read(M, sim_map_file)     
 end if               
 
  call TestExactlike(H,Mask, NoiseMap, M, fid_cls_file, cls_file)
  call HealpixMap_Free(Mask)    

  call mpiStop()               
end if
call HealpixFree(H)

 
 if (get_covariance .and.  (noise_filename_format/='') ) then
 
   if (MpiID==0) then
     print *,'Getting matrices for covariance'
    !Note we don't need all if nchanells >1 because noise is assumed uncorrelated between channels
    !Note covariance is calculated for combined year maps
      
    ncovpower = (nchannels+1)*ncrossweights * ((nchannels+1)*ncrossweights+1)/2
   
    allocate(CovP(nCovPower))
    call CrossPowersToHealpixPowerArray(CovPowers, CovP, dofree=.true.)
     !get all for simplicity, but not all entries are used (noise indep between channels)

    allocate(XiMatrices(nCovPower))
   else
     XiMatrices => dummycoupler
   end if
   call PseudoCl_GetCouplingMatrixArr(XiMatrices, CovP, lmin, lmax, want_pol,ncovpower) !parallelized
   
 end if
 
 if (MpiID==0) then

   print *,'Getting coupling matrix'
   allocate(MaskP(ncrossweights))
   call CrossPowersToHealpixPowerArray(MaskPowers, MaskP, dofree=.true.)
   allocate(Coupler(ncrossweights))
 else
     Coupler => dummycoupler
  end if

  call PseudoCl_GetCouplingMatrixArr(Coupler, MaskP, lmin, lmax, want_pol, ncrossweights) !parallelized

  StTime = GeteTime()  
  call PseudoCl_GetCouplingInversesArr(Coupler,ncrossweights)
  if (mpiId==0) print *,'Coupling inversion time',  GeteTime()   - StTime

! if (cache_couple .and. FileExists(trim(data_dir)//trim(l_stem)//'.couple')) then
!        if (MpiID ==0) then !if we are main thread
!            print *,'Reading cached coupling matrix'
!            call  TCouplingMatrix_Read(Coupler, trim(data_dir)//trim(l_stem)//'.couple',1)
!        end if
! else
! end if
!    if (cache_couple) call  TCouplingMatrix_Write(Coupler, trim(data_dir)//trim(l_stem)//'.couple',1)
   
 if (action==pixtest) then
  call TestPix(H,Mask,cls_file)
  call MpiStop('')
 end if
   
 call HealpixInit(H,nside, lmax,.true., w8dir=w8name,method=division_balanced) 

 if (H%MpiID ==0) then !if we are main thread
   
 if (get_covariance) then
  !Again only for combined-year maps

        do channel = 1, nchannels
         print *,'Get Chat Noise_l: '//trim(Channels(channel)%Name)
        
         allocate(Channels(channel)%NoiseP(ncrossweights))
         call PseudoCl_GetCHatNoise(Channels(channel)%NoiseP, Coupler, WeightMaps, nweights, Channels(channel)%NoiseMap)
        
         do i=1,ncrossweights
           call TBeam_PowerSmooth(Channels(channel)%Beam,Channels(channel)%NoiseP(i),+1)
           call HealpixPower_Write(Channels(channel)%NoiseP(i), &
            trim(file_stem)//Trim(IntToStr(i))//trim(Channels(channel)%Name)//'_noise.dat')
         end do
    
        end do
       
       print *,'reading fiducial power'
       call HealpixPower_Nullify(PFid)
       call HealpixPower_ReadFromTextFile(PFid,fid_cls_file,lmax,pol=want_pol,dolens = .false.)

       print *,'Getting covariance'
    
       call PseudoCl_GetFullCovariance(Coupler, XiMatrices, CovArr(1), PFid,vec_size, Channels, nweights, ENoiseFac, &
             .true.,.true.,point_source_A_frac_error)
       if (get_signal_covariance) then
        call PseudoCl_GetFullCovariance(Coupler, XiMatrices, CovArr(2), PFid,vec_size, Channels, nweights, ENoiseFac, &
            .false.,.true.,0.0)
       end if
       if (get_noise_covariance) then
        call PseudoCl_GetFullCovariance(Coupler, XiMatrices, CovArr(3), PFid,vec_size, Channels, nweights, ENoiseFac, &
            .true.,.false.,0.0)
       end if

       if (.not. get_sim .and. .not. get_analysis) then
             call TCouplingMatrix_ArrayFree(Coupler)
             deallocate(Coupler)
       end if
       call TCouplingMatrix_ArrayFree(XiMatrices)
       deallocate(XiMatrices)
   
       call HealpixPower_Init(SigNoiseCov,lmax, want_pol)
       call HealpixPower_Init(PBinnedDiag,lmax, want_pol)

       if (get_signal_covariance) call HealpixPower_Init(SignalCov,lmax, want_pol)
       if (get_noise_covariance) call HealpixPower_Init(NoiseCov,lmax, want_pol)
       
       if (debug_files) then
       print *,'writing diagonals'
       do justsignal=1,3
      
       if (justsignal==2 .and. .not. get_signal_covariance) cycle 
       if (justsignal==3 .and. .not. get_noise_covariance) cycle 

        do cov_x = 1, ncl_tot
         do cov_y = 1, cov_x 
         
         do polx=1, vec_size
   
          ACov => CovArr(justsignal)%Pol(polx,polx)%Cov(cov_x,cov_y)
          if (.not. associated(ACov%C)) call MpiStop(concat('Not associated',IntToStr(cov_x),' ',IntToStr(cov_y)))
         
           do i=lmin,lmax
            fac=1 !1e6*i*(i+1)/(2*3.1415)
            if (polx==1) then
            ix=C_T
            else if (polx==2) then
            ix=C_C
            else if (polx==3) then
            ix=C_E
            else if (polx==4) then
            ix = C_B
            end if   
            SigNoiseCov%Cl(i,ix) = sqrt(fac*ACov%C(i,i))
            if (i>lmin+5 .and. i< lmax-5) then
            PBinnedDiag%Cl(i,ix) = sqrt(sum(fac*ACov%C(i-5:i+5,i-5:i+5))/11.)
            else
            PBinnedDiag%Cl(i,ix) =SigNoiseCov%Cl(i,ix)        
            end if

          end do
        end do
        covstem = concat(trim(file_stem)//'_vec', vec_size,'_',cov_x,'_',cov_y)
        if (justsignal==2) covstem = trim(covstem)//'_signal'
        if (justsignal==3) covstem = trim(covstem)//'_noise'
        call HealpixPower_Write(SigNoiseCov,trim(covstem)// '_diagcov.dat')
        call HealpixPower_Write(PBinnedDiag,trim(covstem)//'_diagcovbinned.dat')
      
       end do
       end do
      end do
      end if !debug files
          
      print *,'Inverting Chat covariances'

     !CovArr arrays now from lmin_hybrid_calc -> lmax

      allocate(HybridMix(vec_size))

 !Take hybrid mixing only between same polarization spectra
       print *,'Getting hybrid coupling'
       do polx=1, vec_size
       
            StTime = GeteTime()  
            if (polx>1) lmin_hybrid_calc=lmin
            call PseudoCl_InverseCovmatArr(CovArr(1)%Pol(polx,polx),lmin_hybrid_calc,lmax,CovInv)
            print *,'Time for PseudoCl_InverseCovmatArr', GeteTime() - StTime
            StTime = geteTime()

            HybridMix(polx)%n = ncl_tot
            HybridMix(polx)%lmin = lmin
            HybridMix(polx)%lmax = lmax
            allocate(HybridMix(polx)%Cov(ncl_tot))
               
            do cov_x = 1, ncl_tot

                allocate(HybridMix(polx)%Cov(cov_x)%C(lmin:lmax,lmin:lmax))
                HybridMix(polx)%Cov(cov_x)%C=0
                 do cov_y = 1, ncl_tot 
                 if (cov_x>=cov_y) then
                  HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax) &
                       =  HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax) + CovInv%Cov(cov_x,cov_y)%C
                 else
                  HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax) &
                      =  HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax) + &
                        transpose(CovInv%Cov(cov_y,cov_x)%C)
                 end if
                end do
            end do
            
            call TCovMatArray_Free(CovInv)
           
        print *,'Getting joint cov for polix', polx

        allocate(ACovHybrid%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax))
        ACovHybrid%C=0
        do cov_x = 1, ncl_tot
            ACovHybrid%C=ACovHybrid%C + HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax)
        end do  

        print *,'Time to get ACovHybrid', GeteTime() - StTime
       
        print *,'inverting for covariance'
        StTime = geteTime()
        call Matrix_Inverse(ACovHybrid%C)
        print *,'Time to inverse ACovHybrid', GeteTime() - StTime
        StTime = geteTime()
          
        print *,'getting final hybrid matrix'
        
        allocate(tmp(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax))
        allocate(tmp2(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax))

        do cov_x = 1, ncl_tot
            tmp = HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax)
            call Matrix_Mult(ACovHybrid%C,tmp,tmp2)
            do i=lmin_hybrid_calc,lmax
              HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,i) = tmp2(:,i)
            end do
!            HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax) =  &
!              matmul(ACovHybrid%C,HybridMix(polx)%Cov(cov_x)%C(lmin_hybrid_calc:lmax,lmin_hybrid_calc:lmax))
        end do
        deallocate(tmp,tmp2)
           
           
       if (polx==1) then    
        print *,'fixing to uniform for TT at low l'
        !still not understood why noise goes wacky at low l

         do cov_x = 1, ncl_tot
          do l= lmin, lmin_hybrid-1
               HybridMix(polx)%Cov(cov_x)%C(l,:)=0 
               if (cov_x==1) HybridMix(polx)%Cov(cov_x)%C(l,l)=1 
                    !Diagonal up to limin_hybrid in uniform weighting in first channel

          end do    
         end do 
         
       else if(polx==4) then

         print *,'fixing to uniform for BB at l < 120'
         !Getting negative BB estimators for l<100, and hybrid noise going negative
         do cov_x = 1, ncl_tot
          do l= lmin, 120-1
               HybridMix(polx)%Cov(cov_x)%C(l,:)=0 
               if (cov_x==1) HybridMix(polx)%Cov(cov_x)%C(l,l)=1 !Diagonal up to limin_hybrid in uniform weighting

          end do    
         end do 

       end if    

       deallocate(ACovHybrid%C)

        print *,'Time to finish pol', GeteTime() - StTime
    
       call PseudoCl_CovmatArr1DWrite(HybridMix(polx),  concat(trim(data_dir)//trim(out_file_base)//'_hybridmat_pol' &
                   ,polx,'_chan',nchannels,'_w',nweights,'_lmax',lmax,'.dat'))

       end do
       
       allocate(tmp(lmin:lmax,lmin:lmax))
       allocate(tmpcov(lmin:lmax,lmin:lmax))
       
       
      print *,'Getting hybrid covariance'      
      
      do justsignal=1,3
      
       if (justsignal==2 .and. .not. get_signal_covariance) cycle 
       if (justsignal==3 .and. .not. get_noise_covariance) cycle 

      !At low l TT covariance is just uniform weighting
      allocate(HyCov(justsignal)%C(nl*vec_size, nl*vec_size))
      HyCov(justsignal)%C=0
    
      StTime = geteTime()
      
       do polx=1,vec_size
        do poly=1,polx
          if (polx==4 .and. poly<3) cycle
          print *,'polx, poly = ',polx, poly
           
           do cov_x=1, ncl_tot
            do cov_y=1, ncl_tot
              if (polx==poly .and. cov_y>cov_x) then
               call Matrix_Mult_NT(HybridMix(polx)%Cov(cov_x)%C,CovArr(justsignal)%Pol(polx,poly)%cov(cov_y,cov_x)%C,tmp)
              !! deallocate(CovArr%Pol(polx,poly)%cov(cov_y,cov_x)%C)
               call Matrix_Mult_NT(tmp,HybridMix(poly)%Cov(cov_y)%C,tmpcov)
              else
               call Matrix_Mult(HybridMix(polx)%Cov(cov_x)%C,CovArr(justsignal)%Pol(polx,poly)%cov(cov_x,cov_y)%C,tmp)
             !!  deallocate(CovArr%Pol(polx,poly)%cov(cov_x,cov_y)%C)
               call Matrix_Mult_NT(tmp,HybridMix(poly)%Cov(cov_y)%C,tmpcov)
              end if
              HyCov(justsignal)%C((polx-1)*nl+1:polx*nl, (poly-1)*nl+1:poly*nl) = &
                   HyCov(justsignal)%C((polx-1)*nl+1:polx*nl, (poly-1)*nl+1:poly*nl) + tmpcov
!             HyCov(justsignal)%C((polx-1)*nl+1:polx*nl, (poly-1)*nl+1:poly*nl) = &
!HyCov(justsignal)%C((polx-1)*nl+1:polx*nl, (poly-1)*nl+1:poly*nl) + &
!              matmul(HybridMix(polx)%Cov(cov_x)%C,  &
!matmul(CovArr%Pol(polx,poly)%cov(cov_x,cov_y)%C,transpose(HybridMix(poly)%Cov(cov_y)%C)))  
            end do
           end do
          if (polx/=poly) HyCov(justsignal)%C((poly-1)*nl+1:poly*nl, (polx-1)*nl+1:polx*nl) &
                   = transpose(HyCov(justsignal)%C((polx-1)*nl+1:polx*nl, (poly-1)*nl+1:poly*nl))
        end do
       end do 
       
       print *,'freeing CovArr'
       call TCovMatPolArray_Free(CovArr(justsignal))
       print *,'Time for hybrid covariance', getetime() - StTime
        
       end do

       deallocate(tmp, tmpcov)
       
       HyCov(1)%C = HyCov(1)%C*mK**4
       print *, 'Cov 22 = ',HyCov(1)%C(2,2)
       call  MatrixSym_Write_Binary_Single(concat( trim(data_dir)//trim(out_file_base)// &
                        '_hybrid_pol',vec_size,'_lmax',lmax,'_',ncl_tot,'.covmat'),HyCov(1)%C)

       HyCov(1)%C = HyCov(1)%C/mK**4
     
       print *,'writing files'
        SigNoiseCov%Cl=0
        SignalCov%Cl=0 !pure signal variance
        NoiseCov%Cl=0 !pure signal variance
        do i=1,nl
        fac=1 !1e6*i*(i+1)/(2*3.1415)
        SigNoiseCov%Cl(i+lmin-1,C_T) = sqrt(fac*HyCov(1)%C(i,i))
        if (get_signal_covariance) SignalCov%Cl(i+lmin-1,C_T) = sqrt(fac*HyCov(2)%C(i,i))
        if (get_noise_covariance) NoiseCov%Cl(i+lmin-1,C_T) = sqrt(fac*HyCov(3)%C(i,i))
        if (i+lmin-1 >lmin+5 .and. i+lmin-1 < lmax -5)  then
         PBinnedDiag%Cl(i+lmin-1,C_T) = sqrt(sum(fac*HyCov(1)%C(i-5:i+5,i-5:i+5))/11.)
        else
         PBinnedDiag%Cl(i+lmin-1,C_T) =SigNoiseCov%Cl(i,C_T)        
        end if
        if (vec_size>1) then  
         SigNoiseCov%Cl(i+lmin-1,C_C) = &
             sqrt(fac*HyCov(1)%C(nl+i,nl+i))
         SigNoiseCov%Cl(i+lmin-1,C_E) = &
            sqrt(fac*HyCov(1)%C(nl*2+i,nl*2+i))
         if (i+lmin-1 >lmin+5 .and. i+lmin-1 < lmax -5) then
            PBinnedDiag%Cl(i+lmin-1,C_C) = sqrt(sum(fac*HyCov(1)%C(nl+i-5:nl+i+5,nl+i-5:nl+i+5))/11.)
            PBinnedDiag%Cl(i+lmin-1,C_E) = sqrt(sum(fac*HyCov(1)%C(2*nl+i-5:2*nl+i+5,2*nl+i-5:2*nl+i+5))/11.)
            else
            PBinnedDiag%Cl(i+lmin-1,C_C) =SigNoiseCov%Cl(i+lmin-1,C_C)        
            PBinnedDiag%Cl(i+lmin-1,C_E) =SigNoiseCov%Cl(i+lmin-1,C_E)        
         end if
         if (vec_size>3) then
          SigNoiseCov%Cl(i+lmin-1,C_B) = &
            sqrt(fac*HyCov(1)%C(nl*3+i,nl*3+i))
          if (i+lmin-1 >lmin+5 .and. i+lmin-1 < lmax -5)  then
            PBinnedDiag%Cl(i+lmin-1,C_B) = sqrt(sum(fac*HyCov(1)%C(3*nl+i-5:3*nl+i+5,3*nl+i-5:3*nl+i+5))/11.)
            else
            PBinnedDiag%Cl(i+lmin-1,C_B) =SigNoiseCov%Cl(i+lmin-1,C_B)        
          end if

         end if
        end if
        
        end do
        print *,'write file'
        call HealpixPower_Write(SigNoiseCov,trim(file_stem)//'_vec'//trim(IntToStr(vec_size))//'_'// &
               'joint_diagcov.dat')
        if (get_signal_covariance) call HealpixPower_Write(SignalCov,trim(file_stem)//'_vec'//trim(IntToStr(vec_size))//'_'// &
               'joint_diagcov_signal.dat')
        if (get_noise_covariance) call HealpixPower_Write(NoiseCov,trim(file_stem)//'_vec'//trim(IntToStr(vec_size))//'_'// &
               'joint_diagcov_noise.dat')

        call HealpixPower_Write(PBinnedDiag,trim(file_stem)//'_vec'//trim(IntToStr(vec_size))//'_'// &
               'joint_diagcovbinned.dat')
     
      if (wmap_like_files) then  
       file_unit = new_file_unit()
       call CreateTxtFile(trim(data_dir)//trim(out_file_base)//'_wmapstyle_cov.dat',file_unit)
       do i=1, min(1000-lmin+1,nl)
        do j= i+1,min(1000-lmin+1,nl)
          write (file_unit, '(2I5,2E17.7)') i+lmin-1, j+lmin-1, 0.d0, HyCov(1)%C(i,j)/sqrt(HyCov(1)%C(i,i)*HyCov(1)%C(j,j))
        end do
       end do
       call CloseFile(file_unit)
      end if  
     
      deallocate(HyCov(1)%C)          
      if (get_signal_covariance) deallocate(HyCov(2)%C)          
     
      print *,'getting joint noise and point source transfers'

      call HealpixPower_Init(HybridNoise,lmax, want_pol)
      call HealpixPower_Init(PointSourceP,lmax, want_pol)
      
       ix = 0
       ix1 = 0
       do chan1=1, nchannels 
        do x=1,nweights
         ix1 = ix1+1
         ix2=0
         do chan2=1, nchannels
          do y=1,nweights
          ix2=ix2+1
          if (ix2 > ix1) cycle
        
          ix=ix+1  
          
          !Use WMAP5 model for point sources
          if (point_source_A /= 0.d0) then
           A_PtSc = sqrt(Channels(Chan1)%PtSrcA*Channels(Chan2)%PtSrcA)
           do i=lmin,lmax
            PointSourceP%Cl(i,C_T)= PointSourceP%Cl(i,C_T) +  A_PtSc*sum(HybridMix(1)%Cov(ix)%C(i,:))
           end do       
          else
            PointSourceP%Cl(i,C_T)= 0
          end if  
          if (chan1==chan2) then
           i = sym_ix(nweights, x,y)
           HybridNoise%Cl(lmin:lmax,C_T) = HybridNoise%Cl(lmin:lmax,C_T) &
                     + matmul(HybridMix(1)%Cov(ix)%C, Channels(chan1)%NoiseP(i)%Cl(lmin:lmax,C_T))
           if (vec_size>=3) then
             HybridNoise%Cl(lmin:lmax,C_E) = HybridNoise%Cl(lmin:lmax,C_E) &
                     + matmul(HybridMix(3)%Cov(ix)%C, Channels(chan1)%NoiseP(i)%Cl(lmin:lmax,C_E))
            if (vec_size>=4) then
                 HybridNoise%Cl(lmin:lmax,C_B) = HybridNoise%Cl(lmin:lmax,C_B) &
                    + matmul(HybridMix(4)%Cov(ix)%C, Channels(chan1)%NoiseP(i)%Cl(lmin:lmax,C_B))
            end if
           end if 
         end if
         end do
        end do
        end do
       end do
   
       call HealpixPower_Write(HybridNoise, trim(file_stem)//'_hybrid_noise.dat')

       call HealpixPower_Write(PointSourceP, trim(file_stem)//'_PointSources.dat')

      if (get_signal_covariance .and. get_noise_covariance) then 
       call HealpixPower_Init(NEff,lmax, want_pol)
       NEff%Cl(2:lmax,C_T) = (SigNoiseCov%Cl(2:lmax,C_T)/SignalCov%Cl(2:lmax,C_T)-1)* &
               (PFid%Cl(2:lmax,C_T)+PointSourceP%Cl(2:lmax,C_T))
!       if (cross_spectra) NEff%Cl = NEff%Cl*(1 +1./(2*( TotYearWeightMaps()/nweights -1))) !Very approximately
       call HealpixPower_Write(NEff,trim(file_stem)//'_vec'//trim(IntToStr(vec_size))//'_NEff.dat')
       end if
  
    if (get_signal_covariance) then
       call HealpixPower_Init(fSkyEff,lmax, want_pol)
      !One way to define degrees of freedom = fskyEff^2 (following WMAP parameterization for diagonal of covariance)
      !Might be better to include point sources in PFid here
       do i=2,lmax
        fSkyEff%Cl(i,1) = sqrt(2*(PFid%Cl(i,1)+PointSourceP%Cl(i,C_T))**2/(2*i+1)) / SignalCov%Cl(i,1)
       end do
       call HealpixPower_Write(fSkyEff,trim(file_stem)//'_vec'//trim(IntToStr(vec_size))//'_fSkyEff.dat')
      end if
  

        print *,'have done covariance'
        
        if (.not. get_sim .and. .not. get_analysis) then
         call HealpixFree(H)
         call MpiStop()
        end if
  else
    print *,'reading hybrid mix'
    if (get_hybrid) then
     allocate(HybridMix(vec_size))
     do polx=1, vec_size
      call PseudoCl_CovmatArr1DRead(HybridMix(polx), &
        concat(trim(data_dir)//trim(out_file_base)//'_hybridmat_pol',polx,'_chan',nchannels,'_w',nweights,'_lmax',lmax,'.dat'))
     end do
     print *,'reading noise'
     call HealpixPower_Nullify(HybridNoise)
     call HealpixPower_ReadFromTextFile(HybridNoise, trim(file_stem)//'_hybrid_noise.dat',lmax, pol = want_pol)
     call HealpixPower_ReadFromTextFile(PointSourceP, trim(file_stem)//'_PointSources.dat',lmax, pol = want_pol)
     if (point_source_A == 0.d0) PointSourceP%Cl=0
    end if  
  end if
 
!###########Analysis########### 

 if (get_analysis) then

  print *,'Data analysis: reading data maps'

  anastem = file_stem
  if (analysis_root/='') anastem = concat(anastem,'_',analysis_root) 
  if (cross_spectra) then
 
   if (action==nulltest) nyears = InYears - 1
 
   allocate(WMaps(TotYearWeightMaps()))
   if (action==nulltest) then
     call GetWeightedDataNulltestMaps(WMaps)
   else    
     call GetWeightedDataMaps(WMaps)
   end if
   call CrossPCls(H, WMaps, DataCl, .true.)
   deallocate(WMaps)
   nyears = InYears
   call HealpixPower_Write(DataCl,trim(anastem)//'_unsubtracted_cross_cls.dat')
   DataCl%Cl(:,1) = DataCl%Cl(:,1) - PointSourceP%Cl(:,1)
   call HealpixPower_Write(DataCl,trim(anastem)//'_cross_cls.dat')

  else !not cross spectra
    
   allocate(MapArray(nchannels))
   do channel = 1, nchannels
      call HealpixMap_Nullify(MapArray(channel))
      call ProcessDatamap(H, Channels(channel), MapArray(Channel))
   end do

   print *,'Analysing maps'
   
   call AnalyseMap(H, MapArray, DataCl) 
   call HealpixPower_Write(DataCl,trim(anastem)//'_unsubtracted_data_cls.dat')
   DataCl%Cl(:,C_T) = DataCl%Cl(:,C_T) - HybridNoise%Cl(:,C_T) - PointSourceP%Cl(:,C_T)
   call HealpixPower_Write(DataCl,trim(anastem)//'_data_cls.dat')
           
   do channel = 1, nchannels
    call healpixMap_Free(MapArray(channel))
   end do
   deallocate(MapArray)  
 
   end if !not cross spectra
   print *,'Map analysis done'

    if (wmap_like_files .and. (get_signal_covariance .and. get_noise_covariance) ) then
     !Note it is not clear to me what WMAP do with point source contribution to cosmic variance
     !Here I'm putting it in Neff so construacted diagonal covariance is correct
     anastem = trim(data_dir)//trim(out_file_base)
     if (analysis_root/='')  anastem= concat(anastem,'_',analysis_root)
     file_unit = new_file_unit()
     call CreateTxtFile(trim(anastem)//'_wmapstyle_inputs_tt.dat',file_unit)
       do i=2, min(lmax,1000)
          write (file_unit, '(1I5,3E17.7)') i, i*(i+1)*mK**2/(HO_twopi)*DataCl%Cl(i,C_T), &
                      i*(i+1)*mK**2/(HO_twopi)*(NEff%Cl(i,C_T)+PointSourceP%Cl(i,C_T)), fSkyEff%Cl(i,C_T)
       end do
      call CloseFile(file_unit)
    end if

   if (.not. get_sim) call MpiStop()
 end if
 
!###########SIMULATION########### 
  
  if (.not. RandInited) call InitRandom(rand_seed)
  RandInited = .true.
  
  print *,'Sim noise:', sim_noise, 'sim signal:', sim_signal
  
  print *,'reading cl'
  if (lensed_sim) then

   call HealpixPower_ReadFromTextFile(PUnlensed,cls_unlensed_sim_file,lmax,pol=want_pol,dolens = .true.)
   if (cls_unlensed_sim_file_tensor/='') then
    call HealpixPower_ReadFromTextFile(PSim,cls_unlensed_sim_file_tensor,lmax,pol=want_pol,dolens = .false.)
    PUnlensed%Cl(lmin:lmax,:) = PUnlensed%Cl(lmin:lmax,:)  + PSim%Cl(lmin:lmax,:)
   end if
   
  else
  
   call HealpixPower_ReadFromTextFile(PUnlensed,cls_file,lmax,pol=.true.,dolens = .false.)
 
  end if

 
  call HealpixPower_ReadFromTextFile(PSim,cls_file,lmax,pol=want_pol,dolens = .false.)
  if (.not. sim_signal) PSim%Cl=0
  
  if (nsim>1) then
   
   call HealpixPower_Init(CAvg,lmax, want_pol)
   call HealpixPower_Init(CSkew,lmax, want_pol)
   call HealpixPower_Init(CAvgTemp,lmax, want_pol)
   call HealpixPower_Init(Cvar,lmax, want_pol)
   if (sim_signal) call HealpixPower_Init(fSkyEff,lmax, want_pol)
  end if

  print *,'simulating'
    
  if (cross_spectra .or. .not. sim_noise) then
    HybridNoise%Cl=0
  end if

  allocate(simcov(nl*vec_size,nl*vec_size))
  allocate(bigvec(nl*vec_size))
  cl_sample_unit = new_file_unit()
  call CreateTxtFile(trim(sim_stem)//'_cls_hybrid_samples.dat',cl_sample_unit)
  simcov=0 
  do sim_no = 1, nsim

            StTime = GeteTime()  
            
            if (cross_spectra) then
                
                print *,'sim channels'
                if (action==nulltest) then
                 allocate(MapArray(TotYearWeightMaps(nyears-1)/nweights))
                else
                 allocate(MapArray(TotYearWeightMaps()/nweights))
                end if
                call SimulateChannelsUnlensed(H, MapArray, PUnlensed, action==nulltest)           
                if(action==nulltest)  nyears = InYears-1
                allocate(WMaps(TotYearWeightMaps()))
                print *,'get weighted maps'
                call GetWeightedMaps(WMaps, MapArray)
                call HealpixMapArray_Free(MapArray)
                deallocate(MapArray)
                print *, 'get cross Pcls'
                call CrossPCls(H, WMaps, HybridP, .false.)
                nyears = InYears
                HybridP%Cl(:,1) = HybridP%Cl(:,1) - PointSourceP%Cl(:,1) 
                deallocate(WMaps)
            
            else !not cross spectra

             if(sim_no==1) allocate(MapArray(nchannels))
            
            if (lensed_sim .and. sim_signal) then
               call MpiStop('not done lensed sim for multi-channel')
!                print *,'doing lensed sim', sim_no
!               call HealpixAlm_Sim(A, PUnlensed, HasPhi=.true., dopol = want_pol)
!
!               call HealpixAlm2Power(A,P)
!               call HealpixPower_Write(P,concat(trim(sim_stem)//'_cls_simfullsky_alm',sim_no,'.dat'))
!            
!               call HealpixMap_Nullify(GradPhi)
!               call HealpixAlm2GradientMap(H,A, GradPhi,npix,'PHI')
!           
!               call HealpixInterpLensedMap_GradPhi(H,A,GradPhi, MapP, (nside/2048)*1.5, interp_cyl)
!               call HealpixMap_Free(GradPhi)
!               call HealpixAlm_Free(A)
!               
!               if (Channels(1)%Beam%beam_transfer) call MpiStop('Not done beams for lensed')
!               
!               call HealpixMap_Smooth(H, MapP, MapP, lmax, Channels(1)%Beam%fwhm)
!               print *,'Lensed simulation time: ', GeteTime()  - StTime
!           
!               call HealpixMap2Alm(H,MapP,A,lmax)
!               call HealpixAlm2Power(A,P)
!               call HealpixPower_Smooth(P,Channels(1)%Beam%fwhm,+1)
!               call HealpixPower_Write(P,concat(trim(sim_stem)//'_cls_simfullsky',sim_no,'.dat'))
!               chisq = 0
!               do i=50,1500
!                  chisq=chisq+ (2*i+1)*( P%Cl(i,C_T)/PSim%Cl(i,C_T) - log(P%Cl(i,C_T)/PSim%Cl(i,C_T)) -1)
!               end do
!               print *,'sim', sim_no, 'chisq TT= ', chisq
               !  if (sim_no==2) call MpiStop()  
            else           
    
                call SimCombinedChannelsUnlensed(H, MapArray, PUnlensed)           
               
            end if
            
            StTime = GeteTime()     
            call AnalyseMap(H, MapArray, HybridP) 
            HybridP%Cl(:,1) = HybridP%Cl(:,1) - PointSourceP%Cl(:,1) 
            do channel=1, nchannels
              call HealpixMap_Free(MapArray(channel)) 
            end do

        end if
 
            print *,'AnalyseMap time:',  GeteTime() -StTime

            bigvec(1:nl) = (HybridP%Cl(lmin:lmax,C_T) - HybridNoise%Cl(lmin:lmax,C_T) - PSim%Cl(lmin:lmax,C_T))
            if (vec_size>=3) then
              bigvec(nl+1:2*nl) = (HybridP%Cl(lmin:lmax,C_C) - HybridNoise%Cl(lmin:lmax,C_C) - PSim%Cl(lmin:lmax,C_C))
              bigvec(2*nl+1:3*nl) = (HybridP%Cl(lmin:lmax,C_E) - HybridNoise%Cl(lmin:lmax,C_E) - PSim%Cl(lmin:lmax,C_E))
             if (vec_size>=4) then
              bigvec(3*nl+1:4*nl) = (HybridP%Cl(lmin:lmax,C_B) - HybridNoise%Cl(lmin:lmax,C_B) - PSim%Cl(lmin:lmax,C_B))
             end if
            end if
            
            do i=1,nl*vec_size
             simcov(:,i) = simcov(:,i) + bigvec(:)*bigvec(i)
            end do

            CAvg%Cl = CAvg%Cl + HybridP%Cl
            CVar%Cl = Cvar%Cl + (HybridP%Cl - HybridNoise%Cl - PSim%Cl)**2
            CSkew%Cl = CSkew%Cl + (HybridP%Cl - HybridNoise%Cl - PSim%Cl)**3
                
            if (.not. cross_spectra .and. sim_noise) &
             call HealpixPower_Write(HybridP,trim(sim_stem)//'_full_cls_hybrid_sim'//trim(IntToStr(sim_no))//'.dat')
            HybridP%Cl = HybridP%Cl - HybridNoise%Cl 
            call HealpixPower_Write(HybridP,trim(sim_stem)//'_cls_hybrid_sim'//trim(IntToStr(sim_no))//'.dat')
            call Matrix_WriteFileRow(cl_sample_unit, real(HybridP%Cl(lmin:lmax,C_T),dm),lmax-lmin+1)
            call HealpixPower_Free(HybridP)
            
            if (nsim>1) then
              !Do every iteration so we can stop when we like
                CAvgTemp%Cl =  (CSkew%Cl/sim_no)
                CAvgTemp%Cl = sign(abs(CAvgTemp%Cl)**(1/3.),CAvgTemp%Cl) !Safe 1/3 power of negtive number
                call HealpixPower_Write( CAvgTemp,trim(sim_stem)//'_cls_hybrid_skew.dat')
                CAvgTemp%Cl =  CAvg%Cl/sim_no  - HybridNoise%Cl
                call HealpixPower_Write(CAvgTemp,trim(sim_stem)//'_cls_hybrid_avg.dat')
                CAvgTemp%Cl =  sqrt(CVar%Cl/sim_no) 
                call HealpixPower_Write(CAvgTemp,trim(sim_stem)//'_cls_hybrid_var.dat')
                
                if (get_signal_covariance .and. get_noise_covariance .and. sim_signal .and. sim_noise) then
                    !One way to define degrees of freedom = fskyEff^2 (following WMAP parameterization for diagonal of covariance)
                    !Note point source variance is going into fsky
                    do i=2,lmax
                        fSkyEff%Cl(i,1) = sqrt(2*(PSim%Cl(i,1)+PointSourceP%Cl(i,1)+Neff%Cl(i,1))**2/(2*i+1) &
                                     / (CVar%Cl(i,1)/sim_no) )
                    end do
                    call HealpixPower_Write( fskyEff,trim(sim_stem)//'_cls_hybrid_fsky.dat')
                end if
                   
            end if
 end do !Sims
 call CloseFile(cl_sample_unit)
 
   if (nsim>1) then
        simcov = (mK**4/nsim)*simcov
        call  MatrixSym_Write_Binary_Single(concat(trim(data_dir)//trim(out_file_base)//'_sim',nsim,'_hybrid_pol', &
                  vec_size,'_lmax',lmax,'_',ncl_tot,'.covmat'),simcov)

!        if (wmap_like_files .and. (get_signal_covariance .and. get_noise_covariance) &
!                .and. sim_signal .and. sim_noise ) then
!        file_unit = new_file_unit()
!        call CreateTxtFile(concat(trim(data_dir)//trim(out_file_base)//'_sim',nsim,'_wmapstyle_inputs_tt.dat'),file_unit)
!        do i=2, min(lmax,1000)
!            write (file_unit, '(1I5,3E17.7)') i, i*(i+1)*mK**2/(HO_twopi)*DataCl%Cl(i,C_T), &
!                        i*(i+1)*mK**2/(HO_twopi)*NEff%Cl(i,C_T), fSkyEff%Cl(i,C_T)
!        end do
!        call CloseFile(file_unit)
!        end if
     
   end if
   print *,'All done'
 end if

#ifdef MPIPIX
    call HealpixFree(H)
    call mpi_finalize(i)
#endif

#ifdef DEBUG
   write (*,*) 'End of program'
#endif
end program WeightMixer
