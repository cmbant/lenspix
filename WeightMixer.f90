    !Code for calculation of hybrid Pseudo-CL estimators
    !For original code see http://cosmologist.info/weightmixer/
    !Updated for polarization. Buch of other routines for misc purposes, not all well tested
    !This version: AL Dec 2010

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
    use SNModes
    implicit none
    character(LEN=*), parameter :: nulltest='nulltest'
    character(LEN=*), parameter :: pixtest='pixtest'
    integer, parameter :: fname_len = 1024

    character(LEN=fname_len) :: data_dir = 'data/'
    character(LEN=fname_len) :: input_data_dir = ''
    character(LEN=fname_len) :: w8name = '../Healpix_2.00/data/'
    character(LEN=fname_len) :: data_var = ''
    logical :: no_pix_window = .false.
    integer, parameter :: max_detectors = 8
    integer:: nside, healpix_res, lmin,lmax
    integer(I_NPIX) :: npix
    logical :: want_pol, est_noise
    real(dp) :: apodize_mask_fwhm
    integer :: pol_maps = 1
    logical :: sim_noise = .true., sim_signal = .true.
    logical :: fullsky_test = .false.
    logical :: uniform_noise = .false.
    logical :: noise_from_hitcounts = .false.
    logical :: cross_spectra = .false.
    logical :: pol_weights = .false.
    logical :: no_cache = .false.
    real(dp) :: ENoiseFac
    real(dp) :: white_NL, white_NL_P !N_l when testing with white noise

    logical :: want_cl
    logical :: get_exact_cl
    logical :: do_lowl_separation = .false. !project lowl map out of high-l one
    integer :: exact_cl_lmin = 2, exact_cl_lmax =3
    logical :: project_l01 = .true.
    real(dp) :: low_l_fwhm = 7./60
    integer :: l_exact = 20
    integer :: l_exact_margin = 20 !50
    real(dp)  :: WellSupported = 0.99_dp
    real(dp) :: exact_SN_cut = 1e-3_dp
    character(LEN=12) :: WellSupported_txt, exact_SN_cut_txt, exact_tag
    real(dp) :: fake_noise = 0._dp
    real :: noise_colour_factor = 0.0 !assume noise goes like 1+ noise_colour_factor/(2l+1)

    integer :: zero_from_l = 0
    character(LEN=fname_len) :: map_unit = 'muK'
    real(dp) :: map_scale = 1._dp !mutiply map to get muK
    logical :: get_mean_likes = .false.
    logical :: test_big_cut = .false.

    real(dp) :: noise_inv_fwhm = 0.d0
    real :: point_source_A = 0.d0
    real :: point_source_A_frac_error = 0.d0
    logical :: unequal_times_only = .false.
    logical :: subtract_monopole = .false.

    real(sp) :: noise
    logical :: cross_beams = .false.
    integer :: beam_file_column = 2
    character(LEN=fname_len) :: file_stem
    character(LEN=fname_len) :: beam_filename_format, channel_filename_format,detector_filename_format, year_filename_format
    character(LEN=fname_len) :: detector_noise_filename_format, year_noise_filename_format
    character(LEN=fname_len) :: fits_mask, fits_mask_pol !0 and 1 foreground mask
    character(LEN=fname_len) :: processing_mask !optional mask says where hit counts are zero and must be zeroed
    integer :: processing_mask_map =2 !Was 2 for WMAP
    character(LEN=fname_len) :: noise_map_for_window !Map to smooth to get 'inverse noise' weight map if noise_inv_fwhm /=0
    logical :: processing_mask_badpix  = .false. !true to get processing mask from bad pixels
    !Testing
    character(LEN=fname_len) :: check_cls_file1, check_cls_file2

    integer vec_size
    integer nchannels, nweights, nyears, InYears
    Type(TChannel), target, dimension(:), allocatable :: Channels
    Type(HealPixMap), dimension(:), target, allocatable :: WeightMaps, WeightMapsPol, WeightMapsAll
    Type (TCouplingMatrix), dimension(:), pointer :: Coupler, XiMatrices
    Type(TCovMatSet), allocatable :: HybridMix(:)
    Type(HealpixMap), save :: SmoothedNoise
    Type(TCrossBeamSet), target :: CrossBeamSet

    contains

    subroutine MapMulPolWeight(InMap,WeightMap,WeightMapPol,OutMap)
    Type(HealpixMap), intent(in) :: InMap,WeightMap,WeightMapPol
    Type(HealpixMap), intent(out) :: OutMap
    integer i, j

    call HealpixMap_ForceRing(InMap)
    if (InMap%npix /= WeightMap%npix) call MpiStop('MapMulWeight: Map size mismatch')
    call HealpixMap_Init(OutMap,InMap%npix,InMap%nmaps)
    outMap%ordering  = ord_ring

    do j=0, InMap%npix -1
        OutMap%TQU(j,1) = InMap%TQU(j,1) * WeightMap%TQU(j,1)
        if (InMap%nmaps>1) then
            OutMap%TQU(j,2) = InMap%TQU(j,2) * WeightMapPol%TQU(j,1)
            OutMap%TQU(j,3) = InMap%TQU(j,3) * WeightMapPol%TQU(j,1)
        end if
    end do

    end subroutine MapMulPolWeight


    subroutine AnalyseMapAuto(H, Maps, HybridP)
    Type(HealpixInfo) :: H
    Type(HealpixMap) :: Maps(:)
    Type(HealpixPower) :: HybridP, CHat
    Type(HealpixCrossPowers) :: CovPowers, PCls

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
            if (pol_weights) then
                call MapMulPolWeight(Maps(Channel),WeightMaps(i),WeightMapsPol(i),WMaps(ix))
            else
                call HealpixMapMulCut(Maps(Channel),WeightMaps(i),WMaps(ix), 1)
            end if
            !  call HealpixVis_Map2ppmfile(WMaps(i), 'outfiles/wmap.ppm')
        end do
    end do

    call HealpixMapSet2CrossPowers(H, WMaps, PCls, nmaps, lmax,.false.)
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
                    !  call HealpixPower_Write(Chat,'outfiles/Chat.dat')
                    call TBeam_PowerSmooth2(Channels(chan1)%Beam,Channels(chan2)%Beam,CHat,+1)
                    !                call HealpixPower_Write(Chat,'outfiles/Chat_beam.dat')

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
    end subroutine AnalyseMapAuto


    function FormatFilenamePair(FString, Channel1, Channel2) result (formatted)
    character(LEN=*), intent(in) :: FString
    character(LEN=*), intent(in)::  Channel1, Channel2
    character(LEN=fname_len) formatted

    formatted = FString

    call StringReplace('%RES%',IntToStr(healpix_res), formatted)
    call StringReplace('%DVAR%',trim(data_var), formatted)
    call StringReplace('%CHANNEL1%',Channel1,formatted)
    call StringReplace('%CHANNEL2%',Channel2,formatted)

    end function FormatFilenamePair


    function FormatFilename(FString, Channel, Detector, Year) result (formatted)
    character(LEN=*), intent(in) :: FString
    character(LEN=*), intent(in), optional ::  Channel, Detector
    integer, intent(in), optional :: Year
    character(LEN=fname_len) formatted

    formatted = FString

    call StringReplace('%RES%',IntToStr(healpix_res), formatted)
    call StringReplace('%DVAR%',trim(data_var), formatted)
    if (present(Channel)) call StringReplace('%CHANNEL%',Channel,formatted)
    if (present(Detector))  then
        call StringReplace('%DA%',Trim(Detector),formatted)
    end if
    if (present(Year))  call StringReplace('%YEAR%',IntToStr(Year),formatted)

    end function FormatFilename

    function CacheName(FName, nostem) result (cache)
    character(LEN=*) :: FName
    character(LEN=fname_len) :: cache
    logical nostem

    cache = data_dir
    if (.not. nostem) cache = trim(cache)//trim(ExtractFileName(file_stem))
    cache = trim(cache)//trim(ExtractFileName(Fname))
    if (ExtractFileExt(cache)/= '.fits') cache = trim(cache)//'.fits'

    end function CacheName

    subroutine ReadCrossBeamSet(name_template)
    character(LEN=*), intent(in) :: name_template
    integer channel1, channel2
    character(LEN=fname_len) :: beamfile
    Type(TBeam), pointer :: Beam

    allocate(CrossBeamSet%Beams(nchannels,nchannels))
    do channel1 = 1, nchannels
        do channel2 = channel1, nchannels
            beam => CrossBeamSet%Beams(channel1,channel2)
            beam%beam_transfer = .true.
            beamfile = FormatFilenamePair(name_template,Channels(channel1)%Name,Channels(channel2)%Name)
            call TBeam_ReadFile(Beam,beamfile,lmax, col=beam_file_column)
        end do
    end do
    if (.not. no_pix_window) call PixSmoothBeams(lmax)

    end subroutine ReadCrossBeamSet

    subroutine ReadBeams(C, beamfile, lmax)
    !Just weight uniformly for now for effective combined map
    Type(TChannel) :: C
    character(LEN=*), intent(in) :: beamfile
    character(LEN=fname_len) :: file
    integer, intent(in) :: lmax
    integer i
    Type(TBeam) B

    allocate(C%DetectorBeams(C%Count))

    do i=1, C%Count
        file = FormatFilename(beamfile,C%Name, DetectorName(C,i))
        C%DetectorBeams(i)%beam_transfer = .true.
        print *, 'reading beam file: '//trim(file)
        call TBeam_ReadFile(C%DetectorBeams(i),file,lmax)
        if (i==1) then
            call TBeam_ReadFile(C%Beam,file,lmax, col=beam_file_column)
        else
            C%Beam%Beam= C%Beam%Beam + C%DetectorBeams(i)%Beam
        end if
    end do
    C%Beam%Beam = C%Beam%Beam/C%Count

    end subroutine ReadBeams

    subroutine SetGaussianBeams(C, lmax)
    Type(TChannel) :: C
    character(LEN=fname_len) :: file
    integer, intent(in) :: lmax
    integer i
    Type(TBeam) B

    allocate(C%DetectorBeams(C%Count))

    do i=1, C%Count
        C%DetectorBeams(i)%beam_transfer = .true.

        C%DetectorBeams(i)%fwhm = C%Beam%fwhm

        call TBeam_SetGaussian(C%DetectorBeams(i),lmax)
        if (i==1) then
            call TBeam_SetGaussian(C%Beam,lmax)
        else
            C%Beam%Beam= C%Beam%Beam + C%DetectorBeams(i)%Beam
        end if
    end do
    C%Beam%Beam = C%Beam%Beam/C%Count

    end subroutine SetGaussianBeams

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
    integer npw, i, lmax, channel,channel1,channel2

    npw = 4*nside + 1
    if (lmax > npw-1) call MpiStop('lmax too large for pixel window')
    allocate(pixlw(0:npw-1,1:3))
    call LoadPixelWindow(pixlw, nside)
    print *,'pixel window l=2000 or lmax: ',pixlw(min(lmax,2000),1)
    if (cross_beams) then
        do channel1 = 1, nchannels
            do channel2 = channel1, nchannels
                associate(beam => CrossBeamSet%Beams(channel1,channel2)%Beam)
                    beam = beam * pixlw(0:lmax,1)
                    end associate
            end do
        end do
    else
        do channel=1, nchannels
            do i=1, Channels(channel)%Count
                Channels(channel)%DetectorBeams(i)%Beam= Channels(channel)%DetectorBeams(i)%Beam * pixlw(0:lmax,1)
            end do
            Channels(channel)%Beam%Beam= Channels(channel)%Beam%Beam * pixlw(0:lmax,1)
        end do
    end if

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

    subroutine CrossPowersToHealpixPowerArray2(CrossPowers,PowerArray,dofree)
    Type(HealpixCrossPowers) :: CrossPowers
    Type(HealpixPower) :: PowerArray(:,:)
    logical, intent(in), optional :: dofree
    integer i,j,ix,n

    if (CrossPowers%npol>1) call MpiStop('CrossPowersToHealpixPowerArray: not done for pol')
    ix=0
    n = CrossPowers%nmaps / 2
    do i=1, n
        do j=1,i
            ix= ix+1
            call HealpixPower_Init(PowerArray(ix,1), CrossPowers%lmax, CrossPowers%npol==3, dolens=.false., nofree=.true.)
            call HealpixPower_Init(PowerArray(ix,2), CrossPowers%lmax, CrossPowers%npol==3, dolens=.false., nofree=.true.)
            call HealpixPower_Init(PowerArray(ix,3), CrossPowers%lmax, CrossPowers%npol==3, dolens=.false., nofree=.true.)
            PowerArray(ix,1)%Cl(:,C_T) = CrossPowers%Ps(i,j)%Cl(:,1,1)
            PowerArray(ix,2)%Cl(:,C_T) = CrossPowers%Ps(i+n,j)%Cl(:,1,1)
            PowerArray(ix,3)%Cl(:,C_T) = CrossPowers%Ps(i+n,j+n)%Cl(:,1,1)

            if (present(dofree)) then
                if (dofree) then
                    deallocate(CrossPowers%Ps(i,j)%Cl)
                    deallocate(CrossPowers%Ps(i+n,j)%Cl)
                    deallocate(CrossPowers%Ps(i+n,j+n)%Cl)
                end if
            end if
        end do
    end do
    if (dofree) deallocate(CrossPowers%Ps)

    end subroutine CrossPowersToHealpixPowerArray2

    function DetectorName(C, i)
    Type(TChannel), intent(in) :: C
    integer, intent(in) :: i
    character(LEN=48) DetectorName

    if (C%DetectorNames%Count>0) then
        DetectorName=TStringList_Item(C%DetectorNames,i)
    else
        DetectorName=IntToStr(i)
    end if

    end function DetectorName

    subroutine ReadYearMaps(C,fname, M, DA, map_limit)
    character(LEN=*), intent(in):: fname
    character(LEN=fname_len) :: aname
    Type(TChannel), intent(in) :: C
    Type(HealpixMap) :: M(nyears)
    integer year, DA, limit
    integer, intent(in), optional :: map_limit

    do year = 1, nyears
        aname =  FormatFilename(fname, '', DetectorName(C,DA), year)
        call HealpixMap_Nullify(M(year))
        if (present(map_limit)) then
            call HealpixMap_Read(M(year),aname, map_limit)
        else
            call HealpixMap_Read(M(year),aname)
        end if
        call HealpixMap_ForceRing(M(year))
    end do

    end subroutine ReadYearMaps

    subroutine GetWeightedDataMaps(WMaps)
    Type(HealpixMap) :: WMaps(:)
    integer detector, ix, weight, year, channel
    character(LEN=fname_len) :: aname
    Type(HealpixMap) :: M
    real mono

    print *,'Reading maps to get weighted map array'
    ix = 0
    do channel = 1, nchannels
        print *, 'channel '//trim(Channels(channel)%Name)

        do Detector=1, Channels(channel)%Count
            print *, 'Detector', Detector

            do year = 1, nyears
                aname = FormatFilename(year_filename_format, Channels(Channel)%Name, &
                DetectorName(Channels(Channel),Detector), year)
                call HealpixMap_Nullify(M)
                call HealPixMap_read(M, aname, pol_maps)
                call HealpixMap_ForceRing(M)

                if (subtract_monopole) then
                    !!!!
                    mono=sum(dble(M%TQU(:,1)),mask=WeightMaps(1)%TQU(:,1)>0 .and. M%TQU(:,1)/=fmissval)/ &
                    count(WeightMaps(1)%TQU(:,1)>0 .and. M%TQU(:,1)/=fmissval)
                    where (M%TQU(:,1)/= fmissval)
                        M%TQU(:,1)=M%TQU(:,1) - mono
                    end where
                end if

                where (M%TQU/= fmissval)
                    M%TQU = M%TQU*map_scale/mK
                end where

                do weight = 1, nweights
                    ix = ix + 1
                    call HealpixMap_Nullify(WMaps(ix))
                    if (pol_weights) then
                        call MapMulPolWeight(M,WeightMaps(weight),WeightMapsPol(weight),WMaps(ix))
                    else
                        call HealpixMapMulCut(M,WeightMaps(weight),WMaps(ix), 1)
                    end if
                    !                    call HealpixVis_Map2ppmfile(WMaps(ix), 'outfiles/maskmap'//trim(IntToStr(ix))//'.ppm')
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
    character(LEN=fname_len) :: aname

    print *,'Making signal-null maps and weighted array'
    ix = 0
    do channel = 1, nchannels
        print *, 'channel '//trim(Channels(channel)%Name)

        do Detector=1, Channels(channel)%Count
            print *, 'Detector', Detector


            do year = 1, nyears+1
                aname = FormatFilename(year_filename_format, Channels(Channel)%Name, &
                DetectorName(Channels(Channel),Detector), year)
                call HealpixMap_Nullify(yearMaps(year))
                call HealPixMap_read(yearMaps(year), aname, pol_maps)
                call HealpixMap_ForceRing(yearMaps(year))
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
                    if (pol_weights) then
                        call MapMulPolWeight(M,WeightMaps(weight),WeightMapsPol(weight),WMaps(ix))
                    else
                        call HealpixMapMulCut(M,WeightMaps(weight),WMaps(ix), 1)
                    end if
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
    integer maxp

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
                if (want_pol .and. Channels(channel)%detector_has_pol(detector)) then
                    maxp=3
                else
                    maxp=1
                end if
                call HealpixAlm_Assign(SmoothA, A,maxp )
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
                        call HealpixMap_Init(YearMaps(Year), npix, pol = want_pol .and. Channels(channel)%detector_has_pol(detector))
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
                        call HealpixMap_Init(MapArray(ix), npix, pol = &
                        want_pol .and. Channels(channel)%detector_has_pol(detector))
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
                    if (pol_weights) then
                        call MapMulPolWeight(MapArray(ix),WeightMaps(weight),WeightMapsPol(weight),WMaps(nmaps))
                    else
                        call HealpixMapMulCut(MapArray(ix),WeightMaps(weight),WMaps(nmaps), 1)
                    end if
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
    real(dp), allocatable :: TotCount(:,:,:)
    real(dp) :: clWeight(0:lmax)
    integer njoint, jointix, jointix2, index, Pix, l
    Type(HealpixAllCl), pointer :: ThisCross

    print *,'getting CrossPCls'

    nmaps = TotYearWeightMaps()

    njoint = nchannels*nweights
    allocate(TotCount(0:lmax,nchannels*nweights*(nchannels*nweights+1)/2,3))
    TotCount = 0
    do i = 1, njoint*(njoint+1)/2
        call HealpixPower_Nullify(TotCl(i))
        call HealpixPower_Init(TotCl(i),lmax, want_pol)
    end do

    call HealpixPower_Nullify(Chat)

    print *, 'Getting cross powers'

    call HealpixMapSet2CrossPowers(H, WMaps, PCls, nmaps, lmax,.false.)
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
                                    !                    if (weight2 > weight) cycle
                                    !!! if (ix2>ix) cycle
                                    !don't miss alternative TE estimator T year 1 E year 2 and T year 2 E year 1
                                    !Can just double count all temperature estimators; bit inefficient
                                    if (channel==channel2 .and. year==year2 .and. detector==detector2 .or. &
                                    unequal_times_only .and. year==year2) cycle

                                    !Put all together as though from just channel and weight; don't attempt optimal detector weighting
                                    jointix = (channel-1)*nweights + weight
                                    jointix2 = (channel2-1)*nweights + weight2
                                    index =    sym_ix(njoint,jointix,jointix2)
                                    if (weight2 <= weight) then
                                        ThisCross => PCls%Ps(ix,ix2)
                                    else
                                        ThisCross => PCls%Ps(ix2,ix)
                                    end if
                                    call PseudoCl_GetCHat(Coupler(sym_ix(nweights,weight,weight2)),ThisCross, Chat)

                                    call TBeam_PowerSmooth2(Channels(channel)%DetectorBeams(Detector),&
                                    Channels(channel2)%DetectorBeams(Detector2),CHat,+1)
                                    !Just use inverse noise weight as though C_l were full sky uni-weighted
                                    ClWeight(0:lmax) = (Channels(channel)%DetectorBeams(Detector)%Beam * &
                                    Channels(channel2)%DetectorBeams(Detector2)%Beam   &
                                    / ( Channels(channel)%sig0(Detector) * Channels(channel2)%sig0(Detector2)))**2
                                    do l=2,lmax
                                        TotCl(index)%Cl(l,C_T) = TotCl(index)%Cl(l,C_T) + Chat%Cl(l,C_T)*ClWeight(l)
                                    end do
                                    TotCount(:,index,1) = TotCount(:,index,1) + ClWeight
                                    if (size(ThisCross%Cl,2)>1) then
                                        do l=2,lmax
                                            TotCl(index)%Cl(l,C_C) = TotCl(index)%Cl(l,C_C) + Chat%Cl(l,C_C)*ClWeight(l)
                                        end do
                                        TotCount(:,index,2) = TotCount(:,index,2) + ClWeight
                                    end if
                                    if (size(ThisCross%Cl,2)+size(ThisCross%Cl,3)==9) then
                                        do l=2,lmax
                                            TotCl(index)%Cl(l,C_E:C_B) = TotCl(index)%Cl(l,C_E:C_B) + Chat%Cl(l,C_E:C_B)*ClWeight(l)
                                        end do
                                        TotCount(:,index,3) = TotCount(:,index,3) + ClWeight
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    print *,'done Pcl'
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
                        if (TotCount(l,Pix,1)/=0) then
                            TotCl(Pix)%Cl(l,C_T) = TotCl(Pix)%Cl(l,C_T)/TotCount(l,Pix,1)
                        end if
                        if (TotCount(l,Pix,2)/=0) then
                            TotCl(Pix)%Cl(l,C_C) = TotCl(Pix)%Cl(l,C_C)/TotCount(l,Pix,2)
                        end if
                        if (TotCount(l,Pix,3)/=0) then
                            TotCl(Pix)%Cl(l,C_E:C_B) = TotCl(Pix)%Cl(l,C_E:C_B)/TotCount(l,Pix,3)
                        end if
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
    character(LEN=fname_len) :: fname
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
        call ReadYearMaps(C, fname, M, Detector, pol_maps)
        do year=1, nyears
            ix = ix+1
            call HealpixMap_Nullify(WMaps(ix))
            M(year)%TQU = M(year)%TQU*map_scale/mK
            if (pol_weights) then
                call MapMulPolWeight(M(year),WeightMaps(weight),WeightMapsPol(Weight),WMaps(ix))
            else
                call HealpixMapMulCut(M(year),WeightMaps(weight),WMaps(ix), 1)
            end if
        end do
        call HealpixMapArray_Free(M)
    end do

    print *, 'Getting cross powers'

    call HealpixMapSet2CrossPowers(H, WMaps, PCls, nmaps, lmax,.false.)
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
    character(LEN=fname_len) :: fname
    integer DA
    integer year, year2 , i
    integer (KIND=8) :: missing,pixcount
    Type(HealpixMap) :: Mask
    Type(HealpixMap) :: M(nyears)
    real(dp) :: sig0(nyears*(nyears+1)/2)
    real sigvar, hit1,hit2
    integer ix
    character(LEN=fname_len) :: cachefile

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
        call HealpixMap_ForceRing(Mask)

        fname = FormatFilename(year_filename_format, C%Name)
        print *,'Getting noise: '//trim(C%Name)
        do DA = 1, C%Count
            call ReadYearMaps(C,fname, M, DA)

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

    end subroutine EstAndSetNoise


    subroutine CombineYears(C)
    Type(TChannel) :: C
    character(LEN=fname_len) :: aname, outname
    integer year, Detector
    Type(HealpixMap) :: M, MTot
    Type(TBeam) :: Beam
    !Assume second column is hit count


    do Detector = 1,C%Count
        outname = FormatFilename(detector_filename_format, C%Name, DetectorName(C,Detector))
        if (.not. FileExists(outname)) then
            print *,'Combining years '//trim(C%Name), Detector
            if (.not. noise_from_hitcounts) then
                if (nyears>1) call MpiStop('generalise combined maps')
                aname = FormatFilename(year_filename_format, C%Name, DetectorName(C,Detector), Year)
                call HealPixMap_read(Mtot, aname, pol_maps)
                call HealpixMap_ForceRing(Mtot)
            else
                do year = 1, nyears
                    aname = FormatFilename(year_filename_format, C%Name, DetectorName(C,Detector), Year)

                    print *,'reading '//trim(aname)
                    if (year==1) then
                        call HealPixMap_read(Mtot, aname)
                        call HealpixMap_ForceRing(Mtot)
                        MTot%TQU(:,1) = MTot%TQU(:,1)*Mtot%TQU(:,2)
                    else
                        call HealPixMap_read(M, aname)
                        call HealpixMap_ForceRing(M)
                        MTot%TQU(:,1) = MTot%TQU(:,1) + M%TQU(:,1)*M%TQU(:,2)
                        MTot%TQU(:,2) = MTot%TQU(:,2) + M%TQU(:,2)
                    end if
                end do
                call HealpixMap_Free(M)
                where (MTot%TQU(:,2)>0)
                    MTot%TQU(:,1) = MTot%TQU(:,1)/MTot%TQU(:,2)
                end where
            end if

            print *,'writing '//trim(outname)
            call HealpixMap_Write(Mtot, outname)
            call HealpixMap_Free(MTot)
        end if
    end do

    end  subroutine CombineYears

    subroutine ProcessChannelDatamap(H, C, M)
    Type(HealpixInfo) :: H
    Type(TChannel), target :: C
    Type(HealpixMap) ::  M
    Type(HealpixMap) :: AMap
    character(LEN=fname_len) :: cache_name, map_fname, fname
    integer i, polcount

    Type (HealpixAlm) :: A
    Type(HealpixPOwer) :: Chat


    call HealpixMap_Nullify(M)
    if (detector_filename_format=='') call MpiStop('No Data file')
    map_fname = FormatFilename(detector_filename_format, C%Name)
    if (C%Count==1 .and. FileExists(map_fname)) then
        call HealpixMap_Read(M,map_fname, pol_maps)
    else
        cache_name = CacheName(map_fname, .true.)
        call StringReplace('%DA%','allDA_signal',cache_name)
        call StringReplace('.fits','_all.fits',cache_name)

        if (FileExists(cache_name)) then
            print *,'reading cached data map: '// trim(C%Name)
            call HealpixMap_Read(M,cache_name)
        else
            call healpixMap_Init(M, npix, pol = want_pol)

            call HealpixMap_ForceRing(M)
            polcount=0
            do i=1, C%Count
                fname = FormatFilename(detector_filename_format, C%Name,DetectorName(C,i))
                call HealpixMap_Read(AMap, fname, pol_maps)

                call HealpixMap_ForceRing(AMap)
                where (M%TQU(:,1) /= fmissval .and. AMap%TQU(:,1) /= fmissval)
                    M%TQU(:,1) = M%TQU(:,1) + AMap%TQU(:,1)
                elsewhere
                    M%TQU(:,1) = fmissval
                end where

                if (want_pol .and. C%detector_has_pol(i)) then
                    polcount=polcount+1
                    where (M%TQU(:,2:3) /= fmissval .and. AMap%TQU(:,2:3) /= fmissval)
                        M%TQU(:,2:3) = M%TQU(:,2:3) + AMap%TQU(:,2:3)
                    elsewhere
                        M%TQU(:,2:3) = fmissval
                    end where
                end if
            end do
            where (M%TQU(:,1) /= fmissval)
                M%TQU(:,1) = M%TQU(:,1) /C%Count
            end where
            if (want_pol) then
                where (M%TQU(:,2:3) /= fmissval)
                    M%TQU(:,2:3) = M%TQU(:,2:3) /PolCount
                end where
            end if
            call HealpixMap_Free(AMap)
            print *,'writing '//trim(cache_name)
            call HealpixMap_Write(M,cache_name)
        end if
    end if
    where (M%TQU /= fmissval)
        M%TQU = M%TQU * ( map_scale/mK )
    end where
    print *, 'data map min/max ', minval(M%TQU(:,1)), maxval(M%TQU(:,1))

    end subroutine ProcessChannelDatamap

    subroutine GetSmoothedNoise(H)
    Type(HealpixInfo) :: H
    Type(HealpixMap) :: Hits
    character(LEN=fname_len):: cache_name
    real MaxN, minnoise

    call HealpixMap_Nullify(SmoothedNoise)

    if (.not. fullsky_test) then
        if (noise_inv_fwhm/=0.d0) then
            cache_name = CacheName('smoothed_noise',.false.)

            if (.not. no_cache .and. FileExists(cache_name)) then
                print *,'Reading smoothed noise'
                call HealpixMap_Read(SmoothedNoise,cache_name)
            else
                print *,'Getting smoothed noise map'
                if (noise_map_for_window/='') then
                    call HealpixMap_Read(Hits,noise_map_for_window)
                    if (noise_from_hitcounts) then
                        call HealpixMap_SetToIndexOnly(Hits, min(Hits%nmaps,2))
                        Hits%TQU = 1/Hits%TQU
                    else
                        call HealpixMap_SetToIndexOnly(Hits, 1) !use TT noise for smoothed window
                        MaxN = maxval(Hits%TQU)
                        where (Hits%TQU <0)
                            Hits%TQU = MaxN
                        end where
                    end if
                    minnoise = minval( Hits%TQU(:,1) )
                    Hits%TQU(:,1) = Hits%TQU(:,1)+minnoise  !don't need this for WMAP

                    call HealpixMap_udgrade(Hits, SmoothedNoise, nside, pessimistic=.false.)
                    call HealpixMap_Smooth(H, SmoothedNoise, SmoothedNoise, 2*lmax, noise_inv_fwhm)
                    call DeleteFile(cache_name)
                    call HealpixMap_Write(SmoothedNoise,cache_name)
                end if
            end if
            print *,'min/max smoothed noise = ', minval(SmoothedNoise%TQU(:,1))/maxval(SmoothedNoise%TQU(:,1))
        else
            print *, 'using non-smoothed noise'
            call HealpixMap_Read(SmoothedNoise,noise_map_for_window)
            if (noise_from_hitcounts) then
                call healpixMap_SetToIndexOnly(SmoothedNoise,Min(smoothedNoise%nmaps,2))
                where (SmoothedNoise%TQU >0)
                    SmoothedNoise%TQU = 1/SmoothedNoise%TQU
                end where
            else
                call healpixMap_SetToIndexOnly(SmoothedNoise,1)
            end if
            call HealpixMap_ForceRing(SmoothedNoise)
        end if
    end if


    end subroutine GetSmoothedNoise

    subroutine ProcessNoiseMaps(H, C)
    Type(HealpixInfo) :: H
    Type(TChannel), target :: C
    Type(HealpixMap), pointer ::  NoiseMap
    Type(HealpixMap) :: AMap
    character(LEN=fname_len) :: cache_name, noise_fname, fname
    integer minnoise
    integer i, year,npol, polcount

    print *, 'ProcessNoiseMaps'
    NoiseMap =>  C%NoiseMap

    call HealpixMap_Nullify(NoiseMap)
    call HealpixMap_Nullify(AMap)

    if (detector_noise_filename_format/='') then
        noise_fname = FormatFilename(detector_noise_filename_format, C%Name)
        cache_name = CacheName(noise_fname, .true.)
        call StringReplace('%DA%','allDA',cache_name)

        if (FileExists(cache_name) .and. .not. cross_spectra) then
            print *,'reading cached noise map: '// trim(C%Name)
            call HealpixMap_Read(NoiseMap,cache_name, pol_maps)
        else
            if (cross_spectra) then
                allocate(C%DetectorYearNoiseMaps(C%Count,nyears))
            end if
            polcount=0
            do i=1, C%Count
                fname = FormatFilename(detector_noise_filename_format, C%Name, DetectorName(C,i))
                call GetNoiseMap(H, AMap, fname, C%sig0(i))
                call HealpixMap_ForceRing(AMap)
                if (i==1) then
                    call HealpixMap_Init(NoiseMap, AMap%npix, pol_maps)
                end if
                npol=1
                if (want_pol .and. C%detector_has_pol(i) ) then
                    npol=3
                    polcount=polcount+1
                end if
                where (NoiseMap%TQU(:,1:npol) /= fmissval .and. AMap%TQU(:,1:npol) /= fmissval)
                    NoiseMap%TQU(:,1:npol) =  NoiseMap%TQU(:,1:npol) + AMap%TQU(:,1:npol)
                elsewhere
                    NoiseMap%TQU(:,1:npol) = fmissval
                end where

                if (cross_spectra) then
                    if (nyears==1) then
                        call HealpixMap_Nullify(C%DetectorYearNoiseMaps(i,1))
                        print *,'assigning 1 year'
                        call HealpixMap_Assign(C%DetectorYearNoiseMaps(i,1),AMap)
                    else
                        do year = 1, nyears
                            call HealpixMap_Nullify(C%DetectorYearNoiseMaps(i,year))
                            fname = FormatFilename(year_noise_filename_format, C%Name, DetectorName(C,i), year)
                            call GetNoiseMap(H, C%DetectorYearNoiseMaps(i,year), fname, C%sig0(i))
                        end do
                    end if
                end if
            end do

            where (NoiseMap%TQU(:,1) /= fmissval)
                NoiseMap%TQU(:,1) = NoiseMap%TQU(:,1)/C%Count**2
            end where
            if (want_pol) then
                where (NoiseMap%TQU(:,2:3) /= fmissval)
                    NoiseMap%TQU(:,2:3) = NoiseMap%TQU(:,2:3)/polcount**2
                end where
            end if

            call HealpixMap_Free(AMap)

            if (.not. cross_spectra) then
                print *,'writing '//trim(cache_name)
                call HealpixMap_Write(NoiseMap,cache_name)
            end if
        end if
        call HealpixMap_ForceRing(NoiseMap)

        print *,trim(C%Name)//' min/max TT noise = ', &
        minval(NoiseMap%TQU(:,1), mask = NoiseMap%TQU(:,1) /= fmissval), &
        maxval(NoiseMap%TQU(:,1), mask = NoiseMap%TQU(:,1) /= fmissval)
        if (want_pol) then
            print *,trim(C%Name)//' min/max Pol noise = ', &
            minval(NoiseMap%TQU(:,2:3), mask = NoiseMap%TQU(:,2:3) /= fmissval), &
            maxval(NoiseMap%TQU(:,2:3), mask = NoiseMap%TQU(:,2:3)/= fmissval)
            C%ENoiseFac = sum(dble(NoiseMap%TQU(:,2)+NoiseMap%TQU(:,3))/dble(NoiseMap%TQU(:,1)), &
            mask = NoiseMap%TQU(:,2)/= fmissval .and. NoiseMap%TQU(:,1)/= fmissval )/ &
            (2* count(NoiseMap%TQU(:,2)/= fmissval .and. NoiseMap%TQU(:,1)/= fmissval))
            ENoiseFac = C%ENoiseFac
            print *, 'Empirical ENoiseFac = ', ENoiseFac
        end if
    else
        call HealpixMap_Init(NoiseMap, npix,pol = want_pol)
        NoiseMap%TQU(:,1) = noise*NoiseMap%npix/(HO_fourpi)
        NoiseMap%TQU(:,2:3) = ENoiseFac*noise*NoiseMap%npix/(HO_fourpi)
    end if

    if (uniform_noise) then
        !Isotropic white noise, no cut
        print *,'Doing uniform white noise'
        NoiseMap%TQU(:,1) = 1/(sum(1/dble(NoiseMap%TQU(:,1)), mask=NoiseMap%TQU(:,1)>0)/count(NoiseMap%TQU(:,1)>0))
        if (want_pol) &
        NoiseMap%TQU(:,2:3) = 1/(sum(1/dble(NoiseMap%TQU(:,2:3)), mask=NoiseMap%TQU(:,2:3)>0)/count(NoiseMap%TQU(:,2:3)>0))

        white_NL = NoiseMap%TQU(1,1)*HO_fourpi/NoiseMap%npix
        white_NL_P = ENoiseFac*white_NL
    else
        white_NL =  1/(sum(1/dble(NoiseMap%TQU(:,1)), mask=NoiseMap%TQU(:,1)>0)/count(NoiseMap%TQU(:,1)>0)) &
        *HO_fourpi/NoiseMap%npix
        if (want_pol) white_NL_P =  1/(sum(1/dble(NoiseMap%TQU(:,2:3)), &
        mask=NoiseMap%TQU(:,2:3)>0)/count(NoiseMap%TQU(:,2:3)>0))  *HO_fourpi/NoiseMap%npix
    end if

    end subroutine ProcessNoiseMaps

    subroutine GetNoiseMap(H, NoiseMap, noise_map, noise_sig0)
    character(LEN=*), intent(in) :: noise_map
    Type(HealpixInfo) :: H
    Type(HealpixMap) :: NoiseMap
    Type(HealpixMap) :: Hits
    real(dp) :: meanhit, booknoise
    real(sp) :: minnoise
    real(dp), intent(in) :: noise_sig0
    integer i
    print *,'converting noise: '  // trim(noise_map)

    call HealpixMap_Nullify(NoiseMap)
    call HealpixMap_Read(NoiseMap, noise_map)
    if (noise_from_hitcounts) then
        if (want_pol) call MpiStop('noise_from_hitcounts only TT at the mo')

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
    else
        if (.not. want_pol) call HealpixMap_SetToIndexOnly(NoiseMap,1)
        print *,'scaling noise ',(map_scale/mK)**2
        where (NoiseMap%TQU/=fmissval)
            NoiseMap%TQU = NoiseMap%TQU * (map_scale/mK)**2
        end where
    end if
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



    subroutine TestExactLike(H,WeightMap, WeightMapPol,NoiseMap, M, fid_cl_file, sim_cl_file, Coupler)
    use CutSkyAsymm
    Type(TCouplingMatrix) :: coupler(:)
    Type(HealpixInfo) :: H
    Type (TSNModeMat) :: SN
    Type(TSNModes) :: Modes
    Type (TSNModeOptions) :: Opts
    Type(HealpixMap) :: ProjMap, WeightMap, WeightMapPol, NoiseMap, M

    character(LEN=*), intent(in) :: fid_cl_file,sim_cl_file
    Type(HealpixPower) :: PFid, P, ProjCl, CHat
    Type(HealpixPower):: CheckCls1, CheckCls2
    integer l,l_low, nmodes
    integer  npol
    Type(AsymmCouplings) TestW
    real(sp), allocatable ::AlmVec(:), AlmVec2(:)
    real(dp) chisq,amp, like, term
    complex(dp) AMode
    character(LEN=fname_len) :: exact_file, exact_stem
    logical :: noise_cut = .false.
    Type(HealpixCrossPowers) :: PCls
    Type(HealpixMap) :: tmpM,WM, MapArr(3)
    Type(AsymmCouplings) :: Asymm
    Type (HealpixAlm) :: MaskA, MaskAP,MapA, MapAProj
    Type(TCovMat) :: FiducialChol
    real(dp), allocatable :: BigModes(:)
    integer j

    real (sp), allocatable :: vec(:)
    asymm_pol = want_pol

    l_low = l_exact+l_exact_margin
    npol = 1
    if (want_pol) npol=3

    print *,'Doing exact with lmax = ',l_low
    exact_file = concat(trim(data_dir)//'Proj_lexact',l_exact,'_margin',l_exact_margin,'_supp'//trim(WellSupported_txt))
    if (want_pol) exact_file=concat(exact_file,'_pol')
    if (uniform_noise) exact_file=concat(exact_file,'_uninoise')
    exact_file = concat(exact_file, ExtractFileName(fits_mask))

    exact_stem = trim(concat('outfiles/'//trim(exact_tag),l_exact,'_llow',l_low,'_supp')) &
    //trim(WellSupported_txt)//'_SNcut'//trim(exact_SN_cut_txt)//'_'

    exact_file = concat(exact_file,'.dat')

    !    call ppm_masked_map(M,WeightMap, 'outfiles/weighted_map.ppm')

    if (test_big_cut) then
        !Test previous with bigger mask
        call HealpixMap_Read(ProjMap,concat(exact_stem,'lowl_map.fits'))

        call HealpixVis_MapMask2ppm(ProjMap,WeightMap, concat(exact_stem,'lowl_map_bigcut_masked.ppm'), 200.)

        ProjMap%TQU = ProjMap%TQU * WeightMap%TQU

        call HealpixMap2Alm(H,ProjMap,MapAProj,Coupler(1)%lmax)
        call HealpixAlm2Power(MapAProj,ProjCl)
        call HealpixPower_Write(ProjCl,concat(exact_stem,'lowl_map_bigcut_pcl.dat'))

        call HealpixPower_Nullify(CHat)
        call HealpixMap_Nullify(MapArr(1))
        call HealpixMap_Assign(MapArr(1), ProjMap)
        call HealpixMap_Read(ProjMap,concat(exact_stem,'highl_map.fits'))
        ProjMap%TQU = ProjMap%TQU * WeightMap%TQU
        call HealpixMap_Nullify(MapArr(2))
        call HealpixMap_Assign(MapArr(2), ProjMap)
        call HealpixMapSet2CrossPowers(H, MapArr, PCls, 2, Coupler(1)%lmax,.false.)

        call PseudoCl_GetCHat(Coupler(1),PCls%Ps(1,1), Chat)
        call HealpixPower_Write(Chat,concat(exact_stem,'lowl_bigcut_hatcl.dat'))
        call PseudoCl_GetCHat(Coupler(1),PCls%Ps(2,2), Chat)
        call HealpixPower_Write(Chat,concat(exact_stem,'highl_bigcut_hatcl.dat'))
        call MpiStop()
    end if

    call HealpixPower_ReadFromTextFile(PFid,fid_cl_file,l_low,pol=want_pol,dolens = .false.)
    if (zero_from_l/=0) PFId%Cl(zero_from_l:lmax,:)=0
    call HealpixPower_Smooth(PFid,low_l_fwhm)

    !$ call OMP_SET_NUM_THREADS(4)

    if (.not. no_cache .and. FileExists(exact_file)) then
        print *, 'Reading cached projection matrix and modes'
        call TSNModeMat_Read(SN, exact_file, .true.)
        if (SN%MakeOptions%l_exact/=l_exact) call MpiStop('wrong l_exact in cache')
        if (SN%MakeOptions%l_low/=l_low) call MpiStop('wrong l_low in cache')
        Opts = SN%MakeOptions
        print *,'Read file'
    else
        print *, 'frac weight < 1e-5', count(WeightMap%TQU(:,1)<1e-5)/real(WeightMap%npix)
        print *, 'frac weight < 0.1', count(WeightMap%TQU(:,1)<0.1)/real(WeightMap%npix)

        Opts%l_exact = l_exact
        Opts%l_low = l_low
        Opts%SN_cut = exact_SN_cut
        Opts%want_pol = want_pol
        Opts%pol_weights = pol_weights
        Opts%want_cov = .true.
        Opts%ProjectMonoAndDipole = project_l01
        Opts%KeepSupport = 0.99d0
        Opts%ENoiseFac = 4
        Opts%keep_pseudo_noise_cov = do_lowl_separation
        Opts%fake_noise = fake_noise
        if (exact_cl_lmin  <2) then
            Opts%l_min_exact = exact_cl_lmin
            PFid%Cl(exact_cl_lmin:1,1) = PFid%Cl(2,1) !Assume looking for residuals of order the quadrupole
        else
            Opts%l_min_exact = 2
        end if
        call SNModes_GetSNModeMats(H, SN, Opts, WeightMap, WeightMapPol, NoiseMap, PFid)
        call TSNModeMat_Write(SN, exact_file)
    end if

    print *,'udgrade'
    call HealpixMap_udgrade(M,WM , nside, pessimistic=.false.)
    call HealpixMap_Assign(M, WM)
    print *,'mult weight'
    WM%TQU(:,1)=M%TQU(:,1)*WeightMap%TQU(:,1)
    if (want_pol) then
        if (pol_weights) then
            WM%TQU(:,2)=M%TQU(:,2)*WeightMapPol%TQU(:,1)
            WM%TQU(:,3)=M%TQU(:,3)*WeightMapPol%TQU(:,1)
        else
            WM%TQU(:,2)=M%TQU(:,2)*WeightMap%TQU(:,1)
            WM%TQU(:,3)=M%TQU(:,3)*WeightMap%TQU(:,1)
        end if
    end if

    print *, 'Get SN modes for map'
    call SNModes_GetSNModes(H, SN, SN%MakeOptions,WM, Modes)
    call HealpixMap_Free(WM)

    if (get_exact_cl) then
        !!!!!!
        print *,'reading ILC'
        call HealpixMap_Read(tmpM,concat(input_data_dir,'ILC.fits'))
        call HealpixMap_ForceRing(tmpM)
        call HealpixVis_Map2ppmfile(tmpM, concat(exact_stem,'Planck.ppm'),symmetric=.true.)

        print *,'diff map'
        M%TQU(:,1) = tmpM%TQU(:,1) - M%TQU(:,1)/1000
        print *,'diff pic'
        call HealpixVis_Map2ppmfile(M, concat(exact_stem,'WMAP_Planck_diff.ppm'),symmetric=.true.)
        print *,'smooth map'
        call HealpixMap_Smooth(H, M, M, 40, 50.d0/60)
        call HealpixVis_Map2ppmfile(M, concat(exact_stem,'WMAP_Planck_diff_smooth.ppm'),symmetric=.true.)
        call MpiStop('done')

        print *,'Doing exact cl with lmin = ', exact_cl_lmin, 'lmax=',exact_cl_lmax
        call GetFullSkyAlmEstimator(SN,Modes, PFid,MapAProj, exact_cl_lmin,exact_cl_lmax) !get quad and octopole

        call HealpixAlm2Map(H, MapAProj, ProjMap, M%npix)
        call HealpixVis_Map2ppmfile(ProjMap, concat(exact_stem,'exact_alm_map.ppm'),symmetric=.true.)
        call HealpixVis_MapMask2ppm(ProjMap,WeightMap, concat(exact_stem,'exact_alm_map_masked.ppm'), 200.)
        call HealpixMap2Alm(H, ProjMap, MapAProj, exact_cl_lmax + 2)
        call HealpixAlm2Power(MapAProj,ProjCl)
        call HealpixPower_Write(ProjCl, concat(exact_stem,'exact_cl.dat'))

        call mpistop()
    end if

    if (do_lowl_separation) then
        if (want_pol) stop 'not done pol high l projection'
        !Get high-l map with rubbish projected

        print *,'getting lowl map'
        !Eq 20 in paper

        !Theory part
        print *,'theory part of cov'
        allocate(AlmVec(SN%TheoryProj%nl))

        print *, SN%TheoryProj%nl,SN%TheoryProj%nr, Modes%nmodes , SN%DataProj%nr

        AlmVec =  matmul(SN%TheoryProj%M,Modes%SNModes(1)%V)
        print *,' mult by C'
        call HealpixVectorMultPower(AlmVec,PFid, l_low)
        call HealpixMap2Alm(H,WeightMap, MaskA, l_low*2)
        if (pol_weights) then
            call HealpixMap2Alm(H,WeightMapPol, MaskAP, l_low*2)
        else
            MaskAP = MaskA
        end if
        print *,' Get pseudo alm coupling matrix'
        call CutSkyAsymm_GetCoupling(TestW, MaskA, MaskAP, l_low, plusonly = .false.)
        AlmVec = matmul(TestW%WAsymm, AlmVec)

        deallocate(TestW%WAsymm)

        !Noise part
        print *,'noise part of cov'
        allocate(AlmVec2(SN%DataProj%nl))
        AlmVec2 =  matmul(SN%DataProj%M,Modes%SNModes(1)%V)
        print *, sum(abs(AlmVec2))
        AlmVec2 =  matmul(SN%PseudoNoiseCov%WAsymm, AlmVec2)
        print *, sum(abs(AlmVec2))
        AlmVec = AlmVec + AlmVec2
        deallocate(ALmVec2)
        print *, SN%DataProj%nl, SN%TheoryProj%nl, (l_low+1)**2, SN%TheoryProj%nr, SN%DataProj%nr
        call HealpixVector2Alm(AlmVec, MapAProj, l_low, polix=1)

        call HealpixAlm2Power(MapAProj,ProjCl)

        call HealpixPower_Write(ProjCl,concat(exact_stem,'lowl_map_cl.dat'))

        call HealpixMap_Nullify(ProjMap)
        call HealpixAlm2Map(H, MapAProj, ProjMap, M%npix)

        call DeleteFile(concat(exact_stem,'lowl_map.fits'))
        !     call HealpixMap_Write(ProjMap,concat(exact_stem,'lowl_map.fits'))

        call HealpixVis_Map2ppmfile(ProjMap, concat(exact_stem,'lowl_map.ppm'),symmetric=.true.)
        call HealpixVis_MapMask2ppm(ProjMap,WeightMap, concat(exact_stem,'lowl_map_masked.ppm'), 200.)
        call HealpixMap_Nullify(MapArr(1))
        call HealpixMap_Nullify(MapArr(2))
        call HealpixMap_Nullify(MapArr(3))

        call HealpixMap_Assign(MapArr(1), ProjMap)

        ProjMap%TQU = M%TQU - ProjMap%TQU
        call DeleteFile(concat(exact_stem,'highl_map.fits'))
        !     call HealpixMap_Write(ProjMap,concat(exact_stem,'highl_map.fits'))

        call HealpixVis_MapMask2ppm(ProjMap,WeightMap, concat(exact_stem,'highl_map.ppm'))
        call HealpixMap_Assign(MapArr(2), ProjMap)
        MapArr(2)%TQU = MapArr(2)%TQU *WeightMap%TQU

        call HealpixMap_Assign(MapArr(3), M)
        MapArr(3)%TQU = MapArr(3)%TQU *WeightMap%TQU


        MapAProj%TEB(:,0:l_exact,:) = 0
        call HealpixAlm2Map(H, MapAProj, ProjMap, M%npix)
        call HealpixVis_Map2ppmfile(ProjMap, concat(exact_stem,'lowl_map_highl.ppm'),symmetric=.true.)
        call HealpixVis_MapMask2ppm(ProjMap,WeightMap, concat(exact_stem,'lowl_map_highl_masked.ppm'),200.)

        call HealpixVector2Alm(AlmVec, MapAProj, l_low, polix=1)

        call HealpixMap2Alm(H,M,MapA,l_low)
        MapAProj%TEB = MapA%TEB - MapAProj%TEB
        MapAProj%TEB(:,l_exact+1:l_low,:) = 0
        call HealpixAlm2Map(H, MapAProj, ProjMap, M%npix)
        call HealpixVis_MapMask2ppm(ProjMap,WeightMap, concat(exact_stem,'lowl_missing_masked.ppm'),200.)

        call HealpixAlm2Power(MapAProj,ProjCl)
        call HealpixPower_Write(ProjCl,concat(exact_stem,'lowl_missing_pcl.dat'))

        call HealpixAlm_Free(MapAProj)
        call HealpixAlm2Power(MapA,ProjCl)
        call HealpixPower_Write(ProjCl,concat(exact_stem,'all_pcl.dat'))

        print *,'getting Chat'
        call HealpixPower_Nullify(CHat)
        call HealpixMapSet2CrossPowers(H, MapArr, PCls, 3, Coupler(1)%lmax,.false.)

        call PseudoCl_GetCHat(Coupler(1),PCls%Ps(1,1), Chat)
        call HealpixPower_Write(Chat,concat(exact_stem,'lowl_hatcl.dat'))
        call PseudoCl_GetCHat(Coupler(1),PCls%Ps(2,2), Chat)
        call HealpixPower_Write(Chat,concat(exact_stem,'highl_hatcl.dat'))
        call PseudoCl_GetCHat(Coupler(1),PCls%Ps(3,3), Chat)
        call HealpixPower_Write(Chat,concat(exact_stem,'all_cut_hatcl.dat'))
        print *, 'done Chat'

        M%TQU = M%TQU * WeightMap%TQU
        call HealpixMap2Alm(H,M, MapA, l_low)
        call HealpixAlm2Power(MapA,ProjCl)
        call HealpixPower_Write(ProjCl,concat(exact_stem,'all_cut_pcl.dat'))

        deallocate(AlmVec)
        deallocate(SN%PseudoNoiseCov%WASymm)
    end if

    !Chop to just range of interest
    SN%TheoryProj%nl =(l_exact+1)**2
    if (want_pol) SN%TheoryProjPol%nl = (l_exact+1)**2-4

    nmodes = SN%TheoryProj%nr
    if (want_pol) then
        nmodes = nmodes + SN%TheoryProjPol%nr*2
    end if

    print *,'testing likelihood'
    !   if (check_cls_file1 /='') then
    !    call HealpixPower_ReadFromTextFile(checkCls1,check_cls_file1,l_low,pol=want_pol,dolens = .false.)
    !    call HealpixPower_ReadFromTextFile(checkCls2,check_cls_file2,l_low,pol=want_pol,dolens = .false.)
    !    call HealpixPower_ReadFromTextFile(PFid,check_cls_file1,l_low,pol=want_pol,dolens = .false.) !dummy allocate
    !   end if
    !
    !   allocate(BigModes(nmodes))
    !
    !   if (get_mean_likes) then
    !        print *,'Getting mean log likelihoods'
    !        call HealpixPower_ReadFromTextFile(PFid,sim_cl_file,l_low,pol=want_pol,dolens = .false.)
    !       !!!
    !    !    PFid%Cl(2:lmax,C_B)=0
    !        !For full sky result
    !        P%Cl = PFid%Cl
    !        P%Cl(2:l_low,C_T) = P%Cl(2:l_low,C_T) + white_NL
    !        if (want_pol) then
    !         P%Cl(2:l_low,C_E) = P%Cl(2:l_low,C_E) + white_NL_P
    !         P%Cl(2:l_low,C_B) = P%Cl(2:l_low,C_B) + white_NL_P
    !        end if
    !        !For modes
    !        StTime = GeteTime()
    !        call CutSkyAsymm_GetFullCovariance(TheoryProj,TheoryProjPol, FiducialChol, PFid, 2,l_exact)
    !        print *, 'full covariance time', GeteTime()  - StTime
    !        do j=1,nmodes
    !         FiducialChol%C(:,j) = FiducialChol%C(:,j)+SN%BigNoiseCov%C(:,j) + SN%OtherLNoiseCov%C(:,j)
    !        end do
    !        call Matrix_Cholesky(FiducialChol%C)  !Note upper triangular is not zeroed
    !        print *,'got mean matrix'
    !   end if
    !
    !   do i=0, 20
    !
    !    if (check_cls_file1 /='') then
    !     !Check interpolation between two models
    !      PFid%Cl = CheckCls2%Cl*(i/20.) + (1-i/20.)*CheckCls1%Cl
    !    else
    !     call HealpixPower_ReadFromTextFile(PFid,sim_cl_file,l_low,pol=want_pol,dolens = .false.)
    !      ! amp = 1+(i-8)/real(l_exact)/2
    !       amp = i/10.
    !
    !      if (amp<0) cycle
    !!      PFid%Cl(2:l_exact,C_B) =  PFid%Cl(2:l_exact,C_B)*amp
    !
    !   !   amp = 1+(i-8)/real(l_exact)/2.
    !      PFid%Cl(2:lmax,C_B) =  PFid%Cl(2:lmax,C_B)*amp
    !     end if
    !
    !    call CutSkyAsymm_GetFullCovariance(TheoryProj,TheoryProjPol, BigCov, PFid, 2,l_exact)
    !
    !     if (get_mean_likes) then
    !         like= 0
    !        call  Matrix_CholeskyRootInverse(BigCov%C)
    !        do j=1, nmodes
    !         like = like - log(BigCov%C(j,j))
    !       end do
    !       call Matrix_MultTri(BigCov%C,FiducialChol%C,'Right')
    !       like = like + sum(BigCov%C**2)/2
    !     else
    !        do j=1, TheoryProj%nr
    !        BigModes(j) = SNModes(1)%V(j)
    !        end do
    !        if (want_pol) then
    !        do j=1, TheoryProjPol%nr
    !        BigModes(TheoryProj%nr+j) = SNModes(2)%V(j)
    !        BigModes(TheoryProj%nr+TheoryProjPol%nr+j) = SNModes(3)%V(j)
    !        end do
    !        end if
    !!         print *,'testsum', sum(BigCov%C**2)
    !        like = Matrix_GaussianLogLike(BigCov%C,BigModes)
    !     end if
    !
    !    deallocate(BigCov%C)
    !    chisq = 0
    !    do l=2, l_exact
    !     PFid%Cl(l,C_T) = PFid%Cl(l,C_T) + white_NL
    !      if (want_pol) then
    !       PFid%Cl(l,C_E) = PFid%Cl(l,C_E) + white_NL_P
    !       PFid%Cl(l,C_B) = PFid%Cl(l,C_B) + white_NL_P
    !
    !       term = PFid%Cl(l,C_T)*PFid%Cl(l,C_E) - PFid%Cl(l,C_C)**2
    !
    !       chisq = chisq + (2*l+1)* ( &
    !       (PFid%Cl(l,C_T)*P%Cl(l,C_E) + PFid%Cl(l,C_E)*P%Cl(l,C_T) - 2 *PFid%Cl(l,C_C)*P%Cl(l,C_C))/term &
    !        + log( term/ (P%Cl(l,C_T)*P%Cl(l,C_E) - P%Cl(l,C_C)**2)) -2)
    !
    !        chisq = chisq + (2*l+1)* (P%Cl(l,C_B)/PFid%Cl(l,C_B) &
    !                +log(PFid%Cl(l,C_B)/P%Cl(l,C_B)) - 1)
    !
    !      else
    !       chisq = chisq + (2*l+1)* (P%Cl(l,C_T)/PFid%Cl(l,C_T) + log(PFid%Cl(l,C_T)))
    !      end if
    !    end do
    !    print *, i, chisq/2, like
    !
    !end do
    !
    !    print *, 'Time for likes', GeteTime()-StTime
    !

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
    character(LEN=fname_len) :: fname
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

    subroutine GetApodizedMask(H,Mask,cache_fname,fits_mask)
    Type(HealpixInfo) :: H
    Type(HealpixMap) :: Mask, BinaryMask
    character(LEN=*), intent(in) :: cache_fname, fits_mask
    Type(HealpixAlm) :: A
    integer process_mask, i

    if (.not. FileExists(cache_fname)) then
        print *,'Reading BinaryMask'
        call HealpixMap_Nullify(BinaryMask)
        call HealpixMap_Read(BinaryMask, fits_mask)

        if (BinaryMask%nside /= nside) print *, 'Upgrading '
        call HealpixMap_udgrade(BinaryMask, Mask, nside, pessimistic=.false.)
        call HealpixMap_Free(BinaryMask)
        call HealpixMap_ForceRing(Mask)

        print *,'generating apodised high-res mask'
        call HealpixMap2Alm(H,Mask, A, lmax*2, map_ix=min(2,Mask%nmaps)) !!!
        call HealpixAlm_Smooth(A,apodize_mask_fwhm)
        call HealpixAlm2Map(H,A,Mask,npix)

        print *,'minmax mask first pass  = ', minval(Mask%TQU), maxval(Mask%TQU)
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
            if (processing_mask_badpix) then
                where (BinaryMask%TQU(:,1)<0)
                    Mask%TQU(:,1)=0
                end where
            else
                if (processing_mask_map > BinaryMask%nmaps) call MpiStop('processing_mask_map > BinaryMask%nmaps')
                process_mask=0
                do i=0, Mask%npix
                    if (BinaryMask%TQU(i,processing_mask_map)==0)  then
                        if (Mask%TQU(i,1)>0.1) process_mask=process_mask+1
                        Mask%TQU(i,1)=0
                    end if
                end do
                print *,'Masking processing mask, additional %: ', process_mask/real(BinaryMask%npix)*100
            end if
            call HealpixMap_Free(BinaryMask)
        end if

        call HealpixVis_Map2ppmfile(Mask, 'outfiles/mask.ppm')

        call HealpixMap_Write(MAsk,cache_fname)
    else
        print *,'Reading cached mask'
        call HealpixMap_Nullify(Mask)
        call HealpixMap_Read(Mask,cache_fname)
    end if

    print *,'minmax mask = ', minval(Mask%TQU), maxval(Mask%TQU)
    print *,'fsky = ', sum(real(Mask%TQU(:,1),dp))/Mask%npix
    end subroutine GetApodizedMask

    subroutine DifferenceMaps
    Type(HealpixMap) :: M, M2
    real(dp) :: norm1, norm2

    call HealpixMap_read(M,'C:\tmp\ctp3\v3.1\100Ghz23.fits')
    call HealpixMap_read(M2,'C:\tmp\ctp3\v3.1\100Ghz14.fits')
    norm1 = sum(M%TQU(:,1), mask = M%TQU(:,1) /= fmissval)/ count( M%TQU(:,1) /= fmissval)
    norm2 = sum(M2%TQU(:,1), mask = M2%TQU(:,1) /= fmissval)/ count( M2%TQU(:,1) /= fmissval)
    print *,norm1, norm2
    where (M%TQU(:,1)/= fmissval)
        m%tqu(:,1) = m%tqu(:,1) - norm1
    elsewhere
        m%tqu(:,1) = 0
    end where
    where (M2%TQU(:,1)/= fmissval )
        m2%tqu(:,1) = m2%tqu(:,1) - norm2
    elsewhere
        m2%tqu(:,1) = 0
    end where

    call HealpixVis_Map2ppmfile(M, 'z:\M.ppm')
    call HealpixVis_Map2ppmfile(M2, 'z:\M2.ppm')


    where (M%TQU(:,1) /= 0 .and. M2%TQU(:,1) /= 0)
        M%TQU(:,1) = M%TQU(:,1) - M2%TQU(:,1)
    elsewhere
        M%TQU(:,1) = 0
    end where

    call HealpixVis_Map2ppmfile(M, 'z:\diffmap.ppm')

    stop

    end subroutine DifferenceMaps


    subroutine CombineSplitNoiseMaps(tag)
    Type(HealpixMap) :: M, MQ,MU,MN
    character(LEN=*), intent(in) :: tag

    call HealpixMap_Read(M,trim(tag)//'_noise_II.fits')
    call HealpixMap_Read(MQ,trim(tag)//'_QQ.fits')
    call HealpixMap_Read(MU,trim(tag)//'_UU.fits')
    call HealpixMap_ForceRIng(M)
    call HealpixMap_ForceRIng(MQ)
    call HealpixMap_ForceRIng(MU)
    call HealpixMap_Init(MN,M%npix, pol=.true.)
    MN%TQU(:,1)=M%TQU(:,1)
    MN%TQU(:,2)=MQ%TQU(:,1)
    MN%TQU(:,3)=MU%TQU(:,1)
    call HealpixMap_WRite(MN,trim(tag)//'_pol_noise.fits')

    call HealpixMap_Read(M,trim(tag)//'.fits')
    call HealpixMap_Read(MQ,trim(tag)//'_Q.fits')
    call HealpixMap_Read(MU,trim(tag)//'_U.fits')
    call HealpixMap_ForceRIng(M)
    call HealpixMap_ForceRIng(MQ)
    call HealpixMap_ForceRIng(MU)
    call HealpixMap_Init(MN,M%npix, pol=.true.)
    MN%TQU(:,1)=M%TQU(:,1)
    MN%TQU(:,2)=MQ%TQU(:,1)
    MN%TQU(:,3)=MU%TQU(:,1)
    call HealpixMap_WRite(MN,trim(tag)//'_pol.fits')


    call HealpixMap_Free(M)
    call HealpixMap_Free(MQ)
    call HealpixMap_Free(MU)
    call HealpixMap_Free(MN)

    end subroutine CombineSplitNoiseMaps

    subroutine AddPlanckSimMaps(H)
    use pix_tools
    Type(HealpixInfo)  :: H
    Type(HealpixAlm) :: A
    Type(HealpixPower):: P, PS
    real(dp) testvar, MaxN, MaxNP
    Type(HealpixMap) :: NoiseSim, Signal, Noise, SigNoiseSim
    character(LEN=fname_len) :: fname, detector_name
    integer i,detector, neighbours(8), nneigh, j, nav, nin
    real fake_noise_fac
    real pixav(9)

    print *,'adding maps'

    call HealpixMap_Nullify(NoiseSim)
    call HealpixMap_Nullify(Noise)
    call HealpixMap_Nullify(Signal)
    call HealpixMap_Nullify(SigNoiseSim)

    !        do i=1,9
    !        call HealpixMap_Read(NoiseSim,concat(input_data_dir,'mc_noise/ctp3.nmc.0000'//trim(IntToStr(i))//'.mm2048.fits'))
    !        NoiseSim%TQU=NoiseSim%TQU * map_scale/mK;
    !        call HealpixMap2Alm(H,NoiseSim, A, lmax, dopol =  .true.)
    !        call HealpixAlm2Power(A, P)
    !        if (i==1) then
    !          call HealpixPower_assign(PS,P)
    !        else
    !          PS%Cl = PS%Cl + P%Cl
    !        end if
    !        call HealpixPower_Write(P,'fullsky_noise_cls_'//trim(sim_no)//'.dat')
    !        call HealpixMap_Free(NoiseSim)
    !        call HealpixAlm_Free(A)
    !
    !        end do
    !        PS%Cl = PS%Cl / 9
    !        call HealpixPower_Write(PS,'fullsky_noise_avcls.dat')
    !
    !        return

    call InitRandom()

    ! call HealpixMap_Read(NoiseSim,concat(input_data_dir,'ctp3.nmc.00000.mm2048.fits'))
    ! call HealpixMap2Alm(H,NoiseSim, A, lmax,dopol =  .true.)
    ! call HealpixAlm2Power(A, P)
    ! call HealpixPower_Write(P,'0000mm_noise_cls.dat')
    ! call MpiStop()
    call HealpixMap_Read(Signal,trim(input_data_dir)//'mask.fits')
    call HealpixVis_Map2ppmfile(Signal, 'outfiles/mask.ppm')

    call HealpixMap_Read(Signal,trim(input_data_dir)//'hitmap.fits')
    call HealpixVis_Map2ppmfile(Signal, 'outfiles/hitmap.ppm')

    call HealpixMap_Read(Signal,trim(input_data_dir)//'noise_variance.fits')
    call HealpixVis_Map2ppmfile(Signal, 'outfiles/noise_variance.ppm')


    call HealpixMap_Read(Signal,trim(input_data_dir)//'signal.fits')
    call HealpixMap_ForceRing(Signal)
    call HealpixMap2Alm(H,Signal, A, lmax,dopol =  .false.)
    call HealpixAlm2Power(A, P)
    call HealpixPower_Write(P,'outfiles/input_cls.dat')

    !  call HealpixMap_Write(Noise,trim(input_data_dir)//'diag_noise_maps.fits', .true.)

    do detector =1,9
        call HealpixMap_Read(NoiseSim,concat(trim(input_data_dir)//'noise_realization_0',detector,'.fits'))
        !        call HealpixMap_ForceNest(NoiseSim)
        if (any(NoiseSim%TQU(:,1)==fmissval))  print * ,detector, 'has missing pix'
        call HealpixMap_ForceRing(NoiseSim)

        call HealpixMap_Assign(SigNoiseSim,Signal);
        SigNoiseSim%TQU(:,1) = Signal%TQU(:,1) + NoiseSim%TQU(:,1)

        !
        !where (NoiseSim%TQU(:,1) /= fmissval)
        ! SigNoiseSim%TQU(:,1) = Signal%TQU(:,1) + NoiseSim%TQU(:,1)
        !elsewhere
        ! SigNoiseSim%TQU(:,1) = Gaussian1()* sqrt(MaxN)
        ! NoiseSim%TQU(:,4) = MaxN
        !end where
        !where (NoiseSim%TQU(:,2) /= fmissval)
        ! SigNoiseSim%TQU(:,2) = Signal%TQU(:,2) + NoiseSim%TQU(:,2)
        !elsewhere
        ! SigNoiseSim%TQU(:,2) = Gaussian1()* sqrt(MaxNP)
        ! NoiseSim%TQU(:,7) = MaxNP
        !end where
        !where (NoiseSim%TQU(:,3) /= fmissval)
        ! SigNoiseSim%TQU(:,3) = Signal%TQU(:,3) + NoiseSim%TQU(:,3)
        !elsewhere
        ! SigNoiseSim%TQU(:,3) = Gaussian1()* sqrt(MaxNP)
        ! NoiseSim%TQU(:,9) = MaxNP
        !end where
        !
        fname = concat(trim(input_data_dir)//'tot_',detector,'.fits')
        call HealpixMap_Write(SigNoiseSim,fname, .true.)
    end do
    call HealpixMap_Free(NoiseSim)
    call HealpixMap_Free(Noise)
    call HealpixMap_Free(Signal)
    call HealpixMap_Free(SigNoiseSim)

    !
    !        call HealpixMap_Read(NoiseSim,concat(input_data_dir,'ctp3.nmc.00000.mm2048.fits'))
    !        print *,'n noise maps = ',NoiseSim%nmaps
    !        call HealpixMap_ForceRing(NoiseSim)
    !
    !       if (.true.) then
    !        NoiseSim%TQU = NoiseSim%TQU*map_scale/mK
    !        Noise%TQU = Noise%TQU*(map_scale/mK)**2
    !        print *, 'avg noise', sqrt(sum(Noise%TQU(:,1))/Noise%npix)
    !        testvar = sum( dble(NoiseSim%TQU(:,1))/sqrt(Noise%TQU(:,1)), mask = Noise%TQU(:,1)>0) / &
    !                     count(Noise%TQU(:,1)>0)
    !        print *, 'mean TT = ', testvar
    !        testvar = sum( dble(NoiseSim%TQU(:,1)), mask = Noise%TQU(:,1)>0) / &
    !                     count(Noise%TQU(:,1)>0)
    !        print *, 'mean TT = ', testvar
    !        NoiseSim%TQU(:,1) = NoiseSim%TQU(:,1) - testvar
    !        testvar = sum( dble(NoiseSim%TQU(:,1))**2/Noise%TQU(:,1), mask = Noise%TQU(:,1)>0) / &
    !                     count(Noise%TQU(:,1)>0)
    !        print *,'TT noise test:', testvar
    !        testvar = sum( dble(NoiseSim%TQU(:,2))**2/Noise%TQU(:,2), mask = Noise%TQU(:,2)>0) / &
    !                     count(Noise%TQU(:,2)>0)
    !        print *,'QQ noise test:', testvar
    !        testvar = sum( dble(NoiseSim%TQU(:,3))**2/Noise%TQU(:,3), mask = Noise%TQU(:,3)>0) / &
    !                     count(Noise%TQU(:,3)>0)
    !        print *,'UU noise test:', testvar
    !        print *, 'UU-QQ', sum( dble(Noise%TQU(:,3)-Noise%TQU(:,2))**2/Noise%TQU(:,2)**2, mask = Noise%TQU(:,2)>0) / &
    !                     count(Noise%TQU(:,2)>0)
    !        call MpiStop()
    !       end if

    ! call HealpixMap2Alm(H,SigNoiseSim, A, lmax*2, dopol =  .true.)
    !! call HealpixAlm2Power(A, P)
    !     call HealpixPower_Write(P,'test_noise_cls.dat')
    end subroutine AddPlanckSimMaps

    subroutine FillHoles(M)
    !Fill missing pixels with average of surrounding pixels
    use pix_tools
    Type(HealpixMap) :: M
    integer i
    integer neighbours(8), nneigh, j, nav
    real pixav(9)
    integer nin

    nin= M%nmaps

    call HealpixMap_ForceNest(M)
    do i=0, M%npix
        if (any(M%TQU(i,:)==fmissval)) then
            call neighbours_nest(M%nside,i,neighbours,nneigh)
            nav = 0
            pixav =0
            do j=1, nneigh
                if (M%TQU(j,1) /= fmissval) then
                    nav = nav+1
                    pixav(1:nin) = pixav(1:nin) + M%TQU(j,1:nin)
                endif
                if (nav >0) then
                    M%TQU(i,:) = pixav(1:nin) / nav
                    !unclear what to do with the noise variance, just average for now - doesn't matter much for cross spectra
                else
                    print *,'all near pixels missing  '
                end if
            end do
        end if
    end do
    call HealpixMap_ForceRing(M)

    end subroutine FillHoles


    subroutine SplineSmoothSampleArray(vec,n)
    use Cubic_Spline_GCV
    !Assume normalized to approximately 1
    Type (HealpixPower) :: P , P2
    integer n,l
    REAL(DP) vec(:), X(n), DF(n), Y(n), C(n-1,3), &
    var, Wk(7*(n + 2)), SE(n)
    INTEGER IC,JOB,IER


    IC = n-1
    Job = 1
    var = -1
    do l=1,n
        x(l) = l
        !   vec(l) = P%Cl(l+1,1)/P2%Cl(l+1,1)
        DF(l) = 1
    end do
    call CUBGCV(X,vec,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
    do l=1,n-1
        vec(l) = Y(l)
    end do

    !        call CreateTxtFile('z:\interp.txt',2)
    !        do l=2,lmax-1
    !           print *,l,Y(l-1)
    !           write (2,'(1I7, 2E15.5)') l,P%Cl(l,1)/P2%Cl(l,1), Y(l-1)
    !        end do
    !        print *,'err = ', IER
    !        pause
    !        stop

    end subroutine SplineSmoothSampleArray

    subroutine SplineFitCl(P, Template)
    Type(HealpixPower) :: P, Template
    integer i,nl
    real(dp) vec(P%lmax-1)

    print *,'spline fit cl'
    nl = P%lmax -1
    do i=C_T, C_C
        print *, 'i=', i
        if (i==C_C .and. .not. all(Template%Cl(2:P%lmax,i)>0)) then
            !not a variance
            cycle
        end if
        vec = P%Cl(2:P%lmax,i)/Template%Cl(2:P%lmax,i)
        call SplineSmoothSampleArray(vec, nl)
        P%Cl(2:P%lmax,i) = vec*Template%Cl(2:P%lmax,i)
        if (.not. P%pol) exit
    end do
    print *,'end spline fit cl'

    end subroutine SplineFitCl


    end module WeightMixing

    !   subroutine DoINit
    !    use MPIstuff
    !    integer i,ix
    !    call mpi_init_thread(MPI_THREAD_FUNNELED,ix,i)
    !     print *,'mpi_init_thread',MPI_THREAD_FUNNELED,ix,i

    !  end subroutine DoINit


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
    Type(HealpixMap)   :: M, CutM, Mask, MaskPol, GradPhi
    Type(HealpixMap), allocatable, target   :: MapArray(:)
    Type(HealpixPower) :: PBinnedDiag, PUnlensed, PFid, P, PSim, &
    CAvg, CSkew, CAvgTemp,CVar,CoffVar,Chat, HybridNoise, HybridP, DataCl, DiagCov
    Type(HealpixPower) :: SigNoiseCov, OffDiagCov,SignalCov, NoiseCov, NEff, fSkyEff
    Type(HealpixPower), dimension(:,:), allocatable :: MaskP
    Type(HealpixPower), dimension(:,:), allocatable :: CovP

    Type(HealpixAlm)   :: A, A2
    Type (TCouplingMatrix), target :: dummycoupler(1)

    Type(HealpixMap), pointer :: Map1,Map2, amap, NoiseMap
    Type(HealpixCrossPowers) :: CovPowers, MaskPowers, MaskPowersAll

    Type(HealpixMap) :: BinaryMask
    integer rand_seed
    character(LEN=32) :: weights(5)
    character(LEN=fname_len)  :: NuMStr,cache_name,mask_fname
    character(LEN=fname_len)  :: l_stem, cls_file, cls_unlensed_sim_file, cls_unlensed_sim_file_tensor, &
    fid_cls_file, out_file_base, out_file_root, sim_map_file, analysis_root, anastem
    character(LEN=fname_len)  :: beamfile,sim_stem, covstem, detectorNames, detectorPols
    character(LEN=fname_len) :: healpixloc, beam_dir
    integer :: pol_vec_size, i, and_seed
    logical :: err, debug_files

    real(dp) ::  meanN, A_PtSc
    Type(TCovMatPolArray), target  :: CovArr(3)
    Type(TCovMat), pointer :: ACov
    Type(TCovMat) HyCov(3)
    Type(TCovMat) :: LensedCov,ACovHybrid
    real(DP) noise_fac, chisq, mask_lowl_fwhm
    real(dp), allocatable :: tmp(:,:), tmp2(:,:),tmpcov(:,:), simcov(:,:)
    real(dp), allocatable :: bigVec(:), bigVec2(:)
    real(SP) fac
    integer l,ix,j
    integer ix1,ix2, chan1,chan2,x,y

    Type (HealpixPower) :: PointSourceP
    logical :: get_covariance, get_sim, get_analysis
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
    logical :: use_openmp = .true.
    integer mkl_threads, max_threads
    character(LEN=4) add_planck_sim
    real(dp) :: noise_adjustment = 1.d0
    integer, parameter :: smooth_w = 2
    real :: smooth_kernel(-smooth_w*2:smooth_w*2)
    integer status
    !$ integer, external :: omp_get_max_threads

#ifdef MPIPIX
    call mpi_init(i)
#endif
    call GetMpiStat(MpiID, MpiSize)

    !$ if (use_openmp .and. MPiID==0) then
    !$ max_threads = min(16,omp_get_max_threads())
    !$ mkl_threads = max_threads
    !$ print *, 'max_threads = ', max_threads
    !$ end if

    !$ call omp_set_nested(0)
    !$ call omp_set_dynamic(1)
    !$ call omp_set_num_threads(1)

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
        read(map_unit,*) map_scale
    end if
    fid_cls_file = Ini_read_String('fid_cls_file')

    out_file_base = Ini_Read_String('out_file_root')
    out_file_root = 'outfiles/' //trim(out_file_base)
    analysis_root = Ini_read_String('analysis_root',.false.)

    want_pol = Ini_Read_Logical('want_pol')
    if (want_pol) then
        pol_maps = 3
    else
        pol_maps = 1
    end if
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
    noise_colour_factor = Ini_read_Real('noise_colour_factor',0.0)

    beam_transfer = Ini_read_Logical ('beam_transfer')

    year_filename_format = trim(input_data_dir)//Ini_Read_String('year_filename_format')

    Ini_Fail_On_Not_Found = .false.

    data_var = Ini_Read_String('data_var')
    noise_adjustment = Ini_read_Double('noise_adjustment',1.d0)

    !combined_filename_format is computed maps for each detector added over years
    detector_filename_format = Ini_read_String('detector_filename_format')
    if (detector_filename_format /= '') then
        detector_filename_format = trim(input_data_dir)//detector_filename_format
    else
        detector_filename_format = trim(data_dir)//trim(ExtractFileName(year_filename_format))
        call StringReplace('%YEAR%',trim(IntToStr(nyears))//'years',detector_filename_format)
    end if

    detector_noise_filename_format = Ini_Read_String('detector_noise_filename_format')
    if (detector_noise_filename_format == '') then
        detector_noise_filename_format = detector_filename_format
    else
        detector_noise_filename_format= trim(input_data_dir)//trim(detector_noise_filename_format)
    end if
    year_noise_filename_format = Ini_Read_String('year_noise_filename_format')
    if (year_noise_filename_format /='') then
        year_noise_filename_format= trim(input_data_dir)//trim(year_noise_filename_format)
    end if
    beam_dir = Ini_Read_String('beam_dir')
    if (beam_dir=='') beam_dir = input_data_dir

    beam_filename_format = trim(beam_dir)//Ini_Read_String('beam_filename_format')
    cross_beams = Ini_read_Logical('cross_beams', cross_beams)
    beam_file_column = Ini_read_Int('beam_file_column',beam_file_column)

    noise_from_hitcounts = Ini_read_Logical('noise_from_hitcounts',.false.)
    nchannels = Ini_read_Int('nchannels')
    allocate(Channels(nchannels))
    do i = 1, nchannels
        Channels(i)%Name = Ini_Read_String_Array('channel_name',i)
        Channels(i)%Count = Ini_read_int_array('channel_count',i)
        DetectorNames = Ini_Read_String_Array('channel_detector_names',i)
        call TStringList_Init(Channels(i)%DetectorNames)
        if (DetectorNames /='') then
            call TStringList_SetFromString(Channels(i)%DetectorNames,DetectorNames)
        end if
        allocate(Channels(i)%detector_has_pol(Channels(i)%Count))
        detectorPols= Ini_Read_String_Array('channel_detector_polarized',i)
        if (want_pol .and. detectorPols/='') then
            read(detectorPols,*) Channels(i)%detector_has_pol
        else
            Channels(i)%detector_has_pol = want_pol
        end if
        Channels(i)%Ghz = Ini_read_Double_Array('channel_Ghz',i)
        allocate(Channels(i)%sig0(Channels(i)%Count))
        if (noise_from_hitcounts) then
            NumStr = Ini_read_String_Array('noise_sig0',i)
            read(NumStr,*) Channels(i)%sig0
        else
            Channels(i)%sig0 = 1
        end if
        Channels(i)%Beam%beam_transfer = beam_transfer
        if (.not. beam_transfer) then
            Channels(i)%Beam%fwhm = Ini_read_real('fwhm')/60  !Read in arcminute, internally in degrees
        end if
    end do

    Ini_Fail_On_Not_Found = .false.

    add_planck_sim = Ini_read_String('add_planck_sim')

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
    if (cross_Spectra) then
        unequal_times_only = Ini_Read_Logical('unequal_times_only',.false.)
        if (MpiID==0) print *,'doing only cross-spectra from different time periods'
    end if
    no_pix_window = Ini_read_Logical('no_pix_window', .false.)

    if (uniform_noise) then
        nweights = 1
        lmin_hybrid_calc = lmin
    end if
    Ini_Fail_On_Not_Found= .false.

    w8name = Ini_Read_String('w8dir')
    if (w8name=='') then
        call get_environment_variable('HEALPIX', healpixloc, status=status)
        if (status==0) then
            w8name = trim(healpixloc)//'/data/'
        else
            call mpiStop('w8dir not found')
        end if
    end if

    mkl_threads = Ini_Read_Int('mkl_threads',1)

    action = Ini_read_String('action')

    debug_files = Ini_read_Logical('debug_files',.false.)

    check_cls_file1 = Ini_read_String('check_cls_file1')
    check_cls_file2 = Ini_read_String('check_cls_file2')

    fits_mask = trim(input_data_dir)//trim(Ini_read_String('unapodized_mask'))
    if (want_pol) then
        fits_mask_pol = trim(input_data_dir)//trim(Ini_read_String('unapodized_mask_pol'))
        pol_weights = fits_mask /= fits_mask_pol
    end if

    processing_mask = trim(input_data_dir)//Ini_read_String('processing_mask')
    processing_mask_badpix = Ini_read_Logical('processing_mask_badpix', .false.)
    if (.not. processing_mask_badpix .and. processing_mask/='') &
    processing_mask_map = Ini_read_Int('processing_mask_map',1)

    subtract_monopole = Ini_Read_Logical('subtract_monopole',.false.)

    noise_map_for_window = Ini_Read_String('noise_map_for_window')
    if (noise_map_for_window /='') noise_map_for_window = trim(input_data_dir)//trim(noise_map_for_window )
    cls_file = Ini_Read_String('cls_file')
    cls_unlensed_sim_file = Ini_read_string('cls_unlensed_sim_file')
    cls_unlensed_sim_file_tensor =  Ini_read_string('cls_unlensed_sim_file_tensor')

    do_exact_like = Ini_Read_Logical('do_exact_like')
    want_cl = .true.

    if (do_exact_like) then
        want_cl  = Ini_read_Logical('do_cl')
        WellSupported_txt  = Ini_read_string('min_support')
        WellSupported  = Ini_read_real('min_support')
        if (MpiID==0) print *,'min_support = ', WellSupported
        exact_SN_cut = Ini_read_Real('exact_SN_cut')
        exact_SN_cut_txt = Ini_read_String('exact_SN_cut')
        if (MpiID==0) print *,'SN cut = ', exact_SN_cut
        test_big_cut = Ini_Read_Logical('test_big_cut')
        l_exact = Ini_read_Int('l_exact')
        l_exact_margin = Ini_read_Int('l_exact_margin')
        exact_tag = Ini_Read_String('exact_tag')
        get_exact_cl = Ini_read_logical('get_exact_cl')
        zero_from_l = Ini_read_Int('zero_from_l')
        project_l01 = Ini_read_logical('project_l01')
        if (get_exact_cl) then
            exact_cl_lmin = Ini_read_int('exact_cl_lmin')
            exact_cl_lmax = Ini_read_int('exact_cl_lmax')
        end if
        if (exact_tag/='') exact_tag = trim(exact_tag)//'_'
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
    fits_mask_pol = FormatFilename(fits_mask_pol)

    processing_mask = FormatFilename(processing_mask)


    !############### Do Stuff #####################

    ! call HealpixMap_Read(BinaryMask, '/home/cosmos-tmp/aml1005/data/planck/mask_gal_CF857GHz_41pc.fits')
    ! call HealpixMap_Read(Mask, '/home/cosmos-tmp/aml1005/data/planck/mask_L3_v40.fits')
    ! BinaryMask%TQU=BinaryMask%TQU*Mask%TQU
    ! call HealpixMap_Write(BinaryMask, '/home/cosmos-tmp/aml1005/data/planck/mask_gal_CF857GHz_41pc_x_mask_L3_v40.fits')
    ! stop

    !  call CombineSplitNoiseMaps('/home/cosmos/users/aml1005/data/planck/MAP_v41_2048_GALACTIC_0240_11400_217Ghz')
    !  call CombineSplitNoiseMaps('/home/cosmos/users/aml1005/data/planck/MAP_v41_2048_GALACTIC_0240_11400_143Ghz')
    !  call CombineSplitNoiseMaps('/home/cosmos/users/aml1005/data/planck/MAP_v41_2048_GALACTIC_0240_11400_100Ghz')
    !  call MpiStop()

#ifdef MKLTHREADS
    call mkl_set_num_threads(1)
#endif

    ncrossweights = nweights*(nweights+1)/2
    ncl_tot = nweights*nchannels*(nweights*nchannels+1)/2
    if (want_pol) then
        vec_size=pol_vec_size
    else
        vec_size=1
    end if
    if (cross_spectra) then
        cross_noise_scale = 1 + 1/(real(TotYearWeightMaps())/nweights-1)
        if (MpiID==0) print *,'naive cross-noise scaling : ', cross_noise_scale
    end if
    !This scaling is by analogy with the exact full sky case, see appendix C of Hammimeche & Lewis

    !  call DifferenceMaps

    call HealpixInit(H,nside, 2*lmax,.true., w8dir=w8name,method=division_balanced)

    if (H%MpiID ==0) then !if we are main thread

    !print *,'start'
    if (add_planck_sim /= '') then
        call AddPlanckSimMaps(H)
        call MpiStop()
    end if

    print *,'Using nside = ', nside
    print *,'Point source A = ', point_source_A, ' +- ',point_source_A_frac_error*100,'%'

    if (cross_beams) then
        call ReadCrossBeamSet(beam_filename_format)
    else
        do channel = 1, nchannels
            if (Channels(channel)%Beam%beam_transfer) then
                beamfile = FormatFilename(beam_filename_format,Channels(channel)%Name)
                call ReadBeams(Channels(channel), beamfile, lmax)
            else
                call SetGaussianBeams(Channels(channel), lmax)
            end if
            Channels(channel)%PtSrcA =  PtScrcA(channel,channel)
            print *,'Channel '//trim(Channels(Channel)%Name)//' point source A = ', Channels(channel)%PtSrcA
        end do
    end if
    if (.not. no_pix_window) call PixSmoothBeams(lmax)
    print *,'done beams'

    if (apodize_mask_fwhm/=0.d0) then

    mask_fname = trim(data_dir)//trim(IntToStr(lmax))//'_nside'//trim(IntTOStr(nside))//'_fwhm' &
    //trim(RealToStr(real(apodize_mask_fwhm),4))
    file_stem = trim(file_stem)//'_apo'
    l_stem = trim(l_stem)//'_apo'

    call GetApodizedMask(H,Mask,trim(mask_fname)//'_'//trim(ExtractFileName(fits_mask)),fits_mask)
    if (pol_weights) then
        call GetApodizedMask(H,MaskPol,trim(mask_fname)//'_'//trim(ExtractFileName(fits_mask_pol)),fits_mask_pol)
    end if
    else
        !Don't apodize mask
        print *,'Reading BinaryMask'
        call HealpixMap_Read(BinaryMask, fits_mask)
        print *, 'Upgrading '
        call HealpixMap_udgrade(BinaryMask, Mask, nside, pessimistic=.true.)
        call HealpixMap_Free(BinaryMask)
        call HealpixMap_SetToIndexOnly(Mask,min(Mask%nmaps,2))
        call HealpixMap_ForceRing(Mask)
        if (pol_weights) then
            call HealpixMap_Read(BinaryMask, fits_mask_pol)
            print *, 'Upgrading '
            call HealpixMap_udgrade(BinaryMask, MaskPol, nside, pessimistic=.true.)
            call HealpixMap_Free(BinaryMask)
            call HealpixMap_SetToIndexOnly(MaskPol,min(MaskPol%nmaps,2))
            call HealpixMap_ForceRing(MaskPol)
        end if
    end if

    do channel= 1, nchannels
        call CombineYears(Channels(channel))
        if (est_noise) call EstAndSetNoise(Channels(channel))
        call ProcessNoiseMaps(H, Channels(Channel))
    end do

    if (want_cl) then
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
        if (pol_weights) then
            allocate(WeightMapsPol(nweights))
        end if
        do i=1,nweights
            call HealpixMap_Nullify(WeightMaps(i))
            call HealpixMap_Assign(WeightMaps(i),Mask)
            if (pol_weights) then
                call HealpixMap_Nullify(WeightMapsPol(i))
                call HealpixMap_Assign(WeightMapsPol(i),MaskPol)
            end if

            if (weights(i)=='uniform') then
            else if (weights(i)=='invnoise') then
                meanN = sum(SmoothedNoise%TQU(:,1), mask=SmoothedNoise%TQU(:,1) > 0)/count(SmoothedNoise%TQU(:,1) > 0)
                print *,'mean smoothed noise = ', meanN
                print *, 'sum mask', sum(WeightMaps(i)%TQU(:,1)) !this stops ifort 10.1 crashing on next line ???
                where (WeightMaps(i)%TQU(:,1)>0 .and. SmoothedNoise%TQU(:,1) > 0)
                    WeightMaps(i)%TQU(:,1) = WeightMaps(i)%TQU(:,1)*meanN/( SmoothedNoise%TQU(:,1) )
                end where
                print *, 'min/max Weightmap', minval(WeightMaps(i)%TQU(:,1)),maxval(WeightMaps(i)%TQU(:,1))

                if (pol_weights) then
                    where (WeightMapsPol(i)%TQU(:,1)>0 .and. SmoothedNoise%TQU(:,1) > 0)
                        WeightMapsPol(i)%TQU(:,1) = WeightMapsPol(i)%TQU(:,1)*meanN/( SmoothedNoise%TQU(:,1) )
                    end where
                end if
            else if (weights(i)=='mixedinvnoise') then
                call MpiStop('mixedinvnoise assumed not used')
                meanN = sum(SmoothedNoise%TQU(:,1), mask=SmoothedNoise%TQU(:,1) > 0)/count(SmoothedNoise%TQU(:,1) > 0)
                WeightMaps(i)%TQU(:,1) = WeightMaps(i)%TQU(:,1)*meanN/(SmoothedNoise%TQU(:,1)+meanN)
            else
                call MpiStop('unknown weight: '//trim(weights(i)))
            end if
        end do

        call HealpixMap_Free(SmoothedNoise)

        if (action/= pixtest .and. .not. do_exact_like) then
            call HealpixMap_Free(Mask)
            if (pol_weights) call HealpixMap_Free(MaskPol)
        end if

        print *,'Getting new weights power'
        if (.not. pol_weights) then
            call HealpixMapSet2CrossPowers(H, WeightMaps, MaskPowers, nweights, lmax*2,.false.)
        else
            allocate(WeightMapsAll(nweights*2))
            do i=1, nweights
                WeightMapsAll(i) = WeightMaps(i)
                WeightMapsAll(i+nweights) = WeightMapsPol(i)
            end do
            call HealpixMapSet2CrossPowers(H, WeightMapsAll, MaskPowersAll, 2*nweights, lmax*2,.false.)
            deallocate(WeightMapsAll)
        end if

        if (get_covariance) then
            print *,'getting weight powers for the covariance'
            call PseudoCl_WeightsToCovPowers(H, WeightMaps, WeightMapsPol, Channels, CovPowers, nweights, lmax*2, pol_weights)
            !       if (pol_weights) then
            !         !Pol T should give ordering with weight T_1 E_2 etc not T_2 E_1
            !!         call PseudoCl_WeightsToCovPowers(H, WeightMapsPol, WeightMaps, Channels, CovPowers(2), nweights, lmax*2)
            !        call PseudoCl_WeightsToCovPowers(H, WeightMapsPol, WeightMapsPol, Channels, CovPowers(3), nweights, lmax*2)
            !      end if
        end if
    end if !want_cl


    if (do_exact_like) then
        !Harmonic low-l exact likelihood (under devel)

        !  call HealpixMap_Init(NoiseMap, npix,pol = .false.)
        !  NoiseMap%TQU(:,1) = noise*NoiseMap%npix/(HO_fourpi)

        print *,'Doing exact with l_exact = ', l_exact, 'margin =', l_exact_margin
        sim_map_file = concat(trim(data_dir)//'sim_map',lmax,'_nside',nside)
        if (want_pol) sim_map_file = concat(sim_map_file,'_pol')
        if (uniform_noise) sim_map_file=concat(sim_map_file,'_uninoise')
        sim_map_file = concat(sim_map_file,'.fits')
        call InitRandom(rand_seed)
        RandInited = .true.

        ! if (get_sim .or. .not. FileExists(sim_map_file)) then
        if (get_sim ) then
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
            !call HealpixMap_Read(M, sim_map_file)

            print *,'reading map channel 1'
            do channel = 1, 1
                call ProcessChannelDatamap(H, Channels(channel), M)
            end do
            print *,'read map'
        end if
        do Channel =1,1
            call TestExactlike(H,Mask,MaskPol,Channels(Channel)%NoiseMap, M, fid_cls_file, cls_file, Coupler)
        end do
        call HealpixMap_Free(Mask)

        call mpiStop()
    end if

    end if  !MpiID=0

    call HealpixFree(H)


    if (MpiID==0) then
        print *,'Getting coupling matrix'
        allocate(MaskP(ncrossweights,3))
        if (.not. pol_weights) then
            call CrossPowersToHealpixPowerArray(MaskPowers, MaskP(1,1), dofree=.true.)
        else
            call CrossPowersToHealpixPowerArray2(MaskPowersAll, MaskP, dofree=.true.)
        end if
        print *,' getting coupling matrix arr'
        allocate(Coupler(ncrossweights))
    else
        Coupler => dummycoupler
    end if

    call PseudoCl_GetCouplingMatrixArr(Coupler, MaskP, lmin, lmax, want_pol, ncrossweights, pol_weights) !parallelized

    StTime = GeteTime()
    call PseudoCl_GetCouplingInversesArr(Coupler,ncrossweights)
    if (mpiId==0) then
        print *,'Coupling inversion time',  GeteTime()   - StTime
    end if
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

    if (get_covariance .and.  (detector_noise_filename_format/='') ) then
        if (MpiID==0) then
            print *,'Getting matrices for covariance'
            !Note we don't need all if nchanells >1 because noise is assumed uncorrelated between channels
            !Note covariance is calculated for combined year maps

            ncovpower = CovPowers%nmaps*(CovPowers%nmaps +1)/2
            print *,'ncovpower', ncovpower, nchannels, ncrossweights

            allocate(CovP(nCovPower,1))
            call CrossPowersToHealpixPowerArray(CovPowers, CovP(:,1), dofree=.true.)
            !if (pol_weights) then
            ! call CrossPowersToHealpixPowerArray(CovPowers(2), CovP(:,2), dofree=.true.)
            ! call CrossPowersToHealpixPowerArray(CovPowers(3), CovP(:,3), dofree=.true.)
            !end if
            !get all for simplicity, but not all entries are used (noise indep between channels)

            allocate(XiMatrices(nCovPower))
        else
            XiMatrices => dummycoupler
        end if

        !   call PseudoCl_GetCouplingMatrixArr(XiMatrices,CovP, lmin, lmax, want_pol,ncovpower, .false.)
        ! Now use only the T form
        call PseudoCl_GetCouplingMatrixArr(XiMatrices,CovP, lmin, lmax, .false.,ncovpower, .false.)
        !parallelized; no pol weights here as including all required variants explicitly in CovP array

        if (MpiID==0) then
            !crashes??
            ! print *,'freeing CovP'
            ! call PseudoCl_FreePowerArray2(CovP)
            ! print *,'deallocate'
            deallocate(CovP)

            print *,'Xi matrices ...',nCovPower, lmin,lmax
            print *,  XiMatrices(1)%T(2,2), XiMatrices(1)%T(3,3), XiMatrices(1)%T(4,4)
        end if
    end if

    call MpiQuietWait !low CPU usage so can use openmp until MPI wanted

    call HealpixInit(H,nside, lmax,.true., w8dir=w8name,method=division_balanced)

    if (H%MpiID ==0) then !if we are main thread

    if (get_covariance) then
        !Again only for combined-year maps

        do channel = 1, nchannels
            print *,'Get Chat Noise_l: '//trim(Channels(channel)%Name)

            allocate(Channels(channel)%NoiseP(ncrossweights))
            if (pol_weights) then
                call PseudoCl_GetCHatNoise(Channels(channel)%NoiseP, Coupler, &
                WeightMaps,WeightMapsPol,nweights, Channels(channel)%NoiseMap,noise_colour_factor)
            else
                call PseudoCl_GetCHatNoise(Channels(channel)%NoiseP, Coupler, &
                WeightMaps,WeightMaps,nweights, Channels(channel)%NoiseMap,noise_colour_factor)
            end if
            do i=1,ncrossweights
                call TBeam_PowerSmooth(Channels(channel)%Beam,Channels(channel)%NoiseP(i),+1)
                call HealpixPower_Write(Channels(channel)%NoiseP(i), &
                trim(file_stem)//Trim(IntToStr(i))//trim(Channels(channel)%Name)//'_noise.dat')
            end do
        end do

        print *,'reading fiducial power'
        call HealpixPower_Nullify(PFid)
        call HealpixPower_ReadFromTextFile(PFid,fid_cls_file,lmax,pol=want_pol,dolens = .false.)
        print *,'PFid...', PFid%cl(2:5,1)
        print *,'Getting covariance'
#ifdef MKLTHREADS
        call mkl_set_num_threads(mkl_threads)
#endif
        !$ if (use_openmp) then
        !$ call omp_set_num_threads(max_threads)
        !$ print *, 'setting threads = ', max_threads
        !$ end if

        call PseudoCl_GetFullCovariance(Coupler, XiMatrices, CovArr(1), PFid,vec_size, Channels, nweights,  &
        .true.,.true.,point_source_A_frac_error, pol_weights,0,1)

        call PseudoCl_GetFullCovariance(Coupler, XiMatrices, CovArr(1), PFid,vec_size, Channels, nweights,  &
        .true.,.true.,point_source_A_frac_error, pol_weights,1,1)

        if (get_signal_covariance) then
            call PseudoCl_GetFullCovariance(Coupler, XiMatrices, CovArr(2), PFid,vec_size, Channels, nweights,  &
            .false.,.true.,0.0, pol_weights,1,1)
        end if
        if (get_noise_covariance) then
            call PseudoCl_GetFullCovariance(Coupler, XiMatrices, CovArr(3), PFid,vec_size, Channels, nweights,  &
            .true.,.false.,0.0, pol_weights,1,1)
        end if

        !$ if (use_openmp) then
        !$ call omp_set_num_threads(1)
        !$ end if

        if (.not. get_sim .and. .not. get_analysis) then
            call TCouplingMatrix_ArrayFree(Coupler)
            deallocate(Coupler)
        end if
        call TCouplingMatrix_ArrayFree(XiMatrices)
        deallocate(XiMatrices)

        call HealpixPower_Init(SigNoiseCov,lmax, want_pol)
        call HealpixPower_Init(OffDiagCov,lmax, want_pol)
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
        !$ if (use_openmp) then
        !$ call omp_set_num_threads(max_threads)
        !$ print *, 'setting threads = ', max_threads
        !$ end if

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
        if (use_openmp) then
            !$ call omp_set_num_threads(1)
        end if


        deallocate(tmp, tmpcov)

        HyCov(1)%C = HyCov(1)%C*mK**4
        print *, 'Cov 22 = ',HyCov(1)%C(2,2)
        call  MatrixSym_Write_Binary_Single(concat( trim(data_dir)//trim(out_file_base)// &
        '_hybrid_pol',vec_size,'_lmax',lmax,'_',ncl_tot,'.covmat'),HyCov(1)%C)

        HyCov(1)%C = HyCov(1)%C/mK**4

        print *,'writing files'
        SigNoiseCov%Cl=0
        SignalCov%Cl=0 !pure signal variance
        NoiseCov%Cl=0 !pure noise variance
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
                SigNoiseCov%Cl(i+lmin-1,C_C) = sqrt(fac*HyCov(1)%C(nl+i,nl+i))
                SigNoiseCov%Cl(i+lmin-1,C_E) = sqrt(fac*HyCov(1)%C(nl*2+i,nl*2+i))
                OffDiagCov%Cl(i+lmin-1,C_T) = fac*HyCov(1)%C(i,nl+i)
                OffDiagCov%Cl(i+lmin-1,C_C) = fac*HyCov(1)%C(i,nl*2+i)
                OffDiagCov%Cl(i+lmin-1,C_E) = fac*HyCov(1)%C(nl+i,nl*2+i)

                if (get_signal_covariance) SignalCov%Cl(i+lmin-1,C_C) = sqrt(fac*HyCov(2)%C(nl+i,nl+i))
                if (get_noise_covariance)  NoiseCov%Cl(i+lmin-1,C_C) = sqrt(fac*HyCov(3)%C(nl+i,nl+i))
                if (get_signal_covariance) SignalCov%Cl(i+lmin-1,C_E) = sqrt(fac*HyCov(2)%C(nl*2+i,nl*2+i))
                if (get_noise_covariance)  NoiseCov%Cl(i+lmin-1,C_E) = sqrt(fac*HyCov(3)%C(nl*2+i,nl*2+i))

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
                    OffDiagCov%Cl(i+lmin-1,C_B) = fac*HyCov(1)%C(3*nl+i,nl*2+i)
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
        call HealpixPower_Write(OffDiagCov,trim(file_stem)//'_vec'//trim(IntToStr(vec_size))//'_'// &
        'joint_offdiagcov.dat')

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

        if (.not. get_sim) deallocate(HyCov(1)%C)
        if (get_signal_covariance) deallocate(HyCov(2)%C)
        if (get_noise_covariance) deallocate(HyCov(3)%C)

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

    call MpiWakeQuietWait !low CPU usage so can use openmp until MPI wanted

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
            call ProcessChannelDatamap(H, Channels(channel), MapArray(Channel))
        end do

        print *,'Analysing maps'

        call AnalyseMapAuto(H, MapArray, DataCl)
        call HealpixPower_Write(DataCl,trim(anastem)//'_unsubtracted_data_cls.dat')
        DataCl%Cl(:,C_T) = DataCl%Cl(:,C_T) - PointSourceP%Cl(:,C_T)
        DataCl%Cl = DataCl%Cl - HybridNoise%Cl*noise_adjustment
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
        call HealpixPower_Init(COffvar,lmax, want_pol)

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
        call AnalyseMapAuto(H, MapArray, HybridP)
        HybridP%Cl(:,1) = HybridP%Cl(:,1) - PointSourceP%Cl(:,1)
        do channel=1, nchannels
            call HealpixMap_Free(MapArray(channel))
        end do
        end if

        print *,'AnalyseMapAuto time:',  GeteTime() -StTime

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
        COffVar%Cl(:,C_T) = COffVar%Cl(:,C_T) + (HybridP%Cl(:,C_T)-HybridNoise%Cl(:,C_T) - PSim%Cl(:,C_T))* &
        (HybridP%Cl(:,C_C)- PSim%Cl(:,C_C))
        COffVar%Cl(:,C_C) = COffVar%Cl(:,C_C) + (HybridP%Cl(:,C_T)-HybridNoise%Cl(:,C_T)- PSim%Cl(:,C_T))* &
        (HybridP%Cl(:,C_E) - HybridNoise%Cl(:,C_E)- PSim%Cl(:,C_E))
        COffVar%Cl(:,C_E) = COffVar%Cl(:,C_E) + (HybridP%Cl(:,C_E)-HybridNoise%Cl(:,C_E)- PSim%Cl(:,C_E))* &
        (HybridP%Cl(:,C_C) - PSim%Cl(:,C_C))
        if (vec_size>3) then
            COffVar%Cl(:,C_B) = COffVar%Cl(:,C_B) + (HybridP%Cl(:,C_E)-HybridNoise%Cl(:,C_E)- PSim%Cl(:,C_E))* &
            (HybridP%Cl(:,C_B)-HybridNoise%Cl(:,C_B)- PSim%Cl(:,C_B))
        end if

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
            CAvgTemp%Cl =  COffVar%Cl/sim_no
            call HealpixPower_Write(CAvgTemp,trim(sim_stem)//'_cls_hybrid_offdiagvar.dat')

            if (get_covariance) then
                call SplineFitCl(CAvgTemp,SigNoiseCov)
                call HealpixPower_Write(CAvgTemp,trim(sim_stem)//'_cls_hybrid_var_fit'//trim(intToStr(sim_no))//'.dat')
                call HealpixPower_Write(CAvgTemp,trim(sim_stem)//'_cls_hybrid_var_fit.dat')
            end if
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

        if (get_covariance) then
            HyCov(1)%C = HyCov(1)%C*mK**4
            do i=1,nl*vec_size
                bigvec(i) = simcov(i,i)/HyCov(1)%C(i,i)
            end do
            do i=1, vec_size
                call SplineSmoothSampleArray(bigvec((i-1)*nl+1:), nl)
            end do

            !        do i = -smooth_w*2, smooth_w*2
            !         smooth_kernel(i) = exp(-i**2/real(2*smooth_w**2))
            !        end do
            !        smooth_kernel= smooth_kernel/sum(smooth_kernel)
            !
            !        allocate(BigVec2(nl))
            !        do j=1,vec_size
            !         bigVec2 = bigvec( (j-1)*nl+1:j*nl )
            !         do i=1+smooth_w*2,nl-smooth_w*2
            !              bigvec((j-1)*nl+i) = sum(bigvec2(i-smooth_w*2:i+smooth_w*2)*smooth_kernel)
            !         end do
            !        end do


            do i=1,nl*vec_size
                HyCov(1)%C(:,i) = HyCov(1)%C(:,i)*sqrt(bigVec(i))
                HyCov(1)%C(i,:) = HyCov(1)%C(i,:)*sqrt(bigVec(i))
            end do
            !        deallocate(BigVec2)

            call  MatrixSym_Write_Binary_Single(concat( trim(data_dir)//trim(out_file_base)// &
            '_hybrid_pol',vec_size,'_lmax',lmax,'_',ncl_tot,'_simcalib.covmat'),HyCov(1)%C)

            HyCov(1)%C = HyCov(1)%C/mK**4
            SigNoiseCov%Cl=0 !pure signal variance
            fac=1
            do i=1,nl
                SigNoiseCov%Cl(i+lmin-1,C_T) = sqrt(fac*HyCov(1)%C(i,i))
                if (vec_size>1) then
                    SigNoiseCov%Cl(i+lmin-1,C_C) = &
                    sqrt(fac*HyCov(1)%C(nl+i,nl+i))
                    SigNoiseCov%Cl(i+lmin-1,C_E) = &
                    sqrt(fac*HyCov(1)%C(nl*2+i,nl*2+i))
                    if (vec_size>3) then
                        SigNoiseCov%Cl(i+lmin-1,C_B) = &
                        sqrt(fac*HyCov(1)%C(nl*3+i,nl*3+i))
                    end if
                end if
            end do
            print *,'write calib rated diag cov'
            call HealpixPower_Write(SigNoiseCov,trim(file_stem)//'_vec'//trim(IntToStr(vec_size))//'_'// &
            'joint_diagcov_simcalib.dat')
        end if

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
