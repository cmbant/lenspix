    !Program to make pure B mode map using supported modes following method of astro-ph/0305545
    !AL April 2013
    program EBsep
    use CutSkyAsymm
    use HealpixObj
    use HealpixVis
    use spinalm_tools
    use IniFile
    use AmlUtils
    implicit none

    real(dp) :: SUPPORT = 0.99d0, filter_support = 0.9d0, EB_support = 0.99d0
    Type (AsymmCouplings) :: Couplings, Couplings2
    character(LEN=1024) :: healpixloc, w8name
    integer :: lmax, status
    Type(HealpixInfo)  :: H
    Type(HealpixMap) Mask, InMap, LoadMap, FiltMap, CutMap
    Type(HealpixAlm) MaskAlm, Alm
    Type(ProjMat) :: Proj, Proj2, FilterProj
    Type(ProjMatPol) ProjPol
    integer unitm, num_threads
    integer :: nside, lmaxfile, unit
    character(LEN=:), allocatable :: matrix_name, mask_file, out_dir, in_map, output_map, in_map_diff
    logical :: B_only = .true.
    logical :: unit_mask = .true.
    logical :: do_E= .false.
    logical :: do_L_filter = .false.
    integer i
    real(dp), pointer :: rotW(:,:), tmpM(:,:)
    real(dp) :: smooth_mask = 0.d0, rotate_angle=0.d0

    call get_environment_variable('HEALPIX', healpixloc, status=status)
    if (status==0) then
        w8name = trim(healpixloc)//'/data/'
    end if

    call Ini_Open(GetParam(1),1)
    support = Ini_read_Real('support')
    filter_support = Ini_read_real('filter_support')
    nside = Ini_read_Int('nside')
    lmax = Ini_read_Int('lmax')

    B_only = Ini_Read_Logical('B_only', B_only)
    mask_file = trim(Ini_read_String('mask_file'))
    smooth_mask = Ini_read_Double('smooth_mask',smooth_mask)
    unit_mask = Ini_read_Logical('unit_mask',unit_mask)
    if(.not. unit_mask) EB_support = Ini_read_Double('EB_support',EB_support)
    do_E= Ini_Read_Logical('do_E',do_E)
    do_L_filter = Ini_Read_Logical('do_L_filter',do_L_filter)
    out_dir = trim(Ini_read_String('out_dir'))
    num_threads=  Ini_read_int('num_threads')
    if (GetParamCount()>2) then
        in_map = trim(GetParam(2))
        output_map = trim(GetParam(3))
    else
        in_map = trim(Ini_Read_String('input_map'))
        output_map = trim(Ini_Read_String('output_map'))
        if (in_map/='' .and. output_map=='') stop 'no output_map specified'
        in_map_diff = trim(Ini_Read_String('input_map_diff'))
    end if
    i = index(output_map,'.fits')
    if (i>0) then
        output_map = output_map(1:i-1)
    end if
    rotate_angle=Ini_read_Double('rotate_angle',rotate_angle)
    call Ini_Close()
    if (num_threads>0) call mkl_set_num_threads(num_threads)
    call HealpixInit(H, nside, lmax*2,.true., w8dir=w8name, method= division_equalrows)

    matrix_name = out_dir//trim(ExtractFileName(mask_file))//'_lmax'//trim(IntToStr(lmax))//'_support'//trim(RealToStr(real(Support)))
    if (.not. B_only) matrix_name = matrix_name //'_EB'
    if (smooth_mask>0) matrix_name = matrix_name //'_smth'//trim(RealToStr(real(smooth_mask)))
    if (.not. unit_mask) matrix_name = matrix_name //'_apo'
    call HealpixMap_Read(Mask, mask_file)
    print *, 'fsky=',sum(Mask%TQU)/Mask%npix
    if (.not. unit_mask) Mask%TQU = Mask%TQU**2
    call HealpixMap2Alm(H,Mask,MaskAlm,lmax*2, dopol=.false.)
    if (smooth_mask >0) then
        call HealpixAlm_Smooth(MaskAlm, smooth_mask)
        if (.not. unit_mask) stop 'err smooth'
        call HealpixAlm2Map(H,MaskAlm,Mask, Mask%npix)
        print *, 'smooth fsky=',sum(Mask%TQU)/Mask%npix
    end if

    !/scratch/maja1/planck/repository/exchanges/dx11/maps/hfi/HFI_SkyMap_353_2048_DX11c_full.fits
    if (.not. FileExists(matrix_name)) then
        !First get alm of mask
        !Calculate cut-sky coupling matrix
        call CutSkyAsymm_GetCoupling(Couplings,MaskAlm, MaskAlm,lmax, plusonly =B_only, WantTemp=.false., WantPol=.true.)
        if (B_only) then
            print *, 'Getting supported pure mode projection matrix'
            call CutSkyAsymm_GetSupportedModes(Couplings%WPAsymm, Proj, SUPPORT)
            deallocate(Couplings%WPAsymm)
            if (.not. unit_mask) then
                !Want simultaneous eigenvalues of W^2 and W^{(2)}
                !Have diagonalized W^{(2)}
                !Now make D^{-1/2} Proj W^2 Proj^T D^{-1/2}, which have unit eigenvalue if eigenvalues of W^(2) and W^2 are shared
                call HealpixMap_Read(Mask, mask_file)
                call HealpixMap2Alm(H,Mask,MaskAlm,lmax*2)
                call HealpixMap_Free(Mask)
                call CutSkyAsymm_GetCoupling(Couplings2,MaskAlm, MaskAlm,lmax, plusonly =.true., WantTemp=.false., WantPol=.true.)
                call HealpixAlm_Free(MaskAlm)
                allocate(tmpM(Proj%nl,Proj%nr))
                print *,'symmleft'
                call Matrix_Mult_SymmLeft(Couplings2%WPAsymm,Proj%M, tmpM)
                deallocate(Couplings2%WPAsymm)
                allocate(RotW(Proj%nr,Proj%nr))
                print *,'mult'
                call Matrix_Mult_TN(tmpM,tmpM, rotW)
                deallocate(tmpM)
                print *,'norm'
                do i=1, Proj%nr
                    RotW(i,:) = RotW(i,:)/Proj%RootDiag(i)
                    RotW(:,i) = RotW(:,i)/Proj%RootDiag(i)
                end do
                call CutSkyAsymm_GetSupportedModes(RotW, Proj2, EB_support)
                deallocate(RotW)
                do i=1, Proj%nr
                    Proj%M(:,i) = Proj%M(:,i)/Proj%RootDiag(i)
                end do
                allocate(tmpM(Proj%nl,Proj2%nr))
                call Matrix_Mult(Proj%M,Proj2%M, tmpM)
                do i=1, Proj%nr
                    Proj%M(:,i) = Proj%M(:,i)*Proj%RootDiag(i)**2
                end do
                allocate(Proj%Minv(Proj%nl,Proj2%nr))
                call Matrix_Mult(Proj%M,Proj2%M, Proj%Minv)
                deallocate(Proj%M, Proj2%M)
                Proj%nr = Proj2%nr
                Proj%M => tmpM
            end if

            !Save big matrix result
            open(newunit=unit, file = matrix_name, form='unformatted',status='replace')
            write(unit) lmax
            write(unit) Proj%nl,Proj%nr
            if (.not. unit_mask) then
                write(unit) Proj%Minv
            else
                write(unit) Proj%RootDiag
            end if
            write(unit) Proj%M
            close(unit)
        else
            print *, 'Getting supported EB mode projection matrix'
            call CutSkyAsymm_SupportedModesPolFromCoupling(Couplings, ProjPol, SUPPORT)
            !Save big matrix result
            open(newunit=unit, file = matrix_name, form='unformatted',status='replace')
            write(unit) lmax
            write(unit) ProjPol%nl,ProjPol%nr
            write(unit) ProjPol%RootDiag
            write(unit) ProjPol%WComp
            close(unit)
        end if

        print *, 'Done: wrote '//matrix_name
    else
        if (in_map=='') stop 'Nothing to do, computed matrix already exists'
        print *, 'reading stored matrix: '//matrix_name
        open(newunit=unit, file = matrix_name, form='unformatted',status='old')
        read(unit) lmaxfile
        if (lmaxfile/=lmax) stop 'lmax mismatch'
        read(unit) Proj%nl,Proj%nr
        allocate(Proj%M(Proj%nl,Proj%nr))
        if (unit_Mask) then
            allocate(Proj%RootDiag(Proj%nr))
            read(unit) Proj%RootDiag
        else
            allocate(Proj%Minv(Proj%nl,Proj%nr))
            read(unit) Proj%Minv
        end if
        read(unit) Proj%M
        close(unit)
        print *, 'Read matrix'
    end if

    if (in_map /= '') then
        call HealpixMap_Read(LoadMap, in_map)
        if (in_map_diff/='') then
            call HealpixMap_Read(InMap, in_map_diff)
            LoadMap%TQU = (LoadMap%TQU - InMap%TQU)/2
        endif
        if (LoadMap%nside /= nside) then
            print *,'Degrading to, nside = ', nside
            call HealpixMap_udgrade(LoadMap, InMap, nside)
            call HealpixMap_Assign(LoadMap,InMap)
        else
            call HealpixMap_Assign(InMap, LoadMap)
        end if
        print *,'InMap degraded min/max Q,U', minval(InMap%TQU(:,2:3)),maxval(InMap%TQU(:,2:3))
        call HealpixMap_ForceRing(InMap)

        if (rotate_angle/=0.d0) then
            call HealpixMap_ForceRing(LoadMap)
            print *,'rotating by', rotate_angle
            InMap%TQU(:,2) = LoadMap%TQU(:,2)*cos(2*rotate_angle/180*HO_PI) - LoadMap%TQU(:,3)*sin(2*rotate_angle/180*HO_PI)
            InMap%TQU(:,3) = LoadMap%TQU(:,3)*cos(2*rotate_angle/180*HO_PI) + LoadMap%TQU(:,2)*sin(2*rotate_angle/180*HO_PI)
        end if
        call HealpixMap_Free(LoadMap)
        call HealpixMapMulCut(InMap,Mask,CutMap,1)
        print *,'cut min/max/sum Q,U', minval(CutMap%TQU(:,2:3)),maxval(CutMap%TQU(:,2:3)),sum(CutMap%TQU(:,2:3))

        call B_filter_map(CutMap, Proj, FiltMap, Alm, polix=3)
        !   print *, 'min, max B map', minval(InMap%TQU(:,2:3)), maxval(InMap%TQU(:,2:3))
        print *, 'min, max filtered Q/U', minval(FiltMap%TQU(:,2:3)), maxval(FiltMap%TQU(:,2:3))
        print *, 'min, max filtered B', minval(FiltMap%TQU(:,1)), maxval(FiltMap%TQU(:,1))
        !        call HealpixMap_Write(InMap, in_map//'_Bonly',  overwrite=.true.)
        call HealpixMap_Write(FiltMap, output_map//'.fits',  overwrite=.true.)
        if (do_E) then
            call B_filter_map(CutMap, Proj, FiltMap, Alm, polix=2)
            print *, 'min, max filtered Q/U', minval(FiltMap%TQU(:,2:3)), maxval(FiltMap%TQU(:,2:3))
            print *, 'min, max filtered E', minval(FiltMap%TQU(:,1)), maxval(FiltMap%TQU(:,1))
            call HealpixMap_Write(FiltMap, output_map//'_E.fits',  overwrite=.true.)
        end if

        if (do_L_filter) then
            call getLfilterProj(Proj, FilterProj)
            print *,'filtering'
            !        call HealpixMapMulCut(InMap,Mask,CutMap,1)
            call B_filter_map(CutMap, FilterProj, FiltMap,Alm, polix=3)
            call HealpixMap_Write(FiltMap, output_map//'_filt.fits',  overwrite=.true.)
            if (do_E) then
                call B_filter_map(CutMap, FilterProj, FiltMap,Alm, polix=2)
                call HealpixMap_Write(FiltMap, output_map//'_Efilt.fits',  overwrite=.true.)
            end if
        end if

    end if
    call HealpixFree(H)

    contains

    subroutine AlmToBVec(Alm, vec, polix)
    real(sp), allocatable :: vec(:)
    Type(HealpixAlm) :: Alm
    integer polix

    allocate(vec((lmax+1)**2-4))
    call HealpixAlm2Vector(Alm, vec, lmax, polix)
    end subroutine AlmToBVec

    subroutine B_filter_map(M, Proj, ProjMap, Alm, polix)
    !polix = 3, B map, polix=2 actually E map
    Type(HealpixMap) M, ProjMap
    real(sp), allocatable :: vec(:), modes(:)
    Type(ProjMat) :: Proj
    integer polix
    Type(HealpixAlm) :: Alm
    integer i

    call HealpixMap2Alm(H,M,Alm,lmax, dopol=.true.)
    call AlmToBVec(Alm, vec, polix)

    allocate(modes(Proj%nr))
    do i=1, Proj%nr
        modes(i) = dot_product(vec,Proj%M(:,i))
    end do
    if (associated(Proj%Minv)) then
        vec = matmul(Proj%Minv, modes)
    else
        vec = matmul(Proj%M, modes)
    end if
    call HealpixAlm_Init(Alm,lmax,3)
    call HealpixVector2Alm(vec, Alm, lmax, polix=polix)
    ALM%TEB(1,:,:) = Alm%TEB(polix,:,:)
    call HealpixAlm2Map(H, Alm, ProjMap, M%npix)

    end subroutine B_filter_map


    subroutine GetLCov(Couple, Cov, lmin, lmax, pow, scale)
    real(dp), allocatable :: Cov(:,:)
    Type(ProjMat) :: Couple
    !    Type(HealpixPower) :: P
    integer, intent(in) :: lmin, lmax
    real(dp) pow, scale
    integer i,j
    real(dp) :: tmp
    integer l, mix
    !Couple is the transpose so fast sums in dot_product

    if (.not. Allocated(Cov)) allocate(Cov(Couple%nr,Couple%nr), source=0.d0)
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(DYNAMIC), Private(mix, tmp ,j, l)
    do i=1, Couple%nr
        do j=1, i
            tmp = 0
            do l=lmin,lmax
                mix = idx_P(l,-l)
                tmp =tmp + dot_product(Couple%M(mix:mix+2*l,i),Couple%M(mix:mix+2*l,j))*(120./l)**pow*scale
            end do
            if (associated(Couple%RootDiag)) tmp = tmp*Couple%RootDiag(i)*Couple%RootDiag(j)
            Cov(i,j) = Cov(i,j) + tmp
            if (i/=j) Cov(j,i) =  Cov(i,j)
        end do
    end do
    !$OMP END PARALLEL DO

    end subroutine GetLCov

    subroutine getLfilterProj(Proj, FilterProj)
    real(dp), allocatable :: rotCov(:,:), Cov(:,:), N(:,:)
    Type(ProjMat) :: Proj, HasSignalProj,FilterProj
    !Type(THealpixPower) P

    !    call HealpixPower_ReadFromTextFile(P,'',lmax,pol=.true.,dolens = .false.)

    print *, 'get covs'
    call GetLCov(Proj, Cov, 50, 120, pow=2.d0, scale=1.d0)
    allocate(N, source=Cov)
    call GetLCov(Proj, N, 121, lmax, 0.d0, scale=1.d0)
    ! N =N*100
    call GetLCov(Proj, N, 2, 49, pow=0.d0, scale=50.d0)

    print *, 'get root inverse'
    call Matrix_CholeskyRootInverse(N, transpose=.true.)
    print *, 'rotate'
    allocate(rotCov(Proj%nr,Proj%nr))
    call Matrix_RotateSymm(Cov, N, Proj%nr, rotCov, triangular = .true.)
    print *,'get supported'
    call CutSkyAsymm_GetSupportedModes(rotCov, HasSignalProj, filter_support)
    print *,'getting filtered projection'
    allocate(FilterProj%M(Proj%nl, HasSignalProj%nr))
    FilterProj%nr = HasSignalProj%nr
    FilterProj%nl = Proj%nl
    FilterProj%M = matmul(Proj%M,HasSignalProj%M)
    If (associated(Proj%Minv)) then
        allocate(FilterProj%Minv(Proj%nl, HasSignalProj%nr))
        FilterProj%Minv = matmul(Proj%Minv,HasSignalProj%M)
    end if
    end subroutine getLfilterProj

    end program EBsep