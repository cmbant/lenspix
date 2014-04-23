    !Functions for calculating coupling matrices and covariances, etc.
    !Desgined for use on a cluster 2+ CPUs using MPI, not so well tested on single CPU
    !This version AL: April 2008 - Dec 2010
    module PseudoCl
    use healpix_types, ONLY: SP,DP,I4B,SPC
    use HealpixObj
    use MpiStuff
    USE spinalm_tools
    implicit none

    Type TBeam
        logical :: beam_transfer
        real(dp) :: fwhm
        real(dp), pointer :: beam(:) => NULL()
    end type TBeam

    Type TCrossBeamSet
        Type(TBeam), allocatable :: Beams(:,:)
    end Type

    Type TChannel
        character(LEN=16) :: Name
        character(LEN=256) :: noise_file
        integer :: Count  !Number of detectors
        Type (TStringList) :: DetectorNames
        logical, pointer :: detector_has_pol(:) => NULL() !true if TQU, false if just T
        real(dp), pointer  :: sig0(:)  => NULL()
        real(dp) :: Ghz
        real(dp) :: PtSrcA !non-beam-smoothed C_l of fiducial point sources
        real(dp) :: ENoiseFac
        Type(TBeam), pointer :: DetectorBeams(:)  => NULL()
        Type(HealpixMap), pointer :: DetectorYearNoiseMaps(:,:) => NULL()
        Type(TBeam) :: Beam
        Type(HealpixMap) :: NoiseMap
        Type(HealpixPower), pointer :: NoiseP(:) => NULL()
    end Type TChannel


    Type TCouplingMatrix
        integer lmin
        integer lmax
        integer nl
        logical has_pol
        real(dp), dimension(:,:), pointer :: T => NULL() , X => NULL() , &
        EE => NULL() , EB => NULL()
        real(dp), dimension(:,:), pointer ::  InvT => NULL() ,InvX => NULL() , &
        InvEE => NULL() , InvEB => NULL()
    end type TCouplingMatrix

    Type TCovMat
        real(dp), dimension(:,:), pointer :: C  => NULL()
    end Type TCovMat

    Type TCovMatArray
        integer ncl
        Type(TCovMat), dimension(:,:), pointer :: Cov  => NULL()
    end Type TCovMatArray

    Type TCovMatPolArray
        integer vec_size
        Type(TCovMatArray), dimension(:,:), pointer :: Pol  => NULL()
    end Type TCovMatPolArray

    Type TCovMatSet
        integer n, lmin, lmax
        Type(TCovMat), dimension(:), pointer :: Cov  => NULL()
    end Type TCovMatSet


    Type TClArray
        Type(HealpixPower), dimension(:), pointer :: P  => NULL()
    end Type TClArray

    Type TCouplingMatrixArray
        Type(TCouplingMatrix), dimension(:), pointer :: M  => NULL()
    end Type TCouplingMatrixArray


    Real(dp):: NoiseScale
    !Avoid very small numbers is noise*w_i*w_j maps;
    !set by WeightsToCovPowers
    !used by GetFullCovariance
    real(dp) :: cross_noise_scale = 1.d0

    contains

    subroutine TBeam_Free(B)
    Type(TBeam) :: B
    integer status

    deallocate(B%Beam,stat = status)
    nullify(B%Beam)

    end subroutine TBeam_Free

    subroutine TBeam_Assign(B, Bin)
    Type(TBeam) B, Bin

    call TBeam_Free(B)
    B%beam_transfer = Bin%beam_transfer
    B%fwhm = Bin%fwhm
    if (B%beam_transfer) then
        allocate(B%beam(0:size(Bin%beam)-1))
        B%beam = Bin%Beam
    end if

    end subroutine TBeam_Assign


    subroutine TBeam_SetGaussian(B, lmax)
    Type(TBeam) :: B
    integer, intent(in) :: lmax
    real xlc, sigma2
    integer l

    B%beam_transfer = .true.
    allocate(B%beam(0:lmax))
    B%beam = 1
    xlc= 180*sqrt(8.*log(2.))/HO_pi
    sigma2 = (B%fwhm/xlc)**2

    do l=1,lmax
        B%beam(l) =  exp(-l*(l+1)*sigma2/2)
    end do

    end subroutine TBeam_SetGaussian


    subroutine TBeam_ReadFile(B, filename, lmax, col)
    Type(TBeam) :: B
    character(LEN=*), intent(in) :: filename
    integer, intent(in), optional :: col
    integer, intent(in) :: lmax
    integer i,l, beamcol
    real(dp) inB
    character(LEN=1024*8) :: inline
    real(dp), allocatable :: cols(:)

    if (present(col)) then
        beamcol=col
    else
        beamcol=2
    end if
    allocate(cols(col-1))
    B%beam_transfer = .true.
    allocate(B%beam(0:lmax))
    call OpenTxtFile(filename,1)
    do
        read (1,'(a)', err = 500, end =500) inline
        if (inline(1:1) /= '#' .and. inline/='') then
            read(inline,*) l, cols
            if (l<=lmax) B%beam(l) = cols(beamcol-1)
        end if
    end do

500 close(1)
    if (l<lmax) call MpiStop('TBeam_ReadFile: not high enough l in file')

    end subroutine TBeam_ReadFile


    subroutine TBeam_PowerSmooth(B, P, sgn)
    Type(TBeam) :: B
    Type (HealpixPower)::P
    integer, intent(in) :: sgn

    if (P%lmax > size(B%Beam)-1) call MpiStop('TBeam_PowerSmooth: lmax mismatch')

    if (B%beam_transfer) then
        call HealpixPower_Smooth_Beam(P,B%beam,sgn)
    else
        call HealpixPower_Smooth(P,B%fwhm,sgn)
    end if
    end subroutine TBeam_PowerSmooth

    subroutine TBeam_PowerSmooth2(B1,B2, P, sgn)
    Type(TBeam) :: B1, B2
    Type (HealpixPower)::P
    integer, intent(in) :: sgn

    if (P%lmax > size(B1%Beam)-1) call MpiStop('TBeam_PowerSmooth2: lmax mismatch')

    if (.not. B1%beam_transfer) call MpiStop('TBeam_PowerSmooth2: want transfers now')

    call HealpixPower_Smooth_Beam2(P,B1%beam,B2%beam,sgn)

    end subroutine TBeam_PowerSmooth2





    subroutine TCouplingMatrix_Free(M)
    Type(TCouplingMatrix) :: M
    integer status

    deallocate(M%T,stat=status)
    deallocate(M%X, M%EE, M%EB, stat=status)
    deallocate(M%InvT,stat=status)
    deallocate(M%InvX, M%InvEE,M%InvEB, stat=status)
    call TCouplingMatrix_Nullify(M)

    end subroutine TCouplingMatrix_Free

    subroutine TCouplingMatrix_ArrayFree(C)
    Type(TCouplingMatrix) :: C(:)
    integer i

    do i=1,size(C)
        call TCouplingMatrix_Free(C(i))
    end do

    end subroutine TCouplingMatrix_ArrayFree

    subroutine TCouplingMatrix_Nullify(M)
    Type(TCouplingMatrix) :: M

    nullify(M%T,M%X,M%EE,M%EB)
    nullify(M%InvT,M%InvX,M%InvEE,M%InvEB)

    end subroutine TCouplingMatrix_Nullify

    subroutine PseudoCl_GetCouplingInversesArr(M, n)
    !Get M^{-1} by inverting symmetric form and putting in l factors
    use MatrixUtils
    integer, intent(in) :: n
    Type(TCouplingMatrix) :: M(n)
    Type(TCouplingMatrix), pointer :: AM
    integer lmin,lmax,l
    real(dm), dimension(:,:), allocatable :: tmp,Inv11
    integer MpiID, MpiSize, i,j
    integer params(3)
    logical haspol
    Type(TMatrixType), dimension(:), allocatable :: Arr
    Type(TMatrixType) dummyArr(1)
    integer nmat

    call GetMpiStat(MpiID, MpiSize)

    if (MpiID==0) then
        nmat = n
        lmin=M(1)%lmin
        lmax=M(1)%lmax

        if (M(1)%has_pol) nmat=n*3
        print *,'nmat = ',nmat
        allocate(Arr(nmat))
        do i=1,n
            allocate(Arr(i)%M(lmin:lmax,lmin:lmax))
            Arr(i)%M = M(i)%T
            if (M(i)%has_pol) then
                print *,'inverting pol'

                allocate(Arr(i+n)%M(lmin:lmax,lmin:lmax))
                Arr(i+n)%M = M(i)%EE
                allocate(Arr(i+2*n)%M(lmin:lmax,lmin:lmax))
                Arr(i+2*n)%M = M(i)%X
            end if
        end do
        call Matrix_InverseArrayMPI(Arr,nmat)
        do i=1,n
            M(i)%InvT => Arr(i)%M
            do l=lmin,lmax
                M(i)%InvT(l,:)=M(i)%InvT(l,:)/(2*l+1)
            end do
            if (M(i)%has_pol) then
                M(i)%InvEE => Arr(i+n)%M
                M(i)%InvX => Arr(i+2*n)%M

                allocate(inv11(lmin:lmax,lmin:lmax))
                allocate(tmp(lmin:lmax,lmin:lmax))

                Inv11= M(i)%InvEE
                !M%InvEE = M%EE - matmul(M%EB,matmul(Inv11,M%EB))

                call Matrix_Mult_SymmRight(Inv11,M(i)%EB,tmp)
                call Matrix_Mult(M(i)%EB,tmp,M(i)%InvEE)
                do l=lmin,lmax
                    M(i)%InvEE(:,l) = M(i)%EE(:,l) - M(i)%InvEE(:,l)
                end do
                call Matrix_inverse(M(i)%InvEE)

                print *,'getting InvEB'
                call Matrix_Mult_SymmRight(M(i)%InvEE,M(i)%EB,tmp)

                allocate(M(i)%InvEB(lmin:lmax,lmin:lmax))

                call Matrix_Mult_SymmRight(tmp,inv11,M(i)%InvEB,-1._dm)

                deallocate(Inv11,tmp,Arr)

                do l=lmin,lmax
                    M(i)%InvX(l,:)=M(i)%InvX(l,:)/(2*l+1)
                    M(i)%InvEE(l,:)=M(i)%InvEE(l,:)/(2*l+1)
                    M(i)%InvEB(l,:)=M(i)%InvEB(l,:)/(2*l+1)
                end do
            end if
        end do
    else
        call Matrix_InverseArrayMPI(DummyArr,nmat)
    end if

    end subroutine PseudoCl_GetCouplingInversesArr

    subroutine TCovMatArray_Free(C)
    Type(TCovmatArray) :: C
    integer i,j, err

    do i=1, C%ncl
        do j=1, C%ncl
            deallocate(C%Cov(i,j)%C, stat=err)
        end do
    end do
    deallocate(C%Cov)

    end  subroutine TCovMatArray_Free

    subroutine TCouplingMatrix_Write(M, fname,i)
    use AMLUtils
    Type(TCouplingMatrix) :: M
    character(LEN=*), intent(in) :: fname
    integer, intent(in) :: i
    logical B1,B2

    call CreateFile(fname, i ,'UNFORMATTED')
    B1 = associated(M%T)
    B2 = associated(M%InvT)

    write (i) M%lmin,M%lmax,M%nl,M%has_pol,B1,B2

    if (B1) then
        write(i) M%T
        if (M%has_pol) then
            write(i) M%X
            write(i) M%EE
            write(i) M%EB
        end if
    end if
    if (B2) then
        write(i) M%InvT
        if (M%has_pol) then
            write(i) M%InvX
            write(i) M%InvEE
            write(i) M%InvEB
        end if
    end if

    Close(i)

    end subroutine TCouplingMatrix_Write

    subroutine TCouplingMatrix_Read(M, fname,i)
    use AMLUtils
    Type(TCouplingMatrix) :: M
    character(LEN=*), intent(in) :: fname
    integer, intent(in) :: i
    integer lmin,lmax
    logical B1,B2

    call TCouplingMatrix_Free(M)
    call OpenFile(fname, i ,'UNFORMATTED')

    read (i) lmin,lmax,M%nl,M%has_pol,B1,B2
    M%lmin=lmin
    M%lmax=lmax
    if (B1) then
        allocate(M%T(lmin:lmax,lmin:lmax))
        read(i) M%T
        if (M%has_pol) then
            allocate(M%X(lmin:lmax,lmin:lmax))
            allocate(M%EE(lmin:lmax,lmin:lmax))
            allocate(M%EB(lmin:lmax,lmin:lmax))

            read(i) M%X
            read(i) M%EE
            read(i) M%EB
        end if
    else
        nullify(M%T,M%X,M%EE,M%EB)
    end if
    if (B2) then
        allocate(M%InvT(lmin:lmax,lmin:lmax))
        read(i) M%InvT
        if (M%has_pol) then
            allocate(M%InvX(lmin:lmax,lmin:lmax))
            allocate(M%InvEE(lmin:lmax,lmin:lmax))
            allocate(M%InvEB(lmin:lmax,lmin:lmax))

            read(i) M%InvX
            read(i) M%InvEE
            read(i) M%InvEB
        end if
    else
        nullify(M%InvT,M%InvEE,M%InvEB,M%InvT)
    end if

    Close(i)

    end subroutine TCouplingMatrix_Read


    subroutine SetArrayInterlacedSym(MMixed,Mout, nthreads)
    real(dp) ::  Mout(:,:)
    real(dp), intent(in) :: MMixed(:,:,:)
    integer, intent(in) :: nthreads
    integer id, i, j,colix

    do id = 0, nthreads-1
        colix=0
        do i = 1+id, size(Mout,DIM=2), nthreads
            colix=colix+1
            Mout(1:i,i) = MMixed(1:i,colix,id+1)
        end do
    end do
    do i = 1, size(Mout,DIM=2)
        do j=1, i-1
            Mout(i,j) = Mout(j,i)
        end do
    end do

    end subroutine SetArrayInterlacedSym


    subroutine AssignSym(Mout, Minp)
    real(dp) :: Mout(:,:), Minp(:,:)
    integer i,j

    do i = 1, size(Mout,DIM=2)
        do j=1, i
            Mout(i,j) = Minp(j,i)
            Mout(j,i) = Minp(j,i)
        end do
    end do
    end subroutine AssignSym


    subroutine PseudoCl_GetCouplingMatrixArr(M_in, P, inlmin, inlmax, indopol, inncl, inpolweights, inpolspin)
    !Get coupling matrix from power spectrum of the weight function
    !For the moment assume weight same for all maps
    !Returns the symmetric form M_{l1 l2}/(2*l2+1)
    use AMLUtils
    use spinalm_tools
    Type(TCouplingMatrix),dimension(:), target :: M_in
    Type(TCouplingMatrix),dimension(:), pointer :: M

    Type(HealpixPower),dimension(:,:) :: P
    logical, intent(in) :: indopol
    integer, intent(in) :: inlmin,inlmax, inncl
    logical, intent(in) :: inpolweights
    integer, intent(in), optional :: inpolspin
    integer ncl,lmax, npolweights
    logical dopol
    real(dp), dimension(:,:,:), allocatable :: W
    real(dp), dimension(:,:), allocatable ::  tmp
    real(dp), dimension(:,:,:), allocatable ::  blocks

    real(dp), dimension(:), allocatable :: threejj0, threejj2, sthreejj0, sthreejj2,sthreejj20
    integer :: lmin, nl, l1, l2, lplus, lminus,  Plmax
    integer clix
    integer MpiID, MpiSize
    integer colix, ierr, ncols, polspin
    integer id, params(6), TPix, PPix
    double precision IniTime
    !Assume loads of memory

    IniTime = GeteTime()

    call GetMpiStat(MpiID, MpiSize)

    params(1) = inlmax
    params(2) = inncl
    if (indopol) then
        params(3)=1
    else
        params(3)=0
    end if
    params(4) = inlmin
    if (inpolweights) then
        params(5) = 3
    else
        params(5)= 1
    end if
    if (.not. present(inpolspin)) then
        params(6)=2
    else
        params(6)=inpolspin
    end if

#ifdef MPIPIX
    call MPI_BCAST(params,SIze(params),MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

    lmax=params(1)
    ncl=params(2)
    dopol = params(3)==1
    lmin=params(4)
    npolweights = params(5)
    polspin = params(6)
    if (npolweights==1) then
        TPix=1
        PPix=1
    else
        TPix=2
        PPix=3
    end if
    if (MpiID>0) then
        allocate(M(ncl))
    else
        M => M_in
    end if

    allocate(W(0:2*lmax, ncl, npolweights))
    W=0
    allocate(threejj0(0:2*lmax), sthreejj0(0:2*lmax))
    if (dopol) then
        allocate(threejj2(0:2*lmax), sthreejj2(0:2*lmax), sthreejj20(0:2*lmax))
    end if

    ncols = (lmax - lmin + MpiSize) / MpiSize !Max number

    do clix = 1, ncl
        if (MpiID==0) then
            ! print *,'free'
            call TCouplingMatrix_Nullify(M(clix)) !Currently seg faults if try to free, usual intel bug
            if (lmax > P(clix,1)%lmax) call MpiStop('PseudoCl_GetCouplingMatrix: lmax in power spectrum too low')
            if (lmax*2 > P(clix,1)%lmax) write(*,*) 'Warning: PseudoCl_GetCouplingMatrix really needs 2*lmax '
        end if
        M(clix)%has_pol = dopol
        M(clix)%lmin = lmin
        M(clix)%lmax = lmax
        M(clix)%nl = lmax-lmin +1
        allocate(M(clix)%T(lmin:lmax,ncols))
        M(clix)%T=0

        if (dopol) then
            allocate(M(clix)%X(lmin:lmax,ncols))
            M(clix)%X=0
            allocate(M(clix)%EE(lmin:lmax,ncols))
            M(clix)%EE=0
            allocate(M(clix)%EB(lmin:lmax,ncols))
            M(clix)%EB=0
        end if

        if (MpiID==0) then
            do l1 = 0,  min(P(clix,1)%lmax,2*lmax)
                W(l1,clix,1) = (2*l1+1)*P(clix,1)%Cl(l1,C_T)/(4*pi)
                if (npolweights>1) then
                    W(l1,clix,2) = (2*l1+1)*P(clix,2)%Cl(l1,C_T)/(4*pi)
                    W(l1,clix,3) = (2*l1+1)*P(clix,3)%Cl(l1,C_T)/(4*pi)
                end if
            end do
        end if
    end do


#ifdef MPIPIX
    call MPI_BCAST(W,SIze(W),MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
    colix =0
    do l1 = lmin + MpiID, lmax, MpiSize
        colix = colix+1
        if (colix > ncols) stop 'colix > ncols'
        !            print *,MpiID, 'colix',colix
        do l2 = lmin, l1
            lplus =  l1+l2
            lminus = abs(l1-l2)

            call GetThreeJs(threejj0(lminus:),l1,l2,0,0)
            sthreejj0(lminus:lplus) = threejj0(lminus:lplus)**2

            if (dopol) then
                !note that lminus is correct, want max(abs(l1-l2),abs(m1)) where m1=0 here
                call GetThreeJs(threejj2(lminus:),l1,l2,-polspin,polspin)
                sthreejj2(lminus:lplus) = threejj2(lminus:lplus)**2
                sthreejj20(lminus:lplus:2) = threejj0(lminus:lplus:2)*threejj2(lminus:lplus:2)
            end if

            do clix = 1, ncl
                M(clix)%T(l2,colix) = sum(W(lminus:lplus,clix,1)*sthreejj0(lminus:lplus))

                if (dopol) then
                    M(clix)%X(l2,colix) = sum(W(lminus:lplus:2,clix,TPix)*sthreejj20(lminus:lplus:2))
                    M(clix)%EE(l2,colix) = sum(W(lminus:lplus:2,clix,PPix)*sthreejj2(lminus:lplus:2))
                    M(clix)%EB(l2,colix) = sum(W(lminus+1:lplus:2,clix,PPix)*sthreejj2(lminus+1:lplus:2))
                end if
            end do
        end do
    end do

    deallocate(W)
    deallocate(threejj0,sthreejj0)
    if (dopol)  deallocate(threejj2,sthreejj2,sthreejj20)

#ifdef MPIPIX
    allocate(Blocks(lmin:lmax, ncols,MpiSize))
#else
    allocate(tmp(lmin:lmax,lmin:lmax))
#endif
    do clix = 1, ncl
#ifdef MPIPIX
        call MPI_GATHER(M(clix)%T,size(M(clix)%T),MPI_DOUBLE_PRECISION,blocks,size(M(clix)%T),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        deallocate(M(clix)%T)
        if (MpiID==0) then
            allocate(M(clix)%T(lmin:lmax,lmin:lmax))
            call SetArrayInterlacedSym(blocks, M(clix)%T,MpiSize)
        end if
        if (dopol) then
            call MPI_GATHER(M(clix)%X,size(M(clix)%X),MPI_DOUBLE_PRECISION,Blocks,size(M(clix)%X),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            deallocate(M(clix)%X)
            if (MpiID==0) then
                allocate(M(clix)%X(lmin:lmax,lmin:lmax))
                call SetArrayInterlacedSym(blocks, M(clix)%X,MpiSize)
            end if
            call MPI_GATHER(M(clix)%EE,size(M(clix)%EE),MPI_DOUBLE_PRECISION,Blocks,size(M(clix)%EE),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            deallocate(M(clix)%EE)
            if (MpiID==0) then
                allocate(M(clix)%EE(lmin:lmax,lmin:lmax))
                call SetArrayInterlacedSym(Blocks, M(clix)%EE,MpiSize)
            end if
            call MPI_GATHER(M(clix)%EB,size(M(clix)%EB),MPI_DOUBLE_PRECISION,Blocks,size(M(clix)%EB),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            deallocate(M(clix)%EB)
            if (MpiID==0) then
                allocate(M(clix)%EB(lmin:lmax,lmin:lmax))
                call SetArrayInterlacedSym(Blocks, M(clix)%EB,MpiSize)
            end if
        end if
#else
        tmp = M(clix)%T
        deallocate(M(clix)%T)
        allocate(M(clix)%T(lmin:lmax,lmin:Lmax))
        call AssignSym(M(clix)%T,tmp)

        !    M(clix)%T = tmp
        if (dopol) then
            tmp = M(clix)%X
            deallocate(M(clix)%X)
            allocate(M(clix)%X(lmin:lmax,lmin:Lmax))
            call AssignSym(M(clix)%X,tmp)

            tmp = M(clix)%EE
            deallocate(M(clix)%EE)
            allocate(M(clix)%EE(lmin:lmax,lmin:Lmax))
            call AssignSym(M(clix)%EE, tmp)

            tmp = M(clix)%EB
            deallocate(M(clix)%EB)
            allocate(M(clix)%EB(lmin:lmax,lmin:Lmax))
            call AssignSym(M(clix)%EB, tmp)
        end if
#endif
    end do

#ifdef MPIPIX
    if (MpiID>0) deallocate(M)
    deallocate(Blocks)
#else
    deallocate(tmp)
#endif
    if (MpiID==0) print *,'coupling time:', GeteTime() -  IniTime

    end subroutine PseudoCl_GetCouplingMatrixArr


    subroutine PseudoCl_GetCHat(M, PCls, Phat)
    !We only get TE where map index of E is >= that of T
    use MatrixUtils
    integer lmax
    Type(TCouplingMatrix) :: M
    Type(HealpixAllCl) :: Pcls
    Type(HealpixPower) :: Phat
    real(dm) tmp1(M%lmin:M%lmax),tmp2(M%lmin:M%lmax),tmpB(M%lmin:M%lmax)

    if (.not. associated(M%InvT)) then
        Stop 'PseudoCl_GetCHat: need inverses'
    end if

    call HealpixPower_Init(Phat,M%lmax, M%has_pol .and. size(PCls%Cl,2)+size(PCls%Cl,3)>2)
    tmp1=PCls%Cl(M%lmin:M%lmax,1,1)
    call Matrix_MulVec(M%InvT,tmp1,tmp2)
    Phat%Cl(M%lmin:M%lmax,C_T) = tmp2

    if (M%has_pol) then
        if (size(PCls%Cl,2)>1) then
            tmp1=PCls%Cl(M%lmin:M%lmax,2,1)
            call Matrix_MulVec(M%InvX,tmp1,tmp2)
            Phat%Cl(M%lmin:M%lmax,C_C)=tmp2
            if (size(PCls%Cl,3)>1) then
                tmp1 = PCls%Cl(M%lmin:M%lmax,2,2)
                tmpB = PCls%Cl(M%lmin:M%lmax,3,3)
                call Matrix_MulVec(M%InvEE,tmp1,tmp2)
                call Matrix_MulVec(M%InvEB,tmpB,tmp2,1._dm,1._dm)
                Phat%Cl(M%lmin:M%lmax,C_E)=tmp2

                call Matrix_MulVec(M%InvEE,tmpB,tmp2)
                call Matrix_MulVec(M%InvEB,tmp1,tmp2,1._dm,1._dm)
                Phat%Cl(M%lmin:M%lmax,C_B)=tmp2
            end if
        end if
    end if

    end subroutine PseudoCl_GetCHat

    subroutine PseudoCl_GetSimpleCovariance(M, Cov, PFid, vec_size, NoiseP, fwhm)
    !Simple 1/0 mask with noise added into the theory Cl
    !Not used
    use MatrixUtils
    Type(TCouplingMatrix) :: M
    Type(HealpixPower) :: PFid, NoiseP
    Type(TCovMat) :: Cov
    integer, intent(in) :: vec_size
    real(dp), intent(in) :: fwhm
    integer nl,i, status
    real fac
    real rootT(M%lmin:M%lmax), rootE(M%lmin:M%lmax), rootB(M%lmin:M%lmax)
    real(dp), dimension(:,:), allocatable :: Minv,tmp

    nl = M%lmax-M%lmin+1

    if (.not. associated(M%InvT)) then
        call MpiStop('Need inverses')
    end if

    PFid%Cl = PFid%Cl + NoiseP%Cl

    deallocate(Cov%C,stat=status)
    allocate(Cov%C(nl*vec_size,nl*vec_size))
    Cov%C = 0

    allocate(Minv(nl*vec_size,nl*vec_size))

    rootT = sqrt(PFid%Cl(M%lmin:M%lmax,C_T))
    rootE = sqrt(PFid%Cl(M%lmin:M%lmax,C_E))
    rootB = sqrt(PFid%Cl(M%lmin:M%lmax,C_B))

    print *,'doing T'

    do i=1,nl
        Cov%C(1:nl,i) =  2*PFid%Cl(M%lmin+i-1,C_T)*M%T(:,M%lmin+i-1)*PFid%Cl(M%lmin:M%lmax,C_T)
    end do

    print *,'doing Minv'

    Minv(1:nl,1:nl) = M%invT

    print *,'doing pol'
    if (M%has_pol) then
        !ordering is T X E B
        if (vec_size /= 3 .and. vec_size /= 4) &
        stop 'PseudoCl_GetCovariance: vec_size must be 3 (TT, TE, EE) or 4 (+BB)'

        do i=1,nl
            !XX
            fac = rootE(M%lmin+i-1)*rootT(M%lmin+i-1)

            Cov%C(nl+1:2*nl,nl+i) = fac*M%X(:,M%lmin+i-1)*rootE(M%lmin:M%lmax)*rootT(M%lmin:M%lmax) &
            + 2*PFid%Cl(M%lmin+i-1,C_C)*M%T(:,M%lmin+i-1)*PFid%Cl(M%lmin:M%lmax,C_C)

            !EE
            Cov%C(2*nl+1:3*nl,2*nl+i) = 2*PFid%Cl(M%lmin+i-1,C_E)*M%EE(:,M%lmin+i-1)*PFid%Cl(M%lmin:M%lmax,C_E) &
            +2*PFid%Cl(M%lmin+i-1,C_B)*M%EB(:,M%lmin+i-1)*PFid%Cl(M%lmin:M%lmax,C_B)

            !TX
            Cov%C(1:nl,nl+i) = rootT(M%lmin+i-1)*rootT(M%lmin:M%lmax)*&
            (PFid%Cl(M%lmin:M%lmax,C_C)+PFid%Cl(M%lmin+i-1,C_C))*M%T(:,M%lmin+i-1)
            Cov%C(nl+i,1:nl) = Cov%C(1:nl,2*nl+i)

            !TE
            Cov%C(1:nl,2*nl+i) = 2*PFid%Cl(M%lmin+i-1,C_C)*M%T(:,M%lmin+i-1)*PFid%Cl(M%lmin:M%lmax,C_C)
            Cov%C(2*nl+i,1:nl) = Cov%C(1:nl,2*nl+i)

            !XE
            Cov%C(nl+1:2*nl,2*nl+i) = rootT(M%lmin+i-1)*rootT(M%lmin:M%lmax)* &
            (PFid%Cl(M%lmin:M%lmax,C_C)+PFid%Cl(M%lmin+i-1,C_C))*M%T(:,M%lmin+i-1)
            Cov%C(2*nl+i,nl+1:2*nl) = Cov%C(2*nl+i,nl+1:2*nl)


            if (vec_size>3) then
                !BB
                Cov%C(3*nl+1:4*nl,3*nl+i) = 2*PFid%Cl(M%lmin+i-1,C_E)*M%EB(:,M%lmin+i-1)*PFid%Cl(M%lmin:M%lmax,C_E)+&
2               *PFid%Cl(M%lmin+i-1,C_B)*M%EE(:,M%lmin+i-1)*PFid%Cl(M%lmin:M%lmax,C_B)
                !EB
                Cov%C(3*nl+1:4*nl,2*nl+i) =(rootE(M%lmin+i-1)*rootE(M%lmin:M%lmax)&
                +rootB(M%lmin+i-1)*rootB(M%lmin:M%lmax))**2/2*M%EB(:,M%lmin+i-1)
                Cov%C(2*nl+i,3*nl+1:4*nl) = Cov%C(3*nl+1:4*nl,2*nl+i)
            end if
        end do

        Minv(nl+1:2*nl,nl+1:2*nl) = M%invX
        Minv(2*nl+1:3*nl,2*nl+1:3*nl) = M%invEE

        if (vec_size>3) then
            Minv(3*nl+1:4*nl,3*nl+1:4*nl) = M%invEE
            Minv(2*nl+1:3*nl,3*nl+1:4*nl) = M%invEB
            Minv(3*nl+1:4*nl,2*nl+1:3*nl) = M%invEB
        end if
    else
        if (vec_size/=1) stop 'PseudoCl_GetCovariance: vec_size must be 1 for no pol'
    end if

    print *,'doing mult'
    allocate(tmp(nl*vec_size,nl*vec_size))
    call Matrix_Mult(Minv,Cov%C,tmp)
    call Matrix_Mult_NT(tmp,Minv,Cov%C)
    deallocate(tmp)

    deallocate(Minv)
    print *,'done cov'


    end subroutine PseudoCl_GetSimpleCovariance

    function sym_ix(n,ax,ay)
    integer sym_ix
    integer, intent(in) :: n,ax,ay
    integer ix, x,y
    integer x2,y2
    !there's easy analytical result for this too
    if (ay>ax) then
        x=ay
        y=ax
    else
        x=ax
        y=ay
    end if
    ix=0
    do x2=1,n
        do y2=1,x2
            ix=ix+1
            if (x2==x .and. y2==y) then
                sym_ix =ix
                return
            end if
        end do
    end do
    write (*,*) 'n,x,y = ', n,x,y
    call MpiStop('sym_ix: error')

    end function sym_ix

    subroutine TCovMatPolArray_Free(CArr)
    Type(TCovMatPolArray) :: CArr
    integer i,j,polx,poly, status

    do polx=1,CArr%vec_size
        do poly =1, polx
            call TCovMatArray_Free(CArr%Pol(polx,poly))
        end do
    end do
    deallocate(CArr%Pol)

    end subroutine TCovMatPolArray_Free

    subroutine EffectiveCovPower(P)
    Type(HealpixPower) :: P
    real(dp), allocatable :: trueP(:)
    real(dp) :: win, norm
    integer l, l2, avrange,lmin,lmax

    print *,'using effective smoothed C_l'
    lmin = 0
    lmax= P%lmax

    ! get effective average C_l for covariance; try just using first one where cosmic var most important
    ! this tries to cure underestimation of cov because C_l falling rapidly, esp by first trough
    !call HealpixPower_Write(SmoothedP(channel),'outfiles/tmpin'//trim(INtToStr(channel)))
    avrange=40
    allocate(trueP(lmin:lmax))
    trueP(lmin:Lmax) = P%Cl(lmin:lmax,C_T)
    do l= max(200,avrange+lmin), lmax-avrange
        norm = 0
        P%Cl(l,C_T)=0
        do l2=l-avrange,min(lmax,l+avrange)
            win = exp(-real(l2-l)**2/(2*25**2))    !l2^2 and P^2 is OKish
            norm = norm + win
            P%Cl(l,C_T)= P%Cl(l,C_T) + trueP(l2)*win
        end do
        P%Cl(l,C_T)=(P%Cl(l,C_T)/norm)
    end do
    ! call HealpixPower_Write(P, 'outfiles/tmp')
    deallocate(trueP)

    end  subroutine EffectiveCovPower

    function sym_ix_check(chk, n,i,j) result(res)
    logical chk(:)
    integer, intent(in) :: n,i,j
    integer res

    res = sym_ix(n,i,j)
    chk(res) = .true.

    end function sym_ix_check

    subroutine PseudoCl_GetFullCovariance(Coupler, Xi, CArr, PFid, vec_size,  Channels, nweights, &
    incnoise,incsignal, frac_pt_src_error, pol_weights, step, nstep)
    !General result
    !At l>30 seems to be very accurate to use just Xi%T rather than Xi%E (and consistent at this approx anyway)

    use MatrixUtils
    Type(TCovMatPolArray), intent(inout), target :: CArr
    logical, intent(in) :: incnoise, incsignal
    integer, intent(in) :: step, nstep !to split up so don't blow memory
    real, intent(in) :: frac_pt_src_error
    logical, intent(in) :: pol_weights
    Type(TCouplingMatrix), intent(in) :: Xi(:)  !2*ncl(2*ncl+1)/2 vector
    Type(TCouplingMatrix), intent(in) :: Coupler(:)
    Type(HealpixPower), intent(in) :: PFid
    type(TChannel) :: Channels(:)
    Type(TCovMat), pointer :: C, EE, EB, BE, BB
    integer, intent(in) :: vec_size
    integer, intent(in) :: nweights
    real(dp) ::  ENoiseFac_11,ENoiseFac_22, ENoiseFac_12,ENoiseFac_21,ENoiseFac2
    integer aix,nl, ncl_tot, ncl, l, j, status
    real fac, crossNoiseFac
    real, dimension(:,:), allocatable :: rootT, rootE, rootB
    real, dimension(:), allocatable :: rootT_11,RootT_12,rootT_21,rootT_22
    real, dimension(:), allocatable :: BeamAvTE, rootE_11,RootE_12,rootE_21,rootE_22
    real, dimension(:), allocatable :: rootB_11,RootB_12,rootB_21,rootB_22

    real(dp), dimension(:,:), allocatable :: tmp, tmp11,tmp12,tmp21,tmp22
    integer, parameter :: w2=1, wTT=2, w2T =3
    integer i1,i2,j1,j2,ix, x,y, ix2, x2,y2, lmin, lmax
    integer ix_couple,ix2_couple
    integer s_11_s_22,s_12_s_21, n_11_s_22,n_12_s_21, s_11_n_22,s_12_n_21, n_11_n_22,n_12_n_21
    integer PP_s_11_s_22, PP_s_12_s_21, PP_n_11_s_22, PP_n_12_s_21, PP_s_11_n_22, PP_s_12_n_21, PP_n_11_n_22, PP_n_12_n_21
    integer PT_s_11_s_22, PT_s_12_s_21, PT_n_11_s_22, PT_n_12_s_21, PT_s_11_n_22, PT_s_12_n_21, PT_n_11_n_22, PT_n_12_n_21
    integer MM_s_11_s_22, MM_s_12_s_21
    integer MT_s_11_n_22, MT_s_12_n_21, MP_s_12_n_21, PM_n_11_s_22
    integer MT_s_11_s_22, MT_s_12_s_21, PM_s_11_s_22, PM_s_12_s_21
    logical dopol
    integer polx,poly, nchannels, channel
    double precision STime
    Type(HealpixPower), allocatable :: DummyNoise(:)
    Type(HealpixPower), allocatable :: SmoothedP(:)
    Type(HealpixPower) :: AvP

    Type(TBeam) :: AvBeam
    integer chanx, chanx2,chany,chany2
    integer nvarmaps, nmaskmaps
    integer ix_11,ix_12, ix_21,ix_22
    logical, allocatable :: chk(:)

    nchannels = size(Channels)
    lmin = Coupler(1)%lmin
    lmax = Coupler(1)%lmax
    ncl = (nweights*(nweights+1))/2
    ncl_tot = (nweights*nchannels*(nweights*nchannels+1))/2

    print *,'PseudoCl_GetFullCovariance: channels = ',nchannels

    dopol= Coupler(1)%has_pol
    nl = lmax-lmin+1

    if (pol_weights) then
        nvarmaps= (2*nchannels+3)*ncl
        nmaskmaps=(nchannels+1)*ncl
    else
        nvarmaps= (nchannels+1)*ncl
        nmaskmaps =0
    end if
    allocate(chk(size(Xi)))
    chk = .false.

    if (step /= 0) then
        allocate(rootT_11(lmin:lmax),rootT_21(lmin:lmax),rootT_12(lmin:lmax),rootT_22(lmin:lmax))
        allocate(rootT(lmin:lmax,nchannels))
        if (dopol) then
            allocate(rootE_11(lmin:lmax),rootE_21(lmin:lmax),rootE_12(lmin:lmax),rootE_22(lmin:lmax))
            allocate(rootB_11(lmin:lmax),rootB_21(lmin:lmax),rootB_12(lmin:lmax),rootB_22(lmin:lmax))
            allocate(RootE(lmin:lmax,nchannels),RootB(lmin:lmax,nchannels))
            allocate(BeamAvTE(lmin:lmax))
        end if
        call HealpixPower_Nullify(AvP)


        print*,'Noise Scale= ', NoiseScale
        crossNoiseFac = NoiseScale**2*cross_noise_scale ! for the N-N terms, the scaling -can be increased for cross-spectra

        if (.not. associated(Coupler(1)%InvT))  call MpiStop('Need inverses')

        print *,'getting smoothed Cls'

        allocate(SmoothedP(nchannels))
        do channel = 1, nchannels
            call healpixPower_Nullify(SmoothedP(channel))
            call HealpixPower_Assign(SmoothedP(channel), PFid)
            if (incsignal) then
                SmoothedP(channel)%Cl(:,C_T) = SmoothedP(channel)%Cl(:,C_T) + Channels(channel)%PtSrcA

                call TBeam_PowerSmooth(Channels(channel)%Beam, SmoothedP(channel), -1)

                call EffectiveCovPower(SmoothedP(channel))
            else
                SmoothedP(channel)%Cl=0
            end if

            rootT(:,channel) = sqrt(SmoothedP(channel)%Cl(lmin:lmax,C_T))
            if (dopol) then
                !         if (nchannels >1) call MpiStop('not done multi-channel pol')
                rootE(:,channel) = sqrt(SmoothedP(channel)%Cl(lmin:lmax,C_E))
                rootB(:,channel) = sqrt(SmoothedP(channel)%Cl(lmin:lmax,C_B))
            end if
        end do


        if (step==1) then
            print *,'allocating'
            CArr%vec_size = vec_size
            allocate(CArr%Pol(vec_size,vec_size))
            do polx=1,vec_size
                do poly =1, polx
                    allocate(CArr%Pol(polx,poly)%Cov(ncl_tot,ncl_tot))
                    CArr%Pol(polx,poly)%ncl = ncl_tot
                end do
            end do
        end if

        allocate(DummyNoise(nchannels))
        do channel = 1, nchannels
            call HealpixPower_Nullify(DummyNoise(channel))
            call HealpixPower_Init(DummyNoise(channel), lmax, .false.)
            DummyNoise(channel)%Cl = 1
            call TBeam_PowerSmooth(Channels(channel)%Beam, DummyNoise(channel),+1)
        end do

        allocate(tmp(lmin:lmax,lmin:lmax))

    end if !step/=0


    STime = GeteTime()
    print *,'Doing covariance main loops'

    ix=0
    i1=0
    do chanx = 1, nchannels
        do x=1,nweights
            i1 = i1+1
            i2=0
            do chany = 1, nchannels
                do y=1,nweights
                    i2 = i2+1
                    if (i2 > i1) cycle

                    ix=ix+1

                    ix_couple= sym_ix(nweights, x,y)

                    ix2=0
                    j1 =0
                    do chanx2 = 1, nchannels
                        do x2=1,nweights
                            j1 = j1+1
                            j2=0

                            do chany2 = 1, nchannels
                                do y2=1,nweights
                                    j2 = j2+1

                                    if (j2 > j1) cycle
                                    ix2 = ix2+1

                                    if (step/=0) then
                                        rootT_11 = sqrt(rootT(:,chanx)*rootT(:,chanx2))
                                        rootT_12 = sqrt(rootT(:,chanx)*rootT(:,chany2))
                                        rootT_21 = sqrt(rootT(:,chany)*rootT(:,chanx2))
                                        rootT_22 = sqrt(rootT(:,chany)*rootT(:,chany2))

                                        !     call HealpixPower_Assign(AvP, PFid)
                                        !     if (incsignal) then
                                        !      AvP%Cl(:,1) = AvP%Cl(:,1) + &
                                        !     ( Channels(chanx)%PtSrcA*Channels(chany)%PtSrcA*Channels(chanx2)%PtSrcA*Channels(chany2)%PtSrcA)**0.25
                                        !     else
                                        !      AvP%Cl=0
                                        !     end if
                                        !     call TBeam_Assign(AvBeam, Channels(chanx)%Beam)
                                        !     AvBeam%Beam = (Channels(chanx)%Beam%beam*Channels(chany)%Beam%beam* &
                                        !                    Channels(chanx2)%Beam%beam*Channels(chany2)%Beam%beam)**0.25_dp
                                        !     call TBeam_PowerSmooth(AvBeam, AvP, -1)
                                        !!call HealpixPower_Assign(AvP, PFid)
                                        call HealpixPower_Assign(AvP, SmoothedP(chanx))

                                        !Beware overflow if C_l numbers small
                                        AvP%Cl(:,C_T) = (real(SmoothedP(chanx)%Cl(:,C_T),dp)*SmoothedP(chany)%Cl(:,C_T)* &
                                        real(SmoothedP(chanx2)%Cl(:,C_T),dp)*SmoothedP(chany2)%Cl(:,C_T) ) ** 0.25_dp


                                        if (dopol) then
                                            !Added AvP Dec 2010
                                            AvP%Cl(:,C_E) = (real(SmoothedP(chanx)%Cl(:,C_E),dp)*SmoothedP(chany)%Cl(:,C_E)* &
                                            real(SmoothedP(chanx2)%Cl(:,C_E),dp)*SmoothedP(chany2)%Cl(:,C_E) ) ** 0.25_dp
                                            AvP%Cl(:,C_B) = (real(SmoothedP(chanx)%Cl(:,C_B),dp)*SmoothedP(chany)%Cl(:,C_B)* &
                                            real(SmoothedP(chanx2)%Cl(:,C_B),dp)*SmoothedP(chany2)%Cl(:,C_B) ) ** 0.25_dp
                                            where (PFid%Cl(:,C_T)>0)
                                                AvP%Cl(:,C_C)= PFid%Cl(:,C_C)* AvP%Cl(:,C_T)/PFid%Cl(:,C_T)
                                            end where

                                            rootE_11 = sqrt(rootE(:,chanx)*rootE(:,chanx2))
                                            rootE_12 = sqrt(rootE(:,chanx)*rootE(:,chany2))
                                            rootE_21 = sqrt(rootE(:,chany)*rootE(:,chanx2))
                                            rootE_22 = sqrt(rootE(:,chany)*rootE(:,chany2))
                                            rootB_11 = sqrt(rootB(:,chanx)*rootB(:,chanx2))
                                            rootB_12 = sqrt(rootB(:,chanx)*rootB(:,chany2))
                                            rootB_21 = sqrt(rootB(:,chany)*rootB(:,chanx2))
                                            rootB_22 = sqrt(rootB(:,chany)*rootB(:,chany2))
                                            BeamAvTE = sqrt(AvP%Cl(lmin:lmax,C_T)*AvP%Cl(lmin:lmax,C_E))
                                        end if
                                    end if !step /= 0

                                    ix2_couple= sym_ix(nweights, x2,y2)

                                    ix_11 = sym_ix_check(chk,nweights, x,x2)
                                    ix_12 = sym_ix_check(chk,nweights, x,y2)
                                    ix_21 = sym_ix_check(chk,nweights, y,x2)
                                    ix_22 = sym_ix_check(chk,nweights, y,y2)

                                    !T-T
                                    s_11_s_22 = sym_ix_check(chk,nvarmaps,  ix_11, ix_22)
                                    s_12_s_21 = sym_ix_check(chk,nvarmaps,  ix_12, ix_21)

                                    PP_s_11_s_22 = sym_ix_check(chk,nvarmaps,  ix_11+nmaskmaps, ix_22+nmaskmaps)
                                    PP_s_12_s_21 = sym_ix_check(chk,nvarmaps,  ix_12+nmaskmaps, ix_21+nmaskmaps)

                                    PT_s_11_s_22 = sym_ix_check(chk,nvarmaps,  ix_11+nmaskmaps, ix_22)
                                    PT_s_12_s_21 = sym_ix_check(chk,nvarmaps,  ix_12+nmaskmaps, ix_21)

                                    MM_s_11_s_22 = sym_ix_check(chk,nvarmaps,  ix_11+2*nmaskmaps, ix_22+2*nmaskmaps)
                                    MM_s_12_s_21 = sym_ix_check(chk,nvarmaps,  ix_12+2*nmaskmaps, ix_21+2*nmaskmaps)

                                    MT_s_11_s_22 = sym_ix_check(chk,nvarmaps,  ix_11+2*nmaskmaps, ix_22+nmaskmaps)
                                    MT_s_12_s_21 = sym_ix_check(chk,nvarmaps,  ix_12+2*nmaskmaps, ix_21+nmaskmaps)

                                    PM_s_11_s_22 = sym_ix_check(chk,nvarmaps,  ix_11+nmaskmaps, ix_22+2*nmaskmaps)
                                    PM_s_12_s_21 = sym_ix_check(chk,nvarmaps,  ix_12+nmaskmaps, ix_21+2*nmaskmaps)


                                    !nweights-T
                                    n_11_s_22=0; n_12_s_21=0;s_11_n_22=0;s_12_n_21=0
                                    PT_n_11_s_22=0; PT_n_12_s_21=0;PT_s_11_n_22=0;PT_s_12_n_21=0
                                    PP_n_11_s_22=0; PP_n_12_s_21=0;PP_s_11_n_22=0;PP_s_12_n_21=0
                                    MT_s_11_n_22=0; MT_s_12_n_21=0;
                                    MP_s_12_n_21=0; PM_n_11_s_22=0;


                                    if (chanx==chanx2) n_11_s_22 = sym_ix_check(chk,nvarmaps, ix_11+ncl*chanx, ix_22)
                                    if (chanx==chany2) n_12_s_21 = sym_ix_check(chk,nvarmaps, ix_12+ncl*chanx, ix_21)
                                    if (chany==chany2) s_11_n_22 = sym_ix_check(chk,nvarmaps, ix_11, ix_22+ncl*chany)
                                    if (chany==chanx2) s_12_n_21 = sym_ix_check(chk,nvarmaps, ix_12, ix_21+ncl*chany)

                                    if (vec_size>1) then
                                        if (chanx==chanx2) PT_n_11_s_22 = sym_ix_check(chk,nvarmaps, ix_11+ncl*chanx+nmaskmaps, ix_22)
                                        if (chanx==chany2) PT_n_12_s_21 = sym_ix_check(chk,nvarmaps, ix_12+ncl*chanx+nmaskmaps, ix_21)
                                        if (chany==chany2) PT_s_11_n_22 = sym_ix_check(chk,nvarmaps, ix_11+nmaskmaps, ix_22+ncl*chany)
                                        if (chany==chanx2) PT_s_12_n_21 = sym_ix_check(chk,nvarmaps, ix_12+nmaskmaps, ix_21+ncl*chany)

                                        if (chanx==chanx2) PP_n_11_s_22 = sym_ix_check(chk,nvarmaps, ix_11+ncl*chanx+nmaskmaps, ix_22+nmaskmaps)
                                        if (chanx==chany2) PP_n_12_s_21 = sym_ix_check(chk,nvarmaps, ix_12+ncl*chanx+nmaskmaps, ix_21)
                                        if (chany==chany2) PP_s_11_n_22 = sym_ix_check(chk,nvarmaps, ix_11+nmaskmaps, ix_22+ncl*chany)
                                        if (chany==chanx2) PP_s_12_n_21 = sym_ix_check(chk,nvarmaps, ix_12+nmaskmaps, ix_21+ncl*chany)

                                        if (chany==chany2) MT_s_11_n_22 = sym_ix_check(chk,nvarmaps, ix_11+2*nmaskmaps, ix_22+ncl*chany)
                                        if (chany==chanx2) MT_s_12_n_21 = sym_ix_check(chk,nvarmaps, ix_12+2*nmaskmaps, ix_21+ncl*chany)

                                        if (chany==chanx2) MP_s_12_n_21 = sym_ix_check(chk,nvarmaps, ix_12+2*nmaskmaps, ix_21+ncl*chany+nmaskmaps)
                                        if (chanx==chanx2) PM_n_11_s_22 = sym_ix_check(chk,nvarmaps, ix_11+ncl*chanx+nmaskmaps, ix_22+2*nmaskmaps)
                                    end if

                                    !nweights-nweights
                                    n_11_n_22=0; n_12_n_21=0
                                    PP_n_11_n_22=0; PP_n_12_n_21=0
                                    PT_n_11_n_22 = 0
                                    if (chanx==chanx2 .and. chany==chany2) then
                                        n_11_n_22 = sym_ix_check(chk,nvarmaps, ix_11+ncl*chanx, ix_22+ncl*chany)
                                        PP_n_11_n_22 = sym_ix_check(chk,nvarmaps, ix_11+ncl*chanx+nmaskmaps, ix_22+ncl*chany+nmaskmaps)
                                        PT_n_11_n_22 = sym_ix_check(chk,nvarmaps, ix_11+ncl*chanx+nmaskmaps, ix_22+ncl*chany)
                                    end if
                                    if (chanx==chany2 .and. chany==chanx2) then
                                        n_12_n_21 = sym_ix_check(chk,nvarmaps, ix_12+ncl*chanx, ix_21+ncl*chany)
                                        PP_n_12_n_21 = sym_ix_check(chk,nvarmaps, ix_12+ncl*chanx+nmaskmaps, ix_21+ncl*chany+nmaskmaps)
                                    end if

                                    ENoiseFac_11=sqrt(Channels(chanx)%ENoiseFac*Channels(chanx2)%ENoiseFac)
                                    ENoiseFac_22=sqrt(Channels(chany)%ENoiseFac*Channels(chany2)%ENoiseFac)
                                    ENoiseFac_21=sqrt(Channels(chany)%ENoiseFac*Channels(chanx2)%ENoiseFac)
                                    ENoiseFac_12=sqrt(Channels(chanx)%ENoiseFac*Channels(chany2)%ENoiseFac)
                                    ENoiseFac2=ENoiseFac_11*ENoiseFac_22

                                    !assume pol has larger noise than T, so always use pol in order of weights e.g. X = E_2 T_1 not T_2 E_1

                                    if (step==0) cycle

                                    do polx=1,vec_size
                                        do poly=1,polx
                                            !           print *, 'polx, poly', polx, poly

                                            C =>  CArr%Pol(polx,poly)%Cov(ix,ix2)

                                            if (ix2>ix .and. polx==poly .or. polx==4 .and. poly<3) then
                                                !Is symmetric  = transpose(CArr%Pol(polx,poly)%Cov(ix2,ix)%C)
                                                nullify(C%C)
                                                cycle
                                            end if

                                            if (step==1) then
                                                allocate(C%C(lmin:lmax,lmin:lmax))
                                                C%C=0
                                            end if
                                            !order T X E B

                                            if (polx==1 .and. poly==1) then
                                                !TT
                                                do l=lmin,lmax
                                                    if (incnoise) then
                                                        if (n_11_s_22/=0) C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + &
                                                        NoiseScale*Xi(n_11_s_22)%T(:,l)*rootT_22(lmin:lmax)*rootT_22(l)
                                                        if (n_12_s_21/=0) C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + &
                                                        NoiseScale*Xi(n_12_s_21)%T(:,l)*rootT_21(lmin:lmax)*rootT_21(l)
                                                        if (s_11_n_22/=0) C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + &
                                                        NoiseScale*Xi(s_11_n_22)%T(:,l)*rootT_11(lmin:lmax)*rootT_11(l)
                                                        if (s_12_n_21/=0) C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + &
                                                        NoiseScale*Xi(s_12_n_21)%T(:,l)*rootT_12(lmin:lmax)*rootT_12(l)

                                                        if (n_12_n_21/=0) C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + crossNoiseFac*Xi(n_12_n_21)%T(:,l)
                                                        if (n_11_n_22/=0) C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + crossNoiseFac*Xi(n_11_n_22)%T(:,l)
                                                    end if
                                                    if (s_11_s_22/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    AvP%Cl(l,C_T)*(Xi(s_11_s_22)%T(:,l))*AvP%Cl(lmin:lmax,C_T)
                                                    if (s_12_s_21/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    AvP%Cl(l,C_T)*(Xi(s_12_s_21)%T(:,l))*AvP%Cl(lmin:lmax,C_T)
                                                    !    print *, l,C%C(l,l),AvP%Cl(l,C_T) ,Xi(s_12_s_21)%T(l,l)
                                                end do

                                                if (step==nstep) then
                                                    call Matrix_Mult(Coupler(ix_couple)%invT,C%C,tmp)
                                                    call Matrix_Mult_NT(tmp,Coupler(ix2_couple)%invT,C%C)
                                                end if
                                            else if (polx==2 .and. poly==2) then
                                                !XX = E_x T_y E_x2 T_y2 only (by assumption)

                                                do l=lmin,lmax
                                                    if (PT_s_11_s_22/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l)&
                                                    + BeamAvTE(l) * (Xi(PT_s_11_s_22)%T(:,l)*BeamAvTE(lmin:lmax))
                                                    if (MM_s_12_s_21/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l)&
                                                    +  AvP%Cl(l,C_C)* AvP%Cl(lmin:lmax,C_C)*Xi(MM_s_12_s_21)%T(:,l)
                                                    if (incnoise) then
                                                        if (PT_n_11_n_22/=0)  C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + &
                                                        crossNoiseFac*ENoiseFac_11*Xi(PT_n_11_n_22)%T(:,l)
                                                        if (PT_n_11_s_22/=0)  C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + &
                                                        NoiseScale*ENoiseFac_11*rootT_22(l)*rootT_22(lmin:lmax) * Xi(PT_n_11_s_22)%T(:,l)
                                                        if (PT_s_11_n_22/=0)  C%C(lmin:lmax,l)= C%C(lmin:lmax,l) + &
                                                        NoiseScale*rootE_11(l)*rootE_11(lmin:lmax) * Xi(PT_s_11_n_22)%T(:,l)
                                                    end if
                                                end do
                                                if (step==nstep) then
                                                    call Matrix_Mult(Coupler(ix_couple)%invX,C%C,tmp)
                                                    call Matrix_Mult_NT(tmp,Coupler(ix2_couple)%invX,C%C)
                                                end if
                                            else if (polx==3 .and. poly==3) then
                                                !EE
                                                do l=lmin,lmax
                                                    if (PP_s_11_s_22/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    AvP%Cl(l,C_E)*(Xi(PP_s_11_s_22)%T(:,l))*AvP%Cl(lmin:lmax,C_E)
                                                    if (PP_s_12_s_21/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    AvP%Cl(l,C_E)*(Xi(PP_s_12_s_21)%T(:,l))*AvP%Cl(lmin:lmax,C_E)
                                                    if (incnoise) then
                                                        if (PP_n_11_n_22/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l) + crossNoiseFac*EnoiseFac2*Xi(PP_n_11_n_22)%T(:,l)
                                                        if (PP_n_12_n_21/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l) + crossNoiseFac*EnoiseFac2*Xi(PP_n_12_n_21)%T(:,l)
                                                        if (PP_n_11_s_22/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l) + &
                                                        NoiseScale*ENoiseFac_11*rootE_22(l)*rootE_22(lmin:lmax)*Xi(PP_n_11_s_22)%T(:,l)
                                                        if (PP_n_12_s_21/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l) + &
                                                        NoiseScale*ENoiseFac_12*rootE_21(l)*rootE_21(lmin:lmax)*Xi(PP_n_12_s_21)%T(:,l)
                                                        if (PP_s_11_n_22/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l) + &
                                                        NoiseScale*ENoiseFac_22*rootE_11(l)*rootE_11(lmin:lmax)*Xi(PP_s_11_n_22)%T(:,l)
                                                        if (PP_s_12_n_21/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l) + &
                                                        NoiseScale*ENoiseFac_21*rootE_12(l)*rootE_12(lmin:lmax)*Xi(PP_s_12_n_21)%T(:,l)
                                                    end if
                                                end do

                                                if (vec_size==3 .and. step==nstep) then
                                                    !Rather good approx to drop E-B terms entirely for l>>1
                                                    call Matrix_Mult(Coupler(ix_couple)%invEE,C%C,tmp)
                                                    call Matrix_Mult_NT(tmp,Coupler(ix2_couple)%invEE,C%C)
                                                end if
                                            else if (polx==4 .and. poly==4) then
                                                !BB

                                                do l=lmin,lmax
                                                    if (PP_s_11_s_22/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    AvP%Cl(l,C_B)*(Xi(PP_s_11_s_22)%T(:,l))*AvP%Cl(lmin:lmax,C_B)
                                                    if (PP_s_12_s_21/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    AvP%Cl(l,C_B)*(Xi(PP_s_12_s_21)%T(:,l))*AvP%Cl(lmin:lmax,C_B)
                                                    !     + PFid%Cl(l,C_E)*(Xi(s_11_s_22)%EB(:,l)+Xi(s_12_s_21)%EB(:,l))*PFid%Cl(lmin:lmax,C_E) &  !Leakage, check
                                                    if (incnoise) then
                                                        if (n_11_n_22/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l)+ &
                                                        crossNoiseFac*EnoiseFac2*Xi(PP_n_11_n_22)%T(:,l)
                                                        if (n_12_n_21/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l)+ &
                                                        crossNoiseFac*EnoiseFac2*Xi(PP_n_12_n_21)%T(:,l)
                                                        if (n_11_s_22/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l)+ &
                                                        NoiseScale*ENoiseFac_11*rootB_22(l)*rootB_22(lmin:lmax)*Xi(PP_n_11_s_22)%T(:,l)
                                                        if (n_12_s_21/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l)+ &
                                                        NoiseScale*ENoiseFac_12*rootB_21(l)*rootB_21(lmin:lmax)*Xi(PP_n_12_s_21)%T(:,l)
                                                        if (s_11_n_22/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l)+ &
                                                        NoiseScale*ENoiseFac_22*rootB_11(l)*rootB_11(lmin:lmax)*Xi(PP_s_11_n_22)%T(:,l)
                                                        if (s_12_n_21/=0) C%C(lmin:lmax,l)=C%C(lmin:lmax,l)+ &
                                                        NoiseScale*ENoiseFac_21*rootB_21(l)*rootB_21(lmin:lmax)*Xi(PP_s_12_n_21)%T(:,l)
                                                    end if
                                                end do
                                            else if (polx==2 .and. poly==1) then
                                                !X T = ET TT
                                                do l=lmin,lmax
                                                    !Should be using correct beam-averaged C_C, C_T here...
                                                    if (MT_s_11_s_22/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    sqrt(AvP%Cl(l,C_T))* sqrt(AvP%Cl(lmin:lmax,C_T))*(AvP%Cl(lmin:lmax,C_C)+AvP%Cl(l,C_C))/2* &
                                                    (Xi(MT_s_12_s_21)%T(:,l))
                                                    if (MT_s_12_s_21/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    sqrt(AvP%Cl(l,C_T))* sqrt(AvP%Cl(lmin:lmax,C_T))*(AvP%Cl(lmin:lmax,C_C)+AvP%Cl(l,C_C))/2* &
                                                    (Xi(MT_s_12_s_21)%T(:,l))
                                                    if (incnoise) then
                                                        if (MT_s_11_n_22/=0)  C%C(lmin:lmax,l) =  C%C(lmin:lmax,l) + &
                                                        NoiseScale*(AvP%Cl(lmin:lmax,C_C)+AvP%Cl(l,C_C))/2* Xi(MT_s_11_n_22)%T(:,l)
                                                        if (MT_s_12_n_21/=0)  C%C(lmin:lmax,l) =  C%C(lmin:lmax,l) + &
                                                        NoiseScale*(AvP%Cl(lmin:lmax,C_C)+AvP%Cl(l,C_C))/2* Xi(MT_s_12_n_21)%T(:,l)
                                                    end if
                                                end do
                                                if (step==nstep) then
                                                    call Matrix_Mult(Coupler(ix_couple)%invX,C%C,tmp)
                                                    call Matrix_Mult_NT(tmp,Coupler(ix2_couple)%invT,C%C)
                                                end if
                                            else if (polx==3 .and. poly==1) then
                                                !ET = EE TT
                                                do l=lmin,lmax
                                                    if (MM_s_11_s_22/=0) C%C(lmin:lmax,l) =  C%C(lmin:lmax,l) + &
                                                    AvP%Cl(l,C_C)* AvP%Cl(lmin:lmax,C_C)*( Xi(MM_s_11_s_22)%T(:,l) )
                                                    if (MM_s_12_s_21/=0) C%C(lmin:lmax,l) =  C%C(lmin:lmax,l) + &
                                                    AvP%Cl(l,C_C)* AvP%Cl(lmin:lmax,C_C)*( Xi(MM_s_12_s_21)%T(:,l))
                                                end do
                                                if (step==nstep) then
                                                    call Matrix_Mult(Coupler(ix_couple)%invEE,C%C,tmp)
                                                    call Matrix_Mult_NT(tmp,Coupler(ix2_couple)%invT,C%C)
                                                end if
                                            else if (polx==3 .and. poly==2) then
                                                !EX = EE ET

                                                do l=lmin,lmax
                                                    if (PM_s_11_s_22/=0) C%C(lmin:lmax,l) =   C%C(lmin:lmax,l) + &
                                                    sqrt(AvP%Cl(l,C_E))*sqrt(AvP%Cl(lmin:lmax,C_E))*(AvP%Cl(lmin:lmax,C_C)+AvP%Cl(l,C_C) )/2* &
                                                    ( Xi(PM_s_11_s_22)%T(:,l))
                                                    if (PM_s_12_s_21/=0) C%C(lmin:lmax,l) =   C%C(lmin:lmax,l) + &
                                                    sqrt(AvP%Cl(l,C_E))*sqrt(AvP%Cl(lmin:lmax,C_E))*(AvP%Cl(lmin:lmax,C_C)+AvP%Cl(l,C_C) )/2* &
                                                    ( Xi(PM_s_12_s_21)%T(:,l))

                                                    if (incnoise) then
                                                        if (PM_n_11_s_22/=0)    C%C(lmin:lmax,l) =   C%C(lmin:lmax,l) + &
                                                        NoiseScale*ENoiseFac_11* (AvP%Cl(lmin:lmax,C_C)+AvP%Cl(l,C_C))/2* Xi(PM_n_11_s_22)%T(:,l)
                                                        if (MP_s_12_n_21/=0)    C%C(lmin:lmax,l) =   C%C(lmin:lmax,l) + &
                                                        NoiseScale*ENoiseFac_21* (AvP%Cl(lmin:lmax,C_C)+AvP%Cl(l,C_C))/2* Xi(MP_s_12_n_21)%T(:,l)
                                                    end if
                                                end do
                                                if (step==nstep) then
                                                    call Matrix_Mult(Coupler(ix_couple)%invEE,C%C,tmp)
                                                    call Matrix_Mult_NT(tmp,Coupler(ix2_couple)%invX,C%C)
                                                end if
                                            else if (polx==4 .and. poly==3) then
                                                call MpiSTop('full covariance, should check EB now using approx that EE=T')
                                                !EB
                                                do l=lmin,lmax
                                                    !Should be using right averages
                                                    if (PP_s_11_s_22/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    (sqrt(AvP%Cl(l,C_E))*sqrt(AvP%Cl(lmin:lmax,C_E)) &
                                                    +  sqrt(AvP%Cl(l,C_B))*sqrt(AvP%Cl(lmin:lmax,C_B)))**2 *&
                                                    ( Xi(PP_s_11_s_22)%EB(:,l) )/2  ! factor of 2 ? check

                                                    if (PP_s_12_s_21/=0) C%C(lmin:lmax,l) = C%C(lmin:lmax,l) + &
                                                    (sqrt(AvP%Cl(l,C_E))*sqrt(AvP%Cl(lmin:lmax,C_E)) &
                                                    +  sqrt(AvP%Cl(l,C_B))*sqrt(AvP%Cl(lmin:lmax,C_B)))**2 *&
                                                    ( Xi(PP_s_12_s_21)%EB(:,l) )/2  ! factor of 2 ? check

                                                    if (incnoise) then
                                                        if (PP_n_11_n_22/=0) C%C(lmin:lmax,l) =  C%C(lmin:lmax,l) + &
                                                        crossNoiseFac*EnoiseFac2*(Xi(PP_n_11_n_22)%EB(:,l))
                                                        if (PP_n_12_n_21/=0)  C%C(lmin:lmax,l) =  C%C(lmin:lmax,l) + &
                                                        crossNoiseFac*EnoiseFac2*(Xi(PP_n_12_n_21)%EB(:,l))
                                                        if (PP_n_11_s_22/=0) C%C(lmin:lmax,l) =  C%C(lmin:lmax,l)  &
                                                        + NoiseScale*EnoiseFac_11*rootE_22(l)*rootE_22(lmin:lmax)*Xi(PP_n_11_s_22)%EB(:,l)
                                                        if (PP_n_12_s_21/=0) C%C(lmin:lmax,l) =  C%C(lmin:lmax,l)  &
                                                        + NoiseScale*EnoiseFac_12*rootE_21(l)*rootE_21(lmin:lmax)*Xi(PP_n_12_s_21)%EB(:,l)
                                                        if (PP_s_11_n_22/=0) C%C(lmin:lmax,l) =  C%C(lmin:lmax,l)  &
                                                        + NoiseScale*EnoiseFac_22*rootE_11(l)*rootE_11(lmin:lmax)*Xi(PP_s_11_n_22)%EB(:,l)
                                                        if (PP_s_12_n_21/=0) C%C(lmin:lmax,l) =  C%C(lmin:lmax,l)  &
                                                        + NoiseScale*EnoiseFac_21*rootE_12(l)*rootE_12(lmin:lmax)*Xi(PP_s_12_n_21)%EB(:,l)
                                                    end if
                                                end do
                                            end if


                                        end do !pol
                                    end do


                                end do !y2
                            end do   !x2
                        end do
                    end do

                end do !y
            end do !x
        end do
    end do

    if (step==nstep) then
        ix=0
        i1=0
        do chanx = 1, nchannels
            do x=1,nweights
                i1 =i1+1
                i2=0
                do chany = 1, nchannels
                    do y=1,nweights
                        i2 = i2 + 1
                        if (i2 > i1) cycle
                        ix=ix+1

                        ix2=0
                        j1=0
                        do chanx2 = 1, nchannels
                            do x2=1,nweights
                                j1=j1+1
                                j2 = 0
                                do chany2 = 1, nchannels
                                    do y2=1,nweights
                                        j2 = j2+1
                                        if (j2 > j1) cycle

                                        ix2=ix2+1

                                        if (vec_size>3 .and. ix2<=ix) then
                                            if (nchannels>1) call MpiStop('Get covariance need to be more careful with symmetries')
                                            !Put in E/B coupling inverses
                                            !  print *,'Doing E/B cov coupling inverses', ix,ix2
                                            EE =>  CArr%Pol(3,3)%Cov(ix,ix2)
                                            BE =>  CArr%Pol(4,3)%Cov(ix,ix2)
                                            EB =>  CArr%Pol(4,3)%Cov(ix2,ix)  !Is the transpose
                                            BB =>  CArr%Pol(4,4)%Cov(ix,ix2)

                                            allocate(tmp11(lmin:lmax,lmin:lmax))
                                            allocate(tmp12(lmin:lmax,lmin:lmax))
                                            allocate(tmp22(lmin:lmax,lmin:lmax))
                                            allocate(tmp21(lmin:lmax,lmin:lmax))

                                            call Matrix_Mult(Coupler(ix)%invEE,EE%C,tmp11)
                                            call Matrix_Mult(Coupler(ix)%invEB,BE%C,tmp)
                                            tmp11 = tmp11 + tmp


                                            call Matrix_Mult_NT(Coupler(ix)%invEE,EB%C,tmp12)
                                            call Matrix_Mult(Coupler(ix)%invEB,BB%C,tmp)
                                            tmp12 = tmp12 + tmp


                                            call Matrix_Mult(Coupler(ix)%invEB,EE%C,tmp21)
                                            call Matrix_Mult(Coupler(ix)%invEE,BE%C,tmp)
                                            tmp21 = tmp21 + tmp


                                            call Matrix_Mult_NT(Coupler(ix)%invEB,EB%C,tmp22)
                                            call Matrix_Mult(Coupler(ix)%invEE,BB%C,tmp)
                                            tmp22 = tmp22 + tmp

                                            !!!Right multiply
                                            call Matrix_Mult_NT(tmp11,Coupler(ix2)%invEE,EE%C)
                                            call Matrix_Mult_NT(tmp12,Coupler(ix2)%invEB,tmp)
                                            EE%C = EE%C + tmp


                                            call Matrix_Mult_NT(tmp21,Coupler(ix2)%invEE,BE%C)
                                            call Matrix_Mult_NT(tmp22,Coupler(ix2)%invEB,tmp)
                                            BE%C = BE%C + tmp


                                            call Matrix_Mult_NT(tmp11,Coupler(ix2)%invEB,EB%C)
                                            call Matrix_Mult_NT(tmp12,Coupler(ix2)%invEE,tmp)
                                            EB%C = transpose(EB%C + tmp)

                                            call Matrix_Mult_NT(tmp21,Coupler(ix2)%invEB,BB%C)
                                            call Matrix_Mult_NT(tmp22,Coupler(ix2)%invEE,tmp)
                                            BB%C = BB%C + tmp


                                            deallocate(tmp11,tmp12,tmp22,tmp21)
                                        end if

                                        do polx=1,vec_size
                                            do poly=1,polx
                                                if (ix2>ix .and. polx==poly .or. polx==4 .and. poly<3) cycle
                                                C =>  CArr%Pol(polx,poly)%Cov(ix,ix2)
                                                do l=lmin,lmax
                                                    C%C(:,l)= C%C(:,l)*sqrt(DUmmyNoise(chanx2)%Cl(l,C_T)*DummyNoise(chany2)%Cl(l,C_T))
                                                    C%C(l,:)= C%C(l,:)*sqrt(DUmmyNoise(chanx)%Cl(l,C_T)*DummyNoise(chany)%Cl(l,C_T))
                                                end do

                                                if (polx==1 .and. poly==1 .and. frac_pt_src_error/=0 .and. incsignal) then
                                                    !Point source uncertainty
                                                    C%C = C%C + frac_pt_src_error**2 * &
                                                    (Channels(chanx)%PtSrcA*Channels(chany)%PtSrcA*Channels(chanx2)%PtSrcA*Channels(chany2)%PtSrcA)**0.5
                                                end if
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end if !step==nstep

    if (step==0) then
        print *,'Count check', count(chk), nvarmaps

        do j=1,size(Xi)
            print *,j,chk(j)
        end do
    end if
    deallocate(chk)
    if (step==0) return

    deallocate(tmp)

    print *,'Done main loops'

    do channel = 1, nchannels
        call HealpixPower_Free(DummyNoise(channel))
    end do
    deallocate(DummyNoise)
    deallocate(RootT)
    deallocate(rootT_11,rootT_12,rootT_21,rootT_22)
    if (dopol) then
        deallocate(RootE,RootB)
        deallocate(rootE_11,rootE_12,rootE_21,rootE_22)
        deallocate(rootB_11,rootB_12,rootB_21,rootB_22)
        deallocate(BeamAvTE)
    end if
    call HealpixPower_Free(AvP)
    call TBeam_Free(AvBeam)
    print *,'Time for Cov:', step, GeteTime()-STime

    end subroutine PseudoCl_GetFullCovariance


    subroutine PseudoCl_GetCHatNoise(P, Coupler, WeightMaps, WeightMapsPol,nweights, NoiseMap, noise_colour_factor)
    Type(HealpixPower) :: P(:)
    Type(TCouplingMatrix) :: Coupler(:)
    integer, intent(in) :: nweights
    Type(HealpixMap) :: NoiseMap, WeightMaps(nweights),WeightMapsPol(nweights)
    real, intent(in) :: noise_colour_factor
    double precision :: noise_fac !Need to be careful summing lots of single precisions
    double precision :: fac
    double precision, allocatable :: noise_colour(:)

    integer i, ix, x,y

    ix=0
    do x=1,nweights
        do y=1,x
            ix=ix+1

            call HealpixPower_Init(P(ix), Coupler(ix)%lmax, pol = Coupler(ix)%has_pol,dolens = .false.,nofree=.true.)

            fac = HO_fourpi/dble(NoiseMap%npix)**2

            print *,'mean w_i w_j', x,y, sum(dble(WeightMaps(x)%TQU(:,1)*WeightMaps(y)%TQU(:,1)))/NoiseMap%npix

            noise_fac = fac*sum(dble(NoiseMap%TQU(:,1))*WeightMaps(x)%TQU(:,1)*WeightMaps(y)%TQU(:,1))
            allocate(noise_colour(Coupler(ix)%lmin:Coupler(ix)%lmax))
            do i=Coupler(ix)%lmin, Coupler(ix)%lmax
                noise_colour(i) = 1+ noise_colour_factor/real(2*i+1)
            end do
            do i=Coupler(ix)%lmin, Coupler(ix)%lmax
                P(ix)%Cl(i,C_T) = sum(Coupler(ix)%InvT(i,:)*noise_colour)*noise_fac
            end do

            print *,'noise l=50  muK^2:', P(ix)%Cl(50,C_T) *mK**2
            print *,'noise l=500  muK^2:', P(ix)%Cl(500,C_T) *mK**2

            if (Coupler(ix)%has_pol) then
                noise_fac = fac/2*sum(dble((NoiseMap%TQU(:,2)+NoiseMap%TQU(:,3))*WeightMapsPol(x)%TQU(:,1)*WeightMapsPol(y)%TQU(:,1)))
                do i=Coupler(ix)%lmin, Coupler(ix)%lmax
                    P(ix)%Cl(i,C_E) = sum((Coupler(ix)%InvEE(i,:)+Coupler(ix)%InvEB(i,:))*noise_colour)*noise_fac
                    P(ix)%Cl(i,C_B) = P(ix)%Cl(i,C_E)
                end do
            end if
            deallocate(noise_colour)
        end do
    end do

    end subroutine PseudoCl_GetCHatNoise

    subroutine PseudoCl_WeightsToCovPowers(H, WeightMaps1, WeightMaps2,Channels, CovPowers, nweights, lmax, pol_weights)
    !Noise power maps factor of (HO_fourpi/dble(NoiseMap%npix))**2 large
    Type(HealpixInfo)  :: H
    integer, intent (in) :: nweights, lmax
    Type(HealpixMap) :: WeightMaps1(:), WeightMaps2(:)
    Type(TChannel) ::  Channels(:)
    logical, intent(in) :: pol_weights

    Type(HealpixCrossPowers) :: CovPowers
    Type(HealPixMap), dimension(:), target, allocatable :: CovMaps
    Type(HealPixMap), pointer :: amap
    integer nnoise, ncross
    integer ix,i,j, channel
    integer ncovmaps
    integer nmaskcomb, maskcomb
    double precision IniTime

    IniTime = GeteTime()

    nnoise = size(channels)
    ncross = nweights*(nweights+1)/2
    if (pol_weights) then
        nmaskcomb = 3
        ncovmaps = (3+2*nnoise)*ncross
    else
        nmaskcomb=1
        ncovmaps = (1+nnoise)*ncross
    end if

    allocate(CovMaps(ncovmaps))
    do i=1,ncovmaps
        call HealpixMap_Nullify(CovMaps(i))
    end do

    print *,'getting mask/noise powers, ncovmaps =', ncovmaps

    NoiseScale = (HO_fourpi/dble(Channels(1)%NoiseMap%npix))
    !Power of w_i*w_j, w_i*w_j sigma^2
    do maskcomb=1,nmaskcomb
        ix= (maskcomb-1)*(1+nnoise)*ncross
        do i=1, nweights
            do j=1,i
                ix= ix+1
                amap => CovMaps(ix)

                if (maskcomb==1) then
                    call HealpixMap_Assign(AMap,WeightMaps1(i))
                    amap%TQU = amap%TQU * WeightMaps1(j)%TQU
                else if (maskcomb==2) then
                    call HealpixMap_Assign(AMap,WeightMaps2(i))
                    amap%TQU = amap%TQU * WeightMaps2(j)%TQU
                else
                    call HealpixMap_Assign(AMap,WeightMaps1(i))
                    amap%TQU = amap%TQU * WeightMaps2(j)%TQU
                    !mixed P-T never have noise
                end if

                if (maskcomb/=3) then
                    do channel = 1, nnoise
                        call HealpixMap_Assign(CovMaps(ix+ncross*channel),Amap)
                        CovMaps(ix+ncross*channel)%TQU(:,1) = CovMaps(ix+ncross*channel)%TQU(:,1)  * Channels(channel)%NoiseMap%TQU(:,1)
                        ! *(HO_fourpi/dble(NoiseMap%npix))

                        !Avoid very small numbers
                        print *, 'covmap min/max '//trim(Channels(Channel)%Name), &
                        minval( CovMaps(ix+ncross*channel)%TQU(:,1)), maxval( CovMaps(ix+ncross*channel)%TQU(:,1))
                    end do
                end if
            end do
        end do
    end do

    call HealpixMapSet2CrossPowers(H, CovMaps, CovPowers, ncovmaps, lmax, .true.)
    deallocate(CovMaps)

    print *,'PseudoCl_WeightsToCovPowers time:', GeteTime() -  IniTime

    end subroutine PseudoCl_WeightsToCovPowers



    subroutine PseudoCl_InverseCovmatArr(CArr,minl_hybrid,lmax, CArrOut)
    !inverses pol-diagonal covmat
    use MatrixUtils
    Type(TCovMatArray), target :: CArr
    Type(TCovMatArray), optional, target :: CArrOut
    Type(TCovMatArray), pointer :: AC
    integer i,j, lmax
    real(dp), dimension(:,:), allocatable :: A
    integer nl,minl_hybrid, polx, poly

    if (present(CArrOut)) then
        AC => CArrOut
    else
        AC => CArr
    end if

    nl = lmax - minl_hybrid + 1

    print *,'Inverse covmat array nl, ncl = ',nl, CArr%ncl

    allocate(A(nl*CArr%ncl,nl*CArr%ncl))

    do i=1,CArr%ncl
        do j =1,i
            A((i-1)*nl+1:i*nl,(j-1)*nl+1:j*nl) = CArr%Cov(i,j)%C(minl_hybrid:lmax,minl_hybrid:lmax)
            if (i/=j) A((j-1)*nl+1:j*nl,(i-1)*nl+1:i*nl) = transpose(CArr%Cov(i,j)%C(minl_hybrid:lmax,minl_hybrid:lmax))
            if (.not. present(CArrOut)) deallocate(CArr%Cov(i,j)%C)
        end do
    end do

    !Regularize inverse; result only used for coupling so can do what we like without bias
    !Better would be to add diagonal elements proportional to the noise at low l where noise small
    do i=1,nl*CArr%ncl
        A(i,i) = A(i,i)*1.01
    end do

    call Matrix_Inverse(A)

    if (present(CArrOut)) then
        CArrOut%ncl = CArr%ncl
        allocate(CArrOut%Cov(CArr%ncl,CArr%ncl))
    end if

    do i=1,CArr%ncl
        do j =1,i
            allocate(AC%Cov(i,j)%C(minl_hybrid:lmax,minl_hybrid:lmax))
            AC%Cov(i,j)%C = A((i-1)*nl+1:i*nl,(j-1)*nl+1:j*nl)
            if (i/=j) nullify(AC%Cov(j,i)%C)
        end do
    end do

    deallocate(A)

    end subroutine PseudoCl_InverseCovmatArr

    subroutine PseudoCl_CovmatArrWriteold(CArr, fname)
    use MatrixUtils
    Type(TCovMatArray), target :: CArr
    character(LEN=*), intent(in) :: fname
    integer i,j
    integer nl, ncl, sz, file_unit

    file_unit = new_file_unit()
    call CreateFile(fname, file_unit,'unformatted')
    sz= size(CArr%Cov(1,1)%C,DIM=1)
    ncl = CArr%ncl
    print *,'sz,ncl = ',sz,ncl
    write (file_unit) ncl, sz
    do i=1,CArr%ncl
        do j =1,i
            write(file_unit) CArr%Cov(i,j)%C
        end do
    end do
    call CloseFile(file_unit)

    end subroutine PseudoCl_CovmatArrWriteold

    subroutine PseudoCl_CovmatArrReadold(CArr, fname)
    use MatrixUtils
    Type(TCovMatArray), target :: CArr
    character(LEN=*), intent(in) :: fname
    integer i,j
    integer nl, ncl, sz, file_unit

    file_unit = new_file_unit()
    call OpenFile(fname, file_unit,'unformatted')

    read  (file_unit) ncl, sz
    CArr%ncl = ncl
    allocate(Carr%Cov(ncl,ncl))
    do i=1,ncl
        do j =1,i
            allocate(Carr%Cov(i,j)%C(sz,sz))
            read(file_unit) CArr%Cov(i,j)%C
            if (i/=j)  nullify(CArr%Cov(j,i)%C)
        end do
    end do
    call CloseFile(file_unit)

    end subroutine PseudoCl_CovmatArrReadold

    subroutine PseudoCl_CovmatArr1DWrite(CArr, fname)
    use MatrixUtils
    Type(TCovMatSet), target :: CArr
    character(LEN=*), intent(in) :: fname
    integer i,j
    integer nl, ncl, sz, file_unit

    file_unit = new_file_unit()
    call CreateFile(fname, file_unit,'unformatted')
    ncl = CArr%n
    write (file_unit) ncl, CArr%lmin, CArr%lmax
    do i=1,CArr%n
        write(file_unit) CArr%Cov(i)%C
    end do
    call CloseFile(file_unit)

    end subroutine PseudoCl_CovmatArr1DWrite

    subroutine PseudoCl_CovmatArr1DRead(CArr, fname)
    use MatrixUtils
    Type(TCovMatSet), target :: CArr
    character(LEN=*), intent(in) :: fname
    integer i
    integer nl, ncl, file_unit, lmin, lmax

    file_unit = new_file_unit()
    call OpenFile(fname, file_unit,'unformatted')

    read  (file_unit) ncl, lmin, lmax
    CArr%n = ncl
    allocate(Carr%Cov(ncl))
    do i=1,ncl
        allocate(Carr%Cov(i)%C(lmin:lmax,lmin:lmax))
        read(file_unit) CArr%Cov(i)%C
    end do
    call CloseFile(file_unit)

    end subroutine PseudoCl_CovmatArr1DRead

    subroutine PseudoCl_FreePowerArray2(P)
    Type(HealpixPower) :: P(:,:)
    integer i,j

    do i=1,size(P,1)
        do j=1,size(P,2)
            print *,i,j
            call HealpixPower_Free(P(i,j))
        end do
    end do
    end subroutine PseudoCl_FreePowerArray2

    end module PseudoCl
