    !Module for computing coupling matrices over non-symmetric sky cuts

    !Call GetCouplingAsymm(cut_Alm_file, lmax, s, plusonly, almaxr)
    !Fills in the WAsymm array (if s=0), or WPAsymm and WMAsymm arrays (s/=0).
    !WAsymm = W_{(lm)(lm)'}, WPAsymm = W_{+,(lm)(lm)'}, WMAsymm = W_{-,(lm)(lm)'}
    !The mapping from (lm) to an index is given by the idx(s,l,m) function
    !If plusonly = T, only computes WPAsymm, if almaxr is present, uses up to almaxr on right matrix indices
    !The cut it taken from a text file of alms, cut_Alm_file, e.g. as converted from HEALPIX output

    module CutSkyAsymm
    use AMLUtils
    use HealpixObj
    use PseudoCl
    use MatrixUtils
    implicit none

    Type AsymmCouplings
        integer lmaxl, lmaxr
        real(dp), pointer :: WAsymm(:,:), WPAsymm(:,:), WMAsymm(:,:)
    end  Type AsymmCouplings

    Type VecArray
        real(DP), dimension(:), pointer ::  V
    end Type VecArray

    Type TComplexCovMat
        complex(dp), dimension(:,:), pointer :: C
    end Type TComplexCovMat


    Type ProjMat
        integer nl,nr
        real(dp), pointer :: RootDiag(:) => null()
        real(dp), pointer :: M(:,:) => null()
        real(dp), pointer :: Minv(:,:) => null()
    end  Type ProjMat

    Type ProjMatPol
        integer nr,nl
        real(dp), pointer :: RootDiag(:)
        complex(dp), pointer :: WComp(:,:) !complex one
        real(dp), pointer :: EProj(:,:), BProj(:,:)
    end  Type ProjMatPol


    logical :: asymm_pol = .false.
    logical :: EBOnly = .false.
    real(dp), parameter :: NoiseProject = 1e7_dp

    contains

    function idx(s,l,m)
    !Conversion from l,m into vectors
    integer, intent(in) :: l,m,s
    integer idx

    idx = l*(l+1) + m +1 - s**2

    end function


    function idx_T(l,m)
    !Conversion from l,m into vectors
    integer, intent(in) :: l,m
    integer idx_T

    idx_T = l*(l+1) + m +1

    end function

    function idx_P(l,m)
    !Conversion from l,m into vectors
    integer, intent(in) :: l,m
    integer idx_P

    idx_P = l*(l+1) + m +1 -4

    end function


    subroutine HealpixAlm2Vector(Alm, v, lmax, polix)
    Type (HealpixAlm) :: Alm
    integer, intent(in) :: lmax, polix
    real(sp) :: v(:)
    integer off,l,m

    if (polix==1) then
        do l=0, lmax
            v(idx_T(l,0))=real(Alm%TEB(1,l,0))
        end do

        do m=1, lmax
            do l=m, lmax
                v(idx_T(l,m))=root2*real(Alm%TEB(1,l,m))
                v(idx_T(l,-m))=root2*aimag(Alm%TEB(1,l,m))
            end do
        end do
    else
        if (Alm%npol<3) stop 'HealpixAlm2Vector: npol < 3'
        do l=2, lmax
            v(idx_P(l,0))=real(Alm%TEB(polix,l,0))
        end do

        do m=1, lmax
            do l=max(2,m), lmax
                v(idx_P(l,m))=root2*real(Alm%TEB(polix,l,m))
                v(idx_P(l,-m))=root2*aimag(Alm%TEB(polix,l,m))
            end do
        end do
    end if

    end subroutine HealpixAlm2Vector


    subroutine HealpixVector2Alm(v,Alm, lmax, polix)
    Type (HealpixAlm) :: Alm
    integer, intent(in) :: lmax
    integer, intent(in) :: polix
    real(sp) :: v(:)
    integer l,m

    if (polix==1) then
        do l=0, lmax
            Alm%TEB(1,l,0) = v(idx_T(l,0))
        end do

        do m=1, lmax
            do l=m, lmax
                Alm%TEB(1,l,m) = cmplx(v(idx_T(l,m)),v(idx_T(l,-m)))/root2
            end do
        end do
    else
        if (Alm%npol<3) stop 'HealpixVector2Alm: npol < 3'
        do l=2, lmax
            Alm%TEB(polix,l,0) = v(idx_P(l,0))
        end do

        do m=1, lmax
            do l=max(2,m), lmax
                Alm%TEB(polix,l,m) = cmplx(v(idx_P(l,m)),v(idx_P(l,-m)))/root2
            end do
        end do
    end if

    end subroutine HealpixVector2Alm

    subroutine HealpixVectorMultPower(AlmVec,P, lmax)
    real(sp) :: AlmVec(:)
    Type(HealpixPower) :: P
    integer, intent(in) :: lmax
    integer l,m

    do l=0, lmax
        print *,'l=',l
        do m=-l,l
            AlmVec(idx_T(l,m)) = AlmVec(idx_T(l,m)) *P%Cl(l,C_T)
        end do
    end do

    end  subroutine HealpixVectorMultPower



    subroutine HealpixAlm2BigVec(Alm, v, lmax, dopol)
    Type (HealpixAlm) :: Alm
    logical, intent(in) :: dopol
    integer, intent(in) :: lmax
    real(sp) :: v(:)
    integer off,l,m

    do l=0, lmax
        v(idx(0,l,0))=real(Alm%TEB(1,l,0))
    end do

    do m=1, lmax
        do l=m, lmax
            v(idx(0,l,m))=root2*real(Alm%TEB(1,l,m))
            v(idx(0,l,-m))=root2*aimag(Alm%TEB(1,l,m))
        end do
    end do

    if (dopol) then
        off = (lmax+1)**2
        do l=2, lmax
            v(off+idx_P(l,0))=real(Alm%TEB(2,l,0))
        end do

        do m=1, lmax
            do l=max(2,m), lmax
                v(off+idx_P(l,m))=root2*real(Alm%TEB(2,l,m))
                v(off+idx_P(l,-m))=root2*aimag(Alm%TEB(2,l,m))
            end do
        end do

        off = (lmax+1)**2*2-4

        do l=2, lmax
            v(off+idx_P(l,0))=real(Alm%TEB(3,l,0))
        end do
        do m=1, lmax
            do l=max(2,m), lmax
                v(off+idx_P(l,m))=root2*real(Alm%TEB(3,l,m))
                v(off+idx_P(l,-m))=root2*aimag(Alm%TEB(3,l,m))
            end do
        end do
    end if

    end subroutine HealpixAlm2BigVec

    subroutine CutSkyAsymm_GetCoupling(M,WinA, WinAP,lmax, plusonly, almaxr, WantTemp, WantPol)
    Type (HealpixAlm), intent(in) :: WinA, WinAP
    Type (AsymmCouplings) :: M
    logical, intent(in), optional :: WantTemp, WantPol
    integer, intent(in) :: lmax
    integer, intent(in), optional :: almaxr
    logical, intent(in), optional :: plusonly
    integer i,j, l1,l2,m1,m2,m3,lmin,lmx,offs
    real(dp) threejs(lmax*2+1),theejs_minusm(lmax*2+1),threejs0(lmax*2+1),threejsm(lmax*2+1)
    real(dp) threejs2(lmax*2+1),threejsm2(lmax*2+1)
    real(dp) start,fac
    complex(dp), dimension(:,:), allocatable :: Wlm, WlmP
    complex(dp) Wpp,Wpm,tmp2,tmp
    complex(dp) Wplus_pp,Wplus_pm,iWminus_pp,iWminus_pm
    logical dominus
    integer lmaxl,lmaxr,lmaxcut
    logical doTemp, doPol

    if (present(WantTemp)) then
        DoTemp=WantTemp
    else
        DoTemp = .true.
    end if
    if (present(WantPol)) then
        doPol = WantPol
    else
        doPol = asymm_pol
    end if

    lmaxl = lmax
    if (present(almaxr)) then
        lmaxr=almaxr
    else
        lmaxr=lmax
    end if

    dominus = .true.
    if (present(plusonly)) dominus = .not. plusonly

    !Have to go to 2*lmax
    if (doTemp) then
        lmaxcut = min(WinA%lmax,lmaxr*2)
        allocate(Wlm(0:lmaxcut,0:lmaxcut))
    end if
    if (doPol) then
        lmaxcut = min(WinAP%lmax,lmaxr*2)
        allocate(WlmP(0:lmaxcut,0:lmaxcut))
    end if
    if (lmaxcut < lmaxr*2) write(*,*) 'WARNING: getting coupling with lmaxcut < lmax*2'
    do l1=0,lmaxcut
        if (doTemp) Wlm(l1,0:l1) = WinA%TEB(1,l1,0:l1)*sqrt((2*l1+1)/fourpi)
        if (doPol) WlmP(l1,0:l1) = WinAP%TEB(1,l1,0:l1)*sqrt((2*l1+1)/fourpi)
    end do

    M%lmaxl = lmaxl
    M%lmaxr = lmaxr

    if (lmaxl/=lmaxr) call MpiStop('currently assume lmaxl =lmaxr')

    if (doTemp) then
        allocate(M%WAsymm((lmaxl+1)**2,(lmaxr+1)**2))
        M%WAsymm = 0
    end if

    if (doPol) then
        allocate(M%WPAsymm((lmaxl+1)**2-4,(lmaxr+1)**2-4))
        M%WPAsymm = 0
        if (dominus) then
            allocate(M%WMAsymm((lmaxl+1)**2-4,(lmaxr+1)**2-4))
            M%WMAsymm = 0
        end if
    end if

    print *,'Getting coupling matrices'

    !$OMP PARALLEL DO DEFAULT(PRIVATE), SCHEDULE(DYNAMIC), SHARED(lmaxl,lmaxcut,doPol,DoTemp,Wlm,M,WlmP,dominus)
    do l1=0, lmaxl
        do l2=0,l1
            call GetThreeJs(threejs0,l1,l2,0,0)

            if (doPol .and. l1>=2 .and. l2>=2) then
                call GetThreeJs(threejs2,l1,l2,2,-2)
                call GetThreeJs(threejsm2,l1,l2,-2,2)
            end if

            lmx = min(lmaxcut,l1+l2)

            do m1=0,l1
                fac = sqrt(real((2*l1+1)*(2*l2+1),dp))
                if (mod(m1,2)/=0) fac=-fac

                do m2 = 0,l2
                    ! (-1)^m1 factor put in fac above

                    !m1 m2
                    m3 = m1-m2
                    lmin = max(abs(l1-l2),abs(m3))
                    offs = lmin - abs(l1-l2)

                    call GetThreeJs(threejs,l1,l2,-m1,m2)
                    if (DoTemp) then
                        Wpp = fac*sum(Wlm(lmin:lmx,abs(m3))*threejs(1:lmx-lmin+1)*threejs0(1+offs:lmx-lmin+1+offs))
                        if(m3<0) Wpp = conjg(Wpp)*(-1)**m3
                    end if
                    !m1 -m2
                    if (m2/=0) then
                        m3 = m1+m2
                        lmin = max(abs(l1-l2),abs(m3))
                        offs = lmin - abs(l1-l2)
                        call GetThreeJs(theejs_minusm,l1,l2,-m1,-m2)
                        if (DoTemp) then
                            Wpm = fac*sum(Wlm(lmin:lmx,abs(m3))*theejs_minusm(1:lmx-lmin+1)*threejs0(1+offs:lmx-lmin+1+offs))
                            if(m3<0) Wpm = conjg(Wpm)*(-1)**m3
                            Wpm  = (-1)**m2*Wpm !only enters with this
                        end if
                    else
                        Wpm=0
                    end if

                    !we are calculating W(l1,m1, l2,m2)
                    if (DoTemp) then
                        if (m2==0) then
                            if (m1==0) then
                                M%WAsymm(idx_T(l1,0),idx_T(l2,0)) =  real(Wpp)
                            else
                                M%WAsymm(idx_T(l1,m1),idx_T(l2,0)) =  root2*real(Wpp)
                                M%WAsymm(idx_T(l1,-m1),idx_T(l2,0)) =  root2*aimag(Wpp)
                            end if
                        else
                            if (m1==0) then
                                M%WAsymm(idx_T(l1,0),idx_T(l2,m2)) =  root2*real(Wpp)
                                M%WAsymm(idx_T(l1,0),idx_T(l2,-m2)) =  -root2*aimag(Wpp)
                            else
                                M%WAsymm(idx_T(l1,m1),idx_T(l2,m2)) =  real(Wpp) + real(Wpm)
                                M%WAsymm(idx_T(l1,-m1),idx_T(l2,m2)) =  aimag(Wpp) + aimag(Wpm)

                                M%WAsymm(idx_T(l1,m1),idx_T(l2,-m2)) =  -aimag(Wpp) + aimag(Wpm)
                                M%WAsymm(idx_T(l1,-m1),idx_T(l2,-m2)) =  real(Wpp) - real(Wpm)
                            end if
                        end if
                    end if

                    if (.not. doPol .or. l1<2 .or. l2<2) cycle

                    m3 = m1-m2
                    lmin = max(abs(l1-l2),abs(m3))
                    offs = lmin - abs(l1-l2)
                    tmp =sum(WlmP(lmin:lmx,abs(m3))*threejs(1:lmx-lmin+1)*threejs2(1+offs:lmx-lmin+1+offs))
                    if(m3<0) tmp = conjg(tmp)*(-1)**m3
                    tmp = tmp*fac

                    tmp2 =sum(WlmP(lmin:lmx,abs(m3))*threejs(1:lmx-lmin+1)*threejsm2(1+offs:lmx-lmin+1+offs))
                    if(m3<0) tmp2 = conjg(tmp2)*(-1)**m3
                    tmp2 = tmp2*fac

                    Wplus_pp=(tmp+tmp2)/2
                    iWminus_pp=cmplx(0,1)*(tmp-tmp2)/2

                    !m1 -m2
                    if (m2/=0) then
                        m3 = m1+m2
                        lmin = max(abs(l1-l2),abs(m3))
                        offs = lmin - abs(l1-l2)
                        tmp = fac*sum(WlmP(lmin:lmx,abs(m3))*theejs_minusm(1:lmx-lmin+1)*threejs2(1+offs:lmx-lmin+1+offs))
                        if(m3<0) tmp = conjg(tmp)*(-1)**m3
                        tmp  = (-1)**m2*tmp !only enters with this

                        tmp2 = fac*sum(WlmP(lmin:lmx,abs(m3))*theejs_minusm(1:lmx-lmin+1)*threejsm2(1+offs:lmx-lmin+1+offs))
                        if(m3<0) tmp2 = conjg(tmp2)*(-1)**m3
                        tmp2  = (-1)**m2*tmp2 !only enters with this

                        Wplus_pm=(tmp+tmp2)/2
                        iWminus_pm=cmplx(0,1)*(tmp-tmp2)/2
                    else
                        Wplus_pm=0
                        iWminus_pm=0
                    end if

                    if (m2==0) then
                        if (m1==0) then
                            M%WPAsymm(idx_P(l1,0),idx_P(l2,0)) =  real(Wplus_pp)
                            if (dominus) then
                                !zero..
                                M%WMAsymm(idx_P(l1,0),idx_P(l2,0)) =  0
                                if (real(iWminus_pp) /=0) print *,'test not zero WM m1=0,m2=0'
                            end if
                        else
                            M%WPAsymm(idx_P(l1,m1),idx_P(l2,0)) =  root2*real(Wplus_pp)
                            M%WPAsymm(idx_P(l1,-m1),idx_P(l2,0)) =  root2*aimag(Wplus_pp)
                            if (dominus) then
                                M%WMAsymm(idx_P(l1,m1),idx_P(l2,0)) =  root2*real(iWminus_pp)
                                M%WMAsymm(idx_P(l1,-m1),idx_P(l2,0)) =  root2*aimag(iWminus_pp)
                            end if
                        end if
                    else
                        if (m1==0) then
                            M%WPAsymm(idx_P(l1,0),idx_P(l2,m2)) =  root2*real(Wplus_pp)
                            M%WPAsymm(idx_P(l1,0),idx_P(l2,-m2)) =  -root2*aimag(Wplus_pp)
                            if (dominus) then
                                M%WMAsymm(idx_P(l1,0),idx_P(l2,m2)) =  root2*real(iWminus_pp)
                                M%WMAsymm(idx_P(l1,0),idx_P(l2,-m2)) =  -root2*aimag(iWminus_pp)
                            end if
                        else
                            M%WPAsymm(idx_P(l1,m1),idx_P(l2,m2)) =  real(Wplus_pp) + real(Wplus_pm)
                            M%WPAsymm(idx_P(l1,-m1),idx_P(l2,m2)) =  aimag(Wplus_pp) + aimag(Wplus_pm)

                            M%WPAsymm(idx_P(l1,m1),idx_P(l2,-m2)) =  -aimag(Wplus_pp) + aimag(Wplus_pm)
                            M%WPAsymm(idx_P(l1,-m1),idx_P(l2,-m2)) =  real(Wplus_pp) - real(Wplus_pm)

                            if (dominus) then
                                M%WMAsymm(idx_P(l1,m1),idx_P(l2,m2)) =  real(iWminus_pp) + real(iWminus_pm)
                                M%WMAsymm(idx_P(l1,-m1),idx_P(l2,m2)) =  aimag(iWminus_pp) + aimag(iWminus_pm)

                                M%WMAsymm(idx_P(l1,m1),idx_P(l2,-m2)) =  -aimag(iWminus_pp) + aimag(iWminus_pm)
                                M%WMAsymm(idx_P(l1,-m1),idx_P(l2,-m2)) =  real(iWminus_pp) - real(iWminus_pm)
                            end if
                        end if
                    end if
                end do
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    print *,'Got coupling matrices'

    if (DoTemp) then
        do i=1, (lmaxr+1)**2
            do j=1,i-1
                M%WAsymm(j,i) =M%WAsymm(i,j)
            end do
        end do
    end if

    if (doPol) then
        do i=1, (lmaxr+1)**2-4
            do j=1,i-1
                M%WPAsymm(j,i) =M%WPAsymm(i,j)
                if (dominus) M%WMAsymm(j,i) =-M%WMAsymm(i,j)  !Antisymmetric
            end do
        end do
    end if

    if (doTemp) deallocate(Wlm)
    if (doPol) deallocate(WlmP)

    end subroutine CutSkyAsymm_GetCoupling


    subroutine CutSkyAsymm_GetUnmixedPolCovariance(Couple, Cov, P, lmin, lmax)
    !Just take C_E = C_B
    Type (TComplexCovMat) :: Cov
    Type(ProjMatPol) :: Couple
    Type(HealpixPower) :: P
    integer, intent(in) :: lmin, lmax
    integer i,j
    complex(dp) :: tmp
    integer l, mix
    !Couple is the transpose so fast sums in dot_product

    allocate(Cov%C(Couple%nr,Couple%nr))

    do i=1, Couple%nr
        do j=1, i
            tmp = 0
            do l=lmin,lmax
                mix = idx_P(l,-l)
                tmp =tmp + dot_product(Couple%WComp(mix:mix+2*l,i),Couple%WComp(mix:mix+2*l,j))*P%Cl(l,C_E)
            end do
            if (associated(Couple%RootDiag)) tmp = tmp*Couple%RootDiag(i)*Couple%RootDiag(j)
            Cov%C(i,j) = tmp
            if (i/=j) Cov%C(j,i) =  conjg(tmp)
        end do
    end do

    end subroutine CutSkyAsymm_GetUnmixedPolCovariance


    subroutine CutSkyAsymm_GetCovariance(Couple, Cov, P, lmin, lmax, dipole_monopole_add)
    Type (TCovMat) :: Cov
    Type(ProjMat) :: Couple
    Type(HealpixPower) :: P
    logical, intent(in) :: dipole_monopole_add
    integer, intent(in) :: lmin, lmax
    integer i,j
    real(dp) :: tmp
    integer l, mix
    !Couple is the transpose so fast sums in dot_product

    allocate(Cov%C(Couple%nr,Couple%nr))

    do i=1, Couple%nr
        do j=1, i
            tmp = 0

            if (dipole_monopole_add) then
                do l=0,1
                    !Large noise in dipole/monopole
                    mix = idx(0,l,-l)
                    tmp =tmp + dot_product(Couple%M(mix:mix+2*l,i),Couple%M(mix:mix+2*l,j))*NoiseProject
                end do
            end if

            do l=lmin,lmax
                mix = idx(0,l,-l)
                tmp =tmp + dot_product(Couple%M(mix:mix+2*l,i),Couple%M(mix:mix+2*l,j))*P%Cl(l,C_T)
            end do
            if (associated(Couple%RootDiag)) tmp = tmp*Couple%RootDiag(i)*Couple%RootDiag(j)
            Cov%C(i,j) = tmp
            if (i/=j) Cov%C(j,i) =  tmp
        end do
    end do

    end subroutine CutSkyAsymm_GetCovariance

    subroutine CutSkyAsymm_GetFullCovariance(Couple, PolCouple,  Cov, P, lmin, lmax)
    Type (TCovMat) :: Cov
    Type(ProjMat) :: Couple
    Type(ProjMatPol) :: PolCouple
    Type(HealpixPower) :: P
    integer, intent(in) :: lmin, lmax
    integer i,j
    real(dp) :: sum1,sum2,tmp,tmpEE,tmpBB, tmpEB, tmpBE
    integer l, mix, mixT, nmodes
    !Couple is the transpose so fast sums in dot_product
    real(dm) ddot
    external ddot

    nmodes = Couple%nr
    if (asymm_pol) then
        nmodes = nmodes + 2*PolCouple%nr
    end if
    allocate(Cov%C(nmodes,nmodes))
    Cov%C=0

    do i=1, Couple%nr
        do j=1, i
            tmp = 0

            do l=lmin,lmax
                mix = idx_T(l,-l)
                tmp =tmp + dot_product(Couple%M(mix:mix+2*l,i),Couple%M(mix:mix+2*l,j))*P%Cl(l,C_T)
                !     tmp = tmp + ddot(2*l+1, Couple%M(mix,i), 1, Couple%M(mix,j), 1)*P%Cl(l,C_T)
            end do
            Cov%C(i,j) = tmp
            if (i/=j) Cov%C(j,i) =  tmp
        end do
    end do

    if (asymm_pol) then
        do i=1, Couple%nr
            do j=1, PolCouple%nr
                !TE
                tmp = 0
                do l=lmin,lmax
                    mix = idx(2,l,-l)
                    mixT = idx(0,l,-l)
                    tmp =tmp + dot_product(Couple%M(mixT:mixT+2*l,i),&
                    PolCouple%EProj(mix:mix+2*l,j))*P%Cl(l,C_C)
                end do
                Cov%C(i,Couple%nr+j) = tmp
                Cov%C(Couple%nr+j,i) =  tmp

                !TB
                tmp = 0
                do l=lmin,lmax
                    mix = idx(2,l,-l)
                    mixT = idx_T(l,-l)
                    !            tmp =tmp + dot_product(Couple%M(mixT:mixT+2*l,i),&
                    !                -PolCouple%BProj(mix:mix+2*l,j))*P%Cl(l,C_C)
                    tmp = tmp - ddot(2*l+1, Couple%M(mixT,i), 1, PolCouple%BProj(mix,j), 1)*P%Cl(l,C_C)
                end do
                Cov%C(i,Couple%nr+PolCouple%nr+j) = tmp
                Cov%C(Couple%nr+PolCouple%nr+j,i) =  tmp
            end do
        end do

        do i=1, PolCouple%nr
            do j=1, i
                !EE
                tmpEE = 0
                tmpBB = 0
                do l=lmin,lmax
                    mix = idx_P(l,-l)
                    sum1 = ddot(2*l+1, PolCouple%EProj(mix,i), 1, PolCouple%EProj(mix,j), 1)
                    sum2 = ddot(2*l+1, PolCouple%BProj(mix,i), 1, PolCouple%BProj(mix,j), 1)
                    !            sum1=dot_product(PolCouple%EProj(mix:mix+2*l,i),PolCouple%EProj(mix:mix+2*l,j))
                    !           sum2=dot_product(PolCouple%BProj(mix:mix+2*l,i),PolCouple%BProj(mix:mix+2*l,j))
                    tmpEE =tmpEE + sum1*P%Cl(l,C_E)  + sum2*P%Cl(l,C_B)
                    tmpBB =tmpBB + sum1*P%Cl(l,C_B)  + sum2*P%Cl(l,C_E)
                end do
                !if (i==j) print *, i, tmpEE, tmpBB
                ! if (i/=j .and. abs(tmpEE) > 1e-5) print *,i,j, tmpEE, tmpBB
                Cov%C(Couple%nr+i,Couple%nr+j) = tmpEE
                if (i/=j) Cov%C(Couple%nr+j,Couple%nr+i) =  tmpEE
                Cov%C(Couple%nr+PolCouple%nr+i,Couple%nr+PolCouple%nr+j) = tmpBB
                if (i/=j) Cov%C(Couple%nr+PolCouple%nr+j,Couple%nr+PolCouple%nr+i) =  tmpBB

                !EB/BE
                tmpEB = 0
                tmpBE=0
                do l=lmin,lmax
                    mix = idx_P(l,-l)
                    sum1 = ddot(2*l+1, PolCouple%EProj(mix,i), 1, PolCouple%BProj(mix,j), 1)
                    sum2 = ddot(2*l+1, PolCouple%BProj(mix,i), 1, PolCouple%EProj(mix,j), 1)
                    !          sum1=dot_product(PolCouple%EProj(mix:mix+2*l,i),PolCouple%BProj(mix:mix+2*l,j))
                    !          sum2=dot_product(PolCouple%BProj(mix:mix+2*l,i),PolCouple%EProj(mix:mix+2*l,j))
                    tmpEB =tmpEB - sum1*P%Cl(l,C_E) + sum2*P%Cl(l,C_B)
                    tmpBE =tmpBE + sum1*P%Cl(l,C_B) - sum2*P%Cl(l,C_E)
                end do
                Cov%C(Couple%nr+i,Couple%nr+PolCouple%nr+j) = tmpEB
                Cov%C(Couple%nr+PolCouple%nr+j,Couple%nr+i) =  tmpEB
                Cov%C(Couple%nr+PolCouple%nr+i,Couple%nr+j) = tmpBE
                Cov%C(Couple%nr+j,Couple%nr+PolCouple%nr+i) =  tmpBE
            end do
        end do
    end if

    end subroutine CutSkyAsymm_GetFullCovariance


    subroutine CutSkyAsymm_GetNoiseCovariance(Proj, ProjPol,Cov, CovPol,  NoiseA,NoiseAP, lmax, NoiseCov)
    integer, intent(in):: lmax
    Type(HealpixAlm) :: NoiseA, NoiseAP
    Type(ProjMat) :: Proj
    Type(ProjMatPol) :: ProjPol
    Type (TCovMat) :: Cov
    Type(AsymmCouplings) NoiseCov
    Type (TComplexCovmat) :: CovPol
    integer i,j,n
    complex(dp), allocatable :: CNoise(:,:)
    !Assume Q and U uncorrelated and equal variance for the moment

    call CutSkyAsymm_GetCoupling(NoiseCov, NoiseA, NoiseAP, lmax, plusonly = .false.)
    allocate(Cov%C(Proj%nr,Proj%nr))
    call Matrix_RotateSymm(NoiseCov%WAsymm, Proj%M, Proj%nr,  Cov%C)

    do i=1, Proj%nr
        do j=1, Proj%nr
            Cov%C(j,i) =  Cov%C(j,i)/(Proj%RootDiag(i)*Proj%RootDiag(j))
        end do
    end do

    !deallocate(NoiseCov%WAsymm)

    if (asymm_pol) then
        n = size(NoiseCov%WPAsymm,DIM=1)
        allocate(CNoise(n,n))
        CNoise= NoiseCov%WPAsymm - cmplx(0,1)*NoiseCov%WMAsymm
        !CNoise = 1/2*<P P^dag>
        deallocate(NoiseCov%WPAsymm)
        deallocate(NoiseCov%WMAsymm)
        allocate(CovPol%C(ProjPol%nr,ProjPol%nr))


        !Gets U^dag Mat U
        call Matrix_CRotateSymm(CNoise, ProjPol%WComp, ProjPol%nr,  CovPol%C)
        deallocate(CNoise)

        if (associated(ProjPol%RootDiag)) then
            do i=1,ProjPol%nr
                do j=1,ProjPol%nr
                    CovPol%C(j,i) =  CovPol%C(j,i)/(ProjPol%RootDiag(i)*ProjPol%RootDiag(j))
                end do
            end do
        end if

        !
        !        allocate(CovPol%C(ProjPol%nr,ProjPol%nr))
        !
        !        call Matrix_RotateSymm(NoiseCov%WPAsymm, ProjPol%WP_Proj%M, ProjPol%nr,  CovPol%C)
        !        deallocate(NoiseCov%WPAsymm)
        !
        !        if (associated(ProjPol%RootDiag)) then
        !        do i=1,ProjPol%nr
        !        do j=1,ProjPol%nr
        !            CovPol%C(j,i) =  CovPol%C(j,i)/(ProjPol%RootDiag(i)*ProjPol%RootDiag(j))
        !        end do
        !        end do
        !        end if
        !
        !        allocate(CovPolM%C(ProjPol%nr,ProjPol%nr))
        !        call Matrix_RotateAntiSymm(NoiseCov%WMAsymm, ProjPol%WP_Proj%M, ProjPol%nr,  CovPolM%C)
        !        deallocate(NoiseCov%WMAsymm)
        !
        !        if (associated(ProjPol%RootDiag)) then
        !        do i=1,ProjPol%nr
        !        do j=1,ProjPol%nr
        !            CovPolM%C(j,i) =  CovPolM%C(j,i)  /(ProjPol%RootDiag(i)*ProjPol%RootDiag(j))
        !        end do
        !        end do
        !        end if

        !        if (.not. associated(ProjPol%WP_dataProj%M,ProjPol%WP_Proj%M)) then
        !         deallocate(ProjPol%WP_dataProj%M)
        !         deallocate(ProjPol%WP_dataProj%RootDiag)
        !         nullify(ProjPol%WP_dataProj%M)
        !        end if
    end if

    end  subroutine CutSkyAsymm_GetNoiseCovariance


    subroutine CutSkyAsymm_GetSupportedModes(M, Proj, SUPPORT)
    real(dp) :: M(:,:)
    Type(ProjMat) :: Proj
    integer n, nfound
    real(dp), allocatable :: diag(:)
    real(dp), intent(in), optional :: SUPPORT
    real(dp) :: SUP
    real(dp), parameter :: MINSUPPORT = 1e-4_dp

    if (present(SUPPORT)) then
        SUP = SUPPORT
    else
        SUP = MINSUPPORT
    end if
    n=size(M,DIM=1)
    if (n/=size(M,DIM=2)) call MpiStop('CutSkyAsymm_GetSupportedModes: assumed square')

    allocate(diag(n))
    nfound = n
    call Matrix_Diagonalize_Partial(M, diag, n, SUP,1e100_dp, nfound)
    print *,'n modes = ', nfound, ' of ', n
    print *,'n modes > 0.9', count(diag(1:nfound)>0.9)
    print *,'n modes > 0.99', count(diag(1:nfound)>0.99)
    print *,'n modes > 0.999', count(diag(1:nfound)>0.999)
    print *, 'min/max diag = ', minval(diag(1:nfound)),maxval(diag(1:nfound))
    print *,'sum diag^2 = ', sum(diag(1:nfound)**2)

    allocate(Proj%M(n,nfound))
    Proj%nr = nfound
    Proj%nl = n
    Proj%M = M(:,1:nfound)
    allocate(proj%RootDiag(nfound))
    proj%RootDiag = sqrt(diag(1:nfound))
    deallocate(diag)

    end subroutine CutSkyAsymm_GetSupportedModes

    subroutine CutSkyAsymm_SupportedModesPolFromCoupling(Couplings, Proj, SUPPORT)
    Type (AsymmCouplings) :: Couplings
    Type(ProjMatPol) :: Proj
    real(dp), intent(in), optional :: SUPPORT
    complex(dp), allocatable :: Mat(:,:)
    integer n,m
    
    !Diagonalize (W_+ W_-, -W_- W_+) by diagonalizing hermitian W_+ + iW_i

    m = size(Couplings%WPAsymm,DIM=1)
    n = size(Couplings%WPAsymm,DIM=2)
    if (n/=m) call MpiStop('not checked m/=n')

    allocate(Mat(m,n))
    Mat = Couplings%WPAsymm - cmplx(0,1)*Couplings%WMAsymm  !W_+ + iW_-

    deallocate(Couplings%WMAsymm)
    deallocate(Couplings%WPAsymm)

    call CutSkyAsymm_GetSupportedModesPol(Mat, Proj, SUPPORT)
    deallocate(Mat)

    end subroutine CutSkyAsymm_SupportedModesPolFromCoupling

    subroutine CutSkyAsymm_GetSupportedModesPol(M, Proj, SUPPORT)
    complex(dp) :: M(:,:)
    Type(ProjMatPol) :: Proj
    integer n, nfound
    real(dp), allocatable :: diag(:)
    real(dp), intent(in), optional :: SUPPORT
    real(dp) :: SUP
    real(dp), parameter :: MINSUPPORT = 1e-4_dp

    if (present(SUPPORT)) then
        SUP = SUPPORT
    else
        SUP = MINSUPPORT
    end if
    n=size(M,DIM=1)
    if (n/=size(M,DIM=2)) call MpiStop('CutSkyAsymm_GetSupportedModesPol: assumed square')

    allocate(diag(n))
    nfound = n
    call Matrix_CDiagonalize_Partial(M, diag, n, SUP,1e100_dp, nfound)
    print *,'n modes = ', nfound, ' of ', n
    print *,'n modes > 0.999', count(diag(1:nfound)>0.999)
    print *,'n modes > 0.99', count(diag(1:nfound)>0.99)
    print *,'n modes > 0.9', count(diag(1:nfound)>0.9)
    print *, 'min/max diag = ', minval(diag(1:nfound)),maxval(diag(1:nfound))
    print *,'sum diag^2 = ', sum(diag(1:nfound)**2)

    allocate(Proj%WComp(n,nfound))
    Proj%nr = nfound
    Proj%nl = n
    Proj%WComp = M(:,1:nfound)
    allocate(proj%RootDiag(nfound))
    proj%RootDiag = sqrt(diag(1:nfound))
    deallocate(diag)

    end subroutine CutSkyAsymm_GetSupportedModesPol


    subroutine CutSkyAsymm_GetModeMatrixFromMapAlm(A, AP, Proj, ProjPol,lmax, MINSUPPORT, MapModes, MapA)
    Type(healpixAlm):: A, AP
    integer, intent(in) :: lmax
    real(dp), intent(in) :: MINSUPPORT
    Type(AsymmCouplings) :: Asymm
    Type(ProjMat) :: Proj
    Type(ProjMatPol) :: ProjPol
    real(dp), allocatable :: diag(:)
    real(dp), allocatable :: VT(:,:)
    integer i,m,n, npol, nfound
    Type(HealpixAlm),optional :: MapA
    Type (VecArray),optional :: MapModes(3)
    real(SP), dimension(:), allocatable :: Vec, VecE, VecB
    complex(dp), allocatable :: WComp(:,:)
    complex(dp) AMode

    call CutSkyAsymm_GetCoupling(Asymm, A, AP, lmax, plusonly = .false.)

    print *,'get temp supported'
    call CutSkyAsymm_GetSupportedModes(Asymm%WAsymm, Proj, MINSUPPORT)
    deallocate(Asymm%WAsymm)

    if (present(MapModes)) then
        allocate(vec((lmax+1)**2))
        call HealpixAlm2Vector(MapA, vec, lmax, 1)

        allocate(MapModes(1)%V(Proj%nr))
        do i=1, Proj%nr
            MapModes(1)%V(i) = dot_product(vec,Proj%M(:,i))/Proj%RootDiag(i)
        end do
        deallocate(vec)
    end if

    if (asymm_pol) then
        print *,'get pol supported'

        if (present(MapModes)) then
            allocate(vecE((lmax+1)**2-4))
            call HealpixAlm2Vector(MapA, vecE, lmax, 2)
            allocate(vecB((lmax+1)**2-4))
            call HealpixAlm2Vector(MapA, vecB, lmax, 3)
        end if

        if (.not. EBOnly) then
            print *,'doing not EBOnly'
            !Diagonalize (W_+ W_-, -W_- W_+) by diagonalizing hermitian W_+ + iW_i
            call CutSkyAsymm_SupportedModesPolFromCoupling(ASymm, ProjPol, MINSUPPORT)

            print *,'do diagonalize'

            if (present(MapModes)) then
                allocate(MapModes(2)%V(ProjPol%nr))
                allocate(MapModes(3)%V(ProjPol%nr))

                do i=1, ProjPol%nr
                    !AMode = (E+iB)^\dag
                    AMode = dot_product(cmplx(vecE,vecB),ProjPol%WComp(:,i))/ProjPol%RootDiag(i)
                    !Dot_product does SUM (CONJG (vector_a)*vector_b).

                    MapModes(2)%V(i)=real(AMode)
                    MapModes(3)%V(i)=-aimag(AMode)
                    !          MapModes(2)%V(i) = dot_product(vecE,ProjPol%WP_Proj%M(:,i))/ProjPol%RootDiag(i) &
                    !                          + dot_product(vecB,ProjPol%WM_Proj(:,i))/ProjPol%RootDiag(i)
                    !          MapModes(3)%V(i) = dot_product(vecB,ProjPol%WP_Proj%M(:,i))/ProjPol%RootDiag(i) &
                    !                          - dot_product(vecE,ProjPol%WM_Proj(:,i))/ProjPol%RootDiag(i)
                end do
                deallocate(vecE,vecB)
            end if


            !   !SVD of (W_+ W_-)
            !    allocate(VT(m,2*n))
            !    VT(1:m,1:n) = Asymm%WPAsymm
            !    VT(1:m,n+1:2*n) = Asymm%WMAsymm
            !    deallocate(Asymm%WMAsymm)
            !
            !    call Matrix_SVD_VT(VT,m,2*n,Diag,Asymm%WPAsymm )
            !    npol = count(diag>WellSupportedpol)
            !    print *,'num pol modes = ', npol
            !    ProjPol%nl = m
            !    ProjPol%nr = npol
            !    allocate(ProjPol%RootDiag(npol))
            !    allocate(ProjPol%WP_proj%M( ProjPol%nl,ProjPol%nr))
            !    ProjPol%WP_proj%M = transpose(VT(1:npol,1:m))
            !    ProjPol%RootDiag = 1
            !
            !    allocate(ProjPol%WM_proj(ProjPol%nl,ProjPol%nr))
            !    ProjPol%WM_proj = transpose(VT(1:npol,n+1:2*n))
            !    deallocate(VT)
            !
            !    ProjPol%WP_dataProj%nl = m
            !    ProjPol%WP_dataProj%nr = npol
            !    allocate(ProjPol%WP_dataProj%M(ProjPol%nl,ProjPol%nr))
            !    ProjPol%WP_dataProj%M = Asymm%WPAsymm(:,1:npol)
            !    deallocate(Asymm%WPAsymm)
            !    allocate(ProjPol%WP_dataProj%RootDiag(npol))
            !    ProjPol%WP_dataProj%RootDiag = diag(1:npol)

            !    allocate(MapModes(2)%V(ProjPol%nr))
            !    allocate(MapModes(3)%V(ProjPol%nr))
            !    do i=1, ProjPol%nr
            !          MapModes(2)%V(i) = dot_product(vecE,ProjPol%WP_dataProj%M(:,i))/diag(i)
            !          MapModes(3)%V(i) = dot_product(vecB,ProjPol%WP_dataProj%M(:,i))/diag(i)
            !    end do
            !    deallocate(vecE,vecB)
            !    deallocate(diag)
        else
            !Method if want E/B only
            stop 'removed'

            !           print *,'doing EB supported only'
            !            call CutSkyAsymm_GetSupportedModes(Asymm%WPAsymm, ProjPol%WP_proj, WellSupportedpol)
            !            deallocate(Asymm%WPAsymm)
            !            allocate(ProjPol%WM_proj( ProjPol%nl,ProjPol%nr))
            !
            !            print *,'get E/B coupling (WM_proj)'
            !
            !            ! D^{-1/2} U^T (iW_-) (all transpose)
            !            call Matrix_Mult(Asymm%WMAsymm,ProjPol%WP_proj%M, ProjPol%WM_proj)
            !            do i = 1, ProjPol%nr
            !            !minus because WMAsymm is antisymmtric
            !            ProjPol%WM_proj(:,i) = -ProjPol%WM_proj(:,i)/ProjPol%RootDiag(i)
            !            end do
            !           deallocate(Asymm%WMAsymm)
            !
            !
            !         allocate(MapModes(2)%V(ProjPol%nr))
            !         allocate(MapModes(3)%V(ProjPol%nr))
            !         do i=1, ProjPol%nr
            !          MapModes(2)%V(i) = dot_product(vecE,ProjPol%WP_proj%M(:,i))/ProjPol%RootDiag(i)
            !          MapModes(3)%V(i) = dot_product(vecB,ProjPol%WP_proj%M(:,i))/ProjPol%RootDiag(i)
            !         end do
            !         deallocate(vecE,vecB)
            !         ProjPol%WP_dataProj=ProjPol%WP_proj
        end if
    end if

    end subroutine CutSkyAsymm_GetModeMatrixFromMapAlm

    end  module CutSkyAsymm
