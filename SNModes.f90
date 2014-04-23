!Decomposition of T and pol maps into S/N modes
!Note some functions are currently on T
!Al, tidied up to this module March 2009

module SNModes
 use CutSkyAsymm
 use AMLutils
 use HealpixObj
 use spinalm_tools
 use MatrixUtils
 implicit none

 Type TSNModeOptions
   integer :: l_exact, l_low
   real(dP) :: SN_cut  !e.g.  1e-3, where we no longer care about losing signal at l < l_exact
   logical :: want_pol, want_cov
   logical :: pol_weights
   logical :: ProjectMonoAndDipole
   real(dp) :: KeepSupport  !e.g. 0.99, degree of support required for modes kept at l < l_low
   real(dp) :: ENoiseFac !factor by which pol noise assumed to be larger than T   
   logical :: keep_pseudo_noise_cov
   real(dp) :: fake_noise !for adding to modes when calculating S/N projection
   integer :: l_min_exact !normally 2, but if want to recover diople can set power spectrum for l=1 and set l_min_exact=1
  end Type TSNModeOptions

 
 Type TSNModes
   integer :: nmodes , nmodespol     
   Type(VecArray) :: SNModes(3) 
 end Type TSNModes   
    
 Type TSNModeMat
   logical :: has_pol
   integer :: totmodes
   !this and above a dagger of the actual coupling matrix, 
   ! e.g. X_s = DataProjPol^\dagger X_tilde or 
   ! or X_s^\dagger = X_tilde^\dagger DataProjPol
   Type(ProjMat) :: TheoryProj, DataProj    
   Type(ProjMatPol) :: TheoryProjPol, DataProjPol  
   Type(TSNModeOptions) :: MakeOptions

   logical :: has_cov
   !Information for covariance if requested
   !Note that currently calculate low l cov using input fid Cl; so should be using beam-smoothed C_l when then used in likelihood code
   real(dp) :: highlScaleT, highlScaleE, highlScaleB, highlScaleC
   Type(TCovMat) :: BigNoiseCov !big matrix '(T, E, B)' with noise covariance of all modes
   Type(TCovMat) :: OtherlNoiseCov !noise from high l

  !Optional extra data
   Type(AsymmCouplings) :: PseudoNoiseCov  !Temp only, need to back project modes onto sky modes (e.g. pictures)

 end Type TSNModeMat
 
 
contains


  subroutine SNModes_GetSNModeMats(H, SN, Opts, WeightMap,WeightMapPol, NoiseMap, PFid)
    !PFid is fiducial theory spectrum with any smoothing already applied
         Type (HealpixInfo) :: H
         Type (TSNModeOptions), intent(in) :: Opts
         Type (TSNModeMat) :: SN
         Type (HealpixMap) :: WeightMap, WeightMapPol, NoiseMap
         Type (HealpixPower) :: PFid
         
         Type (HealpixAlm) MaskA, MaskAP, NoiseA, NoiseAP
         Type (HealpixMap) :: NWMap
         
         Type(TCovMat) :: JustNoiseCov, NoiseCov, RotCov, Cov, HighlNoiseCov, ModeNoise
         Type(ProjMat) :: SNProj, Proj, HasSignalProj
         real(dp), allocatable :: Linv(:,:), tmpSN(:,:), TmpDataProj(:,:)

         Type(TComplexCovMat) :: JustPolNoiseCov, PolNoiseCov, PolCov, PolHighlNoiseCov, PolModeNoise
         Type(ProjMatPol) :: SNProjPol, ProjPol, HasSignalProjPol
         complex(dp), allocatable :: PolLinv(:,:), PoltmpSN(:,:) , TmpDataProjPol(:,:)        
         real(dp) :: NoiseScale
         real(dp) StTime
         integer i
         
         SN%MakeOptions = Opts
         print *,'SNModes_GetSNModeMats'
                    
         call HealpixMap2Alm(H,WeightMap, MaskA, Opts%l_low*2)
         call HealpixMap_Assign(NWMap,WeightMap)
  
         NWMap%TQU(:,1) = WeightMap%TQU(:,1)**2*NoiseMap%TQU(:,1)
         NoiseScale = (HO_fourpi/dble(NoiseMap%npix))  !Normalization for noise calculated like coupling matrix
         call HealpixMap2Alm(H,NWMap, NoiseA, Opts%l_low*2)
         call HealpixMap_Free(NWMap)
         
         if (Opts%want_pol) then
          if (Opts%pol_weights) then
             print *,'SNModes_GetSNModeMats doing pol'
             call HealpixMap2Alm(H,WeightMapPol, MaskAP, Opts%l_low*2)
             call HealpixMap_Assign(NWMap,WeightMapPol)      
             NWMap%TQU(:,1) = WeightMapPol%TQU(:,1)**2*NoiseMap%TQU(:,1)
             call HealpixMap2Alm(H,NWMap, NoiseAP, Opts%l_low*2)
             call HealpixMap_Free(NWMap)
           else
            NoiseAP = NoiseA
            MaskAP = MaskA
           endif                 
         end if

         print *,'SNModes_GetModeMatrixFromMapAlm'
         StTime = GeteTime()  
         call healpix_sleepMPI(H)
#ifdef MKLTHREADS
        call mkl_set_num_threads(H%MPISize)
#endif
         !Proj%M(nin, nsupported) is U for modes with nearly full support
         call CutSkyAsymm_GetModeMatrixFromMapAlm(MaskA, MaskAP, Proj, ProjPol, Opts%l_low, Opts%KeepSupport)
         print *,'CutSkyAsymm_GetModeMatrixFromMapAlm time',  GeteTime()   - StTime
         call HealpixAlm_Free(MaskA)
         if (opts%pol_weights)  call HealpixAlm_Free(MaskAP)
    !Get Noise    
        StTime = GeteTime()  
        print *,'Get noise covariance'
        call CutSkyAsymm_GetNoiseCovariance(Proj,ProjPol, JustNoiseCov, JustPolNoiseCov, &
                    NoiseA, NoiseAP, Opts%l_low, SN%PseudoNoiseCov)
        call HealpixAlm_Free(NoiseA)
        if (opts%pol_weights)  call HealpixAlm_Free(NoiseAP)
        
        JustNoiseCov%C = JustNoiseCov%C * NoiseScale        
        JustPolNoiseCov%C = JustPolNoiseCov%C * (Opts%ENoiseFac*NoiseScale)
        if (Opts%keep_pseudo_noise_cov) then
         SN%PseudoNoiseCov%WAsymm =  SN%PseudoNoiseCov%WAsymm * NoiseScale 
         if (Opts%want_pol) call mpiStop('not done pol noise cov')
        else 
         deallocate(SN%PseudoNoiseCov%WASymm)
         nullify(SN%PseudoNoiseCov%WASymm)
        end if
        
        print *,'CutSkyAsymm_GetNoiseCovariance time',  GeteTime()   - StTime

        if (Opts%fake_noise/=0) then
          print *,'adding fake noise to diagonal noise var'
          !Note not used for variance, just for projection
          do i=1,Proj%nr
            JustNoiseCov%C(i,i) = JustNoiseCov%C(i,i) + Opts%fake_noise
          end do         
         end if

    !add noise from high l
        print *,'get high l covariance'
        StTime = GeteTime()
        !HighlNoiseCov is only used for projecting modes; to calculate actual variance don't include dipole/monopole large contribution
        call CutSkyAsymm_GetCovariance(Proj, HighlNoiseCov, PFid, Opts%l_exact+1, Opts%l_low, Opts%ProjectMonoAndDipole)
        allocate(noiseCov%C(Proj%nr,Proj%nr))
        NoiseCov%C = JustNoiseCov%C + HighlNoiseCov%C
        print *,'CutSkyAsymm_GetCovariance high l time',  GeteTime()   - StTime
        deallocate(HighlNoiseCov%C)

        StTime = GeteTime() 
        call CutSkyAsymm_GetCovariance(Proj, Cov, PFid, Opts%l_min_exact,Opts%l_exact,.false.)
        print *,'CutSkyAsymm_GetCovariance low l time',  GeteTime()   - StTime

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

        allocate(rotCov%C(Proj%nr,Proj%nr))
        StTime = GeteTime() 
        call Matrix_RotateSymm(Cov%C, LInv, Proj%nr,  rotCov%C, triangular = .true.)
        deallocate(Cov%C)
        print *,'Rotate time',  GeteTime()   - StTime
        StTime = GeteTime() 
        !Get modes of Signal/(signal+ noise)
        call CutSkyAsymm_GetSupportedModes(rotCov%C, HasSignalProj, Opts%SN_cut) 
        print *,'CutSkyAsymm_GetSupportedModes sigal/(signal+noise)',  GeteTime()   - StTime
        deallocate(rotCov%C)
      
        deallocate(HasSignalProj%RootDiag) !Don't need the diagonal since just doing projecting into right basis
        nullify(HasSignalProj%RootDiag)
        
        StTime = GeteTime() 
        SNProj%nl = Proj%nr
        SNProj%nr=HasSignalProj%nr
        allocate(SNProj%M(SNProj%nl,SNProj%nr))
        !SNProj(polix)%M = Proj^T L^{-1}
        call Matrix_Mult(LInv,HasSignalProj%M,SNProj%M) 
        deallocate(LInv)
        deallocate(HasSignalProj%M)
 
! Get projected noise

        if (Opts%want_cov) then
         StTime = GeteTime() 
         allocate(ModeNoise%C(HasSignalProj%nr,HasSignalProj%nr))
!        !tmpSN = <nn^T> = [ Proj^T L^{-1}] NoiseCov [ Proj^T L^{-1}]^T
         call Matrix_RotateSymm(JustNoiseCov%C, SNProj%M, HasSignalProj%nr,  ModeNoise%C)
         print *,'noise part of S/N mode cov',  GeteTime()   - StTime
       end if
       deallocate(JustNoiseCov%C)

!Get theory coupling matrix
        StTime = GeteTime() 
        allocate(tmpSN(HasSignalProj%nl,HasSignalProj%nr))
        do i = 1, Proj%nr
           tmpSN(i,:) = SNProj%M(i,:)*Proj%RootDiag(i)
        end do

!Get data coupling
        allocate(TmpDataProj(HasSignalProj%nl,HasSignalProj%nr))
        do i=1, Proj%nr
            TmpDataProj(i,:) = SNProj%M(i,:)/Proj%RootDiag(i)
        end do
        deallocate(SNProj%M)

        print *,'getting full theory proj'
        SN%TheoryProj%nl = Proj%nl
        SN%TheoryProj%nr = HasSignalProj%nr
        allocate(SN%TheoryProj%M(SN%TheoryProj%nl,SN%TheoryProj%nr))
        call Matrix_Mult(Proj%M, tmpSN,SN%TheoryProj%M)
        deallocate(tmpSN)
  
        print *,'getting big data proj'
        SN%DataProj%nl = Proj%nl
        SN%DataProj%nr = SN%TheoryProj%nr
        allocate(SN%DataProj%M(SN%DataProj%nl,SN%DataProj%nr))
        call Matrix_Mult(Proj%M, tmpDataProj,SN%DataProj%M)
        deallocate(tmpDataProj)
     
        deallocate(Proj%M)
       
 !!!!Doing polarization       
       if (Opts%want_pol) then      
    
                print *,'Doing Polarization'        

                !for projection into S+N eigenmodes just approx pol cov as with no mixing from EE
                StTime = GeteTime()  
                !Get real part
                call CutSkyAsymm_GetUnmixedPolCovariance(ProjPol, PolHighlNoiseCov, PFid, Opts%l_exact+1, Opts%l_low)
                allocate(PolnoiseCov%C(ProjPol%nr,ProjPol%nr))
                PolNoiseCov%C = JustPolNoiseCov%C + PolHighlNoiseCov%C
                print *,'CutSkyAsymm_GetUnmixedPolCovariance high l time',  GeteTime()   - StTime
                deallocate(PolHighlNoiseCov%C)

                StTime = GeteTime() 
                call CutSkyAsymm_GetUnmixedPolCovariance(ProjPol, PolCov, PFid, 2, Opts%l_exact)
                print *,'CutSkyAsymm_GetUnmixedPolCovariance low l time',  GeteTime()   - StTime

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
                call CutSkyAsymm_GetSupportedModesPol(PoltmpSN, HasSignalProjPol,Opts%SN_cut) 
                print *,'CutSkyAsymm_GetSupportedModesPol',  GeteTime()   - StTime
                
                deallocate(PoltmpSN)
                deallocate(HasSignalProjPol%RootDiag)
                nullify(HasSignalProjPol%RootDiag)

                StTime = GeteTime() 

                SNProjPol%nl = ProjPol%nr
                SNProjPol%nr= HasSignalProjPol%nr
                allocate(SNProjPol%WComp(SNProjpol%nl,SNProjPol%nr))
                !SNProj(polix)%M = Proj^\dag L^{-1}
                call Matrix_CMult(PolLInv,HasSignalProjPol%WComp,SNProjPol%WComp) 
                deallocate(PolLInv)
                deallocate(HasSignalProjPol%WComp)
                
            ! Get projected noise

                if (Opts%want_cov) then
                 StTime = GeteTime() 
                 print *,'Getting pol modenoise'
                 allocate(PolModeNoise%C(HasSignalProjPol%nr,HasSignalProjPol%nr))
            !        !tmpSN = <nn^T> = [ Proj^T L^{-1}] NoiseCov [ Proj^T L^{-1}]^T
                 call Matrix_CRotateSymm(JustPolNoiseCov%C, SNProjPol%WComp, SNProjPol%nr,  PolModeNoise%C)
               
                 print *,'pol noise cov',  GeteTime()   - StTime
                end if
                deallocate(JustPolNoiseCov%C) 
           
            !Get theory coupling matrix
                StTime = GeteTime() 
                allocate(PoltmpSN(HasSignalProjPol%nl,HasSignalProjPol%nr))
                do i = 1, ProjPol%nr
                    PoltmpSN(i,:) = SNProjPol%WComp(i,:) * ProjPol%RootDiag(i)
                end do

!Get data coupling
                allocate(TmpDataProjPol(HasSignalProjPol%nl,HasSignalProjPol%nr))
                do i=1, ProjPol%nr
                    TmpDataProjPol(i,:) = SNProjPol%WComp(i,:) / ProjPol%RootDiag(i)
                end do
                deallocate(SNProjPol%WComp)

                print *,'getting big theory proj'
                SN%TheoryProjPol%nl = ProjPol%nl
                SN%TheoryProjPol%nr = HasSignalProjPol%nr
                allocate(SN%TheoryProjPol%WComp(SN%TheoryProjPol%nl,SN%TheoryProjPol%nr))
                call Matrix_CMult(ProjPol%WComp, PoltmpSN,SN%TheoryProjPol%WComp)
                deallocate(PoltmpSN)

                allocate(SN%TheoryProjPol%EProj(SN%TheoryProjPol%nl,SN%TheoryProjPol%nr))
                allocate(SN%TheoryProjPol%BProj(SN%TheoryProjPol%nl,SN%TheoryProjPol%nr))
                SN%TheoryProjPol%EProj = real(SN%TheoryProjPol%WComp)  !'E'^T = E^T EProj + B^T BProj
                SN%TheoryProjPol%BProj = aimag(SN%TheoryProjPol%WComp)
                deallocate(SN%TheoryProjPol%WComp)
                nullify(SN%TheoryProjPol%WComp)
             !I've not tested this yet
                print *,'getting big pol data proj'
                SN%DataProjPol%nl = ProjPol%nl
                SN%DataProjPol%nr = SN%TheoryProjPol%nr
                allocate(SN%DataProjPol%WComp(SN%DataProjPol%nl,SN%DataProjPol%nr))
                call Matrix_CMult(ProjPol%WComp, tmpDataProjPol,SN%DataProjPol%WComp)
                deallocate(ProjPol%WComp)
                nullify(ProjPol%WComp)
                deallocate(tmpDataProjPol)
      end if

!!!Big vectors
       SN%totmodes = SN%TheoryProj%nr
       if (Opts%want_pol) then
        SN%totmodes= SN%totmodes+ SN%TheoryProjPol%nr*2
       end if

       SN%has_cov = Opts%want_cov
       
       if (Opts%want_cov) then
  
           SN%highlScaleT = PFid%Cl(Opts%l_exact+1,C_T)

           print *,'get full high l covariance'         
           call CutSkyAsymm_GetFullCovariance(SN%TheoryProj,SN%TheoryProjPol, SN%OtherLNoiseCov, PFid, Opts%l_exact+1,Opts%l_low)
       
           print *,'Get full covariance'
           
           allocate(SN%BigNoiseCov%C(SN%totmodes,SN%totmodes))
           SN%BigNoiseCov%C=0
           SN%BigNoiseCov%C(1:SN%TheoryProj%nr,1:SN%TheoryProj%nr) = ModeNoise%C
           deallocate(ModeNoise%C)
           if (Opts%want_pol) then
           
            SN%highlScaleB = PFid%Cl(Opts%l_exact+1,C_B)
            SN%highlScaleE = PFid%Cl(Opts%l_exact+1,C_E)
            SN%highlScaleC = PFid%Cl(Opts%l_exact+1,C_C)
           
            SN%BigNoiseCov%C(SN%TheoryProj%nr+1:SN%TheoryProj%nr+SN%TheoryProjPol%nr, &
                             SN%TheoryProj%nr+1:SN%TheoryProj%nr+SN%TheoryProjPol%nr) =  real(PolModeNoise%C)
           
            SN%BigNoiseCov%C(SN%TheoryProj%nr+SN%TheoryProjPol%nr+1:SN%TheoryProj%nr+2*SN%TheoryProjPol%nr,&
                  SN%TheoryProj%nr+SN%TheoryProjPol%nr+1:SN%TheoryProj%nr+2*SN%TheoryProjPol%nr) = real(PolModeNoise%C)

            !BE           
            SN%BigNoiseCov%C(SN%TheoryProj%nr+SN%TheoryProjPol%nr+1:SN%TheoryProj%nr+2*SN%TheoryProjPol%nr,&
                  SN%TheoryProj%nr+1:SN%TheoryProj%nr+SN%TheoryProjPol%nr) = aimag(PolModeNoise%C)
           
           !EB    
            SN%BigNoiseCov%C(SN%TheoryProj%nr+1:SN%TheoryProj%nr+SN%TheoryProjPol%nr, &
                    SN%TheoryProj%nr+SN%TheoryProjPol%nr+1:SN%TheoryProj%nr+2*SN%TheoryProjPol%nr) = - aimag(PolModeNoise%C)

            deallocate(PolModeNoise%C)
           end if
           print *,'Get full high l covariance time',  GeteTime()   - StTime

       end if
#ifdef MKLTHREADS
        call mkl_set_num_threads(1)
#endif
       call healpix_wakeMPI

  end subroutine SNModes_GetSNModeMats
  
  subroutine SNModes_GetSNModes(H, SN, Opts, M, Modes)
  !M is the masked map
    Type (HealpixInfo):: H
    Type (TSNModeMat) :: SN
    Type (TSNModeOptions)  :: Opts
    Type (HealpixMap) :: M
    Type (TSNModes)   :: Modes
    Type (HealpixAlm) :: MapA
    real (sp) , dimension(:), allocatable :: vec, vecE, vecB
    complex (SP) :: amode
    integer i
             
    call HealpixMap2Alm(H,M,MapA,Opts%l_low)  
    
    allocate(vec((Opts%l_low+1)**2))
    call HealpixAlm2Vector(MapA, vec, Opts%l_low, 1)
    
    Modes%nmodes = SN%DataProj%nr    
    allocate(Modes%SNModes(1)%V(SN%DataProj%nr))
    do i=1, SN%DataProj%nr
         Modes%SNModes(1)%V(i) = dot_product(vec,SN%DataProj%M(:,i))
    end do
    deallocate(vec)
   
    if (Opts%want_pol) then
            
            Modes%nmodespol = SN%DataProjPol%nr    
       
            allocate(vecE((Opts%l_low+1)**2-4))
            call HealpixAlm2Vector(MapA, vecE, opts%l_low, 2)
            allocate(vecB((Opts%l_low+1)**2-4))
            call HealpixAlm2Vector(MapA, vecB, opts%l_low, 3)
    
            allocate(Modes%SNModes(2)%V(Modes%nmodespol))
            allocate(Modes%SNModes(3)%V(Modes%nmodespol))

            do i=1, Modes%nmodespol
            !AMode = (E+iB)^\dag
            !get conjugate of mode
              AMode = dot_product(cmplx(vecE,vecB),SN%DataProjPol%WComp(:,i))
               !Dot_product does SUM (CONJG (vector_a)*vector_b).

              Modes%SNModes(2)%V(i)=real(AMode)
              Modes%SNModes(3)%V(i)=-aimag(AMode)
            end do
            deallocate(vecE,vecB)
    
   else 
    Modes%nmodespol = 0             
   end if

   call HealpixAlm_Free(MapA)

  end subroutine SNModes_GetSNModes
  
  subroutine SNModes_WriteLikelihoodData(SN, Modes, exact_file)
  !For reading in by likelihood code
   Type(TSNModes) :: Modes
   Type(TSNModeMat) :: SN
   character(LEN=*), intent(in) :: exact_file
   integer unit
   integer i
       
       print *,'writing files'

        unit = new_file_unit()
        call CreateFile(exact_file,unit,'unformatted')
        if (SN%TheoryProj%nr/= Modes%nmodes) call mpiStop('SNModes_WriteLikelihoodData: error in no of modes')
        write(unit) SN%TheoryProj%nl, SN%TheoryProj%nr
        write(unit) Modes%SNModes(1)%V 
        do i=1, SN%TheoryProj%nr
         write(unit) SN%TheoryProj%M(1:SN%TheoryProj%nl,i)
        end do

       write(unit) SN%highlScaleT
        
       if (Modes%nmodespol /=0) then
            
            write(unit) SN%highlScaleE, SN%highlScaleC, SN%highlScaleB

            write(unit) SN%TheoryProjPol%nl, SN%TheoryProjPol%nr
            write(unit) Modes%SNModes(2)%V,Modes%SNModes(3)%V
 !           write(unit) DiagNoise(2)%V 
            do i=1, SN%TheoryProjPol%nr 
             write(unit) SN%TheoryProjPol%EProj(1:SN%TheoryProjPol%nl,i)
             write(unit) SN%TheoryProjPol%BProj(1:SN%TheoryProjPol%nl,i)
            end do
        end if


        do i=1,SN%totmodes
          write(unit) SN%BigNoiseCov%C(1:i,i)
          write(unit) SN%OtherLNoiseCov%C(1:i,i)
        end do
       i=252353
       write(unit) i
      close(unit)
    
  end subroutine SNModes_WriteLikelihoodData

  subroutine TSNModes_Read(Modes, exact_file)
   Type(TSNModes) :: Modes
   character(LEN=*), intent(in) :: exact_file
   integer unit
   
    unit = new_file_unit()
    call OpenFile(exact_file,unit,'unformatted')
    read(unit) Modes%nmodes, Modes%nmodespol
    allocate( Modes%SNModes(1)%V(Modes%nmodes))
    read(unit) Modes%SNModes(1)%V
    if (Modes%nmodespol /=0) then
     read(unit) Modes%nmodespol
     allocate(Modes%SNModes(2)%V(Modes%nmodespol))
     allocate(Modes%SNModes(3)%V(Modes%nmodespol))
     read(unit) Modes%SNModes(2)%V,Modes%SNModes(3)%V
    end if
    close(unit)

  end subroutine TSNModes_Read

   subroutine TSNModeMat_Read(SN, exact_file, feedback)
    Type(TSNModeMat) :: SN
    character(LEN=*), intent(in) :: exact_file
    logical :: feedback
    integer unit
    integer i, j
    
       unit = new_file_unit()
    
       call OpenFile(exact_file,unit,'unformatted')

        read(unit) SN%MakeOptions
        read(unit) SN%has_pol, SN%has_cov, SN%totmodes
        read(unit) SN%TheoryProj%nl, SN%TheoryProj%nr
        
        if (feedback) then
         print *,'Reading: '//trim(exact_file)
         print *, 'pol = ',SN%has_pol, 'has_cov = ',SN%has_cov, 'totmodes = ', SN%totmodes
        end if
        allocate(SN%TheoryProj%M(SN%TheoryProj%nl,SN%TheoryProj%nr))
        nullify(SN%TheoryProj%RootDiag)
        do i=1, SN%TheoryProj%nr
         read(unit) SN%TheoryProj%M(:,i)
        end do

        read(unit) SN%DataProj%nl, SN%DataProj%nr
        allocate(SN%DataProj%M(SN%DataProj%nl,SN%DataProj%nr))
        nullify(SN%DataProj%RootDiag)
        do i=1, SN%DataProj%nr
         read(unit) SN%DataProj%M(:,i)
        end do


        if (SN%has_pol) then
            read(unit) SN%TheoryProjPol%nl, SN%TheoryProjPol%nr
            nullify(SN%TheoryProjPol%RootDiag)
            allocate(SN%TheoryProjPol%EProj(SN%TheoryProjPol%nl,SN%TheoryProjPol%nr))
            allocate(SN%TheoryProjPol%BProj(SN%TheoryProjPol%nl,SN%TheoryProjPol%nr))
            do i=1, SN%TheoryProjPol%nr
             read(unit) SN%TheoryProjPol%EProj(:,i)
             read(unit) SN%TheoryProjPol%BProj(:,i)
            end do
            
            read(unit) SN%DataProjPol%nl, SN%DataProjPol%nr
            nullify(SN%DataProjPol%RootDiag)
            allocate(SN%DataProjPol%EProj(SN%DataProjPol%nl,SN%DataProjPol%nr))
            allocate(SN%DataProjPol%BProj(SN%DataProjPol%nl,SN%DataProjPol%nr))
            do i=1, SN%DataProjPol%nr
             read(unit) SN%DataProjPol%EProj(:,i)
             read(unit) SN%DataProjPol%BProj(:,i)
            end do
            
        end if
        
        if (SN%has_cov) then
            
            read(unit) SN%highlScaleT, SN%highlScaleE, SN%highlScaleB, SN%highlScaleC

            allocate(SN%bigNoiseCov%C(SN%totmodes,SN%totmodes))
            allocate(SN%OtherlNoiseCov%C(SN%totmodes,SN%totmodes))
            
            do i=1,SN%totmodes
              read(unit) SN%BigNoiseCov%C(1:i,i)
              read(unit) SN%OtherlNoiseCov%C(1:i,i)          
            end do
            do i=1,SN%totmodes
             do j=i+1, SN%totmodes
               SN%BigNoiseCov%C(j,i) = SN%BigNoiseCov%C(i,j) 
               SN%OtherlNoiseCov%C(j,i) = SN%OtherlNoiseCov%C(i,j) 
             end do
            end do 
            
        end if    

        if (SN%MakeOptions%Keep_pseudo_noise_cov) then
          read (unit) SN%PseudoNoiseCov%lmaxl, SN%PseudoNoiseCov%lmaxr
          allocate(SN%PseudoNoiseCov%WASymm((SN%PseudoNoiseCov%lmaxl+1)**2,(SN%PseudoNoiseCov%lmaxr+1)**2))
          read (unit) SN%PseudoNoiseCov%WASymm
        end if
          
        read (unit) i
        if (i/=252353) call MpiStop('TSNModeMat_Read: Bad file') 
        close(unit)
    
    end subroutine TSNModeMat_Read
    
 subroutine TSNModeMat_Write(SN, exact_file)
    Type(TSNModeMat) :: SN
    character(LEN=*), intent(in) :: exact_file
    integer unit
    integer i
    
        unit = new_file_unit()
        call CreateFile(exact_file,unit,'unformatted')
       
        write(unit) SN%MakeOptions
        write(unit) SN%has_pol, SN%has_cov, SN%totmodes
        write(unit) SN%TheoryProj%nl, SN%TheoryProj%nr
        do i=1, SN%TheoryProj%nr
         write(unit) SN%TheoryProj%M(:,i)
        end do
        
        write(unit) SN%DataProj%nl, SN%DataProj%nr
        do i=1, SN%DataProj%nr
         write(unit) SN%DataProj%M(:,i)
        end do

        if (SN%has_pol) then
            write(unit) SN%TheoryProjPol%nl, SN%TheoryProjPol%nr
            do i=1, SN%TheoryProjPol%nr
             write(unit) SN%TheoryProjPol%EProj(:,i)
             write(unit) SN%TheoryProjPol%BProj(:,i)
            end do

            write(unit) SN%DataProjPol%nl, SN%DataProjPol%nr
            do i=1, SN%DataProjPol%nr
             write(unit) SN%DataProjPol%EProj(:,i)
             write(unit) SN%DataProjPol%BProj(:,i)
            end do

        end if
        
        if (SN%has_cov) then
            write(unit) SN%highlScaleT, SN%highlScaleE, SN%highlScaleB, SN%highlScaleC
            do i=1,SN%totmodes
              write(unit) SN%BigNoiseCov%C(1:i,i)
              write(unit) SN%OtherlNoiseCov%C(1:i,i)          
            end do
        end if    
        
        if (SN%MakeOptions%Keep_pseudo_noise_cov) then
          write (unit) SN%PseudoNoiseCov%lmaxl, SN%PseudoNoiseCov%lmaxr
          write (unit) SN%PseudoNoiseCov%WASymm
        end if
          
        i=252353
        write(unit) i
        close(unit)
    
    end subroutine TSNModeMat_Write
    

    subroutine GetFullSkyAlmEstimator(SN,Modes, PFid, OutAlm, lmin_want, lmax_want) 
    !No attempt to optimize this, temp only at the mo; presumably only works for low lmax where invertible
     Type (TSNModes) :: Modes
     Type(TSNModeMat) :: SN
     Type(HealpixPower) :: PFid
     Type(HealpixAlm) :: OutAlm
     integer, intent(in) :: lmin_want, lmax_want
     Type (TCovMat) :: LowLFidCov     
     integer n_wanted     
     real, allocatable ::  AlmVec(:)
     real(dp), allocatable :: tmpSN(:,:), tmpvec(:)
  
     print *,'get full mode covariance'
     if (SN%MakeOptions%want_pol) call MpiStop('GetFullSkyAlmEstimator:not done pol exact cl')

     call CutSkyAsymm_GetFullCovariance(SN%TheoryProj,SN%TheoryProjPol, LowLFidCov, PFid, 2,SN%MakeOptions%l_exact)
       
     LowLFidCov%C =  LowLFidCov%C + SN%BigNoiseCov%C + SN%OtherLNoiseCov%C
    ! Get Linv = [L^{-1}]^T, so Linv is upper triangular
     call Matrix_Inverse(LowLFidCov%C)
     !Just get up to l_exact
     
     if (lmax_want > SN%MakeOptions%l_exact) call MpiStop('GetFullSkyAlmEstimator: wanting ell > l constructed to l_exact')
     n_wanted =(lmax_want+1)**2 - (lmin_want)**2
     
     allocate(tmpvec(n_wanted))
     print *,'get tmpvec'
     tmpvec = matmul(SN%TheoryProj%M((lmin_want)**2+1:(lmax_want+1)**2,:), matmul(LowLFidCov%C,Modes%SNModes(1)%V) )
     print *,'got tmpvec'
     allocate(tmpSN(n_wanted,n_wanted))
     tmpSN =  matmul(matmul( SN%TheoryProj%M((lmin_want)**2+1:(lmax_want+1)**2,:), LowLFidCov%C), &
            transpose( SN%TheoryProj%M((lmin_want)**2+1:(lmax_want+1)**2,:) ))
     print *, 'got tmpSN'
     call Matrix_Inverse(tmpSN)
     allocate(AlmVec( (lmax_want+1)**2 ))
     AlmVec = 0
     AlmVec((lmin_want)**2+1:(lmax_want+1)**2) = matmul(tmpSN, tmpvec) 
     deallocate(LowLFidCov%C)
     deallocate(tmpSN)
     deallocate(tmpvec)
     call HealpixAlm_Init(OutAlm,lmax_want,1)
     call HealpixVector2Alm(AlmVec, OutAlm, lmax_want,1)
     deallocate(AlmVec)
     
     end  subroutine GetFullSkyAlmEstimator   

end module SNModes