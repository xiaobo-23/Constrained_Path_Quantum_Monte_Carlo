!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the measurements for the FT-CPMC simulations, at \tau=BetaT or 
!                      \tau\in[0, BetaT].
! COMMENT: Static Measurements of CPMC simulations.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!   There are two different kinds of static measurements in Ft-CPMC:
!          (0) After choosing all the auxiliary fields, measure at \tau = BetaT point;
!          (1) Sweep from \tau=BetaT to \tau=0 and measure several times.
!             
!   PhyMeaStatBetaT --> Perform the static measurements at \tau=BetaT point;
!
!   PhyMeaStatM2One --> Perform the static measurements during the sweep from \tau=BetaT to \tau=0;
!   ProcMeaM2One --> Subroutine used to perform the process of measured data during the sweep from \tau=BetaT to \tau=0.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!########################################################################################################################
!########################################################################################################################
!############################### The STATIC measurements at tau==BetaT point ############################################
!############################### The STATIC measurements at tau==BetaT point ############################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine PhyMeaStatBetaT(NB, NSW) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PhyMeaStatBetaT(NB, NSW) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to measure the static physical quantities at \tau == BetaT point.
! KEYWORDS: Measure static physical quantities.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Measure the static physical quantities during the sweep of CPMC, at \tau == BetaT point.
!
!     During the measurement, we assign
!          GrnFRTmp00(i, j) = <c_i c_j^+>; GrnFRTmp11(i, j) = <c_i^+ c_j>
!
!     Input: NB  --> Integer index for BIN simulations;
!            NSW --> Integer index for sweeps;
!
!     Outpt: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
      use Observable
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB, NSW
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      integer I0, I1, I2, WalkIndx, Iwalk, NDim, SpnInd
      real(rp) ConfgConst, Rtp0, nTOccpHere
      complex(rp) Ztp0
      real(rp)   , allocatable :: Collect0(:, :, :)
      complex(rp), allocatable :: Collect1(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for static measurements ___________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      TimsStaMe = TimsStaMe + 1
      call system_clock(time1)
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for all the static quantities ________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Initializations for this subroutine _______________________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Initialize the Summation of weights of walkers ___________________
!________________________________________________________________________________________
      ConfgConst = 0.0_rp; WeightSumSwp = 0.0_rp
!________________________________________________________________________________________      
!_________________ (1) For IfMuTqmcNt == T or IfFixnT == T cases ________________________
!________________________________________________________________________________________
      nTOccpHere = 0.0_rp
!**************************************************************************************************     
!___________________ 1. Perform the static measurements for all random walkers ____________________
!**************************************************************************************************
      do WalkIndx = 1, NWkBt, +1
!________________________________________________________________________________________      
!_________________ (0) The integer index of the present random walker ___________________
!________________________________________________________________________________________
         Iwalk = IdptWkIndx(WalkIndx)
!________________________________________________________________________________________      
!_________________ (1) Compute the static Green's function matrix _______________________
!_____________________ exp(-dt*H_T/2) * GrnFRTmp00 * [exp(-dt*H_T/2)]^{-1} ______________
!________________________________________________________________________________________
         call dMat_Copy_QMC(GrnFunct(1, 1, 1, Iwalk), GrnFRTmp00(1, 1, 1))
         if(SymTrotDcp) call GrnFSymTrotR("HT", GrnFRTmp00(1, 1, 1))
!________________________________________________________________________________________      
!_________________ (2) Accumulate weights of walkers in single Process __________________
!________________________________________________________________________________________
         !!!!!!!!!! The sign/phase/weight related factor for measurement --> Walker weight
         ConfgConst = WghtProc(Iwalk)
         !!!!!!!!!! Accumulate the walker weight
         WeightSumSwp = WeightSumSwp + WghtProc(Iwalk)
!________________________________________________________________________________________      
!_________________ (3) Only measure the total density for fixed n_Occ ___________________
!_____________________ For IfMuTqmcNt == .true. or IfFixnT == .true. ____________________
!________________________________________________________________________________________
         if(IfMuTqmcNt .or. IfFixnT) then
            Rtp0 = 0.0_rp
            do I0 = 1, NumNS, +1
               Rtp0 = Rtp0 + 2.0_rp - GrnFRTmp00(I0, I0, 1) - GrnFRTmp00(I0, I0, 2)
            enddo
            nTOccpHere = nTOccpHere + ConfgConst * Rtp0 / dble(NumNS)
!________________________________________________________________________________________      
!_________________ (4) Measure energies, fillings and correlation functions _____________
!_____________________ For IfFixnT == .false. ___________________________________________
!________________________________________________________________________________________ 
         else
!____________________________________________________________________________      
!________________ [0] Process the static Green's functions __________________
!____________________________________________________________________________
            !!!!!!!!!! For IfFftEnPar == T, Compute GrFKCTmp00(k, q) = <c_k c_q^+>
            if(IfFftEnPar) then
               GrFKCTmp11 = rp_Zzero; GrFKCTmp00 = rp_Zzero
               call RghtWvfc_r2k_R(GrnFRTmp00(1, 1, 1), NumNS, NumNS, GrFKCTmp11(1, 1, 1))
               call LeftWvfc_r2k_C(NumNS, NumNS, GrFKCTmp11(1, 1, 1), GrFKCTmp00(1, 1, 1))
            end if
            !!!!!!!!!! Obtain GrnFRTmp11(i, j) = <c_i^+ c_j> = I - GrnFRTmp00^T
            !!!!!!!!!! Obtain GrFKCTmp11(k, q) = <c_k^+ c_q> = I - GrFKCTmp00^T
            if(IfFftEnPar) then
               GrnFRTmp11 = 0.0_rp; GrFKCTmp11 = rp_Zzero
               do SpnInd = 1, NmSpn, +1
                  do I2 = 1, NumNS, +1
                     do I1 = 1, NumNS, +1
                        GrnFRTmp11(I1, I2, SpnInd) = - GrnFRTmp00(I2, I1, SpnInd)
                        GrFKCTmp11(I1, I2, SpnInd) = - GrFKCTmp00(I2, I1, SpnInd)
                     enddo
                     GrnFRTmp11(I2, I2, SpnInd) = 1.0_rp   + GrnFRTmp11(I2, I2, SpnInd)
                     GrFKCTmp11(I2, I2, SpnInd) = rp_Z_One + GrFKCTmp11(I2, I2, SpnInd)
                  enddo
               enddo   
            else
               GrnFRTmp11 = 0.0_rp
               do SpnInd = 1, NmSpn, +1
                  do I2 = 1, NumNS, +1
                     do I1 = 1, NumNS, +1
                        GrnFRTmp11(I1, I2, SpnInd) = - GrnFRTmp00(I2, I1, SpnInd)
                     enddo
                     GrnFRTmp11(I2, I2, SpnInd) = 1.0_rp + GrnFRTmp11(I2, I2, SpnInd)
                  enddo
               enddo   
            end if
!____________________________________________________________________________      
!________________ [1] Measure energies, fillings and GrF matrices ___________
!____________________________________________________________________________
            if(IfFftEnPar) then
               !!!!!!!!!! The G(i, j) in r-space and G(k, q) matrix in k-space
               Rtp0 = ConfgConst; Ztp0 = cmplx(ConfgConst, 0.0_rp, rp)
               call DAXPY(NumNS*NumNS*2, Rtp0, GrnFRTmp11(1, 1, 1), 1, RlSpGrnFtSwp(1, 1, 1, 0), 1)
               call ZAXPY(NumNS*NumNS*2, Ztp0, GrFKCTmp11(1, 1, 1), 1, KSpGreenFSwp(1, 1, 1, 0), 1)
               !!!!!!!!!! Measure the energies and densities
               call ObStaEnrgy_H0FFTW(ConfgConst, GrnFRTmp11, GrFKCTmp11, EngOccCrFSwp(1, 0)) 
            else
               !!!!!!!!!! The G(i, j) in r-space
               call DAXPY(NumNS*NumNS*2, ConfgConst, GrnFRTmp11(1, 1, 1), 1, RlSpGrnFtSwp(1, 1, 1, 0), 1)
               !!!!!!!!!! Measure the energies and densities
               call ObStaEnrgy_SpnDcp(ConfgConst, GrnFRTmp00, GrnFRTmp11, EngOccCrFSwp(1, 0))
            end if
!____________________________________________________________________________      
!________________ [2] Measure n(k) and the paring matrices __________________
!____________________________________________________________________________
            if(IfFftEnPar) then
               call NkPairWvfc_H0FFTW(ConfgConst, GrFKCTmp00, GrFKCTmp11, NkSgleSwp(1, 1, 0), PairMtSwp(1, 1, 0))
            end if
!____________________________________________________________________________      
!________________ [3] Measure r-space correlations (PBC or OBC) _____________
!____________________________________________________________________________
            if(abs(PinSz) < rp_Eps) then
               call ObStaCrFct_SpnDcp(ConfgConst, GrnFRTmp00, GrnFRTmp11, RealSpCrFSwp(1, 1, 0))
            end if
         end if
      enddo
!**************************************************************************************************     
!___________________ 2. Combine all the random walkers from different processes ___________________
!**************************************************************************************************
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)
!________________________________________________________________________________________      
!_________________ (0) Process the weights of all random walkers ________________________
!________________________________________________________________________________________
      Rtp0 = 0.0_rp
      call MPI_ALLREDUCE(WeightSumSwp, Rtp0, 1, rp_MPI_REAL, MPI_SUM, acomm, ierr)
      WeightSumSwp = Rtp0
!________________________________________________________________________________________      
!_________________ (1) For IfMuTqmcNt == T or IfFixnT == T cases ________________________
!_____________________ only measure the number of fermions ______________________________
!________________________________________________________________________________________ 
      if(IfMuTqmcNt .or. IfFixnT) then
         Rtp0 = 0.0_rp
         call MPI_ALLREDUCE(nTOccpHere, Rtp0, 1, rp_MPI_REAL, MPI_SUM, acomm, ierr)
         nTOccpHere = Rtp0
!________________________________________________________________________________________      
!_________________ (2) For other cases, Measure all quantities __________________________
!________________________________________________________________________________________    
      else
!____________________________________________________________________________      
!________________ [0] Collect data of energies and GrF matrices _____________
!____________________________________________________________________________   
         !!!!!!!!!! The energies, fillings and density correlations
         NDim = 40
         allocate(Collect0(NDim, 1, 1)); Collect0 = 0.0_rp
         call MPI_ALLREDUCE(EngOccCrFSwp(1, 0), Collect0(1, 1, 1), NDim, rp_MPI_REAL, MPI_SUM, acomm, ierr)
         call dcopy(NDim, Collect0(1, 1, 1), 1, EngOccCrFSwp(1, 0), 1)
         if(allocated(Collect0)) deallocate(Collect0)
         !!!!!!!!!! The r-space single-particle G(i, j) = <c_i^+ c_j> matrix
         NDim = NumNC * NumNC * 2
         allocate(Collect0(NumNC, NumNC, 2)); Collect0 = 0.0_rp
         call MPI_REDUCE(RlSpGrnFtSwp(1, 1, 1, 0), Collect0(1, 1, 1), NDim, rp_MPI_REAL, MPI_SUM, amstr, acomm, ierr)
         call dcopy(NDim, Collect0(1, 1, 1), 1, RlSpGrnFtSwp(1, 1, 1, 0), 1)
         if(allocated(Collect0)) deallocate(Collect0)
         !!!!!!!!!! The k-space single-particle G(k, q) = <c_k^+ c_q> matrix
         if(IfFftEnPar) then
            NDim = NumNC * NumNC * 2
            allocate(Collect1(NumNC, NumNC, 2)); Collect1 = rp_Zzero
            call MPI_REDUCE(KSpGreenFSwp(1, 1, 1, 0), Collect1(1, 1, 1), NDim, rp_MPI_COMPLEX, MPI_SUM, amstr, acomm, ierr)
            call zcopy(NDim, Collect1(1, 1, 1), 1, KSpGreenFSwp(1, 1, 1, 0), 1)
            if(allocated(Collect1)) deallocate(Collect1)
         end if
!____________________________________________________________________________      
!________________ [1] Collect data of n(k) and the paring matrices __________
!____________________________________________________________________________
         if(IfFftEnPar) then
            !!!!!!!!!! The n(k) as momentum distribution
            NDim = NumNC*2
            allocate(Collect1(NumNC, 2, 1)); Collect1 = rp_Zzero
            call MPI_REDUCE(NkSgleSwp(1, 1, 0), Collect1(1, 1, 1), NDim, rp_MPI_COMPLEX, MPI_SUM, amstr, acomm, ierr)
            call zcopy(NDim, Collect1(1, 1, 1), 1, NkSgleSwp(1, 1, 0), 1)
            if(allocated(Collect1)) deallocate(Collect1)
            !!!!!!!!!! The paring matrices
            NDim = NumNC*NumNC
            allocate(Collect1(NumNC, NumNC, 1)); Collect1 = rp_Zzero
            call MPI_REDUCE(PairMtSwp(1, 1, 0), Collect1(1, 1, 1), NDim, rp_MPI_COMPLEX, MPI_SUM, amstr, acomm, ierr)
            call zcopy(NDim, Collect1(1, 1, 1), 1, PairMtSwp(1, 1, 0), 1)
            if(allocated(Collect1)) deallocate(Collect1)
         end if
!____________________________________________________________________________      
!________________ [2] Collect data of r-space correlations (PBC or OBC) _____
!____________________________________________________________________________
         if(abs(PinSz) < rp_Eps) then
            NDim = NmSitePair*40
            allocate(Collect0(NmSitePair, 40, 1)); Collect0 = 0.0_rp
            call MPI_REDUCE(RealSpCrFSwp(1, 1, 0), Collect0(1, 1, 1), NDim, rp_MPI_REAL, MPI_SUM, amstr, acomm, ierr)
            call dcopy(NDim, Collect0(1, 1, 1), 1, RealSpCrFSwp(1, 1, 0), 1)
            if(allocated(Collect0)) deallocate(Collect0)
         end if
      end if
      call MPI_Barrier(acomm, ierr)
#endif
!**************************************************************************************************     
!___________________ 3. Accumulate the results of physical observables ____________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Process the weights of all random walkers ________________________
!________________________________________________________________________________________
      !!!!!!!!!! Weight average (over all walkers) of the present sweep
      WeightAvgSwp = merge(WeightSumSwp/dble(NmWalkAllP), 1.0_rp, IfSwepReWt)
      !!!!!!!!!! Accumulate weight average of sweeps in a BIN
      WtMeanSumBIN = WtMeanSumBIN + WeightAvgSwp
!________________________________________________________________________________________      
!_________________ (1) For IfMuTqmcNt == T or IfFixnT == T cases ________________________
!_____________________ only the number of fermions ______________________________________
!________________________________________________________________________________________ 
      if(IfMuTqmcNt .or. IfFixnT) then
         WeightList(NObsStat) = merge(WeightAvgSwp, WeightSumSwp/dble(NmWalkAllP), IfMuTqmcNt)
         n_Occ_List(NObsStat) = nTOccpHere / WeightSumSwp
!________________________________________________________________________________________      
!_________________ (2) For other cases, for all quantities ______________________________
!________________________________________________________________________________________    
      else 
!____________________________________________________________________________      
!________________ [0] Results of energies and the GrF matrices ______________
!____________________________________________________________________________
         !!!!!!!!!! The energies, fillings and density correlations
         EngOccCrFSwp(:, 0) = EngOccCrFSwp(:, 0) / WeightSumSwp
         !!!!!!!!!! The r-space single-particle GrF matrix
         RlSpGrnFtSwp(:, :, :, 0) = RlSpGrnFtSwp(:, :, :, 0) / WeightSumSwp
         !!!!!!!!!! The k-space single-particle GrF matrix
         if(IfFftEnPar) then
            KSpGreenFSwp(:, :, :, 0) = KSpGreenFSwp(:, :, :, 0) / WeightSumSwp
         end if
!____________________________________________________________________________      
!________________ [1] Results of n(k) and the paring matrices _______________
!____________________________________________________________________________
         if(IfFftEnPar) then
            NkSgleSwp(:, :, 0) = NkSgleSwp(:, :, 0) / WeightSumSwp
            PairMtSwp(:, :, 0) = PairMtSwp(:, :, 0) / WeightSumSwp
         end if
!____________________________________________________________________________      
!________________ [2] Results of r-space correlations (PBC or OBC) __________
!____________________________________________________________________________ 
         if(abs(PinSz) < rp_Eps) then
            RealSpCrFSwp(:, :, 0) = RealSpCrFSwp(:, :, 0) / WeightSumSwp
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for static measurements ___________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeStaMe = TimeStaMe + TimeIntrvl(time1, time2)
      
   end subroutine PhyMeaStatBetaT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!############################### The STATIC measurements for sweep of tau=0 and tau=BetaT ###############################
!############################### The STATIC measurements for sweep of tau=0 and tau=BetaT ###############################
!########################################################################################################################
!########################################################################################################################

      
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine PhyMeaStatM2One(Iwalk, NB, NSW, NT) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PhyMeaStatM2One(Iwalk, NB, NSW, NT) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to measure the static physical quantities for sweep of tau=0 and tau=BetaT.
! KEYWORDS: Measure static physical quantities.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Measure the static physical quantities during the sweep of CPMC, for the [BetaT, 0]. sweep
!
!     Input: Iwalk --> Integer index for the present random walker.
!            NB    --> Integer index for BIN simulations;
!            NSW   --> Integer index for sweeps;
!            NT    --> Integer index for the time slice;
!
!     Outpt: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
      use Observable
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer Iwalk, NB, NSW, NT
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      integer I1, I2, SpnInd
      real(rp) Rtp0, ConfgConst
      complex(rp) Ztp0
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for static measurements ___________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      TimsB0Mea = TimsB0Mea + 1
      call system_clock(time1)
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for all the static quantities ________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Perform the static measurements for Iwalk-th random walker ________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Prepare the static Green's function matrices _____________________
!________________________________________________________________________________________
!____________________________________________________________________________      
!________________ [0] GrnFRTmp00 = <c_i c_j^+> from propagation _____________
!_______ GrnFRTmp00 --> exp(-dt*K/2) * GrnFRTmp00 * [exp(-dt*K/2)]^{-1} _____
!____________________________________________________________________________
      call dMat_Copy_QMC(GrnFunct(1, 1, 1, Iwalk), GrnFRTmp00(1, 1, 1))
      if(SymTrotDcp) call GrnFSymTrotR("H0", GrnFRTmp00(1, 1, 1))
!____________________________________________________________________________      
!________________ [1] For IfFftEnPar == T case, compute _____________________
!___________________  GrFKCTmp00(k, q) = <c_k c_q^+> ________________________
!____________________________________________________________________________
      if(IfFftEnPar) then
         GrFKCTmp11 = rp_Zzero; GrFKCTmp00 = rp_Zzero
         call RghtWvfc_r2k_R(GrnFRTmp00(1, 1, 1), NumNS, NumNS, GrFKCTmp11(1, 1, 1))
         call LeftWvfc_r2k_C(NumNS, NumNS, GrFKCTmp11(1, 1, 1), GrFKCTmp00(1, 1, 1))
      end if
!____________________________________________________________________________      
!________________ [2] GrnFRTmp11 = <c_i^+ c_j> = I - GrnFRTmp00^T ___________
!____________________________________________________________________________
      if(IfFftEnPar) then
         GrnFRTmp11 = 0.0_rp; GrFKCTmp11 = rp_Zzero
         do SpnInd = 1, NmSpn, +1
            do I2 = 1, NumNS, +1
               do I1 = 1, NumNS, +1
                  GrnFRTmp11(I1, I2, SpnInd) = - GrnFRTmp00(I2, I1, SpnInd)
                  GrFKCTmp11(I1, I2, SpnInd) = - GrFKCTmp00(I2, I1, SpnInd)
               enddo
               GrnFRTmp11(I2, I2, SpnInd) = 1.0_rp   + GrnFRTmp11(I2, I2, SpnInd)
               GrFKCTmp11(I2, I2, SpnInd) = rp_Z_One + GrFKCTmp11(I2, I2, SpnInd)
            enddo
         enddo   
      else
         GrnFRTmp11 = 0.0_rp
         do SpnInd = 1, NmSpn, +1
            do I2 = 1, NumNS, +1
               do I1 = 1, NumNS, +1
                  GrnFRTmp11(I1, I2, SpnInd) = - GrnFRTmp00(I2, I1, SpnInd)
               enddo
               GrnFRTmp11(I2, I2, SpnInd) = 1.0_rp + GrnFRTmp11(I2, I2, SpnInd)
            enddo
         enddo   
      end if
!________________________________________________________________________________________      
!_________________ (1) Weight of the present random walker ______________________________
!________________________________________________________________________________________ 
      ConfgConst = WghtProc(Iwalk)
!________________________________________________________________________________________      
!_________________ (2) Measure energies, fillings and correlation functions _____________
!________________________________________________________________________________________
!____________________________________________________________________________      
!________________ [0] Measure energies, fillings and GrF matrices ___________
!____________________________________________________________________________
      if(IfFftEnPar) then
         !!!!!!!!!! The G(i, j) in r-space and G(k, q) matrix in k-space
         Rtp0 = ConfgConst; Ztp0 = cmplx(ConfgConst, 0.0_rp, rp)
         call DAXPY(NumNS*NumNS*2, Rtp0, GrnFRTmp11(1, 1, 1), 1, RlSpGrnFtSwp(1, 1, 1, 1), 1)
         call ZAXPY(NumNS*NumNS*2, Ztp0, GrFKCTmp11(1, 1, 1), 1, KSpGreenFSwp(1, 1, 1, 1), 1)
         !!!!!!!!!! Measure the energies and densities
         call ObStaEnrgy_H0FFTW(ConfgConst, GrnFRTmp11, GrFKCTmp11, EngOccCrFSwp(1, 1)) 
      else
         !!!!!!!!!! The G(i, j) in r-space
         call DAXPY(NumNS*NumNS*2, ConfgConst, GrnFRTmp11(1, 1, 1), 1, RlSpGrnFtSwp(1, 1, 1, 1), 1)
         !!!!!!!!!! Measure the energies and densities
         call ObStaEnrgy_SpnDcp(ConfgConst, GrnFRTmp00, GrnFRTmp11, EngOccCrFSwp(1, 1))
      end if
!____________________________________________________________________________      
!________________ [1] Measure n(k) and the paring matrices __________________
!____________________________________________________________________________
      if(IfFftEnPar) then
         call NkPairWvfc_H0FFTW(ConfgConst, GrFKCTmp00, GrFKCTmp11, NkSgleSwp(1, 1, 1), PairMtSwp(1, 1, 1))
      end if
!____________________________________________________________________________      
!________________ [2] Measure r-space correlations (PBC or OBC) _____________
!____________________________________________________________________________
      if(abs(PinSz) < rp_Eps) then
         call ObStaCrFct_SpnDcp(ConfgConst, GrnFRTmp00, GrnFRTmp11, RealSpCrFSwp(1, 1, 1))
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for static measurements ___________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeB0Mea = TimeB0Mea + TimeIntrvl(time1, time2)
      
   end subroutine PhyMeaStatM2One
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ProcMeaM2One() 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ProcMeaM2One() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to process the measured results during the sweep of tau==Beta to tau==0.
! KEYWORDS: Process the measured static and dynamic results for Beta --> 0 sweep.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Process the measured static and dynamic results for the [BetaT, 0] sweep.
!
!     Input: (none); Outpt: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
      use Observable
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      integer NDim
      real(rp)   , allocatable :: Collect0(:, :, :)
      complex(rp), allocatable :: Collect1(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for static measurements ___________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      TimsDatPr = TimsDatPr + 1
      call system_clock(time1)
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for process the measured data ________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Combine all the random walkers from different processes ___________________
!**************************************************************************************************
#ifdef MPIPROCESS   
      call MPI_Barrier(acomm, ierr)
!________________________________________________________________________________________      
!_________________ (0) Collect data of energies and GrF matrices ________________________
!________________________________________________________________________________________   
      !!!!!!!!!! The energies, fillings and density correlations
      NDim = 40; 
      allocate(Collect0(NDim, 1, 1)); Collect0 = 0.0_rp
      call MPI_ALLREDUCE(EngOccCrFSwp(1, 1), Collect0(1, 1, 1), NDim, rp_MPI_REAL, MPI_SUM, acomm, ierr)
      call dcopy(NDim, Collect0(1, 1, 1), 1, EngOccCrFSwp(1, 1), 1)
      if(allocated(Collect0)) deallocate(Collect0)
      !!!!!!!!!! The r-space single-particle G(i, j) = <c_i^+ c_j> matrix
      NDim = NumNC * NumNC * 2
      allocate(Collect0(NumNC, NumNC, 2)); Collect0 = 0.0_rp
      call MPI_REDUCE(RlSpGrnFtSwp(1, 1, 1, 1), Collect0(1, 1, 1), NDim, rp_MPI_REAL, MPI_SUM, amstr, acomm, ierr)
      call dcopy(NDim, Collect0(1, 1, 1), 1, RlSpGrnFtSwp(1, 1, 1, 1), 1)
      if(allocated(Collect0)) deallocate(Collect0)
      !!!!!!!!!! The k-space single-particle G(k, q) = <c_k^+ c_q> matrix
      if(IfFftEnPar) then
         NDim = NumNC * NumNC * 2
         allocate(Collect1(NumNC, NumNC, 2)); Collect1 = rp_Zzero
         call MPI_REDUCE(KSpGreenFSwp(1, 1, 1, 1), Collect1(1, 1, 1), NDim, rp_MPI_COMPLEX, MPI_SUM, amstr, acomm, ierr)
         call zcopy(NDim, Collect1(1, 1, 1), 1, KSpGreenFSwp(1, 1, 1, 1), 1)
         if(allocated(Collect1)) deallocate(Collect1)
      end if
!________________________________________________________________________________________      
!_________________ (1) Collect data of n(k) and the paring matrices _____________________
!________________________________________________________________________________________  
      if(IfFftEnPar) then
         !!!!!!!!!! n(k) as momentum distribution
         NDim = NumNC * 2
         allocate(Collect1(NumNC, 2, 1)); Collect1 = rp_Zzero
         call MPI_REDUCE(NkSgleSwp(1, 1, 1), Collect1(1, 1, 1), NDim, rp_MPI_COMPLEX, MPI_SUM, amstr, acomm, ierr)
         call zcopy(NDim, Collect1(1, 1, 1), 1, NkSgleSwp(1, 1, 1), 1)
         if(allocated(Collect1)) deallocate(Collect1)
         !!!!!!!!!! The paring matrices
         NDim = NumNC * NumNC
         allocate(Collect1(NumNC, NumNC, 1)); Collect1 = 0.0_rp
         call MPI_REDUCE(PairMtSwp(1, 1, 1), Collect1(1, 1, 1), NDim, rp_MPI_COMPLEX, MPI_SUM, amstr, acomm, ierr)
         call zcopy(NDim, Collect1(1, 1, 1), 1, PairMtSwp(1, 1, 1), 1)
         if(allocated(Collect1)) deallocate(Collect1)
      end if
!________________________________________________________________________________________      
!_________________ (2) Collect data of r-space correlations (PBC or OBC) ________________
!________________________________________________________________________________________
      if(abs(PinSz) < rp_Eps) then
         NDim = NmSitePair * 40
         allocate(Collect0(NmSitePair, 40, 1)); Collect0 = 0.0_rp
         call MPI_REDUCE(RealSpCrFSwp(1, 1, 1), Collect0(1, 1, 1), NDim, rp_MPI_REAL, MPI_SUM, amstr, acomm, ierr)
         call dcopy(NDim, Collect0(1, 1, 1), 1, RealSpCrFSwp(1, 1, 1), 1)
         if(allocated(Collect0)) deallocate(Collect0)
      end if
!________________________________________________________________________________________      
!_________________ (3) Collect data of dynamic correlations (PBC or OBC) ________________
!________________________________________________________________________________________ 
      if(IfTAU) then
         NDim = (NumTauPnt+1) * NmSitePairTau * 40
         allocate(Collect0(0:NumTauPnt, NmSitePairTau, 40)); Collect0 = 0.0_rp
         call MPI_REDUCE(RealSpCrFTauSwp(0, 1, 1), Collect0(0, 1, 1), NDim, rp_MPI_REAL, MPI_SUM, amstr, acomm, ierr)
         call dcopy(NDim, Collect0(0, 1, 1), 1, RealSpCrFTauSwp(0, 1, 1), 1)
         if(allocated(Collect0)) deallocate(Collect0)
      end if
      call MPI_Barrier(acomm, ierr)
#endif
!**************************************************************************************************     
!___________________ 1. Accumulate the results of physical observables ____________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Results of energies and the GrF matrices _________________________
!________________________________________________________________________________________
      !!!!!!!!!! The energies, fillings and density correlations
      EngOccCrFSwp(:, 1) = EngOccCrFSwp(:, 1) / WeightSumSwp / dble(MeaM2One)
      !!!!!!!!!! The r-space single-particle GrF matrix
      RlSpGrnFtSwp(:, :, :, 1) = RlSpGrnFtSwp(:, :, :, 1) / WeightSumSwp / dble(MeaM2One)
      !!!!!!!!!! The k-space single-particle GrF matrix
      if(IfFftEnPar) then
         KSpGreenFSwp(:, :, :, 1) = KSpGreenFSwp(:, :, :, 1) / WeightSumSwp / dble(MeaM2One)
      end if      
!________________________________________________________________________________________      
!_________________ (1) Results of n(k) and the paring matrices __________________________
!________________________________________________________________________________________
      if(IfFftEnPar) then
         NkSgleSwp(:, :, 1) = NkSgleSwp(:, :, 1) / WeightSumSwp / dble(MeaM2One)
         PairMtSwp(:, :, 1) = PairMtSwp(:, :, 1) / WeightSumSwp / dble(MeaM2One)
      end if
!________________________________________________________________________________________      
!_________________ (2) Results of r-space static correlations (PBC or OBC) ______________
!________________________________________________________________________________________
      if(abs(PinSz) < rp_Eps) then
         RealSpCrFSwp(:, :, 1) = RealSpCrFSwp(:, :, 1) / WeightSumSwp / dble(MeaM2One)  
      end if
!________________________________________________________________________________________      
!_________________ (4) Results of r-space dynamic correlations (PBC or OBC) _____________
!________________________________________________________________________________________
      if(IfTAU) then
         RealSpCrFTauSwp(:, :, :) = RealSpCrFTauSwp(:, :, :) / WeightSumSwp
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for static measurements ___________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeDatPr = TimeDatPr + TimeIntrvl(time1, time2)
      
   end subroutine ProcMeaM2One
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$