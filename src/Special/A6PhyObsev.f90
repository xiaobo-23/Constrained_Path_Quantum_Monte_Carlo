!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL SPIN PINNING FIELDS; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform allocation of the observables before the whole CPMC simulations.
! COMMENT: CPMC Initialization process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   InitPhyObsev --> Subroutine to allocate the observables vectors and matrices used in CPMC;
!
!   CPMCObservInit --> Give initial values to the physical observables.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine InitPhyObsev()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitPhyObsev()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the initialization for observables calculated in CPMC simulations, 
!                mainly initiate some calculation parameters and allocate some vectors and matrices.
! KEYWORDS: Initialization of Observable module.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Initialization of Observable module, including:
!             (0) Initialization of some measurement parameters;
!             (1) Allocate the observable matrices and initializations.
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use Observable
      use StdInOutSt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer MeaTypeCnt
!______________________________________________________________________________________________________________     
!____________________________ Main calculations of Observables Initialization _________________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of initialization process _______________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'InitPhyObsev: Initialization of measured physical observables!')")
      end if
!**************************************************************************************************     
!___________________ 0. Initializations for this subroutine _______________________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) The number of measurements for BetaT and [B, 0] sweep ____________
!________________________________________________________________________________________
      MeaTypeCnt = merge(1, 0, IfM2OneMea)
!________________________________________________________________________________________  
!_________________ (1) Settings for the output folders for results ______________________
!________________________________________________________________________________________
      FileOutAdd(0) = "Output/"; FileOutAdd(1) = "Add_Output/"
      FileAddTxt(0) = ".txt"   ; FileAddTxt(1) = "_Add.txt"
!**************************************************************************************************  
!___________________ 1. Measure energies, occupations and local densities _________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Accumulate Energies and Occupations for every single BIN _________
!________________________________________________________________________________________
      allocate(EngOccCrFSwp(40, 0:MeaTypeCnt)); EngOccCrFSwp = 0.0_rp    
      allocate(EngOccCrFBIN(40, 0:MeaTypeCnt)); EngOccCrFBIN = 0.0_rp
!________________________________________________________________________________________  
!_________________ (1) Expectation energies for every term in the model _________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         allocate(EnHopt1(NmBIN, 0:MeaTypeCnt)); EnHopt1 = 0.0_rp
         allocate(EnHopt2(NmBIN, 0:MeaTypeCnt)); EnHopt2 = 0.0_rp
         allocate(EnHopt3(NmBIN, 0:MeaTypeCnt)); EnHopt3 = 0.0_rp
         allocate(EnZmFld(NmBIN, 0:MeaTypeCnt)); EnZmFld = 0.0_rp
         allocate(EnPinSz(NmBIN, 0:MeaTypeCnt)); EnPinSz = 0.0_rp
         allocate(EnDopCh(NmBIN, 0:MeaTypeCnt)); EnDopCh = 0.0_rp
         allocate(EnHubbU(NmBIN, 0:MeaTypeCnt)); EnHubbU = 0.0_rp
         allocate(EnTotal(NmBIN, 0:MeaTypeCnt)); EnTotal = 0.0_rp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         allocate(EnSinusoidalPinSz(NmBin, 0:MeaTypeCnt)); EnSinusoidalPinSz = 0.0_rp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if(EkDispType >= 1) then
            allocate(EpOvEfg(NmBIN, 2, 0:MeaTypeCnt)); EpOvEfg = 0.0_rp
         end if
      end if
!________________________________________________________________________________________  
!_________________ (2) First-order derivative of model parameters _______________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         allocate(EnHbUCh(NmBIN, 0:MeaTypeCnt)); EnHbUCh = 0.0_rp
         allocate(EnTotCh(NmBIN, 0:MeaTypeCnt)); EnTotCh = 0.0_rp
         allocate(HmtOvt2(NmBIN, 0:MeaTypeCnt)); HmtOvt2 = 0.0_rp
         allocate(HmtOvt3(NmBIN, 0:MeaTypeCnt)); HmtOvt3 = 0.0_rp
         allocate(HmtOvZm(NmBIN, 0:MeaTypeCnt)); HmtOvZm = 0.0_rp
         allocate(HmtOvbU(NmBIN, 0:MeaTypeCnt)); HmtOvbU = 0.0_rp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         allocate(HmOvSinusoidalPinSz(NmBin, 0:MeaTypeCnt)); HmOvSinusoidalPinSz = 0.0_rp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end if
!________________________________________________________________________________________  
!_________________ (3) Occupation, DouOcc, Some correlations ____________________________
!________________________________________________________________________________________   
      if(amyid == amstr) then      
         allocate(nOccpUp(NmBIN, 0:MeaTypeCnt)); nOccpUp = 0.0_rp 
         allocate(nOccpDw(NmBIN, 0:MeaTypeCnt)); nOccpDw = 0.0_rp
         allocate(nOccTot(NmBIN, 0:MeaTypeCnt)); nOccTot = 0.0_rp
         allocate(NmNeTot(NmBIN, 0:MeaTypeCnt)); NmNeTot = 0.0_rp
         allocate(RDouOcc(NmBIN, 0:MeaTypeCnt)); RDouOcc = 0.0_rp
         allocate(RNNDCrF(NmBIN, 0:MeaTypeCnt)); RNNDCrF = 0.0_rp
      end if
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!_________________ (4) Measure total electron density for fixed n_Occ case ______________
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(IfMuTqmcNt .or. IfFixnT) then
         VectorLength = 2*NSwepFix
         allocate(WeightList(VectorLength)); WeightList = + 0.000_rp
         allocate(n_Occ_List(VectorLength)); n_Occ_List = - 100.0_rp
      end if
!________________________________________________________________________________________  
!_________________ (5) The single-particle GrFs in r- and k-space _______________________
!________________________________________________________________________________________
      !!!!!!!!!! The r-space single-particle GrF
      allocate(RlSpGrnFtSwp(NumNC, NumNC, 2, 0:MeaTypeCnt)); RlSpGrnFtSwp = 0.0_rp
      allocate(RlSpGrnFtBIN(NumNC, NumNC, 2, 0:MeaTypeCnt)); RlSpGrnFtBIN = 0.0_rp
      if(amyid == amstr) then
         allocate(RlSpGrnFtAll(NmBIN, (NumNC+1)*NumNC/2, 2, 0:MeaTypeCnt)); RlSpGrnFtAll = 0.0_rp
      end if
      !!!!!!!!!! The k-space single-particle GrF
      if(IfFftEnPar) then
         allocate(KSpGreenFSwp(NumNC, NumNC, 2, 0:MeaTypeCnt)); KSpGreenFSwp = rp_Zzero
         allocate(KSpGreenFBIN(NumNC, NumNC, 2, 0:MeaTypeCnt)); KSpGreenFBIN = rp_Zzero
         if(amyid == amstr) then
            allocate(KSpGreenFAll(NmBIN, (NumNC+1)*NumNC/2, 4, 0:MeaTypeCnt)); KSpGreenFAll = 0.0_rp
         end if
      end if
!________________________________________________________________________________________  
!_________________ (6) The local spin and charge densities ______________________________
!________________________________________________________________________________________
      !!!!!!!!!! Accumulate local densities for an single BIN
      allocate(RSpLclOrdBIN(NumNS, 2, 0:MeaTypeCnt)); RSpLclOrdBIN = 0.0_rp
      !!!!!!!!!! Local densities and local spin moments Sz for all BINS
      if(amyid == amstr) then
         allocate(RSpLclOrdAll(NmBIN, NumNS, 5, 0:MeaTypeCnt)); RSpLclOrdAll = 0.0_rp
         allocate(TotMagMmtAll(NmBIN,           0:MeaTypeCnt)); TotMagMmtAll = 0.0_rp
      end if
!**************************************************************************************************    
!___________________ 2. Compute n(k) and pairing wavefunction from k-space GrF matrix _____________
!**************************************************************************************************
      if(IfFftEnPar) then
!________________________________________________________________________________________
!_________________ (0) For the momentum distribution n(k) _______________________________
!________________________________________________________________________________________
         !!!!!!!!!! The momentum distribution n(k) for an single BIN
         allocate(NkSgleSwp(NumNC, 2, 0:MeaTypeCnt)); NkSgleSwp = rp_Zzero
         allocate(NkSgleBIN(NumNC, 2, 0:MeaTypeCnt)); NkSgleBIN = rp_Zzero
         !!!!!!!!!! The momentum distribution n(k) for all BINS
         if(amyid == amstr) then
            allocate(NkDistrib(NmBIN, NumNC, 2, 0:MeaTypeCnt)); NkDistrib = 0.0_rp
         end if
!________________________________________________________________________________________
!_________________ (1) For the pairing wavefunction and condensate fractions ____________
!________________________________________________________________________________________
         !!!!!!!!!! The pairing matrix for an single BIN
         allocate(PairMtSwp(NumNC, NumNC, 0:MeaTypeCnt)); PairMtSwp = rp_Zzero
         allocate(PairMtBIN(NumNC, NumNC, 0:MeaTypeCnt)); PairMtBIN = rp_Zzero
         allocate(PairMtAcm(NumNC, NumNC, 0:MeaTypeCnt)); PairMtAcm = rp_Zzero
         !!!!!!!!!! The pairing matrix for all BINS
         if(amyid == amstr) then
            allocate(PairWvfct(NmBIN+1, NumNC, 6, 0:MeaTypeCnt)); PairWvfct = 0.0_rp
            allocate(Condfract(NmBIN+1,        2, 0:MeaTypeCnt)); Condfract = 0.0_rp
         end if
      end if
!**************************************************************************************************    
!___________________ 3. Static Correlation functions for both PBC and OBC _________________________
!**************************************************************************************************
      if(abs(PinSz) < rp_Eps) then
!________________________________________________________________________________________  
!_________________ (0) Accumulate r-space correlations for an single BIN ________________
!________________________________________________________________________________________
         !!!!!!!!!! Number of independent site pairs for correlations
         NmSitePair = merge(NumNC, (NumNC+1)*NumNC/2, IfPyObsPBC)
         !!!!!!!!!! The r-space correlation functions in various channels
         allocate(RealSpCrFSwp(NmSitePair, 40, 0:MeaTypeCnt)); RealSpCrFSwp = 0.0_rp
         allocate(RealSpCrFBIN(NmSitePair, 40, 0:MeaTypeCnt)); RealSpCrFBIN = 0.0_rp
!________________________________________________________________________________________  
!_________________ (1) R-space and k-space correlation functions for all BINs ___________
!________________________________________________________________________________________
         if(amyid == amstr) then
            !!!!!!!!!! R-space correlation functions for both OPEN and PERIODIC boundary conditions
            allocate(RealSpCrFAll(NmBIN, NmSitePair, 40, 0:MeaTypeCnt)); RealSpCrFAll = 0.0_rp
            !!!!!!!!!! k-space structure factors for only PERIODIC boundary conditions
            if(IfPyObsPBC) then
               allocate(KSpaceCrFAll(NmBIN, NumNC, 30, 0:MeaTypeCnt)); KSpaceCrFAll = rp_Zzero
            end if
         end if
!________________________________________________________________________________________  
!_________________ (2) Structure factors and Correlation ratios at specific _____________
!_____________________ k points for all BINs, for IfPyObsPBC == T case __________________
!________________________________________________________________________________________
         if( (amyid == amstr) .and. (IfPyObsPBC) ) then
            !!!!!!!!!! Structure factors at specific k points for all BINs
            allocate(KStructFactGamm(NmBIN, 20, 0:MeaTypeCnt)); KStructFactGamm = 0.0_rp
            allocate(KStructFactXPnt(NmBIN, 20, 0:MeaTypeCnt)); KStructFactXPnt = 0.0_rp
            allocate(KStructFactPntM(NmBIN, 20, 0:MeaTypeCnt)); KStructFactPntM = 0.0_rp
            !!!!!!!!!! Correlation ratios at specific k points for all BINs
            allocate(CorlatRatioGamm(NmBIN, 50, 0:MeaTypeCnt)); CorlatRatioGamm = 0.0_rp
            allocate(CorlatRatioXPnt(NmBIN, 50, 0:MeaTypeCnt)); CorlatRatioXPnt = 0.0_rp
            allocate(CorlatRatioPntM(NmBIN, 50, 0:MeaTypeCnt)); CorlatRatioPntM = 0.0_rp
         end if
!________________________________________________________________________________________  
!_________________ (3) The momentum distribution n(k) for IfPyObsPBC == T case __________
!________________________________________________________________________________________
         if( (amyid == amstr) .and. (IfPyObsPBC) .and. (.not. IfFftEnPar) ) then
            allocate(NkDistrib(NmBIN, NumNC, 2, 0:MeaTypeCnt)); NkDistrib = 0.0_rp
         end if
      end if
!**************************************************************************************************    
!___________________ 4. Dynamic Correlations for PERIODIC and OPEN boundary conditions ____________
!**************************************************************************************************
      if(IfTAU) then
!________________________________________________________________________________________  
!_________________ (0) Number of imaginary time points in measurements __________________
!________________________________________________________________________________________
         call GntTauGrid_Dyn()
!________________________________________________________________________________________  
!_________________ (1) Accumulate dynamic correlation results for an single BIN _________
!________________________________________________________________________________________
         !!!!!!!!!! Number of bond pairs that we compute the correlations
         NmSitePairTau = merge(NumNC, 3, IfPyObsPBC)
         !!!!!!!!!! The pair index integers for IfPyObsPBC==F case
         if(.not. IfPyObsPBC) then
            allocate(TauBondPairId(NmSitePairTau, 2))
            TauBondPairId(1, 1) = 1; TauBondPairId(1, 2) = 1
            TauBondPairId(2, 1) = 2; TauBondPairId(2, 2) = 2
            TauBondPairId(3, 1) = 1; TauBondPairId(3, 2) = 2
         end if
         !!!!!!!!!! Array for accumulating the dynamic correlations
         allocate(RealSpCrFTauSwp(0:NumTauPnt, NmSitePairTau, 40)); RealSpCrFTauSwp = 0.0_rp
         allocate(RealSpCrFTauBIN(0:NumTauPnt, NmSitePairTau, 40)); RealSpCrFTauBIN = 0.0_rp
         allocate(KSpCrFtTauforFT(0:NumTauPnt, NmSitePairTau, 34)); KSpCrFtTauforFT = rp_Zzero
!________________________________________________________________________________________  
!_________________ (2) Accumulate dynamic correlation results for all BINs ______________
!________________________________________________________________________________________
         if(amyid == amstr) then
            !!!!!!!!!! For IfPyObsPBC==T case, specific k points; For IfPyObsPBC==F case, the r-space correlations
            allocate(RorKspCrFTauAll(NmBIN, 0:NumTauPnt, merge(3, NmSitePairTau, IfPyObsPBC), 30))
            RorKspCrFTauAll = 0.0_rp
            !!!!!!!!!! Some specific results only for IfPyObsPBC==T case
            if(IfPyObsPBC) then
               !!!!!!!! Store the results at r-space local correlations for all BINS
               allocate(RSpCrFLclTauAll(NmBIN, 0:NumTauPnt, 10)); RSpCrFLclTauAll = 0.0_rp
               !!!!!!!! Store the results of \beta * [ S_{13}(\beta/2) + S_{12}(\beta/2) ]
               if(NmTDM >= LTrot/2) then
                  allocate(Spn13pS12HalfBt(NmBIN)); Spn13pS12HalfBt = 0.0_rp
               end if
               !!!!!!!! Store the results of G(k,beta/2) used to compute Fermi surface for all BINS
               if(NmTDM >= LTrot/2) then
                  allocate(GkTauBetaOv2All(NmBIN, NumNC, 3)); GkTauBetaOv2All = 0.0_rp
               end if
               !!!!!!!! Store the results of Current(k=0,beta/2) as the DC conductivity
               if(NmTDM >= LTrot/2) then
                  allocate(Dc_Conductivity(NmBIN)); Dc_Conductivity = 0.0_rp
               end if
            end if
         end if
!&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#
!_________________ (3) Dynamic correlation function in iw_n axis ________________________
!&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#
         if( (amyid == amstr) .and. (IfPyObsPBC) .and. (NmTDMType == 0 .or. NmTDMType == 1) ) then
            !!!!!!!!!! Number of iw_n to be computed
            NmFrqFermi = 2
            NmFrqBoson = 2
            !!!!!!!!!! Fermionic channel
            allocate(NonIntG0kIwn(       NumNC, 4*NmFrqFermi)); NonIntG0kIwn = 0.0_rp
            allocate(FermiGrF_Iwn(NmBIN, NumNC, 4*NmFrqFermi)); FermiGrF_Iwn = 0.0_rp
            allocate(SelfEnrg_Iwn(NmBIN, NumNC, 4*NmFrqFermi)); SelfEnrg_Iwn = 0.0_rp
            !!!!!!!!!! Bosonic channels
            allocate(BosonCrF_Iwn(NmBIN, NumNC, NmFrqBoson, 30)); BosonCrF_Iwn = 0.0_rp
            !!!!!!!!!! The Quasi-particle weights
            allocate(QuasiParWhgt(NmBIN, NumNC)); QuasiParWhgt = 0.0_rp
            !!!!!!!!!! Drude conductivity and superfluid density
            allocate(DrudeSupfdWt(NmBIN, 6)); DrudeSupfdWt = 0.0_rp
         end if
      end if
      
   end subroutine InitPhyObsev
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine GntTauGrid_Dyn()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GntTauGrid_Dyn()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to generate the tau grid for the dynamic measurements. The tau grid can be 
!                 equal distance or with nonequal distance of tau.
! KEYWORDS: Generate Tau grid for dynamic measurements.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-10
! DESCRIPTION: Generate Tau grid for dynamic measurements.
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use Observable
      implicit none
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________      
      integer I0, NT, NT_Last, TauDiscrt, NmTauCount
      integer TauBgn, TauEnd, TauLeft
      integer TauPntVal_Tmp(0:NmTDM)
!______________________________________________________________________________________________________________     
!____________________________ Main calculations of Generating the tau grid ____________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Some initializations before the discretization of Tau _____________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Determine the length of tau to be discretized ____________________
!________________________________________________________________________________________
      if(NmTDM < LTrot/2) then
         TauDiscrt = NmTDM
      else
         TauDiscrt = LTrot/2
      end if
!________________________________________________________________________________________      
!_________________ (1) Initialization of NumTauPnt and TauPntVal_Tmp ____________________
!________________________________________________________________________________________
      NumTauPnt = 0
      TauPntVal_Tmp = 0
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!___________________ 1. For IfEqDistDt == T case __________________________________________________
!______________________ The tau grid is with equal distance of tau ________________________________
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(IfEqDistDt) then
!________________________________________________________________________________________      
!_________________ (0) The number of tau points in measurements _________________________
!________________________________________________________________________________________
         NumTauPnt = TauDiscrt / NvMeaDyn
         if(mod(TauDiscrt, NvMeaDyn) /= 0) then
            NumTauPnt = TauDiscrt / NvMeaDyn + 1
         end if
!________________________________________________________________________________________      
!_________________ (1) Record values of tau for all the tau points ______________________
!________________________________________________________________________________________
         do I0 = 0, NumTauPnt-1, +1
            TauPntVal_Tmp(I0) = I0 * NvMeaDyn
         enddo
         TauPntVal_Tmp(NumTauPnt) = TauDiscrt
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!___________________ 2. For IfEqDistDt == F case __________________________________________________
!______________________ The tau grid is with nonequal distance of tau _____________________________
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      else
!________________________________________________________________________________________      
!_________________ (0) The Tau == 0 point as the zero-th point __________________________
!________________________________________________________________________________________
         NumTauPnt = 0; TauPntVal_Tmp(NumTauPnt) = 0
         TauLeft = 0
!________________________________________________________________________________________      
!_________________ (1) In [(I0-1)*TauDiscrt/5, I0*TauDiscrt/5], NvMeaDyn = I0 ___________
!________________________________________________________________________________________
         do I0 = 1, 5, +1
            NvMeaDyn = I0
            TauBgn = (I0-1) * TauDiscrt/5 - TauLeft + NvMeaDyn
            TauEnd =  I0    * TauDiscrt/5
            do NT = TauBgn, TauEnd, +1
               if( mod(NT-TauBgn, NvMeaDyn) == 0 ) then
                  NumTauPnt = NumTauPnt + 1
                  TauPntVal_Tmp(NumTauPnt) = NT
               end if
            enddo
            TauLeft = TauEnd - TauPntVal_Tmp(NumTauPnt)
         enddo
!________________________________________________________________________________________      
!_________________ (2) Take the Tau == TauDiscrt point as the last point ________________
!________________________________________________________________________________________
         if( TauPntVal_Tmp(NumTauPnt) < TauDiscrt ) then
            NumTauPnt = NumTauPnt + 1
            TauPntVal_Tmp(NumTauPnt) = TauDiscrt
         end if
      end if
!**************************************************************************************************     
!___________________ 3. Obtain final NumTauPnt integer and TauPntVal vector _______________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Record the integer index for Tau == BetaT/2 ______________________
!________________________________________________________________________________________  
      if(NmTDM >= LTrot/2) then
         TauHalfBetaT = NumTauPnt
      end if
!________________________________________________________________________________________      
!_________________ (1) For NmTDM >= LTrot/2, obtain all tau point in [0, \beta] _________
!________________________________________________________________________________________
      if(NmTDM >= LTrot/2) then
         allocate(TauPntVal0B(0:2*NumTauPnt))
         TauPntVal0B(0:NumTauPnt) = TauPntVal_Tmp(0:NumTauPnt)
         do I0 = NumTauPnt-1, 0, -1
            NT = LTrot - TauPntVal0B(I0)
            TauPntVal0B(2*NumTauPnt-I0) = NT
         enddo
      end if
!________________________________________________________________________________________      
!_________________ (2) For the NmTDM > LTrot/2 case --> Mirror Symmetry _________________
!________________________________________________________________________________________      
      if(NmTDM > LTrot/2) then
         NmTauCount = 0
         do I0 = NumTauPnt-1, 0, -1
            NT = LTrot - TauPntVal_Tmp(I0)
            if(NT <= NmTDM) then
               NmTauCount = NmTauCount + 1
               TauPntVal_Tmp(NumTauPnt+NmTauCount) = NT
               NT_Last = NT
            else
               exit
            end if
         enddo
         if(NT_Last < NmTDM) then
            NmTauCount = NmTauCount + 1
            TauPntVal_Tmp(NumTauPnt+NmTauCount) = NmTDM
         end if
         NumTauPnt = NumTauPnt + NmTauCount
      end if
      allocate(TauPntVal(0:NumTauPnt))
      TauPntVal(0:NumTauPnt) = TauPntVal_Tmp(0:NumTauPnt)

   end subroutine GntTauGrid_Dyn
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine CPMCObservInit()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCObservInit()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the initialization for observables allocated in subroutine 
!                  InitAllocObs, just give zero values for the matrices.
! KEYWORDS: Initialization of matrices in Observable module.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Initialization of constants in CPMC, including:
!             (0) Initialization of some measurement parameters;
!             (1) Allocate the observable matrices and initializations.
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use MPISetting
      use Observable
      implicit none
!______________________________________________________________________________________________________________     
!____________________________ Main calculations of Observables Initialization _________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Initializations for all quantities related to measurements ________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) For summation of weight average of all sweeps ____________________
!________________________________________________________________________________________
      WtMeanSumBIN = 0.0_rp
!________________________________________________________________________________________      
!_________________ (1) For Static physical observables __________________________________
!________________________________________________________________________________________
      !!!!!!!!!! Number of measurements
      NObsStat = 0
      !!!!!!!!!! The energies, fillings, density correlations and single-particle GrFs
      EngOccCrFBIN = 0.0_rp
      RlSpGrnFtBIN = 0.0_rp
      if(IfFftEnPar) KSpGreenFBIN = rp_Zzero
      !!!!!!!!!! The spin and charge densities
      RSpLclOrdBIN = 0.0_rp
      !!!!!!!!!! n(k) and the paring matrices
      if(IfFftEnPar) then
         NkSgleBIN = rp_Zzero; PairMtBIN = rp_Zzero
      end if
      !!!!!!!!!! The r-space correlation functions depending on OBC or PBC
      if(abs(PinSz) < rp_Eps) RealSpCrFBIN = rp_Zzero
!________________________________________________________________________________________      
!_________________ (2) Dynamic observables ______________________________________________
!________________________________________________________________________________________
      if(IfTAU) then
         NObsDynm = 0
         RealSpCrFTauBIN = rp_Zzero
      end if
      
   end subroutine CPMCObservInit
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$