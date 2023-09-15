!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL SPIN PINNING FIELDS; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A few subroutines for performing the CPMC finialization after the whole CPMC simulation, including:
!                  (1) Finalization of vectors and matrices in CoreParamt module;
!                  (2) Finalization of vectors and matrices in Observable module.
! COMMENT: CPMC Finalization process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   CPMCFinalize   --> Subroutine to call subroutines to perform the total CPMC finalization;
!   
!   CoreParamtFinl --> Subroutine to perform the finalization for CoreParamt module;
!   ObservableFinl --> Subroutine to perform the finalization for Observable module.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine CPMCFinalize()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCFinalize()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the CPMC Finalization after the whole CPMC simulation.
! KEYWORDS: CPMC initialization.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: CPMC Finalization process.
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use CoreParamt
      use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of CPMC Finalization _______________________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCFinalize process _________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "(A)") "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         write(*, "(2x, 'CPMCFinalize: Finalization for the QMC simulation!')")
      end if
!**************************************************************************************************  
!___________________ 0. Finalization of CoreParamt module _________________________________________
!**************************************************************************************************
		call CoreParamtFinl()  
!**************************************************************************************************  
!___________________ 1. Finalization of Observable module _________________________________________
!**************************************************************************************************
		call ObservableFinl()	
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCFinalize process _________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(      
      if(amyid == amstr) then
         write(*, "(A)") "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      end if
		
   end subroutine CPMCFinalize
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine CoreParamtFinl()  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CoreParamtFinl()  
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the finalization for the CoreParamt module, mainly deallocate arrays.
! KEYWORDS: Finalization of module.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Finalization for module CoreParamt. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!______________________________ Deallocate the allocated Arrays _______________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************  
!___________________ 0. Define the finite lattice size ____________________________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (2) Nearest-Neighbor unit cells of all unit cells ____________________
!________________________________________________________________________________________ 
		if(allocated(   KpList)) deallocate(   KpList)
		if(allocated(InvKpList)) deallocate(InvKpList)
		if(allocated(   StList)) deallocate(   StList)
		if(allocated(InvStList)) deallocate(InvStList)
!________________________________________________________________________________________  
!_________________ (3) Distance indexes for two unit cells in the finite lattice ________
!________________________________________________________________________________________      
      if(allocated(IminusJ)) deallocate(IminusJ)
!________________________________________________________________________________________  
!_________________ (4) NN and 2NN lattice sites information _____________________________
!________________________________________________________________________________________      
      if(allocated(FNNBond)) deallocate(FNNBond)
      if(allocated(SNNBond)) deallocate(SNNBond)
      if(allocated(TNNBond)) deallocate(TNNBond)
!________________________________________________________________________________________  
!_________________ (5) Open or periodic boundary conditions in a1 and a2 directions _____
!________________________________________________________________________________________      
      if(allocated(FNNStBnd)) deallocate(FNNStBnd) 
      if(allocated(SNNStBnd)) deallocate(SNNStBnd) 
      if(allocated(TNNStBnd)) deallocate(TNNStBnd)
!**************************************************************************************************  
!___________________ 2. Exp(-\Delta\tau\hat{V}) and Exp(-\Delta\tau\hat{K}) _______________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (1) For the kinetic propagator of H_0 part in many-body H ____________
!________________________________________________________________________________________
!____________________________________________________________________________	  
!________________ [0] Eigenvalues and eigenvectors for H_0 Hamiltonian ______
!____________________________________________________________________________
      if(allocated(H0_EgValue)) deallocate(H0_EgValue)
      if(allocated(H0_EgVctor)) deallocate(H0_EgVctor)
!____________________________________________________________________________	  
!________________ [1] Use FFT method to calculate ___________________________
!____________________________________________________________________________
      if(allocated(ExpdtEofH0)) deallocate(ExpdtEofH0)
!____________________________________________________________________________	  
!________________ [2] Use matrix product of full matrices ___________________
!____________________________________________________________________________
      if(allocated(ExpdtH0Mat)) deallocate(ExpdtH0Mat)
!________________________________________________________________________________________  
!_________________ (2) Exp(-\Delta\tau\hat{V}) related quantities _______________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [2] Auxiliary fields for HubbU interaction ________________
!____________________________________________________________________________       
      if(allocated(IsingbU)) deallocate(IsingbU)
!____________________________________________________________________________	  
!________________ [3] Two-component auxiliary fields and ____________________
!____________________ Corresonding \Delta_{ii} and phase ____________________
!____________________________________________________________________________ 
      if(allocated(ExpobU       )) deallocate(ExpobU       )
      if(allocated(ExpobUInv    )) deallocate(ExpobUInv    )
      if(allocated(ExpobU_H0T   )) deallocate(ExpobU_H0T   )
      if(allocated(ExpobUInv_H0T)) deallocate(ExpobUInv_H0T)
      if(allocated(DeltbU_H0T   )) deallocate(DeltbU_H0T   )
!________________________________________________________________________________________  
!_________________ (3) The exponential kinetic matrices for H_0 in B_x and B_T __________
!________________________________________________________________________________________   
!____________________________________________________________________________	  
!________________ [0] Eigenvalues, eigenvectors for H_0 without Mu term _____
!____________________________________________________________________________  
      if(allocated(HT_EigValu)) deallocate(HT_EigValu)
      if(allocated(HT_EigVect)) deallocate(HT_EigVect)
!____________________________________________________________________________	  
!________________ [1] Use FFT method to calculate ___________________________
!____________________________________________________________________________ 
      if(allocated(ExpdtOfHTe)) deallocate(ExpdtOfHTe)
      if(allocated(ExpdtTryHT)) deallocate(ExpdtTryHT)
      if(FFTEXPDTH0 .or. FFTEXPDTHT) then
         call FFTSettingEnd()
      end if
!____________________________________________________________________________	  
!________________ [3] For the UHF kind of trial Hamiltonian H_T _____________
!____________________________________________________________________________
      if(allocated(UHF_Const)) deallocate(UHF_Const)
      if(allocated(MagMoment)) deallocate(MagMoment)
!**************************************************************************************************  
!___________________ 3. The Core matrices used in CPMC simulations ________________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) UDV matrices at left and right sides for propagation _____________
!________________________________________________________________________________________
      if(allocated(ULeft   )) deallocate(ULeft   )
      if(allocated(DLeftVec)) deallocate(DLeftVec)
      if(allocated(VLeft   )) deallocate(VLeft   )
      if(allocated(URght   )) deallocate(URght   )
      if(allocated(DRghtVec)) deallocate(DRghtVec)
      if(allocated(VRght   )) deallocate(VRght   )
      if(allocated(LogScaleRght)) deallocate(LogScaleRght)
      if(allocated(LogScaleLeft)) deallocate(LogScaleLeft)
      if(allocated(LogScaleOfHT)) deallocate(LogScaleOfHT)
!________________________________________________________________________________________  
!_________________ (1) The Green's function matrix for propagation ______________________
!________________________________________________________________________________________
      if(allocated(GrnFunct)) deallocate(GrnFunct)
!________________________________________________________________________________________  
!_________________ (2) Weight and their logs, and all weights vector ____________________
!________________________________________________________________________________________    
      if(allocated(IdptWkIndx)) deallocate(IdptWkIndx)  
      if(allocated(WghtProc)) deallocate(WghtProc)
      if(allocated(Log_Wght)) deallocate(Log_Wght)
      if(allocated(WghtTotl)) deallocate(WghtTotl)
!________________________________________________________________________________________  
!_________________ (3) The growth estimator and the ancestry link _______________________
!________________________________________________________________________________________
      if(allocated(GrowthCoefft)) deallocate(GrowthCoefft)
!________________________________________________________________________________________  
!_________________ (4) The ancestry link and number of ancestry walkers _________________
!________________________________________________________________________________________
      if(allocated(AncestryLink)) deallocate(AncestryLink)
!________________________________________________________________________________________	  
!_________________ (5) Store the UDV matrices for measurement in [BetaT, 0] sweep _______
!________________________________________________________________________________________
      if(allocated(UStMt   )) deallocate(UStMt   )
      if(allocated(DVecStMt)) deallocate(DVecStMt)
      if(allocated(DScalLog)) deallocate(DScalLog)
      if(allocated(VStMt   )) deallocate(VStMt   )
!________________________________________________________________________________________  
!_________________ (6) Temporary Green's Function matrices for Static and Dynamic _______
!________________________________________________________________________________________  
      if(allocated(GrnFRTmp00)) deallocate(GrnFRTmp00)
      if(allocated(GrnFRTmp11)) deallocate(GrnFRTmp11)
      if(allocated(GrnFRTmp22)) deallocate(GrnFRTmp22)
      if(allocated(GrFKCTmp00)) deallocate(GrFKCTmp00)
      if(allocated(GrFKCTmp11)) deallocate(GrFKCTmp11)
!________________________________________________________________________________________  
!_________________ (7) Temporary Green's Function matrices for Dynamic __________________
!________________________________________________________________________________________ 
      if(allocated(GrF00)) deallocate(GrF00)
      if(allocated(GrFT0)) deallocate(GrFT0)
      if(allocated(GrF0T)) deallocate(GrF0T)
      if(allocated(GrFTT)) deallocate(GrFTT)
!________________________________________________________________________________________  
!_________________ (8) Identity matrix, used for initializations ________________________
!________________________________________________________________________________________        
      if(allocated(IdMtR)) deallocate(IdMtR)
!**************************************************************************************************	  
!___________________ 4. Some other quantities used in CPMC simulations ____________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (3) Whether to save the auxiliary-field configurations _______________
!________________________________________________________________________________________
      if(allocated(RdWrtField)) deallocate(RdWrtField)
      if(allocated(RdWrtIntgr)) deallocate(RdWrtIntgr)
      
	end subroutine CoreParamtFinl
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine ObservableFinl()  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ObservableFinl()  
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the finalization for the Observable module, mainly deallocate arrays.
! KEYWORDS: Finalization of Observable module.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Finalization for module Observables. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		use Observable
            use MPISetting
            use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!____________________________ Deallocate allocated arrays in the simulations __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************  
!___________________ 1. Measure energies, occupations and local densities _________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Accumulate Energies and Occupations for every single BIN _________
!________________________________________________________________________________________
      if(allocated(EngOccCrFSwp)) deallocate(EngOccCrFSwp)
      if(allocated(EngOccCrFBIN)) deallocate(EngOccCrFBIN)
!________________________________________________________________________________________  
!_________________ (1) Expectation energies for every term in the model _________________
!________________________________________________________________________________________
      if(allocated(EnHopt1)) deallocate(EnHopt1)
      if(allocated(EnHopt2)) deallocate(EnHopt2)
      if(allocated(EnHopt3)) deallocate(EnHopt3)
      if(allocated(EnZmFld)) deallocate(EnZmFld)
      if(allocated(EnPinSz)) deallocate(EnPinSz)
      if(allocated(EnDopCh)) deallocate(EnDopCh)
      if(allocated(EnHubbU)) deallocate(EnHubbU)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(allocated(EnSinusoidalPinSz)) deallocate(EnSinusoidalPinSz)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(allocated(EnTotal)) deallocate(EnTotal)
      if(allocated(EpOvEfg)) deallocate(EpOvEfg)
!________________________________________________________________________________________  
!_________________ (2) First-order derivative of model parameters _______________________
!________________________________________________________________________________________
      if(allocated(EnHbUCh)) deallocate(EnHbUCh)
      if(allocated(EnTotCh)) deallocate(EnTotCh)
      if(allocated(HmtOvt2)) deallocate(HmtOvt2)
      if(allocated(HmtOvt3)) deallocate(HmtOvt3)
      if(allocated(HmtOvZm)) deallocate(HmtOvZm)
      if(allocated(HmtOvbU)) deallocate(HmtOvbU)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (allocated(HmOvSinusoidalPinSz)) deallocate(HmOvSinusoidalPinSz)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!________________________________________________________________________________________  
!_________________ (3) Occupation, DouOcc, Some correlations ____________________________
!________________________________________________________________________________________   
      if(allocated(nOccpUp)) deallocate(nOccpUp)
      if(allocated(nOccpDw)) deallocate(nOccpDw)
      if(allocated(nOccTot)) deallocate(nOccTot)
      if(allocated(NmNeTot)) deallocate(NmNeTot)
      if(allocated(RDouOcc)) deallocate(RDouOcc)
      if(allocated(RNNDCrF)) deallocate(RNNDCrF)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!_________________ (4) Measure total electron density for fixed n_Occ case ______________
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(allocated(WeightList)) deallocate(WeightList)
      if(allocated(n_Occ_List)) deallocate(n_Occ_List)
!________________________________________________________________________________________  
!_________________ (5) The single-particle GrFs in r- and k-space _______________________
!________________________________________________________________________________________
      if(allocated(RlSpGrnFtSwp)) deallocate(RlSpGrnFtSwp)
      if(allocated(RlSpGrnFtBIN)) deallocate(RlSpGrnFtBIN)
      if(allocated(RlSpGrnFtAll)) deallocate(RlSpGrnFtAll)
      if(allocated(KSpGreenFSwp)) deallocate(KSpGreenFSwp)
      if(allocated(KSpGreenFBIN)) deallocate(KSpGreenFBIN)
      if(allocated(KSpGreenFAll)) deallocate(KSpGreenFAll)
!________________________________________________________________________________________  
!_________________ (6) The local spin and charge densities ______________________________
!________________________________________________________________________________________
      if(allocated(RSpLclOrdBIN)) deallocate(RSpLclOrdBIN)
      if(allocated(RSpLclOrdAll)) deallocate(RSpLclOrdAll)
      if(allocated(TotMagMmtAll)) deallocate(TotMagMmtAll)
!**************************************************************************************************	 
!___________________ 2. Compute n(k) and pairing wavefunction from k-space GrF matrix _____________
!**************************************************************************************************
!________________________________________________________________________________________
!_________________ (0) For the momentum distribution n(k) _______________________________
!________________________________________________________________________________________
      if(allocated(NkSgleSwp)) deallocate(NkSgleSwp)
      if(allocated(NkSgleBIN)) deallocate(NkSgleBIN)
      if(allocated(NkDistrib)) deallocate(NkDistrib)
!________________________________________________________________________________________
!_________________ (1) For the pairing wavefunction and condensate fractions ____________
!________________________________________________________________________________________
      if(allocated(PairMtSwp)) deallocate(PairMtSwp)
      if(allocated(PairMtBIN)) deallocate(PairMtBIN)
      if(allocated(PairMtAcm)) deallocate(PairMtAcm)
      if(allocated(PairWvfct)) deallocate(PairWvfct)
      if(allocated(Condfract)) deallocate(Condfract)
!**************************************************************************************************	 
!___________________ 3. Static Correlation functions for both PBC and OBC _________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Accumulate r-space correlations for an single BIN ________________
!________________________________________________________________________________________
      if(allocated(RealSpCrFSwp)) deallocate(RealSpCrFSwp)
      if(allocated(RealSpCrFBIN)) deallocate(RealSpCrFBIN)
!________________________________________________________________________________________  
!_________________ (1) R-space and k-space correlation functions for all BINs ___________
!________________________________________________________________________________________
      if(allocated(RealSpCrFAll)) deallocate(RealSpCrFAll)
      if(allocated(KSpaceCrFAll)) deallocate(KSpaceCrFAll)
!________________________________________________________________________________________  
!_________________ (2) Structure factors and Correlation ratios at specific _____________
!_____________________ k points for all BINs, for IfPyObsPBC == T case __________________
!________________________________________________________________________________________
      if(allocated(KStructFactGamm)) deallocate(KStructFactGamm)
      if(allocated(KStructFactXPnt)) deallocate(KStructFactXPnt)
      if(allocated(KStructFactPntM)) deallocate(KStructFactPntM)
      if(allocated(CorlatRatioGamm)) deallocate(CorlatRatioGamm)
      if(allocated(CorlatRatioXPnt)) deallocate(CorlatRatioXPnt)
      if(allocated(CorlatRatioPntM)) deallocate(CorlatRatioPntM)
!**************************************************************************************************	 
!___________________ 4. Dynamic Correlations for PERIODIC and OPEN boundary conditions ____________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Number of imaginary time points in measurements __________________
!________________________________________________________________________________________  
      if(allocated(TauPntVal  )) deallocate(TauPntVal  ) 
      if(allocated(TauPntVal0B)) deallocate(TauPntVal0B)
!________________________________________________________________________________________  
!_________________ (1) Accumulate r-space dynamic correlation for single BIN ____________
!_____________________ For both PERIODIC and OPEN boundary conditions ___________________
!________________________________________________________________________________________
      if(allocated(TauBondPairId  )) deallocate(TauBondPairId  )
      if(allocated(RealSpCrFTauSwp)) deallocate(RealSpCrFTauSwp)
      if(allocated(RealSpCrFTauBIN)) deallocate(RealSpCrFTauBIN)
      if(allocated(KSpCrFtTauforFT)) deallocate(KSpCrFtTauforFT)
!________________________________________________________________________________________  
!_________________ (2) Local dynamic correlations and G(k,beta/2) for all BINS __________
!_____________________ Only for PERIODIC boundary conditions ____________________________
!________________________________________________________________________________________
      if(allocated(RorKspCrFTauAll)) deallocate(RorKspCrFTauAll)
      if(allocated(RSpCrFLclTauAll)) deallocate(RSpCrFLclTauAll)
      if(allocated(Spn13pS12HalfBt)) deallocate(Spn13pS12HalfBt)
      if(allocated(GkTauBetaOv2All)) deallocate(GkTauBetaOv2All)
      if(allocated(Dc_Conductivity)) deallocate(Dc_Conductivity)
!&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#
!_________________ (3) Dynamic correlation function in iw_n axis ________________________
!&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#
      if(allocated(NonIntG0kIwn)) deallocate(NonIntG0kIwn)
      if(allocated(FermiGrF_Iwn)) deallocate(FermiGrF_Iwn)
      if(allocated(SelfEnrg_Iwn)) deallocate(SelfEnrg_Iwn)
      if(allocated(BosonCrF_Iwn)) deallocate(BosonCrF_Iwn)
      if(allocated(QuasiParWhgt)) deallocate(QuasiParWhgt)
      if(allocated(DrudeSupfdWt)) deallocate(DrudeSupfdWt)
		
   end subroutine ObservableFinl
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$