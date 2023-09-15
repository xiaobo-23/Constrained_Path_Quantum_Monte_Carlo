!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL SPIN PINNING FIELDS; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A module and a few subroutines used for defining the observables which needs to be calculated in the CPMC 
!            Calculations. These observables includes the equal-time quantities (static) and time-displayed 
!            quantities (dynamic). Including: 
!             (0) Integer quantities used in equal-time measurement;
!             (1) Integer quantities used in time-displayed measurement;
!             (2) Expectation energies for every term in model Hamiltonian and double occupancy;
!             (2) Equal-time Real space correlation functions for some chosen lattice sites pairs in every BIN simulation;
!             (3) Equal-time Real space correlation functions for all lattice sites in every BIN simulation;
!             (4) Time-displayed Real space correlation functions for all lattice sites in every BIN simulation;
!             (5) Average Expectation energies of all terms for all BIN simulations;
!             (6) Average Equal-time Real space correlation functions for some chosen lattice sites pairs for 
!                        all BIN simulations;
!             (7) Average Equal-time reciprocal space correlation functions for all lattice sites for all BIN simulations;
! COMMENT: Physical Observables.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!    Observable     --> Subroutine to define the related quantities for calculations;
!    ObservableInit --> Subroutine to perform the initializations of the quantities defined in this module.
!    ObservableFinl --> Subroutine to perform finalization of Observable module.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Module ______________________________________________________________________
!________________________________________________________________________________________________________________________
	module Observable
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________    
		use RealPrecsn
		implicit none
!**************************************************************************************************  
!___________________ 0. Some parameters used for measurements and outputs _________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Settings for the Measurements ____________________________________
!________________________________________________________________________________________
            logical IfPyObsPBC    ! Boundary conditions for measurements
            logical IfFftEnPar    ! Whether only to measure energies and pairing by k-space GrF
            logical IfSwepReWt    ! Reweighting for results of sweeps in a BIN
            logical IfCrfDfBkg    ! Whether to deduct the back ground in correlation functions
!________________________________________________________________________________________  
!_________________ (1) Accumulate walker weights for measurements _______________________
!________________________________________________________________________________________  
            real(rp) WeightSumSwp ! The sum of weights of present sweep
            real(rp) WeightAvgSwp ! The weight average of present sweep
            real(rp) WtMeanSumBIN ! The summation of weight average of all sweeps for a BIN
!________________________________________________________________________________________  
!_________________ (2) Settings for the static measurements _____________________________
!________________________________________________________________________________________  
            integer NObsStat      ! Number of static measurements at Tau=BetaT
            logical IfM2OneMea    ! Whether to measure in Beta --> 0 path
            integer MeaM2One      ! Number of measurements in single [BetaT, 0] sweep
!________________________________________________________________________________________  
!_________________ (3) Settings for the dynamic measurements ____________________________
!________________________________________________________________________________________
            !!!!!!!!!! Setting for the dynamic measurements
            logical IfTAU         ! Whether to calculate the time-displayed quantities
            logical IfTau0Rand    ! Whether to use to random NT as Tau==0 point for dynamic measurements 
            integer NmTDMType     ! ==0, NmTDM=LTrot/2; ==1, NmTDM=LTrot; ==2, use input NmTDM
            integer NmTDM         ! Total length of the time-displaced measurements (+1)
            logical IfEqDistDt    ! Whether to apply the equal interval \Delta\tau for S(k, tau)
            integer NvMeaDyn      ! Number of interval time slices for Dynamic measurements
            integer TauHalfBetaT  ! Integer index of Tau == BetaT/2, for  NmTDM >= LTrot/2 case
            !!!!!!!!!! Number of measurements 
            integer NObsDynm      ! Number of Dynamic measurements in single BIN
            !!!!!!!!!! Whether to calculate and output Dynamic correlation at all k points
            logical IfDyGrFOut, IfDySpnOut, IfDyDenOut
            logical IfDyPstOut, IfDyDWvOut, IfDyCurrnt
!________________________________________________________________________________________  
!_________________ (4) Settings for the output folders for results ______________________
!________________________________________________________________________________________
            character(20) FileOutAdd(0:1), FileAddTxt(0:1)
!**************************************************************************************************  
!___________________ 1. Measure energies, occupations and local densities _________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Accumulate Energies and Occupations for every single BIN _________
!________________________________________________________________________________________ 
      real(rp), allocatable :: EngOccCrFBIN(:, :) ! Energies and occupations for single BIN
      real(rp), allocatable :: EngOccCrFSwp(:, :) ! Single sweep --> Measurement at tau = BetaT and [B, 0] sweep
!________________________________________________________________________________________  
!_________________ (1) Expectation energies for every term in the model for all BINs ____
!________________________________________________________________________________________
      real(rp), allocatable :: EnHopt1(:, :)   ! Expectation energy for Hopt1 term
      real(rp), allocatable :: EnHopt2(:, :)   ! Expectation energy for Hopt2 term
      real(rp), allocatable :: EnHopt3(:, :)   ! Expectation energy for Hopt3 term
      real(rp), allocatable :: EnZmFld(:, :)   ! Expectation energy for Zeeman field term
      real(rp), allocatable :: EnPinSz(:, :)   ! Expectation energy for Pinning Sz term
      real(rp), allocatable :: EnDopCh(:, :)   ! Expectation energy for dopping chemical potential term
      real(rp), allocatable :: EnHubbU(:, :)   ! Expectation energy for HubbU term

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real(rp), allocatable :: EnSinusoidalPinSz(:, :) ! Expectation energy for sinusoidal spin pinning fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      real(rp), allocatable :: EnTotal(:, :)   ! Expectation energy for Total model Hamiltonian
      real(rp), allocatable :: EpOvEfg(:, :, :)! Expectation energy per particle in unit of E_{FG}
!________________________________________________________________________________________  
!_________________ (2) First-order derivative of model parameters for all BINs __________
!________________________________________________________________________________________
      real(rp), allocatable :: EnHbUCh(:, :)   ! Energy of HubbU term without chemical potential
      real(rp), allocatable :: EnTotCh(:, :)   ! Total energy without chemical potential
      real(rp), allocatable :: HmtOvt2(:, :)   ! First-order derivative of Hopt2
      real(rp), allocatable :: HmtOvt3(:, :)   ! First-order derivative of Hopt3
      real(rp), allocatable :: HmtOvZm(:, :)   ! First-order derivative of ZmFdz
      real(rp), allocatable :: HmtOvbU(:, :)   ! First-order derivative of HubbU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real(rp), allocatable :: HmOvSinusoidalPinSz(:, :)  ! First-order derivative of SinusoidalPinSz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!________________________________________________________________________________________  
!_________________ (3) Occupation, DouOcc, Some correlations for all BINs _______________
!________________________________________________________________________________________      
      real(rp), allocatable :: nOccpUp(:, :)   ! Expectation value for spin-up   occupation
      real(rp), allocatable :: nOccpDw(:, :)   ! Expectation value for spin-down occupation
      real(rp), allocatable :: nOccTot(:, :)   ! Expectation value for total occupation
      real(rp), allocatable :: NmNeTot(:, :)   ! Total number of electrons in the system
      real(rp), allocatable :: RDouOcc(:, :)   ! Expectation value for double occupancy
      real(rp), allocatable :: RNNDCrF(:, :)   ! Expectation value for n_{i\sigma}*n_{i+1,-\sigma}
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
!_________________ (4) Measure total electron density for fixed n_Occ case ______________
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      integer VectorLength
      real(rp), allocatable :: WeightList(:)
      real(rp), allocatable :: n_Occ_List(:)
!________________________________________________________________________________________  
!_________________ (5) The single-particle GrFs in r- and k-space _______________________
!________________________________________________________________________________________
      !!!!!!!!!! The r-space single-particle GrF
      real(rp), allocatable :: RlSpGrnFtSwp(:, :, :, :)     ! Single swap
      real(rp), allocatable :: RlSpGrnFtBIN(:, :, :, :)     ! Single BIN
      real(rp), allocatable :: RlSpGrnFtAll(:, :, :, :)
      !!!!!!!!!! The k-space single-particle GrF
      complex(rp), allocatable :: KSpGreenFSwp(:, :, :, :) ! For every single BIN
      complex(rp), allocatable :: KSpGreenFBIN(:, :, :, :) ! For every single BIN
      real(rp)   , allocatable :: KSpGreenFAll(:, :, :, :) ! For all BINs
!________________________________________________________________________________________  
!_________________ (6) The local spin and charge densities ______________________________
!________________________________________________________________________________________
      real(rp), allocatable :: RSpLclOrdBIN(:, :, :   )   ! Accumulate local densities for an single BIN
      real(rp), allocatable :: RSpLclOrdAll(:, :, :, :)   ! Local densities for all BINS
      real(rp), allocatable :: TotMagMmtAll(:, :      )   ! local spin moments Sz for all BINS
!**************************************************************************************************	 
!___________________ 2. Compute n(k) and pairing wavefunction from k-space GrF matrix _____________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) For the momentum distribution n(k) _______________________________
!________________________________________________________________________________________
      complex(rp), allocatable :: NkSgleSwp(:, :, :)    ! For every single sweep
      complex(rp), allocatable :: NkSgleBIN(:, :, :)    ! For every single BIN
      real(rp)   , allocatable :: NkDistrib(:, :, :, :) ! For all BINs
!________________________________________________________________________________________  
!_________________ (1) For the pairing wavefunction and condensate fractions ____________
!________________________________________________________________________________________
      complex(rp), allocatable :: PairMtSwp(:, :, :)    ! For every single sweep
      complex(rp), allocatable :: PairMtBIN(:, :, :)    ! For every single BIN
      complex(rp), allocatable :: PairMtAcm(:, :, :)    ! Accumulation for BINs
      real(rp)   , allocatable :: PairWvfct(:, :, :, :) ! For all BINs
      real(rp)   , allocatable :: Condfract(:,    :, :) ! For all BINs
!**************************************************************************************************	 
!___________________ 3. Static Correlation functions for both PBC and OBC _________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Accumulate r-space correlations for an single BIN ________________
!________________________________________________________________________________________
      !!!!!!!!!! Number of independent site pairs for correlations
      integer NmSitePair
      !!!!!!!!!! The r-space correlation functions in various channels
      real(rp), allocatable :: RealSpCrFSwp(:, :, :)
      real(rp), allocatable :: RealSpCrFBIN(:, :, :)
!________________________________________________________________________________________  
!_________________ (1) R-space and k-space correlation functions for all BINs ___________
!________________________________________________________________________________________
      !!!!!!!!!! R-space correlation functions for both OPEN and PERIODIC boundary conditions
      real(rp), allocatable :: RealSpCrFAll(:, :, :, :)   
      !!!!!!!!!! k-space structure factors for only PERIODIC boundary conditions
      complex(rp), allocatable :: KSpaceCrFAll(:, :, :, :)
!________________________________________________________________________________________  
!_________________ (2) Structure factors and Correlation ratios at specific _____________
!_____________________ k points for all BINs, for IfPyObsPBC == T case __________________
!________________________________________________________________________________________
      !!!!!!!!!! Structure factors at Gamma, X, M points
      real(rp), allocatable :: KStructFactGamm(:, :, :)
      real(rp), allocatable :: KStructFactXPnt(:, :, :)
      real(rp), allocatable :: KStructFactPntM(:, :, :)
      !!!!!!!!!! Correlation ratios    
      real(rp), allocatable :: CorlatRatioGamm(:, :, :)
      real(rp), allocatable :: CorlatRatioXPnt(:, :, :)
      real(rp), allocatable :: CorlatRatioPntM(:, :, :)
!**************************************************************************************************	 
!___________________ 4. Dynamic Correlations for PERIODIC and OPEN boundary conditions ____________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Number of imaginary time points in measurements __________________
!________________________________________________________________________________________  
      integer NumTauPnt
      integer, allocatable :: TauPntVal(:)      
      integer, allocatable :: TauPntVal0B(:)
!________________________________________________________________________________________  
!_________________ (1) Accumulate r-space dynamic correlation for single BIN ____________
!_____________________ For both PERIODIC and OPEN boundary conditions ___________________
!________________________________________________________________________________________
      integer NmSitePairTau
      integer, allocatable :: TauBondPairId(:, :)
      real(rp)   , allocatable :: RealSpCrFTauSwp(:, :, :)  ! The time-displaced quantities in single sweep
      real(rp)   , allocatable :: RealSpCrFTauBIN(:, :, :)  ! The time-displaced quantities in single BIN
      complex(rp), allocatable :: KSpCrFtTauforFT(:, :, :)  ! Temporary array for Fourier Transformation
!________________________________________________________________________________________  
!_________________ (2) Local dynamic correlations and G(k,beta/2) for all BINS __________
!_____________________ For both PERIODIC and OPEN boundary conditions ___________________
!________________________________________________________________________________________
      real(rp), allocatable :: RorKspCrFTauAll(:, :, :, :)  ! r-space or k-space dynamic correlations
      real(rp), allocatable :: RSpCrFLclTauAll(:, :, :)     ! The local correlation functions
      real(rp), allocatable :: Spn13pS12HalfBt(:      )     ! \beta * [ S_{13}(\beta/2) + S_{12}(\beta/2) ]
      real(rp), allocatable :: GkTauBetaOv2All(:, :, :)     ! The G(k,beta/2) for fermi surface
      real(rp), allocatable :: Dc_Conductivity(:      )     ! The direct current conductivity
!&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#
!_________________ (3) Dynamic correlation function in iw_n axis ________________________
!&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#&%@#
      !!!!!!!!!! Number of iw_n to be computed
      integer NmFrqFermi
      integer NmFrqBoson
      !!!!!!!!!! Fermionic channel
      real(rp), allocatable :: NonIntG0kIwn(   :, :)        ! The G_0(k,iw_n)
      real(rp), allocatable :: FermiGrF_Iwn(:, :, :)        ! The G(k,iw_n)
      real(rp), allocatable :: SelfEnrg_Iwn(:, :, :)        ! The \Sigma(k,iw_n)
      !!!!!!!!!! Bosonic channels
      real(rp), allocatable :: BosonCrF_Iwn(:, :, :, :)     ! The C(k,iw_n) bosonic channel
      !!!!!!!!!! The Quasi-particle weights
      real(rp), allocatable :: QuasiParWhgt(:, :)
      !!!!!!!!!! Drude conductivity and superfluid density
      real(rp), allocatable :: DrudeSupfdWt(:, :) 

   end module Observable
!________________________________________________________________________________________________________________________  
!____________________________________ End Module ________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$