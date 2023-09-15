!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL SPIN PINNING FIELDS; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A few subroutines for performing the data processing after the whole CPMC simulations, including:
!                  (0) Calculate the avgerage and errorbar for expectation energies of all the energy terms;
!                  (1) Calculate the avgerage for accepting ratio of three different interaction terms;
!                  (2) Calculate the avgerage and errorbar for real space correlation functions in only a few 
!                                         lattice site pairs;
!                  (3) Calculate the avgerage and errorbar for reciprocal space correlation functions.
! COMMENT: CPMC data process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   The subroutine used to process all STATIC and DYNAMIC physical observables.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   DQMCDataProc --> Subroutine to call subroutines to perform data processing and output the processed data;
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   The subroutines used to process all STATIC physical observables.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   FinlEngLcDenProc --> Process the energies and occupations;
!   FinlRKStaCrFProc --> Process the correlation functions under periodic boundary conditions;
!   FinlRSpCrFOBProc --> Process the correlation functions under open boundary conditions;
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   The subroutine used to process all DYNAMIC physical observables.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   FinlPBRSpCrFTauProc --> Process the dynamic correlation functions under periodic boundary conditions;
!   FinlRSpCrFTauOBProc --> Process the dynamic correlation functions under open boundary conditions;
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   The subroutine used to perform simple statistics.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   MeanErrBar     --> Compute average and errorbar for an array;   
!   MeanErrBar_Var --> Compute average and variance for an array.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!________________________________________ Begin subroutine ______________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine CPMCDataProc()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCDataProc()
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the data process for part of the results from CPMC.
! KEYWORDS: Data process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Data process after the whole CPMC simulation.
!
!     Input:  (none).   Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!______________________________________________________________________________________________________________   
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
        use RealPrecsn
        use CoreParamt      
        use Observable
        use MPISetting
      implicit none
!______________________________________________________________________________________________________________   
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer MeaInd
!______________________________________________________________________________________________________________   
!_________________________________ Main calculations for Data process _________________________________________
!______________________________________________________________________________________________________________
        if(amyid == amstr) then
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCDataProc process _________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
         write(*, "(A)") "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         write(*, "(2x, 'CPMCDataProc: The final data process!')")
!**************************************************************************************************   
!___________________ 0. Process data for Static Measure only at BetaT point _______________________
!**************************************************************************************************
         do MeaInd = 0, merge(1, 0, IfM2OneMea), +1
!________________________________________________________________________________________     
!_________________ (0) Process results of energies, occupations and local densities _____
!________________________________________________________________________________________
            call FinlEngLcDenProc(MeaInd)
!________________________________________________________________________________________     
!_________________ (1) Process results of Results of momentum distribution ______________
!_____________________ and Condensation fractions, pairing wavefunctions ________________
!________________________________________________________________________________________
            if(IfFftEnPar) then
               call FinlNkPairWfProc(MeaInd)
            end if
!________________________________________________________________________________________     
!_________________ (2) Process results of real-space correlations _______________________
!_____________________ for PERIODIC and OPEN boundary conditions ________________________
!________________________________________________________________________________________
            if(abs(PinSz) < rp_Eps) then
               call FinlRKStaCrFProc(MeaInd)
            end if
!**************************************************************************************************   
!___________________ 1. Process data for Dynamic Measure at BetaT + [BetaT, 0] ____________________
!**************************************************************************************************
            if( (MeaInd == 1) .and. (IfTAU) ) then
!________________________________________________________________________________________     
!_________________ (0) Compute average and errorbars for dynamic observables ____________
!________________________________________________________________________________________         
               call FinlRKDynCrFTauProc()
            end if
         enddo
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCDataProc process _________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(         
         write(*, "(A)") "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         write(*, "()")
         write(*, "()")
      end if
      
   end subroutine CPMCDataProc
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!################################## Average and Errorbars for STATIC observables ########################################
!################################## Average and Errorbars for STATIC observables ########################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine FinlEngLcDenProc(MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FinlEngLcDenProc(MeaInd)
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the final data process for static energies and occupations.
! KEYWORDS: Data process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input: MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
!
!     Outpt: (none).
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
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer MeaInd
!______________________________________________________________________________________________________________   
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, Ix, Iy
      real(rp) Avgr, ErrB
      real(rp) EngAvg(40, 2)          ! Mean values and error bar for the energies and Double Occupancy
!______________________________________________________________________________________________________________   
!_________________________________ Main calculations for Data process _________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************   
!___________________ 0. The QMC average of energies, occupations and local densities ______________
!**************************************************************************************************
!________________________________________________________________________________________     
!_________________ (0) Obtain the mean values and error bars ____________________________
!________________________________________________________________________________________
      call EngDocMean(MeaInd, EngAvg)
!________________________________________________________________________________________     
!_________________ (1) Output the energies, derivatives, occupations ____________________
!________________________________________________________________________________________  
!____________________________________________________________________________     
!________________ [0] First of all, Save all simulation parameters __________
!____________________________________________________________________________
      open( 289, file = Trim(FileOutAdd(MeaInd)) // "05_QMCSimultParamt" // Trim(FileAddTxt(MeaInd)))
      write(289, "(I4, A, I4, A, I6, A )", advance = "no") NumL1, char(9), NumL2, char(9), NumNS, char(9)
      write(289, "(es17.8, A, es17.8, A, es17.8, A)", advance = "no") BetaT, char(9), TempT, char(9), Dltau, char(9)
      write(289, "(es17.8, A, es17.8, A, I4, A)", advance = "no") Hopt1Up, char(9), Hopt1Dn, char(9), EkDispType, char(9)
      write(289, "(es17.8, A, es17.8, A)", advance = "no") Hopt2, char(9), Hopt3, char(9)
      write(289, "(es17.8, A)", advance = "no") ZmFdz, char(9)
      write(289, "(es17.8, A, I4, A)", advance = "no") PinSz, char(9), PinSzType, char(9)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Add sinusoidal spin pinning fileds
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(289, "(es17.8, A, es17.8, A)", advance = "no") SinusoidalPinSz, char(9), LambdaSz, char(9)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      if(EkDispType == 0) then
         write(289, "(es17.8, A)", advance = "no") ChemP, char(9)
      else if(EkDispType >= 1) then
         write(289, "(es17.8, A)", advance = "no") 4.0_rp+HubbU/2.0_rp-ChemP, char(9)
      end if
      write(289, "(es17.8, A           )", advance = "no") ChemP_BT, char(9)
      write(289, "(es17.8, A, es17.8   )", advance = "no") HubbU, char(9), HubbU_UHF
      write(289, "()")
      close(289)
!____________________________________________________________________________     
!________________ [1] Expectation Energies for all terms ____________________
!____________________________________________________________________________
      open(291, file = Trim(FileOutAdd(MeaInd)) // "05_QMCEnergyAverag" // Trim(FileAddTxt(MeaInd)))
      do I1 = 1, 8, +1
         write(291, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(I1, 1), char(9), EngAvg(I1, 2)
      enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Store energy that corresponds to the added spin pinning fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(291, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(22, 1), char(9), EngAvg(22, 2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(291, "()")
      close(291)
!____________________________________________________________________________     
!________________ [2] First-order derivative of model parameters ____________
!____________________________________________________________________________
      open(293, file = Trim(FileOutAdd(MeaInd)) // "05_QMCFstOrdDifPam" // Trim(FileAddTxt(MeaInd)))
      do I1 = 16, 17, +1
         write(293, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(I1, 1), char(9), EngAvg(I1, 2)
      enddo
      write(293, "(A)", advance = "no") char(9)
      do I1 = 18, 21, +1
         write(293, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(I1, 1), char(9), EngAvg(I1, 2)
      enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Store the first-order derivative of the term that corresponds to the added spin pinning fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      write(293, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(23, 1), char(9), EngAvg(23, 2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      write(293, "()")
      close(293)
!____________________________________________________________________________     
!________________ [3] Occupation, DouOcc, Some correlations _________________
!____________________________________________________________________________
      open(295, file = Trim(FileOutAdd(MeaInd)) // "05_QMCOccDouSpnCrF" // Trim(FileAddTxt(MeaInd)))
      do I1 = 31, 34, +1
         write(295, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(I1, 1), char(9), EngAvg(I1, 2)
      enddo
      write(295, "(A)", advance = "no") char(9)
      do I1 = 35, 36, +1
         write(295, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(I1, 1), char(9), EngAvg(I1, 2)
      enddo
      write(295, "()")
      close(295)   
!____________________________________________________________________________     
!________________ [4] Total energy per particle in unit of E_{FG} ___________
!____________________________________________________________________________
      if(EkDispType >= 1) then
         open( 293, file = Trim(FileOutAdd(MeaInd)) // "05_QMCAvgEfermiEfg" // Trim(FileAddTxt(MeaInd)))
         write(293, "(I4)", advance = "no") NumL1
         write(293, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(09, 1), char(9), EngAvg(09, 2)
         write(293, "(A, es17.8, A, es17.8)", advance = "no") char(9), EngAvg(10, 1), char(9), EngAvg(10, 2)
         write(293, "()")
         close(293)
      end if     
!&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$
!________________________ Store the total electron density ______________________________
!&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$&%^$
      open( 899, file = Trim(FileOutAdd(MeaInd)) // "99_TotalDensity_BX" // Trim(FileAddTxt(MeaInd)))
      write(899, "(A, es25.16, A, es25.16)", advance = "no") char(9), EngAvg(31, 1), char(9), EngAvg(31, 2)
      write(899, "(A, es25.16, A, es25.16)", advance = "no") char(9), EngAvg(32, 1), char(9), EngAvg(32, 2)
      write(899, "(A, es25.16, A, es25.16)", advance = "no") char(9), EngAvg(33, 1), char(9), EngAvg(33, 2)
      write(899, "()")
      close(899)
!**************************************************************************************************   
!___________________ 1. Results of single-particle GrF matrices in r- and k-space _________________
!**************************************************************************************************
!________________________________________________________________________________________     
!_________________ (0) The r-space single-particle GrF matrix ___________________________
!________________________________________________________________________________________
      open(303, file = Trim(FileOutAdd(MeaInd)) // "66_RSpDnstyMatFinl" // Trim(FileAddTxt(MeaInd)))
      do I2 = 1, NumNC, +1
         do I1 = 1, I2, +1
            I0 = (I2-1)*I2/2 + I1
            write(303, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
            call MeanErrBar(Bin_Data, RlSpGrnFtAll(BinStart, I0, 1, MeaInd), Avgr, ErrB)
            write(303, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            write(303, "(A)", advance = "no") char(9)
            call MeanErrBar(Bin_Data, RlSpGrnFtAll(BinStart, I0, 2, MeaInd), Avgr, ErrB)
            write(303, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            write(303, "()")
         enddo
      enddo
      close(303)
!________________________________________________________________________________________     
!_________________ (1) The k-space single-particle GrF matrix ___________________________
!________________________________________________________________________________________
      if(IfFftEnPar) then
         open(383, file = Trim(FileOutAdd(MeaInd)) // "66_KDnstyMatrxFinl" // Trim(FileAddTxt(MeaInd)))
         do I2 = 1, NumNC, +1
            do I1 = 1, I2, +1
               I0 = (I2-1)*I2/2 + I1
               write(383, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
               call MeanErrBar(Bin_Data, KSpGreenFAll(BinStart, I0, 1, MeaInd), Avgr, ErrB)
               write(383, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, KSpGreenFAll(BinStart, I0, 2, MeaInd), Avgr, ErrB)
               write(383, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               write(383, "(A)", advance = "no") char(9)
               call MeanErrBar(Bin_Data, KSpGreenFAll(BinStart, I0, 3, MeaInd), Avgr, ErrB)
               write(383, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, KSpGreenFAll(BinStart, I0, 4, MeaInd), Avgr, ErrB)
               write(383, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               write(383, "()")
            enddo
         enddo
         close(383)
      end if
!**************************************************************************************************   
!___________________ 2. Results of local charge and spin densities ________________________________
!**************************************************************************************************
!________________________________________________________________________________________     
!_________________ (0) Local charge and spin densities for all sites ____________________
!________________________________________________________________________________________
      open(311, file = Trim(FileOutAdd(MeaInd)) // "06_RnUpnDw_AvgErrB" // Trim(FileAddTxt(MeaInd)))
      open(312, file = Trim(FileOutAdd(MeaInd)) // "06_RSpzSpx_AvgErrB" // Trim(FileAddTxt(MeaInd)))
      do Iy = 1, NumL2, +1
         do Ix = 1, NumL1, +1
            I1 = (Iy-1)*NumL1 + Ix
            !!!!!! local density for spin-up, down and total
            write(311, "(I6, A, I4, A, I4, A)", advance = "no") I1, char(9), Ix, char(9), Iy, char(9)
            call MeanErrBar(Bin_Data, RSpLclOrdAll(BinStart, I1, 1, MeaInd), Avgr, ErrB)
            write(311, "(A, es25.16, A, es25.16)", advance = "no") char(9), Avgr, char(9), ErrB
            call MeanErrBar(Bin_Data, RSpLclOrdAll(BinStart, I1, 2, MeaInd), Avgr, ErrB)
            write(311, "(A, es25.16, A, es25.16)", advance = "no") char(9), Avgr, char(9), ErrB
            call MeanErrBar(Bin_Data, RSpLclOrdAll(BinStart, I1, 3, MeaInd), Avgr, ErrB)
            write(311, "(A, es25.16, A, es25.16)", advance = "no") char(9), Avgr, char(9), ErrB
            write(311, "()")
            !!!!!! local spin density
            write(312, "(I6, A, I4, A, I4, A)", advance = "no") I1, char(9), Ix, char(9), Iy, char(9)
            call MeanErrBar(Bin_Data, RSpLclOrdAll(BinStart, I1, 4, MeaInd), Avgr, ErrB)
            write(312, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            write(312, "()")
         enddo
      enddo
      close(311)
      close(312)
!________________________________________________________________________________________     
!_________________ (1) Local charge and spin densities for all rows _____________________
!________________________________________________________________________________________
      if(mod(NumNS, 2) == 0) then
         open(313, file = Trim(FileOutAdd(MeaInd)) // "06_RLocalnTotalRow" // Trim(FileAddTxt(MeaInd)))
         open(314, file = Trim(FileOutAdd(MeaInd)) // "06_RLocalSpinSzRow" // Trim(FileAddTxt(MeaInd)))
         do Ix = 1, NumL1, +1
            write(313, "(I4, A)", advance = "no") Ix, char(9)
            write(314, "(I4, A)", advance = "no") Ix, char(9)
            do Iy = 1, NumL2, +1
               I1 = (Iy-1)*NumL1 + Ix
               call MeanErrBar(Bin_Data, RSpLclOrdAll(BinStart, I1, 3, MeaInd), Avgr, ErrB)
               write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, RSpLclOrdAll(BinStart, I1, 4, MeaInd), Avgr, ErrB)
               write(314, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            enddo
            write(313, "()")
            write(314, "()")
         enddo
         close(313)
         close(314)
      end if
!________________________________________________________________________________________     
!_________________ (2) Average Local spin Sz moment _____________________________________
!________________________________________________________________________________________
      if(mod(NumNS, 2) == 0) then
         open( 313, file = Trim(FileOutAdd(MeaInd)) // "06_TotalMagntMomnt" // Trim(FileAddTxt(MeaInd)))
         call MeanErrBar(Bin_Data, TotMagMmtAll(BinStart, MeaInd), Avgr, ErrB)
         write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         close(313)
      end if

   end subroutine FinlEngLcDenProc
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine EngDocMean(MeaInd, EngAvg)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  EngDocMean(MeaInd, EngAvg)
! TYPE:     subroutine
! PURPOSE:  This Subroutine calculates the mean values and error bar for the energies and double occupancy.
! KEYWORDS: Data process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: QMC average for energies, occupations and s-s correlation.
!
!     Input:  EngAvg --> The output average and errorbar. 
!             MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
!  
!     Output: (none).
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
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer MeaInd
        real(rp) EngAvg(40, 2)
!______________________________________________________________________________________________________________   
!_________________________________ Main calculations for Data process _________________________________________
!______________________________________________________________________________________________________________
      call MeanErrBar(Bin_Data, EnHopt1(BinStart, MeaInd), EngAvg(01, 1), EngAvg( 1, 2))
      call MeanErrBar(Bin_Data, EnHopt2(BinStart, MeaInd), EngAvg(02, 1), EngAvg( 2, 2))
      call MeanErrBar(Bin_Data, EnHopt3(BinStart, MeaInd), EngAvg(03, 1), EngAvg( 3, 2))
      call MeanErrBar(Bin_Data, EnZmFld(BinStart, MeaInd), EngAvg(04, 1), EngAvg( 4, 2))
      call MeanErrBar(Bin_Data, EnPinSz(BinStart, MeaInd), EngAvg(05, 1), EngAvg( 5, 2))
      call MeanErrBar(Bin_Data, EnDopCh(BinStart, MeaInd), EngAvg(06, 1), EngAvg( 6, 2))
      call MeanErrBar(Bin_Data, EnHubbU(BinStart, MeaInd), EngAvg(07, 1), EngAvg( 7, 2))
      call MeanErrBar(Bin_Data, EnTotal(BinStart, MeaInd), EngAvg(08, 1), EngAvg( 8, 2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Take average of the energy that corresponds to the added spin pinning fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call MeanErrBar(Bin_Data, EnSinusoidalPinSz(BinStart, MeaInd), EngAvg(22, 1), EngAvg(22, 2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if(EkDispType >= 1) then
         call MeanErrBar(Bin_Data, EpOvEfg(BinStart, 1, MeaInd), EngAvg(09, 1), EngAvg(09, 2))
         call MeanErrBar(Bin_Data, EpOvEfg(BinStart, 2, MeaInd), EngAvg(10, 1), EngAvg(10, 2))
      end if
      
      call MeanErrBar(Bin_Data, EnHbUCh(BinStart, MeaInd), EngAvg(16, 1), EngAvg(16, 2))
      call MeanErrBar(Bin_Data, EnTotCh(BinStart, MeaInd), EngAvg(17, 1), EngAvg(17, 2))
      call MeanErrBar(Bin_Data, HmtOvt2(BinStart, MeaInd), EngAvg(18, 1), EngAvg(18, 2))
      call MeanErrBar(Bin_Data, HmtOvt3(BinStart, MeaInd), EngAvg(19, 1), EngAvg(19, 2))
      call MeanErrBar(Bin_Data, HmtOvZm(BinStart, MeaInd), EngAvg(20, 1), EngAvg(20, 2))
      call MeanErrBar(Bin_Data, HmtOvbU(BinStart, MeaInd), EngAvg(21, 1), EngAvg(21, 2)) 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Take average of the energy that corresponds to the added spin pinning fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call MeanErrBar(Bin_Data, HmOvSinusoidalPinSz(BinStart, MeaInd), EngAvg(23, 1), EngAvg(23, 2))
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      call MeanErrBar(Bin_Data, nOccpUp(BinStart, MeaInd), EngAvg(31, 1), EngAvg(31, 2))
      call MeanErrBar(Bin_Data, nOccpDw(BinStart, MeaInd), EngAvg(32, 1), EngAvg(32, 2))
      call MeanErrBar(Bin_Data, nOccTot(BinStart, MeaInd), EngAvg(33, 1), EngAvg(33, 2))
      call MeanErrBar(Bin_Data, NmNeTot(BinStart, MeaInd), EngAvg(34, 1), EngAvg(34, 2))
      call MeanErrBar(Bin_Data, RDouOcc(BinStart, MeaInd), EngAvg(35, 1), EngAvg(35, 2))
      call MeanErrBar(Bin_Data, RNNDCrF(BinStart, MeaInd), EngAvg(36, 1), EngAvg(36, 2))

   end subroutine EngDocMean
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine FinlNkPairWfProc(MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FinlNkPairWfProc(MeaInd)
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the final data process for n(k) and pairing wavefunction for FFT. 
! KEYWORDS: Data process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION:
!
!     Input: MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
!
!     Outpt: (none).
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
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer MeaInd
!______________________________________________________________________________________________________________   
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, Ix, Iy, Nk
      real(rp) Avgr, ErrB
!______________________________________________________________________________________________________________   
!_________________________________ Main calculations for Data process _________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************   
!___________________ 0. Process the results of momentum distribution ______________________________
!**************************************************************************************************
      open(385, file = Trim(FileOutAdd(MeaInd)) // "MomtumDistFinl" // Trim(FileAddTxt(MeaInd)))
      do Nk = 1, NumNC, +1
         write(385, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
         call MeanErrBar(Bin_Data, NkDistrib(BinStart, Nk, 1, MeaInd), Avgr, ErrB)
         write(385, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, NkDistrib(BinStart, Nk, 2, MeaInd), Avgr, ErrB)
         write(385, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(385, "()")
      enddo
      close(385)
!**************************************************************************************************   
!___________________ 1. Process the results of pairing wavefunctions ______________________________
!**************************************************************************************************
      !!!!!!!!!! The condensate fraction
      open( 383, file = Trim(FileOutAdd(MeaInd)) // "CondstFracFinl" // Trim(FileAddTxt(MeaInd)))
      call MeanErrBar(Bin_Data, Condfract(BinStart, 1, MeaInd), Avgr, ErrB)
      write(383, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
      call MeanErrBar(Bin_Data, Condfract(BinStart, 2, MeaInd), Avgr, ErrB)
      write(383, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
      write(383, "(A, es17.8, A, es17.8)") char(9), Condfract(NmBin+1, 1, MeaInd), char(9), Condfract(NmBin+1, 2, MeaInd)
      close(383)
      !!!!!!!!!! The pairing wavefunction in k and r spaces
      open(385, file = Trim(FileOutAdd(MeaInd)) // "PairWvfcMtFinl" // Trim(FileAddTxt(MeaInd)))
      open(387, file = Trim(FileOutAdd(MeaInd)) // "PairWvfcRlFinl" // Trim(FileAddTxt(MeaInd)))
      do I0 = 1, NumNC, +1
         !!!!! Pairing wavefunction in k space
         write(385, "(I4, A, I4, A)", advance = "no") KpList(I0, 1), char(9), KpList(I0, 2), char(9)
         call MeanErrBar(Bin_Data, PairWvfct(BinStart, I0, 1, MeaInd), Avgr, ErrB)
         write(385, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, PairWvfct(BinStart, I0, 2, MeaInd), Avgr, ErrB)
         write(385, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, PairWvfct(BinStart, I0, 3, MeaInd), Avgr, ErrB)
         write(385, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(385, "(A, es17.8, A, es17.8, A, es17.8)") char(9), PairWvfct(NmBin+1, I0, 1, MeaInd), &
            & char(9), PairWvfct(NmBin+1, I0, 2, MeaInd), char(9), PairWvfct(NmBin+1, I0, 3, MeaInd)
         !!!!! Pairing wavefunction in r space
         Ix = mod(I0-1, NumL1); Iy = (I0-1)/NumL1
         write(387, "(I4, A, I4, A)", advance = "no") Ix-NumL1/2, char(9), Iy-NumL2/2, char(9)
         call MeanErrBar(Bin_Data, PairWvfct(BinStart, I0, 4, MeaInd), Avgr, ErrB)
         write(387, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, PairWvfct(BinStart, I0, 5, MeaInd), Avgr, ErrB)
         write(387, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, PairWvfct(BinStart, I0, 6, MeaInd), Avgr, ErrB)
         write(387, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(387, "(A, es17.8, A, es17.8, A, es17.8)") char(9), PairWvfct(NmBin+1, I0, 4, MeaInd), &
            & char(9), PairWvfct(NmBin+1, I0, 5, MeaInd), char(9), PairWvfct(NmBin+1, I0, 6, MeaInd)
      enddo
      close(385)
      close(387)
   
   end subroutine FinlNkPairWfProc
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine FinlRKStaCrFProc(MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FinlRKStaCrFProc(MeaInd)
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the final data process for r-space correlation functions for PERIODIC
!                 boundary conditions.
! KEYWORDS: Data process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input: MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
!
!     Outpt: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!______________________________________________________________________________________________________________   
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
        use RealPrecsn
        use CoreParamt      
        use Observable
        use MPISetting
        implicit none
!______________________________________________________________________________________________________________   
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer MeaInd
!______________________________________________________________________________________________________________   
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
        integer Nk, I1, I2, Idimj
        real(rp) Avgr, ErrB
      real(rp) KCrFAvg(2)             ! Average values used in momentum space correlation function
      real(rp) KCrFErr(2)             ! Errorbar used in momentum space correlation function
!______________________________________________________________________________________________________________   
!_________________________________ Main calculations for Data process _________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************   
!___________________ 0. Output the results for the present BIN simulation _________________________
!______________________ for both PERIODIC and OPEN boundary conditions ____________________________
!**************************************************************************************************
!________________________________________________________________________________________     
!_________________ (0) Output all the r-space correlation functions _____________________
!_____________________ for both IfPyObsPBC==T and IfPyObsPBC==F cases ___________________
!________________________________________________________________________________________
      !!!!!!!!!! Open the files for the output
      open(304, file = Trim(FileOutAdd(MeaInd)) // "66_RGrnFctCrFc_All" // Trim(FileAddTxt(MeaInd)))
      open(305, file = Trim(FileOutAdd(MeaInd)) // "66_RSpinZZCrFc_All" // Trim(FileAddTxt(MeaInd)))
      open(306, file = Trim(FileOutAdd(MeaInd)) // "66_RSpinPMCrFc_All" // Trim(FileAddTxt(MeaInd)))
      open(307, file = Trim(FileOutAdd(MeaInd)) // "66_RDenDenCrFc_All" // Trim(FileAddTxt(MeaInd)))
      open(308, file = Trim(FileOutAdd(MeaInd)) // "66_RPairStCrFc_All" // Trim(FileAddTxt(MeaInd)))
      open(309, file = Trim(FileOutAdd(MeaInd)) // "66_REdSParCrFc_All" // Trim(FileAddTxt(MeaInd)))
      open(310, file = Trim(FileOutAdd(MeaInd)) // "66_RDWvParCrFc_All" // Trim(FileAddTxt(MeaInd)))
      open(311, file = Trim(FileOutAdd(MeaInd)) // "66_RBndPairCrF_All" // Trim(FileAddTxt(MeaInd)))
      !!!!!!!!!! Initialization of (I1, I2) for IfPyObsPBC == F case
      I2 = 1; I1 = 0
      !!!!!!!!!! Output all the r-space correlations
      do Idimj = 1, NmSitePair, +1
         !!!!!!!! The (I1, I2) index for site pairs
         if(IfPyObsPBC) then
            I1 = StList(Idimj, 1); I2 = StList(Idimj, 2)
            I1 = I1 - NumL1/2 - mod(NumL1, 2); I2 = I2 - NumL2/2 - mod(NumL2, 2)
         else
            I1 = I1 + 1
            if(I1 > I2) then
               I2 = I2 + 1; I1 = 1
            end if
         end if
         !!!!!!!! Compute and Write the results to files
         !!!! Single-particle Green's functions
         write(304, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 01, MeaInd), Avgr, ErrB)
         write(304, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(304, "(A)", advance = "no") char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 02, MeaInd), Avgr, ErrB)
         write(304, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(304, "()")
         !!!! Spin Sz-Sz correlations
         write(305, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
         !! The standard correlation
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 03, MeaInd), Avgr, ErrB)
         write(305, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         !! The Vertex contribution
         write(305, "(A)", advance = "no") char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 17, MeaInd), Avgr, ErrB)
         write(305, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(305, "()")
         !!!! Spin Spm-Spm correlations
         write(306, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
         !! The standard correlation
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 04, MeaInd), Avgr, ErrB)
         write(306, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         !! The Vertex contribution
         write(306, "(A)", advance = "no") char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 18, MeaInd), Avgr, ErrB)
         write(306, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(306, "()")
         !!!! Density-density correlations
         write(307, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
         !! The standard correlation
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 05, MeaInd), Avgr, ErrB)
         write(307, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         !! The Vertex contribution
         write(307, "(A)", advance = "no") char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 19, MeaInd), Avgr, ErrB)
         write(307, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(307, "()")
         !!!! on-site spin-singlet s-wave pairing
         write(308, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
         !! The standard correlation
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 06, MeaInd), Avgr, ErrB)
         write(308, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         !! The Vertex contribution
         write(308, "(A)", advance = "no") char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 20, MeaInd), Avgr, ErrB)
         write(308, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(308, "()")
         !!!! NN extended spin-singlet s-wave pairing
         write(309, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
         !! The standard correlation
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 07, MeaInd), Avgr, ErrB)
         write(309, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         !! The Vertex contribution
         write(309, "(A)", advance = "no") char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 21, MeaInd), Avgr, ErrB)
         write(309, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(309, "()")
         !!!! NN spin-singlet d-wave pairing
         write(310, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
         !! The standard correlation
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 08, MeaInd), Avgr, ErrB)
         write(310, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         !! The Vertex contribution
         write(310, "(A)", advance = "no") char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 22, MeaInd), Avgr, ErrB)
         write(310, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(310, "()")
         !!!! Bond-Bond spin-singlet pairing
         write(311, "(I4, A, I4, A)", advance = "no") I1, char(9), I2, char(9)
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 31, MeaInd), Avgr, ErrB)
         write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 32, MeaInd), Avgr, ErrB)
         write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 33, MeaInd), Avgr, ErrB)
         write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, RealSpCrFAll(BinStart, Idimj, 34, MeaInd), Avgr, ErrB)
         write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         write(311, "()")
      enddo
      close(304); close(305); close(306); close(307); close(308); close(309); close(310); close(311)
!________________________________________________________________________________________     
!_________________ (1) Output the structure factors at all k points _____________________
!_____________________ only for IfPyObsPBC==T case ______________________________________
!________________________________________________________________________________________
      if(IfPyObsPBC) then
         !!!!!!!!!! Open txt files to for storing the results
         open(313, file = Trim(FileOutAdd(MeaInd)) // "07_KGrnFct_AvgErrB" // Trim(FileAddTxt(MeaInd)))
         open(323, file = Trim(FileOutAdd(MeaInd)) // "07_KSpinZZ_AvgErrB" // Trim(FileAddTxt(MeaInd)))
         open(333, file = Trim(FileOutAdd(MeaInd)) // "07_KSpinPM_AvgErrB" // Trim(FileAddTxt(MeaInd)))
         open(343, file = Trim(FileOutAdd(MeaInd)) // "07_KDenDen_AvgErrB" // Trim(FileAddTxt(MeaInd)))
         open(353, file = Trim(FileOutAdd(MeaInd)) // "07_KPairSt_AvgErrB" // Trim(FileAddTxt(MeaInd)))
         open(363, file = Trim(FileOutAdd(MeaInd)) // "07_KEdSPar_AvgErrB" // Trim(FileAddTxt(MeaInd)))
         open(373, file = Trim(FileOutAdd(MeaInd)) // "07_KDWvPar_AvgErrB" // Trim(FileAddTxt(MeaInd)))
         !!!!!!!!!! Iteration of all k points for all correlation functions
         do Nk = 1, NumNC, +1
            !!!!!!!! KGFUpUp and KGFDwDw
            write(313, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9) 
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 01, MeaInd), KCrFAvg, KCrFErr)
            write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            write(313, "(A)", advance = "no") char(9)
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 02, MeaInd), KCrFAvg, KCrFErr)
            write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            write(313, "()")
            !!!!!!!! KSpinZZ
            write(323, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            !!!!!! The standard correlation
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 03, MeaInd), KCrFAvg, KCrFErr)
            write(323, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(323, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            !!!!!! The Vertex contribution
            write(323, "(A)", advance = "no") char(9)
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 17, MeaInd), KCrFAvg, KCrFErr)
            write(323, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(323, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            write(323, "()")        
            !!!!!!!! KSpinPM
            write(333, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            !!!!!! The standard correlation
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 04, MeaInd), KCrFAvg, KCrFErr)
            write(333, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(333, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            !!!!!! The Vertex contribution
            write(333, "(A)", advance = "no") char(9)
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 18, MeaInd), KCrFAvg, KCrFErr)
            write(333, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(333, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            write(333, "()")
            !!!!!!!! KDenDen
            write(343, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            !!!!!! The standard correlation
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 05, MeaInd), KCrFAvg, KCrFErr)
            write(343, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(343, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            !!!!!! The Vertex contribution
            write(343, "(A)", advance = "no") char(9)
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 19, MeaInd), KCrFAvg, KCrFErr)
            write(343, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(343, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            write(343, "()")
            !!!!!!!! KPairSt
            write(353, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            !!!!!! The standard correlation
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 06, MeaInd), KCrFAvg, KCrFErr)
            write(353, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(353, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            !!!!!! The Vertex contribution
            write(353, "(A)", advance = "no") char(9)
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 20, MeaInd), KCrFAvg, KCrFErr)
            write(353, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(353, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            write(353, "()")
            !!!!!!!! KEdSPar
            write(363, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            !!!!!! The standard correlation
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 07, MeaInd), KCrFAvg, KCrFErr)
            write(363, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(363, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            !!!!!! The Vertex contribution
            write(363, "(A)", advance = "no") char(9)
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 21, MeaInd), KCrFAvg, KCrFErr)
            write(363, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(363, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            write(363, "()")
            !!!!!!!! KDWvPar
            write(373, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            !!!!!! The standard correlation
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 08, MeaInd), KCrFAvg, KCrFErr)
            write(373, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(373, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            !!!!!! The Vertex contribution
            write(373, "(A)", advance = "no") char(9)
            call KCorrFMeanA(KSpaceCrFAll(1, Nk, 22, MeaInd), KCrFAvg, KCrFErr)
            write(373, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(1), char(9), KCrFErr(1)
            write(373, "(A, es17.8, A, es17.8)", advance = "no") char(9), KCrFAvg(2), char(9), KCrFErr(2)
            write(373, "()")
         enddo
         !!!!!!!!!! Close the files
         close(313); close(323); close(333); close(343); close(353); close(363); close(373)
      end if
!________________________________________________________________________________________     
!_________________ (2) Output the structure factors and correlation ratios ______________
!_____________________ only for IfPyObsPBC==T case ______________________________________
!________________________________________________________________________________________
      if(IfPyObsPBC) then
!____________________________________________________________________________     
!________________ [0] The structure factors --> Consider symmtry ____________
!____________________________________________________________________________
         !!!!!!!!!! At Gamma point
         open(311, file = Trim(FileOutAdd(MeaInd)) // "09_KStructFactGamm" // Trim(FileAddTxt(MeaInd)))
         do I1 = 3, 8, +1
            call MeanErrBar(Bin_Data, KStructFactGamm(BinStart, I1, MeaInd), Avgr, ErrB)
            write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         enddo
         write(311, "()")
         close(311)
         !!!!!!!!!! At X point
         open(311, file = Trim(FileOutAdd(MeaInd)) // "10_KStructFactXPnt" // Trim(FileAddTxt(MeaInd)))
         do I1 = 3, 8, +1
            call MeanErrBar(Bin_Data, KStructFactXPnt(BinStart, I1, MeaInd), Avgr, ErrB)
            write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         enddo
         write(311, "()")
         close(311)
         !!!!!!!!!! At M point
         open(311, file = Trim(FileOutAdd(MeaInd)) // "11_KStructFactPntM" // Trim(FileAddTxt(MeaInd)))
         do I1 = 3, 8, +1
            call MeanErrBar(Bin_Data, KStructFactPntM(BinStart, I1, MeaInd), Avgr, ErrB)
            write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         enddo
         write(311, "()")
         close(311)   
!____________________________________________________________________________     
!________________ [1] The correlation ratios --> Consider symmtry ___________
!____________________________________________________________________________
         !!!!!!!!!! At Gamma point
         open(311, file = Trim(FileOutAdd(MeaInd)) // "09_CorlatRatioGamm" // Trim(FileAddTxt(MeaInd)))
         do I1 = 3, 8, +1
            call MeanErrBar(Bin_Data, CorlatRatioGamm(BinStart, I1, MeaInd), Avgr, ErrB)
            write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         enddo
         write(311, "()")
         close(311)
         !!!!!!!!!! At X point
         open(311, file = Trim(FileOutAdd(MeaInd)) // "10_CorlatRatioXPnt" // Trim(FileAddTxt(MeaInd)))
         do I1 = 3, 8, +1
            call MeanErrBar(Bin_Data, CorlatRatioXPnt(BinStart, I1, MeaInd), Avgr, ErrB)
            write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         enddo
         write(311, "()")
         close(311)
         !!!!!!!!!! At M point
         open(311, file = Trim(FileOutAdd(MeaInd)) // "11_CorlatRatioPntM" // Trim(FileAddTxt(MeaInd)))
         do I1 = 3, 8, +1
            call MeanErrBar(Bin_Data, CorlatRatioPntM(BinStart, I1, MeaInd), Avgr, ErrB)
            write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         enddo
         write(311, "()")
         close(311)
!________________________________________________________________________________________     
!_________________ (2) Momentum distribution For IfFftEnPar == .false. case _____________
!_____________________ only for IfPyObsPBC==T case ______________________________________
!________________________________________________________________________________________
         if( (IfPyObsPBC) .and. (.not. IfFftEnPar) ) then
            open(385, file = Trim(FileOutAdd(MeaInd)) // "MomtumDistFinl" // Trim(FileAddTxt(MeaInd)))
            do Nk = 1, NumNC, +1
               write(385, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
               call MeanErrBar(Bin_Data, NkDistrib(BinStart, Nk, 1, MeaInd), Avgr, ErrB)
               write(385, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, NkDistrib(BinStart, Nk, 2, MeaInd), Avgr, ErrB)
               write(385, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               write(385, "()")
            enddo
            close(385)
         end if
      end if

   end subroutine FinlRKStaCrFProc
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine KCorrFMeanA(CVector, KCrFAvg, KCrFErr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  KCorrFMeanA(CVector, KCrFAvg, KCrFErr)
! TYPE:     subroutine
! PURPOSE:  This Subroutine calculates the mean values and error bar for the correlation functions in reciprocal
!               space including:
!                 (1) The <Sz*Sz> correlation function;
!                 (2) The <S+*S- + S-*S+> correlation function;
!                 (3) The <n_i * n_j> correlation function;
!                 (4) Momentum distribution.
! KEYWORDS: Data process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input:  NDim    --> Input integer as dimension of CMatrix(N, 4) matrix;
!             CMatrix --> Input complex matrix for calculating average and errorbar;
!
!     Output: KCrFAvg --> The output average values for real and complex parts of CMatrix;
!             KCrFErr --> The output errorbar values for real and complex parts of CMatrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!______________________________________________________________________________________________________________   
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
        use RealPrecsn
      use CoreParamt
        implicit none
!______________________________________________________________________________________________________________   
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
        complex(rp) CVector(NmBIN)
        real(rp) KCrFAvg(2)
        real(rp) KCrFErr(2)
      real(rp) RealVec(NmBIN), ImagVec(NmBIN)
!______________________________________________________________________________________________________________   
!_________________________________ Main calculations for Data process _________________________________________
!______________________________________________________________________________________________________________
      RealVec(1:NmBIN) =  real(CVector(1:NmBIN))
      ImagVec(1:NmBIN) = aimag(CVector(1:NmBIN))
        call MeanErrBar(Bin_Data, RealVec(BinStart), KCrFAvg(1), KCrFErr(1))
        call MeanErrBar(Bin_Data, ImagVec(BinStart), KCrFAvg(2), KCrFErr(2))

    end subroutine KCorrFMeanA
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine KCorrFMeanB(NOrbtSquare, CMatrix, KCrFAvg, KCrFErr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  KCorrFMeanB(NOrbtSquare, CMatrix, KCrFAvg, KCrFErr)
! TYPE:     subroutine
! PURPOSE:  This Subroutine calculates the mean values and error bar for the correlation functions in reciprocal 
!                   space including:
!                 (1) The <Sz*Sz> correlation function;
!                 (2) The <S+*S- + S-*S+> correlation function;
!                 (3) The <n_i * n_j> correlation function;
!                 (4) Momentum distribution.
! KEYWORDS: Data process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     V interaction updating.
!
!     Input:  NDim1   --> Input integer as dimension of CMatrix(NDim1, NDim2) matrix;
!             NDim2   --> Input integer as dimension of CMatrix(NDim1, NDim2) matrix;
!             CMatrix --> Input complex matrix for calculating average and errorbar;
!
!     Output: KCrFAvg --> The output average values for real and complex parts of CMatrix;
!             KCrFErr --> The output errorbar values for real and complex parts of CMatrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!______________________________________________________________________________________________________________   
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
        use RealPrecsn
      use CoreParamt
        implicit none
!______________________________________________________________________________________________________________   
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
        integer NOrbtSquare
        complex(rp) CMatrix(NmBIN, NOrbtSquare)
        real(rp) KCrFAvg(2*NOrbtSquare)
        real(rp) KCrFErr(2*NOrbtSquare)
!______________________________________________________________________________________________________________   
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
        integer Id
      real(rp) RealVec(NmBIN), ImagVec(NmBIN)
!______________________________________________________________________________________________________________   
!_________________________________ Main calculations for Data process _________________________________________
!______________________________________________________________________________________________________________
        do Id = 1, NOrbtSquare, +1
         RealVec(1:NmBIN) = real(CMatrix(1:NmBIN, Id))
         ImagVec(1:NmBIN) = aimag(CMatrix(1:NmBIN, Id))
            call MeanErrBar(Bin_Data, RealVec(BinStart), KCrFAvg(2*Id-1), KCrFErr(2*Id-1))
            call MeanErrBar(Bin_Data, ImagVec(BinStart), KCrFAvg(2*Id  ), KCrFErr(2*Id  ))
        enddo

    end subroutine KCorrFMeanB
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!################################## Average and Errorbars for DYNAMIC observables #######################################
!################################## Average and Errorbars for DYNAMIC observables #######################################
!########################################################################################################################
!######################################################################################################################## 


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine FinlRKDynCrFTauProc()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FinlRKDynCrFTauProc()
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the final data process for dynamic correlations with both OPEN and PERIODIC
!                     boundary conditions.
! KEYWORDS: Data process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION:
!
!     Input:  (none).   Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!______________________________________________________________________________________________________________   
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
        use RealPrecsn
        use CoreParamt      
        use Observable
        use MPISetting
        implicit none
!______________________________________________________________________________________________________________   
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
        integer NTInd, Nk, FrqInd, I0, I1
      real(rp) Avgr, ErrB
!______________________________________________________________________________________________________________   
!_________________________________ Main calculations for Data process _________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************   
!___________________ 0. Dynamic correlations for OPEN and PERIODIC boundary conditions ____________
!**************************************************************************************************    
!________________________________________________________________________________________     
!_________________ (0) Output local dynamic correlations ________________________________
!_____________________ only for IfPyObsPBC == T case ____________________________________
!________________________________________________________________________________________
      if(IfPyObsPBC) then
         !!!!!!!!!! The local correlation functions versus. \tau
         open(342, file = "Add_Output/08_LocalGrnFctcF_Tau.txt")
         open(343, file = "Add_Output/08_LocalSpzDencF_Tau.txt")
         open(344, file = "Add_Output/08_FNNSNNSpinCrF_Tau.txt")
         do NTInd = 0, NumTauPnt, +1
            !!!!!!!! Local single-particle Green's function
            write(342, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
            do I0 = 1, 4, +1
               call MeanErrBar(Bin_Data, RSpCrFLclTauAll(BinStart, NTInd, I0), Avgr, ErrB)
               write(342, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            enddo
            write(342, "()")
            !!!!!!!! Local Sz-Sz and Density-density correlation function
            write(343, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
            do I0 = 5, 8, +1
               call MeanErrBar(Bin_Data, RSpCrFLclTauAll(BinStart, NTInd, I0), Avgr, ErrB)
               write(343, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            enddo
            write(343, "()")
            !!!!!!!! The Sz-Sz correlations for FNN and SNN
            write(344, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
            do I0 = 9, 10, +1
               call MeanErrBar(Bin_Data, RSpCrFLclTauAll(BinStart, NTInd, I0), Avgr, ErrB)
               write(344, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            enddo
            write(344, "()")
         enddo
         close(342); close(343); close(344)
         !!!!!!!!!! \beta * [ S_{13}(\beta/2) + S_{12}(\beta/2) ]
         if( NmTDM >= LTrot/2 ) then
            open( 345, file = "Add_Output/S13pS12HfBFinl.txt")
            call MeanErrBar(Bin_Data, Spn13pS12HalfBt(BinStart), Avgr, ErrB)
            write(345, "(A, es17.8, A, es17.8)") char(9), Avgr, char(9), ErrB
            close(345)
         end if
      end if
!________________________________________________________________________________________     
!_________________ (1) IfPyObsPBC==T --> k-point Dynamic correlations ___________________
!_____________________ IfPyObsPBC==F --> r-space Dynamic correlations ___________________
!________________________________________________________________________________________
      open(341, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "GFUpUpAvgErr_Tau.txt")
      open(342, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "GFDnDnAvgErr_Tau.txt")
      open(343, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "SpinZZAvgErr_Tau.txt")
      open(344, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "SpinPMAvgErr_Tau.txt")
      open(345, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "DenDenAvgErr_Tau.txt")
      open(346, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "CurrntAvgErr_Tau.txt")
      open(347, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "PairStAvgErr_Tau.txt")
      open(348, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "EdSParAvgErr_Tau.txt")
      open(349, file = "Add_Output/08_" // merge("K", "R", IfPyObsPBC) // "DWvParAvgErr_Tau.txt")
      do NTInd = 0, NumTauPnt, +1
         !!!!!!!!!! The dynamic single-particle Green's functions
         do I0 = 1, 2, +1
            write(340+I0, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
            do I1 = 1, merge(2, NmSitePairTau, IfPyObsPBC), +1
               call MeanErrBar(Bin_Data, RorKspCrFTauAll(BinStart, NTInd, I1, 2*I0-1), Avgr, ErrB)
               write(340+I0, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, RorKspCrFTauAll(BinStart, NTInd, I1, 2*I0  ), Avgr, ErrB)
               write(340+I0, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            enddo
         enddo
         !!!!!!!!!! The dynamic correlation functions in bosonic channels
         do I0 = 3, 9, +1
            write(340+I0, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
            do I1 = 1, merge(3, NmSitePairTau, IfPyObsPBC), +1
               call MeanErrBar(Bin_Data, RorKspCrFTauAll(BinStart, NTInd, I1, 2+I0), Avgr, ErrB)
               write(340+I0, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            enddo
         enddo
         !!!!!!!!!! Write to a new line
         do I0 = 1, 9, +1
            write(340+I0, "()")
         enddo
      enddo
      close(341); close(342); close(343); close(344); close(345); close(346); close(347); close(348); close(349)
!________________________________________________________________________________________     
!_________________ (2) Store the G(k, Beta/2) for fermi surface _________________________
!_____________________ Store Current(k=0, Beta/2) as dc conductivity ____________________
!_____________________ Only for IfPyObsPBC == T case ____________________________________
!________________________________________________________________________________________
      if( (IfPyObsPBC) .and. (NmTDM >= LTrot/2) ) then
         !!!!!!!!!! Approximated spectral function A(w=e_F, k) = \beta * G(k, Beta/2)
         open(311, file = Trim("Add_Output/GkHalfBetaFinl.txt"))
         do Nk = 1, NumNC, +1
            write(311, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            call MeanErrBar(Bin_Data, GkTauBetaOv2All(BinStart, Nk, 1), Avgr, ErrB)
            if(Avgr < 0.0_rp) Avgr = rp_Eps
            write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            call MeanErrBar(Bin_Data, GkTauBetaOv2All(BinStart, Nk, 2), Avgr, ErrB)
            if(Avgr < 0.0_rp) Avgr = rp_Eps
            write(311, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            write(311, "()")
         enddo
         close(311)
         !!!!!!!!!! Approximated dc conductivity \sigma_{dc} = \beta^2/\pi * Current(k=0, Beta/2)
         open( 311, file = "Add_Output/HdcConductFinl.txt") 
         call MeanErrBar(Bin_Data, Dc_Conductivity(BinStart), Avgr, ErrB)
         write(311, "(A, es17.8, A, es17.8)") char(9), Avgr, char(9), ErrB
         close(311)
      end if
!**************************************************************************************************   
!___________________ 1. For IfPyObsPBC==T case, Compute Matsubara Frequency quantities ____________
!**************************************************************************************************
      if( (IfPyObsPBC) .and. (NmTDMType == 0 .or. NmTDMType == 1) ) then
!________________________________________________________________________________________     
!_________________ (0) For fermion channel --> Green's function and self-energy _________
!________________________________________________________________________________________
         open(313, file = "Add_Output/TT_GrFCrFIwnDyData/GrnFctkIwnFinl.txt")
         open(314, file = "Add_Output/TT_GrFCrFIwnDyData/SlfEngkIwnFinl.txt")
         open(315, file = "Add_Output/TT_GrFCrFIwnDyData/QusParWghtFinl.txt")
         do Nk = 1, NumNC, +1
            write(313, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            write(314, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            write(315, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
            !!!!!!!! G(k,iw_n) and \Sigma(k,iw_n) in spin-up channel
            do FrqInd = 1, NmFrqFermi, +1
               !!!!!! G(k,iw_n)
               call MeanErrBar(Bin_Data, FermiGrF_Iwn(BinStart, Nk, 2*FrqInd-1), Avgr, ErrB)
               write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, FermiGrF_Iwn(BinStart, Nk, 2*FrqInd  ), Avgr, ErrB)
               write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               !!!!!! \Sigma(k,iw_n)
               call MeanErrBar(Bin_Data, SelfEnrg_Iwn(BinStart, Nk, 2*FrqInd-1), Avgr, ErrB)
               write(314, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, SelfEnrg_Iwn(BinStart, Nk, 2*FrqInd  ), Avgr, ErrB)
               write(314, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            enddo
            write(313, "(A)", advance = "no") char(9)
            write(314, "(A)", advance = "no") char(9)
            !!!!!!!! G(k,iw_n) and \Sigma(k,iw_n) in spin-down channel
            do FrqInd = 1, NmFrqFermi, +1
               !!!!!! G(k,iw_n)
               call MeanErrBar(Bin_Data, FermiGrF_Iwn(BinStart, Nk, 2*NmFrqFermi+2*FrqInd-1), Avgr, ErrB)
               write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, FermiGrF_Iwn(BinStart, Nk, 2*NmFrqFermi+2*FrqInd  ), Avgr, ErrB)
               write(313, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               !!!!!! \Sigma(k,iw_n)
               call MeanErrBar(Bin_Data, SelfEnrg_Iwn(BinStart, Nk, 2*NmFrqFermi+2*FrqInd-1), Avgr, ErrB)
               write(314, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               call MeanErrBar(Bin_Data, SelfEnrg_Iwn(BinStart, Nk, 2*NmFrqFermi+2*FrqInd  ), Avgr, ErrB)
               write(314, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            enddo
            !!!!!!!! The quasi-particle weights
            call MeanErrBar(Bin_Data, QuasiParWhgt(BinStart, Nk), Avgr, ErrB)
            write(315, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
            !!!!!!!! New line for output
            write(313, "()")
            write(314, "()")
            write(315, "()")
         enddo
         close(313); close(314); close(315)
!________________________________________________________________________________________     
!_________________ (1) For bosonic channels --> Correlation in various channels _________
!________________________________________________________________________________________
         !!!!!!!!!! S(k, iw_m) for various bosonic channels
         open(311, file = "Add_Output/TT_GrFCrFIwnDyData/SpinZZkIwnFinl.txt")
         open(312, file = "Add_Output/TT_GrFCrFIwnDyData/SpinPMkIwnFinl.txt")
         open(313, file = "Add_Output/TT_GrFCrFIwnDyData/DenDenkIwnFinl.txt")
         open(314, file = "Add_Output/TT_GrFCrFIwnDyData/CurrntkIwnFinl.txt")
         open(315, file = "Add_Output/TT_GrFCrFIwnDyData/PairStkIwnFinl.txt")
         open(316, file = "Add_Output/TT_GrFCrFIwnDyData/EdSParkIwnFinl.txt")
         open(317, file = "Add_Output/TT_GrFCrFIwnDyData/DWvParkIwnFinl.txt")
         do Nk = 1, NumNC, +1
            !!!!!!!! For the standard correlation
            do I0 = 1, 7, +1
               write(310+I0, "(I4, A, I4, A)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9)
               do FrqInd = 1, NmFrqBoson, +1
                  call MeanErrBar(Bin_Data, BosonCrF_Iwn(BinStart, Nk, FrqInd, I0), Avgr, ErrB)
                  write(310+I0, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               enddo
            enddo
            !!!!!!!! For the Vertex Contribution
            do I0 = 1, 7, +1
               write(310+I0, "(A)", advance = "no") char(9)
               do FrqInd = 1, NmFrqBoson, +1
                  call MeanErrBar(Bin_Data, BosonCrF_Iwn(BinStart, Nk, FrqInd, 15+I0), Avgr, ErrB)
                  write(310+I0, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
               enddo
            enddo
            !!!!!!!! Write to a new line
            do I0 = 1, 7, +1
               write(310+I0, "()")
            enddo
         enddo
         close(311); close(312); close(313); close(314); close(315); close(316); close(317)
         !!!!!!!!!! Drude conductivity and superfluid density
         open( 299, file = "Add_Output/TT_GrFCrFIwnDyData/DrudeSpfldFinl.txt")
         !!!!!!!! For the standard correlation
         call MeanErrBar(Bin_Data, DrudeSupfdWt(BinStart, 1), Avgr, ErrB)
         write(299, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, DrudeSupfdWt(BinStart, 2), Avgr, ErrB)
         write(299, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, DrudeSupfdWt(BinStart, 3), Avgr, ErrB)
         write(299, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         !!!!!!!! For the Vertex Contribution
         write(299, "(A)", advance = "no") char(9)
         call MeanErrBar(Bin_Data, DrudeSupfdWt(BinStart, 4), Avgr, ErrB)
         write(299, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, DrudeSupfdWt(BinStart, 5), Avgr, ErrB)
         write(299, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         call MeanErrBar(Bin_Data, DrudeSupfdWt(BinStart, 6), Avgr, ErrB)
         write(299, "(A, es17.8, A, es17.8)", advance = "no") char(9), Avgr, char(9), ErrB
         close(299)
      end if

   end subroutine FinlRKDynCrFTauProc
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!################################## For the simple statistics ###########################################################
!################################## For the simple statistics ###########################################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine MeanErrBar(NDim, VecA, Avgr, ErrB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  MeanErrBar(NDim, VecA, Avgr, ErrB)
! TYPE:     subroutine
! PURPOSE:  This Subroutine calculates the average and errorbar value for an array of data VecA(1:N).
! KEYWORDS: Mean and ErrorBar.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-10
! DESCRIPTION: Calculate the mean and errorbar.
!
!     Input:  NDim --> Length of the VecA vector;
!             VecA --> The data vector;
!
!     Output: Avgr --> The mean value of VecA(StartIndex:NDim).
! 
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!______________________________________________________________________________________________________________   
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
        use RealPrecsn
      use CoreParamt
        implicit none
!______________________________________________________________________________________________________________   
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
        integer NDim            ! Ending integer index for data processing
        real(rp) VecA(NDim)     ! Data vector
        real(rp) Avgr           ! Average value
        real(rp) ErrB           ! Error Bar
!______________________________________________________________________________________________________________   
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
        integer I1
      integer BgnInd
      integer EndInd
      integer NmData
!______________________________________________________________________________________________________________   
!_____________________________ Main calculations Mean values and error bar ____________________________________
!______________________________________________________________________________________________________________
      if(IfCutSmLg) then
!**************************************************************************************************   
!_________________ 0. Sort the vector VecA(1:N) by Ascending order ________________________________
!**************************************************************************************************
         call QuckSortR(1, NDim, VecA(1))
         BgnInd = 2
         EndInd = NDim - 1
         NmData = NDim - 2
      else
         BgnInd = 1
         EndInd = NDim
         NmData = NDim
      end if
!**************************************************************************************************   
!_________________ 1. Calculate the average and error bar for VecA(2:NDim-1) ______________________
!**************************************************************************************************      
      Avgr = 0.0_rp
      do I1 = BgnInd, EndInd, +1
         Avgr = Avgr + VecA(I1)
      enddo
      Avgr = Avgr / dble(NmData)
         
      ErrB = 0.0_rp
      do I1 = BgnInd, EndInd, +1
         ErrB = ErrB + (VecA(I1) - Avgr)**2
      enddo
      ErrB = ErrB / dble(NmData) / dble(NmData-1)
      ErrB = sqrt(ErrB)

   end subroutine MeanErrBar
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
    subroutine MeanErrBar_Var(NDim, VecA, Avgr, ErrB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  MeanErrBar_Var(NDim, VecA, Avgr, ErrB)
! TYPE:     subroutine
! PURPOSE:  This Subroutine calculates the average and errorbar value for an array of data VecA(1:N).
! KEYWORDS: Mean and ErrorBar.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate the mean and errorbar.
!
!     Input:  NDim --> Length of the VecA vector;
!             VecA --> The data vector;
!
!     Output: Avgr --> The mean value of VecA(StartIndex:NDim).
! 
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!______________________________________________________________________________________________________________   
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
        use RealPrecsn
      use CoreParamt
        implicit none
!______________________________________________________________________________________________________________   
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
        integer NDim            ! Ending integer index for data processing
        real(rp) VecA(NDim)     ! Data vector
        real(rp) Avgr           ! Average value
        real(rp) ErrB           ! Error Bar
!______________________________________________________________________________________________________________   
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
        integer I1
      integer BgnInd
      integer EndInd
      integer NmData
!______________________________________________________________________________________________________________   
!_____________________________ Main calculations Mean values and error bar ____________________________________
!______________________________________________________________________________________________________________
      if(IfCutSmLg) then
!**************************************************************************************************   
!_________________ 0. Sort the vector VecA(1:N) by Ascending order ________________________________
!**************************************************************************************************
         call QuckSortR(1, NDim, VecA(1))
         BgnInd = 2
         EndInd = NDim - 1
         NmData = NDim - 2
      else
         BgnInd = 1
         EndInd = NDim
         NmData = NDim
      end if
!**************************************************************************************************   
!_________________ 1. Calculate the average and error bar for VecA(2:NDim-1) ______________________
!**************************************************************************************************      
      Avgr = 0.0_rp
      do I1 = BgnInd, EndInd, +1
         Avgr = Avgr + VecA(I1)
      enddo
      Avgr = Avgr / dble(NmData)
         
      ErrB = 0.0_rp
      do I1 = BgnInd, EndInd, +1
         ErrB = ErrB + (VecA(I1) - Avgr)**2
      enddo
      ErrB = ErrB / dble(NmData)

   end subroutine MeanErrBar_Var
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$