!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL SPIN PINNING FIELDS; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the data process in every BIN simulation after the simulation process
!              in this BIN. This data process only calculates the average values for static observables.
! COMMENT: Post process for static observables in every BIN simulation.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   PostStatic   --> Subroutine to calculate average values for the static observables;
!
!   ProcEngLcDen --> Subroutine to process static data of the weight, energies and occupations;
!   ProcRKStaCrF --> Subroutine to process static data of r-space and k-space correlation functions;
!
!   DeductCrFBkgSta --> Subroutine to deduct the background for the correlation functions;
!   VertexContrbSta --> Subroutine to compute the vertex contribution for all the correlation functions;
!   FourierTransSta --> Subroutine to perform Fourier Transformation for correlation functions.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine PostStatic(NB, MeaInd) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PostStatic(NB, MeaInd) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the post process for the CPMC calculated values of observables, 
!                namely calculates the average values for the observables in a single BIN for NSwep calculations.
! KEYWORDS: Post Process of observables.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Post Process of static observables.
!
!     Input: NB     --> The iteration number of NmBin for the CPMC simulation;  
!            MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
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
      use TimeRecord
      use QMCTimeRec
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB, MeaInd
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in data post process ______________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&         
      TimsDatPr = TimsDatPr + 1
      call system_clock(time1)
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for Average values of observables ____________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. BIN process of all STATIC physical observables ____________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Results of the energies and local densities ______________________
!________________________________________________________________________________________
      call ProcEngLcDen(NB, MeaInd)
!________________________________________________________________________________________      
!_________________ (1) Results of momentum distribution and pairing wavefunction ________
!________________________________________________________________________________________
      if(IfFftEnPar) then
         call ProcNkPairWf(NB, MeaInd)
      end if
!________________________________________________________________________________________      
!_________________ (2) Results of r-space correlations for PERIODIC and OPEN BCs ________
!________________________________________________________________________________________
      if(abs(PinSz) < rp_Eps) then
         call ProcRKStaCrF(NB, MeaInd)
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for data process __________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&     
      call system_clock(time2)
      TimeDatPr = TimeDatPr + TimeIntrvl(time1, time2)

   end subroutine PostStatic
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!############################### Process results of energies, occupations and local densities ###########################
!############################### Process results of energies, occupations and local densities ###########################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ProcEngLcDen(NB, MeaInd) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ProcEngLcDen(NB, MeaInd) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform BIN data process for static energies and occupations.
! KEYWORDS: Post Process of observables.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Process energies and density-density correlations.
!
!     Input: NB     --> The iteration number of NmBin for the CPMC simulation;  
!            MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
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
      use TimeRecord
      use QMCTimeRec
      use MPISetting
      use StdInOutSt
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB, MeaInd
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, Ix, Iy, SiteParity
      complex(rp) Ztp0
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for Average values of observables ____________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Results of the energies and occupations ___________________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) The average of all results for the present BIN ___________________
!________________________________________________________________________________________
      EngOccCrFBIN(:, MeaInd) = EngOccCrFBIN(:, MeaInd) / WtMeanSumBIN
      if(MeaInd == 0) then
         if( (mod(amyid, merge(anprc/4, 1, anprc>4)) == 0) .or. (amyid == anprc-1) ) then
            write(*, "(28x, 'PostStatic: ', I4.4, 2x, I4.4, 3es20.11)") NB, amyid, EngOccCrFBIN(31:32, MeaInd), &
               & EngOccCrFBIN(17, MeaInd)
         end if
         if(amyid == amstr) then
            open( 999, file = Trim(FMnt), access = "append")
            write(999, "('PostStatic: ', I4.4, 2x, I4.4, 3es20.11)") NB, amyid, EngOccCrFBIN(31:32, MeaInd), &
               & EngOccCrFBIN(17, MeaInd)
            write(999, "()")
            write(999, "()")
            close(999)
         end if
      end if
!________________________________________________________________________________________      
!_________________ (1) Store and output the results for the present BIN _________________
!________________________________________________________________________________________
      if(amyid == amstr) then 
!____________________________________________________________________________      
!________________ [0] The energies of all terms in Hamiltonian ______________
!____________________________________________________________________________  
         !!!!!!!!!! Record results of the present BIN   
         EnHopt1(NB, MeaInd) = EngOccCrFBIN(01, MeaInd)
         EnHopt2(NB, MeaInd) = EngOccCrFBIN(02, MeaInd)
         EnHopt3(NB, MeaInd) = EngOccCrFBIN(03, MeaInd)
         EnZmFld(NB, MeaInd) = EngOccCrFBIN(04, MeaInd)
         EnPinSz(NB, MeaInd) = EngOccCrFBIN(05, MeaInd)
         EnDopCh(NB, MeaInd) = EngOccCrFBIN(06, MeaInd)
         EnHubbU(NB, MeaInd) = EngOccCrFBIN(07, MeaInd)
         EnTotal(NB, MeaInd) = EngOccCrFBIN(08, MeaInd)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         EnSinusoidalPinSz(NB, MeaInd) = EngOccCrFBIN(22, MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         !!!!!!!!!! Output results of the present BIN 
         open( 291, file = Trim(FileOutAdd(MeaInd)) // "01_QMCExpectEnergy" // Trim(FileAddTxt(MeaInd)), access = "append") 
         write(291, "(I4.4)", advance = "no") NB
         do I0 = 1, 8, +1
            write(291, "(A, es25.16)", advance = "no") char(9), EngOccCrFBIN(I0, MeaInd)
         enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         write(291, "(A, es25.16)", advance = "no") char(9), EngOccCrFBIN(22, MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         write(291, "()")
         close(291)
!____________________________________________________________________________      
!________________ [1] First-order derivative of model parameters ____________
!____________________________________________________________________________
         !!!!!!!!!! Record results of the present BIN
         EnHbUCh(NB, MeaInd) = EngOccCrFBIN(16, MeaInd)
         EnTotCh(NB, MeaInd) = EngOccCrFBIN(17, MeaInd)
         HmtOvt2(NB, MeaInd) = EngOccCrFBIN(18, MeaInd)
         HmtOvt3(NB, MeaInd) = EngOccCrFBIN(19, MeaInd)
         HmtOvZm(NB, MeaInd) = EngOccCrFBIN(20, MeaInd)
         HmtOvbU(NB, MeaInd) = EngOccCrFBIN(21, MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         HmOvSinusoidalPinSz(NB, MeaInd) = EngOccCrFBIN(23, MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !!!!!!!!!! Output results of the present BIN 
         open( 293, file = Trim(FileOutAdd(MeaInd)) // "01_FstOrdDiffParam" // Trim(FileAddTxt(MeaInd)), access = "append")
         write(293, "(I4.4)", advance = "no") NB
         do I0 = 16, 17, +1
            write(293, "(A, es25.16)", advance = "no") char(9), EngOccCrFBIN(I0, MeaInd)
         enddo
         write(293, "(A)", advance = "no") char(9)
         do I0 = 18, 21, +1
            write(293, "(A, es25.16)", advance = "no") char(9), EngOccCrFBIN(I0, MeaInd)
         enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         write(293, "(A, es25.16)", advance = "no") char(9), EngOccCrFBIN(23, MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         write(293, "()")
         close(293)
!____________________________________________________________________________      
!________________ [2] Fillings, DouOcc and density correlations _____________
!____________________________________________________________________________
         !!!!!!!!!! Record results of the present BIN
         nOccpUp(NB, MeaInd) = EngOccCrFBIN(31, MeaInd)
         nOccpDw(NB, MeaInd) = EngOccCrFBIN(32, MeaInd)
         nOccTot(NB, MeaInd) = EngOccCrFBIN(33, MeaInd)
         NmNeTot(NB, MeaInd) = EngOccCrFBIN(34, MeaInd)
         RDouOcc(NB, MeaInd) = EngOccCrFBIN(35, MeaInd)
         RNNDCrF(NB, MeaInd) = EngOccCrFBIN(36, MeaInd)
         !!!!!!!!!! Output results of the present BIN 
         open( 295, file = Trim(FileOutAdd(MeaInd)) // "01_OccDouSpnCorrFc" // Trim(FileAddTxt(MeaInd)), access = "append")
         write(295, "(I4.4)", advance = "no") NB
         do I0 = 31, 34, +1
            write(295, "(A, es25.16)", advance = "no") char(9), EngOccCrFBIN(I0, MeaInd)
         enddo
         write(295, "(A)", advance = "no") char(9)
         do I0 = 35, 36, +1
            write(295, "(A, es25.16)", advance = "no") char(9), EngOccCrFBIN(I0, MeaInd)
         enddo
         write(295, "()")
         close(295)
!____________________________________________________________________________      
!________________ [3] Total energy per particle in unit of E_{FG} ___________
!____________________________________________________________________________
         if(EkDispType >= 1) then
            !!!!!!!!!! Record results of the present BIN
            EpOvEfg(NB, 1, MeaInd) = EnTotCh(NB, MeaInd) * dble(NumNS)/NumNe   * dble(NumNS)/NumNe   / rp_pi
            EpOvEfg(NB, 2, MeaInd) = EnTotCh(NB, MeaInd) / nOccTot(NB, MeaInd) / nOccTot(NB, MeaInd) / rp_pi
            !!!!!!!!!! Output results of the present BIN 
            open( 293, file = Trim(FileOutAdd(MeaInd)) // "01_EnPerFermionEfg" // Trim(FileAddTxt(MeaInd)), access = "append")
            write(293, "(I4, A, es25.16, A, es25.16)") NB, char(9), EpOvEfg(NB, 1, MeaInd), &
               & char(9), EpOvEfg(NB, 2, MeaInd)
            close(293)
         end if
      end if
!**************************************************************************************************     
!___________________ 1. Results of single-particle GrF matrices in r- and k-space _________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) The r-space single-particle GrF matrix ___________________________
!________________________________________________________________________________________
!____________________________________________________________________________      
!________________ [0] Collect the results of the present BIN ________________
!____________________________________________________________________________
      RlSpGrnFtBIN(:, :, :, MeaInd) = RlSpGrnFtBIN(:, :, :, MeaInd) / WtMeanSumBIN
!____________________________________________________________________________      
!________________ [1] Store the results of the present BIN __________________
!____________________________________________________________________________
      if(amyid == amstr) then
         do I2 = 1, NumNC, +1
            do I1 = 1, I2, +1
               I0 = (I2-1)*I2/2 + I1
               RlSpGrnFtAll(NB, I0, 1, MeaInd) = RlSpGrnFtBIN(I1, I2, 1, MeaInd) + RlSpGrnFtBIN(I2, I1, 1, MeaInd)
               RlSpGrnFtAll(NB, I0, 2, MeaInd) = RlSpGrnFtBIN(I1, I2, 2, MeaInd) + RlSpGrnFtBIN(I2, I1, 2, MeaInd)
            enddo
         enddo
         RlSpGrnFtAll(NB, :, :, MeaInd) = RlSpGrnFtAll(NB, :, :, MeaInd) / 2.0_rp
      end if
!________________________________________________________________________________________      
!_________________ (1) The k-space single-particle GrF matrix ___________________________
!________________________________________________________________________________________
      if(IfFftEnPar) then
!____________________________________________________________________________      
!________________ [0] Collect the results of the present BIN ________________
!____________________________________________________________________________
         KSpGreenFBIN(:, :, :, MeaInd) = KSpGreenFBIN(:, :, :, MeaInd) / WtMeanSumBIN
!____________________________________________________________________________      
!________________ [1] Store the results of the present BIN __________________
!____________________________________________________________________________
         if(amyid == amstr) then
            do I2 = 1, NumNC, +1
               do I1 = 1, I2, +1
                  I0 = (I2-1)*I2/2 + I1
                  !!!!!!!! Record the results of the present BIN simulation
                  Ztp0 = ( KSpGreenFBIN(I1, I2, 1, MeaInd) + conjg(KSpGreenFBIN(I2, I1, 1, MeaInd)) ) / 2.0_rp
                  KSpGreenFAll(NB, I0, 1, MeaInd) =  real(Ztp0)
                  KSpGreenFAll(NB, I0, 2, MeaInd) = aimag(Ztp0)
                  Ztp0 = ( KSpGreenFBIN(I1, I2, 2, MeaInd) + conjg(KSpGreenFBIN(I2, I1, 2, MeaInd)) ) / 2.0_rp
                  KSpGreenFAll(NB, I0, 3, MeaInd) =  real(Ztp0)
                  KSpGreenFAll(NB, I0, 4, MeaInd) = aimag(Ztp0)
               enddo
            enddo
         end if
      end if
!**************************************************************************************************     
!___________________ 2. Results of local charge and spin density __________________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Compute densities of the present BIN from RlSpGrnFtBIN ___________
!________________________________________________________________________________________
      do I0 = 1, NumNS, +1
         RSpLclOrdBIN(I0, 1, MeaInd) = RlSpGrnFtBIN(I0, I0, 1, MeaInd)
         RSpLclOrdBIN(I0, 2, MeaInd) = RlSpGrnFtBIN(I0, I0, 2, MeaInd)
      enddo
!________________________________________________________________________________________      
!_________________ (1) Store and output results of every BIN simulation _________________
!________________________________________________________________________________________
      if(amyid == amstr) then
!____________________________________________________________________________      
!________________ [0] Record results of the present BIN _____________________
!____________________________________________________________________________     
         !!!!!!!! Store the results for Average and Errorbar
         do I0 = 1, NumNS, +1
            !!!!!! Local charge density
            RSpLclOrdAll(NB, I0, 1, MeaInd) = RSpLclOrdBIN(I0, 1, MeaInd)
            RSpLclOrdAll(NB, I0, 2, MeaInd) = RSpLclOrdBIN(I0, 2, MeaInd)
            RSpLclOrdAll(NB, I0, 3, MeaInd) = RSpLclOrdBIN(I0, 1, MeaInd) + RSpLclOrdBIN(I0, 2, MeaInd)
            !!!!!! Local spin density
            RSpLclOrdAll(NB, I0, 4, MeaInd) = ( RSpLclOrdBIN(I0, 1, MeaInd) - RSpLclOrdBIN(I0, 2, MeaInd) ) / 2.0_rp
         enddo
         !!!!!!!! The Average local spin moment Sz        
         TotMagMmtAll(NB, MeaInd) = 0.0_rp
         do I0 = 1, NumNS, +1
            Ix = StList(I0, 1); Iy = StList(I0, 2)
            SiteParity = (-1)**(mod(Ix+Iy+1, 2))
            TotMagMmtAll(NB, MeaInd) = TotMagMmtAll(NB, MeaInd) + dble(SiteParity)*RSpLclOrdAll(NB, I0, 4, MeaInd)
         enddo
         TotMagMmtAll(NB, MeaInd) = TotMagMmtAll(NB, MeaInd) / dble(NumNS)
!____________________________________________________________________________      
!________________ [1] Output results of the present BIN _____________________
!____________________________________________________________________________ 
         !!!!!!!! The r-space local charge and spin densities
         open( 301, file = Trim(FileOutAdd(MeaInd)) // "02_RLocalNeDensity" // Trim(FileAddTxt(MeaInd)), access = "append")
         open( 302, file = Trim(FileOutAdd(MeaInd)) // "02_RLocalSpzSpxOrd" // Trim(FileAddTxt(MeaInd)), access = "append")
         write(301, "(I6)") NB
         write(302, "(I6)") NB
         do Iy = 1, NumL2, +1
            do Ix = 1, NumL1, +1
               I0 = (Iy-1)*NumL1 + Ix
               !!!!!! local density for spin-up, down and total
               write(301, "(I6, A, I4, A, I4)", advance = "no") I0, char(9), Ix, char(9), Iy
               write(301, "(A, es25.16, A, es25.16, A, es25.16)", advance = "no") &
                  & char(9), RSpLclOrdAll(NB, I0, 1, MeaInd), char(9), RSpLclOrdAll(NB, I0, 2, MeaInd), &
                  & char(9), RSpLclOrdAll(NB, I0, 3, MeaInd)
               write(301, "()")
               !!!!!! local spin density
               write(302, "(I6, A, I4, A, I4)", advance = "no") I0, char(9), Ix, char(9), Iy
               write(302, "(A, es25.16      )", advance = "no") char(9), RSpLclOrdAll(NB, I0, 4, MeaInd)
               write(302, "()")
            enddo
         enddo
         close(301); close(302)
         !!!!!!!! The Average local spin moment Sz
         if(mod(NumNS, 2) == 0) then
            open( 313, file = Trim(FileOutAdd(MeaInd)) // "02_TotalMagntMomnt" // Trim(FileAddTxt(MeaInd)), access = "append")
            write(313, "(I4, A, es25.16)") NB, char(9), TotMagMmtAll(NB, MeaInd)
            close(313)
         end if
      end if

   end subroutine ProcEngLcDen
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!############################### IfFftEnPar == T --> Process n(k) and pairing wavefunction ##############################
!############################### IfFftEnPar == T --> Process n(k) and pairing wavefunction ##############################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ProcNkPairWf(NB, MeaInd) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ProcNkPairWf(NB, MeaInd) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to process the results of momentum distribution and pairing wavefunction for
!                   IfFftEnPar == T case, after every single BIN simulation.
! KEYWORDS: Post Process of n(k) and pairing wavefunction.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Process pairing wavefunctions and condensation fraction.
!
!     Input: NB     --> The iteration number of NmBin for the CPMC simulation;  
!            MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
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
      use StdInOutSt
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB, MeaInd
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, Ix, Iy, Nk
      integer VecR(2)
      real(rp) Rtp1, VecK(2)
      complex(rp) Ztp1
      real(rp)   , allocatable :: Collect0(:   )
      complex(rp), allocatable :: Collect1(:, :)
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for data process after BIN simulation ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Process the results of momentum distribution n(k) _________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Collect the results of the present BIN ___________________________
!________________________________________________________________________________________
      NkSgleBIN(:, :, MeaInd) = NkSgleBIN(:, :, MeaInd) / WtMeanSumBIN
!________________________________________________________________________________________      
!_________________ (1) Output the average results of the present BIN ____________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         open(291, file = Trim(FileOutAdd(MeaInd)) // "MomtumDistBINs" // Trim(FileAddTxt(MeaInd)), access = "append")
         do Nk = 1, NumNC, +1
            NkDistrib(NB, Nk, 1, MeaInd) = dble(NkSgleBIN(Nk, 1, MeaInd))
            NkDistrib(NB, Nk, 2, MeaInd) = dble(NkSgleBIN(Nk, 2, MeaInd))
            write(291, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
            write(291, "(A, es25.16, A, es25.16)", advance = "no") char(9), NkDistrib(NB, Nk, 1, MeaInd), &
                                                                 & char(9), NkDistrib(NB, Nk, 2, MeaInd)
            write(291, "()")
         enddo
         close(291)
      end if
!**************************************************************************************************     
!___________________ 1. Process the results of pairing wavefunctions ______________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Collect the results of the present BIN ___________________________
!________________________________________________________________________________________
      !!!!!!!!!! Average of the measurements
      PairMtBIN(:, :, MeaInd) = PairMtBIN(:, :, MeaInd) / WtMeanSumBIN
      !!!!!!!!!! Deduct the diagonal background term for PairMtBIN
      if(HubbU <= 0.0_rp) then
         do I0 = 1, NumNC, +1
            Ix = - mod(I0-1, NumL1); Iy = - (I0-1)/NumL1
            if(Ix < 0) Ix = Ix + NumL1
            if(Iy < 0) Iy = Iy + NumL2
            I1 = Iy*NumL1 + Ix + 1
            PairMtBIN(I0, I0, MeaInd) = PairMtBIN(I0, I0, MeaInd) - &
                              & cmplx(real(NkSgleBIN(I0, 1, MeaInd))*real(NkSgleBIN(I1, 2, MeaInd)), 0.0_rp, rp)
         enddo
      end if
!________________________________________________________________________________________      
!_________________ (1) Diagonalize PairMtBIN matrix and store results ___________________
!________________________________________________________________________________________
      if(amyid == amstr) then
!____________________________________________________________________________      
!________________ [0] Diagonalize PairMtBIN and store results _______________
!____________________________________________________________________________
         !!!!!!!!!! Make PairMtBIN matrix Hermitian
         PairMtBIN(:, :, MeaInd) = ( PairMtBIN(:, :, MeaInd) + conjg(transpose(PairMtBIN(:, :, MeaInd))) ) / 2.0_rp
         !!!!!!!!!! Diagonalize the PairMtBIN matrix
         allocate(Collect0(NumNC), Collect1(NumNC, NumNC))
         Collect0 = 0.0_rp; Collect1 = rp_Zzero
         call MatrDiagZ2(NumNC, NumNC, PairMtBIN(1, 1, MeaInd), Collect0(1), Collect1(1, 1))
         !!!!!!!!!! Obtain the condensation fraction
         Condfract(NB, 1, MeaInd) = Collect0(NumNC) * 2.0_rp / NumNe
         Condfract(NB, 2, MeaInd) = Collect0(NumNC) * 2.0_rp / NmNeTot(NB, MeaInd)
         !!!!!!!!!! Deal with the global phase of wavefunction
         if( real(Collect1(1, NumNC)) < 0.0_rp ) then
            do I1 = 1, NumNC, +1
               Collect1(I1, NumNC) = - Collect1(I1, NumNC)
            enddo
         end if
         !!!!!!!!!! The k space wavefunction for pairing
         do I0 = 1, NumNC, +1
            PairWvfct(NB, I0, 1, MeaInd) =  dble( Collect1(I0, NumNC) )
            PairWvfct(NB, I0, 2, MeaInd) = aimag( Collect1(I0, NumNC) )
            PairWvfct(NB, I0, 3, MeaInd) =  sqrt( dble( Collect1(I0, NumNC)*conjg(Collect1(I0, NumNC)) ) )
         enddo
         !!!!!!!!!! The r space wavefunction for pairing
         do I0 = 1, NumNC, +1
            Ix = mod(I0-1, NumL1); Iy = (I0-1)/NumL1
            VecR(1) = -NumL1/2 + Ix
            VecR(2) = -NumL2/2 + Iy
            Ztp1 = rp_Zzero
            do I1 = 1, NumNC, +1
               Ix = mod(I1-1, NumL1); Iy = (I1-1)/NumL1
               VecK(1) = dble(Ix) * 2.0_rp * rp_pi / dble(NumL1)
               VecK(2) = dble(Iy) * 2.0_rp * rp_pi / dble(NumL2)
               Rtp1 = dot_product(VecK, dble(VecR))
               Ztp1 = Ztp1 + exp(cmplx(0.0_rp, -Rtp1, rp))*Collect1(I1, NumNC)
            enddo
            Ztp1 = Ztp1 / dble(NumNC)
            PairWvfct(NB, I0, 4, MeaInd) =  dble( Ztp1 )
            PairWvfct(NB, I0, 5, MeaInd) = aimag( Ztp1 )
            PairWvfct(NB, I0, 6, MeaInd) =  sqrt( dble(Ztp1*conjg(Ztp1)) )
         enddo
!____________________________________________________________________________      
!________________ [1] Diagonalize PairMtAcm and store results _______________
!____________________________________________________________________________
         !!!!!!!!!! Accumulate the pairing matrix
         PairMtAcm(:, :, MeaInd) = PairMtAcm(:, :, MeaInd) + PairMtBIN(:, :, MeaInd)
         PairMtBIN(:, :, MeaInd) = PairMtAcm(:, :, MeaInd) / dble(NB)
         !!!!!!!!!! Diagonalize the PairMtAcm matrix
         Collect0 = 0.0_rp; Collect1 = rp_Zzero
         if(NB == NmBin) then
            call MatrDiagZ2(NumNC, NumNC, PairMtBIN(1, 1, MeaInd), Collect0(1), Collect1(1, 1))
         else
            call MatrDiagZ1(NumNC, NumNC, PairMtBIN(1, 1, MeaInd), Collect0(1))
         end if
         !!!!!!!!!! Obtain the condensation fraction
         Condfract(NmBin+1, 1, MeaInd) = Collect0(NumNC) * 2.0_rp / NumNe
         Condfract(NmBin+1, 2, MeaInd) = Collect0(NumNC) * 2.0_rp / (sum(NmNeTot(1:NB, MeaInd))/dble(NB))
         !!!!!!!!!! For NB == NmBin, obtain the pairing wavefunction
         if(NB == NmBin) then
            !!!!!!!! Fix the gauge
            if( real(Collect1(1, NumNC)) < 0.0_rp ) then
               do I1 = 1, NumNC, +1
                  Collect1(I1, NumNC) = - Collect1(I1, NumNC)
               enddo
            end if
            !!!!!!!! The pairing wavefunction in k space
            do I0 = 1, NumNC, +1
               PairWvfct(NmBin+1, I0, 1, MeaInd) =  dble( Collect1(I0, NumNC) )
               PairWvfct(NmBin+1, I0, 2, MeaInd) = aimag( Collect1(I0, NumNC) )
               PairWvfct(NmBin+1, I0, 3, MeaInd) =  sqrt( dble( Collect1(I0, NumNC)*conjg(Collect1(I0, NumNC)) ) )
            enddo
            !!!!!!!! The pairing wavefunction in r space
            do I0 = 1, NumNC, +1
               Ix = mod(I0-1, NumL1); Iy = (I0-1)/NumL1
               VecR(1) = -NumL1/2 + Ix
               VecR(2) = -NumL2/2 + Iy
               Ztp1 = rp_Zzero
               do I1 = 1, NumNC, +1
                  Ix = mod(I1-1, NumL1); Iy = (I1-1)/NumL1
                  VecK(1) = dble(Ix) * 2.0_rp * rp_pi / dble(NumL1)
                  VecK(2) = dble(Iy) * 2.0_rp * rp_pi / dble(NumL2)
                  Rtp1 = dot_product(VecK, dble(VecR))
                  Ztp1 = Ztp1 + exp(cmplx(0.0_rp, -Rtp1, rp))*Collect1(I1, NumNC)
               enddo
               Ztp1 = Ztp1 / dble(NumNC)
               PairWvfct(NmBin+1, I0, 4, MeaInd) =  dble( Ztp1 )
               PairWvfct(NmBin+1, I0, 5, MeaInd) = aimag( Ztp1 )
               PairWvfct(NmBin+1, I0, 6, MeaInd) =  sqrt( dble(Ztp1*conjg(Ztp1)) )
            enddo
         end if
         !!!!!!!!!! Deallocate the matrices used above
         if(allocated(Collect0)) deallocate(Collect0)
         if(allocated(Collect1)) deallocate(Collect1)
      end if
!________________________________________________________________________________________      
!_________________ (2) Output the average results of the present BIN ____________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         !!!!!!!! The condensate fraction
         open( 289, file = Trim(FileOutAdd(MeaInd)) // "CondstFracBINs" // Trim(FileAddTxt(MeaInd)), access = "append")
         write(289, "(I4, A, es25.16, A, es25.16, A, es25.16, A, es25.16)") NB, char(9), Condfract(NB, 1, MeaInd), &
            & char(9), Condfract(NB, 2, MeaInd), char(9), Condfract(NmBin+1, 1, MeaInd), &
            & char(9), Condfract(NmBin+1, 2, MeaInd)
         close(289)
         !!!!!!!! The pairing wavefunction in k and r spaces
         open(291, file = Trim(FileOutAdd(MeaInd)) // "PairWvfcMtBINs" // Trim(FileAddTxt(MeaInd)), access = "append")
         open(293, file = Trim(FileOutAdd(MeaInd)) // "PairWvfcRlBINs" // Trim(FileAddTxt(MeaInd)), access = "append")
         do I0 = 1, NumNC, +1
            !!!!!! The k space pairing wavefunction
            write(291, "(I4, A, I4, A, es25.16, A, es25.16, A, es25.16)") KpList(I0, 1), char(9), KpList(I0, 2), &
               & char(9), PairWvfct(NB, I0, 1, MeaInd), char(9), PairWvfct(NB, I0, 2, MeaInd), &
               & char(9), PairWvfct(NB, I0, 3, MeaInd)
            !!!!!! The r space pairing wavefunction
            Ix = mod(I0-1, NumL1); Iy = (I0-1)/NumL1
            write(293, "(I4, A, I4, A, es25.16, A, es25.16, A, es25.16)") Ix-NumL1/2, char(9), Iy-NumL2/2, &
               & char(9), PairWvfct(NB, I0, 4, MeaInd), char(9), PairWvfct(NB, I0, 5, MeaInd), &
               & char(9), PairWvfct(NB, I0, 6, MeaInd)
         enddo
         close(291); close(293)
      end if

   end subroutine ProcNkPairWf
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!############################ Process static data of r- and k-space correlation functions ###############################
!############################ Process static data of r- and k-space correlation functions ###############################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ProcRKStaCrF(NB, MeaInd) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ProcRKStaCrF(NB, MeaInd) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform BIN data process for r-space correlation functions for PERIODIC
!                 boundary conditions.
! KEYWORDS: Post Process of observables.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Post Process of r-space correlation functions.
!
!     Input: NB     --> The iteration number of NmBin for the CPMC simulation;  
!            MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
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
      use TimeRecord
      use QMCTimeRec
      use MPISetting
      use StdInOutSt
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB           ! Number index of the BIN 
      integer MeaInd       ! Index for measure kinds
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, Idimj, Nk
      integer KPntIndx(10)
      real(rp) RlCrFTmp(15)
      complex(rp) ZTempA(60), ZTempB(60)
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for Average values of observables ____________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Results of the r-space correlation functions ______________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Collect the results of the present BIN ___________________________
!________________________________________________________________________________________
      RealSpCrFBIN(:, :, MeaInd) = RealSpCrFBIN(:, :, MeaInd) / WtMeanSumBIN
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
!_________________ (1) The Vertex Contribution for pairing correlations _________________
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
      if(amyid == amstr) then
         call VertexContrbSta(MeaInd)
      end if
!________________________________________________________________________________________      
!_________________ (2) Deduct the background for correlation functions __________________
!________________________________________________________________________________________
      if( (amyid == amstr) .and. (IfCrfDfBkg) ) then
         call DeductCrFBkgSta(MeaInd)
      end if
!**************************************************************************************************     
!___________________ 1. Postprocess of r- and k-space correlations for both _______________________
!_____________________________ PERIODIC and OPEN boundary conditions ______________________________
!**************************************************************************************************
      if(amyid == amstr) then
!________________________________________________________________________________________      
!_________________ (0) Store results of r-space correlations for the presnt BIN _________
!_____________________ for both IfPyObsPBC==T and IfPyObsPBC==F cases ___________________
!________________________________________________________________________________________
         do I0 = 1, 40, +1
            do Idimj = 1, NmSitePair, +1
               RealSpCrFAll(NB, Idimj, I0, MeaInd) = RealSpCrFBIN(Idimj, I0, MeaInd)
            enddo
         enddo
!________________________________________________________________________________________      
!_________________ (1) Compute, store k-space structure factors for the presnt BIN ______
!_____________________ only for the IfPyObsPBC==T case __________________________________
!_____________________ From RealSpCrFBIN to KSpaceCrFAll ________________________________
!________________________________________________________________________________________
         if(IfPyObsPBC) then
            call FourierTransSta(NB, MeaInd)
         end if
!________________________________________________________________________________________      
!_________________ (2) Obtain structure factors and correlation ratios at _______________
!_____________________ special k points only for the IfPyObsPBC==T case _________________
!________________________________________________________________________________________
         if(IfPyObsPBC) then
!____________________________________________________________________________      
!________________ [0] The structure factors --> Consider symmtry ____________
!____________________________________________________________________________
            !!!!!!!!!! Initialization for the results
            ZTempA = cmplx(1.0E+64_rp, 0.0_rp, rp)
            !!!!!!!!!! For Gamma=(0,0) Point
            KPntIndx(1) = InvKpList(0, 0)
            ZTempA(01:16) = KSpaceCrFAll(NB, KPntIndx(1), 01:16, MeaInd)
            !!!!!!!!!! For X1=(pi,0) and X2=(0,pi) Points
            KPntIndx(2) = InvKpList(+NumL1/2,        0)
            KPntIndx(3) = InvKpList(       0, +NumL2/2)
            do I1 = 17, 32, +1
               ZTempA(I1) = rp_Zzero
               do I0 = 2, merge(3, 2, NumL1==NumL2), +1
                  ZTempA(I1) = ZTempA(I1) + KSpaceCrFAll(NB, KPntIndx(I0), I1-16, MeaInd)
               enddo
               ZTempA(I1) = ZTempA(I1) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
            enddo
            !!!!!!!!!! For M=(pi/2,pi/2) Point
            KPntIndx(4) = InvKpList(+NumL1/2, +NumL2/2)
            ZTempA(33:48) = KSpaceCrFAll(NB, KPntIndx(4), 01:16, MeaInd)
!____________________________________________________________________________      
!________________ [1] The correlation ratios --> Consider symmtry ___________
!____________________________________________________________________________
            !!!!!!!!!! Initialization for the results
            ZTempB = rp_Zzero
            !!!!!!!!!! For Gamma=(0,0) Point
            KPntIndx(1) = InvKpList(+1,  0)
            KPntIndx(2) = InvKpList( 0, +1)
            do I1 = 01, 16, +1
               ZTempB(I1) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  ZTempB(I1) = ZTempB(I1) + KSpaceCrFAll(NB, KPntIndx(I0), I1, MeaInd)
               enddo
               ZTempB(I1) = ZTempB(I1) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               if(abs(ZTempA(I1)) > rp_Eps) then
                  ZTempB(I1) = rp_Z_One - ZTempB(I1) / ZTempA(I1)
               else
                  ZTempB(I1) = rp_Zzero
               end if
            enddo
            !!!!!!!!!! For X1=(pi,0) and X2=(0,pi) Points
            KPntIndx(5) = InvKpList(+NumL1/2-1,          0)
            KPntIndx(6) = InvKpList(+NumL1/2  ,         +1)
            KPntIndx(7) = InvKpList(         0, +NumL2/2-1)
            KPntIndx(8) = InvKpList(        +1, +NumL2/2  )
            do I1 = 17, 32, +1
               ZTempB(I1) = rp_Zzero
               do I0 = 5, merge(8, 5, NumL1==NumL2), +1
                  ZTempB(I1) = ZTempB(I1) + KSpaceCrFAll(NB, KPntIndx(I0), I1-16, MeaInd)
               enddo
               ZTempB(I1) = ZTempB(I1) / merge(4.0_rp, 1.0_rp, NumL1==NumL2)
               if(abs(ZTempA(I1)) > rp_Eps) then
                  ZTempB(I1) = rp_Z_One - ZTempB(I1) / ZTempA(I1)
               else
                  ZTempB(I1) = rp_Zzero
               end if
            enddo
            !!!!!!!!!! For M=(pi,pi) Point
            KPntIndx(09) = InvKpList(+NumL1/2-1, +NumL2/2  )
            KPntIndx(10) = InvKpList(+NumL1/2  , +NumL2/2-1) 
            do I1 = 33, 48, +1
               ZTempB(I1) = rp_Zzero
               do I0 = 09, merge(10, 09, NumL1==NumL2), +1
                  ZTempB(I1) = ZTempB(I1) + KSpaceCrFAll(NB, KPntIndx(I0), I1-32, MeaInd)
               enddo
               ZTempB(I1) = ZTempB(I1) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               if(abs(ZTempA(I1)) > rp_Eps) then
                  ZTempB(I1) = rp_Z_One - ZTempB(I1) / ZTempA(I1)
               else
                  ZTempB(I1) = rp_Zzero
               end if
            enddo
!____________________________________________________________________________      
!________________ [2] Store the results for the present BIN _________________
!____________________________________________________________________________
            !!!!!!!!!! Scaling for structure factors by 1/NumNC
            ZTempA = ZTempA / dble(NumNC)
            !!!!!!!!!! Store results at Gamma=(0,0) Point for Average and Errorbar
            KStructFactGamm(NB, 01:16, MeaInd) = real(ZTempA(01:16))
            CorlatRatioGamm(NB, 01:16, MeaInd) = real(ZTempB(01:16))
            !!!!!!!!!! Store results at X1=(pi,0) and X2=(0,pi) Points for Average and Errorbar
            KStructFactXPnt(NB, 01:16, MeaInd) = real(ZTempA(17:32))
            CorlatRatioXPnt(NB, 01:16, MeaInd) = real(ZTempB(17:32))
            !!!!!!!!!! Store results at M=(pi/2,pi/2) Points for Average and Errorbar
            KStructFactPntM(NB, 01:16, MeaInd) = real(ZTempA(33:48))
            CorlatRatioPntM(NB, 01:16, MeaInd) = real(ZTempB(33:48))
         end if
      end if
!**************************************************************************************************     
!___________________ 2. Output the results for the present BIN simulation _________________________
!______________________ for both PERIODIC and OPEN boundary conditions ____________________________
!**************************************************************************************************
      if(amyid == amstr) then
!________________________________________________________________________________________      
!_________________ (0) Output all the r-space correlation functions _____________________
!_____________________ for both IfPyObsPBC==T and IfPyObsPBC==F cases ___________________
!________________________________________________________________________________________
         !!!!!!!!!! Open the files for the output
         open(304, file = Trim(FileOutAdd(MeaInd)) // "55_RGrnFctCrFc_All" // Trim(FileAddTxt(MeaInd)), access = "append")
         open(305, file = Trim(FileOutAdd(MeaInd)) // "55_RSpinZZCrFc_All" // Trim(FileAddTxt(MeaInd)), access = "append")
         open(306, file = Trim(FileOutAdd(MeaInd)) // "55_RSpinPMCrFc_All" // Trim(FileAddTxt(MeaInd)), access = "append")
         open(307, file = Trim(FileOutAdd(MeaInd)) // "55_RDenDenCrFc_All" // Trim(FileAddTxt(MeaInd)), access = "append")
         open(308, file = Trim(FileOutAdd(MeaInd)) // "55_RPairStCrFc_All" // Trim(FileAddTxt(MeaInd)), access = "append")
         open(309, file = Trim(FileOutAdd(MeaInd)) // "55_REdSParCrFc_All" // Trim(FileAddTxt(MeaInd)), access = "append")
         open(310, file = Trim(FileOutAdd(MeaInd)) // "55_RDWvParCrFc_All" // Trim(FileAddTxt(MeaInd)), access = "append")
         open(311, file = Trim(FileOutAdd(MeaInd)) // "55_RBndPairCrF_All" // Trim(FileAddTxt(MeaInd)), access = "append")
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
            !!!!!!!! Write the results to files
            write(304, "(I4, A, I4, A, es25.16, A, es25.16)") I1, char(9), I2, char(9), &
                                    & RealSpCrFAll(NB, Idimj, 01, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 02, MeaInd)
            write(305, "(I4, A, I4, A, es25.16, A, es25.16)") I1, char(9), I2, char(9), &
                                    & RealSpCrFAll(NB, Idimj, 03, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 17, MeaInd)
            write(306, "(I4, A, I4, A, es25.16, A, es25.16)") I1, char(9), I2, char(9), &
                                    & RealSpCrFAll(NB, Idimj, 04, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 18, MeaInd)
            write(307, "(I4, A, I4, A, es25.16, A, es25.16)") I1, char(9), I2, char(9), &
                                    & RealSpCrFAll(NB, Idimj, 05, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 19, MeaInd)
            write(308, "(I4, A, I4, A, es25.16, A, es25.16)") I1, char(9), I2, char(9), &
                                    & RealSpCrFAll(NB, Idimj, 06, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 20, MeaInd)
            write(309, "(I4, A, I4, A, es25.16, A, es25.16)") I1, char(9), I2, char(9), &
                                    & RealSpCrFAll(NB, Idimj, 07, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 21, MeaInd)
            write(310, "(I4, A, I4, A, es25.16, A, es25.16)") I1, char(9), I2, char(9), &
                                    & RealSpCrFAll(NB, Idimj, 08, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 22, MeaInd)
            write(311, "(I4, A, I4, A, es25.16, A, es25.16, A, es25.16, A, es25.16)") I1, char(9), I2, char(9), &
                        & RealSpCrFAll(NB, Idimj, 31, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 32, MeaInd), char(9), &
                        & RealSpCrFAll(NB, Idimj, 33, MeaInd), char(9), RealSpCrFAll(NB, Idimj, 34, MeaInd)
         enddo
         close(304); close(305); close(306); close(307); close(308); close(309); close(310); close(311)
!________________________________________________________________________________________      
!_________________ (1) Output the structure factors and correlation ratios ______________
!_____________________ only for IfPyObsPBC==T case ______________________________________
!________________________________________________________________________________________
         if(IfPyObsPBC) then
            open(331, file = Trim(FileOutAdd(MeaInd)) // "03_KGrnFctCorrFunc" // Trim(FileAddTxt(MeaInd)), access = "append")
            open(332, file = Trim(FileOutAdd(MeaInd)) // "03_KSpinZZCorrFunc" // Trim(FileAddTxt(MeaInd)), access = "append")
            open(333, file = Trim(FileOutAdd(MeaInd)) // "03_KSpinPMCorrFunc" // Trim(FileAddTxt(MeaInd)), access = "append")
            open(334, file = Trim(FileOutAdd(MeaInd)) // "03_KDenDenCorrFunc" // Trim(FileAddTxt(MeaInd)), access = "append")
            open(335, file = Trim(FileOutAdd(MeaInd)) // "03_KPairStCorrFunc" // Trim(FileAddTxt(MeaInd)), access = "append")
            open(336, file = Trim(FileOutAdd(MeaInd)) // "03_KEdSParCorrFunc" // Trim(FileAddTxt(MeaInd)), access = "append")
            open(337, file = Trim(FileOutAdd(MeaInd)) // "03_KDWvParCorrFunc" // Trim(FileAddTxt(MeaInd)), access = "append")
            do I0 = 331, 337, +1
               write(I0, "(I4.4, A)", advance = "no") NB, char(9)
            enddo
            write(331, "(es25.16, A, es25.16)") KStructFactXPnt(NB, 01, MeaInd)*dble(NumNC), char(9), &
               & KStructFactXPnt(NB, 02, MeaInd)*dble(NumNC)
            write(332, "(es25.16, A, es25.16, A, es25.16, A, es25.16)") &
               & KStructFactXPnt(NB, 03, MeaInd), char(9), KStructFactPntM(NB, 03, MeaInd), char(9), &
               & CorlatRatioXPnt(NB, 03, MeaInd), char(9), CorlatRatioPntM(NB, 03, MeaInd)
            write(333, "(es25.16, A, es25.16, A, es25.16, A, es25.16)") &
               & KStructFactXPnt(NB, 04, MeaInd), char(9), KStructFactPntM(NB, 04, MeaInd), char(9), &
               & CorlatRatioXPnt(NB, 04, MeaInd), char(9), CorlatRatioPntM(NB, 04, MeaInd)
            write(334, "(es25.16, A, es25.16, A, es25.16, A, es25.16)") &
               & KStructFactXPnt(NB, 05, MeaInd), char(9), KStructFactPntM(NB, 05, MeaInd), char(9), &
               & CorlatRatioXPnt(NB, 05, MeaInd), char(9), CorlatRatioPntM(NB, 05, MeaInd)
            write(335, "(es25.16, A, es25.16, A, es25.16, A, es25.16)") &
               & KStructFactGamm(NB, 06, MeaInd), char(9), KStructFactXPnt(NB, 06, MeaInd), char(9), &
               & CorlatRatioGamm(NB, 06, MeaInd), char(9), CorlatRatioXPnt(NB, 06, MeaInd)
            write(336, "(es25.16, A, es25.16, A, es25.16, A, es25.16)") &
               & KStructFactGamm(NB, 07, MeaInd), char(9), KStructFactXPnt(NB, 07, MeaInd), char(9), &
               & CorlatRatioGamm(NB, 07, MeaInd), char(9), CorlatRatioXPnt(NB, 07, MeaInd)
            write(337, "(es25.16, A, es25.16, A, es25.16, A, es25.16)") &
               & KStructFactGamm(NB, 08, MeaInd), char(9), KStructFactXPnt(NB, 08, MeaInd), char(9), &
               & CorlatRatioGamm(NB, 08, MeaInd), char(9), CorlatRatioXPnt(NB, 08, MeaInd)
            close(331); close(332); close(333); close(334); close(335); close(336); close(337)
         end if
!________________________________________________________________________________________      
!_________________ (2) Momentum distribution For IfFftEnPar == .false. case _____________
!_____________________ only for IfPyObsPBC==T case ______________________________________
!________________________________________________________________________________________
         if( (IfPyObsPBC) .and. (.not. IfFftEnPar) ) then
            open(291, file = Trim(FileOutAdd(MeaInd)) // "MomtumDistBINs" // Trim(FileAddTxt(MeaInd)), access = "append")
            do Nk = 1, NumNC, +1
               NkDistrib(NB, Nk, 1, MeaInd) = real( KSpaceCrFAll(NB, Nk, 01, MeaInd) )
               NkDistrib(NB, Nk, 2, MeaInd) = real( KSpaceCrFAll(NB, Nk, 02, MeaInd) )
               write(291, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
               write(291, "(A, es25.16, A, es25.16)") char(9), NkDistrib(NB, Nk, 1, MeaInd), &
                                                    & char(9), NkDistrib(NB, Nk, 2, MeaInd)
            enddo
            close(291)
         end if
      end if

   end subroutine ProcRKStaCrF
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine VertexContrbSta(MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  VertexContrbSta(MeaInd) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compute the vertex contribution for various bosonic correlation functions. 
! KEYWORDS: Compute the vertex contribution. For the static correlation functions.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Fourier Transformation.
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
      integer Idimj, I0, I1, I2, I3, I4, Kd, Hd
      real(rp) RlCrFTmp(15), Rtp0
      real(rp), allocatable :: Collect0(:, :)
!______________________________________________________________________________________________________________     
!_____________________________ Compute the vertex contribution for correlations _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Vertex contribution = C_{QMC} - C_{uncorrelated} __________________________
!************************************************************************************************** 
!________________________________________________________________________________________      
!_________________ (0) Allocate the temporary matrice ___________________________________
!________________________________________________________________________________________ 
      allocate(Collect0(NmSitePair, 14)); Collect0 = rp_Zzero
!________________________________________________________________________________________      
!_________________ (1) Compute the uncorrelated part of the correlations ________________
!________________________________________________________________________________________
      !!!!!!!!!! Computation applying translational symmetry
   !$OMP PARALLEL &
   !$OMP PRIVATE(Idimj, I0, I1, I2, I3, I4, Kd, Hd, Rtp0, RlCrFTmp)
   !$OMP DO REDUCTION(+ : Collect0)
      do I2 = 1, NumNC, +1
         do I1 = 1, merge(NumNC, I2, IfPyObsPBC), +1
            !!!!!!!! The integer index for the present term 
            Idimj = merge(IminusJ(I1, I2), (I2-1)*I2/2+I1, IfPyObsPBC)
            !!!!!!!! The real constant, 1 for I1==I2 and 0 for I1/=I2
            Rtp0 = merge(1.0_rp, 0.0_rp, I1==I2)
            !!!!!!!! For the spin correlations, <SzSz> and <S+S- + S-S+>/2
            !!!!!! The <SzSz> correlation
            RlCrFTmp(1) =   + ( RlSpGrnFtBIN(I1, I1, 1, MeaInd) - RlSpGrnFtBIN(I1, I1, 2, MeaInd) ) * &
                                             & ( RlSpGrnFtBIN(I2, I2, 1, MeaInd) - RlSpGrnFtBIN(I2, I2, 2, MeaInd) ) &
                          & + RlSpGrnFtBIN(I1, I2, 1, MeaInd)*( Rtp0 - RlSpGrnFtBIN(I2, I1, 1, MeaInd) ) &
                          & + RlSpGrnFtBIN(I1, I2, 2, MeaInd)*( Rtp0 - RlSpGrnFtBIN(I2, I1, 2, MeaInd) )
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            Collect0(Idimj, 01) = Collect0(Idimj, 01) + RlCrFTmp(1)
            !!!!!! The Sxy-Sxy correlation, <S+S- + S-S+>/2
            RlCrFTmp(2) =   + RlSpGrnFtBIN(I1, I2, 1, MeaInd)*( Rtp0 - RlSpGrnFtBIN(I2, I1, 2, MeaInd) ) &
                          & + RlSpGrnFtBIN(I1, I2, 2, MeaInd)*( Rtp0 - RlSpGrnFtBIN(I2, I1, 1, MeaInd) )
            RlCrFTmp(2) = RlCrFTmp(2) / 2.0_rp
            Collect0(Idimj, 02) = Collect0(Idimj, 02) + RlCrFTmp(2)
            !!!!!!!! For the density-density correlation
            RlCrFTmp(1) =   + ( RlSpGrnFtBIN(I1, I1, 1, MeaInd) + RlSpGrnFtBIN(I1, I1, 2, MeaInd) ) * &
                                             & ( RlSpGrnFtBIN(I2, I2, 1, MeaInd) + RlSpGrnFtBIN(I2, I2, 2, MeaInd) ) &
                          & + RlSpGrnFtBIN(I1, I2, 1, MeaInd)*( Rtp0 - RlSpGrnFtBIN(I2, I1, 1, MeaInd) ) &
                          & + RlSpGrnFtBIN(I1, I2, 2, MeaInd)*( Rtp0 - RlSpGrnFtBIN(I2, I1, 2, MeaInd) )
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            Collect0(Idimj, 03) = Collect0(Idimj, 03) + RlCrFTmp(1)
            !!!!!!!! For the on-site spin-singlet pairing
            RlCrFTmp(1) =   + RlSpGrnFtBIN(I1, I2, 1, MeaInd)*RlSpGrnFtBIN(I1, I2, 2, MeaInd) &
                          & + ( Rtp0 - RlSpGrnFtBIN(I2, I1, 1, MeaInd) ) * ( Rtp0 - RlSpGrnFtBIN(I2, I1, 2, MeaInd) )
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            Collect0(Idimj, 04) = Collect0(Idimj, 04) + RlCrFTmp(1)
            !!!!!!!! For all the nearest-neighbor pairing correlations
            RlCrFTmp = 0.0_rp
            do Kd = 1, 4, +1
               I3 = FNNBond(I1, Kd)
               do Hd = 1, 4, +1
                  I4 = FNNBond(I2, Hd)
                  RlCrFTmp(15) =   + RlSpGrnFtBIN(I1, I2, 1, MeaInd)*RlSpGrnFtBIN(I3, I4, 2, MeaInd) &
                                 & + RlSpGrnFtBIN(I1, I4, 1, MeaInd)*RlSpGrnFtBIN(I3, I2, 2, MeaInd) &
                                 & + RlSpGrnFtBIN(I3, I2, 1, MeaInd)*RlSpGrnFtBIN(I1, I4, 2, MeaInd) &
                                 & + RlSpGrnFtBIN(I3, I4, 1, MeaInd)*RlSpGrnFtBIN(I1, I2, 2, MeaInd)
                  !!!!!! For the extended s-wave pairing
                  RlCrFTmp(1) = RlCrFTmp(1) + RlCrFTmp(15)
                  !!!!!! For the d-wave pairing
                  RlCrFTmp(2) = RlCrFTmp(2) + dble(1-2*mod(Kd, 2)) * dble(1-2*mod(Hd, 2)) * RlCrFTmp(15)
               enddo
            enddo
            !!!!!! For the extended s-wave pairing
            RlCrFTmp(1) = RlCrFTmp(1) / 32.0_rp
            Collect0(Idimj, 05) = Collect0(Idimj, 05) + RlCrFTmp(1)
            !!!!!! For the NN extended s-wave pairing
            RlCrFTmp(2) = RlCrFTmp(2) / 32.0_rp
            Collect0(Idimj, 06) = Collect0(Idimj, 06) + RlCrFTmp(2)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
      !!!!!!!!!! The size-related constant
      if(IfPyObsPBC) Collect0 = Collect0 / dble(NumNC)
!________________________________________________________________________________________      
!_________________ (2) Compute the vertex contributions _________________________________
!________________________________________________________________________________________
      do Idimj = 1, NmSitePair, +1
         RealSpCrFBIN(Idimj, 17:30, MeaInd) = RealSpCrFBIN(Idimj, 03:16, MeaInd) - Collect0(Idimj, 01:14)
      enddo
!________________________________________________________________________________________      
!_________________ (3) Deallocate the temporary matrix __________________________________
!________________________________________________________________________________________
      if(allocated(Collect0)) deallocate(Collect0)

   end subroutine VertexContrbSta
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine DeductCrFBkgSta(MeaInd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  DeductCrFBkgSta(MeaInd) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to deduct the background for the correlation functions in bosonic channels using
!                   RlSpGrnFtBIN results. For the static correlation functions.
! KEYWORDS: Deduct the background for correlations.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Fourier Transformation.
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
      integer I1, I2, I3, I4, Idimj
      real(rp) RlCrFTmp(15)
      real(rp), allocatable :: Collect0(:, :)
!______________________________________________________________________________________________________________     
!_____________________________ Deduct the background for the correlation functions ____________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Computation for Deducting the background for correlations _________________
!************************************************************************************************** 
!________________________________________________________________________________________      
!_________________ (0) Allocate the temporary matrice ___________________________________
!________________________________________________________________________________________ 
      allocate(Collect0(NmSitePair, 5)); Collect0 = 0.0_rp
!________________________________________________________________________________________      
!_________________ (1) Compute the background for correlations __________________________
!________________________________________________________________________________________
      !!!!!!!!!! Computation applying translational symmetry
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2, I3, I4, Idimj, RlCrFTmp)
   !$OMP DO REDUCTION(+ : Collect0)
      do I2 = 1, NumNC, +1
         do I1 = 1, merge(NumNC, I2, IfPyObsPBC), +1
            !!!!!!!! The integer index for the present term 
            Idimj = merge(IminusJ(I1, I2), (I2-1)*I2/2+I1, IfPyObsPBC)
            !!!!!!!! For the Sz-Sz correlation function
            RlCrFTmp(1) =     ( RlSpGrnFtBIN(I1, I1, 1, MeaInd) - RlSpGrnFtBIN(I1, I1, 2, MeaInd) ) &
                          & * ( RlSpGrnFtBIN(I2, I2, 1, MeaInd) - RlSpGrnFtBIN(I2, I2, 2, MeaInd) ) / 4.0_rp
            Collect0(Idimj, 01) = Collect0(Idimj, 01) + RlCrFTmp(1)
            !!!!!!!! For the density-density correlation function
            RlCrFTmp(1) =     ( RlSpGrnFtBIN(I1, I1, 1, MeaInd) + RlSpGrnFtBIN(I1, I1, 2, MeaInd) ) &
                          & * ( RlSpGrnFtBIN(I2, I2, 1, MeaInd) + RlSpGrnFtBIN(I2, I2, 2, MeaInd) ) / 4.0_rp
            Collect0(Idimj, 02) = Collect0(Idimj, 02) + RlCrFTmp(1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL         
      !!!!!!!!!! The size-related constant
      if(IfPyObsPBC) Collect0 = Collect0 / dble(NumNC)
!________________________________________________________________________________________      
!_________________ (2) Deduct the background in the correlations ________________________
!________________________________________________________________________________________
      RealSpCrFBIN(:, 03, MeaInd) = RealSpCrFBIN(:, 03, MeaInd) - Collect0(:, 01)    ! Sz-Sz correlation
      RealSpCrFBIN(:, 05, MeaInd) = RealSpCrFBIN(:, 05, MeaInd) - Collect0(:, 02)    ! Density-density correlation
!________________________________________________________________________________________      
!_________________ (3) Deallocate the temporary matrix __________________________________
!________________________________________________________________________________________
      if(allocated(Collect0)) deallocate(Collect0)

   end subroutine DeductCrFBkgSta
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine FourierTransSta(NB, MeaInd) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FourierTransSta(NB, MeaInd) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to Perform Fourier transformation for SpinZZ, SpinPM, GreenF, DenDen into 
!                   reciprocal space.
! KEYWORDS: Fourier transformation for equal-time quantities.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Fourier Transformation.
!
!     Input: NB     --> Integer index for the present BIN for KSpaceCrFAll array;
!            MeaInd --> The measure kind; == 0, only measure at BetaTl; == 1, measure BetaT + [BetaT, 0];
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
      integer NB, MeaInd
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Nk, imj, I0, I1, I2
      real(rp) VecK(2), VecR(2)
      real(rp) Rtp1
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for Fourier Transformation ___________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Fourier Transformation for all the data ___________________________________
!************************************************************************************************** 
!________________________________________________________________________________________      
!_________________ (0) Re-Initialization of KSpaceCrFAll results ________________________
!________________________________________________________________________________________          
      KSpaceCrFAll(NB, :, :, MeaInd) = rp_Zzero
!________________________________________________________________________________________      
!_________________ (1) Data process for all NTInd _______________________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(Nk, VecK, imj, I1, I2, VecR, Rtp1)
   !$OMP DO
      do Nk = 1, NumNC, +1
         VecK = KpList(Nk, 1) * B1Vec + KpList(Nk, 2) * B2Vec
         do imj = 1, NumNC, +1
            I1 = StList(imj, 1)
            I2 = StList(imj, 2)
            I1 = I1 - NumL1/2 - mod(NumL1, 2)
            I2 = I2 - NumL2/2 - mod(NumL2, 2)
            VecR = I1 * A1Vec + I2 * A2Vec
            Rtp1 = dot_product(VecK, VecR)
            do I0 = 1, 30, +1
               KSpaceCrFAll(NB, Nk, I0, MeaInd) = KSpaceCrFAll(NB, Nk, I0, MeaInd) + &
                                          & exp(cmplx(0.0_rp, Rtp1, rp)) * RealSpCrFBIN(imj, I0, MeaInd)
            enddo
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
      
   end subroutine FourierTransSta
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$