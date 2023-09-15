!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to Calculate the growth estimator, namely the constants accumulated during the
!               propagation, to avoid too large or too small weights.
! COMMENT: Calculate the growth estimator.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   GrwthCoeff --> Used for tunning the constants in the weights;
!   RdWtGrowth --> Read and write the growth estimators from the QMC simulations.
!           
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine GrwthCoeff(NB, NSW, NT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrwthCoeff(NB, NSW, NT)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to adjust the constants GrowthCoefft appearing in the weights of random walkers.
!                  GrowthCoefft is used to avoid the appearance of too small or too large weights. 
! KEYWORDS: Adjust the constant in weights.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Adjust the constant in weights.
!
!     Input:  NB  --> Integer index for BIN;
!             NSW --> Integer index for sweep;
!             NT  --> The time slice.
!     
!     Output: (none)
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
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB, NSW, NT
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      integer I0, I1, I2
      integer NTInd, NT_Init, NT_Finl
      real(rp) WeightSum, Rtp0, Rtp1, Rtp2
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for Calculating growth estimator __________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      TimsConst = TimsConst + 1
      call system_clock(time1)
!______________________________________________________________________________________________________________     
!___________________________ Main calculations of adjusting the constant ______________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Initializations for tunning the growth estimator constant _________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Determine time slice interval of tunning constants _______________
!________________________________________________________________________________________ 
      !!!!!!!!!! For NT_Init
      if( (NT == LTrot) .and. (mod(LTrot, NvPop) /= 0) ) then
         NT_Init = LTrot - mod(LTrot, NvPop) + 1
      else
         NT_Init = NT - NvPop + 1
      end if
      !!!!!!!!!! For NT_Finl
      if( (NB == 0) .and. (NSW == 1) ) then
         NT_Finl = LTrot
      else
         NT_Finl = NT
      end if
!________________________________________________________________________________________      
!_________________ (1) Accumulate all weights at the present NT point ___________________
!________________________________________________________________________________________ 
      Sum_WghtNew = sum(WghtProc)
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)
      WeightSum = 0.0_rp
      call MPI_ALLREDUCE(Sum_WghtNew, WeightSum, 1, rp_MPI_REAL, MPI_SUM, acomm, ierr)
      Sum_WghtNew = WeightSum
#endif
!**************************************************************************************************     
!___________________ 1. Tune growth estimator constant after population control ___________________
!______________________ Sum_WghtNew --> Summation of all weights at NT         point ______________
!______________________ Sum_WghtOld --> Summation of all weights at NT-NvPop+1 point ______________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) For the case of NT/NvPop == 1 ____________________________________
!________________________________________________________________________________________ 
      if( NT/NvPop == 1 ) then
!____________________________________________________________________________      
!________________ [0] [exp(facNew)/exp(fac)]^{NvPop} = ______________________
!____________________ Sum_WghtOld/(Sum_WghtNew*MeanWghtPop) _________________
!____________________________________________________________________________
         Rtp1 = ( log(Sum_WghtOld) - log(Sum_WghtNew) - log(MeanWghtPop) ) / dble(NvPop)
!____________________________________________________________________________      
!________________ [1] GrowthCoefft(NT-NvPop+1 : NT) = ... / Const ___________
!____________________________________________________________________________
         do NTInd = NT_Init, NT_Finl, +1
            GrowthCoefft(NTInd) = GrowthCoefft(NTInd) + Rtp1
         enddo
!________________________________________________________________________________________      
!_________________ (2) For the case of NT/NvPop > 1 _____________________________________
!________________________________________________________________________________________ 
      else if( NT/NvPop > 1 ) then
!____________________________________________________________________________      
!________________ [0] [exp(facNew)/exp(fac)]^{NvPop} = ______________________
!____________________ Sum_WghtOld/(Sum_WghtNew*MeanWghtPop) _________________
!____________________________________________________________________________
         if( (mod(LTrot, NvPop) /= 0) .and. (NT == LTrot) ) then
            Rtp1 = ( log(Sum_WghtOld) - log(Sum_WghtNew) - log(MeanWghtPop) ) / dble(mod(LTrot, NvPop))
         else
            Rtp1 = ( log(Sum_WghtOld) - log(Sum_WghtNew) - log(MeanWghtPop) ) / dble(NvPop)
         end if
!____________________________________________________________________________      
!________________ [1] GrowthCoefft(NT-NvPop+1 : NT) = ... / Const ___________
!____________________________________________________________________________
         do NTInd = NT_Init, NT_Finl, +1
            GrowthCoefft(NTInd) = GrowthCoefft(NTInd) + Rtp1
         enddo
      end if
!________________________________________________________________________________________      
!_________________ (3) Reset Sum_WghtOld and Sum_WghtNew for next run ___________________
!________________________________________________________________________________________
      Sum_WghtOld = Sum_WghtNew
      Sum_WghtNew = 0.0_rp
!**************************************************************************************************     
!___________________ 2. Output the information of tunning constants here __________________________
!**************************************************************************************************
      if(amyid == amstr) then
         !!!!!!!!!! Print head of output for the present BIN
         if(NB == 0 .and. NSW == PoptOutput .and. NT == NvPop) then
            open( 291, file = "Output/00_Adjust_Constant.txt", access = "append")
            write(291, "()")
            write(291, "()")
            if(IfMuTqmcNt .or. IfFixnT) then
               write(291, "()")
               write(291, "('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')")
               if(IfMuTqmcNt) then
                  write(291, "('@@@@@@@@@@@ Adjust ChemP_BT to reach the QMC density @@@@@@@@@@@@@@@@')")
               else if(IfFixnT) then
                  write(291, "('@@@@@@@@@@@ Adjust ChemP to reach fixed electron density @@@@@@@@@@@@')")
               end if
               write(291, "('@@@@@@@@@@@@@@@@@@@@@@@@@@ Iteration ', I3.3, ' @@@@@@@@@@@@@@@@@@@@@@@@@@@@')") &
                  & Fix_Iterate
               write(291, "('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')")
               write(291, "()")
               write(291, "('________________________________________________________________')")
               write(291, "('__________________________ CPMCWarmUp __________________________')") 
               write(291, "('________________________________________________________________')")
            else
               write(291, "()")
               write(291, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
               write(291, "('&&&&&&&&&&&&&&&& CPMC Simulation for Data Statistics &&&&&&&&&&&&&&&&')")
               write(291, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
               write(291, "()")
               write(291, "('################################################################')")
               write(291, "('######################### CPMC WarmUp ##########################')") 
               write(291, "('################################################################')")
            end if
            close(291)
         else if(NB > 0 .and. NSW == NSwep .and. NT == NvPop) then
            open( 291, file = "Output/00_Adjust_Constant.txt", access = "append")
            write(291, "()")
            write(291, "()")
            if(IfMuTqmcNt) then
               write(291, "('________________________________________________________________')")
               write(291, "('__________________________ CPMCOccMuT __________________________')") 
               write(291, "('________________________________________________________________')")
            else if(IfFixnT) then
               write(291, "('________________________________________________________________')")
               write(291, "('__________________________ CPMCMeaFix __________________________')") 
               write(291, "('________________________________________________________________')")
            else
               write(291, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
               write(291, "('++++++++++++++++++++++++++ BIN ', I4.4,' ++++++++++++++++++++++++++++')") NB
               write(291, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
            end if
            close(291)
         end if
         !!!!!!!!!! Print Information of growth estimator for the PoptOutput sweep
         if( (NB == 0 .and. NSW == PoptOutput) .or. (NB > 0 .and. NSW == NSwep) ) then
            open(291, file = "Output/00_Adjust_Constant.txt", access = "append")
            if(NT == NvPop) then
               write(291, "()")
               write(291, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
               write(291, "('&&&&&&&&&&&&&&&&&&&&&&&& Sweep ', I4.4,' &&&&&&&&&&&&&&&&&&&&&&&')") NSW
               write(291, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
            end if
            write(291, "('NT, log(Diff), Diff, log[Const(NT)], Const(NT) = ', I5.5, 4es15.6)") NT, Rtp1, &
               & exp(Rtp1), GrowthCoefft(NT), exp(GrowthCoefft(NT))
            close(291)
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for Calculating growth estimator __________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeConst = TimeConst + TimeIntrvl(time1, time2)
      
   end subroutine GrwthCoeff
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine RdWtGrowth(SaveOrRead)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RdWtGrowth(SaveOrRead)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to store the growth estimator in the AFQMC simulations.
! KEYWORDS: Number convertion.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-25
! DESCRIPTION: Save the growth estimator to the output file.
!
!     Input: SaveOrRead --> == "Save" for saving configurations; == "Read" for reading configurations.
!                     
!     Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      character(4) SaveOrRead
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Stat, NTInd, LTrotTempt
      real(rp) Rtp0
!______________________________________________________________________________________________________________     
!___________________________ Main calculations of storing growth estimator ____________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************     
!___________________ 0. Store the growth estimator from the simulation ____________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Save the growth estimator results ________________________________
!________________________________________________________________________________________
      if(SaveOrRead == "Save") then
         if(amyid == amstr) then
            open(299, file = "Output/GrowthEstimator.txt")
            do NTInd = 1, LTrot, +1
               write(299, "(es25.16)") GrowthCoefft(NTInd)
            enddo
            close(299)
         end if
!________________________________________________________________________________________      
!_________________ (1) Read the growth estimator for following simulation _______________
!________________________________________________________________________________________
      else if(SaveOrRead == "Read") then
         !!!!!!!!!! Read the growth estimator
         !!!!!!!! Determine the LTrot of the input file
         open(299, err = 58, file = "Output/GrowthEstimator.txt", status = "old")
         LTrotTempt = 0
         do while(.true.)
            read(299, *, iostat = Stat) Rtp0
            if(Stat /= 0) exit
            LTrotTempt = LTrotTempt + 1
         enddo
         close(299)
         !!!!!!!! Read the growth estimators
         open(299, err = 58, file = "Output/GrowthEstimator.txt", status = "old")
         do NTInd = 1, merge(LTrot, LTrotTempt, LTrot<=LTrotTempt), +1
            read(299, *) GrowthCoefft(NTInd)
         enddo
         !!!!!!!! For LTrot > LTrotTempt case, fill the rest
         if(LTrot > LTrotTempt) then
            GrowthCoefft(LTrotTempt+1:LTrot) = GrowthCoefft(LTrotTempt)
         end if
         close(299)
         if(amyid == amstr) write(*, "(28x, 'Read the growth estimators from the input file!')")
         go to 59
         !!!!!!!!!! If the file doesn't exist, initialize GrowthCoefft
58       continue
         GrowthCoefft(1:LTrot) = log(1.0_rp)
59       continue
      end if
      
   end subroutine RdWtGrowth
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$