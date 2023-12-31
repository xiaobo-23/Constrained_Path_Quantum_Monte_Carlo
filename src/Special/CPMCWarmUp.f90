!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine CPMCWarmUp()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCWarmUp()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the warm up of the CPMC simulation.
! KEYWORDS: CPMC Simulation process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Perform the warm up of the CPMC simulation. 
!
!     Input:  (none)   Output: (none)
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
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      logical IfMeasure
      integer NB, NSW, Iwalk, nWarmHere
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for warm up part process __________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of CPMC Warm up ____________________________________________
!______________________________________________________________________________________________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%      
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)	
#endif
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCWarmUp process ___________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         if(IfMuTqmcNt .or. IfFixnT) then
            write(*, "(33x, 'Perform warm up part for present iteration!')")
         else
            write(*, "()")
            write(*, "(16x, 'CPMCWarmUp: Perform Warm up part of CPMC simulation!')")
            write(*, "(28x, 'ChemP    == ', sp, es23.16)") ChemP
            write(*, "(28x, 'ChemP_BT == ', sp, es23.16)") ChemP_BT
            if(ReadFldMea) then
               write(*, "(28x, 'ReadFldMea == T. No need of warm up process!')")
            end if
         end if
      end if
!**************************************************************************************************	  
!___________________ 0. Some initializations for warm up process __________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Initializations of time-counting module for present BIN __________
!________________________________________________________________________________________			
      call QMCTimeRecInit()
!________________________________________________________________________________________ 	  
!_________________ (1) Set Warm up with NB == 0 and no measurement ______________________
!________________________________________________________________________________________
      NB = 0; IfMeasure = .false.
!________________________________________________________________________________________ 	  
!_________________ (2) Setting for outputting information during simulation _____________
!________________________________________________________________________________________ 
      nWarmHere = nWarm
      if(Fix_Iterate >= 1) nWarmHere = 2
      StabOutput = nWarmHere
      PoptOutput = nWarmHere
!**************************************************************************************************	  
!___________________ 1. Perform the warm up process for the simulation ____________________________
!************************************************************************************************** 
      do NSW = 1, nWarmHere, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Propagate from \tau=0 to \tau=BetaT to sample paths ______________
!________________________________________________________________________________________
         call SweepOne2M(IfMeasure, NB, NSW)
      enddo
!**************************************************************************************************	  
!___________________ 2. Store the sampled paths (configurations) and Growth estimator _____________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Store the paths (configurations) from CP construction ____________
!________________________________________________________________________________________
      if(IfSaveFlds) call SaveFewCfg()
!________________________________________________________________________________________ 	  
!_________________ (1) Store the growth estimator _______________________________________
!________________________________________________________________________________________
      call RdWtGrowth("Save")
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for warm up part process __________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&     
		call system_clock(time2)
      TimeSgBIN = TimeIntrvl(time1, time2)
!**************************************************************************************************	  
!___________________ 3. Output consumed time of this starting part of CPMC ________________________
!**************************************************************************************************
      call SaveCalTim(NB)
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%      
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)	
#endif
		
   end subroutine CPMCWarmUp
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$