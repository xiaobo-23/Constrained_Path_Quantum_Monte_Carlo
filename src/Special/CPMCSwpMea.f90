!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine CPMCSwpMea() 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCSwpMea() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the sweep for the LTrot time slices (and flipping), during the sweep 
!              we also perform the meansurement for both equal-time and unequal-time physical quantities.
! KEYWORDS: Sweep and Measurement.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Totally, we perform NmBin different sweep and measurement.
!
!     In every Bin, we have NSwep different, full sweep for the space-time lattices;
!     A full sweep means a sweep from t=0 to t=M and, then from t=M back to t=1.
!
!     Input: (none);  Outpt: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use TimeRecord
      use QMCTimeRec
      use RandomNumb
      use Observable
      use CoreParamt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      logical IfMeasure
      integer Iwalk, NB, NSW
      character(40) FlBIN
!______________________________________________________________________________________________________________     
!____________________________ Main calculations of Sweeps and Measures ________________________________________
!______________________________________________________________________________________________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%      
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr) 
#endif
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCSwpMea process ___________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'CPMCSwpMea: Perform Sweep and Measure part of CPMC simulation!')")
         if(SaveFldMea) then
            write(*, "(28x, 'SaveFldMea == T. Save all the sampled paths for measurements!')")
         end if
      end if
!**************************************************************************************************     
!___________________ 0. Iterations of all the BINs for sweeps and measurements ____________________
!**************************************************************************************************      
      do NB = 1, NmBin, +1
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v% 
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of the number of present BIN ____________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(28x, 'BIN = ', I4.4)") NB
         end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!__________________________________ Counting time consumed for present BIN ____________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
         write(FlBIN, "('BIN_', I4.4)") NB
         call SubroutineBgn(Trim(FlBIN), 15, time1)
!**************************************************************************************************     
!___________________ 1. Some initializations for the Present BIN simulation _______________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Open Tag for measurements in sweeps and measures _________________
!_____________________ and the StabOutput parameter for sweep and measure _______________
!________________________________________________________________________________________          
         IfMeasure = .true.
!________________________________________________________________________________________      
!_________________ (1) Initializations of observables for present BIN ___________________
!________________________________________________________________________________________
         call CPMCObservInit()
!________________________________________________________________________________________      
!_________________ (2) Initializations of time-counting module for present BIN __________
!________________________________________________________________________________________       
         call QMCTimeRecInit()
!________________________________________________________________________________________      
!_________________ (3) Setting for outputting information during simulation _____________
!________________________________________________________________________________________    
         !!!!!!!!!! The sweep index for output NumStab and population control  
         StabOutput = NvStbOut
         PoptOutput = NvPopOut
         !!!!!!!!!! Randomly choose a sweep to output weight information
         WghtOutNSW = -10000
         if(IfWghtMxMn) then
            if(amyid == amstr) WghtOutNSW = ceiling( spring_sfmt_stream()*NSwep )
#ifdef MPIPROCESS
            call MPI_BCAST(WghtOutNSW, 1, MPI_INTEGER, amstr, acomm, ierr)
#endif
         end if
!**************************************************************************************************     
!___________________ 2. Iteration of all the sweeps in present BIN simulation _____________________
!**************************************************************************************************         
         do NSW = 1, NSwep, +1
!________________________________________________________________________________________      
!_________________ (0) Propagate from \tau=0 to \tau=BetaT to sample paths ______________
!_____________________ And perform measurements at \tau=BetaT point _____________________
!________________________________________________________________________________________
            call SweepOne2M(IfMeasure, NB, NSW)
!________________________________________________________________________________________      
!_________________ (1) Propagate from \tau=BetaT to \tau=0 for additonal ________________
!_____________________ static and dynamic measurements __________________________________
!________________________________________________________________________________________
            if(IfM2OneMea) call SweepM2One(NB, NSW)
!________________________________________________________________________________________      
!_________________ (2) Process results of physical observables for an single sweep ______
!________________________________________________________________________________________
            call SwpDatProc()
         enddo
!**************************************************************************************************     
!___________________ 3. Store the sampled paths (configurations) and Growth estimator _____________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Store the paths (configurations) from CP construction ____________
!________________________________________________________________________________________
         if(IfSaveFlds) call SaveFewCfg()
!________________________________________________________________________________________      
!_________________ (1) Store the growth estimator _______________________________________
!________________________________________________________________________________________
         call RdWtGrowth("Save")
!**************************************************************************************************     
!___________________ 4. PostProcess for both the static and Dynamic quantities ____________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Process the results of static measurements _______________________
!________________________________________________________________________________________
!____________________________________________________________________________      
!________________ [0] Results for measurements only at BetaT ________________
!____________________________________________________________________________
         call PostStatic(NB, 0)
!____________________________________________________________________________      
!________________ [1] Results for measurements at BetaT + [BetaT, 0] ________
!____________________________________________________________________________
         if(IfM2OneMea) call PostStatic(NB, 1)
!________________________________________________________________________________________      
!_________________ (1) Process the results of dynamic measurements ______________________
!________________________________________________________________________________________
         if(IfTAU) call PostDynamc(NB)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!__________________________________ Counting time consumed for present BIN ____________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
         call SubroutineEnd(Trim(FlBIN), 15, time1, time2)
         TimeSgBIN = TimeIntrvl(time1, time2)
!**************************************************************************************************   
!___________________ 5. Output consumed time of present BIN of CPMC simulation ____________________
!**************************************************************************************************
         call SaveCalTim(NB)
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%      
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr) 
#endif   
      enddo
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%      
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr) 
#endif
      
   end subroutine CPMCSwpMea
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$