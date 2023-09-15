!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to compute the ChemP_BT of H_T Hamiltonian for the case of fixed \mu simulations, 
!              by self-consistent QMC calculations.
! COMMENT: Determine ChemP_BT to make n_{H_T} = n.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!  
!   CPMCGetMuT --> Subroutine used to perform all the calculations to determine ChemP_BT.
!             
!   CPMCOccMuT --> Subroutine used to perform the sweep as [0, Beta] to compute electron density;
!   ComptnTMuT --> Subroutine used to calculate the density average and error bar;
!          
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine CPMCGetMuT()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCGetMuT()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform calculations to determine ChemP_BT of H_T Hamiltonian to make 
!                 n_T = n. 
! KEYWORDS: CPMC measurement of total electron density.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Self-consistent iterations to determine ChemP_BT.
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
      use Observable
      use CoreParamt
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2, time3, time4 
      real(rp) nT_BT
      character(40) FlBIN  ! Write the integer NB into this character 
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of CPMC Simulation _________________________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of the process of present subroutine ____________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'CPMCGetMuT: Warm up and measurement for simulations to determine ChemP_BT!')")
         write(*, "()")
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	  
!___________________ 0. Iterations of all the iterations for the measurements _____________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do Fix_Iterate = 0, NItrGetMuT-1, +1
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v% 
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)		
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!__________________________________ Counting time consumed for present BIN ____________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
         write(FlBIN, "(I3.3, ' Iteration')") Fix_Iterate
         call SubroutineBgn(Trim(FlBIN), 15, time3)
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of the number of present BIN ____________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
         if(amyid == amstr) then
            write(*, "(28x, 'Iteration = ', I2.2)") Fix_Iterate
            write(*, "(33x, 'Input ChemP_BT = ', sp, es19.12)") ChemP_BT
         end if
!________________________________________________________________________________________ 	  
!_________________ (0) The warm up process for the present Fix_Iterate __________________
!________________________________________________________________________________________
         call SubroutineBgn("CPMCWarmUp", 19, time1)
         call CPMCWarmUp() 
         call SubroutineEnd("CPMCWarmUp", 19, time1, time2)
!________________________________________________________________________________________ 	  
!_________________ (1) The sweep for measuring total density ____________________________
!________________________________________________________________________________________
         call SubroutineBgn("CPMCOccMuT", 19, time1)
         call CPMCOccMuT() 
         call SubroutineEnd("CPMCOccMuT", 19, time1, time2)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!__________________________________ Counting time consumed for present BIN ____________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
         call SubroutineEnd(Trim(FlBIN), 15, time3, time4)
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%      
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)	
#endif
      enddo
!&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$$$
!___________________ 1. Calculate GreenF and static observables for H_T Hamiltonian _______________
!&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$$$
!________________________________________________________________________________________  
!_________________ (0) Calculate some Physical Quantities for H_T Hmlt __________________
!________________________________________________________________________________________
      call PhyMeaStatHT()
!________________________________________________________________________________________  
!_________________ (1) Filling of HT Hamiltonian and output _____________________________
!________________________________________________________________________________________
      call Calculate_nTOfHT(nT_BT)     
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)
      if( (mod(amyid, merge(anprc/4, 1, anprc>4)) == 0) .or. (amyid == anprc-1) ) then
         write(*, "(28x, 'For Process ', I4.4, ', ChemP_BT, nT_BT = ', sp, es19.12, ss, es20.12)") amyid, ChemP_BT, nT_BT
      end if
      call MPI_Barrier(acomm, ierr)
#else
         write(*, "(28x, 'For Process ', I4.4, ', ChemP_BT, nT_BT = ', sp, es19.12, ss, es20.12)") amyid, ChemP_BT, nT_BT
#endif
!&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$$$
!___________________ 2. Reset IfMuTqmcNt to be IfMuTqmcNt = .false. _______________________________
!&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$&&&&####$$$$$$
      IfMuTqmcNt = .false.
		
   end subroutine CPMCGetMuT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine CPMCOccMuT()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCOccMuT()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the measurements of total electron density for the many-body system
!                 and then compute the ChemP_BT to make nT = n.
! KEYWORDS: Compute ChemP_BT.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Only measure the total electron density.
!
!     Input: (none); Output: (none)
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
      logical IfMeasure
      integer NB, NSW, Iwalk
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
         write(*, "(33x, 'Measure the total density and compute ChemP_BT!')")
      end if
!**************************************************************************************************	  
!___________________ 0. Some initializations for warm up process __________________________________
!**************************************************************************************************			
!________________________________________________________________________________________ 	  
!_________________ (0) Close Tag for measurements in sweeps and measures ________________
!________________________________________________________________________________________	         
      IfMeasure = .true.; NObsStat = 0
      WeightList = + 00.00_rp; n_Occ_List = - 100.0_rp
      NB = 8000 + Fix_Iterate
!________________________________________________________________________________________ 	  
!_________________ (1) Initializations of time-counting module for present BIN __________
!________________________________________________________________________________________			
      call QMCTimeRecInit()
!________________________________________________________________________________________ 	  
!_________________ (2) Setting for outputting information during simulation _____________
!________________________________________________________________________________________      
      StabOutput = NvStbOut
      PoptOutput = NvPopOut
!**************************************************************************************************	  
!___________________ 1. Carry out the warm up process for the simulation __________________________
!************************************************************************************************** 
      do NSW = 1, NSwep, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Propagate from \tau=0 to \tau=BetaT to sample paths ______________
!_____________________ And perform measurements at \tau=BetaT point _____________________
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
!**************************************************************************************************	  
!___________________ 3. Compute the total density and determine ChemP_BT __________________________
!______________________ And compute the initial GreenF matrix _____________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Compute the total density and determine ChemP_BT _________________
!________________________________________________________________________________________
      call ComptnTMuT()
!________________________________________________________________________________________ 	  
!_________________ (1) The initial GreenF matrix for H_T Hamiltonian ____________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         write(*, "(33x, 'Recompute the GreenF matrix for the H_T Hamiltonian!')")
      end if
      call InitGrFctMat()
!________________________________________________________________________________________ 	  
!_________________ (2) The initial GreenF matrix for H_T Hamiltonian ____________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         write(*, "(33x, 'Reset the \Delta matrix used for updating auxiliary fields!')")
      end if
      call ResetDeltabU()
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for warm up part process __________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&     
		call system_clock(time2)
      TimeSgBIN = TimeIntrvl(time1, time2)
!**************************************************************************************************	  
!___________________ 4. Output consumed time of this starting part of CPMC ________________________
!**************************************************************************************************
      call SaveCalTim(NB)
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%
!________________________ Call MPI_Barrier to make all the processes synchronous ______________________________
!v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%&v%      
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)	
#endif
		
   end subroutine CPMCOccMuT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ComptnTMuT()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ComptnTMuT()
! TYPE:     subroutine
! PURPOSE:  This subroutine first measures the total electron density and then use it to determine the ChemP_BT.
! KEYWORDS: Measure nT and get ChemP_BT.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: First measures the total electron density, and then determine the ChemP_BT as n_T = n.
!
!     Input: (none); Output: (none)
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
      real(rp) nT_BT, Rtp0
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations electron density and sign ____________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Compute the total electron density ________________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Compute the weight summation of all measurements _________________
!________________________________________________________________________________________
      WtMeanSumBIN = sum(WeightList(1:NObsStat))
!________________________________________________________________________________________ 	  
!_________________ (1) Reweighting of all measurements __________________________________
!________________________________________________________________________________________
      nT_Now = dot_product(WeightList(1:NObsStat), n_Occ_List(1:NObsStat)) / WtMeanSumBIN
!**************************************************************************************************	  
!___________________ 1. Determine the ChemP_BT parameter for H_T Hamiltonian ______________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Compute ChemP_BT by bisection method _____________________________
!________________________________________________________________________________________
      call SelfTuneChemP_BT(nT_Now, 1.0E-5_rp, nT_BT)
!________________________________________________________________________________________ 	  
!_________________ (1) Get the average of ChemP_BT from all processes ___________________
!________________________________________________________________________________________
#ifdef MPIPROCESS
      !!!!!!!!!! Average of ChemP_BT
      call MPI_Barrier(acomm, ierr)
      call MPI_ALLREDUCE(ChemP_BT, Rtp0, 1, rp_MPI_REAL, MPI_SUM, acomm, ierr)
      ChemP_BT = Rtp0 / dble(anprc)
      !!!!!!!!!! Recompute nT_BT
      call Calculate_nTOfHT(nT_BT)
#endif
!________________________________________________________________________________________ 	  
!_________________ (2) Output the results of ChemP_BT and nT_BT _________________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         write(*, "(38x, 'ChemP   , nT_Bx = ', sp, es19.12, ss, es20.12)") ChemP   , nT_Now
         write(*, "(38x, 'ChemP_BT, nT_BT = ', sp, es19.12, ss, es20.12)") ChemP_BT, nT_BT
      end if

   end subroutine ComptnTMuT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$