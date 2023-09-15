!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform some simple initializations for CPMC simulations.
! COMMENT: Simple initializations.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!   Includes: (0) Initialization of random number generator;
!             (1) List all of the defined preprocessor conditions.
!             
!   InitRandPrep --> Subroutine used to perform the simple initializations for CPMC simulations;
!
!   RandNumberSeeds --> Subroutine used to generate the random number seeds for all processes;
!   PreprocessCheck --> Subroutine used to check preprocessor definitions and output them.   
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine InitRandPrep()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitRandPrep() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the initialization for the random number generator, as well as 
!                    some preprocessors.
! KEYWORDS: Initialization of Rand and Preprocessors.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Initialization for spring RAND, defining new MPI data type for population control, and preprocessors. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use RandomNumb
      use CoreParamt
      use Observable
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of First Initialization ____________________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of initialization process _______________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'InitRandPrep: Initializations of Rand, Folder and preprocessors!')")
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	  
!_________________ 0. Initialize the Random number generator ______________________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call RandNumberSeeds(.false.)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	  
!_________________ 1. Initialize all the preprocessors used in the simulations ____________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call PreprocessCheck()

   end subroutine InitRandPrep
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine RandNumberSeeds(IfReload)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RandNumberSeeds(IfReload) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to generate the seeds for random numbers in the program, and initiate the 
!              random number generator.
! KEYWORDS: Initialization of Random numbers seeds.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: The seeds of random number should be different for different time, as well as for different processes. 
!              And we should also allow input of these random seeds, and we should also store them for future use. 
!
!     Input:  IfReload --> == T fo reloading seeds and initiate the random number generator again.
!
!     Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use RandomNumb
		use MPISetting
      implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfReload
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer time
      integer ISeed
      character(20) ProcessID
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Generating random numebr seeds __________________________
!______________________________________________________________________________________________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	  
!___________________ 0. If there is input files, simply read the seeds ____________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!________________________________________________________________________________________ 	  
!_________________ (0) Read the seeds of random number generator ________________________
!________________________________________________________________________________________
      write(ProcessID, "(I3.3)") amyid
      open(157, err = 377, file = "Output/RandNumbSeeds/RandSeed_" // Trim(ProcessID) // ".txt", status = "old")
      read(157, *) ISeed
      close(157)
!________________________________________________________________________________________ 	  
!_________________ (1) Initiate random number generator for all processes _______________
!________________________________________________________________________________________
      call spring_sfmt_init(ISeed)
!________________________________________________________________________________________ 	  
!_________________ (2) Output information of reading seeds and go to 457 ________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         if(.not. IfReload) then
            write(*, "(30x, 'Read the random number seeds from input files!')")
         else
            write(*, "(28x, 'Reload random number seeds from InitRandPrep subroutine!')")
         end if
      end if
      go to 457
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	  
!___________________ 1. If there is no input files, generate the seeds ____________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
377   continue
!________________________________________________________________________________________ 	  
!_________________ (0) Generate the seed with dependence on time and processId __________
!________________________________________________________________________________________
      call SYSTEM_CLOCK(time)
      ISeed = nint( abs( time - ( (10*amyid+2020) * 1991 + 2019 ) * (95049+amyid) ) * (sqrt(5.0_rp)-1.0_rp)/2.0_rp )
!________________________________________________________________________________________ 	  
!_________________ (1) Initiate random number generator for all processes _______________
!________________________________________________________________________________________
      call spring_sfmt_init(ISeed)
!________________________________________________________________________________________ 	  
!_________________ (2) Output information of reading seeds ______________________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         write(*, "(30x, 'Generate different random number seeds for processes!')")
      end if
!________________________________________________________________________________________ 	  
!_________________ (3) Store the random seeds for all processes _________________________
!________________________________________________________________________________________ 
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)	
#endif
      write(ProcessID, "(I3.3)") amyid
      open( 157, file = "Output/RandNumbSeeds/RandSeed_" // Trim(ProcessID) // ".txt")
      write(157, *) ISeed
      close(157)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&	  
!___________________ 2. If successfully read seeds from input files, get here from 0 ______________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
457   continue
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)	
#endif

   end subroutine RandNumberSeeds
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine PreprocessCheck()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PreprocessCheck() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to check the definitions of preprocessors and output them.
! KEYWORDS: Initialization of Preprocessors.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Check the definitions of preprocessors and output them. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		use CoreParamt
      use Observable
      use MPISetting
#ifdef OMPTHREADS
      use OMP_LIB
#endif
		implicit none
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of preprocessor check ______________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________ 	  
!_________________ (00) Output whether the code is sequential or parallel _______________
!________________________________________________________________________________________      
      if(amyid == amstr) then
#ifdef MPIPROCESS
         write(*, "(30x, 'MPIPROCESS = T --> MPI parallel finite-T CPMC simulation with anprc = ', I4)") anprc
#else
         write(*, "(30x, 'MPIPROCESS = F --> Sequential finite-T CPMC simulation!')")
#endif   
      end if
!________________________________________________________________________________________ 	  
!_________________ (01) Output whether OpenMP has been activated ________________________
!________________________________________________________________________________________ 
#ifdef OMPTHREADS
   !$OMP PARALLEL
      OMPNumThread = OMP_GET_NUM_THREADS()
   !$OMP END PARALLEL
#endif

      if(amyid == amstr) then
#ifdef OMPTHREADS
         write(*, "(30x, 'OMPTHREADS = T --> OpenMP has been opened with OMPThreads = ', I4)") OMPNumThread
#else
         write(*, "(30x, 'OMPTHREADS = F --> OpenMP has not been activated!')")
#endif
      end if
!________________________________________________________________________________________ 	  
!_________________ (02) Output whether to use the symmetric Trotter Decomps _____________
!________________________________________________________________________________________
      if(amyid == amstr) then
         if(SymTrotDcp) then
            write(*, "(30x, 'SymTrotDcp = T --> Symmetric Trotter decomposition for exp(-dt*H)!')")
         else
            write(*, "(30x, 'SymTrotDcp = F --> Asymmetric Trotter decomposition for exp(-dt*H)!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (03) Output methods to multiply exp(+/-dt*H_0) or exp(+/-dt*H_T) _____
!________________________________________________________________________________________
      if(amyid == amstr) then
         if(FFTEXPDTH0) then
            write(*, "(30x, 'FFTEXPDTH0 = T --> FFT method to calculate exp(+/-dt*H_0)*Matrix!')")
         else
            write(*, "(30x, 'FFTEXPDTH0 = F --> DGEMM method to calculate exp(+/-dt*H_0)*Matrix!')")
         end if
         if(FFTEXPDTHT) then
            write(*, "(30x, 'FFTEXPDTHT = T --> FFT method to calculate exp(+/-dt*H_T)*Matrix!')")
         else
            write(*, "(30x, 'FFTEXPDTHT = F --> DGEMM method to calculate exp(+/-dt*H_T)*Matrix!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (04) Output the measurement methods for energies _____________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         if(IfFftEnPar) then
            write(*, "(30x, 'IfFftEnPar = T --> Measure energies and pairing matrix using k-space GreenF!')")
         else
            write(*, "(30x, 'IfFftEnPar = F --> Measure energies and correlations using r-space GreenF!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (05) Output the dispersion type used in this simulation ______________
!________________________________________________________________________________________
      if(amyid == amstr) then
         if(EkDispType == 0) then
            write(*, "(30x, 'EkDispType = 0 --> The standard Hubbard band dispersion is used!')")
         else if(EkDispType == 1) then
            write(*, "(30x, 'EkDispType = 1 --> The Hubbard band dispersion with 4t shift is used!')")
         else if(EkDispType == 2) then
            write(*, "(30x, 'EkDispType = 2 --> The standard Quadratic band dispersion is used!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (06) Output the method of numerical stablization _____________________
!________________________________________________________________________________________
      if(amyid == amstr) then
#ifdef QRSTABLIZE
         write(*, "(30x, 'QRSTABLIZE = T --> QR algorithm for numerical stablization!')")
#else
         write(*, "(30x, 'QRSTABLIZE = F --> SVD algorithm for numerical stablization!')")
#endif
      end if
!________________________________________________________________________________________ 	  
!_________________ (07) Output the updating method for auxiliary fields _________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         if(UpdtMethod == 0) then
            write(*, "(30x, 'UpdtMethod = 0 --> Local update method for choosing auxiliary fields!')")
         else if(UpdtMethod == 1) then
            write(*, "(30x, 'UpdtMethod = 1 --> Delayed update method for choosing auxiliary fields!')")
         else if(UpdtMethod == 2) then
            write(*, "(30x, 'UpdtMethod = 2 --> Force-bias update method for choosing auxiliary fields!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (08) Canonical or Grand Canonical ensemble calculation _______________
!________________________________________________________________________________________      
      if(amyid == amstr) then
         if(IfFixnT) then
            write(*, "(30x, 'IfFixnT    = T --> QMC simulations with fixed total electron density!')")
         else
            write(*, "(30x, 'IfFixnT    = F --> QMC simulations with fixed chemical potential!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (09) Perform measurements with types of OPEN or PERIODIC _____________
!______________________ Boundary conditions _____________________________________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         if(IfPyObsPBC) then
            write(*, "(30x, 'IfPyObsPBC = T --> QMC measurements with PERIODIC boundary conditions!')")
         else
            write(*, "(30x, 'IfPyObsPBC = F --> QMC measurements with OPEN     boundary conditions!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (10) Output whether to perform additional measure ____________________
!______________________ in [BetaT, 0] sweep _____________________________________________
!________________________________________________________________________________________  
      if(amyid == amstr) then
         if(IfM2OneMea) then
            write(*, "(30x, 'IfM2OneMea = T --> Perform additional measurements in [BetaT, 0] sweep!')")
         else
            write(*, "(30x, 'IfM2OneMea = F --> No sweep and measurements in [BetaT, 0] path!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (11) Output whether to perform additional weighting __________________
!______________________ for all the sweep measures in a single BIN ______________________
!________________________________________________________________________________________  
      if(amyid == amstr) then
         if(IfSwepReWt) then
            write(*, "(30x, 'IfSwepReWt = T --> Perform additional weighting for sweep measures!')")
         else
            write(*, "(30x, 'IfSwepReWt = F --> No additional weighting for sweep measures!')")
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (12) Output whether to perform the dynamic measurements ______________
!________________________________________________________________________________________  
      if(amyid == amstr) then
         if(IfTAU) then
            write(*, "(30x, 'IfTAU      = T --> Do calculate Dynamic correlation functions!')")
            if(IfEqDistDt) then
               write(*, "(30x, '               --> Equal distances between adjacent Tau points!')")
            else
               write(*, "(30x, '               --> Unequal distances between adjacent Tau points!')")
            end if
            if(IfTau0Rand) then
               write(*, "(30x, '               --> Randomly choose NT\in[0, LTrot-1] as Tau==0 point!')")
            else
               write(*, "(30x, '               --> Simply choose NT==0 as Tau==0 point!')")
            end if
            if(IfDyGrFOut) then
               write(*, "(30x, '               --> Output dynamic single-particle GreenF       at all k points!')")
            end if
            if(IfDySpnOut) then
               write(*, "(30x, '               --> Output dynamic spin-spin       correlations at all k points!')")
            end if
            if(IfDyDenOut) then
               write(*, "(30x, '               --> Output dynamic density-density correlations at all k points!')")
            end if
            if(IfDyPstOut) then
               write(*, "(30x, '               --> Output dynamic pair-pair       correlations at all k points!')")
            end if
            if(IfDyDWvOut) then
               write(*, "(30x, '               --> Output dynamic dwave-dwave     correlations at all k points!')")
            end if
            if(IfDyCurrnt) then
               write(*, "(30x, '               --> Output dynamic current-current correlations at all k points!')")
            end if
         else
            write(*, "(30x, 'IfTAU      = F --> Do not calculate Dynamic correlation functions!')")
         end if
      end if

   end subroutine PreprocessCheck
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$