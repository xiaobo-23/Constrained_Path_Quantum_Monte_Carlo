!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: This program performs the Finite-temperature CPMC calculations of Hubbard Model on Square lattice including: 
!               (0) Nearest-Neighbor hopping terms  --> Hopt1 (with plaquette-staggered, same for up and down spins)
!               (1) Nearest-Neighbor hopping terms  --> Hopt2 (opposite for spin-up and spin-down)
!               (2) On-site Hubbard Repulsion term  --> HubbU
!               (3) Doping chemical potential
!          Notice: In the calculations, we introduce the on-site 2-component Ising fields for the U interaction term.
!                  Including four kinds of HS transformations: (0) CDW channel; (1) Sz channel; (2) Sx channel; 
!                         (3) Sy channel.
!
!          This program has taken care of general phase problem, which is the corresponding thing in complex model 
!                Hamiltonian system as the extension of sign problem in real model Hamiltonian system. 
!
!          We can now simulate the Hubbard model away from half-filling. And this program applies the symmetric 
!                 Trotter Decomposition as   
!             exp(-\Delta\tau\hat{H}) = exp(-\Delta\tau\hat{K}/2) * exp(-\Delta\tau\hat{V}) * exp(-\Delta\tau\hat{K}/2)
!  
!          In this program, the H_T in B_T matrix uses the restricted Hartree-Fock mean-field, which is actually
!               just an effective chemical potential. 
!
!          Here, we write the model Hamiltonian as following:
!            (0) For B_T, it's H_T = H_0 + \mu_T \sum_i (n_{i\up} + n_{i\dw})
!            (1) For B_X, it's H   = H_0 + H_U = H_T + H_U + (H_0-H_T) = H_T + H_U - \mu_T \sum_i (n_{i\up} + n_{i\dw})
!
!          In this program, we want to fix the doping chemical potential ChemP. Then we adjust the \mu_T 
!               self-consistently to make the total electron density of the simulations converge.
!
! COMMENT: Used for total CPMC calculations.   
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! TARGET: The following physical properties of 2D Hubbard model on Square Lattice are calculated:
!         (0) Energies; (1) Static correlation functions; (2) Dynamic correlation functions.
!
!         Containing following programs or subroutines/functions:
!
!    main --> The main program of this CPMC code.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Main Program ________________________________________________________________
!________________________________________________________________________________________________________________________
      program Main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  Main
! TYPE:     Main Program
! PURPOSE:  This Subroutine performs the whole CPMC simulation.
! KEYWORDS: Total CPMC calculations.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     By Monte Carlo sampling of the auxiliary fields, all the physical properties are calculated.
!
!     Input: (none);   Outpt: (none). 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this Module ________________________________________
!______________________________________________________________________________________________________________
#ifdef OMPTHREADS
      use OMP_LIB
#endif
            implicit none       
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
            integer(8) time1 
            integer(8) time2
!**************************************************************************************************************         
!________________________________ Set MKL_NUM_THREADS and OMP_NUM_THREADS _____________________________________
!**************************************************************************************************************
#ifdef OMPTHREADS
      integer OMPThreads
      integer MKLThreads
      open(57, err = 58, file = "Input/NumThreads.txt")
      read(57, *) OMPThreads
      read(57, *) MKLThreads
      close(57)
      call OMP_set_num_threads(OMPThreads)
      call MKL_set_num_threads(MKLThreads)
      go to 59
      
58    continue
      call OMP_set_num_threads(4)
      call MKL_set_num_threads(4)
      
59    continue
#endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
!_________________________ 0. The initialization process for some modules ________________________
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            call MPISettingInit()
            call TimeRecordInit()
            call StdInOutStInit()
            call RealPrecsnInit()      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
!_________________________ 1. The Whole Constraint-Path Monte Carlo simulation ___________________
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!___________________________________________________________________________________        
!_________________ (0) Time for starting point of CPMC pragram ___________________
!___________________________________________________________________________________
            call SubroutineBgn("CPMC", 1, time1)
!___________________________________________________________________________________        
!_________________ (1) Perform the CPMC simulation _______________________________
!___________________________________________________________________________________   
            call CPMC()
!___________________________________________________________________________________        
!_________________ (2) Time for Ending point of CPMC pragram _____________________
!___________________________________________________________________________________              
            call SubroutineEnd("CPMC", 1, time1, time2)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
!_________________________ 2. Finalization of MPI setting in the program _________________________
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call MPISettingFinl()
            
   end program Main
!________________________________________________________________________________________________________________________  
!____________________________________ End Main Program __________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
