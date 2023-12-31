!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A few subroutines for performing the main part of CPMC simulation as updating process, including:
!                  (0) First calculate the Green's function as a start;
!                  (1) Perform the standard Sweep and measurement in CPMC simulation.
! COMMENT: CPMC Simulation process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   CPMCSimulate --> Subroutine to call subroutines to perform the total CPMC simulations;
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine CPMCSimulate()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCSimulate()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform main CPMC simulation.
! KEYWORDS: CPMC Simulation process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: CPMC Simulation:
!                  (0) Calculate the Green Function for time slice=M;
!                  (1) Perform the warm up for the CPMC siumulations;
!                  (2) Perform the standard Sweep and measurement in CPMC simulation.
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
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
!______________________________________________________________________________________________________________     
!_______________________________ Main calculations of CPMC Simulation _________________________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCSimulate process _________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "(A)") "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         write(*, "(2x, 'CPMCSimulate: The core QMC simulation!')")
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     
!___________________ 0. Start CPMC simulation by Getting H_T properties ___________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call SubroutineBgn("CPMCSmuBgn", 11, time1)
      call CPMCSmuBgn()
      call SubroutineEnd("CPMCSmuBgn", 11, time1, time2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     
!___________________ 1. Determine ChemP_BT or ChemP self-consistently _____________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!________________________________________________________________________________________         
!_________________ (0) For IfMuTqmcNt == T case, determine ChemP_BT _____________________
!________________________________________________________________________________________
      if(IfMuTqmcNt) then
         call SubroutineBgn("CPMCGetMuT", 11, time1)
         call CPMCGetMuT()
         call SubroutineEnd("CPMCGetMuT", 11, time1, time2)
!________________________________________________________________________________________         
!_________________ (1) For IfFixnT == T case, determine ChemP ___________________________
!________________________________________________________________________________________
      else if(IfFixnT) then
         call SubroutineBgn("CPMCFixdnT", 11, time1)
         call CPMCFixdnT() 
         call SubroutineEnd("CPMCFixdnT", 11, time1, time2)
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     
!___________________ 2. Warm up to achieve the equilibrium of the system __________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&         
      call SubroutineBgn("CPMCWarmUp", 11, time1)
      call CPMCWarmUp() 
      call SubroutineEnd("CPMCWarmUp", 11, time1, time2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&     
!___________________ 3. Sweep and measure the physical observables ________________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call SubroutineBgn("CPMCSwpMea", 11, time1)
            call CPMCSwpMea() 
      call SubroutineEnd("CPMCSwpMea", 11, time1, time2)
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCSimulate process _________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(      
      if(amyid == amstr) then
         write(*, "(A)") "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         write(*, "()")
         write(*, "()")
      end if
            
      end subroutine CPMCSimulate
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$