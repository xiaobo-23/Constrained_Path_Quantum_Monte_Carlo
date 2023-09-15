!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL SPIN PINNING FIELDS; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform Initialization of CPMC output files before the whole CPMC simulations.
! COMMENT: CPMC Initialization process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   InitMCOutput --> Subroutine to Initialization of CPMC output files used in CPMC;
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine InitMCOutput()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitMCOutput()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the initialization for all output files handles for CPMC data. 
! KEYWORDS: Initialization of CPMC output file handles.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Initialization of CPMC output file handles, including:
!             (0) Initializations of The output files Handles integer;
!             (1) Open the output files for the CPMC calculations;
!             (2) Output the CPMC calculation information.
!
!     Input:  (none)   Output: (none)
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
!_____________________________ Main calculations of CPMC output file handles __________________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of initialization process _______________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'InitMCOutput: Some output settings for output files!')")
      end if
!**************************************************************************************************  
!______________________________ 0. Output the information for the CPMC calculations _______________
!************************************************************************************************** 
      if(amyid == amstr) then
!________________________________________________________________________________________  
!_________________ (0) The energies and energy derivatives over parameters ______________
!________________________________________________________________________________________
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         open( 283, file = "Output/01_QMCExpectEnergy.txt", access = "append") 
         write(283, "('BIN', 14x, 'EHopt1', 21x, 'EHopt2', 21x, 'EHopt3', 21x, 'EZmFld', 21x, 'EPinSz', 21x, &
            & 'EDopCh', 21x, 'EHubbU', 21x, 'ESinusoidalPinSz', 21x, 'ETotal')")
         close(283)

         open( 283, file = "Output/01_FstOrdDiffParam.txt", access = "append")
         write(283, "('BIN', 14x, 'EHbUCh', 21x, 'ETotCh', 24x, 'HmOvt2', 21x, 'HmOvt3', 21x, 'HmOvZm', 21x, &
            & 'HmOvbU', 21x, 'HmOvSinusoidalPinSz')")
         close(283)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         open( 283, file = "Output/01_OccDouSpnCorrFc.txt", access = "append")
         write(283, "('BIN', 14x, 'nOcpUp', 21x, 'nOcpDw', 21x, 'nOcpTt', 21x, 'Num_Ne', 24x, 'DouOcc', 21x, &
            & 'NNDCrF')")
         close(283)
         
         if(IfM2OneMea) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            open( 283, file = "Add_Output/01_QMCExpectEnergy_Add.txt", access = "append") 
            write(283, "('BIN', 14x, 'EHopt1', 21x, 'EHopt2', 21x, 'EHopt3', 21x, 'EZmFld', 21x, 'EPinSz', 21x, &
               & 'EDopCh', 21x, 'EHubbU', 21x, 'ESinusoidalPinSz', 21x, 'ETotal')")
            close(283)

            open( 283, file = "Add_Output/01_FstOrdDiffParam_Add.txt", access = "append")
            write(283, "('BIN', 14x, 'EHbUCh', 21x, 'ETotCh', 24x, 'HmOvt2', 21x, 'HmOvt3', 21x, 'HmOvZm', 21x, &
               & 'HmOvbU', 21x, 'HmOvSinusoidalPinSz')")
            close(283)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            open( 283, file = "Add_Output/01_OccDouSpnCorrFc_Add.txt", access = "append")
            write(283, "('BIN', 14x, 'nOcpUp', 21x, 'nOcpDw', 21x, 'nOcpTt', 21x, 'Num_Ne', 24x, 'DouOcc', 21x, &
               & 'NNDCrF')")
            close(283)
         end if
!________________________________________________________________________________________  
!_________________ (1) Information of consumed time of the simulation ___________________
!________________________________________________________________________________________         
         open( 291, file = "Output/00_CPMC_Time_Count.txt", access = "append")
         write(291, "('Output the information of consumed time of the simulation!')")
         close(291)
!________________________________________________________________________________________  
!_________________ (2) Information of numerical stablization process ____________________
!________________________________________________________________________________________
         open( 291 , file = "Output/00_NmStbGrnFStaDyn.txt", access = "append")
         write(291, "('Output the GrF difference of the numerical stablizarion process!')")
         close(291)
!________________________________________________________________________________________  
!_________________ (3) Information of population control process ________________________
!________________________________________________________________________________________
         open( 291 , file = "Output/00_PopulationCntrl.txt", access = "append")
         write(291, "('Output information of the population control process!')")
         close(291)
!________________________________________________________________________________________  
!_________________ (4) Information of adjusting growth estimator ________________________
!________________________________________________________________________________________
         open( 291 , file = "Output/00_Adjust_Constant.txt", access = "append")
         write(291, "('Output information of adjusting the constant during propagation!')")
         close(291)
      end if
		
   end subroutine InitMCOutput
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$