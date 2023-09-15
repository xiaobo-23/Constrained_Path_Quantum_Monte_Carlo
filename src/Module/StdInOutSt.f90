!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A module and a few subroutines used for setting Standard Input/Output process.
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He, Rong-Qiang He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   StdInOutSt        --> Module to define the related characters used to write the Error and Warning Messages;
!   StdInOutStInit()  --> A small subroutine to perform Initializations of quantities in Module StdInOutSt;
!   ErrOut(AMsg) --> Subroutine to write the Error Message into the Error output file;
!   WrnOut(AMsg) --> Subroutine to write the Warning Message into the Warning output file.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Module ______________________________________________________________________
!________________________________________________________________________________________________________________________
	module StdInOutSt
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________    
		implicit none
!______________________________________________________________________________________________________________	  
!_____________________________ The integer handle controling the Monitor outputing ____________________________
!______________________________________________________________________________________________________________		
		integer MonitorOt       ! The monitor output process
!______________________________________________________________________________________________________________	  
!_____________________________________ Output Folder Related Quantities _______________________________________
!______________________________________________________________________________________________________________       
		character(64)   FPGI		 ! Prefix for general input file name
		character(64)   FPGO    ! Prefix for output file name
		character(64)   FPIO		 ! Prefix for output file name for specific iteration
		character(64)   GLPI    ! The input file folder name contanining Gauss-Legendre Integral points
!______________________________________________________________________________________________________________	  
!_____________________________ Output file containing all the calculation details _____________________________
!______________________________________________________________________________________________________________       
		character(64)   FLog		 ! Calculation information output file name
!______________________________________________________________________________________________________________	  
!_____________________________ Output file containing Error, Warning Messages__________________________________
!______________________________________________________________________________________________________________       
		character(64)   FErr		 !   error output file name
		character(64)   FWrn		 ! warning output file name
		character(64)   FMnt		 ! The video output
	  
		character(8192) ErrMsg	 !   error output massage to file ErrMsg
		character(8192) WrnMsg	 ! warning output massage to file WrnMsg
		character(8192) MntMsg	 ! monitor output massage to the Monitor	  
	
	end module StdInOutSt
!________________________________________________________________________________________________________________________  
!____________________________________ End Module ________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine StdInOutStInit()  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  StdInOutStInit() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine Give Some Initializations to the StdInOutSt module.
! KEYWORDS: Initialization of the StdInOutSt module
! AUTHOR:   Yuan-Yao He, Rong-Qiang He
! TIME:     2020-02-27
! DESCRIPTION: We give the values to the quanties defined in the StdInOutSt RealPrecsn. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use StdInOutSt
		use TimeRecord
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________ Output Folder name Related Quantities ______________________________________
!______________________________________________________________________________________________________________ 
		FPGI = "Input/I"                                        ! Prefix for general input file name
		FPGO = "Output/O"                                       ! Prefix for output file name
		FPIO = "N/A"                                            ! Prefix for output file name for specific iteration
		GLPI = "Input/Gauss-Legendre-Points/"                   ! The input file folder name contanining Gauss-Legendre Integral points on my computer
!      GLPI = "/home/yyhe/work/Tables/Gauss-Legendre-Points/" ! The input file folder name contanining Gauss-Legendre Integral points on my Server
!______________________________________________________________________________________________________________	  
!____________________________ Error, Warning and Calculation output file name _________________________________
!______________________________________________________________________________________________________________       
		call date_and_time(ccyymmdd, hhmmss)                                                     ! Present Time
		write(FLog, "(A16, A8, A1, A10, A4)") "Output/CPMC.Log.", ccyymmdd, ".", hhmmss, ".txt"  ! File name  --> Total calculations
		write(FWrn, "(A16, A8, A1, A10, A4)") "Output/CPMC.Wrn.", ccyymmdd, ".", hhmmss, ".txt"  ! File name  --> Warning message
		write(FErr, "(A16, A8, A1, A10, A4)") "Output/CPMC.Err.", ccyymmdd, ".", hhmmss, ".txt"  ! File name  --> Error Message
		write(FMnt, "(A16, A8, A1, A10, A4)") "Output/CPMC.Mnt.", ccyymmdd, ".", hhmmss, ".txt"  ! File name  --> Normal output message
!________________________________________________________________________________________	  
!_______ (0) Monitor Output for this whole program ______________________________________
!________________________________________________________________________________________
		MonitorOt = 229

	end subroutine StdInOutStInit
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$