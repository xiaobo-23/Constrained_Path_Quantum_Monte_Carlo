!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A few subroutines used for calculating the time used for the subroutine calculations.
! COMMENT: Time Recording for (C)DMFT calculations.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   SubroutineSkp --> Subroutine to record the information if some subroutine is skipped in the calculations;
!   SubroutineBgn --> Subroutine to record the information if some subroutine begins in the calculations;
!   SubroutineEnd --> Subroutine to record the information if some subroutine ends in the calculations.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________                            
	subroutine SubroutineSkp(CaptionName, EmptyLen)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SubroutineSkp(CaptionName, EmptyLen)
! TYPE:     subroutine
! PURPOSE:  This Subroutine store the information if some subroutine in the whole calculations is skipped.
! KEYWORDS: Calculation information recording.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input: CaptionName --> The name of the subroutine that has been skipped;
!            EmptyLen    --> The length of the spacing for the output of the information.
!
!     Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this Module ________________________________________
!______________________________________________________________________________________________________________
		use StdInOutSt
		use TimeRecord
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		character(*) CaptionName  ! The name of the subroutine that has been skipped
		integer EmptyLen          ! The length of the spacing for the output of the information
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      character(100) FileName
!______________________________________________________________________________________________________________	  
!______________________________________ Storing the information _______________________________________________
!______________________________________________________________________________________________________________
		if(amyid == amstr) then
			call PresentMoment()
			open( 50, file = trim(FLog), access = "append")
         write(FileName, *) EmptyLen
         write(50, "(A, " // adjustl(FileName) // "x, '  > ', A, ' Skipped')") DatTim, CaptionName
			close(50)
		end if

	end subroutine SubroutineSkp 
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________                            
	subroutine SubroutineBgn(CaptionName, EmptyLen, time1)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SubroutineBgn(CaptionName, EmptyLen, time1)
! TYPE:     subroutine
! PURPOSE:  This Subroutine store the information if some subroutine begins in the calculations.
! KEYWORDS: Calculation information recording.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input: CaptionName --> The name of the subroutine that has bagan to calculate;
!            EmptyLen    --> The length of the spacing for the output of the information;
!            time1       --> The exact time that the subroutine begins to calculate.
!
!     Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this Module ________________________________________
!______________________________________________________________________________________________________________
		use StdInOutSt
		use TimeRecord
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		character(*) CaptionName  ! The name of the subroutine that has been skipped
		integer EmptyLen          ! The length of the spacing for the output of the information
		integer(8) time1          ! The exact time that the subroutine begins to calculate.
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      character(100) FileName
!______________________________________________________________________________________________________________	  
!______________________________________ Storing the information _______________________________________________
!______________________________________________________________________________________________________________
		if(amyid == amstr) then
			call PresentMoment()
			open( 50, file = Trim(FLog), access = "append")
         write(FileName, *) EmptyLen
         write(50, "(A, " // adjustl(FileName) // "x, '  + ', A, ' Bgn')") DatTim, CaptionName
			close(50)
		end if
		
		call system_clock(time1)

	end subroutine SubroutineBgn
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________                            
	subroutine SubroutineEnd(CaptionName, EmptyLen, time1, time2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SubroutineEnd(CaptionName, EmptyLen, time1, time2)
! TYPE:     subroutine
! PURPOSE:  This Subroutine store the information if some subroutine ends in the calculations and it calculate the 
!              time that it costs  in the calculation as TimeIntrvl(time1, time2)  -->  TimeIntrvl(time1, time2) / 3600 .
! KEYWORDS: Calculation information recording.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input: CaptionName --> The name of the subroutine that has bagan to calculate;
!            EmptyLen    --> The length of the spacing for the output of the information;
!            time1       --> The exact time that the subroutine begins to calculate.
!            time2       --> The exact time that the subroutine ends to calculate.
!
!     Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this Module ________________________________________
!______________________________________________________________________________________________________________
		use StdInOutSt
		use TimeRecord
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		character(*) CaptionName  ! The name of the subroutine that has been skipped
		integer EmptyLen          ! The length of the spacing for the output of the information
		integer(8) time1          ! The exact time that the subroutine begins to calculate.
		integer(8) time2          ! The exact time that the subroutine ends to calculate.
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      character(100) FileName
!______________________________________________________________________________________________________________	  
!______________________________________ Storing the information _______________________________________________
!______________________________________________________________________________________________________________
		call system_clock(time2)
		
		if(amyid == amstr) then
			call PresentMoment()
			open( 50, file = trim(FLog), access = "append")
         write(FileName, *) EmptyLen
         write(50, &
    & "(A, " // adjustl(FileName) // "x, '  - ', A, ' End, Elapsed Time = ', F14.5, 's = ', F12.5, 'm = ', F10.5, 'h')") &
    & DatTim, CaptionName, TimeIntrvl(time1, time2), TimeIntrvl(time1, time2)/60, TimeIntrvl(time1, time2)/3600
			close(50)
		end if

   end subroutine SubroutineEnd	 
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$