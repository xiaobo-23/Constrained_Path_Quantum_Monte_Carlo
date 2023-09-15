!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A module and a few subroutines used for Recording the calculation time consumed by the program. 
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!          TimeRecord--Module to define the related quantities uded in the recording of time;
!          TimeRecordInit()--A small subroutine used to perform the Initializations of the TimeRecord Module.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Module ______________________________________________________________________
!________________________________________________________________________________________________________________________
   module TimeRecord
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this Module ________________________________________
!______________________________________________________________________________________________________________   
      implicit none
!______________________________________________________________________________________________________________	  
!_________________________________ Some Quantities used in the program ________________________________________
!______________________________________________________________________________________________________________      
      integer(8) TimRat      ! The Computer counting rate For the Integer*8 kind of number
      integer(8) TimMax      ! The Maximum counting rate For the Integer*8 kind of number, then it overflowed
      real(8) TimThh         ! The time used from zero to overflowed (Unit : second)

      integer, parameter :: DatTimLen = 32  ! The length of the time recording
   
      character(8) ccyymmdd  ! Date, ccyy --> year(1998); mm --> Month(1-12); dd --> day(0-31).
      character(10) hhmmss   ! hhmmss.sss: hh --> hour(0-23); mm --> minute(0-59); ss.sss --> second and millisecond.
      character(19) DatTim   ! Used to store the Present time as : ccyymmdd hhmmss and 5 spacing
!______________________________________________________________________________________________________________	  
!_____________________ Module Containing A Function used to calculate the time Interval  ______________________
!______________________________________________________________________________________________________________ 	  
      contains
!______________________________________________________________________________________________________________	  
!__________ Module Containing A Subroutine used to Calculate the time Interval between two times  _____________
!______________________________________________________________________________________________________________ 
         function TimeIntrvl(time1, time2)
     
            implicit none
   
            integer(8) time1, time2, time0
            real(8) TimeIntrvl
      
            time0 = time2 - time1                         ! Computer Counting times
            if(time0 < 0) time0 = time0 + TimMax + 1      ! if smaller
            TimeIntrvl = time0 * 1._8 / TimRat            ! Unit in second

         end function TimeIntrvl          
!______________________________________________________________________________________________________________	  
!_____________________ Module Containing A Subroutine used to present Date and Time  __________________________
!______________________________________________________________________________________________________________ 
         subroutine PresentMoment()
     
            implicit none
      
            call date_and_time(ccyymmdd, hhmmss)                     ! Catch the present moment time
            write(DatTim, 11) ccyymmdd(1 : 4), ccyymmdd(5 : 6), ccyymmdd(7 : 8), hhmmss(1 : 2), hhmmss(3 : 4), hhmmss(5 : 6)
11		      format(A4, "-", A2, "-", A2, " ", A2, ":", A2, ":", A2)  ! Write the time into Character DatTim

         end subroutine PresentMoment
!______________________________________________________________________________________________________________	  
!________________ Module Containing A Function used to Transform the real second into Integer  ________________
!______________________________________________________________________________________________________________
         function IntSec(Second)
            
            implicit none
            
            integer IntSec
            real(8) Second
            
            IntSec = int(Second)
            
         end function IntSec
!______________________________________________________________________________________________________________	  
!________________ Module Containing A Function used to Transform Real Second into minute  _____________________
!______________________________________________________________________________________________________________
         function Sec2Min(Second)
            
            implicit none
            
            real(8) Second, Sec2Min
            
            Sec2Min = Second / 60
            
         end function Sec2Min                   
!______________________________________________________________________________________________________________	  
!________________ Module Containing A Function used to Transform Real Second into Hour  _______________________
!______________________________________________________________________________________________________________
         function Sec2Hour(Second)
            
            implicit none
            
            real(8) Second, Sec2Hour
            
            Sec2Hour = Second / 3600
            
         end function Sec2Hour         
   
   end module TimeRecord
!________________________________________________________________________________________________________________________  
!____________________________________ End Module ________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine TimeRecordInit()  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  TimeRecordInit() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine Give Some Initializations to the TimeRecord module.
! KEYWORDS: Initialization of the TimeRecord module
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: We give the values to the quanties defined in the TimeRecord RealPrecsn. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this Module ________________________________________
!______________________________________________________________________________________________________________      
      use TimeRecord
      implicit none
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________      
      integer(8) time0
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________      
      call system_clock(time0, TimRat, TimMax)
      
      TimThh = TimMax * 1.0_8 / TimRat

   end subroutine TimeRecordInit
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$