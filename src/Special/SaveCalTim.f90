!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine SaveCalTim(NB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SaveCalTim(NB)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to output the calling times and consumed calculation time by all different 
!              parts in a single BIN QMC simulations, including the following parts:
!
!             0. The main calculations in QMC
!               (0) H_0 propagation                        --> TimsH0Prp, TimeH0Prp
!               (1) H_U propagation                        --> TimsbUPrp, TimebUPrp
!               (2) HubbU updating                         --> TimsUptbU, TimeUptbU
!               (3) Numerical stablization                 --> TimsNmStb, TimeNmStb  
!               (4) Population control                     --> TimsPopCt, TimePopCt
!               (5) Growth estimator                       --> TimsConst, TimsConst
!               (7) Static  measurements                   --> TimsStaMe, TimeStaMe
!               (8) Dynamic measurements                   --> TimsDynMe, TimeDynMe
!               (9) Data process in every BIN              --> TimsDatPr, TimeDatPr
!  
!             1. All the matrix operations in QMC
!               (0) UDV decomposition                      --> TimsUDVOt, TimeUDVOt
!               (1) Matrix multiplication                  --> TimsMtprd, TimeMtprd
!               (2) Matrix Inverse                         --> TimsMtInv, TimeMtInv
!               (3) Matrix determinant                     --> TimsMtDet, TimeMtDet
!  
!             2. Time consumed for Beta --> 0 sweep and measure
!               (0) Beta --> 0 sweep and measure           --> TimsM2One, TimeM2One
! 
!             3. Summation of time
!               (0) All the matrix operations              --> TimeMtSum
!               (1) Main calculation time                  --> TimeTotal
!               (2) Time consumed in single BIN            --> TimeSgBIN
! 
! KEYWORDS: Output calculation time.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Output all the calculation time consumed in CPMC simulation in a single BIN.
!
!     Input: (none);  Outpt: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		use RealPrecsn
		use CoreParamt
		use TimeRecord
		use QMCTimeRec
		use Observable
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer NB           ! Loop integer for NmBin calculations
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of Sweep and Measure __________________________________________
!______________________________________________________________________________________________________________
      if(amyid == amstr) then
         !!!!!!!!!! Process the time consumed for every single part
         TimeStaMe = TimeStaMe + TimeB0Mea
         TimeNmStb = TimeNmStb + TimeB0Stb
         TimeTotal = TimeH0Prp + TimebUPrp + TimeUptbU + TimeNmStb + TimePopCt + TimeConst + TimeStaMe + TimeDyMea &
                             & + TimeDatPr
         TimeMtSum = TimeMtFFT + TimeUDVOt + TimeMtprd + TimeMtInv + TimeMtDet + TimeEqSet
         !!!!!!!!!! Open the output file
         open( 281, file = "Output/00_CPMC_Time_Count.txt", access = "append")
         write(281, "()")
         write(281, "()")
         !!!!!!!!!! Print the head of the file for IfFixnT == T case
         if(NB == -1) then
            write(281, "()")
            write(281, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
            write(281, "('&&&&&&&&&&&&&& CPMC Simulation for Data Statistics &&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
            write(281, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
            write(281, "()")
            write(281, "('##########################################################################')")
            write(281, "('####################### CPMC preparing ###################################')") 
            write(281, "('##########################################################################')")
         else if(NB == 0) then
            if(IfMuTqmcNt .or. IfFixnT) then
               write(281, "()")
               write(281, "('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')")
               if(IfMuTqmcNt) then
                  write(281, "('@@@@@@@@@@@ Adjust ChemP_BT to reach the QMC density @@@@@@@@@@@@@@@@')")
               else if(IfFixnT) then
                  write(281, "('@@@@@@@@@@@ Adjust ChemP to reach fixed electron density @@@@@@@@@@@@')")
               end if
               write(281, "('@@@@@@@@@@@@@@@@@@@@@@@@@@ Iteration ', I3.3, ' @@@@@@@@@@@@@@@@@@@@@@@@@@@@')") &
                  & Fix_Iterate
               write(281, "('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')")
               write(281, "()")
               write(281, "('__________________________________________________________________________')")
               write(281, "('__________________________ CPMCWarmUp ____________________________________')") 
               write(281, "('__________________________________________________________________________')")
            else
               write(281, "()")
               write(281, "('==========================================================================')")
               write(281, "('========================== CPMC Warm up ==================================')") 
               write(281, "('==========================================================================')")
            end if
         else if(NB > 0) then
            if(IfMuTqmcNt) then
               write(281, "()")
               write(281, "('__________________________________________________________________________')")
               write(281, "('__________________________ CPMCOccMuT ____________________________________')") 
               write(281, "('__________________________________________________________________________')")
            else if(IfFixnT) then
               write(281, "()")
               write(281, "('__________________________________________________________________________')")
               write(281, "('__________________________ CPMCMeaFix ____________________________________')") 
               write(281, "('__________________________________________________________________________')")
            else
               write(281, "()")
               write(281, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
               write(281, "('++++++++++++++++++++++++++ BIN ', I4.4,' ++++++++++++++++++++++++++++++++++++++')") NB
               write(281, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
               write(281, "(A12, I12)") "NObsStat = ", NObsStat
               write(281, "(A12, I12)") "MeaM2One = ", MeaM2One
               write(281, "(A12, I12)") "NObsDynm = ", NObsDynm
               write(281, "()")
            end if
         end if
			!!!!!!!!!! Output the number of calling (calculation)
         write(281, "(A12, I12)") "TimsH0Prp = ", TimsH0Prp
         write(281, "(A12, I12)") "TimsbUPrp = ", TimsbUPrp
			write(281, "(A12, I12)") "TimsUptbU = ", TimsUptbU
			write(281, "(A12, I12)") "TimsNmStb = ", TimsNmStb
         write(281, "(A12, I12)") "TimsPopCt = ", TimsPopCt
         write(281, "(A12, I12)") "TimsConst = ", TimsConst
         write(281, "(A12, I12)") "TimsStaMe = ", TimsStaMe
         write(281, "(A12, I12)") "TimsDatPr = ", TimsDatPr
         write(281, "()")
         write(281, "(A12, I12)") "TimsM2One = ", TimsM2One
         write(281, "(A12, I12)") "TimsB0Pgt = ", TimsB0Pgt
         write(281, "(A12, I12)") "TimsB0Stb = ", TimsB0Stb
         write(281, "(A12, I12)") "TimsB0Mea = ", TimsB0Mea
         write(281, "(A12, I12)") "TimsDyMea = ", TimsDyMea
         write(281, "()")
         write(281, "(A12, I12)") "TimsMtFFT = ", TimsMtFFT   
			write(281, "(A12, I12)") "TimsUDVOt = ", TimsUDVOt
			write(281, "(A12, I12)") "TimsMtInv = ", TimsMtInv
         write(281, "(A12, I12)") "TimsMtprd = ", TimsMtprd
         write(281, "(A12, I12)") "TimsMtDet = ", TimsMtDet
         write(281, "(A12, I12)") "TimsEqSet = ", TimsEqSet
         write(281, "()")
         write(281, "(A12, I12)") "TimsGFSta = ", TimsGFSta
         write(281, "(A12, I12)") "TimsGFDyn = ", TimsGFDyn

         write(281, "()")
         
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeH0Prp = ", TimeH0Prp, "s = ", &
            & Sec2Min(TimeH0Prp), "m = ", Sec2Hour(TimeH0Prp), "h = ", 100*TimeH0Prp/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimebUPrp = ", TimebUPrp, "s = ", &
            & Sec2Min(TimebUPrp), "m = ", Sec2Hour(TimebUPrp), "h = ", 100*TimebUPrp/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeUptbU = ", TimeUptbU, "s = ", &
            & Sec2Min(TimeUptbU), "m = ", Sec2Hour(TimeUptbU), "h = ", 100*TimeUptbU/TimeSgBIN, "%"
			write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeNmStb = ", TimeNmStb, "s = ", &
            & Sec2Min(TimeNmStb), "m = ", Sec2Hour(TimeNmStb), "h = ", 100*TimeNmStb/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimePopCt = ", TimePopCt, "s = ", &
            & Sec2Min(TimePopCt), "m = ", Sec2Hour(TimePopCt), "h = ", 100*TimePopCt/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeConst = ", TimeConst, "s = ", &
            & Sec2Min(TimeConst), "m = ", Sec2Hour(TimeConst), "h = ", 100*TimeConst/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeStaMe = ", TimeStaMe, "s = ", &
            & Sec2Min(TimeStaMe), "m = ", Sec2Hour(TimeStaMe), "h = ", 100*TimeStaMe/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeDyMea = ", TimeDyMea, "s = ", &
            & Sec2Min(TimeDyMea), "m = ", Sec2Hour(TimeDyMea), "h = ", 100*TimeDyMea/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeDatPr = ", TimeDatPr, "s = ", &
            & Sec2Min(TimeDatPr), "m = ", Sec2Hour(TimeDatPr), "h = ", 100*TimeDatPr/TimeSgBIN, "%"
         write(281, "()")
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeM2One = ", TimeM2One, "s = ", &
            & Sec2Min(TimeM2One), "m = ", Sec2Hour(TimeM2One), "h = ", 100*TimeM2One/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeB0Pgt = ", TimeB0Pgt, "s = ", &
            & Sec2Min(TimeB0Pgt), "m = ", Sec2Hour(TimeB0Pgt), "h = ", 100*TimeB0Pgt/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeB0Stb = ", TimeB0Stb, "s = ", &
            & Sec2Min(TimeB0Stb), "m = ", Sec2Hour(TimeB0Stb), "h = ", 100*TimeB0Stb/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeB0Mea = ", TimeB0Mea, "s = ", &
            & Sec2Min(TimeB0Mea), "m = ", Sec2Hour(TimeB0Mea), "h = ", 100*TimeB0Mea/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeDyMea = ", TimeDyMea, "s = ", &
            & Sec2Min(TimeDyMea), "m = ", Sec2Hour(TimeDyMea), "h = ", 100*TimeDyMea/TimeSgBIN, "%"
         write(281, "()")
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeMtFFT = ", TimeMtFFT, "s = ", &
            & Sec2Min(TimeMtFFT), "m = ", Sec2Hour(TimeMtFFT), "h = ", 100*TimeMtFFT/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeUDVOt = ", TimeUDVOt, "s = ", &
            & Sec2Min(TimeUDVOt), "m = ", Sec2Hour(TimeUDVOt), "h = ", 100*TimeUDVOt/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeMtInv = ", TimeMtInv, "s = ", &
            & Sec2Min(TimeMtInv), "m = ", Sec2Hour(TimeMtInv), "h = ", 100*TimeMtInv/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeMtprd = ", TimeMtprd, "s = ", &
            & Sec2Min(TimeMtprd), "m = ", Sec2Hour(TimeMtprd), "h = ", 100*TimeMtprd/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeMtDet = ", TimeMtDet, "s = ", &
            & Sec2Min(TimeMtDet), "m = ", Sec2Hour(TimeMtDet), "h = ", 100*TimeMtDet/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeEqSet = ", TimeEqSet, "s = ", &
            & Sec2Min(TimeEqSet), "m = ", Sec2Hour(TimeEqSet), "h = ", 100*TimeEqSet/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeMtSum = ", TimeMtSum, "s = ", &
            & Sec2Min(TimeMtSum), "m = ", Sec2Hour(TimeMtSum), "h = ", 100*TimeMtSum/TimeSgBIN, "%"
         write(281, "()")
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeGFSta = ", TimeGFSta, "s = ", &
            & Sec2Min(TimeGFSta), "m = ", Sec2Hour(TimeGFSta), "h = ", 100*TimeGFSta/TimeSgBIN, "%"
         write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeGFDyn = ", TimeGFDyn, "s = ", &
            & Sec2Min(TimeGFDyn), "m = ", Sec2Hour(TimeGFDyn), "h = ", 100*TimeGFDyn/TimeSgBIN, "%"

         write(281, "()")
         
			write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeTotal = ", TimeTotal, "s = ", &
            & Sec2Min(TimeTotal), "m = ", Sec2Hour(TimeTotal), "h = ", 100*TimeTotal/TimeSgBIN, "%"
			write(281, "(A12, F16.5, A4, F13.5, A4, F10.5, A4, F10.5, A1)") "TimeSgBIN = ", TimeSgBIN, "s = ", &
            & Sec2Min(TimeSgBIN), "m = ", Sec2Hour(TimeSgBIN), "h = ", 100*TimeSgBIN/TimeSgBIN, "%"
			
         if(IfFixnT) then
            write(281, "('__________________________________________________________________________')")
            write(281, "('__________________________________________________________________________')")
            write(281, "('__________________________________________________________________________')")
         else
            if(NB == -1) then
               write(281, "('##########################################################################')")
               write(281, "('##########################################################################')")
               write(281, "('##########################################################################')")
            else if(NB == 0) then 
               write(281, "('==========================================================================')")
               write(281, "('==========================================================================')")
               write(281, "('==========================================================================')")
            else
               write(281, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
               write(281, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
               write(281, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
            end if
         end if
         close(281)
      end if
		
	end subroutine SaveCalTim
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
