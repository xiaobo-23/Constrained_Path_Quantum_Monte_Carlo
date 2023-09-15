!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A few subroutines used for calculating the multiplication of AMat matrix and the exponential Kinetic matrix. 
!               Including 4 different cases as: 
!                     (1) AMat* exp(-/+\Delta\tau*H_0);       (2)  exp(-/+\Delta\tau*H_0)*AMat; 
!                     (3) AMat*[exp(-/+\Delta\tau*H_0)]^{-1}; (4) [exp(-/+\Delta\tau*H_0)]^{-1}*AMat;
!          Includes four different subroutines.
! COMMENT: Exponential kinetic matrix multiplying matrix AMat (at left and right side).
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!     Calculate Exp(-/+\Delta\tau*H_0)*AMat  and  Exp(-/+\Delta\tau*H_0/2)*AMat  and
!               AMat*Exp(-/+\Delta\tau*H_0)  and  AMat*Exp(-/+\Delta\tau*H_0/2).
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   These four subroutines are for exp(-/+\Delta\tau*H_0) and exp(-/+\Delta\tau*H_0/2)
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   RghtMultExpH0    --> Calculate  exp(-\Delta\tau*H_0)      *AMat;
!   RghtMultExpH0Inv --> Calculate [exp(-\Delta\tau*H_0)]^{-1}*AMat;
!   LeftMultExpH0    --> Calculate AMat* exp(-\Delta\tau*H_0)      ;
!   LeftMultExpH0Inv --> Calculate AMat*[exp(-\Delta\tau*H_0)]^{-1}.
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   These four subroutines are for exp(-/+\Delta\tau*H_T) and exp(-/+\Delta\tau*H_T/2)
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   RghtMultExpHT    --> Calculate  exp(-\Delta\tau*H_T)      *AMat;
!   RghtMultExpHTInv --> Calculate [exp(-\Delta\tau*H_T)]^{-1}*AMat;
!   LeftMultExpHT    --> Calculate AMat* exp(-\Delta\tau*H_T)      ;
!   LeftMultExpHTInv --> Calculate AMat*[exp(-\Delta\tau*H_T)]^{-1}.
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   These four subroutines Calculate the matrix products by matrix multiplication.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ExpdtKMultMatFull_R --> Calculate ExpK*AMat simply by matrix multiplication, real version;
!   MatMultExpdtKFull_R --> Calculate AMat*ExpK simply by matrix multiplication, real version;
!   ExpdtKMultMatFull_C --> Calculate ExpK*AMat simply by matrix multiplication, complex version;
!   MatMultExpdtKFull_C --> Calculate AMat*ExpK simply by matrix multiplication, complex version;
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   These four subroutines Calculate the matrix products by FFT method.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ExpdtKMultFFTEgVl_R --> Calculate ExpK*AMat simply by FFT method, real version;
!   FFTEgVlMultExpdtK_R --> Calculate AMat*ExpK simply by FFT method, real version;
!   ExpdtKMultFFTEgVl_C --> Calculate ExpK*AMat simply by FFT method, complex version;
!   FFTEgVlMultExpdtK_C --> Calculate AMat*ExpK simply by FFT method, complex version.
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   These four subroutines transform the wavefunction from r to k space.
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   RghtWvfc_r2k_R --> Transform wavefunction at right side from r to k space, real version;
!   LeftWvfc_r2k_R --> Transform wavefunction at left  side from r to k space, real version;
!   RghtWvfc_r2k_C --> Transform wavefunction at right side from r to k space, complex version;
!   LeftWvfc_r2k_C --> Transform wavefunction at Left  side from r to k space, complex version;
!  
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!########################################################################################################################
!########################################################################################################################
!######################### Exp(-/+\Delta\tau*H_0(/2))*AMat and AMat*Exp(-/+\Delta\tau*H_0(/2)) ##########################
!######################### Exp(-/+\Delta\tau*H_0(/2))*AMat and AMat*Exp(-/+\Delta\tau*H_0(/2)) ##########################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine RghtMultExpH0(IfHalfExpK, AMat, ND2) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RghtMultExpH0(IfHalfExpK, AMat, ND2)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculations as  Exp(-\Delta\tau*H_0)*AMat  or  
!                         Exp(-\Delta\tau*H_0/2)*AMat.
! KEYWORDS: Exp(-\Delta\tau*H_0)*AMat  or  Exp(-\Delta\tau*H_0/2)*AMat.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate Exp(-\Delta\tau*H_0)*AMat  or  Exp(-\Delta\tau*H_0/2)*AMat.
!
!     Input: IfHalfExpK --> Logical quantity to indicate using H_0 or H_0/2;
!            AMat       --> Input matrix;
!            ND2        --> Dimension of AMat(NumNS, ND2, NmSpn) matrix.
! 
!     Outpt: AMat       --> Finally = exp(-\Delta\tau*H_0)*AMat  or  exp(-\Delta\tau*H_0/2)*AMat
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfHalfExpK                 ! == T, use H_0/2; == F, use H_0
      integer ND2                        ! Dimension of AMat
      real(rp) AMat(NumNS, NumNS, NmSpn) ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                   ! Starting time point
      integer(8) time2                   ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsH0Prp = TimsH0Prp + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate exp(-\Delta\tau*H_0)*AMat  or  exp(-\Delta\tau*H_0/2)*AMat _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Two methods to calculate Exp(-\Delta\tau*H_0(/2)))*AMat ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For FFTEXPDTH0 == T case, FFT for matrix product _________________
!________________________________________________________________________________________
      if(FFTEXPDTH0) then
         if(IfHalfExpK) then
            call ExpdtKMultFFTEgVl_R(ExpdtEofH0(1, 1, 3), AMat(1, 1, 1), NumNS, ND2)
         else
            call ExpdtKMultFFTEgVl_R(ExpdtEofH0(1, 1, 1), AMat(1, 1, 1), NumNS, ND2)
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) For FFTEXPDTH0 == F case, direct matrix multiplication ___________
!________________________________________________________________________________________
      else
         if(IfHalfExpK) then
            call ExpdtKMultMatFull_R(ExpdtH0Mat(1, 1, 1, 3), AMat(1, 1, 1), NumNS, ND2)
         else
            call ExpdtKMultMatFull_R(ExpdtH0Mat(1, 1, 1, 1), AMat(1, 1, 1), NumNS, ND2)
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeH0Prp = TimeH0Prp + TimeIntrvl(time1, time2)
		
   end subroutine RghtMultExpH0
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine RghtMultExpH0Inv(IfHalfExpK, AMat, ND2) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RghtMultExpH0Inv(IfHalfExpK, AMat, ND2)  
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculations as  Exp(+\Delta\tau*H_0)*AMat  or  
!                     Exp(+\Delta\tau*H_0/2)*AMat.
! KEYWORDS: Exp(+\Delta\tau*H_0)*AMat  or  Exp(+\Delta\tau*H_0/2)*AMat.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate Exp(+\Delta\tau*H_0)*AMat  or  Exp(+\Delta\tau*H_0/2)*AMat.
!
!     Input: IfHalfExpK --> Logical quantity to indicate using H_0 or H_0/2;
!            AMat       --> Input matrix;
!            ND2        --> Dimension of AMat(NumNS, ND2, NmSpn) matrix.
! 
!     Outpt: AMat       --> Finally = exp(+\Delta\tau*H_0)*AMat  or  exp(+\Delta\tau*H_0/2)*AMat.
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfHalfExpK                  ! == T, use H_0/2; == F, use H_0
      integer ND2                         ! Dimension of AMat matrix
      real(rp) AMat(NumNS, NumNS, NmSpn)  ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsH0Prp = TimsH0Prp + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate exp(+\Delta\tau*H_0)*AMat  or  exp(+\Delta\tau*H_0/2)*AMat _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Two methods to calculate Exp(+\Delta\tau*H_0(/2)))*AMat ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For FFTEXPDTH0 == T case, FFT for matrix product _________________
!________________________________________________________________________________________
      if(FFTEXPDTH0) then
         if(IfHalfExpK) then
            call ExpdtKMultFFTEgVl_R(ExpdtEofH0(1, 1, 4), AMat(1, 1, 1), NumNS, ND2)
         else
            call ExpdtKMultFFTEgVl_R(ExpdtEofH0(1, 1, 2), AMat(1, 1, 1), NumNS, ND2)
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) For FFTEXPDTH0 == F case, direct matrix multiplication ___________
!________________________________________________________________________________________
      else
         if(IfHalfExpK) then
            call ExpdtKMultMatFull_R(ExpdtH0Mat(1, 1, 1, 4), AMat(1, 1, 1), NumNS, ND2)
         else
            call ExpdtKMultMatFull_R(ExpdtH0Mat(1, 1, 1, 2), AMat(1, 1, 1), NumNS, ND2)
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeH0Prp = TimeH0Prp + TimeIntrvl(time1, time2)
		
   end subroutine RghtMultExpH0Inv
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine LeftMultExpH0(IfHalfExpK, ND1, AMat) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  LeftMultExpH0(IfHalfExpK, ND1, AMat)  
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculations as  AMat*Exp(-\Delta\tau*H_0)  or  
!                    AMat*Exp(-\Delta\tau*H_0/2).
! KEYWORDS: AMat*Exp(-\Delta\tau*H_0)  or  AMat*Exp(-\Delta\tau*H_0/2).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate AMat*Exp(-\Delta\tau*H_0)  or  AMat*Exp(-\Delta\tau*H_0/2).
!
!     Input: IfHalfExpK --> Logical quantity to indicate using H_0 or H_0/2;
!            ND1        --> Dimension of AMat matrix;
!            AMat       --> Input matrix;
! 
!     Outpt: AMat       --> Finally = AMat*Exp(-\Delta\tau*H_0)  or  AMat*Exp(-\Delta\tau*H_0/2)
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfHalfExpK                  ! == T, use H_0/2; == F, use H_0
      integer ND1
      real(rp) AMat(NumNS, NumNS, NmSpn)  ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsH0Prp = TimsH0Prp + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate AMat*Exp(-\Delta\tau*H_0)  or  AMat*Exp(-\Delta\tau*H_0/2) _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Two methods to calculate AMat*Exp(-\Delta\tau*H_0(/2))) ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For FFTEXPDTH0 == T case, FFT for matrix product _________________
!________________________________________________________________________________________
      if(FFTEXPDTH0) then
         if(IfHalfExpK) then
            call FFTEgVlMultExpdtK_R(ND1, NumNS, AMat(1, 1, 1), ExpdtEofH0(1, 1, 3))
         else
            call FFTEgVlMultExpdtK_R(ND1, NumNS, AMat(1, 1, 1), ExpdtEofH0(1, 1, 1))
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) For FFTEXPDTH0 == F case, direct matrix multiplication ___________
!________________________________________________________________________________________
      else  
         if(IfHalfExpK) then
            call MatMultExpdtKFull_R(ND1, NumNS, AMat(1, 1, 1), ExpdtH0Mat(1, 1, 1, 3))
         else
            call MatMultExpdtKFull_R(ND1, NumNS, AMat(1, 1, 1), ExpdtH0Mat(1, 1, 1, 1))
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeH0Prp = TimeH0Prp + TimeIntrvl(time1, time2)
		
   end subroutine LeftMultExpH0
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine LeftMultExpH0Inv(IfHalfExpK, ND1, AMat) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  LeftMultExpH0Inv(IfHalfExpK, ND1, AMat)  
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculations as  AMat*Exp(+\Delta\tau*H_0)  or  
!                  AMat*Exp(+\Delta\tau*H_0/2).
! KEYWORDS: AMat*Exp(+\Delta\tau*H_0)  or  AMat*Exp(+\Delta\tau*H_0/2).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate AMat*Exp(+\Delta\tau*H_0)  or  AMat*Exp(+\Delta\tau*H_0/2).
!
!     Input: IfHalfExpK --> Logical quantity to indicate using H_0 or H_0/2;
!            ND1        --> Dimension of AMat matrix;
!            AMat       --> Input matrix;
! 
!     Outpt: AMat      --> Finally = AMat*Exp(+\Delta\tau*H_0)  or  AMat*Exp(+\Delta\tau*H_0/2)
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfHalfExpK                  ! == T, use H_0/2; == F, use H_0
      integer ND1
      real(rp) AMat(NumNS, NumNS, NmSpn)  ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsH0Prp = TimsH0Prp + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate AMat*Exp(+\Delta\tau*H_0)  or  AMat*Exp(+\Delta\tau*H_0/2) _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Two methods to calculate AMat*Exp(+\Delta\tau*H_0(/2))) ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For FFTEXPDTH0 == T case, FFT for matrix product _________________
!________________________________________________________________________________________
      if(FFTEXPDTH0) then
         if(IfHalfExpK) then
            call FFTEgVlMultExpdtK_R(ND1, NumNS, AMat(1, 1, 1), ExpdtEofH0(1, 1, 4))
         else
            call FFTEgVlMultExpdtK_R(ND1, NumNS, AMat(1, 1, 1), ExpdtEofH0(1, 1, 2))
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) For FFTEXPDTH0 == F case, direct matrix multiplication ___________
!________________________________________________________________________________________
      else
         if(IfHalfExpK) then
            call MatMultExpdtKFull_R(ND1, NumNS, AMat(1, 1, 1), ExpdtH0Mat(1, 1, 1, 4))
         else
            call MatMultExpdtKFull_R(ND1, NumNS, AMat(1, 1, 1), ExpdtH0Mat(1, 1, 1, 2))
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeH0Prp = TimeH0Prp + TimeIntrvl(time1, time2)
		
   end subroutine LeftMultExpH0Inv
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!######################### Exp(-/+\Delta\tau*H_T(/2))*AMat and AMat*Exp(-/+\Delta\tau*H_T(/2)) ##########################
!######################### Exp(-/+\Delta\tau*H_T(/2))*AMat and AMat*Exp(-/+\Delta\tau*H_T(/2)) ##########################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine RghtMultExpHT(IfHalfExpK, AMat, ND2) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RghtMultExpHT(IfHalfExpK, AMat, ND2)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculations as  Exp(-\Delta\tau*H_T)*AMat  or  
!                         Exp(-\Delta\tau*H_T/2)*AMat.
! KEYWORDS: Exp(-\Delta\tau*H_T)*AMat  or  Exp(-\Delta\tau*H_T/2)*AMat.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate Exp(-\Delta\tau*H_T)*AMat  or  Exp(-\Delta\tau*H_T/2)*AMat.
!
!     Input: IfHalfExpK --> Logical quantity to indicate using H_T or H_T/2;
!            AMat       --> Input matrix;
!            ND2        --> Dimension of AMat(NumNS, ND2, NmSpn) matrix.
! 
!     Outpt: AMat       --> Finally = exp(-\Delta\tau*H_T)*AMat  or  exp(-\Delta\tau*H_T/2)*AMat
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfHalfExpK                 ! == T, use H_T/2; == F, use H_T
      integer ND2                        ! Dimension of AMat
      real(rp) AMat(NumNS, NumNS, NmSpn) ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                   ! Starting time point
      integer(8) time2                   ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_T propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsH0Prp = TimsH0Prp + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate exp(-\Delta\tau*H_T)*AMat  or  exp(-\Delta\tau*H_T/2)*AMat _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Two methods to calculate Exp(-\Delta\tau*H_T(/2)))*AMat ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For FFTEXPDTHT == T case, FFT for matrix product _________________
!________________________________________________________________________________________
      if(FFTEXPDTHT) then
         if(IfHalfExpK) then
            call ExpdtKMultFFTEgVl_R(ExpdtOfHTe(1, 1, 3), AMat(1, 1, 1), NumNS, ND2)
         else
            call ExpdtKMultFFTEgVl_R(ExpdtOfHTe(1, 1, 1), AMat(1, 1, 1), NumNS, ND2)
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) For FFTEXPDTHT == F case, direct matrix multiplication ___________
!________________________________________________________________________________________
      else
         if(IfHalfExpK) then
            call ExpdtKMultMatFull_R(ExpdtTryHT(1, 1, 1, 3), AMat(1, 1, 1), NumNS, ND2)
         else
            call ExpdtKMultMatFull_R(ExpdtTryHT(1, 1, 1, 1), AMat(1, 1, 1), NumNS, ND2)
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_T propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeH0Prp = TimeH0Prp + TimeIntrvl(time1, time2)
		
   end subroutine RghtMultExpHT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine RghtMultExpHTInv(IfHalfExpK, AMat, ND2) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RghtMultExpHTInv(IfHalfExpK, AMat, ND2)  
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculations as  Exp(+\Delta\tau*H_T)*AMat  or  
!                     Exp(+\Delta\tau*H_T/2)*AMat.
! KEYWORDS: Exp(+\Delta\tau*H_T)*AMat  or  Exp(+\Delta\tau*H_T/2)*AMat.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate Exp(+\Delta\tau*H_T)*AMat  or  Exp(+\Delta\tau*H_T/2)*AMat.
!
!     Input: IfHalfExpK --> Logical quantity to indicate using H_T or H_T/2;
!            AMat       --> Input matrix;
!            ND2        --> Dimension of AMat(NumNS, ND2, NmSpn) matrix.
! 
!     Outpt: AMat       --> Finally = exp(+\Delta\tau*H_T)*AMat  or  exp(+\Delta\tau*H_T/2)*AMat.
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfHalfExpK                  ! == T, use H_T/2; == F, use H_T
      integer ND2                         ! Dimension of AMat matrix
      real(rp) AMat(NumNS, NumNS, NmSpn)  ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_T propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsH0Prp = TimsH0Prp + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate exp(+\Delta\tau*H_T)*AMat  or  exp(+\Delta\tau*H_T/2)*AMat _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Two methods to calculate Exp(+\Delta\tau*H_T(/2)))*AMat ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For FFTEXPDTHT == T case, FFT for matrix product _________________
!________________________________________________________________________________________
      if(FFTEXPDTHT) then
         if(IfHalfExpK) then
            call ExpdtKMultFFTEgVl_R(ExpdtOfHTe(1, 1, 4), AMat(1, 1, 1), NumNS, ND2)
         else
            call ExpdtKMultFFTEgVl_R(ExpdtOfHTe(1, 1, 2), AMat(1, 1, 1), NumNS, ND2)
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) For FFTEXPDTHT == F case, direct matrix multiplication ___________
!________________________________________________________________________________________
      else
         if(IfHalfExpK) then
            call ExpdtKMultMatFull_R(ExpdtTryHT(1, 1, 1, 4), AMat(1, 1, 1), NumNS, ND2)
         else
            call ExpdtKMultMatFull_R(ExpdtTryHT(1, 1, 1, 2), AMat(1, 1, 1), NumNS, ND2)
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_T propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeH0Prp = TimeH0Prp + TimeIntrvl(time1, time2)
		
   end subroutine RghtMultExpHTInv
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine LeftMultExpHT(IfHalfExpK, ND1, AMat) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  LeftMultExpHT(IfHalfExpK, ND1, AMat)  
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculations as  AMat*Exp(-\Delta\tau*H_T)  or  
!                    AMat*Exp(-\Delta\tau*H_T/2).
! KEYWORDS: AMat*Exp(-\Delta\tau*H_T)  or  AMat*Exp(-\Delta\tau*H_T/2).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate AMat*Exp(-\Delta\tau*H_T)  or  AMat*Exp(-\Delta\tau*H_T/2).
!
!     Input: IfHalfExpK --> Logical quantity to indicate using H_T or H_T/2;
!            ND1        --> Dimension of AMat matrix;
!            AMat       --> Input matrix;
! 
!     Outpt: AMat       --> Finally = AMat*Exp(-\Delta\tau*H_T)  or  AMat*Exp(-\Delta\tau*H_T/2)
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfHalfExpK                  ! == T, use H_T/2; == F, use H_T
      integer ND1
      real(rp) AMat(NumNS, NumNS, NmSpn)  ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_T propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsH0Prp = TimsH0Prp + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate AMat*Exp(-\Delta\tau*H_T)  or  AMat*Exp(-\Delta\tau*H_T/2) _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Two methods to calculate AMat*Exp(-\Delta\tau*H_T(/2))) ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For FFTEXPDTHT == T case, FFT for matrix product _________________
!________________________________________________________________________________________
      if(FFTEXPDTHT) then
         if(IfHalfExpK) then
            call FFTEgVlMultExpdtK_R(ND1, NumNS, AMat(1, 1, 1), ExpdtOfHTe(1, 1, 3))
         else
            call FFTEgVlMultExpdtK_R(ND1, NumNS, AMat(1, 1, 1), ExpdtOfHTe(1, 1, 1))
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) For FFTEXPDTHT == F case, direct matrix multiplication ___________
!________________________________________________________________________________________
      else  
         if(IfHalfExpK) then
            call MatMultExpdtKFull_R(ND1, NumNS, AMat(1, 1, 1), ExpdtTryHT(1, 1, 1, 3))
         else
            call MatMultExpdtKFull_R(ND1, NumNS, AMat(1, 1, 1), ExpdtTryHT(1, 1, 1, 1))
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_T propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeH0Prp = TimeH0Prp + TimeIntrvl(time1, time2)
		
   end subroutine LeftMultExpHT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine LeftMultExpHTInv(IfHalfExpK, ND1, AMat) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  LeftMultExpHTInv(IfHalfExpK, ND1, AMat)  
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculations as  AMat*Exp(+\Delta\tau*H_0)  or  
!                  AMat*Exp(+\Delta\tau*H_0/2).
! KEYWORDS: AMat*Exp(+\Delta\tau*H_0)  or  AMat*Exp(+\Delta\tau*H_0/2).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate AMat*Exp(+\Delta\tau*H_0)  or  AMat*Exp(+\Delta\tau*H_0/2).
!
!     Input: IfHalfExpK --> Logical quantity to indicate using H_0 or H_0/2;
!            ND1        --> Dimension of AMat matrix;
!            AMat       --> Input matrix;
! 
!     Outpt: AMat      --> Finally = AMat*Exp(+\Delta\tau*H_0)  or  AMat*Exp(+\Delta\tau*H_0/2)
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfHalfExpK                  ! == T, use H_0/2; == F, use H_0
      integer ND1
      real(rp) AMat(NumNS, NumNS, NmSpn)  ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsH0Prp = TimsH0Prp + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate AMat*Exp(+\Delta\tau*H_0)  or  AMat*Exp(+\Delta\tau*H_0/2) _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Two methods to calculate AMat*Exp(+\Delta\tau*H_0(/2))) ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For FFTEXPDTHT == T case, FFT for matrix product _________________
!________________________________________________________________________________________
      if(FFTEXPDTHT) then
         if(IfHalfExpK) then
            call FFTEgVlMultExpdtK_R(ND1, NumNS, AMat(1, 1, 1), ExpdtOfHTe(1, 1, 4))
         else
            call FFTEgVlMultExpdtK_R(ND1, NumNS, AMat(1, 1, 1), ExpdtOfHTe(1, 1, 2))
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) For FFTEXPDTHT == F case, direct matrix multiplication ___________
!________________________________________________________________________________________
      else
         if(IfHalfExpK) then
            call MatMultExpdtKFull_R(ND1, NumNS, AMat(1, 1, 1), ExpdtTryHT(1, 1, 1, 4))
         else
            call MatMultExpdtKFull_R(ND1, NumNS, AMat(1, 1, 1), ExpdtTryHT(1, 1, 1, 2))
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in H_0 propagation process ________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeH0Prp = TimeH0Prp + TimeIntrvl(time1, time2)
		
   end subroutine LeftMultExpHTInv
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   

!########################################################################################################################
!########################################################################################################################
!######################################## Direct Matrix multiplication in real space ####################################
!######################################## Direct Matrix multiplication in real space ####################################
!######################################## Direct Matrix multiplication in real space ####################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ExpdtKMultMatFull_R(ExpMat, AMat, ND2, ND2Run)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ExpdtKMultMatFull_R(ExpMat, AMat, ND2, ND2Run)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculation as ExpMat*AMat, with spin coupled and decoupled, 
!                     by direct matrix product.
! KEYWORDS: ExpMat*AMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate ExpMat*AMat with spin coupled or decoupled.
!
!     Input: ExpMat --> Input Exp matrix;
!            AMat   --> Input matrix;
!            ND2    --> Dimension of AMat matrix;
!            ND2Run --> Dimension that will be run;
!
!     Outpt: AMat   --> Finally = ExpMat*AMat.
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND2, ND2Run
      real(rp) ExpMat(NumNS, NumNS, NmSpn)    ! The input Exp(-/+dt*H0) or Exp(-/+dt*H0/2) Matrix
      real(rp)   AMat(NumNS, ND2  , NmSpn)    ! Input/Output real matrix     
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                        ! Starting time point
		integer(8) time2                        ! Ending   time point
      integer SpnInd
      real(rp), allocatable :: TmpMat(:, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate ExpMat*AMat with spin coupled and decoupled for H_0 ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Allocate the temporary matrix and initialization __________________________
!**************************************************************************************************
      allocate(TmpMat(NumNS, ND2Run))
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == T case --> NmSpn == 1 and Only spin-up channel __________
!______________________ For DCPDNOSIGN == F case --> NmSpn == 2 and both channels _________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Perform matrix multiplication ____________________________________
!________________________________________________________________________________________ 
         TmpMat = 0.0_rp
         call DGEMM("N", "N", NumNS, ND2Run, NumNS, 1.0_rp, ExpMat(1, 1, SpnInd), NumNS, AMat(1, 1, SpnInd), NumNS, &
            & 0.0_rp, TmpMat(1, 1), NumNS)
!________________________________________________________________________________________ 	  
!_________________ (1) Final result of ExpMat*AMat ______________________________________
!________________________________________________________________________________________
         AMat(1:NumNS, 1:ND2Run, SpnInd) = TmpMat(1:NumNS, 1:ND2Run)
      enddo
!**************************************************************************************************	  
!___________________ 2. Deallocate the temporary matrix ___________________________________________
!**************************************************************************************************
      if(allocated(TmpMat)) deallocate(TmpMat)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine ExpdtKMultMatFull_R
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine MatMultExpdtKFull_R(ND1Run, ND1, AMat, ExpMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  MatMultExpdtKFull_R(ND1Run, ND1, AMat, ExpMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculation as AMat*ExpMat, with spin coupled and decoupled, 
!                   by direct matrix product.
! KEYWORDS: AMat*ExpMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate AMat*ExpMat with spin coupled or decoupled.
!
!     Input: ExpMat --> Input Exp matrix;
!            AMat   --> Input matrix;
!            ND1Run --> Dimension that will be run;
!            ND1    --> Dimension of AMat;
! 
!     Outpt: AMat   --> Finally = AMat*ExpMat.
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1Run, ND1
      real(rp)   AMat(ND1  , NumNS, NmSpn)    ! Input/Output real matrix
      real(rp) ExpMat(NumNS, NumNS, NmSpn)    ! The input Exp(-/+dt*H0) or Exp(-/+dt*H0/2) Matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                    ! Starting time point
      integer(8) time2                    ! Ending   time point
      integer SpnInd
      real(rp), allocatable :: TmpMat(:, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate AMat*ExpMat with spin coupled and decoupled for H_0 ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Allocate the temporary matrix and initialization __________________________
!**************************************************************************************************
      allocate(TmpMat(ND1Run, NumNS))
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == T case --> NmSpn == 1 and Only spin-up channel __________
!______________________ For DCPDNOSIGN == F case --> NmSpn == 2 and both channels _________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Perform matrix multiplication ____________________________________
!________________________________________________________________________________________
         TmpMat = 0.0_rp
         call DGEMM("N", "N", ND1Run, NumNS, NumNS, 1.0_rp, AMat(1, 1, SpnInd), ND1, ExpMat(1, 1, SpnInd), NumNS, &
            & 0.0_rp, TmpMat(1, 1), ND1Run)
!________________________________________________________________________________________ 	  
!_________________ (1) Final result of AMat*ExpMat ______________________________________
!________________________________________________________________________________________
         AMat(1:ND1Run, 1:NumNS, SpnInd) = TmpMat(1:ND1Run, 1:NumNS)
      enddo
!**************************************************************************************************	  
!___________________ 2. Deallocate the temporary matrix ___________________________________________
!**************************************************************************************************
      if(allocated(TmpMat)) deallocate(TmpMat)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine MatMultExpdtKFull_R
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ExpdtKMultMatFull_C(ExpMat, AMat, ND2, ND2Run)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ExpdtKMultMatFull_C(ExpMat, AMat, ND2, ND2Run)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculation as ExpMat*AMat, with spin coupled and decoupled, 
!                     by direct matrix product.
! KEYWORDS: ExpMat*AMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate ExpMat*AMat with spin coupled or decoupled.
!
!     Input: ExpMat --> Input Exp matrix;
!            AMat   --> Input matrix;
!            ND2    --> Dimension of AMat matrix;
!            ND2Run --> Dimension that will be run;
!
!     Outpt: AMat   --> Finally = ExpMat*AMat.
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND2, ND2Run
      complex(rp) ExpMat(NumNS, NumNS, NmSpn)    ! The input Exp(-/+dt*H0) or Exp(-/+dt*H0/2) Matrix
      complex(rp)   AMat(NumNS, ND2  , NmSpn)    ! Input/Output complex matrix     
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                        ! Starting time point
		integer(8) time2                        ! Ending   time point
      integer SpnInd
      complex(rp), allocatable :: TmpMat(:, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate ExpMat*AMat with spin coupled and decoupled for H_0 ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Allocate the temporary matrix and initialization __________________________
!**************************************************************************************************
      allocate(TmpMat(NumNS, ND2Run))
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == T case --> NmSpn == 1 and Only spin-up channel __________
!______________________ For DCPDNOSIGN == F case --> NmSpn == 2 and both channels _________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Perform matrix multiplication ____________________________________
!________________________________________________________________________________________ 
         TmpMat = rp_Zzero
         call ZGEMM("N", "N", NumNS, ND2Run, NumNS, rp_Z_One, ExpMat(1, 1, SpnInd), NumNS, AMat(1, 1, SpnInd), NumNS, &
            & rp_Zzero, TmpMat(1, 1), NumNS)
!________________________________________________________________________________________ 	  
!_________________ (1) Final result of ExpMat*AMat ______________________________________
!________________________________________________________________________________________
            AMat(1:NumNS, 1:ND2Run, SpnInd) = TmpMat(1:NumNS, 1:ND2Run)
      enddo
!**************************************************************************************************	  
!___________________ 2. Deallocate the temporary matrix ___________________________________________
!**************************************************************************************************
      if(allocated(TmpMat)) deallocate(TmpMat)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine ExpdtKMultMatFull_C
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine MatMultExpdtKFull_C(ND1Run, ND1, AMat, ExpMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  MatMultExpdtKFull_C(ND1Run, ND1, AMat, ExpMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculation as AMat*ExpMat, with spin coupled and decoupled, 
!                   by direct matrix product.
! KEYWORDS: AMat*ExpMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate AMat*ExpMat with spin coupled or decoupled.
!
!     Input: ExpMat --> Input Exp matrix;
!            AMat   --> Input matrix;
!            ND1Run --> Dimension that will be run;
!            ND1    --> Dimension of AMat;
! 
!     Outpt: AMat   --> Finally = AMat*ExpMat.
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
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1Run, ND1
      complex(rp)   AMat(ND1  , NumNS, NmSpn)    ! Input/Output complex matrix
      complex(rp) ExpMat(NumNS, NumNS, NmSpn)    ! The input Exp(-/+dt*H0) or Exp(-/+dt*H0/2) Matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                    ! Starting time point
      integer(8) time2                    ! Ending   time point
      integer SpnInd
      complex(rp), allocatable :: TmpMat(:, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate AMat*ExpMat with spin coupled and decoupled for H_0 ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Allocate the temporary matrix and initialization __________________________
!**************************************************************************************************
      allocate(TmpMat(ND1Run, NumNS))
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == T case --> NmSpn == 1 and Only spin-up channel __________
!______________________ For DCPDNOSIGN == F case --> NmSpn == 2 and both channels _________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Perform matrix multiplication ____________________________________
!________________________________________________________________________________________
         TmpMat = rp_Zzero
         call ZGEMM("N", "N", ND1Run, NumNS, NumNS, rp_Z_One, AMat(1, 1, SpnInd), ND1, ExpMat(1, 1, SpnInd), NumNS, &
            & rp_Zzero, TmpMat(1, 1), ND1Run)
!________________________________________________________________________________________ 	  
!_________________ (1) Final result of AMat*ExpMat ______________________________________
!________________________________________________________________________________________
         AMat(1:ND1Run, 1:NumNS, SpnInd) = TmpMat(1:ND1Run, 1:NumNS)
      enddo
!**************************************************************************************************	  
!___________________ 2. Deallocate the temporary matrix ___________________________________________
!**************************************************************************************************
      if(allocated(TmpMat)) deallocate(TmpMat)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine MatMultExpdtKFull_C
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!######################################## Matrix multiplication by FFT method ###########################################
!######################################## Matrix multiplication by FFT method ###########################################
!######################################## Matrix multiplication by FFT method ###########################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ExpdtKMultFFTEgVl_R(ExpdtEk, AMat, ND2, ND2Run)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ExpdtKMultFFTEgVl_R(ExpdtEk, AMat, ND2, ND2Run)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculation as ExpMat*AMat, with spin coupled and decoupled, 
!                     by FFT method.
! KEYWORDS: ExpMat*AMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate ExpMat*AMat with spin coupled or decoupled.
!
!     Input: ExpdtEk --> Exp(-dt*E_k) as the diagonal of the ExpMat matrix;
!            AMat    --> Input matrix;
!            ND2     --> Dimension of AMat;
!            ND2Run  --> Dimension that will be run;
! 
!     Outpt: AMat   --> Finally = ExpMat*AMat.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use RealPrecsn
      use FFTSetting
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND2, ND2Run
      real(rp) ExpdtEk(NumNS, NmSpn)
      real(rp) AMat(NumNS, ND2, NmSpn)       ! Input/Output real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                       ! Starting time point
		integer(8) time2                       ! Ending   time point
      integer SpnInd, I1, I2
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtFFT = TimsMtFFT + 1
      call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate ExpMat*AMat with spin coupled and decoupled for H_0 ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> NmSpn == 1 and Only spin-up channel __________
!______________________ For DCPDNOSIGN == F case --> NmSpn == 2 and both channels _________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
         do I2 = 1, ND2Run, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Perform the X --> K FFT calculation ______________________________
!________________________________________________________________________________________
            FFTVectRNs(1:NumNS) = AMat(1:NumNS, I2, SpnInd)
            call DFFTW_EXECUTE_DFT_R2C(Plan_XK_RC, FFTVectRNs, FFTVectCDm)
!________________________________________________________________________________________ 	  
!_________________ (1) Calculate Exp(\Lambda(k)) * FFTVectCDm ___________________________
!________________________________________________________________________________________
            do I1 = 1, FFTDm, +1
               FFTVectCDm(I1) = ExpdtEk(I1, SpnInd) * FFTVectCDm(I1)
            enddo
!________________________________________________________________________________________ 	  
!_________________ (2) Perform the K --> X FFT calculation ______________________________
!________________________________________________________________________________________
            FFTVectRNs = 0.0_rp
            call DFFTW_EXECUTE_DFT_C2R(Plan_KX_CR, FFTVectCDm, FFTVectRNs)
            AMat(1:NumNS, I2, SpnInd) = FFTVectRNs(1:NumNS)
         enddo
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtFFT = TimeMtFFT + TimeIntrvl(time1, time2)
		
   end subroutine ExpdtKMultFFTEgVl_R
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine FFTEgVlMultExpdtK_R(ND1Run, ND1, AMat, ExpdtEk)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FFTEgVlMultExpdtK_R(ND1Run, ND1, AMat, ExpdtEk)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculation as AMat*ExpMat, with spin coupled and decoupled, 
!                   by FFT method.
! KEYWORDS: AMat*ExpMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate AMat*ExpMat with spin coupled or decoupled.
!
!**************************************************************
!**************** Notice! Notice! Notice! *********************
!**************************************************************
!     Since ExpMat = U * exp(\Lambda) * U^+, then we have AMat*ExpMat = AMat * U * exp(\Lambda) * U^+
!     However, U^+ means the R-->K FFT, while U means the K-->R FFT. So we need to calculate AMat*ExpMat by 
!                FFT as following:
!         (0) First Calculate ExpMat^+ * AMat^+ = U * exp(\Lambda) * U^+ * AMat^+ by the FFT method;
!         (1) Obtain AMat*ExpMat = (ExpMat^+ * AMat^+)^+.
!
!     Input: ND1Run  --> Dimension that will be run;
!            ND1     --> Dimension of AMat;
!            AMat    --> Input matrix;
!            ExpdtEk --> Exp(-dt*E_k) as the diagonal of the ExpMat matrix;
!            
!     Outpt: AMat    --> Finally = AMat*ExpMat.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use RealPrecsn
      use FFTSetting
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1Run, ND1
      real(rp) AMat(ND1, NumNS, NmSpn)       ! Input/Output real matrix
      real(rp) ExpdtEk(NumNS, NmSpn)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                       ! Starting time point
		integer(8) time2                       ! Ending   time point
      integer SpnInd, I1, I2
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtFFT = TimsMtFFT + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate AMat*ExpMat with spin coupled and decoupled for H_0 ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> NmSpn == 1 and Only spin-up channel __________
!______________________ For DCPDNOSIGN == F case --> NmSpn == 2 and both channels _________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, ND1Run, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Perform the X --> K FFT calculation ______________________________
!________________________________________________________________________________________
            FFTVectRNs(1:NumNS) = AMat(I1, 1:NumNS, SpnInd)
            call DFFTW_EXECUTE_DFT_R2C(Plan_XK_RC, FFTVectRNs, FFTVectCDm)
!________________________________________________________________________________________ 	  
!_________________ (1) Calculate FFTVectCDm * Exp(\Lambda(k)) ___________________________
!________________________________________________________________________________________
            do I2 = 1, FFTDm, +1
               FFTVectCDm(I2) = FFTVectCDm(I2) * ExpdtEk(I2, SpnInd)
            enddo
!________________________________________________________________________________________ 	  
!_________________ (2) Perform the K --> X FFT calculation ______________________________
!________________________________________________________________________________________
            FFTVectRNs = 0.0_rp
            call DFFTW_EXECUTE_DFT_C2R(Plan_KX_CR, FFTVectCDm, FFTVectRNs)
            AMat(I1, 1:NumNS, SpnInd) = FFTVectRNs(1:NumNS)
         enddo
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtFFT = TimeMtFFT + TimeIntrvl(time1, time2)
      
   end subroutine FFTEgVlMultExpdtK_R
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ExpdtKMultFFTEgVl_C(ExpdtEk, AMat, ND2, ND2Run)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ExpdtKMultFFTEgVl_C(ExpdtEk, AMat, ND2, ND2Run)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculation as ExpMat*AMat, with spin coupled and decoupled, 
!                     by FFT method.
! KEYWORDS: ExpMat*AMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate ExpMat*AMat with spin coupled or decoupled.
!
!     Input: ExpdtEk --> Exp(-dt*E_k) as the diagonal of the ExpMat matrix;
!            AMat    --> Input matrix;
!            ND2     --> Dimension of AMat;
!            ND2Run  --> Dimension that will be run;
! 
!     Outpt: AMat    --> Finally = ExpMat*AMat.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use RealPrecsn
      use FFTSetting
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND2, ND2Run
      real(rp) ExpdtEk(NumNS, NmSpn)
      complex(rp) AMat(NumNS, ND2, NmSpn)    ! Input/Output complex matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                       ! Starting time point
		integer(8) time2                       ! Ending   time point
      integer SpnInd, I1, I2
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtFFT = TimsMtFFT + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate ExpMat*AMat with spin coupled and decoupled for H_0 ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> NmSpn == 1 and Only spin-up channel __________
!______________________ For DCPDNOSIGN == F case --> NmSpn == 2 and both channels _________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
         do I2 = 1, ND2Run, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Perform the X --> K FFT calculation ______________________________
!________________________________________________________________________________________
            FFTVcApCNs(1:NumNS) = AMat(1:NumNS, I2, SpnInd)
            call DFFTW_EXECUTE_DFT(Plan_XK_CC, FFTVcApCNs, FFTVcBtCNs)
!________________________________________________________________________________________ 	  
!_________________ (1) Calculate Exp(\Lambda(k)) * FFTVcBtCNs ___________________________
!________________________________________________________________________________________
            do I1 = 1, NumNS, +1
               FFTVcBtCNs(I1) = ExpdtEk(I1, SpnInd) * FFTVcBtCNs(I1)
            enddo
!________________________________________________________________________________________ 	  
!_________________ (2) Perform the K --> X FFT calculation ______________________________
!________________________________________________________________________________________
            FFTVcApCNs = rp_Zzero
            call DFFTW_EXECUTE_DFT(Plan_KX_CC, FFTVcBtCNs, FFTVcApCNs)
            AMat(1:NumNS, I2, SpnInd) = FFTVcApCNs(1:NumNS)
         enddo
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtFFT = TimeMtFFT + TimeIntrvl(time1, time2)
		
   end subroutine ExpdtKMultFFTEgVl_C
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine FFTEgVlMultExpdtK_C(ND1Run, ND1, AMat, ExpdtEk)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FFTEgVlMultExpdtK_C(ND1Run, ND1, AMat, ExpdtEk)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the calculation as AMat*ExpMat, with spin coupled and decoupled, 
!                   by FFT method.
! KEYWORDS: AMat*ExpMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate AMat*ExpMat with spin coupled or decoupled.
!
!**************************************************************
!**************** Notice! Notice! Notice! *********************
!**************************************************************
!     Since ExpMat = U * exp(\Lambda) * U^+, then we have AMat*ExpMat = AMat * U * exp(\Lambda) * U^+
!     However, U^+ means the R-->K FFT, while U means the K-->R FFT. So we need to calculate AMat*ExpMat by 
!                FFT as following:
!         (0) First Calculate ExpMat^+ * AMat^+ = U * exp(\Lambda) * U^+ * AMat^+ by the FFT method;
!         (1) Obtain AMat*ExpMat = (ExpMat^+ * AMat^+)^+.
!
!     Input: ND1Run  --> Dimension that will be run;
!            ND1     --> Dimension of AMat;
!            AMat    --> Input matrix;
!            ExpdtEk --> Exp(-dt*E_k) as the diagonal of the ExpMat matrix;
!            
!     Outpt: AMat    --> Finally = AMat*ExpMat.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use RealPrecsn
      use FFTSetting
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1Run, ND1
      complex(rp) AMat(ND1, NumNS, NmSpn)    ! Input/Output complex matrix
      real(rp) ExpdtEk(NumNS, NmSpn)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                       ! Starting time point
		integer(8) time2                       ! Ending   time point
      integer SpnInd, I1, I2
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtFFT = TimsMtFFT + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________ Calculate AMat*ExpMat with spin coupled and decoupled for H_0 ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> NmSpn == 1 and Only spin-up channel __________
!______________________ For DCPDNOSIGN == F case --> NmSpn == 2 and both channels _________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, ND1Run, +1
!________________________________________________________________________________________ 	  
!_________________ (0) Perform the X --> K FFT calculation ______________________________
!________________________________________________________________________________________
            FFTVcApCNs(1:NumNS) = conjg(AMat(I1, 1:NumNS, SpnInd))
            call DFFTW_EXECUTE_DFT(Plan_XK_CC, FFTVcApCNs, FFTVcBtCNs)
!________________________________________________________________________________________ 	  
!_________________ (1) Calculate FFTVcBtCNs * Exp(\Lambda(k)) ___________________________
!________________________________________________________________________________________
            do I2 = 1, NumNS, +1
               FFTVcBtCNs(I2) = FFTVcBtCNs(I2) * ExpdtEk(I2, SpnInd)
            enddo
!________________________________________________________________________________________ 	  
!_________________ (2) Perform the K --> X FFT calculation ______________________________
!________________________________________________________________________________________
            FFTVcApCNs = rp_Zzero
            call DFFTW_EXECUTE_DFT(Plan_KX_CC, FFTVcBtCNs, FFTVcApCNs)
            AMat(I1, 1:NumNS, SpnInd) = conjg(FFTVcApCNs(1:NumNS))
         enddo
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtFFT = TimeMtFFT + TimeIntrvl(time1, time2)
		
   end subroutine FFTEgVlMultExpdtK_C
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
!########################################################################################################################
!########################################################################################################################
!######################################## Transform wavefunctions between real and momentum space #######################
!######################################## Transform wavefunctions between real and momentum space #######################
!######################################## Transform wavefunctions between real and momentum space #######################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine RghtWvfc_r2k_R(AMat, ND2, ND2Run, BMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RghtWvfc_r2k_R(AMat, ND2, ND2Run, BMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to transform the wavefunction at right side from r space to k space.
! KEYWORDS: Transform wavefunction from r to k space, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Right side Wavefunction transformation from r to k space by FFT method. 
!
!     Input: AMat   --> Input matrix under r space basis, real;
!            ND2    --> Dimension of AMat;
!            ND2Run --> Dimension that will be run;
! 
!     Outpt: BMat --> Output matrix under k space basis, complex;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use RealPrecsn
      use FFTSetting
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND2, ND2Run
      real(rp)    AMat(NumNS, ND2   , NmSpn)     ! Input  real    matrix
      complex(rp) BMat(NumNS, ND2Run, NmSpn)     ! Output complex matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                           ! Starting time point
      integer(8) time2                           ! Ending   time point
      integer SpnInd, I1, I1p, I1m, I2, Ix, Iy   ! Iteration Integers
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtFFT = TimsMtFFT + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Transform wavefunction from r to k space _____________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For both DCPDNOSIGN == T and DCPDNOSIGN == F case _________________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
         do I2 = 1, ND2Run, +1
!____________________________________________________________________________ 	  
!________________ [0] Prepare the input vector for FFT ______________________
!____________________________________________________________________________
            FFTVectRNs(1:NumNS) = AMat(1:NumNS, I2, SpnInd)
!____________________________________________________________________________ 	  
!________________ [1] Execute the FFT operation _____________________________
!____________________________________________________________________________
            call DFFTW_EXECUTE_DFT_R2C(Plan_XK_RC, FFTVectRNs, FFTVectCDm)
            FFTVectCDm = FFTVectCDm / sqrt(dble(NumNS))
!____________________________________________________________________________ 	  
!________________ [2] Restore (NumNS) vector from FFTVectCDm(FFTDm) _________
!____________________________________________________________________________
            ! Ix = 0, Iy \in [0, NumL2-1]
            Ix = 0
            do Iy = 0, NumL2-1, +1
               I1  = Iy*(NumL1/2+1) + Ix + 1
               I1p = Iy*NumL1 + Ix + 1
               BMat(I1p, I2, SpnInd) = FFTVectCDm(I1)
            enddo
            ! Iy = 0, Ix \in [1, NumL1-1]
            Iy = 0
            do Ix = 1, NumL1/2, +1
               I1  = Ix + 1
               I1p =         Ix + 1
               I1m = NumL1 - Ix + 1
               BMat(I1p, I2, SpnInd) =       FFTVectCDm(I1)
               BMat(I1m, I2, SpnInd) = conjg(FFTVectCDm(I1))
            enddo
            ! Ix \in [1, NumL1-1], Iy \in [1, NumL2-1]
            do Iy = 1, NumL2-1, +1
               do Ix = 1, NumL1/2, +1
                  I1  = Iy *(NumL1/2+1) + Ix + 1
                  I1p =        Iy * NumL1 +         Ix + 1
                  I1m = (NumL2-Iy)* NumL1 + NumL1 - Ix + 1
                  BMat(I1p, I2, SpnInd) =       FFTVectCDm(I1)
                  BMat(I1m, I2, SpnInd) = conjg(FFTVectCDm(I1))
               enddo
            enddo
         enddo
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtFFT = TimeMtFFT + TimeIntrvl(time1, time2)
		
   end subroutine RghtWvfc_r2k_R
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine LeftWvfc_r2k_R(ND1Run, ND1, AMat, BMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  LeftWvfc_r2k_R(ND1Run, ND1, AMat, BMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to transform the wavefunction at right side from r space to k space.
! KEYWORDS: Transform wavefunction from r to k space, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Left side Wavefunction transformation from r to k space by FFT method. 
!
!     Input: ND1Run --> Dimension that will be run;
!            ND1    --> Dimension of AMat;
!            AMat   --> Input matrix under r space basis, real;
! 
!     Outpt: BMat --> Output matrix under k space basis, complex;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use RealPrecsn
      use FFTSetting
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1, ND1Run
      real(rp)    AMat(ND1   , NumNS, NmSpn)     ! Input  real    matrix
      complex(rp) BMat(ND1Run, NumNS, NmSpn)     ! Output complex matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                           ! Starting time point
      integer(8) time2                           ! Ending   time point
      integer SpnInd, I1, I2, I2p, I2m, Ix, Iy   ! Iteration Integers
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtFFT = TimsMtFFT + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Transform wavefunction from r to k space _____________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For both DCPDNOSIGN == T and DCPDNOSIGN == F case _________________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, ND1Run, +1
!____________________________________________________________________________ 	  
!________________ [0] Prepare the input vector for FFT ______________________
!____________________________________________________________________________
            FFTVectRNs(1:NumNS) = AMat(I1, 1:NumNS, SpnInd)
!____________________________________________________________________________ 	  
!________________ [1] Execute the FFT operation _____________________________
!____________________________________________________________________________
            call DFFTW_EXECUTE_DFT_R2C(Plan_XK_RC, FFTVectRNs, FFTVectCDm)
            FFTVectCDm = FFTVectCDm / sqrt(dble(NumNS))
!____________________________________________________________________________ 	  
!________________ [2] Restore (NumNS) vector from FFTVectCDm(FFTDm) _________
!____________________________________________________________________________
            ! Ix = 0, Iy \in [0, NumL2-1]
            Ix = 0
            do Iy = 0, NumL2-1, +1
               I2  = Iy*(NumL1/2+1) + Ix + 1
               I2p = Iy*NumL1 + Ix + 1
               BMat(I1, I2p, SpnInd) = conjg(FFTVectCDm(I2))
            enddo
            ! Iy = 0, Ix \in [1, NumL1-1]
            Iy = 0
            do Ix = 1, NumL1/2, +1
               I2  = Ix + 1
               I2p =         Ix + 1
               I2m = NumL1 - Ix + 1
               BMat(I1, I2p, SpnInd) = conjg(FFTVectCDm(I2))
               BMat(I1, I2m, SpnInd) =       FFTVectCDm(I2)
            enddo
            ! Ix \in [1, NumL1-1], Iy \in [1, NumL2-1]
            do Iy = 1, NumL2-1, +1
               do Ix = 1, NumL1/2, +1
                  I2  = Iy *(NumL1/2+1) + Ix + 1
                  I2p =        Iy * NumL1 +         Ix + 1
                  I2m = (NumL2-Iy)* NumL1 + NumL1 - Ix + 1
                  BMat(I1, I2p, SpnInd) = conjg(FFTVectCDm(I2))
                  BMat(I1, I2m, SpnInd) =       FFTVectCDm(I2)
               enddo
            enddo
         enddo
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtFFT = TimeMtFFT + TimeIntrvl(time1, time2)
		
   end subroutine LeftWvfc_r2k_R
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine RghtWvfc_r2k_C(AMat, ND2, ND2Run, BMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RghtWvfc_r2k_C(AMat, ND2, ND2Run, BMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to transform the wavefunction at right side from r space to k space.
! KEYWORDS: Transform wavefunction from r to k space, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Right side Wavefunction transformation from r to k space by FFT method. 
!
!     Input: AMat   --> Input matrix under r space basis, complex;
!            ND2    --> Dimension of AMat;
!            ND2Run --> Dimension that will be run;
! 
!     Outpt: BMat --> Output matrix under k space basis, complex;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use RealPrecsn
      use FFTSetting
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND2, ND2Run
      complex(rp) AMat(NumNS, ND2   , NmSpn)     ! Input  complex matrix
      complex(rp) BMat(NumNS, ND2Run, NmSpn)     ! Output complex matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                           ! Starting time point
      integer(8) time2                           ! Ending   time point
      integer SpnInd, I1, I2                     ! Iteration Integers
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtFFT = TimsMtFFT + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Transform wavefunction from r to k space _____________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For both DCPDNOSIGN == T and DCPDNOSIGN == F case _________________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
         do I2 = 1, ND2Run, +1
!____________________________________________________________________________ 	  
!________________ [0] Prepare the input vector for FFT ______________________
!____________________________________________________________________________
            FFTVcApCNs(1:NumNS) = AMat(1:NumNS, I2, SpnInd)
!____________________________________________________________________________ 	  
!________________ [1] Execute the FFT operation _____________________________
!____________________________________________________________________________
            call DFFTW_EXECUTE_DFT(Plan_XK_CC, FFTVcApCNs, FFTVcBtCNs)
            FFTVcBtCNs = FFTVcBtCNs / sqrt(dble(NumNS))
!____________________________________________________________________________ 	  
!________________ [2] Store FFTVcBtCNs to the BMat matrix ______________________
!____________________________________________________________________________
            BMat(1:NumNS, I2, SpnInd) = FFTVcBtCNs(1:NumNS)
         enddo
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtFFT = TimeMtFFT + TimeIntrvl(time1, time2)
		
   end subroutine RghtWvfc_r2k_C
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine LeftWvfc_r2k_C(ND1Run, ND1, AMat, BMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  LeftWvfc_r2k_C(ND1Run, ND1, AMat, BMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to transform the wavefunction at right side from r space to k space.
! KEYWORDS: Transform wavefunction from r to k space, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Left side Wavefunction transformation from r to k space by FFT method. 
!
!     Input: ND1Run --> Dimension that will be run;
!            ND1    --> Dimension of AMat;
!            AMat   --> Input matrix under r space basis, real;
! 
!     Outpt: BMat --> Output matrix under k space basis, complex;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
		use RealPrecsn
      use FFTSetting
      use TimeRecord
      use QMCTimeRec
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1, ND1Run
      complex(rp) AMat(ND1   , NumNS, NmSpn)     ! Input  complex matrix
      complex(rp) BMat(ND1Run, NumNS, NmSpn)     ! Output complex matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1                           ! Starting time point
      integer(8) time2                           ! Ending   time point
      integer SpnInd, I1, I2                     ! Iteration Integers
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtFFT = TimsMtFFT + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Transform wavefunction from r to k space _____________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For both DCPDNOSIGN == T and DCPDNOSIGN == F case _________________________
!**************************************************************************************************
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, ND1Run, +1
!____________________________________________________________________________ 	  
!________________ [0] Prepare the input vector for FFT ______________________
!____________________________________________________________________________
            FFTVcApCNs(1:NumNS) = conjg(AMat(I1, 1:NumNS, SpnInd))
!____________________________________________________________________________ 	  
!________________ [1] Execute the FFT operation _____________________________
!____________________________________________________________________________
            call DFFTW_EXECUTE_DFT(Plan_XK_CC, FFTVcApCNs, FFTVcBtCNs)
            FFTVcBtCNs = FFTVcBtCNs / sqrt(dble(NumNS))
!____________________________________________________________________________ 	  
!________________ [2] Store FFTVcBtCNs to the BMat matrix ___________________
!____________________________________________________________________________
            BMat(I1, 1:NumNS, SpnInd) = conjg(FFTVcBtCNs(1:NumNS))
         enddo
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtFFT = TimeMtFFT + TimeIntrvl(time1, time2)
		
   end subroutine LeftWvfc_r2k_C
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$