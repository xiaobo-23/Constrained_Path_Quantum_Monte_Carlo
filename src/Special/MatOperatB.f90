!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform matrix operations including multiplication, inverse, determinant 
!               and UDV decomposition for QMC simulations, based on spin-coupled and spin-decoupled cases. 
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!    All the following subroutines call the internal subroutines of Intel MKL package. There is no leading 
!         dimension problem for the input matrices and there is spin coupled or decoupled problem.
!
!    zMatDet_LogDet --> Subroutine to calculate the log of determinant of complex square matrix;
!    dMatDet_LogDet --> Subroutine to calculate the log of determinant of real    square matrix;

!    zMatInv_LogDet --> Subroutine to calculate the inverse and log(determinant) of complex square matrix;
!    dMatInv_LogDet --> Subroutine to calculate the inverse and log(determinant) of real    square matrix;
!
!    zMatEqSet_Left_LogDet --> Subroutine to calculate XMat*AMat=BMat (or BMat*AMat^{-1}), complex version;
!    zMatEqSet_Rght_LogDet --> Subroutine to calculate AMat*XMat=BMat (or AMat^{-1}*BMat), complex version;
!    dMatEqSet_Left_LogDet --> Subroutine to calculate XMat*AMat=BMat (or BMat*AMat^{-1}), real    version;
!    dMatEqSet_Rght_LogDet --> Subroutine to calculate AMat*XMat=BMat (or AMat^{-1}*BMat), real    version;
!  
!    zMatUDVdcp --> Subroutine to calculate the UDV decomposition of complex matrix;  
!    dMatUDVdcp --> Subroutine to calculate the UDV decomposition of real    matrix;    
!
!    zMatPrdAll --> Subroutine to calculate matrix product A(ND1, ND2)*B(ND2, ND3) for complex matrices;
!    dMatPrdAll --> Subroutine to calculate matrix product A(ND1, ND2)*B(ND2, ND3) for real    matrices;
!
!    zMatDiagnz --> Subroutine used to diagonalize complex Hermitian matrix;
!    dMatDiagnz --> Subroutine used to diagonalize real symmetric    matrix;
!
!    zMat_Compare --> Subroutine used to compare two complex (NumNS, NumNS) matrices;
!    dMat_Compare --> Subroutine used to compare two real    (NumNS, NumNS) matrices.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   
!########################################################################################################################
!########################################################################################################################
!############################################### Calculate Matrix Determinant ###########################################
!############################################### Calculate Matrix Determinant ###########################################
!############################################### Calculate Matrix Determinant ###########################################
!########################################################################################################################
!########################################################################################################################

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatDet_LogDet(NDim, zMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatDet_LogDet(NDim, zMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate log determinant for complex square matrix, based on 
!                  spin-coupled or spin-decoupled cases.
! KEYWORDS: Calculate Matrix determinant for complex matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix determinant, complex version. Based on spin-coupled or spin-decoupled cases of the system.
!
!     Input: NDim --> Dimension of input A matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: LogzDet --> Log of Complex determinant.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer NDim                         ! Dimension of zMat square matrix
      complex(rp) LogzDet                  ! Log Determinant of zMat matrix
		complex(rp) zMat(NDim, NDim)         ! The input zMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix determinant _____________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&   
		TimsMtDet = TimsMtDet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Calculate log of determinant for the matrix _______________________________
!**************************************************************************************************
      LogzDet = rp_Zzero
      call MatDetermtZ_LogDet(NDim, zMat(1, 1), NDim, LogzDet)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix determinant _____________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&    
		call system_clock(time2)
		TimeMtDet = TimeMtDet + TimeIntrvl(time1, time2)
		
   end subroutine zMatDet_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatDet_LogDet(NDim, dMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatDet_LogDet(NDim, dMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate log determinant for real square matrix, based on 
!                  spin-coupled or spin-decoupled cases.
! KEYWORDS: Calculate Matrix determinant for real matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix determinant, real version. Based on spin-coupled or spin-decoupled cases of the system.
!
!     Input: NDim --> Dimension of input A matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: LogzDet --> Log of Complex determinant.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer NDim                         ! Dimension of dMat square matrix
      complex(rp) LogzDet                  ! Log Determinant of dMat matrix
		real(rp) dMat(NDim, NDim)            ! The input dMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix determinant _____________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&   
		TimsMtDet = TimsMtDet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Calculate log of determinant for the matrix _______________________________
!**************************************************************************************************
      LogzDet = rp_Zzero
      call MatDetermtR_LogDet(NDim, dMat(1, 1), NDim, LogzDet)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix determinant _____________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&    
		call system_clock(time2)
		TimeMtDet = TimeMtDet + TimeIntrvl(time1, time2)
		
   end subroutine dMatDet_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
   
   
   
!########################################################################################################################
!########################################################################################################################
!############################################### Calculate Matrix Inverse ###############################################
!############################################### Calculate Matrix Inverse ###############################################
!############################################### Calculate Matrix Inverse ###############################################
!########################################################################################################################
!########################################################################################################################   

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatInv_LogDet(NDim, zMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatInv_LogDet(NDim, zMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate inverse and log determinant for complex square matrix, based on 
!                  spin-coupled or spin-decoupled cases.
! KEYWORDS: Calculate Matrix inverse and log determinant for complex matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix Inverse and log determinant, complex version. Based on spin-coupled or spin-decoupled 
!             cases of the system.
!
!     Input: NDim --> Dimension of input A matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: LogzDet --> Log of Complex determinant.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer NDim                         ! Dimension of zMat square matrix
      complex(rp) LogzDet                  ! Log Determinant of zMat matrix
		complex(rp) zMat(NDim, NDim)         ! The input zMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix inverse _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsMtInv = TimsMtInv + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Calculate inverse and log of determinant for the matrix ___________________
!**************************************************************************************************
      LogzDet = rp_Zzero
      call MatInversZ_LogDet(NDim, zMat(1, 1), NDim, LogzDet)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix inverse _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%& 
		call system_clock(time2)
		TimeMtInv = TimeMtInv + TimeIntrvl(time1, time2)
		
   end subroutine zMatInv_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatInv_LogDet(NDim, dMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatInv_LogDet(NDim, dMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate inverse and log determinant for real square matrix, based on 
!                  spin-coupled or spin-decoupled cases.
! KEYWORDS: Calculate Matrix inverse and log determinant for real matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix Inverse and log determinant, real version. Based on spin-coupled or spin-decoupled cases 
!                 of the system.
!
!     Input: NDim --> Dimension of input A matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: LogzDet --> Log of Complex determinant.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer NDim, DmUp, DmDw             ! Dimension of dMat square matrix
      complex(rp) LogzDet                  ! Log Determinant of dMat matrix
		real(rp) dMat(NDim, NDim)            ! The input dMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
      
      integer Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetUp, LogzDetDw
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix inverse _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsMtInv = TimsMtInv + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 1. Calculate inverse and log of determinant for the matrix ___________________
!**************************************************************************************************
      LogzDet = rp_Zzero
      call MatInversR_LogDet(DmUp, dMat(1, 1), NDim, LogzDet)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix inverse _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%& 
		call system_clock(time2)
		TimeMtInv = TimeMtInv + TimeIntrvl(time1, time2)
		
   end subroutine dMatInv_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!########################################################################################################################
!########################################################################################################################
!############################################### Solve Linear Equation Set ##############################################
!############################################### Solve Linear Equation Set ##############################################
!############################################### Solve Linear Equation Set ##############################################
!########################################################################################################################
!######################################################################################################################## 
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatEqSet_Left_LogDet(ND1, ND2, AMat, BMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatEqSet_Left_LogDet(ND1, ND2, AMat, BMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve XMat(ND1, ND2)*AMat(ND2, ND2)=BMat(ND1, ND2) or calculate BMat*AMat^{-1}.
! KEYWORDS: Solve XMat*AMat=BMat or Calculate BMat*AMat^{-1}, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat*AMat=BMat or Calculate BMat*AMat^{-1}, by the method of performing LU decomposition and 
!                solving linear equation set.
!
!     Input: ND1  --> Dimension of input matrices AMat(ND2, ND2) and BMat(ND1, ND2);
!            ND2  --> Dimension of input matrices AMat(ND2, ND2) and BMat(ND1, ND2);
!            AMat --> Input AMat(ND2, ND2) matrix;
!            BMat --> Input BMat(ND1, ND2) matrix;
!
!     Outpt: BMat    --> Output solution XMat matrix or result of BMat*AMat^{-1};
!            LogzDet --> Log of Complex determinant of AMat matrix.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1, ND2
      complex(rp) LogzDet
		complex(rp) AMat(ND2, ND2) 
      complex(rp) BMat(ND1, ND2) 
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsEqSet = TimsEqSet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Solve XMat from XMat(ND1, ND2)*AMat(ND2, ND2)=BMat(ND1, ND2) ______________
!**************************************************************************************************
      LogzDet = rp_Zzero
      call SlvLnEqSetZ_Left_LogDet(ND1, ND2, AMat(1, 1), ND2, BMat(1, 1), ND1, LogzDet)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeEqSet = TimeEqSet + TimeIntrvl(time1, time2)
		
   end subroutine zMatEqSet_Left_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatEqSet_Rght_LogDet(ND1, ND2, AMat, BMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatEqSet_Rght_LogDet(ND1, ND2, AMat, BMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve AMat(ND1, ND1)*XMat(ND1, ND2)=BMat(ND1, ND2) or calculate AMat^{-1}*BMat.
! KEYWORDS: Solve AMat*XMat=BMat or Calculate AMat^{-1}*BMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat*XMat=BMat or Calculate AMat^{-1}*BMat, by the method of performing LU decomposition and 
!                  solving linear equation set.
!
!     Input: ND1  --> Dimension of input matrices AMat(ND1, ND1) and BMat(ND1, ND2);
!            ND2  --> Dimension of input matrices AMat(ND1, ND1) and BMat(ND1, ND2);
!            AMat --> Input AMat(ND1, ND1) matrix;
!            BMat --> Input BMat(ND1, ND2) matrix;
!
!     Outpt: BMat    --> Output solution XMat matrix or result of AMat^{-1}*BMat;
!            LogzDet --> Log of Complex determinant of AMat matrix.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1, ND2
      complex(rp) LogzDet
		complex(rp) AMat(ND1, ND1) 
      complex(rp) BMat(ND1, ND2) 
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsEqSet = TimsEqSet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Solve XMat from AMat(ND1, ND1)*XMat(ND1, ND2)=BMat(ND1, ND2) ______________
!**************************************************************************************************
      LogzDet = rp_Zzero
      call SlvLnEqSetZ_Rght_LogDet(ND1, ND2, AMat(1, 1), ND1, BMat(1, 1), ND1, LogzDet)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeEqSet = TimeEqSet + TimeIntrvl(time1, time2)
		
   end subroutine zMatEqSet_Rght_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatEqSet_Left_LogDet(ND1, ND2, AMat, BMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatEqSet_Left_LogDet(ND1, ND2, AMat, BMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve XMat(ND1, ND2)*AMat(ND2, ND2)=BMat(ND1, ND2) or calculate BMat*AMat^{-1}.
! KEYWORDS: Solve XMat*AMat=BMat or Calculate BMat*AMat^{-1}, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat*AMat=BMat or Calculate BMat*AMat^{-1}, by the method of performing LU decomposition and 
!                   solving linear equation set.
!
!     Input: ND1  --> Dimension of input matrices AMat(ND2, ND2) and BMat(ND1, ND2);
!            ND2  --> Dimension of input matrices AMat(ND2, ND2) and BMat(ND1, ND2);
!            AMat --> Input AMat(ND2, ND2) matrix;
!            BMat --> Input BMat(ND1, ND2) matrix;
!
!     Outpt: BMat    --> Output solution XMat matrix or result of BMat*AMat^{-1};
!            LogzDet --> Log of real determinant of AMat matrix.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1, ND2
      complex(rp) LogzDet
		real(rp) AMat(ND2, ND2) 
      real(rp) BMat(ND1, ND2) 
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsEqSet = TimsEqSet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Solve XMat from XMat(ND1, ND2)*AMat(ND2, ND2)=BMat(ND1, ND2) ______________
!**************************************************************************************************
      LogzDet = rp_Zzero
      call SlvLnEqSetR_Left_LogDet(ND1, ND2, AMat(1, 1), ND2, BMat(1, 1), ND1, LogzDet)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeEqSet = TimeEqSet + TimeIntrvl(time1, time2)
		
   end subroutine dMatEqSet_Left_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatEqSet_Rght_LogDet(ND1, ND2, AMat, BMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatEqSet_Rght_LogDet(ND1, ND2, AMat, BMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve AMat(ND1, ND1)*XMat(ND1, ND2)=BMat(ND1, ND2) or calculate AMat^{-1}*BMat.
! KEYWORDS: Solve AMat*XMat=BMat or Calculate AMat^{-1}*BMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat*XMat=BMat or Calculate AMat^{-1}*BMat, by the method of performing LU decomposition and 
!               solving linear equation set.
!
!     Input: ND1  --> Dimension of input matrices AMat(ND1, ND1) and BMat(ND1, ND2);
!            ND2  --> Dimension of input matrices AMat(ND1, ND1) and BMat(ND1, ND2);
!            AMat --> Input AMat(ND1, ND1) matrix;
!            BMat --> Input BMat(ND1, ND2) matrix;
!
!     Outpt: BMat    --> Output solution XMat matrix or result of AMat^{-1}*BMat;
!            LogzDet --> Log of real determinant of AMat matrix.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1, ND2
      complex(rp) LogzDet
		real(rp) AMat(ND1, ND1) 
      real(rp) BMat(ND1, ND2) 
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsEqSet = TimsEqSet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Solve XMat from AMat(ND1, ND1)*XMat(ND1, ND2)=BMat(ND1, ND2) ______________
!**************************************************************************************************
      LogzDet = rp_Zzero
      call SlvLnEqSetR_Rght_LogDet(ND1, ND2, AMat(1, 1), ND1, BMat(1, 1), ND1, LogzDet)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeEqSet = TimeEqSet + TimeIntrvl(time1, time2)
		
   end subroutine dMatEqSet_Rght_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
   
!########################################################################################################################
!########################################################################################################################
!############################################### UDV Decomposition for Matrix ###########################################
!############################################### UDV Decomposition for Matrix ###########################################
!############################################### UDV Decomposition for Matrix ###########################################
!########################################################################################################################
!######################################################################################################################## 
 
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatUDVdcp(NDim, zMat, UMat, DVec, VMat, IfUDVCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatUDVdcp(NDim, zMat, UMat, DVec, VMat, IfUDVCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate perform the UDV decomposition for zMat(NDim, NDim) matrix, based on 
!                  spin-coupled or spin-decoupled cases.
! KEYWORDS: Calculate UDV decomposition for complex matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate UDV decomposition for complex matrix. Based on spin-coupled or spin-decoupled cases of the system.
!
!     Input: NDim --> The second dimension of input zMat matrix;
!            zMat --> Input complex matrix;
!            IfUDVCheck --> Whether to perform the final UDV check.
!
!     Outpt: UMat --> Result U matrix;
!            DVec --> Result D vector;
!            VMat --> Result V matrix.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
      use RealPrecsn
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfUDVCheck
		integer NDim                         ! Dimension of zMat square matrix
		complex(rp) zMat(NDim, NDim)         ! The input zMat matrix
      complex(rp) UMat(NDim, NDim)         ! The input zMat matrix
      real(rp)    DVec(NDim)               ! The input zMat matrix
      complex(rp) VMat(NDim, NDim)         ! The input zMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in calculating UDV decomposition __________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
      TimsUDVOt = TimsUDVOt + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating UDV decomposition _____________________________
!______________________________________________________________________________________________________________ 
#ifdef QRSTABLIZE
      call UDVdcmpsZ_QRPivot(NDim, NDim, zMat(1, 1), NDim, UMat(1, 1), NDim, DVec(1), VMat(1, 1), NDim, IfUDVCheck)
#else
      call UDVdcmpsZ_SVD(    NDim, NDim, zMat(1, 1), NDim, UMat(1, 1), NDim, DVec(1), VMat(1, 1), NDim, IfUDVCheck)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in calculating UDV decomposition __________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeUDVOt = TimeUDVOt + TimeIntrvl(time1, time2)
		
   end subroutine zMatUDVdcp
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatUDVdcp(NDim, dMat, UMat, DVec, VMat, IfUDVCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatUDVdcp(NDim, dMat, UMat, DVec, VMat, IfUDVCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate perform the UDV decomposition for dMat(NDim, NDim) matrix, based on 
!                  spin-coupled or spin-decoupled cases.
! KEYWORDS: Calculate UDV decomposition for real matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate UDV decomposition for real matrix. Based on spin-coupled or spin-decoupled cases of the system.
!
!     Input: NDim --> The second dimension of input dMat matrix;
!            dMat --> Input real matrix;
!            IfUDVCheck --> Whether to perform the final UDV check.
!
!     Outpt: UMat --> Result U matrix;
!            DVec --> Result D vector;
!            VMat --> Result V matrix.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      logical IfUDVCheck
		integer NDim                      ! Dimension of dMat square matrix
		real(rp) dMat(NDim, NDim)         ! The input dMat matrix
      real(rp) UMat(NDim, NDim)         ! The Output UMat matrix
      real(rp) DVec(NDim)               ! The Output D    vector
      real(rp) VMat(NDim, NDim)         ! The Output VMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in calculating UDV decomposition __________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
      TimsUDVOt = TimsUDVOt + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating UDV decomposition _____________________________
!______________________________________________________________________________________________________________ 
#ifdef QRSTABLIZE
      call UDVdcmpsR_QRPivot(NDim, NDim, dMat(1, 1), NDim, UMat(1, 1), NDim, DVec(1), VMat(1, 1), NDim, IfUDVCheck)
#else
      call UDVdcmpsR_SVD(    NDim, NDim, dMat(1, 1), NDim, UMat(1, 1), NDim, DVec(1), VMat(1, 1), NDim, IfUDVCheck)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in calculating UDV decomposition __________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeUDVOt = TimeUDVOt + TimeIntrvl(time1, time2)
		
   end subroutine dMatUDVdcp
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
   
!########################################################################################################################
!########################################################################################################################
!############################################### Calculate Matrix Product ###############################################
!############################################### Calculate Matrix Product ###############################################
!############################################### Calculate Matrix Product ###############################################
!########################################################################################################################
!######################################################################################################################## 
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatPrdAll(TransA, TransB, ND1, ND3, ND2, Ztp1, AMat, LDA, BMat, LDB, Ztp2, CMat, LDC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatPrdAll(TransA, TransB, ND1, ND3, ND2, Ztp1, AMat, LDA, BMat, LDB, Ztp2, CMat, LDC)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of AMat(ND1, ND2) and BMat(ND2, ND3).
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, complex version. 
!
!     CMat(ND1, ND3) = AMat(ND1, ND2) * BMat(ND2, ND3).
!
!     Input: (none); Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      character TransA, TransB
      integer ND1, ND2, ND3
      integer LDA, LDB, LDC
      complex(kind=kind(0.d0)) Ztp1, Ztp2
		complex(kind=kind(0.d0)) AMat(LDA, *)     ! Input AMat(ND1, ND2) matrix
		complex(kind=kind(0.d0)) BMat(LDB, *)     ! Input BMat(ND2, ND3) matrix
		complex(kind=kind(0.d0)) CMat(LDC, *)     ! Ouput CMat(ND1, ND3) matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Calculate Matrix product for AMat and BMat ________________________________
!**************************************************************************************************
      call ZGEMM(TransA, TransB, ND1, ND3, ND2, Ztp1, AMat, LDA, BMat, LDB, Ztp2, CMat, LDC)  
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine zMatPrdAll
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatPrdAll(TransA, TransB, ND1, ND3, ND2, Rtp1, AMat, LDA, BMat, LDB, Rtp2, CMat, LDC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatPrdAll(TransA, TransB, ND1, ND3, ND2, Rtp1, AMat, LDA, BMat, LDB, Rtp2, CMat, LDC)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of AMat(ND1, ND2) and BMat(ND2, ND3).
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, real version. 
!
!     CMat(ND1, ND3) = AMat(ND1, ND2) * BMat(ND2, ND3).
!
!     Input: (none); Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      character TransA, TransB
      integer ND1, ND2, ND3
      integer LDA, LDB, LDC
      real(kind=kind(0.d0)) Rtp1, Rtp2
		real(kind=kind(0.d0)) AMat(LDA, *)     ! Input AMat(ND1, ND2) matrix
		real(kind=kind(0.d0)) BMat(LDB, *)     ! Input BMat(ND2, ND3) matrix
		real(kind=kind(0.d0)) CMat(LDC, *)     ! Ouput CMat(ND1, ND3) matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Calculate Matrix product for AMat and BMat ________________________________
!**************************************************************************************************
      call DGEMM(TransA, TransB, ND1, ND3, ND2, Rtp1, AMat, LDA, BMat, LDB, Rtp2, CMat, LDC)  
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine dMatPrdAll
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMat_VecVecPrd(ND1, ND2, Ztp1, VecA, IncA, VecB, IncB, CMat, LDC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMat_VecVecPrd(ND1, ND2, Ztp1, VecA, IncA, VecB, IncB, CMat, LDC)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate CMat = CMat + Ztp1*VecA*VecB'.
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate CMat = CMat + Ztp1*VecA*VecB', complex version. 
!
!     Input: (none); Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1, ND2, IncA, IncB, LDC
      complex(kind=kind(0.d0)) Ztp1
		complex(kind=kind(0.d0)) VecA(1)            ! Input VecA(ND1     ) matrix
		complex(kind=kind(0.d0)) VecB(1)            ! Input VecB(ND2     ) matrix
		complex(kind=kind(0.d0)) CMat(LDC, *)       ! Ouput CMat(ND1, ND2) matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Calculate Matrix product for AMat and BMat ________________________________
!**************************************************************************************************
      call ZGERU(ND1, ND2, Ztp1, VecA(1), IncA, VecB(1), IncB, CMat(1, 1), LDC)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine zMat_VecVecPrd
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMat_VecVecPrd(ND1, ND2, Rtp1, VecA, IncA, VecB, IncB, CMat, LDC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMat_VecVecPrd(ND1, ND2, Rtp1, VecA, IncA, VecB, IncB, CMat, LDC)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate CMat = CMat + Rtp1*VecA*VecB'.
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate CMat = CMat + Rtp1*VecA*VecB', real version. 
!
!     Input: (none); Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________  
		use TimeRecord
		use QMCTimeRec
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer ND1, ND2, IncA, IncB, LDC
      real(kind=kind(0.d0)) Rtp1
		real(kind=kind(0.d0)) VecA(1)            ! Input VecA(ND1     ) matrix
		real(kind=kind(0.d0)) VecB(1)            ! Input VecB(ND2     ) matrix
		real(kind=kind(0.d0)) CMat(LDC, *)       ! Ouput CMat(ND1, ND2) matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Calculate Matrix product for AMat and BMat ________________________________
!**************************************************************************************************
      call DGER(ND1, ND2, Rtp1, VecA(1), IncA, VecB(1), IncB, CMat(1, 1), LDC)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine dMat_VecVecPrd
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
   
!########################################################################################################################
!########################################################################################################################
!############################################### Diagonalize Hermitian and Symmetric Matrices ###########################
!############################################### Diagonalize Hermitian and Symmetric Matrices ###########################
!############################################### Diagonalize Hermitian and Symmetric Matrices ###########################
!########################################################################################################################
!########################################################################################################################
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatDiagnz(NDim, zMat, EVal, EVec)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatDiagnz(NDim, zMat, EVal, EVec)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to diagonalize zMat(NDim, NDim) matrix which is Hermitian.
! KEYWORDS: Diagonalize Hermitian Matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Diagonalize Hermitian matrix;
!
!     Input:  NDim --> Dimension of the input matrix;
!             zMat --> Input square matrix;
!
!     Output: Eval --> Output eigenvalues;
!             EVec --> Output eigenvector matrix.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NDim
		complex(rp) zMat(NDim, NDim)     ! Input AMat(ND1, ND2) matrix
		real(rp)    EVal(NDim)           ! Input BMat(ND2, ND3) matrix
		complex(rp) EVec(NDim, NDim)     ! Ouput CMat(ND1, ND3) matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer Istat
      integer Info
      integer LWork
      integer LrWork
      
      real(rp), allocatable :: rWork(:)
      complex(rp), allocatable :: Work(:)  
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Diagonalizing Hermitian matrix ____________________________
!______________________________________________________________________________________________________________ 
      LWork  = 2 * NDim - 1
      LrWork = 3 * NDim - 2
      allocate( Work(LWork ), stat=istat)
      allocate(rWork(LrWork), stat=istat)

      Eval = 0.0d0
      call zcopy(NDim*NDim, zMat(1, 1), 1, EVec(1, 1), 1)
				
      Work  = rp_Zzero
      rWork = 0.0d0
      call ZHEEV("V", "U", NDim, EVec, NDim, Eval, Work, LWork, rWork, Info)
      if ( info /= 0 ) then
         write(*, "('zMatDiagnz, error in lapack subroutine ZHEEV!', 2x, 'ierror = ', I4)") info
      end if
      
      if(allocated( Work)) deallocate( Work)
      if(allocated(rWork)) deallocate(rWork)
		
   end subroutine zMatDiagnz
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatDiagnz(NDim, dMat, EVal, EVec)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatDiagnz(NDim, dMat, EVal, EVec)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to diagonalize dMat(NDim, NDim) matrix which is Hermitian.
! KEYWORDS: Diagonalize real symmetric Matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Diagonalize Hermitian matrix;
!
!     Input:  NDim --> Dimension of the input matrix;
!             dMat --> Input square matrix;
!
!     Output: Eval --> Output eigenvalues;
!             EVec --> Output eigenvector matrix.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer DmUp, DmDw, NDim
		real(rp) dMat(NDim, NDim)     ! Input AMat(ND1, ND2) matrix
		real(rp) EVal(NDim)           ! Input BMat(ND2, ND3) matrix
		real(rp) EVec(NDim, NDim)     ! Ouput CMat(ND1, ND3) matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer Istat
      integer Info
      integer LWork
      
      real(rp), allocatable :: Work(:)  
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Diagonalizing Hermitian matrix ____________________________
!______________________________________________________________________________________________________________ 
      LWork = 3 * NDim - 1
      allocate(Work(LWork), stat=istat)

      Eval = 0.0d0
      call dcopy(NDim*NDim, dMat(1, 1), 1, EVec(1, 1), 1)
				
      Work = 0.d0
      call DSYEV("V", "U", NDim, EVec, NDim, Eval, Work, LWork, Info)
      if ( info /= 0 ) then
         write(*, "('dMatDiagnz, error in lapack subroutine DSYEV!', 2x, 'ierror = ', I4)") info
      end if
      
      if(allocated(Work)) deallocate(Work)
		
   end subroutine dMatDiagnz
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
   
!########################################################################################################################
!########################################################################################################################
!############################################### Comparing two matrices #################################################
!############################################### Comparing two matrices #################################################
!############################################### Comparing two matrices #################################################
!########################################################################################################################
!########################################################################################################################

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine dMat_Compare(AMat, BMat, XMaxm, XMean)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  dMat_Compare(AMat, BMat, XMaxm, XMean)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compare two (NumNS, NumNS) real different matrices and calculate their
!             differences in elements.
! KEYWORDS: Compare two matrices, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Compare two matrices, real version. 
!
!     Input: AMat --> Input matrix A(N, M);
!            BMat --> Input matrix A(N, M);
!     
!     Outpt: XMaxm --> The largest absolute difference;
!            XMean --> The average difference in all elements.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		real(rp) XMaxm               ! Maximum difference for matrix elements in A, B matrices
		real(rp) XMean               ! Average difference for matrix elements in A, B matrices
		real(rp) AMat(NumNS, NumNS)  ! Input A matrix
		real(rp) BMat(NumNS, NumNS)  ! Input B matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1      ! Loop integer 
		integer I2      ! Loop integer
		real(rp) Diff   ! Difference between the elements in A and B
      
      real(rp), allocatable :: XMaxmVec(:)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Creating matrix _________________________________________
!______________________________________________________________________________________________________________
		XMaxm = 0.d0
		XMean = 0.d0
      
      allocate(XMaxmVec(NumNS))
      XMaxmVec = 0.d0
      
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2, Diff)
   !$OMP DO REDUCTION(+ : XMean)
		do I1 = 1, NumNS, +1
			do I2 = 1, NumNS, +1
				Diff = dabs(BMat(I1, I2) - AMat(I1, I2))
				if(Diff .gt. XMaxmVec(I1)) then
					XMaxmVec(I1) = Diff
				end if
				XMean = XMean + Diff
			enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
				
		XMean = XMean / dble(NumNS * NumNS)
      XMaxm = maxval(XMaxmVec)
      
      deallocate(XMaxmVec)
		
	end subroutine dMat_Compare
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine zMat_Compare(AMat, BMat, XMaxm, XMean)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  zMat_Compare(AMat, BMat, XMaxm, XMean)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compare two (NumNS, NumNS) complex different matrices and calculate their
!             differences in elements.
! KEYWORDS: Compare two matrices, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Compare two matrices, complex version. 
!
!     Input: AMat --> Input matrix A(N, M);
!            BMat --> Input matrix A(N, M);
!     
!     Outpt: XMaxm --> The largest absolute difference;
!            XMean --> The average difference in all elements.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		complex(rp) XMaxm               ! Maximum difference for matrix elements in A, B matrices
		complex(rp) XMean               ! Average difference for matrix elements in A, B matrices
		complex(rp) AMat(NumNS, NumNS)  ! Input A matrix
		complex(rp) BMat(NumNS, NumNS)  ! Input B matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1      ! Loop integer 
		integer I2      ! Loop integer
		real(rp) Diff   ! Difference between the elements in A and B
      
      real(rp), allocatable :: XMaxmVec(:)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Creating matrix _________________________________________
!______________________________________________________________________________________________________________
		XMaxm = 0.d0
		XMean = 0.d0
      
      allocate(XMaxmVec(NumNS))
      XMaxmVec = 0.d0
      
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2, Diff)
   !$OMP DO REDUCTION(+ : XMean)
		do I1 = 1, NumNS, +1
			do I2 = 1, NumNS, +1
            Diff = dsqrt( dreal( (BMat(I1, I2) - AMat(I1, I2)) * dconjg(BMat(I1, I2) - AMat(I1, I2)) ) )
				if(Diff .gt. XMaxmVec(I1)) then
					XMaxmVec(I1) = Diff
				end if
				XMean = XMean + Diff
			enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
				
		XMean = XMean / dble(NumNS * NumNS)
      XMaxm = maxval(XMaxmVec)
      
      deallocate(XMaxmVec)
		
	end subroutine zMat_Compare
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$