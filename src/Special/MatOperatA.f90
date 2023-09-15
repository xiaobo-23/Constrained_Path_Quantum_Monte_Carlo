!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform matrix operations including multiplication, inverse, determinant and 
!               UDV decomposition for QMC simulations, based on spin-coupled and spin-decoupled cases. 
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!    All the following subroutines call the internal subroutines of Intel MKL package, and they are designed to 
!         perform matrix operations according to whether the spin-up and spin-down are decoupled or not for the
!         simulated model.
!
!    zMatDet_LogDet_QMC --> Subroutine to calculate the log of determinant of complex square matrix; 
!    dMatDet_LogDet_QMC --> Subroutine to calculate the log of determinant of real    square matrix; 
!
!    zMatInv_LogDet_QMC --> Subroutine to calculate the inverse and log(determinant) of complex square matrix;
!    dMatInv_LogDet_QMC --> Subroutine to calculate the inverse and log(determinant) of real    square matrix;
!
!    zMatEqSet_Left_LogDet_QMC --> Subroutine to calculate XMat*AMat=BMat (or BMat*AMat^{-1}), complex version; 
!    zMatEqSet_Rght_LogDet_QMC --> Subroutine to calculate AMat*XMat=BMat (or AMat^{-1}*BMat), complex version; 
!    dMatEqSet_Left_LogDet_QMC --> Subroutine to calculate XMat*AMat=BMat (or BMat*AMat^{-1}), real    version; 
!    dMatEqSet_Rght_LogDet_QMC --> Subroutine to calculate AMat*XMat=BMat (or AMat^{-1}*BMat), real    version; 
!  
!    zMatUDV_QMC --> Subroutine to calculate the UDV decomposition of complex matrix;  
!    dMatUDV_QMC --> Subroutine to calculate the UDV decomposition of real    matrix; 
!
!    zMat_QR_QMC --> Subroutine to calculate the QR decomposition of complex matrix;  
!    dMat_QR_QMC --> Subroutine to calculate the QR decomposition of real    matrix; 
!
!    Matrix product of complex version:
!      zMatPrd_NN_QMC --> Subroutine to calculate matrix product  A(ND1, ND2)   * B(ND2, ND3)   , according to H; 
!      zMatPrd_CN_QMC --> Subroutine to calculate matrix product [A(ND2, ND1)]^+* B(ND2, ND3)   , according to H; 
!      zMatPrd_NC_QMC --> Subroutine to calculate matrix product  A(ND1, ND2)   *[B(ND3, ND2)]^+, according to H; 
!      zMatPrd_CC_QMC --> Subroutine to calculate matrix product [A(ND2, ND1)]^+*[B(ND3, ND2)]^+, according to H; 
!      zMatPrd_TT_QMC --> Subroutine to calculate matrix product [A(ND2, ND1)]^T*[B(ND3, ND2)]^T, according to H; 
!
!    Matrix product of real version:
!      dMatPrd_NN_QMC --> Subroutine to calculate matrix product  A(ND1, ND2)   * B(ND2, ND3)   , according to H; 
!      dMatPrd_TN_QMC --> Subroutine to calculate matrix product [A(ND2, ND1)]^T* B(ND2, ND3)   , according to H; 
!      dMatPrd_NT_QMC --> Subroutine to calculate matrix product  A(ND1, ND2)   *[B(ND3, ND2)]^T, according to H; 
!      dMatPrd_TT_QMC --> Subroutine to calculate matrix product [A(ND2, ND1)]^T*[B(ND3, ND2)]^T, according to H;  
!
!    zMat_Diag_QMC --> Subroutine used to diagonalize complex Hermitian matrix;
!    dMat_Diag_QMC --> Subroutine used to diagonalize real symmetric    matrix;
!
!    zMat_Compare_QMC --> Subroutine used to compare two complex (NumNS, NumNS, 2) matrices;
!    dMat_Compare_QMC --> Subroutine used to compare two real    (NumNS, NumNS, 2) matrices;
!
!    zMat_Copy_QMC --> Subroutine used to copy complex (NumNS, NumNS, 2) matrix;
!    dMat_Copy_QMC --> Subroutine used to copy real    (NumNS, NumNS, 2) matrix.
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
	subroutine zMatDet_LogDet_QMC(NDim, zMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatDet_LogDet_QMC(NDim, zMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate log determinant for complex square matrix, based on 
!                  sign problem or not cases.
! KEYWORDS: Calculate Matrix determinant for complex matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix determinant, complex version. Based on sign problem or not cases of the system.
!
!     Input: NDim --> Dimension of input A matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: LogzDet --> Log of Complex determinant.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
#ifdef DCPDNOSIGN
      complex(rp) zMat(NDim, NDim, 1)      ! The input zMat matrix
#else
		complex(rp) zMat(NDim, NDim, 2)      ! The input zMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
      
      integer Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetUp, LogzDetDw
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix determinant _____________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&   
		TimsMtDet = TimsMtDet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      LogzDet = rp_Zzero
      call MatDetermtZ_LogDet(NDim, zMat(1, 1, 1), NDim, LogzDet)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
      LogzDet = rp_Zzero; LogzDetUp = rp_Zzero; LogzDetDw = rp_Zzero
      call MatDetermtZ_LogDet(NDim, zMat(1, 1, 1), NDim, LogzDetUp)
      call MatDetermtZ_LogDet(NDim, zMat(1, 1, 2), NDim, LogzDetDw)
      LogzDet = LogzDetUp + LogzDetDw
#endif
!**************************************************************************************************	  
!___________________ 2. Reset the imaginary part of LogzDet _______________________________________
!**************************************************************************************************      
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0_rp/rp_pi)
      ImLogzDet = ImLogzDet - 2.0_rp*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix determinant _____________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&    
		call system_clock(time2)
		TimeMtDet = TimeMtDet + TimeIntrvl(time1, time2)
		
   end subroutine zMatDet_LogDet_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatDet_LogDet_QMC(NDim, dMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatDet_LogDet_QMC(NDim, dMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate log determinant for real square matrix, based on 
!                  sign problem or not cases.
! KEYWORDS: Calculate Matrix determinant for real matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix determinant, real version. Based on sign problem or not cases of the system.
!
!     Input: NDim --> Dimension of input A matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: LogzDet --> Log of Complex determinant.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
#ifdef DCPDNOSIGN
      real(rp) dMat(NDim, NDim, 1)         ! The input dMat matrix
#else
		real(rp) dMat(NDim, NDim, 2)         ! The input dMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
      
      integer Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetUp, LogzDetDw
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix determinant _____________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&   
		TimsMtDet = TimsMtDet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      LogzDet = rp_Zzero
      call MatDetermtR_LogDet(NDim, dMat(1, 1, 1), NDim, LogzDet)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
      LogzDet = rp_Zzero; LogzDetUp = rp_Zzero; LogzDetDw = rp_Zzero
      call MatDetermtR_LogDet(NDim, dMat(1, 1, 1), NDim, LogzDetUp)
      call MatDetermtR_LogDet(NDim, dMat(1, 1, 2), NDim, LogzDetDw)
      LogzDet = LogzDetUp + LogzDetDw
#endif
!**************************************************************************************************	  
!___________________ 2. Reset the imaginary part of LogzDet _______________________________________
!**************************************************************************************************      
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0_rp/rp_pi)
      ImLogzDet = ImLogzDet - 2.0_rp*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix determinant _____________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&    
		call system_clock(time2)
		TimeMtDet = TimeMtDet + TimeIntrvl(time1, time2)
		
   end subroutine dMatDet_LogDet_QMC
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
	subroutine zMatInv_LogDet_QMC(NDim, zMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatInv_LogDet_QMC(NDim, zMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate inverse and log determinant for complex square matrix, based on 
!                  sign problem or not cases.
! KEYWORDS: Calculate Matrix inverse and log determinant for complex matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix Inverse and log determinant, complex version. Based on sign problem or not cases 
!        of the system.
!
!     Input: NDim --> Dimension of input A matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: LogzDet --> Log of Complex determinant.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
#ifdef DCPDNOSIGN
      complex(rp) zMat(NDim, NDim, 1)      ! The input zMat matrix
#else
		complex(rp) zMat(NDim, NDim, 2)      ! The input zMat matrix
#endif
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
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      LogzDet = rp_Zzero
      call MatInversZ_LogDet(NDim, zMat(1, 1, 1), NDim, LogzDet)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
      LogzDet = rp_Zzero; LogzDetUp = rp_Zzero; LogzDetDw = rp_Zzero
      call MatInversZ_LogDet(NDim, zMat(1, 1, 1), NDim, LogzDetUp)
      call MatInversZ_LogDet(NDim, zMat(1, 1, 2), NDim, LogzDetDw)
      LogzDet = LogzDetUp + LogzDetDw
#endif
!**************************************************************************************************	  
!___________________ 2. Reset the imaginary part of LogzDet _______________________________________
!**************************************************************************************************      
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0_rp/rp_pi)
      ImLogzDet = ImLogzDet - 2.0_rp*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix inverse _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%& 
		call system_clock(time2)
		TimeMtInv = TimeMtInv + TimeIntrvl(time1, time2)
		
   end subroutine zMatInv_LogDet_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatInv_LogDet_QMC(NDim, dMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatInv_LogDet_QMC(NDim, dMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate inverse and log determinant for real square matrix, based on 
!                  sign problem or not cases.
! KEYWORDS: Calculate Matrix inverse and log determinant for real matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix Inverse and log determinant, real version. Based on sign problem or not cases 
!           of the system.
!
!     Input: NDim --> Dimension of input A matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: LogzDet --> Log of Complex determinant.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
		integer NDim                        ! Dimension of dMat square matrix
      complex(rp) LogzDet                 ! Log Determinant of dMat matrix
#ifdef DCPDNOSIGN
      real(rp) dMat(NDim, NDim, 1)        ! The input dMat matrix
#else
		real(rp) dMat(NDim, NDim, 2)        ! The input dMat matrix
#endif
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
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      LogzDet = rp_Zzero
      call MatInversR_LogDet(NDim, dMat(1, 1, 1), NDim, LogzDet)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
      LogzDet = rp_Zzero; LogzDetUp = rp_Zzero; LogzDetDw = rp_Zzero
      call MatInversR_LogDet(NDim, dMat(1, 1, 1), NDim, LogzDetUp)
      call MatInversR_LogDet(NDim, dMat(1, 1, 2), NDim, LogzDetDw)
      LogzDet = LogzDetUp + LogzDetDw
#endif
!**************************************************************************************************	  
!___________________ 2. Reset the imaginary part of LogzDet _______________________________________
!**************************************************************************************************      
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0_rp/rp_pi)
      ImLogzDet = ImLogzDet - 2.0_rp*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix inverse _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%& 
		call system_clock(time2)
		TimeMtInv = TimeMtInv + TimeIntrvl(time1, time2)
		
   end subroutine dMatInv_LogDet_QMC
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
	subroutine zMatEqSet_Left_LogDet_QMC(ND1, ND2, AMat, BMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatEqSet_Left_LogDet_QMC(ND1, ND2, AMat, BMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve XMat(ND1, ND2)*AMat(ND2, ND2)=BMat(ND1, ND2) or calculate BMat*AMat^{-1}.
! KEYWORDS: Solve XMat*AMat=BMat or Calculate BMat*AMat^{-1}, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat*AMat=BMat or Calculate BMat*AMat^{-1}, by the method of performing LU decomposition and solving 
!             linear equation set.
!
!     Input: ND1   --> Dimension of input matrices AMat(ND2, ND2) and BMat(ND1, ND2);
!            ND2   --> Dimension of input matrices AMat(ND2, ND2) and BMat(ND1, ND2);
!            AMat  --> Input AMat(ND2, ND2) matrix;
!            BMat  --> Input BMat(ND1, ND2) matrix;
!
!     Outpt: BMat    --> Output solution XMat matrix or result of BMat*AMat^{-1};
!            LogzDet --> Log of Complex determinant of AMat matrix.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
#ifdef DCPDNOSIGN
		complex(rp) AMat(ND2, ND2, 1) 
      complex(rp) BMat(ND1, ND2, 1) 
#else
		complex(rp) AMat(ND2, ND2, 2) 
      complex(rp) BMat(ND1, ND2, 2) 
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
      
      integer Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetUp, LogzDetDw
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsEqSet = TimsEqSet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      LogzDet = rp_Zzero
      call SlvLnEqSetZ_Left_LogDet(ND1, ND2, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND1, LogzDet)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
      LogzDet = rp_Zzero; LogzDetUp = rp_Zzero; LogzDetDw = rp_Zzero
      call SlvLnEqSetZ_Left_LogDet(ND1, ND2, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND1, LogzDet)
      call SlvLnEqSetZ_Left_LogDet(ND1, ND2, AMat(1, 1, 2), ND2, BMat(1, 1, 2), ND1, LogzDet)
#endif
!**************************************************************************************************	  
!___________________ 2. Reset the imaginary part of LogzDet _______________________________________
!**************************************************************************************************      
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0_rp/rp_pi)
      ImLogzDet = ImLogzDet - 2.0_rp*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeEqSet = TimeEqSet + TimeIntrvl(time1, time2)
		
   end subroutine zMatEqSet_Left_LogDet_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatEqSet_Rght_LogDet_QMC(ND1, ND2, AMat, BMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatEqSet_Rght_LogDet_QMC(ND1, ND2, AMat, BMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve AMat(ND1, ND1)*XMat(ND1, ND2)=BMat(ND1, ND2) or calculate AMat^{-1}*BMat.
! KEYWORDS: Solve AMat*XMat=BMat or Calculate AMat^{-1}*BMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat*XMat=BMat or Calculate AMat^{-1}*BMat, by the method of performing LU decomposition and solving 
!               linear equation set.
!
!     Input: ND1   --> Dimension of input matrices AMat(ND1, ND1) and BMat(ND1, ND2);
!            ND2   --> Dimension of input matrices AMat(ND1, ND1) and BMat(ND1, ND2);
!            AMat  --> Input AMat(ND1, ND1) matrix;
!            BMat  --> Input BMat(ND1, ND2) matrix;
!
!     Outpt: BMat    --> Output solution XMat matrix or result of AMat^{-1}*BMat;
!            LogzDet --> Log of Complex determinant of AMat matrix.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
#ifdef DCPDNOSIGN
		complex(rp) AMat(ND1, ND1, 1) 
      complex(rp) BMat(ND1, ND2, 1) 
#else
		complex(rp) AMat(ND1, ND1, 2) 
      complex(rp) BMat(ND1, ND2, 2) 
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
      
      integer Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetUp, LogzDetDw
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsEqSet = TimsEqSet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      LogzDet = rp_Zzero
      call SlvLnEqSetZ_Rght_LogDet(ND1, ND2, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND1, LogzDet)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
      LogzDet = rp_Zzero; LogzDetUp = rp_Zzero; LogzDetDw = rp_Zzero
      call SlvLnEqSetZ_Rght_LogDet(ND1, ND2, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND1, LogzDetUp)
      call SlvLnEqSetZ_Rght_LogDet(ND1, ND2, AMat(1, 1, 2), ND1, BMat(1, 1, 2), ND1, LogzDetDw)
      LogzDet = LogzDetUp + LogzDetDw
#endif
!**************************************************************************************************	  
!___________________ 2. Reset the imaginary part of LogzDet _______________________________________
!**************************************************************************************************      
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0_rp/rp_pi)
      ImLogzDet = ImLogzDet - 2.0_rp*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeEqSet = TimeEqSet + TimeIntrvl(time1, time2)
		
   end subroutine zMatEqSet_Rght_LogDet_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatEqSet_Left_LogDet_QMC(ND1, ND2, AMat, BMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatEqSet_Left_LogDet_QMC(ND1, ND2, AMat, BMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve XMat(ND1, ND2)*AMat(ND2, ND2)=BMat(ND1, ND2) or calculate BMat*AMat^{-1}.
! KEYWORDS: Solve XMat*AMat=BMat or Calculate BMat*AMat^{-1}, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat*AMat=BMat or Calculate BMat*AMat^{-1}, by the method of performing LU decomposition and solving 
!            linear equation set.
!
!     Input: ND1   --> Dimension of input matrices AMat(ND2, ND2) and BMat(ND1, ND2);
!            ND2   --> Dimension of input matrices AMat(ND2, ND2) and BMat(ND1, ND2);
!            AMat  --> Input AMat(ND2, ND2) matrix;
!            BMat  --> Input BMat(ND1, ND2) matrix;
!
!     Outpt: BMat    --> Output solution XMat matrix or result of BMat*AMat^{-1};
!            LogzDet --> Log of real determinant of AMat matrix.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
#ifdef DCPDNOSIGN
		real(rp) AMat(ND2, ND2, 1) 
      real(rp) BMat(ND1, ND2, 1) 
#else
		real(rp) AMat(ND2, ND2, 2) 
      real(rp) BMat(ND1, ND2, 2) 
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
      integer Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetUp, LogzDetDw
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsEqSet = TimsEqSet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      LogzDet = rp_Zzero
      call SlvLnEqSetR_Left_LogDet(ND1, ND2, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND1, LogzDet)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
      LogzDet = rp_Zzero; LogzDetUp = rp_Zzero; LogzDetDw = rp_Zzero
      call SlvLnEqSetR_Left_LogDet(ND1, ND2, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND1, LogzDetUp)
      call SlvLnEqSetR_Left_LogDet(ND1, ND2, AMat(1, 1, 2), ND2, BMat(1, 1, 2), ND1, LogzDetDw)
      LogzDet = LogzDetUp + LogzDetDw
#endif
!**************************************************************************************************	  
!___________________ 2. Reset the imaginary part of LogzDet _______________________________________
!**************************************************************************************************      
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0_rp/rp_pi)
      ImLogzDet = ImLogzDet - 2.0_rp*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeEqSet = TimeEqSet + TimeIntrvl(time1, time2)
		
   end subroutine dMatEqSet_Left_LogDet_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatEqSet_Rght_LogDet_QMC(ND1, ND2, AMat, BMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatEqSet_Rght_LogDet_QMC(ND1, ND2, AMat, BMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve AMat(ND1, ND1)*XMat(ND1, ND2)=BMat(ND1, ND2) or calculate AMat^{-1}*BMat.
! KEYWORDS: Solve AMat*XMat=BMat or Calculate AMat^{-1}*BMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat*XMat=BMat or Calculate AMat^{-1}*BMat, by the method of performing LU decomposition and solving 
!               linear equation set.
!
!     Input: ND1   --> Dimension of input matrices AMat(ND1, ND1) and BMat(ND1, ND2);
!            ND2   --> Dimension of input matrices AMat(ND1, ND1) and BMat(ND1, ND2);
!            AMat  --> Input AMat(ND1, ND1) matrix;
!            BMat  --> Input BMat(ND1, ND2) matrix;
!
!     Outpt: BMat    --> Output solution XMat matrix or result of AMat^{-1}*BMat;
!            LogzDet --> Log of real determinant of AMat matrix.     
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
#ifdef DCPDNOSIGN
		real(rp) AMat(ND1, ND1, 1) 
      real(rp) BMat(ND1, ND2, 1)
#else
		real(rp) AMat(ND1, ND1, 2) 
      real(rp) BMat(ND1, ND2, 2) 
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                    ! Starting time point
		integer(8) time2                    ! Ending   time point
      
      integer Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetUp, LogzDetDw
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&      
		TimsEqSet = TimsEqSet + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of log determinant of matrix _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      LogzDet = rp_Zzero
      call SlvLnEqSetR_Rght_LogDet(ND1, ND2, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND1, LogzDet)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
      LogzDet = rp_Zzero; LogzDetUp = rp_Zzero; LogzDetDw = rp_Zzero
      call SlvLnEqSetR_Rght_LogDet(ND1, ND2, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND1, LogzDetUp)
      call SlvLnEqSetR_Rght_LogDet(ND1, ND2, AMat(1, 1, 2), ND1, BMat(1, 1, 2), ND1, LogzDetDw)
      LogzDet = LogzDetUp + LogzDetDw
#endif
!**************************************************************************************************	  
!___________________ 2. Reset the imaginary part of LogzDet _______________________________________
!**************************************************************************************************      
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0_rp/rp_pi)
      ImLogzDet = ImLogzDet - 2.0_rp*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in solving linear equation set ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeEqSet = TimeEqSet + TimeIntrvl(time1, time2)
		
   end subroutine dMatEqSet_Rght_LogDet_QMC
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
	subroutine zMatUDV_QMC(ND1, ND2, zMat, UMat, DVec, VMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatUDV_QMC(ND1, ND2, zMat, UMat, DVec, VMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate perform the UDV decomposition for zMat(2*NumNS, NDim) matrix, 
!                  based on sign problem or not cases.
! KEYWORDS: Calculate UDV decomposition for complex matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate UDV decomposition for complex matrix. Based on sign problem or not cases of the system.
!     This only works for ND1 >= ND2 cases.
!
!     Input: ND1  --> Dimension of zMat(ND1, ND2) matrix;
!            ND2  --> Dimension of zMat(ND1, ND2) matrix;
!            zMat --> Input complex matrix.
!
!     Outpt: UMat --> Result U matrix;
!            DVec --> Result D vector;
!            VMat --> Result V matrix.
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
      integer ND1, ND2
      complex(rp) zMat(ND1, ND2, NmSpn)      ! The input dMat matrix
      complex(rp) UMat(ND1, ND2, NmSpn)      ! The Output UMat matrix
      real(rp)    DVec(ND2,      NmSpn)      ! The Output D    vector
      complex(rp) VMat(ND2, ND2, NmSpn)      ! The Output VMat matrix
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
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
!________________________________________________________________________________________ 	  
!_________________ (0) zMat = UMat * DVec * VMat by SVD or QR methods ___________________
!________________________________________________________________________________________
#ifdef QRSTABLIZE
      call UDVdcmpsZ_QRPivot(ND1, ND2, zMat(1, 1, 1), ND1, UMat(1, 1, 1), ND1, DVec(1, 1), &
         & VMat(1, 1, 1), ND2, IfCheckUDV)
#else
      call UDVdcmpsZ_SVD    (ND1, ND2, zMat(1, 1, 1), ND1, UMat(1, 1, 1), ND1, DVec(1, 1), &
         & VMat(1, 1, 1), ND2, IfCheckUDV)
#endif
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
!________________________________________________________________________________________ 	  
!_________________ (0) zMat = UMat * DVec * VMat by SVD or QR methods ___________________
!________________________________________________________________________________________
#ifdef QRSTABLIZE
      call UDVdcmpsZ_QRPivot(ND1, ND2, zMat(1, 1, 1), ND1, UMat(1, 1, 1), ND1, DVec(1, 1), &
         & VMat(1, 1, 1), ND2, IfCheckUDV)
      call UDVdcmpsZ_QRPivot(ND1, ND2, zMat(1, 1, 2), ND1, UMat(1, 1, 2), ND1, DVec(1, 2), &
         & VMat(1, 1, 2), ND2, IfCheckUDV)
#else
      call UDVdcmpsZ_SVD(ND1, ND2, zMat(1, 1, 1), ND1, UMat(1, 1, 1), ND1, DVec(1, 1), &
         & VMat(1, 1, 1), ND2, IfCheckUDV)
      call UDVdcmpsZ_SVD(ND1, ND2, zMat(1, 1, 2), ND1, UMat(1, 1, 2), ND1, DVec(1, 2), &
         & VMat(1, 1, 2), ND2, IfCheckUDV)
#endif

#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in calculating UDV decomposition __________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeUDVOt = TimeUDVOt + TimeIntrvl(time1, time2)
		
   end subroutine zMatUDV_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatUDV_QMC(ND1, ND2, dMat, UMat, DVec, VMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatUDV_QMC(ND1, ND2, dMat, UMat, DVec, VMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate perform the UDV decomposition for dMat(2*NumNS, NDim) matrix, 
!                  based on sign problem or not cases.
! KEYWORDS: Calculate UDV decomposition for real matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate UDV decomposition for real matrix. Based on sign problem or not cases of the system.
!     This only works for ND1 >= ND2 cases.
!
!     Input: ND1  --> Dimension of zMat(ND1, ND2) matrix;
!            ND2  --> Dimension of zMat(ND1, ND2) matrix;
!            dMat --> Input real matrix.
!
!     Outpt: UMat --> Result U matrix;
!            DVec --> Result D vector;
!            VMat --> Result V matrix.
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
      integer ND1, ND2
      real(rp) dMat(ND1, ND2, NmSpn)      ! The input dMat matrix
      real(rp) UMat(ND1, ND2, NmSpn)      ! The Output UMat matrix
      real(rp) DVec(ND2,      NmSpn)      ! The Output D    vector
      real(rp) VMat(ND2, ND2, NmSpn)      ! The Output VMat matrix
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
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
!________________________________________________________________________________________ 	  
!_________________ (0) zMat = UMat * DVec * VMat by SVD or QR methods ___________________
!________________________________________________________________________________________
#ifdef QRSTABLIZE
      call UDVdcmpsR_QRPivot(ND1, ND2, dMat(1, 1, 1), ND1, UMat(1, 1, 1), ND1, DVec(1, 1), &
         & VMat(1, 1, 1), ND2, IfCheckUDV)
#else
      call UDVdcmpsR_SVD    (ND1, ND2, dMat(1, 1, 1), ND1, UMat(1, 1, 1), ND1, DVec(1, 1), &
         & VMat(1, 1, 1), ND2, IfCheckUDV)
#endif
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
!________________________________________________________________________________________ 	  
!_________________ (0) zMat = UMat * DVec * VMat by SVD or QR methods ___________________
!________________________________________________________________________________________
#ifdef QRSTABLIZE
      call UDVdcmpsR_QRPivot(ND1, ND2, dMat(1, 1, 1), ND1, UMat(1, 1, 1), ND1, DVec(1, 1), &
         & VMat(1, 1, 1), ND2, IfCheckUDV)
      call UDVdcmpsR_QRPivot(ND1, ND2, dMat(1, 1, 2), ND1, UMat(1, 1, 2), ND1, DVec(1, 2), &
         & VMat(1, 1, 2), ND2, IfCheckUDV)
#else
      call UDVdcmpsR_SVD(ND1, ND2, dMat(1, 1, 1), ND1, UMat(1, 1, 1), ND1, DVec(1, 1), &
         & VMat(1, 1, 1), ND2, IfCheckUDV)
      call UDVdcmpsR_SVD(ND1, ND2, dMat(1, 1, 2), ND1, UMat(1, 1, 2), ND1, DVec(1, 2), &
         & VMat(1, 1, 2), ND2, IfCheckUDV)
#endif

#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in calculating UDV decomposition __________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
		call system_clock(time2)
		TimeUDVOt = TimeUDVOt + TimeIntrvl(time1, time2)
		
   end subroutine dMatUDV_QMC
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
	subroutine zMatPrd_NN_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatPrd_NN_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of AMat(ND1, ND2, 2) * BMat(ND2, ND3, 2).
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, complex version.
!
!     CMat(ND1, ND3) = Ztp1*AMat(ND1, ND2)*BMat(ND2, ND3) + Ztp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      complex(rp) Ztp1, Ztp2
#ifdef DCPDNOSIGN
		complex(rp) AMat(ND1, ND2, 1)     ! Input AMat matrix
		complex(rp) BMat(ND2, ND3, 1)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		complex(rp) AMat(ND1, ND2, 2)     ! Input AMat matrix
		complex(rp) BMat(ND2, ND3, 2)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call ZGEMM("N", "N", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND2, Ztp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call ZGEMM("N", "N", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND2, Ztp2, CMat(1, 1, 1), ND1)
      call ZGEMM("N", "N", ND1, ND3, ND2, Ztp1, AMat(1, 1, 2), ND1, BMat(1, 1, 2), ND2, Ztp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine zMatPrd_NN_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatPrd_CN_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatPrd_CN_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of [AMat(ND2, ND1, 2)]^+ * BMat(ND2, ND3, 2).
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, complex version.
!
!     CMat(ND1, ND3) = Ztp1*[AMat(ND2, ND1)]^+*BMat(ND2, ND3) + Ztp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      complex(rp) Ztp1, Ztp2
#ifdef DCPDNOSIGN
		complex(rp) AMat(ND2, ND1, 1)     ! Input AMat matrix
		complex(rp) BMat(ND2, ND3, 1)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		complex(rp) AMat(ND2, ND1, 2)     ! Input AMat matrix
		complex(rp) BMat(ND2, ND3, 2)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call ZGEMM("C", "N", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND2, Ztp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call ZGEMM("C", "N", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND2, Ztp2, CMat(1, 1, 1), ND1)
      call ZGEMM("C", "N", ND1, ND3, ND2, Ztp1, AMat(1, 1, 2), ND2, BMat(1, 1, 2), ND2, Ztp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine zMatPrd_CN_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatPrd_NC_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatPrd_NC_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of AMat(ND1, ND2, 2) * [BMat(ND3, ND2, 2)]^+.
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, complex version.
!
!     CMat(ND1, ND3) = Ztp1*AMat(ND1, ND2)*[BMat(ND3, ND2)]^+ + Ztp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      complex(rp) Ztp1, Ztp2
#ifdef DCPDNOSIGN
		complex(rp) AMat(ND1, ND2, 1)     ! Input AMat matrix
		complex(rp) BMat(ND3, ND2, 1)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		complex(rp) AMat(ND1, ND2, 2)     ! Input AMat matrix
		complex(rp) BMat(ND3, ND2, 2)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call ZGEMM("N", "C", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND3, Ztp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call ZGEMM("N", "C", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND3, Ztp2, CMat(1, 1, 1), ND1)
      call ZGEMM("N", "C", ND1, ND3, ND2, Ztp1, AMat(1, 1, 2), ND1, BMat(1, 1, 2), ND3, Ztp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine zMatPrd_NC_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatPrd_CC_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatPrd_CC_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of [AMat(ND2, ND1, 2)]^+ * [BMat(ND3, ND2, 2)]^+.
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, complex version.
!
!     CMat(ND1, ND3) = Ztp1*[AMat(ND2, ND1)]^+*[BMat(ND3, ND2)]^+  +  Ztp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      complex(rp) Ztp1, Ztp2
#ifdef DCPDNOSIGN
		complex(rp) AMat(ND2, ND1, 1)     ! Input AMat matrix
		complex(rp) BMat(ND3, ND2, 1)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		complex(rp) AMat(ND2, ND1, 2)     ! Input AMat matrix
		complex(rp) BMat(ND3, ND2, 2)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call ZGEMM("C", "C", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND3, Ztp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call ZGEMM("C", "C", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND3, Ztp2, CMat(1, 1, 1), ND1)
      call ZGEMM("C", "C", ND1, ND3, ND2, Ztp1, AMat(1, 1, 2), ND2, BMat(1, 1, 2), ND3, Ztp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine zMatPrd_CC_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMatPrd_TT_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMatPrd_TT_QMC(ND1, ND2, ND3, Ztp1, AMat, BMat, Ztp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of [AMat(ND2, ND1, 2)]^T * [BMat(ND3, ND2, 2)]^T.
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, complex version.
!
!     CMat(ND1, ND3) = Ztp1*[AMat(ND2, ND1)]^T*[BMat(ND3, ND2)]^T  +  Ztp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      complex(rp) Ztp1, Ztp2
#ifdef DCPDNOSIGN
		complex(rp) AMat(ND2, ND1, 1)     ! Input AMat matrix
		complex(rp) BMat(ND3, ND2, 1)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		complex(rp) AMat(ND2, ND1, 2)     ! Input AMat matrix
		complex(rp) BMat(ND3, ND2, 2)     ! Input BMat matrix
		complex(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call ZGEMM("T", "T", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND3, Ztp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call ZGEMM("T", "T", ND1, ND3, ND2, Ztp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND3, Ztp2, CMat(1, 1, 1), ND1)
      call ZGEMM("T", "T", ND1, ND3, ND2, Ztp1, AMat(1, 1, 2), ND2, BMat(1, 1, 2), ND3, Ztp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine zMatPrd_TT_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatPrd_NN_QMC(ND1, ND2, ND3, Rtp1, AMat, BMat, Rtp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatPrd_NN_QMC(ND1, ND2, ND3, Rtp1, AMat, BMat, Rtp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of AMat(ND1, ND2, 2) * BMat(ND2, ND3, 2).
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, real version.
!
!     CMat(ND1, ND3) = Rtp1*AMat(ND1, ND2)*BMat(ND2, ND3) + Rtp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      real(rp) Rtp1, Rtp2
#ifdef DCPDNOSIGN
		real(rp) AMat(ND1, ND2, 1)     ! Input AMat matrix
		real(rp) BMat(ND2, ND3, 1)     ! Input BMat matrix
		real(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		real(rp) AMat(ND1, ND2, 2)     ! Input AMat matrix
		real(rp) BMat(ND2, ND3, 2)     ! Input BMat matrix
		real(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call DGEMM("N", "N", ND1, ND3, ND2, Rtp1, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND2, Rtp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call DGEMM("N", "N", ND1, ND3, ND2, Rtp1, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND2, Rtp2, CMat(1, 1, 1), ND1)
      call DGEMM("N", "N", ND1, ND3, ND2, Rtp1, AMat(1, 1, 2), ND1, BMat(1, 1, 2), ND2, Rtp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine dMatPrd_NN_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatPrd_TN_QMC(ND1, ND2, ND3, Rtp1, AMat, BMat, Rtp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatPrd_TN_QMC(ND1, ND2, ND3, Rtp1, AMat, BMat, Rtp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of [AMat(ND2, ND1, 2)]^T * BMat(ND2, ND3, 2).
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, real version.
!
!     CMat(ND1, ND3) = Rtp1*[AMat(ND2, ND1)]^T*BMat(ND2, ND3) + Rtp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      real(rp) Rtp1, Rtp2
#ifdef DCPDNOSIGN
		real(rp) AMat(ND2, ND1, 1)     ! Input AMat matrix
		real(rp) BMat(ND2, ND3, 1)     ! Input BMat matrix
		real(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		real(rp) AMat(ND2, ND1, 2)     ! Input AMat matrix
		real(rp) BMat(ND2, ND3, 2)     ! Input BMat matrix
		real(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call DGEMM("T", "N", ND1, ND3, ND2, Rtp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND2, Rtp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call DGEMM("T", "N", ND1, ND3, ND2, Rtp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND2, Rtp2, CMat(1, 1, 1), ND1)
      call DGEMM("T", "N", ND1, ND3, ND2, Rtp1, AMat(1, 1, 2), ND2, BMat(1, 1, 2), ND2, Rtp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine dMatPrd_TN_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatPrd_NT_QMC(ND1, ND2, ND3, Rtp1, AMat, BMat, Rtp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatPrd_NT_QMC(ND1, ND2, ND3, Rtp1, AMat, BMat, Rtp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of AMat(ND1, ND2, 2) * [BMat(ND3, ND2, 2)]^T.
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, real version.
!
!     CMat(ND1, ND3) = Rtp1*AMat(ND1, ND2)*[BMat(ND3, ND2)]^T + Rtp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      real(rp) Rtp1, Rtp2
#ifdef DCPDNOSIGN
		real(rp) AMat(ND1, ND2, 1)     ! Input AMat matrix
		real(rp) BMat(ND3, ND2, 1)     ! Input BMat matrix
		real(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		real(rp) AMat(ND1, ND2, 2)     ! Input AMat matrix
		real(rp) BMat(ND3, ND2, 2)     ! Input BMat matrix
		real(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call DGEMM("N", "T", ND1, ND3, ND2, Rtp1, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND3, Rtp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call DGEMM("N", "T", ND1, ND3, ND2, Rtp1, AMat(1, 1, 1), ND1, BMat(1, 1, 1), ND3, Rtp2, CMat(1, 1, 1), ND1)
      call DGEMM("N", "T", ND1, ND3, ND2, Rtp1, AMat(1, 1, 2), ND1, BMat(1, 1, 2), ND3, Rtp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine dMatPrd_NT_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMatPrd_TT_QMC(ND1, ND2, ND3, Rtp1, AMat, BMat, Rtp2, CMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMatPrd_TT_QMC(ND1, ND2, ND3, Rtp1, AMat, BMat, Rtp2, CMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate matrix product of [AMat(ND1, ND2, 2)]^T * [BMat(ND3, ND2, 2)]^T.
! KEYWORDS: Calculate Matrix product.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix product, real version.
!
!     CMat(ND1, ND3) = Rtp1*[AMat(ND1, ND2)]^T*[BMat(ND3, ND2)]^T + Rtp2*CMat(ND1, ND3).
!
!     Input: 
!
!     Outpt: 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
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
      integer ND1, ND2, ND3
      real(rp) Rtp1, Rtp2
#ifdef DCPDNOSIGN
		real(rp) AMat(ND2, ND1, 1)     ! Input AMat matrix
		real(rp) BMat(ND3, ND2, 1)     ! Input BMat matrix
		real(rp) CMat(ND1, ND3, 1)     ! Ouput CMat matrix
#else
		real(rp) AMat(ND2, ND1, 2)     ! Input AMat matrix
		real(rp) BMat(ND3, ND2, 2)     ! Input BMat matrix
		real(rp) CMat(ND1, ND3, 2)     ! Ouput CMat matrix
#endif
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer(8) time1                  ! Starting time point
		integer(8) time2                  ! Ending   time point
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsMtprd = TimsMtprd + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of Calculating matrix product ________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> Only calculate spin-up channel _______________
!**************************************************************************************************
#ifdef DCPDNOSIGN
      call DGEMM("T", "T", ND1, ND3, ND2, Rtp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND3, Rtp2, CMat(1, 1, 1), ND1) 
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Calculate both up and dn channels ____________
!**************************************************************************************************
#else       
      call DGEMM("T", "T", ND1, ND3, ND2, Rtp1, AMat(1, 1, 1), ND2, BMat(1, 1, 1), ND3, Rtp2, CMat(1, 1, 1), ND1)
      call DGEMM("T", "T", ND1, ND3, ND2, Rtp1, AMat(1, 1, 2), ND2, BMat(1, 1, 2), ND3, Rtp2, CMat(1, 1, 2), ND1)
#endif
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in matrix product _________________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
		TimeMtprd = TimeMtprd + TimeIntrvl(time1, time2)
		
   end subroutine dMatPrd_TT_QMC
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
	subroutine zMat_Diag_QMC(NDim, zMat, EVal, EVec)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMat_Diag_QMC(NDim, zMat, EVal, EVec)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to diagonalize zMat(NDim*NDim) matrix which is Hermitian.
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NDim
#ifdef DCPDNOSIGN
		complex(rp) zMat(NDim, NDim, 1)     ! Input AMat matrix
		real(rp)    EVal(NDim,       1)     ! Ouput Eigenvalues matrix
		complex(rp) EVec(NDim, NDim, 1)     ! Ouput Eigenvectors matrix
#else
		complex(rp) zMat(NDim, NDim, 2)     ! Input AMat matrix
		real(rp)    EVal(NDim,       2)     ! Ouput Eigenvalues matrix
		complex(rp) EVec(NDim, NDim, 2)     ! Ouput Eigenvectors matrix
#endif
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
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> No sign problem ______________________________
!**************************************************************************************************
#ifdef DCPDNOSIGN
!________________________________________________________________________________________ 	  
!_________________ (0) Allocate temperary vectors and initilizations ____________________
!________________________________________________________________________________________
      LWork  = 2 * NDim - 1
      LrWork = 3 * NDim - 2
      allocate( Work(LWork ), stat=istat)
      allocate(rWork(LrWork), stat=istat)
      
      Eval = 0.0d0
      call zcopy(NDim*NDim, zMat(1, 1, 1), 1, EVec(1, 1, 1), 1)
!________________________________________________________________________________________ 	  
!_________________ (1) Call MKL sobroutines to do the diagonalization ___________________
!________________________________________________________________________________________         
      Work  = rp_Zzero
      rWork = 0.0d0
		call ZHEEV("V", "U", NDim, EVec(1, 1, 1), NDim, Eval(1, 1), Work, LWork, rWork, Info)
		if( info /= 0 ) then
         write(*, "('zMat_Diag_QMC: error in lapack subroutine ZHEEV! ierror = ', I4)") info
      end if   
!________________________________________________________________________________________ 	  
!_________________ (2) Deallocate the temperary vectors _________________________________
!________________________________________________________________________________________ 
      if(allocated( Work)) deallocate( Work)
      if(allocated(rWork)) deallocate(rWork)  
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
!________________________________________________________________________________________ 	  
!_________________ (0) Allocate temperary vectors and initilizations ____________________
!________________________________________________________________________________________
      LWork  = 2 * NDim - 1
      LrWork = 3 * NDim - 2
      allocate( Work(LWork ), stat=istat)
      allocate(rWork(LrWork), stat=istat)
      
      Eval = 0.0d0
      call zcopy(NDim*NDim*2, zMat(1, 1, 1), 1, EVec(1, 1, 1), 1)
!________________________________________________________________________________________ 	  
!_________________ (1) Diagonalization for spin-up part _________________________________
!________________________________________________________________________________________      
      Work  = rp_Zzero; rWork = 0.0d0
		call ZHEEV("V", "U", NDim, EVec(1, 1, 1), NDim, Eval(1, 1), Work, LWork, rWork, Info)
		if( info /= 0 ) then
         write(*, "('zMat_Diag_QMC UP: error in lapack subroutine ZHEEV! ierror = ', I4)") info
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) Diagonalization for spin-up part _________________________________
!________________________________________________________________________________________         
      Work  = rp_Zzero; rWork = 0.0d0
		call ZHEEV("V", "U", NDim, EVec(1, 1, 2), NDim, Eval(1, 2), Work, LWork, rWork, Info)
		if( info /= 0 ) then
         write(*, "('zMat_Diag_QMC UP: error in lapack subroutine ZHEEV! ierror = ', I4)") info
      end if
!________________________________________________________________________________________ 	  
!_________________ (3) Deallocate the temperary vectors _________________________________
!________________________________________________________________________________________ 
      if(allocated( Work)) deallocate( Work)
      if(allocated(rWork)) deallocate(rWork)  
#endif
		
   end subroutine zMat_Diag_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMat_Diag_QMC(NDim, dMat, EVal, EVec)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMat_Diag_QMC(NDim, dMat, EVal, EVec)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to diagonalize dMat(NDim*NDim) matrix which is Hermitian.
! KEYWORDS: Diagonalize real symmetric Matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Diagonalize real symmetric Matrix.
!
!     Input:  NDim --> Dimension of the input matrix;
!             dMat --> Input square matrix;
!
!     Output: Eval --> Output eigenvalues;
!             EVec --> Output eigenvector matrix.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NDim
#ifdef DCPDNOSIGN
		real(rp) dMat(NDim, NDim, 1)     ! Input AMat matrix
		real(rp) EVal(NDim,       1)     ! Ouput Eigenvalues matrix
		real(rp) EVec(NDim, NDim, 1)     ! Ouput Eigenvectors matrix
#else
		real(rp) dMat(NDim, NDim, 2)     ! Input AMat matrix
		real(rp) EVal(NDim,       2)     ! Ouput Eigenvalues matrix
		real(rp) EVec(NDim, NDim, 2)     ! Ouput Eigenvectors matrix
#endif
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
!**************************************************************************************************	  
!___________________ 0. For DCPDNOSIGN == T case --> No sign problem ______________________________
!**************************************************************************************************
#ifdef DCPDNOSIGN
!________________________________________________________________________________________ 	  
!_________________ (0) Allocate temperary vectors and initilizations ____________________
!________________________________________________________________________________________
      LWork = 3 * NDim - 1
      allocate(Work(LWork ), stat=istat)
      
      Eval = 0.0d0
      call dcopy(NDim*NDim, dMat(1, 1, 1), 1, EVec(1, 1, 1), 1)
!________________________________________________________________________________________ 	  
!_________________ (1) Call MKL sobroutines to do the diagonalization ___________________
!________________________________________________________________________________________         
      Work = 0.0d0
		call DSYEV("V", "U", NDim, EVec(1, 1, 1), NDim, Eval(1, 1), Work, LWork, Info)
		if( info /= 0 ) then
         write(*, "('dMat_Diag_QMC: error in lapack subroutine DSYEV! ierror = ', I4)") info
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) Deallocate the temperary vectors _________________________________
!________________________________________________________________________________________ 
      if(allocated(Work)) deallocate(Work)
!**************************************************************************************************	  
!___________________ 1. For DCPDNOSIGN == F case --> Sign problem or coupled ______________________
!**************************************************************************************************
#else
!________________________________________________________________________________________ 	  
!_________________ (0) Allocate temperary vectors and initilizations ____________________
!________________________________________________________________________________________
      LWork = 3 * NDim - 1
      allocate(Work(LWork ), stat=istat)
      
      Eval = 0.0d0
      call dcopy(NDim*NDim*2, dMat(1, 1, 1), 1, EVec(1, 1, 1), 1)
!________________________________________________________________________________________ 	  
!_________________ (1) Diagonalization for the spin-up part _____________________________
!________________________________________________________________________________________
      Work = 0.0d0
		call DSYEV("V", "U", NDim, EVec(1, 1, 1), NDim, Eval(1, 1), Work, LWork, Info)
		if( info /= 0 ) then
         write(*, "('dMat_Diag_QMC UP: error in lapack subroutine DSYEV! ierror = ', I4)") info
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) Diagonalization for the spin-down part ___________________________
!________________________________________________________________________________________      
      Work = 0.0d0
		call DSYEV("V", "U", NDim, EVec(1, 1, 2), NDim, Eval(1, 2), Work, LWork, Info)
		if( info /= 0 ) then
         write(*, "('dMat_Diag_QMC DW: error in lapack subroutine DSYEV! ierror = ', I4)") info
      end if
!________________________________________________________________________________________ 	  
!_________________ (3) Deallocate the temperary vectors _________________________________
!________________________________________________________________________________________ 
      if(allocated(Work)) deallocate(Work)
#endif
		
   end subroutine dMat_Diag_QMC
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
	subroutine zMat_Compare_QMC(AMat, BMat, XMaxm, XMean)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMat_Compare_QMC(AMat, BMat, XMaxm, XMean)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compare two complex (NumNS, NumNS, 2) matrix.
! KEYWORDS: Compare complex (NumNS, NumNS, 2) matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Compare complex (NumNS, NumNS, 2) matrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) XMaxm, XMean
      complex(rp) AMat(NumNS, NumNS, NmSpn)
      complex(rp) BMat(NumNS, NumNS, NmSpn)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________       
      real(rp) XMaxmUp, XMeanUp
      real(rp) XMaxmDw, XMeanDw
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of copy the matrix ___________________________________________
!______________________________________________________________________________________________________________ 
#ifdef DCPDNOSIGN
      call MatrCmpr_C(NumNS, NumNS, AMat(1, 1, 1), NumNS, BMat(1, 1, 1), NumNS, XMaxmUp, XMeanUp)
      XMaxm = XMaxmUp
      XMean = XMeanUp
#else
      call MatrCmpr_C(NumNS, NumNS, AMat(1, 1, 1), NumNS, BMat(1, 1, 1), NumNS, XMaxmUp, XMeanUp)
      call MatrCmpr_C(NumNS, NumNS, AMat(1, 1, 2), NumNS, BMat(1, 1, 2), NumNS, XMaxmDw, XMeanDw)
      XMaxm = max(XMaxmUp, XMaxmDw)
      XMean = ( XMeanUp + XMeanDw ) / 2.0_rp
#endif
		
   end subroutine zMat_Compare_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMat_Compare_QMC(AMat, BMat, XMaxm, XMean)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMat_Compare_QMC(AMat, BMat, XMaxm, XMean)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compare two real (NumNS, NumNS, 2) matrix.
! KEYWORDS: Compare real (NumNS, NumNS, 2) matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Compare real (NumNS, NumNS, 2) matrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) XMaxm, XMean
      real(rp) AMat(NumNS, NumNS, NmSpn)
      real(rp) BMat(NumNS, NumNS, NmSpn)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________       
      real(rp) XMaxmUp, XMeanUp
      real(rp) XMaxmDw, XMeanDw
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of copy the matrix ___________________________________________
!______________________________________________________________________________________________________________ 
#ifdef DCPDNOSIGN
      call MatrCmpr_R(NumNS, NumNS, AMat(1, 1, 1), NumNS, BMat(1, 1, 1), NumNS, XMaxmUp, XMeanUp)
      XMaxm = XMaxmUp
      XMean = XMeanUp
#else
      call MatrCmpr_R(NumNS, NumNS, AMat(1, 1, 1), NumNS, BMat(1, 1, 1), NumNS, XMaxmUp, XMeanUp)
      call MatrCmpr_R(NumNS, NumNS, AMat(1, 1, 2), NumNS, BMat(1, 1, 2), NumNS, XMaxmDw, XMeanDw)
      XMaxm = max(XMaxmUp, XMaxmDw)
      XMean = ( XMeanUp + XMeanDw ) / 2.0_rp
#endif
		
   end subroutine dMat_Compare_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
!########################################################################################################################
!########################################################################################################################
!############################################### Copy (NumNS, NumNS, 2) matrix ##########################################
!############################################### Copy (NumNS, NumNS, 2) matrix ##########################################
!############################################### Copy (NumNS, NumNS, 2) matrix ##########################################
!########################################################################################################################
!######################################################################################################################## 
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine zMat_Copy_QMC(AMat, BMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  zMat_Copy_QMC(AMat, BMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to copy complex (NumNS, NumNS, 2) matrix.
! KEYWORDS: Copy complex (NumNS, NumNS, 2) matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Copy complex (NumNS, NumNS, 2) matrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      complex(rp) AMat(NumNS, NumNS, NmSpn)
      complex(rp) BMat(NumNS, NumNS, NmSpn)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of copy the matrix ___________________________________________
!______________________________________________________________________________________________________________  
      call zcopy(NmSpn*NumNS*NumNS, AMat(1, 1, 1), 1, BMat(1, 1, 1), 1)
		
   end subroutine zMat_Copy_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine dMat_Copy_QMC(AMat, BMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  dMat_Copy_QMC(AMat, BMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to copy real (NumNS, NumNS, 2) matrix.
! KEYWORDS: Copy real (NumNS, NumNS, 2) matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Copy real (NumNS, NumNS, 2) matrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use RealPrecsn
      use CoreParamt
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) AMat(NumNS, NumNS, NmSpn)
      real(rp) BMat(NumNS, NumNS, NmSpn)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations of copy the matrix ___________________________________________
!______________________________________________________________________________________________________________ 
      call dcopy(NmSpn*NumNS*NumNS, AMat(1, 1, 1), 1, BMat(1, 1, 1), 1)
		
   end subroutine dMat_Copy_QMC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
