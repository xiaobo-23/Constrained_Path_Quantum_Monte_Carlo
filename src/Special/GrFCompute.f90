!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to calculate the equal-time Green Function G(\beta, \beta), G(0, 0) and G(\tau, \tau)
!              and Dynamic Green Function matrices G(\tau, 0) and G(0, \tau), with both real and complex versions. 
!
!          These subroutines are used when the Numerical Stabalization of DQMC are carried out by QR and SVD algorithm
!              with and without column pivoting.
!
! COMMENT: Calculate the equal-time Green Function.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!   GrFStatcCB    --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, complex version;
!   GrFStatcC0    --> Subroutine to calculate equal-time GreenF matrix at \tau=0 , complex version; 
!   GrFStaticC_LR --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, complex version, UL and UR are unitary;
!   GrFStaticC_UR --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, complex version, UR is unitary;
!   GrFStaticC_UL --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, complex version, UL is unitary;
!   GrFStaticC_NO --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, complex version, UL and UR are general;
!   GrFDynamcC_LR --> Subroutine to calculate Dynamic GreenF matrix G(tau, 0) and G(0, tau), complex version;
!
!   GrFStatcRB    --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, real version;
!   GrFStatcR0    --> Subroutine to calculate equal-time GreenF matrix at \tau=0 , real version; 
!   GrFStaticR_LR --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, real version, UL and UR are unitary;
!   GrFStaticR_UR --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, real version, UR is unitary;
!   GrFStaticR_UL --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, real version, UL is unitary;
!   GrFStaticR_NO --> Subroutine to calculate equal-time GreenF matrix at \tau=NT, real version, UL and UR are general;
!   GrFDynamcR_LR --> Subroutine to calculate Dynamic GreenF matrix G(tau, 0) and G(0, tau), real version;
!
!   SplitDMat_SgleD --> Subroutine used to split the DMat matrix and obtain logdet of DMax matrix;
!   SplitDMat_RtLtD --> Subroutine used to split the DRMat and DLMat matrices and obtain logdet of DLMax*DRMax matrix;
!
!   GrnFSymTrotR --> Subroutine used to obtain the GreenF within symmetric Trotter decomposition,    real version;
!   GrnFSymTrotC --> Subroutine used to obtain the GreenF within symmetric Trotter decomposition, complex version.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
!################################################# For Complex Version ##################################################
!################################################# For Complex Version ##################################################
!################################################# For Complex Version ##################################################
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStatcCB(NT, LogScaleRt, UR, DRVec, VR, LogzDet, GreenF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStatcCB(NT, LogScaleRt, UR, DRVec, VR, LogzDet, GreenF) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, 
!                 VR as G = (1 + UR*DRVec*VR)^{-1}.  Partial single-particle Green's function. 
! KEYWORDS: Calculate static Green's function, Complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     GreenF = (1 + UR*DRVec*VR)^{-1} = I - UR * (1 + DRVec*VR*UR)^{-1} * DRVec*VR
!                                     = I - UR * (1 + DRmax*DRmin*VR*UR)^{-1} * DRmax*DRmin*VR
!                                     = I - UR * [(DRmax)^{-1} + DRmin*VR*UR]^{-1} * DRmin*VR
!
!     zDet = det(1 + UR*DRVec*VR) = det(1 + DRVec*VR*UR) = det(1 + DRmax*DRmin*VR*UR)
!          = det(DRmax) * det[(DRmax)^{-1} + DRmin*VR*UR]
!
!     LogzDet = log(zDet) = log[det(DRmax)] + log{det[(DRmax)^{-1} + DRmin*VR*UR]}
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      complex(rp) LogzDet                      ! Log of the determinant of (1 + UR*DRVec*VR)
      complex(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)     DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      complex(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      complex(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)      
      complex(rp), allocatable :: VMatTmpt(:, :, :)
      complex(rp), allocatable :: UVTmpMat(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!__________________ Main calculations of G = I - UR * [(DRmax)^{-1} + DRmin*VR*UR]^{-1} * DRmin*VR ____________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
            allocate(VMatTmpt(NumNS, NumNS, NmSpn)); VMatTmpt = rp_Zzero
      allocate(UVTmpMat(NumNS, NumNS, NmSpn)); UVTmpMat = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Split the rescaled DRVec into DRMin and DRMax ____________________
!_____________________ Calculate Logdet of DRMax matrix _________________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_SgleD(NT, DRVec, LogScaleRt, DRMinMat, DRMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green's function using the following formula ________________
!___________ GreenF = 1_{Ns} - UR * [(DRmax)^{-1} + DRmin*VR*UR]^{-1} * DRmin*VR __________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate UVTmpMat = VR * UR _____________________________________
!________________________________________________________________________________________
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, VR, UR, rp_Zzero, UVTmpMat)
!________________________________________________________________________________________         
!_________________ (1) UVTmpMat = (DRmax)^{-1} + DRmin * UVTmpMat _______________________
!_____________________ VMatTmpt = DRmin * VR ____________________________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               VMatTmpt(I1, I2, SpnInd) = DRMinMat(I1, SpnInd) *       VR(I1, I2, SpnInd)
               UVTmpMat(I1, I2, SpnInd) = DRMinMat(I1, SpnInd) * UVTmpMat(I1, I2, SpnInd)
            enddo
            UVTmpMat(I2, I2, SpnInd) = UVTmpMat(I2, I2, SpnInd) + DRMaxInv(I2, SpnInd)
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) VMatTmpt = UVTmpMat^{-1} * VMatTmpt ______________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call zMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVTmpMat, VMatTmpt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = I - UR * VMatTmpt _______________________________________
!________________________________________________________________________________________
      GreenF = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, -rp_Z_One, UR, VMatTmpt, rp_Zzero, GreenF)
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, NumNS, +1
            GreenF(I1, I1, SpnInd) = rp_Z_One + GreenF(I1, I1, SpnInd)
         enddo
      enddo
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DRMinMat)) deallocate(DRMinMat)
      if(allocated(DRMaxInv)) deallocate(DRMaxInv)
      if(allocated(VMatTmpt)) deallocate(VMatTmpt)
            if(allocated(UVTmpMat)) deallocate(UVTmpMat)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStatcCB
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStatcC0(LogScaleLt, VL, DLVec, UL, LogzDet, GreenF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStatcC0(LogScaleLt, VL, DLVec, UL, LogzDet, GreenF) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of VL, DLVec, 
!                 UL as G = (1 + VL*DLVec*UL)^{-1}.  Single-particle Green's function at Tau = 0. 
! KEYWORDS: Calculate static Green's function, Complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     GreenF = (1 + VL*DLVec*UL)^{-1} = I - VL * (1 + DLVec*UL*VL)^{-1} * DLVec*UL
!                                     = I - VL * (1 + DLmax*DLmin*UL*VL)^{-1} * DLmax*DLmin*UL
!                                     = I - VL * [(DLmax)^{-1} + DLmin*UL*VL]^{-1} * DLmin*UL
!
!     zDet = det(1 + VL*DLVec*UL) = (1 + DLVec*UL*VL) = det(1 + DLmax*DLmin*UL*VL)
!          = det(DLmax) * det[(DLmax)^{-1} + DLmin*UL*VL]
!
!     LogzDet = log(zDet) = log[det(DLmax)] + log{det[(DLmax)^{-1} + DLmin*UL*VL]}
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
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet                     ! The log of Det[1+B(\beta, 0)]
      complex(rp)     VL(NumNS, NumNS, NmSpn) ! The VL matrix for B(\beta, tau)
      real(rp)     DLVec(NumNS       , NmSpn) ! The DL vector for B(\beta, tau) 
      complex(rp)     UL(NumNS, NumNS, NmSpn) ! The UL matrix for B(\beta, tau)
      complex(rp) GreenF(NumNS, NumNS, NmSpn) ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      complex(rp), allocatable :: UMatTmpt(:, :, :)
      complex(rp), allocatable :: UVTmpMat(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!__________________ Main calculations of G = I - VL * [(DLmax)^{-1} + DLmin*UL*VL]^{-1} * DLmin*UL ____________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
            allocate(UMatTmpt(NumNS, NumNS, NmSpn)); UMatTmpt = rp_Zzero
      allocate(UVTmpMat(NumNS, NumNS, NmSpn)); UVTmpMat = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Split the rescaled DLVec into DLMin and DLMax ____________________
!_____________________ Calculate Logdet of DLMax matrix _________________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_SgleD(LTrot, DLVec, LogScaleLt, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green's function using the following formula ________________
!___________ GreenF = 1_{Ns} - VL * [(DLmax)^{-1} + DLmin*UL*VL]^{-1} * DLmin*UL __________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate UVTmpMat = UL * VL _____________________________________
!________________________________________________________________________________________
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, UL, VL, rp_Zzero, UVTmpMat)
!________________________________________________________________________________________         
!_________________ (1) UVTmpMat = (DLmax)^{-1} + DLmin * UVTmpMat _______________________
!_____________________ UMatTmpt = DLmin * UL ____________________________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               UMatTmpt(I1, I2, SpnInd) = DLMinMat(I1, SpnInd) *       UL(I1, I2, SpnInd)
               UVTmpMat(I1, I2, SpnInd) = DLMinMat(I1, SpnInd) * UVTmpMat(I1, I2, SpnInd)
            enddo
            UVTmpMat(I2, I2, SpnInd) = UVTmpMat(I2, I2, SpnInd) + DLMaxInv(I2, SpnInd)
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) UMatTmpt = UVTmpMat^{-1} * UMatTmpt ______________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call zMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVTmpMat, UMatTmpt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = I - VL * UMatTmpt _______________________________________
!________________________________________________________________________________________
      GreenF = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, -rp_Z_One, VL, UMatTmpt, rp_Zzero, GreenF)   
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, NumNS, +1
            GreenF(I1, I1, SpnInd) = rp_Z_One + GreenF(I1, I1, SpnInd)
         enddo
      enddo
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat)) deallocate(DLMinMat)
      if(allocated(DLMaxInv)) deallocate(DLMaxInv)
      if(allocated(UMatTmpt)) deallocate(UMatTmpt)
            if(allocated(UVTmpMat)) deallocate(UVTmpMat)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStatcC0
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStaticC_LR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStaticC_LR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, VR, 
!                 VL, DLVec, UL as G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. 
! KEYWORDS: Calculate static Green's function, Complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. In this calculation, we have applied some technique to 
!            stablize the numerical calculations.
!
!     Here, UR, UL are unitary matrices, but VR, VL are not.
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!       = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!       = UL^{-1}*DLMax^{-1} * [DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]^{-1} 
!                            * DRMax^{-1}*UR^{-1}
!
!     dDet = det[1 + B(beta, 0)] = det(1 + VL*DLVec*UL*UR*DRVec*VR) = det(1 + UR*DRVec*VR * VL*DLVec*UL)
!          = det(UR*DRMax) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                         * det(DLMax*UL)
!          = det(UL*UR) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                  * det(DRMax*DLMax)
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet
      complex(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)     DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      complex(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      complex(rp)  VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)     DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      complex(rp)  UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      complex(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      complex(rp), allocatable :: ULInvTmpMt(:, :, :)
      complex(rp), allocatable :: URInvTmpMt(:, :, :)
      complex(rp), allocatable :: VRVLTmpMat(:, :, :)
      complex(rp), allocatable :: ULURTmpMat(:, :, :)
      complex(rp), allocatable :: UVInvTmpMt(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
            allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = rp_Zzero
            allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = rp_Zzero
      allocate(VRVLTmpMat(NumNS, NumNS, NmSpn)); VRVLTmpMat = rp_Zzero
      allocate(ULURTmpMat(NumNS, NumNS, NmSpn)); ULURTmpMat = rp_Zzero
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green Function matrix _______________________________________
!_______________ G = UL^{-1} * DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!___________________ + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1} * UR^{-1} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate ULURTmpMat = UL * UR ___________________________________
!_______________________________ VRVLTmpMat = VR * VL ___________________________________
!________________________________________________________________________________________
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, UL, UR, rp_Zzero, ULURTmpMat)
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, VR, VL, rp_Zzero, VRVLTmpMat)
!________________________________________________________________________________________         
!_________________ (1) ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                + DRMin*(VR*VL)*DLMin _____________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULInvTmpMt(I1, I2, SpnInd) = conjg(UL(I2, I1, SpnInd)) * DLMaxInv(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * conjg(UR(I2, I1, SpnInd))
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * conjg(ULURTmpMat(I2, I1, SpnInd))* DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * VRVLTmpMat(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call zMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = ULInvTmpMt * URInvTmpMt _________________________________
!________________________________________________________________________________________
      GreenF = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, ULInvTmpMt, URInvTmpMt, rp_Zzero, GreenF)
!**************************************************************************************************     
!___________________ 3. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate log[det(UL*UR)] ________________________________________
!________________________________________________________________________________________ 
      LogzDetTmp = rp_Zzero
      call zMatDet_LogDet_QMC(NumNS, ULURTmpMat, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (1) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 4. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
            if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
            if(allocated(VRVLTmpMat)) deallocate(VRVLTmpMat)
      if(allocated(ULURTmpMat)) deallocate(ULURTmpMat)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStaticC_LR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStaticC_UR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStaticC_UR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, VR, 
!                 VL, DLVec, UL as G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. 
! KEYWORDS: Calculate static Green's function, Complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. In this calculation, we have applied some technique to 
!            stablize the numerical calculations.
!
!     Here, UR is unitary matrix, and UL is general matrix.
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!       = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!       = UL^{-1}*DLMax^{-1} * [DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]^{-1} 
!                            * DRMax^{-1}*UR^{-1}
!
!     dDet = det[1 + B(beta, 0)] = det(1 + VL*DLVec*UL*UR*DRVec*VR) = det(1 + UR*DRVec*VR * VL*DLVec*UL)
!          = det(UR*DRMax) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                         * det(DLMax*UL)
!          = det(UL*UR) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                  * det(DRMax*DLMax)
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet
      complex(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)     DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      complex(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      complex(rp)  VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)     DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      complex(rp)  UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      complex(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      complex(rp), allocatable :: VRVLTmpMat(:, :, :)
      complex(rp), allocatable :: ULInvTmpMt(:, :, :)
      complex(rp), allocatable :: URInvTmpMt(:, :, :)
      complex(rp), allocatable :: UVInvTmpMt(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
      allocate(VRVLTmpMat(NumNS, NumNS, NmSpn)); VRVLTmpMat = rp_Zzero
            allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = rp_Zzero
            allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = rp_Zzero
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green Function matrix _______________________________________
!_______________ G = UL^{-1} * DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!___________________ + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1} * UR^{-1} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate ULInvTmpMt = UL^{-1} ___________________________________
!_______________________________ Determinant of UR matrix _______________________________
!_______________________________ UVInvTmpMt = (UL*UR)^{-1} = UR^+ * UL^{-1} _____________
!_______________________________ VRVLTmpMat = VR * VL ___________________________________
!________________________________________________________________________________________      
      call zMat_Copy_QMC(UL(1, 1, 1), ULInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call zMatInv_LogDet_QMC(NumNS, ULInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp

      call zMat_Copy_QMC(UR(1, 1, 1), URInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call zMatDet_LogDet_QMC(NumNS, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
      
      call zMatPrd_CN_QMC(NumNS, NumNS, NumNS, rp_Z_One, UR, ULInvTmpMt, rp_Zzero, UVInvTmpMt)

      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, VR, VL, rp_Zzero, VRVLTmpMat)
!________________________________________________________________________________________         
!_________________ (1) ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                 + DRMin*(VR*VL)*DLMin ____________________________
!________________________________________________________________________________________
      URInvTmpMt = rp_Zzero
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULInvTmpMt(I1, I2, SpnInd) = ULInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * conjg(UR(I2, I1, SpnInd))
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UVInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * VRVLTmpMat(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call zMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = ULInvTmpMt * URInvTmpMt _________________________________
!________________________________________________________________________________________
      GreenF = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, ULInvTmpMt, URInvTmpMt, rp_Zzero, GreenF)
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
      if(allocated(VRVLTmpMat)) deallocate(VRVLTmpMat)
            if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStaticC_UR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStaticC_UL(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStaticC_UL(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, VR, 
!                 VL, DLVec, UL as G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. 
! KEYWORDS: Calculate static Green's function, Complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. In this calculation, we have applied some technique to 
!            stablize the numerical calculations.
!
!     Here, UL is unitary matrix, and UR is general matrix.
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!       = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!       = UL^{-1}*DLMax^{-1} * [DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]^{-1} 
!                            * DRMax^{-1}*UR^{-1}
!
!     dDet = det[1 + B(beta, 0)] = det(1 + VL*DLVec*UL*UR*DRVec*VR) = det(1 + UR*DRVec*VR * VL*DLVec*UL)
!          = det(UR*DRMax) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                         * det(DLMax*UL)
!          = det(UL*UR) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                  * det(DRMax*DLMax)
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet
      complex(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)     DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      complex(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      complex(rp)  VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)     DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      complex(rp)  UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      complex(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      complex(rp), allocatable :: VRVLTmpMat(:, :, :)
      complex(rp), allocatable :: ULInvTmpMt(:, :, :)
      complex(rp), allocatable :: URInvTmpMt(:, :, :)
      complex(rp), allocatable :: UVInvTmpMt(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
      allocate(VRVLTmpMat(NumNS, NumNS, NmSpn)); VRVLTmpMat = rp_Zzero
            allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = rp_Zzero
            allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = rp_Zzero
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green Function matrix _______________________________________
!_______________ G = UL^{-1} * DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!___________________ + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1} * UR^{-1} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate Determinant of UL matrix _______________________________
!_______________________________ URInvTmpMt = UR^{-1} ___________________________________
!_______________________________ UVInvTmpMt = (UL*UR)^{-1} = UR^{-1} * UL^+ _____________
!_______________________________ VRVLTmpMat = VR * VL ___________________________________
!________________________________________________________________________________________
      call zMat_Copy_QMC(UL(1, 1, 1), ULInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call zMatDet_LogDet_QMC(NumNS, ULInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp

      call zMat_Copy_QMC(UR(1, 1, 1), URInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call zMatInv_LogDet_QMC(NumNS, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
      
      call zMatPrd_NC_QMC(NumNS, NumNS, NumNS, rp_Z_One, URInvTmpMt, UL, rp_Zzero, UVInvTmpMt)
      
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, VR, VL, rp_Zzero, VRVLTmpMat)
!________________________________________________________________________________________         
!_________________ (1) ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                 + DRMin*(VR*VL)*DLMin ____________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULInvTmpMt(I1, I2, SpnInd) = conjg(UL(I2, I1, SpnInd)) * DLMaxInv(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * URInvTmpMt(I1, I2, SpnInd)
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UVInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * VRVLTmpMat(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call zMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = ULInvTmpMt * URInvTmpMt _________________________________
!________________________________________________________________________________________
      GreenF = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, ULInvTmpMt, URInvTmpMt, rp_Zzero, GreenF)
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
      if(allocated(VRVLTmpMat)) deallocate(VRVLTmpMat)
            if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStaticC_UL
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStaticC_NO(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStaticC_NO(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, VR, 
!                 VL, DLVec, UL as G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. 
! KEYWORDS: Calculate static Green's function, Complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. In this calculation, we have applied some technique to 
!            stablize the numerical calculations.
!
!     Here, UL and UR are general matrices, rather than unitary matrices.
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!       = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!       = UL^{-1}*DLMax^{-1} * [DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]^{-1} 
!                            * DRMax^{-1}*UR^{-1}
!
!     dDet = det[1 + B(beta, 0)] = det(1 + VL*DLVec*UL*UR*DRVec*VR) = det(1 + UR*DRVec*VR * VL*DLVec*UL)
!          = det(UR*DRMax) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                         * det(DLMax*UL)
!          = det(UL*UR) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                  * det(DRMax*DLMax)
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet
      complex(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)     DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      complex(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      complex(rp)  VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)     DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      complex(rp)  UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      complex(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      complex(rp), allocatable :: VRVLTmpMat(:, :, :)
      complex(rp), allocatable :: ULInvTmpMt(:, :, :)
      complex(rp), allocatable :: URInvTmpMt(:, :, :)
      complex(rp), allocatable :: UVInvTmpMt(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
      allocate(VRVLTmpMat(NumNS, NumNS, NmSpn)); VRVLTmpMat = rp_Zzero
            allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = rp_Zzero
            allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = rp_Zzero
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green Function matrix _______________________________________
!_______________ G = UL^{-1} * DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!___________________ + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1} * UR^{-1} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate ULInvTmpMt = UL^{-1} ___________________________________
!_______________________________ URInvTmpMt = UR^{-1} ___________________________________
!_______________________________ UVInvTmpMt = (UL*UR)^{-1} = UR^{-1} * UL^{-1} __________
!_______________________________ VRVLTmpMat = VR * VL ___________________________________
!________________________________________________________________________________________
      call zMat_Copy_QMC(UL(1, 1, 1), ULInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call zMatInv_LogDet_QMC(NumNS, ULInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp

      call zMat_Copy_QMC(UR(1, 1, 1), URInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call zMatInv_LogDet_QMC(NumNS, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
      
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, URInvTmpMt, ULInvTmpMt, rp_Zzero, UVInvTmpMt)
      
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, VR, VL, rp_Zzero, VRVLTmpMat)
!________________________________________________________________________________________         
!_________________ (1) ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                 + DRMin*(VR*VL)*DLMin ____________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULInvTmpMt(I1, I2, SpnInd) = ULInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * URInvTmpMt(I1, I2, SpnInd)
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UVInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * VRVLTmpMat(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call zMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = ULInvTmpMt * URInvTmpMt _________________________________
!________________________________________________________________________________________
      GreenF = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, ULInvTmpMt, URInvTmpMt, rp_Zzero, GreenF)
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
      if(allocated(VRVLTmpMat)) deallocate(VRVLTmpMat)
            if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStaticC_NO
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFDynamcC_LR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GrnFTT, GrnFT0, GrnF0T) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFDynamcC_LR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GrnFTT, GrnFT0, GrnF0T) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to recompute the dynamic single-particle Green's function as G(tau, 0) 
!                  and G(0, tau) matrices. 
! KEYWORDS: Recompute Dynamic Green's function, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     By QRP for numerical stabilization, UR, UL are unitary matrices, but VR, VL are not.
!     By SVD for numerical stabilization, UR, UL, VR, VL are all unitary matrices.
!     This subroutine can be used for both cases of QRP and SVD decompositions for numerical stabilization, since
!           we only use the property that UR, UL are unitary matrices, during the calculations.
!
!     G(tau, tau) = [I + B(tau, 0)*B(beta, tau)]^{-1} = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} 
!                 = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!                 = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!                 = UL^{-1}*DLMax^{-1} 
!                           * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin ]^{-1} 
!                           * DRMax^{-1}*UR^{-1}
!
!     G(tau, 0) = [I + B(tau, 0)*B(beta, tau)]^{-1} * B(tau, 0)
!               = (I + UR*DRVec*VR * VL*DLVec*UL)^{-1} * UR*DRVec*VR
!               = UL^{-1}*DLMax^{-1} 
!                           * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin ]^{-1} 
!                           * DRMin*VR
!
!     G(0, tau) = [B(tau, 0)]^{-1} * { I - [I + B(tau, 0)*B(beta, tau)]^{-1} }
!               = [I + B(beta, tau)*B(tau, 0)]^{-1} * B(beta, tau)
!               = B(beta, tau) * [I + B(tau, 0)*B(beta, tau)]^{-1}
!               = VL*DLVec*UL * (I + UR*DRVec*VR * VL*DLVec*UL)^{-1}
!               = VL*DLMin 
!                           * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin ]^{-1} 
!                           * DRMax^{-1}*UR^{-1}
!
!     Since there are some common terms in calculating these two dynamic single-particle Green's functions, we
!               perform the calculations in following steps.
!       (0) Calculate (UL*UR) and (VR*VL).
!       (1) Calculate DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin and
!                          UL^{-1}*DLMax^{-1}, (DRMin*VR), obtain G(tau, 0).
!       (2) Calculate (VL*DLMin) and (DRMax^{-1}*UR^{-1}), obtain G(0, tau).
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
      integer NT
      real(rp) LogScaleRt(NmSpn), LogScaleLt(NmSpn)
      complex(rp) LogzDet
      complex(rp) UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)    DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      complex(rp) VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      complex(rp) VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)    DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      complex(rp) UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      complex(rp) GrnFTT(NumNS, NumNS, NmSpn)  ! The output G(tau, tau) matrix
      complex(rp) GrnFT0(NumNS, NumNS, NmSpn)  ! The output G(tau,   0) matrix
      complex(rp) GrnF0T(NumNS, NumNS, NmSpn)  ! The output G(  0, tau) matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      complex(rp), allocatable :: ULURTmpMat(:, :, :)
      complex(rp), allocatable :: UVInvTmpMt(:, :, :)
      complex(rp), allocatable :: ULInvTmpMt(:, :, :)
      complex(rp), allocatable :: VLTmpMatrx(:, :, :)
      complex(rp), allocatable :: URInvTmpMt(:, :, :)
      complex(rp), allocatable :: VRTmpMatrx(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFDyn = TimsGFDyn + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G(tau, 0) and G(0, tau) matrices ________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
      allocate(ULURTmpMat(NumNS, NumNS, NmSpn)); ULURTmpMat = rp_Zzero
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = rp_Zzero
      allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = rp_Zzero
      allocate(VLTmpMatrx(NumNS, NumNS, NmSpn)); VLTmpMatrx = rp_Zzero
      allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = rp_Zzero
            allocate(VRTmpMatrx(NumNS, NumNS, NmSpn)); VRTmpMatrx = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate all the pieces for Green's Function matrices ____________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate ULURTmpMat = UL * UR ___________________________________
!_______________________________ UVInvTmpMt = VR * VL ___________________________________
!________________________________________________________________________________________
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, UL, UR, rp_Zzero, ULURTmpMat)
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, VR, VL, rp_Zzero, UVInvTmpMt)
!________________________________________________________________________________________         
!_________________ (1) UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                + DRMin*(VR*VL)*DLMin _____________________________
!_____________________ ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ VLTmpMatrx = VL * DLMin __________________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ VRTmpMatrx = DRMin * VR __________________________________________
!________________________________________________________________________________________      
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd)*conjg(ULURTmpMat(I2, I1, SpnInd))*DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * UVInvTmpMt(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
               ULInvTmpMt(I1, I2, SpnInd) = conjg(UL(I2, I1, SpnInd)) * DLMaxInv(I2, SpnInd)
               VLTmpMatrx(I1, I2, SpnInd) =       VL(I1, I2, SpnInd)  * DLMinMat(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * conjg(UR(I2, I1, SpnInd))
               VRTmpMatrx(I1, I2, SpnInd) = DRMinMat(I1, SpnInd) *       VR(I1, I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) Calculate the log determinant log[det(UL*UR)] ____________________
!________________________________________________________________________________________ 
      LogzDetTmp = rp_Zzero
      call zMatDet_LogDet_QMC(NumNS, ULURTmpMat, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!**************************************************************************************************     
!___________________ 2. Calculate G(tau, tau) and G(0, tau) matrices ______________________________
!________ G(tau, tau) = UL^{-1}*DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} _________________
!______________________      + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1}*UR^{-1} ____________________
!________ G(  0, tau) = VL*DLMin * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} ___________________________
!______________________      + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1}*UR^{-1} ____________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      ULURTmpMat = rp_Zzero
      call zMat_Copy_QMC(UVInvTmpMt(1, 1, 1), ULURTmpMat(1, 1, 1))
      LogzDetTmp = rp_Zzero
      call zMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, ULURTmpMat, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (1) GrnFTT = ULInvTmpMt * UVInvTmpMt^{-1} * URInvTmpMt ________________
!________________________________________________________________________________________
      GrnFTT = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, ULInvTmpMt, URInvTmpMt, rp_Zzero, GrnFTT)
!________________________________________________________________________________________         
!_________________ (2) GrnF0T = VLTmpMatrx * UVInvTmpMt^{-1} * URInvTmpMt ________________
!________________________________________________________________________________________
      GrnF0T = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, VLTmpMatrx, URInvTmpMt, rp_Zzero, GrnF0T)
!**************************************************************************************************     
!___________________ 3. Calculate G(tau, 0) matrice _______________________________________________
!_____________ G(tau, 0) = UL^{-1}*DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} ______________
!_________________________     + DRMin*(VR*VL)*DLMin ]^{-1} * DRMin*VR ____________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) VRTmpMatrx = UVInvTmpMt^{-1} * VRTmpMatrx ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call zMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, VRTmpMatrx, LogzDetTmp)
!________________________________________________________________________________________         
!_________________ (1) GrnFT0 = ULInvTmpMt * UVInvTmpMt^{-1} * VRTmpMatrx ________________
!________________________________________________________________________________________
      GrnFT0 = rp_Zzero
      call zMatPrd_NN_QMC(NumNS, NumNS, NumNS, rp_Z_One, ULInvTmpMt, VRTmpMatrx, rp_Zzero, GrnFT0)
!**************************************************************************************************     
!___________________ 4. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet = real(LogzDet); ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 5. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
      if(allocated(ULURTmpMat)) deallocate(ULURTmpMat)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
            if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
      if(allocated(VLTmpMatrx)) deallocate(VLTmpMatrx)
      if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(VRTmpMatrx)) deallocate(VRTmpMatrx)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFDyn = TimeGFDyn + TimeIntrvl(time1, time2)
      
   end subroutine GrFDynamcC_LR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
!################################################# For Real Version #####################################################
!################################################# For Real Version #####################################################
!################################################# For Real Version #####################################################
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStatcRB(NT, LogScaleRt, UR, DRVec, VR, LogzDet, GreenF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStatcRB(NT, LogScaleRt, UR, DRVec, VR, LogzDet, GreenF) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, 
!                 VR as G = (1 + UR*DRVec*VR)^{-1}.  Partial single-particle Green's function. 
! KEYWORDS: Calculate static Green's function, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     GreenF = (1 + UR*DRVec*VR)^{-1} = I - UR * (1 + DRVec*VR*UR)^{-1} * DRVec*VR
!                                     = I - UR * (1 + DRmax*DRmin*VR*UR)^{-1} * DRmax*DRmin*VR
!                                     = I - UR * [(DRmax)^{-1} + DRmin*VR*UR]^{-1} * DRmin*VR
!
!     zDet = det(1 + UR*DRVec*VR) = det(1 + DRVec*VR*UR) = det(1 + DRmax*DRmin*VR*UR)
!          = det(DRmax) * det[(DRmax)^{-1} + DRmin*VR*UR]
!
!     LogzDet = log(zDet) = log[det(DRmax)] + log{det[(DRmax)^{-1} + DRmin*VR*UR]}
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      complex(rp) LogzDet                   ! Log of the determinant of (1 + UR*DRVec*VR)
      real(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)  DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      real(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      real(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: VMatTmpt(:, :, :)
      real(rp), allocatable :: UVTmpMat(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!__________________ Main calculations of G = I - UR * [(DRmax)^{-1} + DRmin*VR*UR]^{-1} * DRmin*VR ____________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
            allocate(VMatTmpt(NumNS, NumNS, NmSpn)); VMatTmpt = 0.0_rp
      allocate(UVTmpMat(NumNS, NumNS, NmSpn)); UVTmpMat = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) Split the rescaled DRVec into DRMin and DRMax ____________________
!_____________________ Calculate Logdet of DRMax matrix _________________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_SgleD(NT, DRVec, LogScaleRt, DRMinMat, DRMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green's function using the following formula ________________
!___________ GreenF = 1_{Ns} - UR * [(DRmax)^{-1} + DRmin*VR*UR]^{-1} * DRmin*VR __________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate UVTmpMat = VR * UR _____________________________________
!________________________________________________________________________________________
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, VR, UR, 0.0_rp, UVTmpMat)
!________________________________________________________________________________________         
!_________________ (1) UVTmpMat = (DRmax)^{-1} + DRmin * UVTmpMat _______________________
!_____________________ VMatTmpt = DRmin * VR ____________________________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               VMatTmpt(I1, I2, SpnInd) = DRMinMat(I1, SpnInd) *       VR(I1, I2, SpnInd)
               UVTmpMat(I1, I2, SpnInd) = DRMinMat(I1, SpnInd) * UVTmpMat(I1, I2, SpnInd)
            enddo
            UVTmpMat(I2, I2, SpnInd) = UVTmpMat(I2, I2, SpnInd) + DRMaxInv(I2, SpnInd)
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) VMatTmpt = UVTmpMat^{-1} * VMatTmpt ______________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call dMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVTmpMat, VMatTmpt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = I - UR * VMatTmpt _______________________________________
!________________________________________________________________________________________
      GreenF = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, -1.0_rp, UR, VMatTmpt, 0.0_rp, GreenF)
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, NumNS, +1
            GreenF(I1, I1, SpnInd) = 1.0_rp + GreenF(I1, I1, SpnInd)
         enddo
      enddo
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DRMinMat)) deallocate(DRMinMat)
      if(allocated(DRMaxInv)) deallocate(DRMaxInv)
      if(allocated(VMatTmpt)) deallocate(VMatTmpt)
            if(allocated(UVTmpMat)) deallocate(UVTmpMat)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStatcRB
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStatcR0(LogScaleLt, VL, DLVec, UL, LogzDet, GreenF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStatcR0(LogScaleLt, VL, DLVec, UL, LogzDet, GreenF) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of VL, DLVec, 
!                 UL as G = (1 + VL*DLVec*UL)^{-1}.  Single-particle Green's function at Tau = 0. 
! KEYWORDS: Calculate static Green's function, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     GreenF = (1 + VL*DLVec*UL)^{-1} = I - VL * (1 + DLVec*UL*VL)^{-1} * DLVec*UL
!                                     = I - VL * (1 + DLmax*DLmin*UL*VL)^{-1} * DLmax*DLmin*UL
!                                     = I - VL * [(DLmax)^{-1} + DLmin*UL*VL]^{-1} * DLmin*UL
!
!     zDet = det(1 + VL*DLVec*UL) = (1 + DLVec*UL*VL) = det(1 + DLmax*DLmin*UL*VL)
!          = det(DLmax) * det[(DLmax)^{-1} + DLmin*UL*VL]
!
!     LogzDet = log(zDet) = log[det(DLmax)] + log{det[(DLmax)^{-1} + DLmin*UL*VL]}
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
      real(rp) LogScaleLt(NmSpn)
      real(rp)     VL(NumNS, NumNS, NmSpn) ! The VL matrix for B(\beta, tau)
      real(rp)  DLVec(NumNS       , NmSpn) ! The DL vector for B(\beta, tau) 
      real(rp)     UL(NumNS, NumNS, NmSpn) ! The UL matrix for B(\beta, tau)
      real(rp) GreenF(NumNS, NumNS, NmSpn) ! The output static single-particle Green's function matrix
      complex(rp) LogzDet                  ! The log of Det[1+B(\beta, 0)]
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      real(rp), allocatable :: UMatTmpt(:, :, :)
      real(rp), allocatable :: UVTmpMat(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!__________________ Main calculations of G = I - VL * [(DLmax)^{-1} + DLmin*UL*VL]^{-1} * DLmin*UL ____________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
            allocate(UMatTmpt(NumNS, NumNS, NmSpn)); UMatTmpt = 0.0_rp
      allocate(UVTmpMat(NumNS, NumNS, NmSpn)); UVTmpMat = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) Split the rescaled DLVec into DLMin and DLMax ____________________
!_____________________ Calculate Logdet of DLMax matrix _________________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_SgleD(LTrot, DLVec, LogScaleLt, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green's function using the following formula ________________
!___________ GreenF = 1_{Ns} - VL * [(DLmax)^{-1} + DLmin*UL*VL]^{-1} * DLmin*UL __________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate UVTmpMat = UL * VL _____________________________________
!________________________________________________________________________________________
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, UL, VL, 0.0_rp, UVTmpMat)
!________________________________________________________________________________________         
!_________________ (1) UVTmpMat = (DLmax)^{-1} + DLmin * UVTmpMat _______________________
!_____________________ UMatTmpt = DLmin * UL ____________________________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               UMatTmpt(I1, I2, SpnInd) = DLMinMat(I1, SpnInd) *       UL(I1, I2, SpnInd)
               UVTmpMat(I1, I2, SpnInd) = DLMinMat(I1, SpnInd) * UVTmpMat(I1, I2, SpnInd)
            enddo
            UVTmpMat(I2, I2, SpnInd) = UVTmpMat(I2, I2, SpnInd) + DLMaxInv(I2, SpnInd)
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) UMatTmpt = UVTmpMat^{-1} * UMatTmpt ______________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call dMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVTmpMat, UMatTmpt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = I - VL * UMatTmpt _______________________________________
!________________________________________________________________________________________
      GreenF = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, -1.0_rp, VL, UMatTmpt, 0.0_rp, GreenF) 
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, NumNS, +1
            GreenF(I1, I1, SpnInd) = 1.0_rp + GreenF(I1, I1, SpnInd)
         enddo
      enddo
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat)) deallocate(DLMinMat)
      if(allocated(DLMaxInv)) deallocate(DLMaxInv)
      if(allocated(UMatTmpt)) deallocate(UMatTmpt)
            if(allocated(UVTmpMat)) deallocate(UVTmpMat)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStatcR0
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStaticR_LR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStaticR_LR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF)  
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, VR, 
!                 VL, DLVec, UL as G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. 
! KEYWORDS: Calculate static Green's function, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. In this calculation, we have applied some technique to 
!            stablize the numerical calculations.
!
!     Here, UR, UL are unitary matrices, but VR, VL are not.
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!       = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!       = UL^{-1}*DLMax^{-1} * [DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]^{-1} 
!                            * DRMax^{-1}*UR^{-1}
!
!     dDet = det[1 + B(beta, 0)] = det(1 + VL*DLVec*UL*UR*DRVec*VR) = det(1 + UR*DRVec*VR * VL*DLVec*UL)
!          = det(UR*DRMax) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                         * det(DLMax*UL)
!          = det(UL*UR) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                  * det(DRMax*DLMax)
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet
      real(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)  DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      real(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      real(rp)  VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)  DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      real(rp)  UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      real(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      real(rp), allocatable :: ULInvTmpMt(:, :, :)
      real(rp), allocatable :: URInvTmpMt(:, :, :)
      real(rp), allocatable :: VRVLTmpMat(:, :, :)
      real(rp), allocatable :: ULURTmpMat(:, :, :)
      real(rp), allocatable :: UVInvTmpMt(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
            allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = 0.0_rp
            allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = 0.0_rp
      allocate(VRVLTmpMat(NumNS, NumNS, NmSpn)); VRVLTmpMat = 0.0_rp
      allocate(ULURTmpMat(NumNS, NumNS, NmSpn)); ULURTmpMat = 0.0_rp
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green Function matrix _______________________________________
!_______________ G = UL^{-1} * DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!___________________ + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1} * UR^{-1} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate ULURTmpMat = UL * UR ___________________________________
!_______________________________ VRVLTmpMat = VR * VL ___________________________________
!________________________________________________________________________________________
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, UL, UR, 0.0_rp, ULURTmpMat)
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, VR, VL, 0.0_rp, VRVLTmpMat)
!________________________________________________________________________________________         
!_________________ (1) ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                + DRMin*(VR*VL)*DLMin _____________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULInvTmpMt(I1, I2, SpnInd) = UL(I2, I1, SpnInd) * DLMaxInv(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UR(I2, I1, SpnInd)
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * ULURTmpMat(I2, I1, SpnInd) * DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * VRVLTmpMat(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call dMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = ULInvTmpMt * URInvTmpMt _________________________________
!________________________________________________________________________________________
      GreenF = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, ULInvTmpMt, URInvTmpMt, 0.0_rp, GreenF)
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate log[det(UL*UR)] ________________________________________
!________________________________________________________________________________________ 
      LogzDetTmp = rp_Zzero
      call dMatDet_LogDet_QMC(NumNS, ULURTmpMat, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (1) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
            if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
            if(allocated(VRVLTmpMat)) deallocate(VRVLTmpMat)
      if(allocated(ULURTmpMat)) deallocate(ULURTmpMat)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStaticR_LR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStaticR_UR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStaticR_UR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF)
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, VR, 
!                 VL, DLVec, UL as G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. 
! KEYWORDS: Calculate static Green's function, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. In this calculation, we have applied some technique to 
!            stablize the numerical calculations.
!
!     Here, UR is unitary matrix, and UL is general matrix.
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!       = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!       = UL^{-1}*DLMax^{-1} * [DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]^{-1} 
!                            * DRMax^{-1}*UR^{-1}
!
!     dDet = det[1 + B(beta, 0)] = det(1 + VL*DLVec*UL*UR*DRVec*VR) = det(1 + UR*DRVec*VR * VL*DLVec*UL)
!          = det(UR*DRMax) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                         * det(DLMax*UL)
!          = det(UL*UR) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                  * det(DRMax*DLMax)
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet
      real(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)  DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      real(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      real(rp)  VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)  DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      real(rp)  UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      real(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      real(rp), allocatable :: VRVLTmpMat(:, :, :)
      real(rp), allocatable :: ULInvTmpMt(:, :, :)
      real(rp), allocatable :: URInvTmpMt(:, :, :)
      real(rp), allocatable :: UVInvTmpMt(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
      allocate(VRVLTmpMat(NumNS, NumNS, NmSpn)); VRVLTmpMat = 0.0_rp
            allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = 0.0_rp
            allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = 0.0_rp
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green Function matrix _______________________________________
!_______________ G = UL^{-1} * DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!___________________ + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1} * UR^{-1} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate ULInvTmpMt = UL^{-1} ___________________________________
!_______________________________ Determinant of UR matrix _______________________________
!_______________________________ UVInvTmpMt = (UL*UR)^{-1} = UR^+ * UL^{-1} _____________
!_______________________________ VRVLTmpMat = VR * VL ___________________________________
!________________________________________________________________________________________
      call dMat_Copy_QMC(UL(1, 1, 1), ULInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call dMatInv_LogDet_QMC(NumNS, ULInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp

      call dMat_Copy_QMC(UR(1, 1, 1), URInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call dMatDet_LogDet_QMC(NumNS, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
      
      call dMatPrd_TN_QMC(NumNS, NumNS, NumNS, 1.0_rp, UR, ULInvTmpMt, 0.0_rp, UVInvTmpMt)

      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, VR, VL, 0.0_rp, VRVLTmpMat)
!________________________________________________________________________________________         
!_________________ (1) ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                 + DRMin*(VR*VL)*DLMin ____________________________
!________________________________________________________________________________________
      URInvTmpMt = 0.0_rp
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULInvTmpMt(I1, I2, SpnInd) = ULInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UR(I2, I1, SpnInd)
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UVInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * VRVLTmpMat(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call dMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = ULInvTmpMt * URInvTmpMt _________________________________
!________________________________________________________________________________________
      GreenF = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, ULInvTmpMt, URInvTmpMt, 0.0_rp, GreenF)
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
      if(allocated(VRVLTmpMat)) deallocate(VRVLTmpMat)
            if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStaticR_UR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStaticR_UL(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStaticR_UL(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF)
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, VR, 
!                 VL, DLVec, UL as G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. 
! KEYWORDS: Calculate static Green's function, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. In this calculation, we have applied some technique to 
!            stablize the numerical calculations.
!
!     Here, UL is unitary matrix, and UR is general matrix.
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!       = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!       = UL^{-1}*DLMax^{-1} * [DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]^{-1} 
!                            * DRMax^{-1}*UR^{-1}
!
!     dDet = det[1 + B(beta, 0)] = det(1 + VL*DLVec*UL*UR*DRVec*VR) = det(1 + UR*DRVec*VR * VL*DLVec*UL)
!          = det(UR*DRMax) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                         * det(DLMax*UL)
!          = det(UL*UR) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                  * det(DRMax*DLMax)
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet
      real(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)  DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      real(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      real(rp)  VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)  DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      real(rp)  UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      real(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      real(rp), allocatable :: VRVLTmpMat(:, :, :)
      real(rp), allocatable :: ULInvTmpMt(:, :, :)
      real(rp), allocatable :: URInvTmpMt(:, :, :)
      real(rp), allocatable :: UVInvTmpMt(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
      allocate(VRVLTmpMat(NumNS, NumNS, NmSpn)); VRVLTmpMat = 0.0_rp
            allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = 0.0_rp
            allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = 0.0_rp
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green Function matrix _______________________________________
!_______________ G = UL^{-1} * DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!___________________ + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1} * UR^{-1} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate Determinant of UL matrix _______________________________
!_______________________________ URInvTmpMt = UR^{-1} ___________________________________
!_______________________________ UVInvTmpMt = (UL*UR)^{-1} = UR^{-1} * UL^+ _____________
!_______________________________ VRVLTmpMat = VR * VL ___________________________________
!________________________________________________________________________________________      
      call dMat_Copy_QMC(UL(1, 1, 1), ULInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call dMatDet_LogDet_QMC(NumNS, ULInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp

      call dMat_Copy_QMC(UR(1, 1, 1), URInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call dMatInv_LogDet_QMC(NumNS, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
      
      call dMatPrd_NT_QMC(NumNS, NumNS, NumNS, 1.0_rp, URInvTmpMt, UL, 0.0_rp, UVInvTmpMt)

      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, VR, VL, 0.0_rp, VRVLTmpMat)
!________________________________________________________________________________________         
!_________________ (1) ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                 + DRMin*(VR*VL)*DLMin ____________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULInvTmpMt(I1, I2, SpnInd) = UL(I2, I1, SpnInd) * DLMaxInv(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * URInvTmpMt(I1, I2, SpnInd)
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UVInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * VRVLTmpMat(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call dMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = ULInvTmpMt * URInvTmpMt _________________________________
!________________________________________________________________________________________
      GreenF = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, ULInvTmpMt, URInvTmpMt, 0.0_rp, GreenF)
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
      if(allocated(VRVLTmpMat)) deallocate(VRVLTmpMat)
            if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStaticR_UL
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFStaticR_No(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFStaticR_No(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GreenF)
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to calculate the static Green's function from the quantities of UR, DRVec, VR, 
!                 VL, DLVec, UL as G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. 
! KEYWORDS: Calculate static Green's function, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1}. In this calculation, we have applied some technique to 
!            stablize the numerical calculations.
!
!     Here, UL and UR are general matrices, rather than unitary matrices.
!
!     G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!       = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!       = UL^{-1}*DLMax^{-1} * [DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]^{-1} 
!                            * DRMax^{-1}*UR^{-1}
!
!     dDet = det[1 + B(beta, 0)] = det(1 + VL*DLVec*UL*UR*DRVec*VR) = det(1 + UR*DRVec*VR * VL*DLVec*UL)
!          = det(UR*DRMax) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                         * det(DLMax*UL)
!          = det(UL*UR) * det[DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin]
!                                  * det(DRMax*DLMax)
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
      integer NT
      real(rp) LogScaleRt(NmSpn)
      real(rp) LogScaleLt(NmSpn)
      complex(rp) LogzDet
      real(rp)  UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp)  DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      real(rp)  VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      real(rp)  VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp)  DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      real(rp)  UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      real(rp) GreenF(NumNS, NumNS, NmSpn)  ! The output static single-particle Green's function matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________   
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      real(rp), allocatable :: VRVLTmpMat(:, :, :)
      real(rp), allocatable :: ULInvTmpMt(:, :, :)
      real(rp), allocatable :: URInvTmpMt(:, :, :)
      real(rp), allocatable :: UVInvTmpMt(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFSta = TimsGFSta + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
      allocate(VRVLTmpMat(NumNS, NumNS, NmSpn)); VRVLTmpMat = 0.0_rp
            allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = 0.0_rp
            allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = 0.0_rp
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate the Green Function matrix _______________________________________
!_______________ G = UL^{-1} * DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!___________________ + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1} * UR^{-1} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate ULInvTmpMt = UL^{-1} ___________________________________
!_______________________________ URInvTmpMt = UR^{-1} ___________________________________
!_______________________________ UVInvTmpMt = (UL*UR)^{-1} = UR^{-1} * UL^{-1} __________
!_______________________________ VRVLTmpMat = VR * VL ___________________________________
!________________________________________________________________________________________
      call dMat_Copy_QMC(UL(1, 1, 1), ULInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call dMatInv_LogDet_QMC(NumNS, ULInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp

      call dMat_Copy_QMC(UR(1, 1, 1), URInvTmpMt(1, 1, 1)) 
      LogzDetTmp = rp_Zzero
      call dMatInv_LogDet_QMC(NumNS, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
      
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, URInvTmpMt, ULInvTmpMt, 0.0_rp, UVInvTmpMt)

      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, VR, VL, 0.0_rp, VRVLTmpMat)
!________________________________________________________________________________________         
!_________________ (1) ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                 + DRMin*(VR*VL)*DLMin ____________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULInvTmpMt(I1, I2, SpnInd) = ULInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * URInvTmpMt(I1, I2, SpnInd)
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UVInvTmpMt(I1, I2, SpnInd) * DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * VRVLTmpMat(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call dMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (3) GreenF = ULInvTmpMt * URInvTmpMt _________________________________
!________________________________________________________________________________________
      GreenF = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, ULInvTmpMt, URInvTmpMt, 0.0_rp, GreenF)
!**************************************************************************************************     
!___________________ 2. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet =  real(LogzDet)
      ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 3. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
      if(allocated(VRVLTmpMat)) deallocate(VRVLTmpMat)
            if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFSta = TimeGFSta + TimeIntrvl(time1, time2)
      
   end subroutine GrFStaticR_No
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrFDynamcR_LR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GrnFTT, GrnFT0, GrnF0T) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrFDynamcR_LR(NT, LogScaleRt, LogScaleLt, UR, DRVec, VR, VL, DLVec, UL, LogzDet, GrnFTT, GrnFT0, GrnF0T) 
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to recompute the dynamic single-particle Green's function as G(tau, 0) 
!                  and G(0, tau) matrices. 
! KEYWORDS: Recompute Dynamic Green's function, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     By QRP for numerical stabilization, UR, UL are unitary matrices, but VR, VL are not.
!     By SVD for numerical stabilization, UR, UL, VR, VL are all unitary matrices.
!     This subroutine can be used for both cases of QRP and SVD decompositions for numerical stabilization, since
!           we only use the property that UR, UL are unitary matrices, during the calculations.
!
!     G(tau, tau) = [I + B(tau, 0)*B(beta, tau)]^{-1} = (1 + UR*DRVec*VR * VL*DLVec*UL)^{-1} 
!                 = UL^{-1} * [(UL*UR)^{-1} + DRVec*(VR*VL)*DLVec]^{-1} * UR^{-1}
!                 = UL^{-1} * [(UL*UR)^{-1} + DRMax*DRMin*(VR*VL)*DLMin*DLMax]^{-1} * UR^{-1}
!                 = UL^{-1}*DLMax^{-1} 
!                           * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin ]^{-1} 
!                           * DRMax^{-1}*UR^{-1}
!
!     G(tau, 0) = [I + B(tau, 0)*B(beta, tau)]^{-1} * B(tau, 0)
!               = (I + UR*DRVec*VR * VL*DLVec*UL)^{-1} * UR*DRVec*VR
!               = UL^{-1}*DLMax^{-1} 
!                           * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin ]^{-1} 
!                           * DRMin*VR
!
!     G(0, tau) = [B(tau, 0)]^{-1} * { I - [I + B(tau, 0)*B(beta, tau)]^{-1} }
!               = [I + B(beta, tau)*B(tau, 0)]^{-1} * B(beta, tau)
!               = B(beta, tau) * [I + B(tau, 0)*B(beta, tau)]^{-1}
!               = VL*DLVec*UL * (I + UR*DRVec*VR * VL*DLVec*UL)^{-1}
!               = VL*DLMin 
!                           * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin ]^{-1} 
!                           * DRMax^{-1}*UR^{-1}
!
!     Since there are some common terms in calculating these two dynamic single-particle Green's functions, we
!               perform the calculations in following steps.
!       (0) Calculate (UL*UR) and (VR*VL).
!       (1) Calculate DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} + DRMin*(VR*VL)*DLMin and
!                          UL^{-1}*DLMax^{-1}, (DRMin*VR), obtain G(tau, 0).
!       (2) Calculate (VL*DLMin) and (DRMax^{-1}*UR^{-1}), obtain G(0, tau).
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
      integer NT
      complex(rp) LogzDet
      real(rp) LogScaleRt(NmSpn), LogScaleLt(NmSpn)
      real(rp) UR   (NumNS, NumNS, NmSpn)  ! The UR matrix for B(tau, 0)
      real(rp) DRVec(NumNS       , NmSpn)  ! The DR vector for B(tau, 0) 
      real(rp) VR   (NumNS, NumNS, NmSpn)  ! The VR matrix for B(tau, 0)
      real(rp) VL   (NumNS, NumNS, NmSpn)  ! The VL matrix for B(\beta, tau)
      real(rp) DLVec(NumNS       , NmSpn)  ! The DL vector for B(\beta, tau) 
      real(rp) UL   (NumNS, NumNS, NmSpn)  ! The UL matrix for B(\beta, tau)
      real(rp) GrnFTT(NumNS, NumNS, NmSpn)  ! The output G(tau, tau) matrix
      real(rp) GrnFT0(NumNS, NumNS, NmSpn)  ! The output G(tau,   0) matrix
      real(rp) GrnF0T(NumNS, NumNS, NmSpn)  ! The output G(  0, tau) matrix
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
            integer I1, I2, SpnInd, Itp0
      real(rp) ReLogzDet, ImLogzDet
      complex(rp) LogzDetTmp
      real(rp), allocatable :: DRMinMat(:, :)
      real(rp), allocatable :: DRMaxInv(:, :)
      real(rp), allocatable :: DLMinMat(:, :)
      real(rp), allocatable :: DLMaxInv(:, :)
      real(rp), allocatable :: ULURTmpMat(:, :, :)
      real(rp), allocatable :: UVInvTmpMt(:, :, :)
      real(rp), allocatable :: ULInvTmpMt(:, :, :)
      real(rp), allocatable :: VLTmpMatrx(:, :, :)
      real(rp), allocatable :: URInvTmpMt(:, :, :)
      real(rp), allocatable :: VRTmpMatrx(:, :, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
            TimsGFDyn = TimsGFDyn + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of G(tau, 0) and G(0, tau) matrices ________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate some necessary matrices and initializations ______________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate necessary matrices for calculations _____________________
!________________________________________________________________________________________
      allocate(DRMinMat(NumNS, NmSpn)); DRMinMat = 0.0_rp
      allocate(DRMaxInv(NumNS, NmSpn)); DRMaxInv = 0.0_rp
      allocate(DLMinMat(NumNS, NmSpn)); DLMinMat = 0.0_rp
      allocate(DLMaxInv(NumNS, NmSpn)); DLMaxInv = 0.0_rp
      allocate(ULURTmpMat(NumNS, NumNS, NmSpn)); ULURTmpMat = 0.0_rp
      allocate(UVInvTmpMt(NumNS, NumNS, NmSpn)); UVInvTmpMt = 0.0_rp
      allocate(ULInvTmpMt(NumNS, NumNS, NmSpn)); ULInvTmpMt = 0.0_rp
      allocate(VLTmpMatrx(NumNS, NumNS, NmSpn)); VLTmpMatrx = 0.0_rp
      allocate(URInvTmpMt(NumNS, NumNS, NmSpn)); URInvTmpMt = 0.0_rp
            allocate(VRTmpMatrx(NumNS, NumNS, NmSpn)); VRTmpMatrx = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) Split DRVec = DRMax*DRMin, DLVec = DLMin*DLMax ___________________
!_____________________ Calculate Logdet of DRMax*DLMax matrix ___________________________
!________________________________________________________________________________________
      LogzDet = rp_Zzero
      call SplitDMat_RtLtD(NT, DRVec, DLVec, LogScaleRt, LogScaleLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogzDet)
!**************************************************************************************************     
!___________________ 1. Calculate all the pieces for Green's Function matrices ____________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Calculate ULURTmpMat = UL * UR ___________________________________
!_______________________________ UVInvTmpMt = VR * VL ___________________________________
!________________________________________________________________________________________
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, UL, UR, 0.0_rp, ULURTmpMat)
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, VR, VL, 0.0_rp, UVInvTmpMt)
!________________________________________________________________________________________         
!_________________ (1) UVInvTmpMt = DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} __________________
!_____________________                + DRMin*(VR*VL)*DLMin _____________________________
!_____________________ ULInvTmpMt = UL^{-1} * DLMax^{-1} ________________________________
!_____________________ VLTmpMatrx = VL * DLMin __________________________________________
!_____________________ URInvTmpMt = DRMax^{-1} * UR^{-1} ________________________________
!_____________________ VRTmpMatrx = DRMin * VR __________________________________________
!________________________________________________________________________________________         
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               UVInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * ULURTmpMat(I2, I1, SpnInd) *DLMaxInv(I2, SpnInd) &
                     & + DRMinMat(I1, SpnInd) * UVInvTmpMt(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
               ULInvTmpMt(I1, I2, SpnInd) = UL(I2, I1, SpnInd) * DLMaxInv(I2, SpnInd)
               VLTmpMatrx(I1, I2, SpnInd) = VL(I1, I2, SpnInd) * DLMinMat(I2, SpnInd)
               URInvTmpMt(I1, I2, SpnInd) = DRMaxInv(I1, SpnInd) * UR(I2, I1, SpnInd)
               VRTmpMatrx(I1, I2, SpnInd) = DRMinMat(I1, SpnInd) * VR(I1, I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________         
!_________________ (2) Calculate the log determinant log[det(UL*UR)] ____________________
!________________________________________________________________________________________ 
      LogzDetTmp = rp_Zzero
      call dMatDet_LogDet_QMC(NumNS, ULURTmpMat, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!**************************************************************************************************     
!___________________ 2. Calculate G(tau, tau) and G(0, tau) matrices ______________________________
!________ G(tau, tau) = UL^{-1}*DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} _________________
!______________________      + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1}*UR^{-1} ____________________
!________ G(  0, tau) = VL*DLMin * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} ___________________________
!______________________      + DRMin*(VR*VL)*DLMin ]^{-1} * DRMax^{-1}*UR^{-1} ____________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) URInvTmpMt = UVInvTmpMt^{-1} * URInvTmpMt ________________________
!________________________________________________________________________________________
      ULURTmpMat = 0.0_rp
      call dMat_Copy_QMC(UVInvTmpMt(1, 1, 1), ULURTmpMat(1, 1, 1))
      LogzDetTmp = rp_Zzero
      call dMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, ULURTmpMat, URInvTmpMt, LogzDetTmp)
      LogzDet = LogzDet + LogzDetTmp
!________________________________________________________________________________________         
!_________________ (1) GrnFTT = ULInvTmpMt * UVInvTmpMt^{-1} * URInvTmpMt ________________
!________________________________________________________________________________________
      GrnFTT = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, ULInvTmpMt, URInvTmpMt, 0.0_rp, GrnFTT)
!________________________________________________________________________________________         
!_________________ (2) GrnF0T = VLTmpMatrx * UVInvTmpMt^{-1} * URInvTmpMt ________________
!________________________________________________________________________________________
      GrnF0T = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, VLTmpMatrx, URInvTmpMt, 0.0_rp, GrnF0T)
!**************************************************************************************************     
!___________________ 3. Calculate G(tau, 0) matrice _______________________________________________
!_____________ G(tau, 0) = UL^{-1}*DLMax^{-1} * [ DRMax^{-1}*(UL*UR)^{-1}*DLMax^{-1} ______________
!_________________________     + DRMin*(VR*VL)*DLMin ]^{-1} * DRMin*VR ____________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) VRTmpMatrx = UVInvTmpMt^{-1} * VRTmpMatrx ________________________
!________________________________________________________________________________________
      LogzDetTmp = rp_Zzero
      call dMatEqSet_Rght_LogDet_QMC(NumNS, NumNS, UVInvTmpMt, VRTmpMatrx, LogzDetTmp)
!________________________________________________________________________________________         
!_________________ (1) GrnFT0 = ULInvTmpMt * UVInvTmpMt^{-1} * VRTmpMatrx ________________
!________________________________________________________________________________________
      GrnFT0 = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, ULInvTmpMt, VRTmpMatrx, 0.0_rp, GrnFT0)
!**************************************************************************************************     
!___________________ 4. Recalculate the log of the determinant ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Reset the imaginary part of LogzDet ______________________________
!________________________________________________________________________________________ 
      ReLogzDet = real(LogzDet); ImLogzDet = aimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = cmplx(ReLogzDet, ImLogzDet, rp)
!**************************************************************************************************     
!___________________ 5. Deallocate all the used matrices in this subroutine _______________________
!**************************************************************************************************
      if(allocated(DLMinMat  )) deallocate(DLMinMat  )
      if(allocated(DLMaxInv  )) deallocate(DLMaxInv  )
      if(allocated(DRMinMat  )) deallocate(DRMinMat  )
      if(allocated(DRMaxInv  )) deallocate(DRMaxInv  )
      if(allocated(ULURTmpMat)) deallocate(ULURTmpMat)
      if(allocated(UVInvTmpMt)) deallocate(UVInvTmpMt)
            if(allocated(ULInvTmpMt)) deallocate(ULInvTmpMt)
      if(allocated(VLTmpMatrx)) deallocate(VLTmpMatrx)
      if(allocated(URInvTmpMt)) deallocate(URInvTmpMt)
      if(allocated(VRTmpMatrx)) deallocate(VRTmpMatrx)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed in Calculating Static GrF Matrices ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
      TimeGFDyn = TimeGFDyn + TimeIntrvl(time1, time2)
      
   end subroutine GrFDynamcR_LR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
!######################################## Split diagonal D matrix D = DMax*DMin #########################################
!######################################## Split diagonal D matrix D = DMax*DMin #########################################
!######################################## Split diagonal D matrix D = DMax*DMin #########################################
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine SplitDMat_SgleD(NT, DVec, LogScale, DMinMat, DMaxInv, LogDetDMax)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SplitDMat_SgleD(NT, DVec, LogScale, DMinMat, DMaxInv, LogDetDMax)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to split the DVec matrix as DVec = DVecMax*DVecMin, and obtain inverse of
!                 of DVecMax and the log of its determinant, in a numerically stable way.
! KEYWORDS: Split DVec as DVec = DVecMax*DVecMin.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: The procedure used as follows is to ensure that the largest and smallest number in DVec do not surpass
!                 the range of real(rp) formulation.
!
!     Input:  NT       --> Time slice index;
!             DVec     --> The input diagonal matrix;
!             LogScale --> The input log of scals used for DVec;
!
!     Output: DMinMat    --> The output of DVecMin matrix (diagonal elements);
!             DMaxInv    --> The output inverse of DVecMax matrix (diagonal elements);
!             LogDetDMax --> The log of determinant of DVecMax matrix.
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
      integer NT
      real(rp) DVec(NumNS, NmSpn)
      real(rp) LogScale(NmSpn)
      real(rp) DMinMat(NumNS, NmSpn)
      real(rp) DMaxInv(NumNS, NmSpn)
      complex(rp) LogDetDMax
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, SpnInd
      real(rp) LogDiagVal
!______________________________________________________________________________________________________________     
!___________________________ Calculate the splitting DVec = DVecMax * DVecMin _________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Split DVec matrix as  DVec = DVecMax * DVecMin ____________________________
!______________________ And calculate logdet of DVevMax matrix ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Initialization of output quantities ______________________________
!________________________________________________________________________________________
      DMinMat = 0.0_rp
      DMaxInv = 0.0_rp
      LogDetDMax = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Choose DVec and MuCoefft according to LRIndex ____________________
!________________________________________________________________________________________
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, NumNS, +1
            LogDiagVal = log(DVec(I1, SpnInd)) + LogScale(SpnInd) - NT*Dltau*ChemP_Sgle
            if( LogDiagVal > 0.0_rp ) then
               DMinMat(I1, SpnInd) = 1.0_rp
               if( LogDiagVal >= +GrFTrcDhld*log(10.0_rp) ) then
                  DMaxInv(I1, SpnInd) = 0.0_rp
               else
                  DMaxInv(I1, SpnInd) = exp(-LogDiagVal)
               end if
               LogDetDMax = LogDetDMax + cmplx(LogDiagVal, 0.0_rp, rp)
            else
               if( LogDiagVal <= -GrFTrcDhld*log(10.0_rp) ) then
                  DMinMat(I1, SpnInd) = 0.0_rp
               else
                  DMinMat(I1, SpnInd) = exp(+LogDiagVal)
               end if
               DMaxInv(I1, SpnInd) = 1.0_rp
            end if
         enddo
      enddo
      
   end subroutine SplitDMat_SgleD
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine SplitDMat_RtLtD(NT, DRVec, DLVec, LogDRt, LogDLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogDetDMax)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SplitDMat_RtLtD(NT, DRVec, DLVec, LogDRt, LogDLt, DRMinMat, DRMaxInv, DLMinMat, DLMaxInv, LogDetDMax)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to split the DVec matrix as DVec = DVecMax*DVecMin, and obtain inverse of
!                 of DVecMax and the log of its determinant, in a numerically stable way.
! KEYWORDS: Split DRVec, DLVec as DRVec = DRVecMax*DRVecMin, DLVec = DLVecMin*DLVecMax.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: The procedure used as follows is to ensure that the largest and smallest number in DVec do not surpass
!                  the range of real(rp) formulation.
!
!     Input:  NT     --> Time slice index;
!             DRVec  --> The input diagonal matrix for right side;
!             DLVec  --> The input diagonal matrix for left side;
!             LogDRt --> The input log of scals used for DRVec;
!             LogDLt --> The input log of scals used for DLVec;
!
!     Output: DRMinMat   --> The output of DRMin matrix (diagonal elements);
!             DRMaxInv   --> The output inverse of DRMax matrix (diagonal elements);
!             DLMinMat   --> The output of DLMin matrix (diagonal elements);
!             DLMaxInv   --> The output inverse of DLMax matrix (diagonal elements);
!             LogDetDMax --> The log of determinant of DRMax*DLMax matrix.
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
      integer NT
      real(rp) DRVec(NumNS, NmSpn)
      real(rp) DLVec(NumNS, NmSpn)
      real(rp) LogDRt(NmSpn), LogDLt(NmSpn)
      real(rp) DRMinMat(NumNS, NmSpn)
      real(rp) DRMaxInv(NumNS, NmSpn)
      real(rp) DLMinMat(NumNS, NmSpn)
      real(rp) DLMaxInv(NumNS, NmSpn)
      complex(rp) LogDetDMax
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, SpnInd
      real(rp) LogDiagVal
!______________________________________________________________________________________________________________     
!___________________________ Calculate the splitting DVec = DVecMax * DVecMin _________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Split DVec matrix as  DVec = DVecMax * DVecMin ____________________________
!______________________ And calculate logdet of DVevMax matrix ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Initialization of output quantities ______________________________
!________________________________________________________________________________________
      DRMinMat = 0.0_rp; DRMaxInv = 0.0_rp
      DLMinMat = 0.0_rp; DLMaxInv = 0.0_rp
      LogDetDMax = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Choose DVec and MuCoefft according to LRIndex ____________________
!________________________________________________________________________________________
      do SpnInd = 1, NmSpn, +1
         do I1 = 1, NumNS, +1
!____________________________________________________________________________         
!________________ [0] For DRVec matrix of the right side ____________________
!____________________________________________________________________________
            LogDiagVal = log(DRVec(I1, SpnInd)) + LogDRt(SpnInd) - NT*Dltau*ChemP_Rght
            if( LogDiagVal > 0.0_rp ) then
               DRMinMat(I1, SpnInd) = 1.0_rp
               if( LogDiagVal >= +GrFTrcDhld*log(10.0_rp) ) then
                  DRMaxInv(I1, SpnInd) = 0.0_rp
               else
                  DRMaxInv(I1, SpnInd) = exp(-LogDiagVal)
               end if
               LogDetDMax = LogDetDMax + cmplx(LogDiagVal, 0.0_rp, rp)
            else
               if( LogDiagVal <= -GrFTrcDhld*log(10.0_rp) ) then
                  DRMinMat(I1, SpnInd) = 0.0_rp
               else
                  DRMinMat(I1, SpnInd) = exp(+LogDiagVal)
               end if
               DRMaxInv(I1, SpnInd) = 1.0_rp
            end if
!____________________________________________________________________________         
!________________ [1] For DLVec matrix of the Left side _____________________
!____________________________________________________________________________
            LogDiagVal = log(DLVec(I1, SpnInd)) + LogDLt(SpnInd) - (LTrot-NT)*Dltau*ChemP_Left
            if( LogDiagVal > 0.0_rp ) then
               DLMinMat(I1, SpnInd) = 1.0_rp
               if( LogDiagVal >= +GrFTrcDhld*log(10.0_rp) ) then
                  DLMaxInv(I1, SpnInd) = 0.0_rp
               else
                  DLMaxInv(I1, SpnInd) = exp(-LogDiagVal)
               end if
               LogDetDMax = LogDetDMax + cmplx(LogDiagVal, 0.0_rp, rp)
            else
               if( LogDiagVal <= -GrFTrcDhld*log(10.0_rp) ) then
                  DLMinMat(I1, SpnInd) = 0.0_rp
               else
                  DLMinMat(I1, SpnInd) = exp(+LogDiagVal)
               end if
               DLMaxInv(I1, SpnInd) = 1.0_rp
            end if
         enddo
      enddo
      
   end subroutine SplitDMat_RtLtD
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
!################################### GreenF with Symmetric Trotter Decomposition ########################################
!################################### GreenF with Symmetric Trotter Decomposition ########################################
!################################### GreenF with Symmetric Trotter Decomposition ########################################
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrnFSymTrotR(H0orHT, GreenF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrnFSymTrotR(H0orHT, GreenF)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to obtain the GreenF matrix within symmetric Trotter decomposition. 
! KEYWORDS: Only measure the total electron density.
! AUTHOR:   Yuan-Yao He
! TIME:     2018-11-07
! DESCRIPTION:
!
!     Calculate GreenF --> exp(-dt*K/2) * GreenF * [exp(-dt*K/2)]^{-1}. Real version.
!
!     Input:  H0orHT --> =="HT" --> exp(-dt*H_T/2) and =="H0" --> exp(-dt*H_0/2);
!             GreenF --> Green's Function matrice <c_i * c_j^+> within asymmetric Trotter decomposition.
!     
!     Output: GreenF --> Green's Function matrice <c_i * c_j^+> within  Symmetric Trotter decomposition.
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
      character(2) H0orHT
      real(rp) GreenF(NumNS, NumNS, NmSpn)
!______________________________________________________________________________________________________________     
!__________________________ Obtain GreenF within symmetric Trotter decomposition ______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. For H0orHT == "H0" case, use exp(-dt*H_0/2) matrix ________________________
!**************************************************************************************************
      if(H0orHT == "H0") then
!________________________________________________________________________________________         
!_________________ (0) For FFTEXPDTH0 == T case --> Use FFT method ______________________
!________________________________________________________________________________________ 
         if(FFTEXPDTH0) then
            call ExpdtKMultFFTEgVl_R(ExpdtEofH0(1, 1, 3), GreenF(1, 1, 1), NumNS, NumNS)
            call FFTEgVlMultExpdtK_R(NumNS, NumNS, GreenF(1, 1, 1), ExpdtEofH0(1, 1, 4))
!________________________________________________________________________________________         
!_________________ (1) For FFTEXPDTH0 == F case --> call DGEMM subroutine _______________
!________________________________________________________________________________________ 
         else
            call ExpdtKMultMatFull_R(ExpdtH0Mat(1, 1, 1, 3), GreenF(1, 1, 1), NumNS, NumNS)
            call MatMultExpdtKFull_R(NumNS, NumNS, GreenF(1, 1, 1), ExpdtH0Mat(1, 1, 1, 4))
         end if
!**************************************************************************************************     
!___________________ 1. For H0orHT == "HT" case, use exp(-dt*H_T/2) matrix ________________________
!**************************************************************************************************
      else if(H0orHT == "HT") then
!________________________________________________________________________________________         
!_________________ (0) For FFTEXPDTHT == T case --> Use FFT method ______________________
!________________________________________________________________________________________ 
         if(FFTEXPDTHT) then
            call ExpdtKMultFFTEgVl_R(ExpdtOfHTe(1, 1, 3), GreenF(1, 1, 1), NumNS, NumNS)
            call FFTEgVlMultExpdtK_R(NumNS, NumNS, GreenF(1, 1, 1), ExpdtOfHTe(1, 1, 4))
!________________________________________________________________________________________         
!_________________ (1) For FFTEXPDTHT == F case --> call DGEMM subroutine _______________
!________________________________________________________________________________________ 
         else
            call ExpdtKMultMatFull_R(ExpdtTryHT(1, 1, 1, 3), GreenF(1, 1, 1), NumNS, NumNS)
            call MatMultExpdtKFull_R(NumNS, NumNS, GreenF(1, 1, 1), ExpdtTryHT(1, 1, 1, 4))
         end if         
      end if
            
   end subroutine GrnFSymTrotR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine GrnFSymTrotC(H0orHT, GreenF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GrnFSymTrotC(H0orHT, GreenF)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to obtain the GreenF matrix within symmetric Trotter decomposition. 
! KEYWORDS: Only measure the total electron density.
! AUTHOR:   Yuan-Yao He
! TIME:     2018-11-07
! DESCRIPTION:
!
!     Input:  H0orHT --> =="HT" --> exp(-dt*H_T/2) and =="H0" --> exp(-dt*H_0/2);
!             GreenF --> Green's Function matrice <c_i * c_j^+> within asymmetric Trotter decomposition.
!     
!     Output: GreenF --> Green's Function matrice <c_i * c_j^+> within  Symmetric Trotter decomposition.
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
      character(2) H0orHT
      complex(rp) GreenF(NumNS, NumNS, NmSpn)
!______________________________________________________________________________________________________________     
!__________________________ Obtain GreenF within symmetric Trotter decomposition ______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. For H0orHT == "H0" case, use exp(-dt*H_0/2) matrix ________________________
!**************************************************************************************************
      if(H0orHT == "H0") then
!________________________________________________________________________________________         
!_________________ (0) For FFTEXPDTH0 == T case --> Use FFT method ______________________
!________________________________________________________________________________________ 
         if(FFTEXPDTH0) then
            call ExpdtKMultFFTEgVl_C(ExpdtEofH0(1, 1, 3), GreenF(1, 1, 1), NumNS, NumNS)
            call FFTEgVlMultExpdtK_C(NumNS, NumNS, GreenF(1, 1, 1), ExpdtEofH0(1, 1, 4))
!________________________________________________________________________________________         
!_________________ (1) For FFTEXPDTH0 == F case --> call ZGEMM subroutine _______________
!________________________________________________________________________________________
         else
            call ExpdtKMultMatFull_C(ExpdtH0Mat(1, 1, 1, 3), GreenF(1, 1, 1), NumNS, NumNS)
            call MatMultExpdtKFull_C(NumNS, NumNS, GreenF(1, 1, 1), ExpdtH0Mat(1, 1, 1, 4))
         end if
!**************************************************************************************************     
!___________________ 1. For H0orHT == "HT" case, use exp(-dt*H_T/2) matrix ________________________
!**************************************************************************************************
      else if(H0orHT == "HT") then
!________________________________________________________________________________________         
!_________________ (0) For FFTEXPDTHT == T case --> Use FFT method ______________________
!________________________________________________________________________________________
         if(FFTEXPDTHT) then
            call ExpdtKMultFFTEgVl_C(ExpdtOfHTe(1, 1, 3), GreenF(1, 1, 1), NumNS, NumNS)
            call FFTEgVlMultExpdtK_C(NumNS, NumNS, GreenF(1, 1, 1), ExpdtOfHTe(1, 1, 4))
!________________________________________________________________________________________         
!_________________ (1) For FFTEXPDTHT == F case --> call ZGEMM subroutine _______________
!________________________________________________________________________________________ 
         else
            call ExpdtKMultMatFull_C(ExpdtTryHT(1, 1, 1, 3), GreenF(1, 1, 1), NumNS, NumNS)
            call MatMultExpdtKFull_C(NumNS, NumNS, GreenF(1, 1, 1), ExpdtTryHT(1, 1, 1, 4))
         end if
      end if
            
  end subroutine GrnFSymTrotC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$