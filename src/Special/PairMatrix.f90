!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to compute the n(k) and pairing matrix using the k-space single-particle Green's
!              function matrix.
! COMMENT: Measure n(k) and Pairing matrix.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   NkPairWvfc_H0FFTW --> Subroutine used to the n(k) and pairing matrix using G(k, k');
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine NkPairWvfc_H0FFTW(CfgConst, GFk_All, GFk_AllC, NkDistr, PairMat) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  NkPairWvfc_H0FFTW(CfgConst, GFk_All, GFk_AllC, NkDistr, PairMat) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to the n(k) and pairing matrix using G(k, k').
! KEYWORDS: Measure n(k) and Pairing matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Measure the pairing matrix M(k, s; k', s') = < \Delta_s(k)^+ \Delta_s'(k') >.
!
!         For spin-singlet pairing: 
!                    \Delta_s(k)^+ = 1/sqrt(2)*(c_{k\up}^+ c_{-k\down}^+ - c_{-k\down}^+ c_{k\up}^+)
!         For spin-triplet pairing: 
!                    \Delta_{t,u}(k)^+ =            c_{k\up}^+ c_{-k\up}^+
!                    \Delta_{t,0}(k)^+ = 1/sqrt(2)*(c_{k\up}^+ c_{-k\down}^+ + c_{-k\down}^+ c_{k\up}^+)
!                    \Delta_{t,d}(k)^+ =            c_{k\down}^+ c_{-k\down}^+
!
!     For the symmetry reason, we only study the \Delta_s(k)^+, \Delta_{t,u}(k)^+, \Delta_{t,d}(k)^+ channels.
!
!     Input:  CfgConst --> The configuration/walker related constant (like sign/phase or walker weight);
!             GFk_All  --> Single-particle Green's function in k space as GFk_All (k, q) = <c_k c_q^+>;
!             GFk_AllC --> Single-particle Green's function in k space as GFk_AllC(k, q) = <c_k^+ c_q>;
!             NkDistr  --> Accumulated measuring results of n(k) as momentum distribution;
!             PairMat  --> Accumulated measuring results of the pairing matrix;
!     
!     Output: NkDistr --> Accumulated measuring results of n(k) as momentum distribution;
!             PairMat --> Accumulated measuring results of the pairing matrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		use RealPrecsn
		use CoreParamt
		use Observable
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) CfgConst
      complex(rp) GFk_All(NumNS, NumNS, 2), GFk_AllC(NumNS, NumNS, 2)
      complex(rp) NkDistr(NumNC, 2), PairMat(NumNC, NumNC)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, I3, I4, Ix, Iy
      complex(rp), allocatable :: NkDistrHere(:, :)
      complex(rp), allocatable :: PairMatHere(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Initializations for this subroutine _______________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Allocate temporary matrices for measurements _____________________
!________________________________________________________________________________________
      allocate(NkDistrHere(NumNC, 2))
      allocate(PairMatHere(NumNC, NumNC))
      NkDistrHere = rp_Zzero
      PairMatHere = rp_Zzero
!**************************************************************************************************	  
!___________________ 1. Measure momentum distribution n(k) using G(k, k') _________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Measure the Momentum Distribution ________________________________
!________________________________________________________________________________________
      do I0 = 1, NumNC, +1
         NkDistrHere(I0, 1) = GFk_AllC(I0, I0, 1)
         NkDistrHere(I0, 2) = GFk_AllC(I0, I0, 2)
      enddo
!________________________________________________________________________________________ 	  
!_________________ (1) Accumulate the results and consider the sign _____________________
!________________________________________________________________________________________
      NkDistr = NkDistr + NkDistrHere * CfgConst
!**************************************************************************************************	  
!___________________ 2. Measure the pairing matrix of the system __________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Measure the pairing matrix using GFk_AllC matrix _________________
!________________________________________________________________________________________
      do I1 = 1, NumNC, +1
         Ix = - mod(I1-1, NumL1); Iy = - (I1-1)/NumL1
         if(Ix < 0) Ix = Ix + NumL1
         if(Iy < 0) Iy = Iy + NumL2
         I3 = Iy*NumL1 + Ix + 1
         do I2 = 1, NumNC, +1
            Ix = - mod(I2-1, NumL1); Iy = - (I2-1)/NumL1
            if(Ix < 0) Ix = Ix + NumL1
            if(Iy < 0) Iy = Iy + NumL2
            I4 = Iy*NumL1 + Ix + 1
            if(HubbU <= 0.0_rp) then
               PairMatHere(I1, I2) = GFk_AllC(I1, I2, 1) * GFk_AllC(I3, I4, 2)
            else
               PairMatHere(I1, I2) =   ( GFk_AllC(I1, I2, 1) * GFk_AllC(I3, I4, 2) &
                                     & + GFk_AllC(I1, I2, 2) * GFk_AllC(I3, I4, 1) &
                                     & + GFk_AllC(I1, I4, 1) * GFk_AllC(I3, I2, 2) &
                                     & + GFk_AllC(I1, I4, 2) * GFk_AllC(I3, I2, 1) ) / 2.0_rp
            end if
         enddo
      enddo
!________________________________________________________________________________________ 	  
!_________________ (1) Accumulate the results and consider the sign _____________________
!________________________________________________________________________________________
      PairMat = PairMat + PairMatHere * CfgConst
!**************************************************************************************************	  
!___________________ 3. Finalization for this measurement _________________________________________
!**************************************************************************************************
      if(allocated(NkDistrHere)) deallocate(NkDistrHere) 
      if(allocated(PairMatHere)) deallocate(PairMatHere) 
      
   end subroutine NkPairWvfc_H0FFTW
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$