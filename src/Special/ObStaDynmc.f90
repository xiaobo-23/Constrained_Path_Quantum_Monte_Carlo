!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the dynamic measurements of single-particle Green's function and 
!                correlation functions in the CPMC simulations. For both periodic and open boundary conditions.
! COMMENT: Measurements of real-space correlations.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   ObDynCrFct_SpnDcp --> Subroutine used to perform the dynamic measurements for both PERIODIC and OPEN BCs.
!               
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine ObDynCrFct_SpnDcp(NTInd, CfgConst, RealSpCrFTau) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ObDynCrFct_SpnDcp(NTInd, CfgConst, RealSpCrFTau) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to measure the time-displaced physical quantities.
! KEYWORDS: Calculate the Dynamic physical quantities.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Measure the time-displaced physical quantities.
!
!                 GT0_{i,j} = + <c_i(\tau)c_j^+>
!                 G0T_{i,j} = + <c_j^+(\tau)c_i>
!                 G00_{i,j} = + <c_ic_j^+>
!                 GTT_{i,j} = + <c_i(\tau)c_j^+(\tau)>
!
!     For the index in RealSpCrFTau(:, :, 01:40): 
!                                      01~04 for single-particle Green's functions;
!                                      05~19 for regular correlations; 
!                                      20~34 for vertex contributions of the pairing correlations;
!                                      35~40 for computing the background for correlations.
!
!     Input: NTInd    --> Imaginary-time \tau in analytic expression;  
!            CfgConst --> The configuration/walker related constant (like sign/phase or walker weight);       
!
!     Outpt: RealSpCrFTau --> The accumulated results.
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
		use Observable
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NTInd
      real(rp) CfgConst
      real(rp) RealSpCrFTau(0:NumTauPnt, NmSitePairTau, 40)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      integer I0, I1, I2, I3, I4, Idimj, Kd, Hd, I0ItrMax
      real(rp) RlCrFTmp(15)
      real(rp), allocatable :: RSpCrFTauHere(:, :)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for dynamic measurements __________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      TimsDyMea = TimsDyMea + 1
      call system_clock(time1)
!______________________________________________________________________________________________________________	  
!_____________________________ Main calculations for all the dynamic quantities _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Initializations for this subroutine _______________________________________
!************************************************************************************************** 
!________________________________________________________________________________________ 	  
!_________________ (0) Temporary matrix for r-space correlation functions _______________
!________________________________________________________________________________________
      allocate(RSpCrFTauHere(NmSitePairTau, 40)); RSpCrFTauHere = 0.0_rp
!**************************************************************************************************	  
!___________________ 1. Calculate all the dynamic correlation functions ___________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Measure all the dynamic r-space correlation functions ____________
!________________________________________________________________________________________
      I0ItrMax = merge(NumNC*NumNC, NmSitePairTau, IfPyObsPBC)
   !$OMP PARALLEL &
   !$OMP PRIVATE(I0, I1, I2, I3, I4, Kd, Hd, Idimj, RlCrFTmp)
   !$OMP DO REDUCTION(+ : RSpCrFTauHere)
		do I0 = 1, I0ItrMax, +1
!____________________________________________________________________________ 	  
!________________ The integer index for the present term ____________________
!____________________________________________________________________________
         if(IfPyObsPBC) then
            I1 = (I0-1)/NumNC + 1; I2 = mod(I0-1, NumNC) + 1
         else
            I1 = TauBondPairId(I0, 1); I2 = TauBondPairId(I0, 2)
         end if
         Idimj = merge(IminusJ(I1, I2), I0, IfPyObsPBC)
!____________________________________________________________________________ 	  
!________________ Spin-up and down single-particle Green's functions ________
!________________________ with both Tau>0 and Tau<0 _________________________
!____________________________________________________________________________
         !!!!!!!!!! Spin-up channel
         RlCrFTmp(1) = + GrFT0(I1, I2, 1)
         RlCrFTmp(2) = + GrF0T(I1, I2, 1)
         RSpCrFTauHere(Idimj, 01) = RSpCrFTauHere(Idimj, 01) + RlCrFTmp(1)
         RSpCrFTauHere(Idimj, 02) = RSpCrFTauHere(Idimj, 02) + RlCrFTmp(2)
         !!!!!!!!!! Spin-down channel
         RlCrFTmp(1) = + GrFT0(I1, I2, 2)
         RlCrFTmp(2) = + GrF0T(I1, I2, 2)
         RSpCrFTauHere(Idimj, 03) = RSpCrFTauHere(Idimj, 03) + RlCrFTmp(1)
         RSpCrFTauHere(Idimj, 04) = RSpCrFTauHere(Idimj, 04) + RlCrFTmp(2)
!____________________________________________________________________________ 	  
!________________ Spin correlations, <SzSz> and <S+S- + S-S+>/2 _____________
!____________________________________________________________________________
         !!!!!!!!!! The <SzSz> correlation
         RlCrFTmp(1) =     ( GrFTT(I1, I1, 1) - GrFTT(I1, I1, 2) ) * ( GrF00(I2, I2, 1) - GrF00(I2, I2, 2) ) &
                       & +   GrF0T(I2, I1, 1) * GrFT0(I1, I2, 1) + GrF0T(I2, I1, 2) * GrFT0(I1, I2, 2)
         RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
         RSpCrFTauHere(Idimj, 05) = RSpCrFTauHere(Idimj, 05) + RlCrFTmp(1)
         !!!!!!!!!! The <S+S- + S-S+>/2 correlation
         RlCrFTmp(2) = GrF0T(I2, I1, 1) * GrFT0(I1, I2, 2) + GrF0T(I2, I1, 2) * GrFT0(I1, I2, 1)
         RlCrFTmp(2) = RlCrFTmp(2) / 2.0_rp
         RSpCrFTauHere(Idimj, 06) = RSpCrFTauHere(Idimj, 06) + RlCrFTmp(2)
!____________________________________________________________________________ 	  
!________________ Density-density correlation function ______________________
!____________________________________________________________________________
         RlCrFTmp(1) =   + GrF0T(I2, I1, 1) * GrFT0(I1, I2, 1) + GrF0T(I2, I1, 2) * GrFT0(I1, I2, 2) &
                       & + ( 2.0_rp - GrFTT(I1, I1, 1) - GrFTT(I1, I1, 2) ) * & 
                       &   ( 2.0_rp - GrF00(I2, I2, 1) - GrF00(I2, I2, 2) )
         RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
         RSpCrFTauHere(Idimj, 07) = RSpCrFTauHere(Idimj, 07) + RlCrFTmp(1)
!____________________________________________________________________________ 	  
!________________ Current-Current correlations in x-direction _______________
!____________________________________________________________________________ 
         I3 = FNNBond(I1, 2); I4 = FNNBond(I2, 2)
         RlCrFTmp(1) =   ( - GrFTT(I3, I1, 1) + GrFTT(I1, I3, 1) - GrFTT(I3, I1, 2) + GrFTT(I1, I3, 2) ) * &
                       & ( - GrF00(I4, I2, 1) + GrF00(I2, I4, 1) - GrF00(I4, I2, 2) + GrF00(I2, I4, 2) )
         RlCrFTmp(1) = RlCrFTmp(1) + GrF0T(I4, I1, 1)*GrFT0(I3, I2, 1) - GrF0T(I4, I3, 1)*GrFT0(I1, I2, 1) &
                                 & - GrF0T(I2, I1, 1)*GrFT0(I3, I4, 1) + GrF0T(I2, I3, 1)*GrFT0(I1, I4, 1) &
                                 & + GrF0T(I4, I1, 2)*GrFT0(I3, I2, 2) - GrF0T(I4, I3, 2)*GrFT0(I1, I2, 2) &
                                 & - GrF0T(I2, I1, 2)*GrFT0(I3, I4, 2) + GrF0T(I2, I3, 2)*GrFT0(I1, I4, 2)
         RSpCrFTauHere(Idimj, 08) = RSpCrFTauHere(Idimj, 08) - RlCrFTmp(1)
!____________________________________________________________________________ 	  
!________________ On-site spin-singlet s-wave pairing correlation ___________
!____________________________________________________________________________   
         RlCrFTmp(1) = GrFT0(I1, I2, 1) * GrFT0(I1, I2, 2) + GrF0T(I2, I1, 1) * GrF0T(I2, I1, 2)
         RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
         RSpCrFTauHere(Idimj, 09) = RSpCrFTauHere(Idimj, 09) + RlCrFTmp(1)
!____________________________________________________________________________ 	  
!________________ Various off-site pairing correlation functions ____________
!____________________________________________________________________________
         !!!!!!!!!! Initialization for the accumulation
         RlCrFTmp = 0.0_rp
         !!!!!!!!!! Summation over all four nearest neighbors
         do Kd = 1, 4, +1
            I3 = FNNBond(I1, Kd)
            do Hd = 1, 4, +1
               I4 = FNNBond(I2, Hd)
               !!!!!!!! The spin-singlet pairing channels
               RlCrFTmp(15) =   + GrF0T(I2, I1, 1)*GrF0T(I4, I3, 2) + GrF0T(I2, I3, 1)*GrF0T(I4, I1, 2) &
                              & + GrF0T(I4, I1, 1)*GrF0T(I2, I3, 2) + GrF0T(I4, I3, 1)*GrF0T(I2, I1, 2)
               !!!!!!!! Accumulate the pairing results
               !!!!!! For the NN extended s-wave pairing
               RlCrFTmp(1) = RlCrFTmp(1) + RlCrFTmp(15)
               !!!!!! For the NN d-wave pairing
               RlCrFTmp(2) = RlCrFTmp(2) + dble(1-2*mod(Kd, 2)) * dble(1-2*mod(Hd, 2)) * RlCrFTmp(15)
            enddo
         enddo
         !!!!!!!!!! For the extended s-wave pairing
         RlCrFTmp(1) = RlCrFTmp(1) / 32.0_rp
         RSpCrFTauHere(Idimj, 10) = RSpCrFTauHere(Idimj, 10) + RlCrFTmp(1)
         !!!!!!!!!! For the d-wave pairing
         RlCrFTmp(2) = RlCrFTmp(2) / 32.0_rp
         RSpCrFTauHere(Idimj, 11) = RSpCrFTauHere(Idimj, 11) + RlCrFTmp(2)
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%	  
!________________ Recording the background for the correlations _____________
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
         !!!!!!!!!! For IfPyObsPBC==F case, background for all correlations; For IfPyObsPBC==T case, only 1~NumNC
         if(IfPyObsPBC .and. I0 <= NumNC) then
            !!!!!!!! For Sz-Sz correlation functions
            RSpCrFTauHere(I2, 35) = ( GrFTT(I2, I2, 2) - GrFTT(I2, I2, 1) ) / 2.0_rp            ! S_{I2}^z(\tau)
            RSpCrFTauHere(I2, 36) = ( GrF00(I2, I2, 2) - GrF00(I2, I2, 1) ) / 2.0_rp            ! S_{I2}^z(0)
            !!!!!!!! For density-density correlation functions
            RSpCrFTauHere(I2, 37) = ( 2.0_rp - GrFTT(I2, I2, 1) - GrFTT(I2, I2, 2) ) / 2.0_rp   ! n_{I2}(\tau)
            RSpCrFTauHere(I2, 38) = ( 2.0_rp - GrF00(I2, I2, 1) - GrF00(I2, I2, 2) ) / 2.0_rp   ! n_{I2}(0)
            !!!!!!!! For Current-Current correlation functions
            I4 = FNNBond(I2, 2)
            RSpCrFTauHere(I2, 39) = - GrFTT(I4, I2, 1) + GrFTT(I2, I4, 1) - GrFTT(I4, I2, 2) + GrFTT(I2, I4, 2)
            RSpCrFTauHere(I2, 40) = - GrF00(I4, I2, 1) + GrF00(I2, I4, 1) - GrF00(I4, I2, 2) + GrF00(I2, I4, 2)
         else if(.not. IfPyObsPBC) then
            !!!!!!!! For Sz-Sz correlation functions
            RSpCrFTauHere(I0, 35) = ( GrFTT(I1, I1, 2) - GrFTT(I1, I1, 1) ) / 2.0_rp
            RSpCrFTauHere(I0, 36) = ( GrF00(I2, I2, 2) - GrF00(I2, I2, 1) ) / 2.0_rp
            !!!!!!!! For density-density correlation functions
            RSpCrFTauHere(I0, 37) = ( 2.0_rp - GrFTT(I1, I1, 1) - GrFTT(I1, I1, 2) ) / 2.0_rp
            RSpCrFTauHere(I0, 38) = ( 2.0_rp - GrF00(I2, I2, 1) - GrF00(I2, I2, 2) ) / 2.0_rp
            !!!!!!!! For Current-Current correlation functions
            I3 = FNNBond(I1, 2); I4 = FNNBond(I2, 2)
            RSpCrFTauHere(I0, 39) = - GrFTT(I3, I1, 1) + GrFTT(I1, I3, 1) - GrFTT(I3, I1, 2) + GrFTT(I1, I3, 2)
            RSpCrFTauHere(I0, 40) = - GrF00(I4, I2, 1) + GrF00(I2, I4, 1) - GrF00(I4, I2, 2) + GrF00(I2, I4, 2)
         end if
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!________________________________________________________________________________________ 	  
!_________________ (1) Constant related to applying the translational symmetry __________
!________________________________________________________________________________________
      if(IfPyObsPBC) RSpCrFTauHere = RSpCrFTauHere / dble(NumNC)
!________________________________________________________________________________________ 	  
!_________________ (2) Pick up the constant and Accumulate all the results ______________
!________________________________________________________________________________________
      do I0 = 1, 40, +1 
         do Idimj = 1, NmSitePairTau, +1
            RealSpCrFTau(NTInd, Idimj, I0) = RealSpCrFTau(NTInd, Idimj, I0) + RSpCrFTauHere(Idimj, I0)*CfgConst
         enddo
      enddo
!**************************************************************************************************	  
!___________________ 2. Finalizations for this subroutine _________________________________________
!**************************************************************************************************
      if(allocated(RSpCrFTauHere)) deallocate(RSpCrFTauHere)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for ddynamic measurements _________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
      call system_clock(time2)
		TimeDyMea = TimeDyMea + TimeIntrvl(time1, time2)

   end subroutine ObDynCrFct_SpnDcp
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$