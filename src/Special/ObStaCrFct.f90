!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the measurements of real-space correlation functions in the 
!                     DQMC simulations. For both periodic and open boundary conditions.
! COMMENT: Measurements of real-space correlations.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-20
! PURPOSE: Different subroutines are introduced as following:
!             
!   ObStaCrFct_SpnDcp --> Subroutine to measure correlations for spin decoupled case with both PBC and OBC.
!               
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ObStaCrFct_SpnDcp(CfgConst, GrF, GrFC, RealSpCrF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ObStaCrFct_SpnDcp(CfgConst, GrF, GrFC, RealSpCrF)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate r-space correlation functions for both PBC and OBC cases.
! KEYWORDS: Calculate r-space correlation functions.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-11-04
! DESCRIPTION: Measure real space correlation function under PBC.
!
!     Notice: We basically measure <C_i C_j> in the soubroutine as correlation functions.
!                  For the PERIODIC boundary conditions, <C_i C_j> only depends on (i-j), so here we impose the 
!                        translational symmetry during the measurements;
!                  For the OPEN     boundary conditions, since the <C_i C_j> matrix must be Hermitian matrix,
!                        we measure the (j\in[1,NumNC], i\in[j, NumNC]) (upper triangular part).
!
!     For the index in RealSpCrF(:, 01:40), 01~02 for single-particle Green's functions;
!                                           03~16 for regular correlations; 
!                                           17~30 for vertex contributions of all correlation functions; 
!                                           31~34 for bond-bond spin-singlet correlations.
! 
!     Input:  CfgConst --> The configuration/walker related constant (like sign/phase or walker weight);
!             GrF      --> Green's Function matrice <c_i * c_j^+>;
!             GrFC     --> Green's Function matrice <c_i^+ * c_j>; 
!     
!     Output: RealSpCrF --> Accumulated measuring results of r-space correlation functions.
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
      real(rp) GrF (NumNS, NumNS, 2)
      real(rp) GrFC(NumNS, NumNS, 2)
      real(rp) RealSpCrF(NmSitePair, 40)
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, I3, I4, Kd, Hd, Idimj
      real(rp) RlCrFTmp(15)
      real(rp), allocatable :: RealSpCrFHere(:, :)
!______________________________________________________________________________________________________________     
!__________________________ Main calculations of real space correlation functions _____________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Initializations for this subroutine _______________________________________
!************************************************************************************************** 
!________________________________________________________________________________________      
!_________________ (0) Temporary matrix for r-space correlation functions _______________
!________________________________________________________________________________________
      allocate(RealSpCrFHere(NmSitePair, 40)); RealSpCrFHere = 0.0_rp
!**************************************************************************************************     
!___________________ 1. Measure r-space correlation functions and post-process ____________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Measure all the r-space correlation functions ____________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(Idimj, I0, I1, I2, I3, I4, Kd, Hd, RlCrFTmp)
   !$OMP DO REDUCTION(+ : RealSpCrFHere)
      do I2 = 1, NumNC, +1
         do I1 = 1, merge(NumNC, I2, IfPyObsPBC), +1
!____________________________________________________________________________      
!________________ The integer index for the present term ____________________
!____________________________________________________________________________
            Idimj = merge(IminusJ(I1, I2), (I2-1)*I2/2+I1, IfPyObsPBC)
!____________________________________________________________________________      
!________________ Spin-up and down single-particle Green's function _________
!____________________________________________________________________________
            !!!!!!!!!! Spin-up channel
            RlCrFTmp(1) = GrFC(I1, I2, 1)
            RealSpCrFHere(Idimj, 01) = RealSpCrFHere(Idimj, 01) + RlCrFTmp(1)
            !!!!!!!!!! Spin-down channel
            RlCrFTmp(2) = GrFC(I1, I2, 2)
            RealSpCrFHere(Idimj, 02) = RealSpCrFHere(Idimj, 02) + RlCrFTmp(2)   
!____________________________________________________________________________      
!________________ Spin correlations, <SzSz> and <S+S- + S-S+>/2 _____________
!____________________________________________________________________________
            !!!!!!!!!! The <SzSz> correlation
            RlCrFTmp(1) =   + ( GrFC(I1, I1, 1) - GrFC(I1, I1, 2) )*( GrFC(I2, I2, 1) - GrFC(I2, I2, 2) ) &
                          & +   GrFC(I1, I2, 1)*GrF(I1, I2, 1) + GrFC(I1, I2, 2)*GrF(I1, I2, 2)
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            RealSpCrFHere(Idimj, 03) = RealSpCrFHere(Idimj, 03) + RlCrFTmp(1)
            !!!!!!!!!! The <S+S- + S-S+>/2 correlation
            RlCrFTmp(2) = GrFC(I1, I2, 1)*GrF(I1, I2, 2) + GrFC(I1, I2, 2)*GrF(I1, I2, 1)
            RlCrFTmp(2) = RlCrFTmp(2) / 2.0_rp
            RealSpCrFHere(Idimj, 04) = RealSpCrFHere(Idimj, 04) + RlCrFTmp(2)           
!____________________________________________________________________________      
!________________ Density-density correlation function ______________________
!____________________________________________________________________________
            RlCrFTmp(1) =   + ( GrFC(I1, I1, 1) + GrFC(I1, I1, 2) )*( GrFC(I2, I2, 1) + GrFC(I2, I2, 2) ) &
                          & +   GrFC(I1, I2, 1)*GrF(I1, I2, 1) + GrFC(I1, I2, 2)*GrF(I1, I2, 2)
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            RealSpCrFHere(Idimj, 05) = RealSpCrFHere(Idimj, 05) + RlCrFTmp(1)
!____________________________________________________________________________      
!________________ On-site spin-singlet s-wave pairing correlation ___________
!____________________________________________________________________________      
            RlCrFTmp(1) = GrFC(I1, I2, 1)*GrFC(I1, I2, 2) + GrF(I1, I2, 1)*GrF(I1, I2, 2)
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            RealSpCrFHere(Idimj, 06) = RealSpCrFHere(Idimj, 06) + RlCrFTmp(1)
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
                  RlCrFTmp(15) =   + GrFC(I1, I2, 1)*GrFC(I3, I4, 2) + GrFC(I1, I4, 1)*GrFC(I3, I2, 2) &
                                 & + GrFC(I3, I2, 1)*GrFC(I1, I4, 2) + GrFC(I3, I4, 1)*GrFC(I1, I2, 2)
                  !!!!!!!! Accumulate the pairing results
                  !!!!!! For the NN extended s-wave pairing
                  RlCrFTmp(1) = RlCrFTmp(1) + RlCrFTmp(15)
                  !!!!!! For the NN d-wave pairing
                  RlCrFTmp(2) = RlCrFTmp(2) + dble(1-2*mod(Kd, 2)) * dble(1-2*mod(Hd, 2)) * RlCrFTmp(15)
               enddo
            enddo
            !!!!!!!!!! For the extended s-wave pairing
            RlCrFTmp(1) = RlCrFTmp(1) / 32.0_rp
            RealSpCrFHere(Idimj, 07) = RealSpCrFHere(Idimj, 07) + RlCrFTmp(1)
            !!!!!!!!!! For the d-wave pairing
            RlCrFTmp(2) = RlCrFTmp(2) / 32.0_rp
            RealSpCrFHere(Idimj, 08) = RealSpCrFHere(Idimj, 08) + RlCrFTmp(2)
!____________________________________________________________________________      
!________________ (NN-bond)-(NN-bond) spin-singlet correlation ______________
!____________________________________________________________________________
            !!!!!!!!!! Initialization for the accumulation
            RlCrFTmp = 0.0_rp
            !!!!!!!!!! v-v, v-h, h-v, h-h, v-->vertical bond, h-->horizontal bond; 1 for vertical, 2 for horizontal
            do Hd = 1, 2, +1
               I4 = FNNBond(I2, Hd)
               do Kd = 1, 2, +1
                  I3 = FNNBond(I1, Kd)
                  I0 = (Hd-1)*2 + Kd
                  RlCrFTmp(I0) =   + GrFC(I1, I2, 1)*GrFC(I3, I4, 2) + GrFC(I1, I4, 1)*GrFC(I3, I2, 2) &
                                 & + GrFC(I3, I2, 1)*GrFC(I1, I4, 2) + GrFC(I3, I4, 1)*GrFC(I1, I2, 2)
                  RlCrFTmp(I0) = RlCrFTmp(I0) / 2.0_rp
               enddo
            enddo
            !!!!!!!!!! Accumulate the results
            RealSpCrFHere(Idimj, 31:34) = RealSpCrFHere(Idimj, 31:34) + RlCrFTmp(1:4)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!________________________________________________________________________________________      
!_________________ (1) Constant related to applying the translational symmetry __________
!________________________________________________________________________________________
      if(IfPyObsPBC) RealSpCrFHere = RealSpCrFHere / dble(NumNC)
!________________________________________________________________________________________      
!_________________ (2) Pick up the constant and Accumulate all the results ______________
!________________________________________________________________________________________
      RealSpCrF = RealSpCrF + RealSpCrFHere * CfgConst
!**************************************************************************************************     
!___________________ 2. Finalizations for this subroutine _________________________________________
!**************************************************************************************************
      if(allocated(RealSpCrFHere)) deallocate(RealSpCrFHere)
      
   end subroutine ObStaCrFct_SpnDcp
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$