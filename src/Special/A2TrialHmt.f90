!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL PINNING FIELDS INTO THE TRIAL DENSITY MATRIX; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to prepare the trial Hamiltonian to constructing the path.
! COMMENT: Set up trial Hamiltonian H_T.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!       In this code, we write the many-body Hamiltonian H as H = H_0 + H_I = H_T + (H_I + H_0 - H_T).
!       So during the construction of paths, we use H_T instead of H_0 for the kinetic terms.
!             
!   InitTrialHmt  --> Subroutine used to set up the trial Hamiltonian H_T matrix for B_T propagator;
!
!   ExpHTkEgVlFFT --> Subroutine used to calculate Exp(-/+dt*H_T) and Exp(-/+dt*H_T/2) by FFT;
!
!   ExpHTkFullMat --> Subroutine used to obtain eigenvalues and eigenvectors of H_T matrix (no ChemP_BT term);
!   UHFReadSzUeff --> Subroutine used to read the U_eff and Sz parameters for UHF simulations;
!   GenerateHTMat --> Subroutine used to calculate the whole real space Hamiltonian matrix;
!
!   CalculateMuBT    --> Subroutine used to obtain ChemP_BT parameter;
!   SelfTuneChemP_BT --> Subroutine used to tune the ChemP_BT parameter self-consistently to arrive at the desired nT;
!   Calculate_nTOfHT --> Subroutine used to obtain electron density of H_T Hamiltonian; 
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine InitTrialHmt()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitTrialHmt()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to set up the trial Hamiltonian H_T for the simulation.
! KEYWORDS: Set up trial Hamiltonian H_T.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Constructing the trial Hamiltonian H_T includes the following steps:
!        (0) Read the effective U and magnetic moments for UHF kind of H_T if there is input files;
!        (1) Set up the H_T matrix and diagonalize it to obtain eigenvalues as well as eigenvectors;
!        (2) Determine the value of ChemP_BT for RHF or canonical ensemble simulation.
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!___________________________ Main calculations of exponential hopping matrices ________________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of initialization process _______________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'InitTrialHmt: Initialization of the trial Hamiltonian H_T!')")
      end if  
!**************************************************************************************************     
!___________________ 0. Exp(-/+dt*H_T) and Exp(-/+dt*H_T/2) matrices by two methods _______________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Some common quantities needed by both methods ____________________
!________________________________________________________________________________________
!____________________________________________________________________________      
!________________ [0] Allocate arrays for UHF-type HT Hamiltonian ___________
!____________________________________________________________________________ 
      HubbU_UHF = 0.0_rp
      allocate(MagMoment(NumNS, 2)); MagMoment = 0.0_rp
      allocate(UHF_Const(NumNS   )); UHF_Const = 0.0_rp
!____________________________________________________________________________      
!________________ [1] Eigenvalues and eigenvectors for HT Hamiltonian _______
!____________________________________________________________________________
      allocate(HT_EigValu(NumNS,        NmSpn)); HT_EigValu = 0.0_rp
      allocate(HT_EigVect(NumNS, NumNS, NmSpn)); HT_EigVect = 0.0_rp
!________________________________________________________________________________________      
!_________________ (1) Apply FFT method method for kinetic propagating __________________
!________________________________________________________________________________________
      if(FFTEXPDTHT) then
!____________________________________________________________________________      
!________________ [0] Output the information ________________________________
!____________________________________________________________________________
         if(amyid == amstr) then
            write(*, "(30x, 'Apply FFT method for exp(-/+dt*H_T) and exp(-/+dt*H_T/2)!')")
            write(*, "(30x, 'Restricted Hartree-Fock for H_T Hamiltonian.')") 
         end if
!____________________________________________________________________________  
!________________ [1] Exp(-dt*e_i) with e_i as eigenvalues of H_T ___________
!____________________________________________________________________________ 
         allocate(ExpdtOfHTe(NumNS, NmSpn, 4)); ExpdtOfHTe = 0.0_rp
!____________________________________________________________________________      
!________________ [2] Prepare FFT for for kinetic propagating _______________
!____________________________________________________________________________         
         call FFTSettingBgn()
         call ExpHTkEgVlFFT()
!________________________________________________________________________________________      
!_________________ (2) Apply matrix product method for kinetic propagating ______________
!________________________________________________________________________________________         
      else   
!____________________________________________________________________________      
!________________ [0] Output the information ________________________________
!____________________________________________________________________________
         if(amyid == amstr) then
            write(*, "(30x, 'Apply DGEMM method for exp(-/+dt*H_T) and exp(-/+dt*H_T/2)!')")
         end if
!____________________________________________________________________________      
!________________ [1] Exp(+/-dt*H_T) or Exp(+/-dt*H_T/2) ____________________
!____________________________________________________________________________
         allocate(ExpdtTryHT(NumNS, NumNS, NmSpn, 4)); ExpdtTryHT = 0.0_rp
!____________________________________________________________________________      
!________________ [2] Prepare ZGEMM for exp(-/+dt*H_T(1/2)) matrix __________
!____________________________________________________________________________ 
         call ExpHTkFullMat()
      end if
!**************************************************************************************************     
!___________________ 1. Obtain the parameter Mu_BT for the simulation _____________________________
!**************************************************************************************************
      call CalculateMuBT()
      
   end subroutine InitTrialHmt
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
!########################################################################################################################
!########################################################################################################################
!######################################### Exp(-/+dt*H_T(/2)) by FFT method #############################################
!######################################### Exp(-/+dt*H_T(/2)) by FFT method #############################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ExpHTkEgVlFFT()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ExpHTkEgVlFFT()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to prepare some quantities for the FFT calculation of multiplying Exp(-Dltau*K) 
!                   by some matrix.
! KEYWORDS: Prepare for FFT calculation of AMat*exp(-/+dt*K).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Mainly Includes to things:
!         (0) Exp(-dt*\Lambda), where Lambda is the diagonal eigenvalue matrix for H_T. Here, the eigen energies 
!                must be in the correct order for FFT calculations;
!         (1) The FFT transforms the real space matrix into reciprocal space, and it can only guarantee the 
!                reciprocal space matrix is diagonal against to k, not necessarily to spin. So if the model 
!                Hamiltonian has spin-flip terms, we also need to find the matrix diagonalizing the matrix
!                in spin space.
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use StdInOutSt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, I1p, I1m, I2, Ix, Iy, IndK
      integer SpnInd
      real(rp) KxUp, KyUp, KxDn, KyDn
      real(rp) E_FNN_Hop_Up, E_SNN_Hop_Up, E_TNN_Hop_Up
      real(rp) E_FNN_Hop_Dn, E_SNN_Hop_Dn, E_TNN_Hop_Dn
      real(rp), allocatable :: AllEigEng(:)
      real(rp), allocatable :: TmpMat(:, :, :)
!______________________________________________________________________________________________________________     
!___________________________  Calculate exponential hopping matrices by FFT method ____________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Calculate exp(-/+dt*H_T) and exp(-/+dt*H_T/2) matrices ____________________
!************************************************************************************************** 
!________________________________________________________________________________________      
!_________________ (0) Get the eigenvalue of H_T according to FFT style _________________
!________________________________________________________________________________________
      IndK = 0
      do I2 = 0, NumL2-1, +1
         KyUp = dble(I2) * 2.0_rp * rp_pi / dble(NumL2)
         KyDn = dble(I2) * 2.0_rp * rp_pi / dble(NumL2)
         if(KyUp > rp_pi) KyUp = KyUp - 2.0_rp*rp_pi
         if(KyDn > rp_pi) KyDn = KyDn - 2.0_rp*rp_pi
         do I1 = 0, NumL1/2, +1
            IndK = IndK + 1
!____________________________________________________________________________      
!________________ [0] The Kx, Ky point \in (-pi, pi] region _________________
!____________________________________________________________________________
            KxUp = dble(I1) * 2.0_rp * rp_pi / dble(NumL1)
            KxDn = dble(I1) * 2.0_rp * rp_pi / dble(NumL1)
            if(KxUp > rp_pi) KxUp = KxUp - 2.0_rp*rp_pi
            if(KxDn > rp_pi) KxDn = KxDn - 2.0_rp*rp_pi
!____________________________________________________________________________      
!________________ [1] The dispersion for NN hopping term ____________________
!____________________________________________________________________________
            if(EkDispType <= 1) then
               E_FNN_Hop_Up = 2.0_rp * Hopt1Up * ( cos(KxUp) + cos(KyUp) )
               E_FNN_Hop_Dn = 2.0_rp * Hopt1Dn * ( cos(KxDn) + cos(KyDn) )
            else if(EkDispType == 2) then
               E_FNN_Hop_Up = Hopt1Up * ( 4.0_rp - KxUp*KxUp - KyUp*KyUp )
               E_FNN_Hop_Dn = Hopt1Dn * ( 4.0_rp - KxDn*KxDn - KyDn*KyDn )
            end if
!____________________________________________________________________________      
!________________ [2] The dispersion for NNN hopping term ___________________
!____________________________________________________________________________
            E_SNN_Hop_Up = 4.0_rp * Hopt2 * cos(KxUp) * cos(KyUp)
            E_SNN_Hop_Dn = 4.0_rp * Hopt2 * cos(KxDn) * cos(KyDn)
!____________________________________________________________________________      
!________________ [3] The dispersion for 3NN hopping term ___________________
!____________________________________________________________________________
            E_TNN_Hop_Up = 2.0_rp * Hopt3 * ( cos(2.0_rp*KxUp) + cos(2.0_rp*KyUp) )
            E_TNN_Hop_Dn = 2.0_rp * Hopt3 * ( cos(2.0_rp*KxDn) + cos(2.0_rp*KyDn) )
!____________________________________________________________________________      
!________________ [4] Obtain the exponential eigenvalues ____________________
!____________________________________________________________________________
            HT_EigValu(IndK, 1) = E_FNN_Hop_Up + E_SNN_Hop_Up + E_TNN_Hop_Up + ZmFdz/2.0_rp
            HT_EigValu(IndK, 2) = E_FNN_Hop_Dn + E_SNN_Hop_Dn + E_TNN_Hop_Dn - ZmFdz/2.0_rp
!____________________________________________________________________________      
!________________ [5] Obtain the exponential EigVal --> ExpdtOfHTe __________
!____________________________________________________________________________
            do SpnInd = 1, NmSpn, +1
               ExpdtOfHTe(IndK, SpnInd, 1) = exp( - Dltau * HT_EigValu(IndK, SpnInd)          )
               ExpdtOfHTe(IndK, SpnInd, 2) = exp( + Dltau * HT_EigValu(IndK, SpnInd)          )
               ExpdtOfHTe(IndK, SpnInd, 3) = exp( - Dltau * HT_EigValu(IndK, SpnInd) / 2.0_rp )
               ExpdtOfHTe(IndK, SpnInd, 4) = exp( + Dltau * HT_EigValu(IndK, SpnInd) / 2.0_rp )     
            enddo
         enddo
      enddo
!________________________________________________________________________________________      
!_________________ (1) Scale ExpdtOfHTe by 1.0_rp/dble(NumNS) for FFT use _______________
!________________________________________________________________________________________
      ExpdtOfHTe = ExpdtOfHTe / dble(NumNS)
!________________________________________________________________________________________      
!_________________ (2) Sort the eigenvalues for both spin-up and spin-down ______________
!________________________________________________________________________________________
      allocate(AllEigEng(NumNS))
      do SpnInd = 1, NmSpn, +1
         ! Ix = 0, Iy \in [0, NumL2-1]
         Ix = 0
         do Iy = 0, NumL2-1, +1
            I1  = Iy*(NumL1/2+1) + Ix + 1
            I1p = Iy*NumL1 + Ix + 1
            AllEigEng(I1p) = HT_EigValu(I1, SpnInd)
         enddo
         ! Iy = 0, Ix \in [1, NumL1-1]
         Iy = 0
         do Ix = 1, NumL1/2, +1
            I1  = Ix + 1
            I1p =         Ix + 1
            I1m = NumL1 - Ix + 1
            AllEigEng(I1p) = HT_EigValu(I1, SpnInd)
            AllEigEng(I1m) = HT_EigValu(I1, SpnInd)
         enddo
         ! Ix \in [1, NumL1-1], Iy \in [1, NumL2-1]
         do Iy = 1, NumL2-1, +1
            do Ix = 1, NumL1/2, +1
               I1  = Iy *(NumL1/2+1) + Ix + 1
               I1p =        Iy * NumL1 +         Ix + 1
               I1m = (NumL2-Iy)* NumL1 + NumL1 - Ix + 1
               AllEigEng(I1p) = HT_EigValu(I1, SpnInd)
               AllEigEng(I1m) = HT_EigValu(I1, SpnInd)
            enddo
         enddo
         HT_EigValu(1:NumNS, SpnInd) = AllEigEng(1:NumNS)
         call QuckSortR(1, NumNS, HT_EigValu(1, SpnInd))
      enddo
      if(allocated(AllEigEng)) deallocate(AllEigEng)
!________________________________________________________________________________________      
!_________________ (3) Output energy of H_T Hamiltonian _________________________________
!________________________________________________________________________________________ 
      if(amyid == amstr) then 
         open( 291, file = "Output/00_NonInteractBand.txt", access = "append")
         write(291, "('_________________________________________________________________________')")
         write(291, "('@@@@@@@@@@@@@@@@@@ Single-particle Energy levels of H_T @@@@@@@@@@@@@@@@@')")
         write(291, "('_________________________________________________________________________')")
!____________________________________________________________________________      
!________________ [0] For The spin-up part of H_T ___________________________
!____________________________________________________________________________
         write(291, "('ExpHTkEgVlFFT --> Decoupled spin-up part in H_T:')")
         write(291, "()")
         do I1 = 1, NumNC, +1
            write(291, "(I6, A, es25.16)") I1, char(9), HT_EigValu(I1, 1)
         enddo
         if(IfFixnT) then
            write(291, "()")
            write(291, "('Total Energy = ', es25.16)") sum(HT_EigValu(1:NumNe/2, 1))
            write(291, "('  Energy gap = ', es25.16)") (HT_EigValu(NumNe/2+1, 1)-HT_EigValu(NumNe/2, 1))/2.0_rp
         end if
         write(291, "()")
         write(291, "()")
!____________________________________________________________________________      
!________________ [1] For The spin-down part of H_T _________________________
!____________________________________________________________________________
         write(291, "('ExpHTkEgVlFFT --> Decoupled spin-down part in H_T:')")
         write(291, "()")
         do I1 = 1, NumNC, +1
            write(291, "(I6, A, es25.16)") I1, char(9), HT_EigValu(I1, 2)
         enddo
         if(IfFixnT) then
            write(291, "()")
            write(291, "('Total Energy = ', es25.16)") sum(HT_EigValu(1:NumNe/2, 2))
            write(291, "('  Energy gap = ', es25.16)") (HT_EigValu(NumNe/2+1, 2)-HT_EigValu(NumNe/2, 2))/2.0_rp
         end if
         close(291)
      end if
!**************************************************************************************************     
!___________________ 1. Obtain eigenvalues and eigenvectors for H_T Hamiltonian ___________________
!______________________ to set the ULeft, DLeftVec and VLeft matrices _____________________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Allocate temporary matrix for the H_T Hamiltonian ________________
!________________________________________________________________________________________
      allocate(TmpMat(NumNS, NumNS, NmSpn)); TmpMat = 0.0_rp
!________________________________________________________________________________________      
!_________________ (1) Obtain the tight-binding Hamiltonian matrix ______________________
!________________________________________________________________________________________
      call GenerateHTMat(TmpMat)
!________________________________________________________________________________________      
!_________________ (2) Diagonalize the tight-binding Hamiltonian matrix _________________
!________________________________________________________________________________________   
      HT_EigValu = 0.0_rp; HT_EigVect = 0.0_rp
      call dMat_Diag_QMC(NumNS, TmpMat(1, 1, 1), HT_EigValu(1, 1), HT_EigVect(1, 1, 1))
!________________________________________________________________________________________      
!_________________ (3) Deallocate the allocated matrix __________________________________
!________________________________________________________________________________________  
      if(allocated(TmpMat)) deallocate(TmpMat)

   end subroutine ExpHTkEgVlFFT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
!########################################################################################################################
!########################################################################################################################
!################################# Exp(-/+dt*H_T(/2)) by Full matrix multiplication #####################################
!################################# Exp(-/+dt*H_T(/2)) by Full matrix multiplication #####################################
!########################################################################################################################
!########################################################################################################################

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ExpHTkFullMat()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ExpHTkFullMat()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the eigenvalues (without \mu_T term) and eigenvector matrix for 
!                     H_T Hamiltonian, and obtain Exp(-/+dt*H_T) and Exp(-/+dt*H_T/2). 
! KEYWORDS: Exp(-/+dt*H_T) and Exp(-/+dt*H_T/2).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input:  (none). Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use StdInOutSt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, Ix, Iy, SpnInd
      real(rp) Rtp0(4)
      real(rp), allocatable :: TmpMat(:, :, :)
!______________________________________________________________________________________________________________     
!___________________________  Calculate exponential hopping matrices without breakup __________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Allocate the matrices used in the following _______________________________
!**************************************************************************************************
      allocate(TmpMat(NumNS, NumNS, NmSpn)); TmpMat = 0.0_rp
!**************************************************************************************************     
!___________________ 1. Construct and Diagonalize the H_T Hamiltonian _____________________________
!**************************************************************************************************
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
!_________________ (0) Read the magnetic moment from UHF mean-field _____________________
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
      call UHFReadSzUeff()
!________________________________________________________________________________________      
!_________________ (1) Obtain the tight-binding Hamiltonian matrix ______________________
!________________________________________________________________________________________
      call GenerateHTMat(TmpMat)
!________________________________________________________________________________________      
!_________________ (2) Diagonalize the tight-binding Hamiltonian matrix _________________
!________________________________________________________________________________________   
      HT_EigValu = 0.0_rp; HT_EigVect = 0.0_rp
      call dMat_Diag_QMC(NumNS, TmpMat(1, 1, 1), HT_EigValu(1, 1), HT_EigVect(1, 1, 1))
!________________________________________________________________________________________      
!_________________ (3) Output energy of H_T Hamiltonian _________________________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         open( 291, file = "Output/00_NonInteractBand.txt", access = "append")
         write(291, "('_________________________________________________________________________')")
         write(291, "('@@@@@@@@@@@@@@@@@@ Single-particle Energy levels of H_T @@@@@@@@@@@@@@@@@')")
         write(291, "('_________________________________________________________________________')")
!____________________________________________________________________________      
!________________ [0] For The spin-up part of H_T ___________________________
!____________________________________________________________________________
         write(291, "('ExpHTkFullMat --> Decoupled spin-up part in H_T:')")
         write(291, "()")
         do I1 = 1, NumNS, +1
            write(291, "(I6, A, es25.16)") I1, char(9), HT_EigValu(I1, 1)
         enddo
         if(IfFixnT) then
            write(291, "()")
            write(291, "('Total Energy = ', es25.16)") sum(HT_EigValu(1:NumNe/2, 1))
            write(291, "('  Energy gap = ', es25.16)") (HT_EigValu(NumNe/2+1, 1)-HT_EigValu(NumNe/2, 1))/2.0_rp
         end if
         write(291, "()")
         write(291, "()")
!____________________________________________________________________________      
!________________ [1] For The spin-up part of H_T ___________________________
!____________________________________________________________________________
         write(291, "('ExpHTkFullMat --> Decoupled spin-down part in H_T:')")
         write(291, "()")
         do I1 = 1, NumNS, +1
            write(291, "(I6, A, es25.16)") I1, char(9), HT_EigValu(I1, 2)
         enddo
         if(IfFixnT) then
            write(291, "()")
            write(291, "('Total Energy = ', es25.16)") sum(HT_EigValu(1:NumNe/2, 2))
            write(291, "('  Energy gap = ', es25.16)") (HT_EigValu(NumNe/2+1, 2)-HT_EigValu(NumNe/2, 2))/2.0_rp
         end if
         close(291)
      end if
!**************************************************************************************************     
!___________________ 2. The exp(-/+dt*H_T) and exp(-/+dt*H_T/2) matrices __________________________
!______________________ exp(-/+dt*H_T  ) = U * exp(-/+dt*\Lambda  ) * U^+ _________________________
!______________________ exp(-/+dt*H_T/2) = U * exp(-/+dt*\Lambda/2) * U^+ _________________________
!************************************************************************************************** 
!________________________________________________________________________________________      
!_________________ (0) exp(-/+dt*\Lambda  ) * U^+ _______________________________________
!_____________________ exp(-/+dt*\Lambda/2) * U^+ _______________________________________
!________________________________________________________________________________________
      ExpdtTryHT = 0.0_rp
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, SpnInd, I2, Rtp0)
   !$OMP DO
      do I1 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            Rtp0(1) = exp( - Dltau * HT_EigValu(I1, SpnInd)          )
            Rtp0(2) = exp( + Dltau * HT_EigValu(I1, SpnInd)          )
            Rtp0(3) = exp( - Dltau * HT_EigValu(I1, SpnInd) / 2.0_rp )
            Rtp0(4) = exp( + Dltau * HT_EigValu(I1, SpnInd) / 2.0_rp )
            do I2 = 1, NumNS, +1
               ExpdtTryHT(I1, I2, SpnInd, 1) = Rtp0(1) * HT_EigVect(I2, I1, SpnInd)
               ExpdtTryHT(I1, I2, SpnInd, 2) = Rtp0(2) * HT_EigVect(I2, I1, SpnInd)
               ExpdtTryHT(I1, I2, SpnInd, 3) = Rtp0(3) * HT_EigVect(I2, I1, SpnInd)
               ExpdtTryHT(I1, I2, SpnInd, 4) = Rtp0(4) * HT_EigVect(I2, I1, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL
!________________________________________________________________________________________      
!_________________ (1) U * [exp(-/+dt*\Lambda  ) * U^+] _________________________________
!_____________________ U * [exp(-/+dt*\Lambda/2) * U^+] _________________________________
!________________________________________________________________________________________
      TmpMat = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, HT_EigVect(1, 1, 1), ExpdtTryHT(1, 1, 1, 1), &
         & 0.0_rp, TmpMat(1, 1, 1))
      call dMat_Copy_QMC(TmpMat(1, 1, 1), ExpdtTryHT(1, 1, 1, 1))
      
      TmpMat = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, HT_EigVect(1, 1, 1), ExpdtTryHT(1, 1, 1, 2), &
         & 0.0_rp, TmpMat(1, 1, 1))
      call dMat_Copy_QMC(TmpMat(1, 1, 1), ExpdtTryHT(1, 1, 1, 2))
      
      TmpMat = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, HT_EigVect(1, 1, 1), ExpdtTryHT(1, 1, 1, 3), &
         & 0.0_rp, TmpMat(1, 1, 1))
      call dMat_Copy_QMC(TmpMat(1, 1, 1), ExpdtTryHT(1, 1, 1, 3))
      
      TmpMat = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, HT_EigVect(1, 1, 1), ExpdtTryHT(1, 1, 1, 4), &
         & 0.0_rp, TmpMat(1, 1, 1))
      call dMat_Copy_QMC(TmpMat(1, 1, 1), ExpdtTryHT(1, 1, 1, 4)) 
!**************************************************************************************************     
!_________________ 2. Deallocate all the used matrices ____________________________________________
!**************************************************************************************************
      if(allocated(TmpMat)) deallocate(TmpMat)  
      
   end subroutine ExpHTkFullMat
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine UHFReadSzUeff()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  UHFReadSzUeff()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to read the effective U and Sz for the UHF kind of trial Hamiltonian H_T. 
! KEYWORDS: Read U_eff and Sz for UHF calculation.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input:  (none). Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use StdInOutSt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, Itp0(3)
      real(rp) Rtp0(6)
!______________________________________________________________________________________________________________     
!___________________________ Read the U_eff and Sz for the UHF calculations ___________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Read the U_eff and Sz parameters for the H_T Hamiltonian __________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) If the file exist, read the U_eff and Sz parameters ______________
!_____________________ And set the UHF energy constant __________________________________
!________________________________________________________________________________________
!@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
!____________ Moments from SCF-HF simulation --> Without Errorbar _________
!@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
      open(799, err = 177, file = "Output/SCF_UHF_ElectronDensity.txt", status = "old")
      read(799, *) HubbU_UHF
      do I0 = 1, NumNS, +1
         read(799, *) Itp0(1), Rtp0(1:2)
         MagMoment(I0, 1) = Rtp0(1)
         MagMoment(I0, 2) = Rtp0(2)
      enddo
      close(799)
      
      do I0 = 1, NumNS, +1
         UHF_Const(I0) = - Dltau * HubbU_UHF * MagMoment(I0, 1) * MagMoment(I0, 2)
      enddo
      
      if(amyid == amstr) then
         write(*, "(30x, 'Unrestricted Hartree-Fock from mean-field for the H_T Hamiltonian.')") 
         write(*, "(30x, '     Read Magnetic moments from SCF_UHF_ElectronDensity.txt!')") 
         write(*, "(30x, '     HubbU_UHF, n_u, n_d = ', f6.3, 2f12.8)") HubbU_UHF, MagMoment(1, 1), MagMoment(1, 2)
      end if
      
      go to 201
!@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
!____________ Moments from QMC    simulation --> With Errorbar ____________
!@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
177   continue
    
      open(799, err = 199, file = "Output/06_RnUpnDw_AvgErrB_Input.txt", status = "old")
      read(799, *) HubbU_UHF
      do I0 = 1, NumNS, +1
         read(799, *) Itp0(1:3), Rtp0(1:6)
         MagMoment(I0, 1) = Rtp0(1)
         MagMoment(I0, 2) = Rtp0(3)
      enddo
      close(799)
      
      do I0 = 1, NumNS, +1
         UHF_Const(I0) = - Dltau * HubbU_UHF * MagMoment(I0, 1) * MagMoment(I0, 2)
      enddo     
      
      if(amyid == amstr) then
         write(*, "(30x, 'Unrestricted Hartree-Fock from QMC for H_T Hamiltonian.')") 
         write(*, "(30x, '     Read Magnetic moments from 06_RnUpnDw_AvgErrB_Input.txt!')") 
         write(*, "(30x, '     HubbU_UHF, n_u, n_d = ', f6.3, 2f12.8)") HubbU_UHF, MagMoment(1, 1), MagMoment(1, 2)
      end if
      
      go to 201
!________________________________________________________________________________________      
!_________________ (1) If the file doesn't exist, just use RHF H_T Hamiltonian __________
!________________________________________________________________________________________      
199   continue
      
      HubbU_UHF = 0.0_rp
      MagMoment = 0.0_rp
      UHF_Const = 0.0_rp
      
      if(amyid == amstr) then
         write(*, "(30x, 'Restricted Hartree-Fock for H_T Hamiltonian.')") 
      end if
      
201   continue

   end subroutine UHFReadSzUeff
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine GenerateHTMat(HmltM)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  TghtBndHmltMt(HmltM)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to construct the real space Hamiltonian matrix for the model on a finite-size
!                 system of the free electron system.
! KEYWORDS: Construct H_T Hamiltonian matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Generate the real space tight-binding Hamiltonian matrix.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use RandomNumb
      use CoreParamt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) HmltM(NumNS, NumNS, NmSpn) 
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, I3, I4
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      integer tmp, tmp1, tmp2, waveGrid
      real(rp) phaseFactor, tmpValue
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      integer Ix, Iy
      integer NF, SiteParity
      real(rp) Rtp1, Rtp2
!______________________________________________________________________________________________________________     
!_____________________ Main calculations of projector matrix and exponential hopping matrices _________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Construct the noninteracting Hamiltonian matrix ___________________________
!************************************************************************************************** 
      HmltM = 0.0_rp
!________________________________________________________________________________________      
!_________________ (0) The NN hopping terms _____________________________________________
!________________________________________________________________________________________
      do I0 = 1, NumNS, +1
!____________________________________________________________________________      
!________________ [0] The bond and lattice sites ____________________________
!____________________________________________________________________________
         Ix = StList(I0, 1); Iy = StList(I0, 2)      
         I1 = FNNBond(I0, 1); I2 = FNNBond(I0, 2)
!____________________________________________________________________________      
!________________ [1] Spin-up part of NN hopping ____________________________
!____________________________________________________________________________
         Rtp1 = Hopt1Up * dble(FNNStBnd(I0, 1))
         HmltM(I0, I1, 1) = HmltM(I0, I1, 1) + Rtp1
         HmltM(I1, I0, 1) = HmltM(I1, I0, 1) + Rtp1
         
         Rtp2 = Hopt1Up * dble(FNNStBnd(I0, 2))
         HmltM(I0, I2, 1) = HmltM(I0, I2, 1) + Rtp2
         HmltM(I2, I0, 1) = HmltM(I2, I0, 1) + Rtp2
!____________________________________________________________________________      
!________________ [2] Spin-down part of NN hopping __________________________
!____________________________________________________________________________
         Rtp1 = Hopt1Dn * dble(FNNStBnd(I0, 1))         
         HmltM(I0, I1, 2) = HmltM(I0, I1, 2) + Rtp1
         HmltM(I1, I0, 2) = HmltM(I1, I0, 2) + Rtp1

         Rtp2 = Hopt1Dn * dble(FNNStBnd(I0, 2))
         HmltM(I0, I2, 2) = HmltM(I0, I2, 2) + Rtp2
         HmltM(I2, I0, 2) = HmltM(I2, I0, 2) + Rtp2
      enddo
!________________________________________________________________________________________      
!_________________ (1) The 2NN hopping terms ____________________________________________
!________________________________________________________________________________________
      if( abs(Hopt2) > rp_Eps ) then
!*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
!________________ For uniform type of t2 hopping --> Fermi surface ________________
!*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
         do I0 = 1, NumNS, +1
!____________________________________________________________________________      
!________________ [0] The bond and lattice sites ____________________________
!____________________________________________________________________________
            Ix = StList(I0, 1); Iy = StList(I0, 2)   
            I1 = SNNBond(I0, 1); I2 = SNNBond(I0, 2)
!____________________________________________________________________________      
!________________ [2] Spin-up part of NNN hopping ___________________________
!____________________________________________________________________________
            Rtp1 = Hopt2 * dble(SNNStBnd(I0, 1))
            HmltM(I0, I1, 1) = HmltM(I0, I1, 1) + Rtp1
            HmltM(I1, I0, 1) = HmltM(I1, I0, 1) + Rtp1
               
            Rtp2 = Hopt2 * dble(SNNStBnd(I0, 2))
            HmltM(I0, I2, 1) = HmltM(I0, I2, 1) + Rtp2
            HmltM(I2, I0, 1) = HmltM(I2, I0, 1) + Rtp2
!____________________________________________________________________________      
!________________ [3] Spin-down part of NNN hopping _________________________
!____________________________________________________________________________  
            Rtp1 = Hopt2 * dble(SNNStBnd(I0, 1))
            HmltM(I0, I1, 2) = HmltM(I0, I1, 2) + Rtp1
            HmltM(I1, I0, 2) = HmltM(I1, I0, 2) + Rtp1
               
            Rtp2 = Hopt2 * dble(SNNStBnd(I0, 2))
            HmltM(I0, I2, 2) = HmltM(I0, I2, 2) + Rtp2
            HmltM(I2, I0, 2) = HmltM(I2, I0, 2) + Rtp2  
         enddo
      end if
!________________________________________________________________________________________      
!_________________ (2) The 3NN hopping terms ____________________________________________
!________________________________________________________________________________________
      if( abs(Hopt3) > rp_Eps ) then
         do I0 = 1, NumNS, +1
!____________________________________________________________________________      
!________________ [0] The bond and lattice sites ____________________________
!____________________________________________________________________________
            Ix = StList(I0, 1); Iy = StList(I0, 2)
            I1 = TNNBond(I0, 1); I2 = TNNBond(I0, 2)
!____________________________________________________________________________      
!________________ [1] Spin-up part of NN hopping ____________________________
!____________________________________________________________________________
            Rtp1 = Hopt3 * dble(TNNStBnd(I0, 1))
            HmltM(I0, I1, 1) = HmltM(I0, I1, 1) + Rtp1
            HmltM(I1, I0, 1) = HmltM(I1, I0, 1) + Rtp1
         
            Rtp2 = Hopt3 * dble(TNNStBnd(I0, 2))
            HmltM(I0, I2, 1) = HmltM(I0, I2, 1) + Rtp2
            HmltM(I2, I0, 1) = HmltM(I2, I0, 1) + Rtp2
!____________________________________________________________________________      
!________________ [2] Spin-down part of 3NN hopping _________________________
!____________________________________________________________________________
            Rtp1 = Hopt3 * dble(TNNStBnd(I0, 1))
            HmltM(I0, I1, 2) = HmltM(I0, I1, 2) + Rtp1
            HmltM(I1, I0, 2) = HmltM(I1, I0, 2) + Rtp1
            
            Rtp2 = Hopt3 * dble(TNNStBnd(I0, 2))
            HmltM(I0, I2, 2) = HmltM(I0, I2, 2) + Rtp2
            HmltM(I2, I0, 2) = HmltM(I2, I0, 2) + Rtp2
         enddo
      end if
!________________________________________________________________________________________      
!_________________ (3) The Zeeman field as z-direction magnetic field ___________________
!________________________________________________________________________________________
      if( abs(ZmFdz) > rp_Eps ) then
         do I0 = 1, NumNS, +1
            HmltM(I0, I0, 1) = HmltM(I0, I0, 1) + ZmFdz/2.0_rp
            HmltM(I0, I0, 2) = HmltM(I0, I0, 2) - ZmFdz/2.0_rp
         enddo
      end if
!________________________________________________________________________________________      
!_________________ (4) The pinning field for local order parameters _____________________
!________________________________________________________________________________________
      if( abs(PinSz) >= rp_Eps ) then
!____________________________________________________________________________      
!________________ [0] Pinning AFM Sz along one edge of ribbon _______________
!____________________________________________________________________________
         if(PinSzType == 0) then
            do Iy = 1, NumL2, +1
               SiteParity = (-1)**(mod(Iy, 2))
               Ix = 1
               I1 = (Iy-1)*NumL1 + Ix
               HmltM(I1, I1, 1) = HmltM(I1, I1, 1) + PinSz/2.0_rp * dble(SiteParity)
               HmltM(I1, I1, 2) = HmltM(I1, I1, 2) - PinSz/2.0_rp * dble(SiteParity)
            enddo
!____________________________________________________________________________      
!________________ [1] Pinning AFM Sz along two edges of ribbon ______________
!____________________________________________________________________________
         else if(PinSzType == 1) then
            do Iy = 1, NumL2, +1
               SiteParity = (-1)**(mod(Iy, 2))
         
               Ix = 1
               I1 = (Iy-1)*NumL1 + Ix
               HmltM(I1, I1, 1) = HmltM(I1, I1, 1) + PinSz/2.0_rp * dble(SiteParity)
               HmltM(I1, I1, 2) = HmltM(I1, I1, 2) - PinSz/2.0_rp * dble(SiteParity)
         
               Ix = NumL1
               I2 = (Iy-1)*NumL1 + Ix
               HmltM(I2, I2, 1) = HmltM(I2, I2, 1) - PinSz/2.0_rp * dble(SiteParity)
               HmltM(I2, I2, 2) = HmltM(I2, I2, 2) + PinSz/2.0_rp * dble(SiteParity)
            enddo
!____________________________________________________________________________      
!________________ [2] Pinning AFM Sz at all lattice sites ___________________
!____________________________________________________________________________
         else if(PinSzType == 2) then
            do Ix = 1, NumL1, +1
               do Iy = 1, NumL2, +1
                  SiteParity = (-1)**(mod(Ix+Iy, 2))
                  I1 = (Iy-1)*NumL1 + Ix
                  HmltM(I1, I1, 1) = HmltM(I1, I1, 1) + PinSz/2.0_rp * dble(SiteParity)
                  HmltM(I1, I1, 2) = HmltM(I1, I1, 2) - PinSz/2.0_rp * dble(SiteParity)
               enddo
            enddo
         end if
      end if
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
!_________________ (5) The UHF mean-field term __________________________________________
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if( abs(HubbU_UHF) > rp_Eps ) then
         do I0 = 1, NumNS, +1
            HmltM(I0, I0, 1) = HmltM(I0, I0, 1) + HubbU_UHF*MagMoment(I0, 2)
            HmltM(I0, I0, 2) = HmltM(I0, I0, 2) + HubbU_UHF*MagMoment(I0, 1)
         enddo
      end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (ifSinusoidalPinning) then
         do tmp1 = 1, NumL1, +1
            waveGrid = mod((tmp1 - 1), NumL1)

            do tmp2 = 1, NumL2, +1
               tmp = (tmp2 - 1) * NumL1 + tmp1
               phaseFactor = (-1.0_rp) ** (tmp1 + tmp2)
               HmltM(tmp, tmp, 1) = HmltM(tmp, tmp, 1) &
                  &+ phaseFactor * SinusoidalPinSz/2.0_rp * cos((2.0_rp * rp_pi / LambdaSz) * waveGrid + rp_pi)
               HmltM(tmp, tmp, 2) = HmltM(tmp, tmp, 2) &
                  &- phaseFactor * SinusoidalPinSz/2.0_rp * cos((2.0_rp * rp_pi / LambdaSz) * waveGrid + rp_pi)
               tmpValue = 0 &
                  &+ phaseFactor * SinusoidalPinSz * cos((2.0_rp * rp_pi / LambdaSz) * waveGrid + rp_pi)
               write(*, *) tmp, phaseFactor, tmpValue
            end do
         end do
      end if

      
      ! Print amplitude and wavelength of the sinusoidal pinning fields to screen
      if (amyid == amstr) then
         write(*, "(35x, 'ifSinusoidalPinning == True --> applying sinusoidal spin pinning fields')")
         write(*, "(35x, 'SinusoidalPinSz == ', es23.16)") SinusoidalPinSz
         write(*, "(35x, 'LambdaSz == ', es23.16)") LambdaSz
      end if

      
      ! Check the sinusoidal wave and save the configuration into a file
      if (amyid == amstr) then
         open(300, file = "Add_Output/01_Sinusoidal_Spin_Pinning_Add.txt", access = "append")
         do tmp1 = 1, NumL1, +1
            waveGrid = mod((tmp1 - 1), NumL1)
            do tmp2 = 1, NumL2, +1 
               tmp = (tmp2 - 1) * NumL1 + tmp1
               phaseFactor = (-1.0_rp) ** (tmp1 + tmp2)
               tmpValue = phaseFactor * SinusoidalPinSz &
                  &* cos((2.0_rp * rp_pi / LambdaSz) * waveGrid + rp_pi)
               write(300, "(3I5, f12.8)") tmp, tmp1, tmp2, tmpValue
            end do
         end do 
         close(300)
      end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end subroutine GenerateHTMat
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!############################################# Determine the ChemP_BT parameter in trial H_T ############################
!############################################# Determine the ChemP_BT parameter in trial H_T ############################
!########################################################################################################################
!########################################################################################################################
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine CalculateMuBT()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CalculateMuBT()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to obtain the ChemP_BT before the simulations.
! KEYWORDS: Obtain ChemP_BT parameter.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Input:  (none). Output: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use StdInOutSt
      use RandomNumb
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      real(rp) Rtp0(6)
      real(rp) Dnsty_BT, Dnsty_BT_Err
      real(rp) nT_BT
!______________________________________________________________________________________________________________     
!___________________________ Obtain the ChemP_BT model parameter ______________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. MuBTType == 0  -->  directly use input value of ChemP _____________________
!**************************************************************************************************
      if(MuBTType == 0) then
!________________________________________________________________________________________      
!_________________ (0) Set ChemP_BT = ChemP and compute the nT_BT _______________________
!________________________________________________________________________________________
         ChemP_BT = ChemP
         call Calculate_nTOfHT(nT_BT)
!________________________________________________________________________________________      
!_________________ (1) Output the information of MuBTType, ChemP_BT and nT_BT ___________
!________________________________________________________________________________________
         if(amyid == amstr) then
            write(*, "(30x, 'MuBTType == 0 ---> Use ChemP_BT == ChemP for simulation!')")
            write(*, "(35x, 'ChemP_BT == ', es23.16)") ChemP_BT
            write(*, "(35x, '   nT_BT == ', es23.16)") nT_BT
         end if
!**************************************************************************************************     
!___________________ 1. MuBTType == 1  -->  directly use input value of ChemP_BT __________________
!**************************************************************************************************
      else if(MuBTType == 1) then
!________________________________________________________________________________________      
!_________________ (0) Use input ChemP_BT and compute the nT_BT _________________________
!________________________________________________________________________________________
         call Calculate_nTOfHT(nT_BT)
!________________________________________________________________________________________      
!_________________ (1) Output the information of MuBTType, ChemP_BT and nT_BT ___________
!________________________________________________________________________________________
         if(amyid == amstr) then
            write(*, "(30x, 'MuBTType == 1 ---> Read ChemP_BT from the input!')")
            write(*, "(35x, 'ChemP_BT == ', es23.16)") ChemP_BT
            write(*, "(35x, '   nT_BT == ', es23.16)") nT_BT
         end if
!**************************************************************************************************     
!___________________ 2. MuBTType == 2  -->  calculate ChemP_BT from input Dnsty_BX ________________
!**************************************************************************************************
      else if(MuBTType == 2) then
!________________________________________________________________________________________      
!_________________ (0) First read the input total density for B_X _______________________
!________________________________________________________________________________________
         open(57, err = 191, file = "Output/99_TotalDensity_BX_Input.txt", status = "old")
         read(57, *) Rtp0(1:6)
         close(57)
         Dnsty_BT     = Rtp0(5)
         Dnsty_BT_Err = Rtp0(6)
!________________________________________________________________________________________      
!_________________ (1) Apply the bisection calculation to obtain ChemP_BT _______________
!________________________________________________________________________________________
         call SelfTuneChemP_BT(Dnsty_BT, Dnsty_BT_Err, nT_BT)
!________________________________________________________________________________________      
!_________________ (2) Output the information of MuBTType, ChemP_BT and nT_BT ___________
!________________________________________________________________________________________
         if(amyid == amstr) then
            write(*, "(30x, 'MuBTType == 2 ---> Calculate ChemP_BT from input file of Dnsty_BX!')")
            write(*, "(35x, 'ChemP_BT == ', es23.16)") ChemP_BT
            write(*, "(35x, '   nT_BT == ', es23.16)") nT_BT
            write(*, "(35x, 'abs(nT_BT-Dnsty_BT) == ', es23.16)") abs(nT_BT-Dnsty_BT)
         end if
         go to 291
!________________________________________________________________________________________      
!_________________ (3) If the input file does not exist, use the input ChemP_BT _________
!________________________________________________________________________________________
191      continue
         call Calculate_nTOfHT(nT_BT)
!________________________________________________________________________________________      
!_________________ (4) Output the information of MuBTType, ChemP_BT and nT_BT ___________
!________________________________________________________________________________________
         if(amyid == amstr) then
            write(*, "(30x, 'MuBTType == 2. But No input file of Dnsty_BX.')")
            write(*, "(35x, '               Simply use input ChemP_BT for simulation!')")
            write(*, "(35x, 'ChemP_BT == ', es23.16)") ChemP_BT
            write(*, "(35x, '   nT_BT == ', es23.16)") nT_BT
         end if
291      continue
!**************************************************************************************************     
!___________________ 3. MuBTType == 3  -->  calculate ChemP_BT from input Fix_nT __________________
!______________________ Used for IfFixnT = .true. and tune ChemP self-consistently ________________
!**************************************************************************************************
      else if(MuBTType == 3) then
!________________________________________________________________________________________      
!_________________ (0) Obtain the Dnsty_BT and Dnsty_BT_Err _____________________________
!________________________________________________________________________________________
         Dnsty_BT     = Fix_nT
         Dnsty_BT_Err = 1.0E-5_rp
!________________________________________________________________________________________      
!_________________ (1) Apply the bisection calculation to obtain ChemP_BT _______________
!________________________________________________________________________________________
         call SelfTuneChemP_BT(Dnsty_BT, Dnsty_BT_Err, nT_BT)
!________________________________________________________________________________________      
!_________________ (2) Output the information of MuBTType, ChemP_BT and nT_BT ___________
!________________________________________________________________________________________
         if(amyid == amstr) then
            write(*, "(30x, 'MuBTType == 3 ---> Calculate ChemP_BT from input Fix_nT!')")
            write(*, "(35x, 'ChemP_BT == ', sp, es19.12)") ChemP_BT
            write(*, "(35x, '   nT_BT == ', sp, es19.12)") nT_BT
            write(*, "(35x, 'abs(nT_BT-Dnsty_BT) == ', es19.12)") abs(nT_BT-Dnsty_BT)
         end if    
!**************************************************************************************************     
!___________________ 4. MuBTType == 4  -->  Determine ChemP_BT self-consistently from _____________
!______________________ computing the density of the system _______________________________________
!**************************************************************************************************
      else if(MuBTType == 4) then
!________________________________________________________________________________________      
!_________________ (0) Set the initial ChemP_BT and compute the nT_BT ___________________
!________________________________________________________________________________________
         ChemP_BT = ChemP / 3.0_rp
         call Calculate_nTOfHT(nT_BT)
!________________________________________________________________________________________      
!_________________ (1) Output the information of MuBTType, ChemP_BT and nT_BT ___________
!________________________________________________________________________________________
         if(amyid == amstr) then
            write(*, "(30x, 'MuBTType == 4 ---> Determine ChemP_BT self-consistently by simulations!')")
            write(*, "(35x, 'The initial ChemP_BT and corresponding nT_BT:')")
            write(*, "(35x, '     ChemP_BT == ', es23.16)") ChemP_BT
            write(*, "(35x, '        nT_BT == ', es23.16)") nT_BT
         end if
      end if

end subroutine CalculateMuBT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine SelfTuneChemP_BT(Dnsty_BT, Dnsty_BT_Err, nT_BT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SelfTuneChemP_BT(Dnsty_BT, Dnsty_BT_Err, nT_BT)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is obtain the ChemP_BT parameter self-consistently from Dnsty_BT and Dnsty_BT_Err by 
!                       Bisection method.
! KEYWORDS: Calculate ChemP_BT from Dnsty_BT and Dnsty_BT_Err.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate ChemP_BT from Dnsty_BT and Dnsty_BT_Err.
!
!     Input:  Dnsty_BT     --> The desired electron density of the H_T Hamiltonian;
!             Dnsty_BT_Err --> The errorbar of input Dnsty_BT;
!
!     Output: nT_BT.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use MPISetting
      use RandomNumb
      use CoreParamt
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) Dnsty_BT, Dnsty_BT_Err
      real(rp) nT_BT
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer IterateCount
      real(rp) Rtp0, nTDiffMaxm
      real(rp) ChemP_Large_Ne, ChemP_Small_Ne
!______________________________________________________________________________________________________________     
!___________________________ Main calculations of nT for the H_T Hamiltonian __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Apply the bisection method to obtain ChemP_BT value _______________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Just initialize the large and small ChemP_BT _____________________
!_____________________ Make it slightly different for different processes _______________
!________________________________________________________________________________________    
      ChemP_BT = 0.0_rp
      Rtp0 = dble(amyid)/dble(anprc) * spring_sfmt_stream() * 0.1_rp
      ChemP_Small_Ne = + 7.0_rp + Rtp0
      ChemP_Large_Ne = - 6.0_rp - Rtp0
!________________________________________________________________________________________      
!_________________ (1) Determine the Threshold for nT_BT ________________________________
!________________________________________________________________________________________
      nTDiffMaxm = Dnsty_BT_Err * 1.0E-6_rp
!________________________________________________________________________________________      
!_________________ (2) Calculate ChemP_BT by the bisection method _______________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         open(291, file = "Output/00_Obtain_ChemP_BT.txt", access = "append")
      end if
         
      IterateCount = 0
      nT_BT = -10.0_rp
      
      do while( abs(nT_BT-Dnsty_BT) >= nTDiffMaxm )
!____________________________________________________________________________      
!________________ [0] Try a new value of ChemP_BT ___________________________
!____________________________________________________________________________
         ChemP_BT = (ChemP_Small_Ne + ChemP_Large_Ne) / 2.0_rp
!____________________________________________________________________________      
!________________ [1] Calculate nT_BT for present ChemP_BT __________________
!____________________________________________________________________________
         call Calculate_nTOfHT(nT_BT)
!____________________________________________________________________________      
!________________ [2] Reset ChemP_Small_Ne and ChemP_Large_Ne _______________
!____________________________________________________________________________
         if(nT_BT > Dnsty_BT) then
            ChemP_Small_Ne = ChemP_Small_Ne
            ChemP_Large_Ne = ChemP_BT
         else
            ChemP_Small_Ne = ChemP_BT
            ChemP_Large_Ne = ChemP_Large_Ne
         end if
!____________________________________________________________________________      
!________________ [3] Output information for every iteration ________________
!____________________________________________________________________________         
         IterateCount = IterateCount + 1
         if(amyid == amstr) then
            write(291, "(I4, 2es23.14, SPes23.14, 2es23.14)") IterateCount, ChemP_BT, nT_BT, abs(nT_BT-Dnsty_BT), &
               & ChemP_Small_Ne, ChemP_Large_Ne
         end if
      enddo 
      
      if(amyid == amstr) then    
         write(291, "()")
         write(291, "()")
         close(291)
      end if
!________________________________________________________________________________________      
!_________________ (3) Take average of results from all processes _______________________
!________________________________________________________________________________________
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)
      Rtp0 = 0.0_rp
      call MPI_ALLREDUCE(ChemP_BT, Rtp0, 1, rp_MPI_REAL, MPI_SUM, acomm, ierr)
      ChemP_BT = Rtp0 / dble(anprc)
#endif
      
   end subroutine SelfTuneChemP_BT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine Calculate_nTOfHT(nT_BT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  Calculate_nTOfHT(nT_BT)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the electron density in the H_T Hamiltonian.
! KEYWORDS: Calculate nT for H_T Hamiltonian using ChemP_BT.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate the electron density in the H_T Hamiltonian using ChemP_BT
!
!     Input:  As above.
!
!     Output: nT_BT.
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
      real(rp) nT_BT
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      integer I0
      real(rp) Rtp1, Rtp2
      real(rp), external :: FermiDiracFunct
!______________________________________________________________________________________________________________     
!___________________________ Main calculations of nT for the H_T Hamiltonian __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Calculate nT_BT for spin decoupled case ___________________________________ 
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Initialization for nT_BT _________________________________________
!________________________________________________________________________________________
      nT_BT = 0.0_rp
!________________________________________________________________________________________      
!_________________ (1) Calculate nT_BT by Fermi-Dirac distribution function _____________
!________________________________________________________________________________________
      do I0 = 1, NumNS, +1
         !!!!! Occupation of Spin-up part
         Rtp1 = HT_EigValu(I0, 1) + ChemP_BT - HubbU_UHF/2.0_rp
         nT_BT = nT_BT + FermiDiracFunct(BetaT, Rtp1)
         !!!!! Occupation of Spin-down part
         Rtp2 = HT_EigValu(I0, 2) + ChemP_BT - HubbU_UHF/2.0_rp
         nT_BT = nT_BT + FermiDiracFunct(BetaT, Rtp2)
      enddo
!________________________________________________________________________________________      
!_________________ (2) The Final result of nT_BT ________________________________________
!________________________________________________________________________________________ 
      nT_BT = nT_BT / dble(NumNS)
      
   end subroutine Calculate_nTOfHT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   real(rp) function FermiDiracFunct(Beta, Energy)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FermiDiracFunct(Beta, Energy)
! TYPE:     function
! PURPOSE:  This function is used to calculate the Fermi-Dirac distribution function in a stable way.
! KEYWORDS: Fermi-Dirac distribution function.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate the Fermi-Dirac distribution function in a stable way. 
!
!     Input: Beta   --> The inverse temperature;
!            Energy --> The energy level.
!
!     Outpt: FermiDiracFunct --> The Fermi-Dirac distribution function value.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use CoreParamt
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) Beta
      real(rp) Energy
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      real(rp) Rtp0, Rtp1
!______________________________________________________________________________________________________________     
!___________________________ Calculate the Fermi-Dirac distribution function __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Calculate the Fermi-Dirac distribution function in stable way _____________
!**************************************************************************************************
      Rtp0 = Beta * Energy
      if( Rtp0 > +200.0_rp*log(10.0_rp) ) then
         Rtp1 = 0.0_rp
      else if( (Rtp0 <= +200.0_rp*log(10.0_rp)) .and. (Rtp0 >= 0.0_rp) ) then
         Rtp1 = exp(-Rtp0) / ( 1.0_rp + exp(-Rtp0) )
      else if( (Rtp0 < 0.0_rp) .and. (Rtp0 >= -200.0_rp*log(10.0_rp)) ) then
         Rtp1 = 1.0_rp / ( 1.0_rp + exp(+Rtp0) )
      else if(Rtp0 < -200.0_rp*log(10.0_rp)) then
         Rtp1 = 1.0_rp
      end if
      FermiDiracFunct = Rtp1

   end function FermiDiracFunct
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   real(rp) function ModfdFermiDirac(Beta, TauVal, Energy)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ModfdFermiDirac(Beta, TauVal, Energy)
! TYPE:     function
! PURPOSE:  This function is used to compute exp[(beta-tau)*E]/[1+exp(beta*E)] in a stable way.
! KEYWORDS: exp[(beta-tau)*E]/[1+exp(beta*E)] function.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate the exp[(beta-tau)*E]/[1+exp(beta*E)] function in a stable way. 
!
!     Input: Beta   --> The inverse temperature;
!            TauVal --> The value of tau \in [0, beta];
!            Energy --> The energy level.
!
!     Outpt: FermiDiracFunct --> The exp[(beta-tau)*E]/[1+exp(beta*E)] function value.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use CoreParamt
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) Beta, TauVal
      real(rp) Energy
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________ 
      real(rp) Rtp0, Rtp1, Rtp2, Rtp3
!______________________________________________________________________________________________________________     
!___________________________ Calculate exp[(beta-tau)*E]/[1+exp(beta*E)] function _____________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Calculate exp[(beta-tau)*E]/[1+exp(beta*E)] function in stable way ________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Store the result of Beta*Energy __________________________________
!________________________________________________________________________________________
      Rtp0 = Beta * Energy
!________________________________________________________________________________________      
!_________________ (1) Compute exp[(beta-tau)*E]/[1+exp(beta*E)] ________________________
!________________________________________________________________________________________
      if(Energy >= 0.0_rp) then
         Rtp2 = merge(0.0_rp, exp(-Rtp0), Rtp0>+200.0_rp*log(10.0_rp))
         Rtp1 = TauVal * Energy
         Rtp3 = merge(0.0_rp, exp(-Rtp1), Rtp1>+200.0_rp*log(10.0_rp))
         ModfdFermiDirac = Rtp3 / (1.0_rp + Rtp2)
      else
         Rtp2 = merge(0.0_rp, exp(+Rtp0), Rtp0<-200.0_rp*log(10.0_rp))
         Rtp1 = (Beta-TauVal) * Energy
         Rtp3 = merge(0.0_rp, exp(+Rtp1), Rtp1<-200.0_rp*log(10.0_rp))
         ModfdFermiDirac = Rtp3 / (1.0_rp + Rtp2)
      end if

   end function ModfdFermiDirac
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$