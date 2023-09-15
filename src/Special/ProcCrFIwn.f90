!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to compute G(k,iw_n) for fermion channel and C(k,iw_m) for various bosonic channels.
! COMMENT: Compute imaginary-frequency observables.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Compute G(k,iw_n) and C(k,iw_m) for every BIN simulation.
!             
!    ComptGrFCrFkIwn --> Subroutine to calculate G(k,iw_n) and C(k,iw_n) from imag-tau data;
!
!    GreenFkIwnFermi --> Subroutine to calculate G(k,iw_n) from imag-tau data for fermionic channel;
!    GreenFkIwnBoson --> Subroutine to calculate C(k,iw_n) from imag-tau data for bosonic channels;
!    ComputeTauToIwn --> Subroutine to calculate G(k,iw_n) or C(k,iw_n) from G(k,\tau) or C(k,\tau);
!    FermiDiracFunct --> Subroutine to calculate Fermi-Dirac distribution.      
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ComptGrFCrFkIwn(NB) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ComptGrFCrFkIwn(NB) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compute the imaginary-frequency quantities from the imaginary-time ones,
!                 by the simple trapezoidal formula.
! KEYWORDS: G(k,iw_n) and C(k,iw_n).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Compute the imaginary-frequency quantities.
!
!     Input: NB --> Integer index for BIN simulations;
!
!     Outpt: (none).
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
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB
!______________________________________________________________________________________________________________     
!______________________________ Main calculations for Fourier Transformation __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. For the fermionic channels --> Green's function and self-energy ___________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Compute G(k,iw_n) from G(k,\tau) for all k pints _________________
!________________________________________________________________________________________
      call GreenFkIwnFermi(NB)
!**************************************************************************************************     
!___________________ 1. For the bosonic channels --> Correlation in various channels ______________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) C(k,iw_n) for the spin Sz-Sz correlation function ________________
!________________________________________________________________________________________
      !!!!!!!!!! For the standard correlation
      call GreenFkIwnBoson(NB, 05)
      !!!!!!!!!! For the Vertex Contribution
      call GreenFkIwnBoson(NB, 20)
!________________________________________________________________________________________      
!_________________ (1) C(k,iw_n) for the spin Sxy-Sxy correlation function ______________
!________________________________________________________________________________________
      !!!!!!!!!! For the standard correlation
      call GreenFkIwnBoson(NB, 06)
      !!!!!!!!!! For the Vertex Contribution
      call GreenFkIwnBoson(NB, 21)
!________________________________________________________________________________________      
!_________________ (2) C(k,iw_n) for the density-density correlation function ___________
!________________________________________________________________________________________
      !!!!!!!!!! For the standard correlation
      call GreenFkIwnBoson(NB, 07)
      !!!!!!!!!! For the Vertex Contribution
      call GreenFkIwnBoson(NB, 22)
!________________________________________________________________________________________      
!_________________ (3) C(k,iw_n) for the current-current correlation function ___________
!________________________________________________________________________________________
      !!!!!!!!!! For the standard correlation
      call GreenFkIwnBoson(NB, 08)
      !!!!!!!!!! For the Vertex Contribution
      call GreenFkIwnBoson(NB, 23)
!________________________________________________________________________________________      
!_________________ (4) C(k,iw_n) for on-site s-wave pair-pair correlation function ______
!________________________________________________________________________________________
      !!!!!!!!!! For the standard correlation
      call GreenFkIwnBoson(NB, 09)
      !!!!!!!!!! For the Vertex Contribution
      call GreenFkIwnBoson(NB, 24)
!________________________________________________________________________________________      
!_________________ (5) C(k,iw_n) for extended s-wave pair-pair correlation function _____
!________________________________________________________________________________________
      !!!!!!!!!! For the standard correlation
      call GreenFkIwnBoson(NB, 10)
      !!!!!!!!!! For the Vertex Contribution
      call GreenFkIwnBoson(NB, 25)
!________________________________________________________________________________________      
!_________________ (6) C(k,iw_n) for d-wave pair-pair correlation function ______________
!________________________________________________________________________________________
      !!!!!!!!!! For the standard correlation
      call GreenFkIwnBoson(NB, 11)
      !!!!!!!!!! For the Vertex Contribution
      call GreenFkIwnBoson(NB, 26)
!**************************************************************************************************     
!___________________ 2. Output the results of G(k,iw_n) and C(k,iw_n) _____________________________
!**************************************************************************************************
      if(amyid == amstr) then
         call OutputGrFCrFIwn(NB)
      end if
      
   end subroutine ComptGrFCrFkIwn
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine GreenFkIwnFermi(NB) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GreenFkIwnFermi(NB) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compute the imaginary-frequency quantities from the imaginary-time ones,
!                 by the simple trapezoidal formula, for single-particle Green's functions.
! KEYWORDS: Compute G(k,iw_n) from G(k,\tau).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Compute the imaginary-frequency G(k,iw_n).
!
!     Input: NB --> Integer index for BIN simulations;
!
!     Outpt: (none).
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
      integer NB
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Nk, FrqInd, TauInd, NumTauProc
      real(rp) TauVal
      real(rp) KxUp, KyUp, KxDn, KyDn, EpsonKUp, EpsonKDn
      complex(rp) Ztp0, Ztp1, Ztp2
      real(rp), external :: FermiDiracFunct
      real(rp), external :: ModfdFermiDirac
      real(rp), allocatable :: GrnFCrFTau(:, :)
      complex(rp), allocatable :: GrnFCrFIwn(:, :)
!______________________________________________________________________________________________________________     
!______________________________ Main calculations for computing G(k,iw_n) from G(k,\tau) ______________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Initialization for this data process ______________________________________
!**************************************************************************************************
      NumTauProc = NumTauPnt
      if(NmTDMType == 0) NumTauProc = 2*NumTauPnt
      allocate(GrnFCrFTau(0:NumTauProc, 0:1)); GrnFCrFTau = 0.0_rp
      allocate(GrnFCrFIwn(0:NmFrqFermi-1, 0:1)); GrnFCrFIwn = rp_Zzero
!**************************************************************************************************     
!___________________ 1. For the fermionic channels --> Green's function and self-energy ___________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Compute the interacting G(k,iw_n) from G(k,\tau) _________________
!________________________________________________________________________________________
      do Nk = 1, NumNC, +1
         !!!!!!!!!! First impose the anti-perodicity Imag-tau data
         if(NmTDMType == 0) then
            do TauInd = 0, NumTauPnt-1, +1
               !!!!!!!! The spin-up channel
               GrnFCrFTau(            TauInd, 0) = real( KSpCrFtTauforFT(TauInd, Nk, 1) )
               GrnFCrFTau(2*NumTauPnt-TauInd, 0) = real( KSpCrFtTauforFT(TauInd, Nk, 2) )
               !!!!!!!! The spin-down channel
               GrnFCrFTau(            TauInd, 1) = real( KSpCrFtTauforFT(TauInd, Nk, 3) )
               GrnFCrFTau(2*NumTauPnt-TauInd, 1) = real( KSpCrFtTauforFT(TauInd, Nk, 4) )
            enddo
            Ztp1 = ( KSpCrFtTauforFT(NumTauPnt, Nk, 1) + KSpCrFtTauforFT(NumTauPnt, Nk, 2) ) / 2.0_rp
            Ztp2 = ( KSpCrFtTauforFT(NumTauPnt, Nk, 3) + KSpCrFtTauforFT(NumTauPnt, Nk, 4) ) / 2.0_rp
            GrnFCrFTau(NumTauPnt, 0) = real(Ztp1)
            GrnFCrFTau(NumTauPnt, 1) = real(Ztp2)
         else if(NmTDMType == 1) then
            do TauInd = 0, NumTauPnt, +1
               Ztp1 = ( KSpCrFtTauforFT(TauInd, Nk, 1) + KSpCrFtTauforFT(NumTauPnt-TauInd, Nk, 2) ) / 2.0_rp
               Ztp2 = ( KSpCrFtTauforFT(TauInd, Nk, 3) + KSpCrFtTauforFT(NumTauPnt-TauInd, Nk, 4) ) / 2.0_rp
               GrnFCrFTau(TauInd, 0) = real(Ztp1)
               GrnFCrFTau(TauInd, 1) = real(Ztp2)
            enddo
         end if
         !!!!!!!!!! The physial single-particle Green's function
         GrnFCrFTau = - GrnFCrFTau 
         !!!!!!!!!! G(k,iw_n) = int_{0}^{\beta}G(k,\tau)e^{iw_n\tau}d\tau, spin-up and spin-down, real and Imag parts
         GrnFCrFIwn = 0.0_rp
         call ComputeTauToIwn("Fermi", NumTauProc, GrnFCrFTau(0, 0), NmFrqFermi, GrnFCrFIwn(0, 0))
         call ComputeTauToIwn("Fermi", NumTauProc, GrnFCrFTau(0, 1), NmFrqFermi, GrnFCrFIwn(0, 1))
         do FrqInd = 1, NmFrqFermi, +1
            !!!!!!!! The spin-up channel
            FermiGrF_Iwn(NB, Nk, 2*FrqInd-1) =  real(GrnFCrFIwn(FrqInd-1, 0))
            FermiGrF_Iwn(NB, Nk, 2*FrqInd  ) = aimag(GrnFCrFIwn(FrqInd-1, 0))
            !!!!!!!! The spin-down channel
            FermiGrF_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd-1) =  real(GrnFCrFIwn(FrqInd-1, 1))
            FermiGrF_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd  ) = aimag(GrnFCrFIwn(FrqInd-1, 1))
         enddo
      enddo
!________________________________________________________________________________________      
!_________________ (1) Compute the free fermion G_0(k,iw_n) from G_0(k,\tau) ____________
!_____________________ Only for NB == 1 _________________________________________________
!________________________________________________________________________________________
      if(NB == 1) then
         do Nk = 1, NumNC, +1
            !!!!!!!!!! The (Kx, Ky) point and dispersion in spin-up and down channels
            !!!!!!!! For spin-up channel
            KxUp = dble(KpList(Nk, 1)) * 2.0_rp*rp_pi/dble(NumL1)
            KyUp = dble(KpList(Nk, 2)) * 2.0_rp*rp_pi/dble(NumL2)
            if(KxUp > rp_pi) KxUp = KxUp - 2.0_rp*rp_pi
            if(KyUp > rp_pi) KyUp = KyUp - 2.0_rp*rp_pi
            if(EkDispType <= 1) then
               EpsonKUp =   + 2.0_rp*Hopt1Up*(cos(KxUp)+cos(KyUp)) + 4.0_rp*Hopt2*cos(KxUp)*cos(KyUp) &
                          & + 2.0_rp*Hopt3*(cos(2.0_rp*KxUp)+cos(2.0_rp*KyUp)) + ZmFdz/2.0_rp + ChemP
            else if(EkDispType == 2) then
               EpsonKUp = Hopt1Up*(- KxUp*KxUp - KyUp*KyUp) + ZmFdz/2.0_rp + (ChemP + 4.0_rp*Hopt1Up)
            end if
            !!!!!!!! For spin-down channel
            KxDn = dble(KpList(Nk, 1)) * 2.0_rp*rp_pi/dble(NumL1)
            KyDn = dble(KpList(Nk, 2)) * 2.0_rp*rp_pi/dble(NumL2)       
            if(KxDn > rp_pi) KxDn = KxDn - 2.0_rp*rp_pi
            if(KyDn > rp_pi) KyDn = KyDn - 2.0_rp*rp_pi  
            if(EkDispType <= 1) then 
               EpsonKDn =   + 2.0_rp*Hopt1Dn*(cos(KxDn)+cos(KyDn)) + 4.0_rp*Hopt2*cos(KxDn)*cos(KyDn) &
                          & + 2.0_rp*Hopt3*(cos(2.0_rp*KxDn)+cos(2.0_rp*KyDn)) - ZmFdz/2.0_rp + ChemP
            else if(EkDispType == 2) then
               EpsonKDn = Hopt1Dn*(- KxDn*KxDn - KyDn*KyDn) - ZmFdz/2.0_rp + (ChemP + 4.0_rp*Hopt1Dn)
            end if
            !!!!!!!!!! Compute G_0(k,iw_n) analytically
            GrnFCrFTau = 0.0_rp
            do TauInd = 0, NumTauProc, +1
               TauVal = TauPntVal0B(TauInd) * Dltau
               ! GrnFCrFTau(TauInd, 0) = - exp(-TauVal*EpsonKUp) * (1.0_rp - FermiDiracFunct(BetaT, EpsonKUp))
               ! GrnFCrFTau(TauInd, 1) = - exp(-TauVal*EpsonKDn) * (1.0_rp - FermiDiracFunct(BetaT, EpsonKDn))
               GrnFCrFTau(TauInd, 0) = - ModfdFermiDirac(BetaT, TauVal, EpsonKUp)
               GrnFCrFTau(TauInd, 1) = - ModfdFermiDirac(BetaT, TauVal, EpsonKDn)
            enddo
            !!!!!!!!!! Compute G_0(k,iw_n)
            GrnFCrFIwn = 0.0_rp
            call ComputeTauToIwn("Fermi", NumTauProc, GrnFCrFTau(0, 0), NmFrqFermi, GrnFCrFIwn(0, 0))
            call ComputeTauToIwn("Fermi", NumTauProc, GrnFCrFTau(0, 1), NmFrqFermi, GrnFCrFIwn(0, 1))
            do FrqInd = 1, NmFrqFermi, +1
               !!!!!!!! The spin-up channel
               NonIntG0kIwn(Nk, 2*FrqInd-1) =  real(GrnFCrFIwn(FrqInd-1, 0))
               NonIntG0kIwn(Nk, 2*FrqInd  ) = aimag(GrnFCrFIwn(FrqInd-1, 0))
               !!!!!!!! The spin-down channel
               NonIntG0kIwn(Nk, 2*NmFrqFermi+2*FrqInd-1) =  real(GrnFCrFIwn(FrqInd-1, 1))
               NonIntG0kIwn(Nk, 2*NmFrqFermi+2*FrqInd  ) = aimag(GrnFCrFIwn(FrqInd-1, 1))
            enddo
         enddo
      end if
!________________________________________________________________________________________      
!_________________ (2) Compute the self-energy by Dyson Equation ________________________
!_____________________ Obtain the quasi-particle weights at all k points ________________
!________________________________________________________________________________________
      do Nk = 1, NumNC, +1
         do FrqInd = 1, NmFrqFermi, +1
            !!!!!!!! The spin-up channel
            Ztp0 = cmplx(NonIntG0kIwn(    Nk, 2*FrqInd-1), NonIntG0kIwn(    Nk, 2*FrqInd), rp)
            Ztp1 = cmplx(FermiGrF_Iwn(NB, Nk, 2*FrqInd-1), FermiGrF_Iwn(NB, Nk, 2*FrqInd), rp)
            Ztp2 = rp_Z_One/Ztp0 - rp_Z_One/Ztp1
            SelfEnrg_Iwn(NB, Nk, 2*FrqInd-1) =  real(Ztp2)
            SelfEnrg_Iwn(NB, Nk, 2*FrqInd  ) = aimag(Ztp2)
            !!!!!!!! The spin-down channel
            Ztp0 = cmplx(NonIntG0kIwn(    Nk, 2*NmFrqFermi+2*FrqInd-1), NonIntG0kIwn(    Nk, 2*NmFrqFermi+2*FrqInd), rp)
            Ztp1 = cmplx(FermiGrF_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd-1), FermiGrF_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd), rp)
            Ztp2 = rp_Z_One/Ztp0 - rp_Z_One/Ztp1
            SelfEnrg_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd-1) =  real(Ztp2)
            SelfEnrg_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd  ) = aimag(Ztp2)
         enddo
         QuasiParWhgt(NB, Nk) = 1.0_rp / ( 1.0_rp - SelfEnrg_Iwn(NB, Nk, 2)/(rp_pi/BetaT) )
      enddo
!**************************************************************************************************     
!___________________ 2. Finilization for this data process ________________________________________
!**************************************************************************************************
      if(allocated(GrnFCrFTau)) deallocate(GrnFCrFTau)
      if(allocated(GrnFCrFIwn)) deallocate(GrnFCrFIwn)
      
   end subroutine GreenFkIwnFermi
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine GreenFkIwnBoson(NB, CrFIndex) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GreenFkIwnBoson(NB, CrFIndex) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compute the imaginary-frequency quantities from the imaginary-time ones,
!                 by the simple trapezoidal formula, for bosonic correlation functions.
! KEYWORDS: Compute C(k,iw_n) from C(k,\tau).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Compute the imaginary-frequency C(k,iw_n).
!
!     Input: NB --> Integer index for BIN simulations;
!            CrFIndex --> Integer index for bosonic channels.
!
!     Outpt: (none).
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
      integer NB
      integer CrFIndex
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Nk, FrqInd, TauInd, NumTauProc
      integer Idtp0, Idtpx, Idtpy
      complex(rp) Ztp0
      real(rp), allocatable :: GrnFCrFTau(:)
      complex(rp), allocatable :: GrnFCrFIwn(:)
!______________________________________________________________________________________________________________     
!______________________________ Main calculations for computing C(k,iw_n) from C(k,\tau) ______________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. For the bosonic channels --> C(k,iw_n) from C(k,\tau) _____________________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Initialization for this data process _____________________________
!________________________________________________________________________________________
      NumTauProc = NumTauPnt/2
      if(NmTDMType == 0) NumTauProc = NumTauPnt
      allocate(GrnFCrFTau(0:NumTauProc)); GrnFCrFTau = 0.0_rp
      allocate(GrnFCrFIwn(0:NmFrqBoson-1)); GrnFCrFIwn = rp_Zzero
!________________________________________________________________________________________      
!_________________ (1) Compute the interacting C(k,iw_n) from C(k,\tau) _________________
!________________________________________________________________________________________
      do Nk = 1, NumNC, +1
         !!!!!!!!!! First impose the perodicity for Imag-tau data
         if(NmTDMType == 0) then
            do TauInd = 0, NumTauPnt, +1
               GrnFCrFTau(TauInd) = real( KSpCrFtTauforFT(TauInd, Nk, CrFIndex) )
            enddo
         else if(NmTDMType == 1) then
            do TauInd = 0, NumTauPnt/2, +1
               Ztp0 = KSpCrFtTauforFT(TauInd, Nk, CrFIndex) + KSpCrFtTauforFT(NumTauPnt-TauInd, Nk, CrFIndex)
               Ztp0 = Ztp0 / 2.0_rp
               GrnFCrFTau(TauInd) = real(Ztp0)
            enddo
         end if
         !!!!!!!!!! C(k,iw_n) = int_{0}^{\beta}C(k,\tau)e^{iw_n\tau}d\tau, only save the real part
         GrnFCrFIwn = 0.0_rp
         call ComputeTauToIwn("Boson", NumTauProc, GrnFCrFTau(0), NmFrqBoson, GrnFCrFIwn(0))
         GrnFCrFIwn = GrnFCrFIwn + conjg(GrnFCrFIwn)
         do FrqInd = 1, NmFrqBoson, +1
            BosonCrF_Iwn(NB, Nk, FrqInd, CrFIndex-4) = real(GrnFCrFIwn(FrqInd-1))
         enddo
      enddo
!________________________________________________________________________________________      
!_________________ (2) Compute Drude conductivity and superfluid density ________________
!________________________________________________________________________________________
!____________________________________________________________________________      
!________________ Current-current correlation --> Whole results _____________
!____________________________________________________________________________
      if(CrFIndex == 08) then
         Idtp0 = InvKpList(0, 0); Idtpx = InvKpList(+1, 0); Idtpy = InvKpList(0, +1)
         DrudeSupfdWt(NB, 1) = (BosonCrF_Iwn(NB, Idtpx, 1, CrFIndex-4) - BosonCrF_Iwn(NB, Idtp0, 2, CrFIndex-4)) / 4.0_rp 
         DrudeSupfdWt(NB, 2) = (BosonCrF_Iwn(NB, Idtpx, 1, CrFIndex-4) - BosonCrF_Iwn(NB, Idtpy, 1, CrFIndex-4)) / 4.0_rp
         DrudeSupfdWt(NB, 3) = DrudeSupfdWt(NB, 2) / (2.0_rp*TempT/rp_pi)
!____________________________________________________________________________      
!________________ Current-current correlation --> Vertex Contribution _______
!____________________________________________________________________________
      else if(CrFIndex == 23) then
         Idtp0 = InvKpList(0, 0); Idtpx = InvKpList(+1, 0); Idtpy = InvKpList(0, +1)
         DrudeSupfdWt(NB, 4) = (BosonCrF_Iwn(NB, Idtpx, 1, CrFIndex-4) - BosonCrF_Iwn(NB, Idtp0, 2, CrFIndex-4)) / 4.0_rp 
         DrudeSupfdWt(NB, 5) = (BosonCrF_Iwn(NB, Idtpx, 1, CrFIndex-4) - BosonCrF_Iwn(NB, Idtpy, 1, CrFIndex-4)) / 4.0_rp
         DrudeSupfdWt(NB, 6) = DrudeSupfdWt(NB, 5) / (2.0_rp*TempT/rp_pi)
      end if
!________________________________________________________________________________________      
!_________________ (3) Finilization for this data process _______________________________
!________________________________________________________________________________________
      if(allocated(GrnFCrFTau)) deallocate(GrnFCrFTau)
      if(allocated(GrnFCrFIwn)) deallocate(GrnFCrFIwn)
      
   end subroutine GreenFkIwnBoson
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine OutputGrFCrFIwn(NB) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  OutputGrFCrFIwn(NB) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to output the G(k,iw_n) and C(k,iw_n) data in every BIN simulation.
! KEYWORDS: Output G(k,iw_n) and C(k,iw_n).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Output G(k,iw_n) and C(k,iw_n) in every BIN.
!
!     Input: NB --> Integer index for BIN simulations;
!
!     Outpt: (none).
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
      integer NB
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Nk, FrqInd
!______________________________________________________________________________________________________________     
!______________________________ Main calculations for output __________________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. For the fermionic channels --> Green's function and self-energy ___________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) Output G(k,iw_n) data ____________________________________________
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/GrnFctkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!!!! The spin-up channel
         do FrqInd = 1, NmFrqFermi, +1
            write(299, "(A, es17.8, A, es17.8)", advance = "no") &
                                             & char(9), FermiGrF_Iwn(NB, Nk, 2*FrqInd-1), &
                                             & char(9), FermiGrF_Iwn(NB, Nk, 2*FrqInd)
         enddo
         !!!!!!!!!! The spin-down channel
         do FrqInd = 1, NmFrqFermi, +1
            write(299, "(A, es17.8, A, es17.8)", advance = "no") &
                                             & char(9), FermiGrF_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd-1), &
                                             & char(9), FermiGrF_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd  )
         enddo
         write(299, "()")
      enddo
      close(299)
!________________________________________________________________________________________      
!_________________ (1) Output self-energy \Sigma(k,iw_n) data ___________________________
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/SlfEngkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!!!! The spin-up channel
         do FrqInd = 1, NmFrqFermi, +1
            write(299, "(A, es17.8, A, es17.8)", advance = "no") &
                                             & char(9), SelfEnrg_Iwn(NB, Nk, 2*FrqInd-1), &
                                             & char(9), SelfEnrg_Iwn(NB, Nk, 2*FrqInd)
         enddo
         !!!!!!!!!! The spin-down channel
         do FrqInd = 1, NmFrqFermi, +1
            write(299, "(A, es17.8, A, es17.8)", advance = "no") &
                                             & char(9), SelfEnrg_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd-1), &
                                             & char(9), SelfEnrg_Iwn(NB, Nk, 2*NmFrqFermi+2*FrqInd  )
         enddo
         write(299, "()")
      enddo
      close(299)
!________________________________________________________________________________________      
!_________________ (2) Output Quasi-particle weights data _______________________________
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/QusParWghtBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4, A, es17.8)") KpList(Nk, 1), char(9), KpList(Nk, 2), char(9), QuasiParWhgt(NB, Nk)
      enddo
      close(299)
!**************************************************************************************************     
!___________________ 1. For the bosonic channels --> Correlation in various channels ______________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) C(k,iw_n) for the spin Sz-Sz correlation function ________________
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/SpinZZkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!! For the standard correlation
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 01)
         enddo
         !!!!!!!! For the Vertex Contribution
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 16)
         enddo
         write(299, "()")
      enddo
      close(299)
!________________________________________________________________________________________      
!_________________ (1) C(k,iw_n) for the spin Sxy-Sxy correlation function ______________
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/SpinPMkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!! For the standard correlation
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 02)
         enddo
         !!!!!!!! For the Vertex Contribution
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 17)
         enddo
         write(299, "()")
      enddo
      close(299)
!________________________________________________________________________________________      
!_________________ (2) C(k,iw_n) for the density-density correlation function ___________
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/DenDenkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!! For the standard correlation
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 03)
         enddo
         !!!!!!!! For the Vertex Contribution
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 18)
         enddo
         write(299, "()")
      enddo
      close(299)
!________________________________________________________________________________________      
!_________________ (3) C(k,iw_n) for the current-current correlation function ___________
!________________________________________________________________________________________
      !!!!!!!!!! The C(k,iw_n) results
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/CurrntkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!! For the standard correlation
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 04)
         enddo
         !!!!!!!! For the Vertex Contribution
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 19)
         enddo
         write(299, "()")
      enddo
      close(299)
      !!!!!!!!!! Drude conductivity and superfluid density
      open( 299, file = "Add_Output/TT_GrFCrFIwnDyData/DrudeSpfldBINs.txt", access = "append")
      write(299, "(I4)", advance = "no") NB
      !!!!!!!! For the standard correlation
      write(299, "(A, es17.8, A, es17.8, A, es17.8)", advance = "no") char(9), DrudeSupfdWt(NB, 1), &
         & char(9), DrudeSupfdWt(NB, 2), char(9), DrudeSupfdWt(NB, 3)
      !!!!!!!! For the Vertex Contribution
      write(299, "(A, es17.8, A, es17.8, A, es17.8)", advance = "no") char(9), DrudeSupfdWt(NB, 4), &
         & char(9), DrudeSupfdWt(NB, 5), char(9), DrudeSupfdWt(NB, 6)
      write(299, "()")
      close(299)
!________________________________________________________________________________________      
!_________________ (4) C(k,iw_n) for on-site s-wave pair-pair correlation function ______
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/PairStkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!! For the standard correlation
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 05)
         enddo
         !!!!!!!! For the Vertex Contribution
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 20)
         enddo
         write(299, "()")
      enddo
      close(299)
!________________________________________________________________________________________      
!_________________ (5) C(k,iw_n) for extended s-wave pair-pair correlation function _____
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/EdSParkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!! For the standard correlation
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 06)
         enddo
         !!!!!!!! For the Vertex Contribution
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 21)
         enddo
         write(299, "()")
      enddo
      close(299)
!________________________________________________________________________________________      
!_________________ (6) C(k,iw_n) for d-wave pair-pair correlation function ______________
!________________________________________________________________________________________
      open(299, file = "Add_Output/TT_GrFCrFIwnDyData/DWvParkIwnBINs.txt", access = "append")
      do Nk = 1, NumNC, +1
         write(299, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
         !!!!!!!! For the standard correlation
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 07)
         enddo
         !!!!!!!! For the Vertex Contribution
         do FrqInd = 1, NmFrqBoson, +1
            write(299, "(A, es17.8)", advance = "no") char(9), BosonCrF_Iwn(NB, Nk, FrqInd, 22)
         enddo
         write(299, "()")
      enddo
      close(299)

   end subroutine OutputGrFCrFIwn
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ComputeTauToIwn(FBType, NmTau, CrFctTau, NmFrq, CrFctIwN) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ComputeTauToIwn(FBType, NmTau, CrFctTau, NmFrq, CrFctIwN) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compute the imaginary-frequency quantities from the imaginary-time ones,
!                 by the simple trapezoidal formula, for both fermions and bosons.
! KEYWORDS: Compute G(k,iw_n) from G(k,\tau).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Compute the imaginary-frequency G(k,iw_n).
!
!     Input: FBType   --> For fermion or boson channels;
!            NmTau    --> Number of tau points to process;
!            CrFctTau --> The imag-tau data as input;
!            NmFrq    --> Number of Imag-iw_n points;
!
!     Outpt: CrFctIwN --> The iw_n data as output.
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
      character(5) FBType
      integer NmTau, NmFrq
      real(rp) CrFctTau(0:NmTau)
      complex(rp) CrFctIwN(0:NmFrq-1)
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Nk, TauInd, FrqInd
      real(rp) TauBgn, TauEnd
      real(rp) CefSlope, CefItcpt
      real(rp) Omega0, OmegaN, Rtp0, Rtp1
      real(rp) RealPart, ImagPart
!______________________________________________________________________________________________________________     
!______________________________ Main calculations for \tau --> iw_n integral __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. For the fermionic channels --> Green's function and self-energy ___________
!**************************************************************************************************
!________________________________________________________________________________________      
!_________________ (0) The first iw_0 for boson or fermion, initialization ______________
!________________________________________________________________________________________
      if(Trim(FBType) == "Fermi") then
         Omega0 = rp_pi / BetaT
      else if(Trim(FBType) == "Boson") then
         Omega0 = 0.0_rp
      end if
      CrFctIwN = rp_Zzero
!________________________________________________________________________________________      
!_________________ (1) The iteration of all k points ____________________________________
!________________________________________________________________________________________
      do TauInd = 0, NmTau-1, +1
         !!!!!!!!!! Value of \tau for TauInd index
         TauBgn = TauPntVal0B(TauInd  ) * Dltau
         TauEnd = TauPntVal0B(TauInd+1) * Dltau
         !!!!!!!!!! Fit G(tau) in [TauInd, TauInd+1]*Dltau region == An*\tau + Bn
         CefSlope = ( CrFctTau(TauInd+1) - CrFctTau(TauInd) ) / ( TauEnd - TauBgn )
         CefItcpt = ( TauEnd*CrFctTau(TauInd) - TauBgn*CrFctTau(TauInd+1) )  / ( TauEnd - TauBgn )
         !!!!!!!!!! Iteration for all the iw_n
         do FrqInd = 0, NmFrq-1, +1
            !!!!!!!! The value of w_n
            OmegaN = Omega0 + dble(2*FrqInd)*rp_pi/BetaT
            if(abs(OmegaN) <= rp_Eps) then
               !!!!!!!! The real part of the numerical integral
               Rtp0 = 0.5_rp*CefSlope*TauBgn*TauBgn + CefItcpt*TauBgn
               Rtp1 = 0.5_rp*CefSlope*TauEnd*TauEnd + CefItcpt*TauEnd
               RealPart = Rtp1 - Rtp0
               !!!!!!!! The imaginary part of the numerical integral
               ImagPart = 0.0_rp
            else
               !!!!!!!! The real part of the numerical integral
               Rtp0 = + CefItcpt/OmegaN * sin(OmegaN*TauBgn) &
                              & + CefSlope/OmegaN**2 * ( +OmegaN*TauBgn*sin(OmegaN*TauBgn) + cos(OmegaN*TauBgn) )
               Rtp1 = + CefItcpt/OmegaN * sin(OmegaN*TauEnd) &
                              & + CefSlope/OmegaN**2 * ( +OmegaN*TauEnd*sin(OmegaN*TauEnd) + cos(OmegaN*TauEnd) )
               RealPart = Rtp1 - Rtp0
               !!!!!!!! The imaginary part of the numerical integral
               Rtp0 = - CefItcpt/OmegaN * cos(OmegaN*TauBgn) &
                              & + CefSlope/OmegaN**2 * ( -OmegaN*TauBgn*cos(OmegaN*TauBgn) + sin(OmegaN*TauBgn) )
               Rtp1 = - CefItcpt/OmegaN * cos(OmegaN*TauEnd) &
                              & + CefSlope/OmegaN**2 * ( -OmegaN*TauEnd*cos(OmegaN*TauEnd) + sin(OmegaN*TauEnd) )
               ImagPart = Rtp1 - Rtp0
            end if
            !!!!!!!! Accumulate GrnFtIwN result by combining real and imaginary parts
            CrFctIwN(FrqInd) = CrFctIwN(FrqInd) + cmplx(RealPart, ImagPart, rp)
         enddo
      enddo

   end subroutine ComputeTauToIwn
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$