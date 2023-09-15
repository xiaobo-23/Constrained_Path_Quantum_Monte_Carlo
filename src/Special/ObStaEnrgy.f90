!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL SPIN PINNING FIELDS; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the measurements of energies in the CPMC simulations.
! COMMENT: Energy measurement.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   ObStaEnrgy_SpnDcp --> Subroutine to perform measurements of energies.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine ObStaEnrgy_SpnDcp(CfgConst, GrF, GrFC, EngOccCrF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ObStaEnrgy_SpnDcp(CfgConst, GrF, GrFC, EngOccCrF)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the expectation energies for all the terms in the model Hamiltonian,
!               as well as the first-order derivatives over model parameters, and phase, sign.
! KEYWORDS: Calculate expectation energy.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate the expectation energies for all the terms in the model Hamiltonian.
!
!     Input:  CfgConst --> The configuration/walker related constant (like sign/phase or walker weight);
!             GrF      --> Green's Function matrice <c_i * c_j^+>;
!             GrFC     --> Green's Function matrice <c_i^+ * c_j>; 
!     
!     Output: EngOccCrF --> Accumulated measuring results of energies and occupations.
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
      real(rp) EngOccCrF(40)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, Ix, Iy, SiteParity
		real(rp) Rtp1, Rtp2, Rtp3
      real(rp)            HmOvHopt2, HmOvHopt3, HmOvZmFld, HmOvPinSz, HmOvHubbU, HmOvHbUCh
      real(rp) EgHopt1  , EgHopt2  , EgHopt3  , EgZmFld  , EgPinSz  , EgDopCh  , EgHubbU  , EgTotal
      real(rp) EgHbUCh, EgTotCh
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      integer tmp, tmp1, tmp2, waveGrid
      real(rp) phaseFactor, tmpHmOvSinusoidalPinSz, EgSinusoidalPinSz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      real(rp) EgnOccpUp, EgnOccpDw, EgRDouOcc, EgRNNDCrF
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Measure the expectation energies for the model ____________________________
!**************************************************************************************************           
!________________________________________________________________________________________ 	  
!_________________ (0) Expectation Energy of NN hopping term --> spin-independent _______
!________________________________________________________________________________________
      EgHopt1 = 0.0_rp
      do I0 = 1, NumNS, +1
!____________________________________________________________________________ 	  
!________________ [0] The two NN lattice sites of I0 site ___________________
!____________________________________________________________________________         
			I1 = FNNBond(I0, 1)
         I2 = FNNBond(I0, 2)
!____________________________________________________________________________ 	  
!________________ [1] Spin-up part of NN hopping ____________________________
!____________________________________________________________________________
         Rtp1 = Hopt1Up * dble(FNNStBnd(I0, 1))
         Rtp2 = Hopt1Up * dble(FNNStBnd(I0, 2))
         EgHopt1 = EgHopt1 + ( GrFC(I0, I1, 1) + GrFC(I1, I0, 1) ) * Hopt1Up * dble(FNNStBnd(I0, 1))
         EgHopt1 = EgHopt1 + ( GrFC(I0, I2, 1) + GrFC(I2, I0, 1) ) * Hopt1Up * dble(FNNStBnd(I0, 2))
!____________________________________________________________________________ 	  
!________________ [2] Spin-down part of NN hopping __________________________
!____________________________________________________________________________
         Rtp1 = Hopt1Dn * dble(FNNStBnd(I0, 1))
         Rtp2 = Hopt1Dn * dble(FNNStBnd(I0, 2))
         EgHopt1 = EgHopt1 + ( GrFC(I0, I1, 2) + GrFC(I1, I0, 2) ) * Hopt1Dn * dble(FNNStBnd(I0, 1))
         EgHopt1 = EgHopt1 + ( GrFC(I0, I2, 2) + GrFC(I2, I0, 2) ) * Hopt1Dn * dble(FNNStBnd(I0, 2))
      enddo
      EgHopt1 = EgHopt1 / dble(NumNS)
!________________________________________________________________________________________ 	  
!_________________ (1) Expectation Energy of 2NN hopping term --> spin-independent ______
!________________________________________________________________________________________
      HmOvHopt2 = 0.0_rp; EgHopt2 = 0.0_rp
!*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
!________________ For uniform type of t2 hopping --> Fermi surface ________________
!*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&  
      do I0 = 1, NumNS, +1
!____________________________________________________________________________ 	  
!________________ [0] The two SNN lattice sites of I0 site __________________
!____________________________________________________________________________         
			I1 = SNNBond(I0, 1)
         I2 = SNNBond(I0, 2)
!____________________________________________________________________________	  
!________________ [1] Spin-up part of NNN hopping ___________________________
!____________________________________________________________________________
         Rtp1 = dble(SNNStBnd(I0, 1))
         Rtp2 = dble(SNNStBnd(I0, 2))
         HmOvHopt2 = HmOvHopt2 + ( GrFC(I0, I1, 1) + GrFC(I1, I0, 1) ) * Rtp1
         HmOvHopt2 = HmOvHopt2 + ( GrFC(I0, I2, 1) + GrFC(I2, I0, 1) ) * Rtp2
!____________________________________________________________________________ 	  
!________________ [2] Spin-down part of NNN hopping _________________________
!____________________________________________________________________________
         Rtp1 = dble(SNNStBnd(I0, 1))
         Rtp2 = dble(SNNStBnd(I0, 2))
         HmOvHopt2 = HmOvHopt2 + ( GrFC(I0, I1, 2) + GrFC(I1, I0, 2) ) * Rtp1
         HmOvHopt2 = HmOvHopt2 + ( GrFC(I0, I2, 2) + GrFC(I2, I0, 2) ) * Rtp2
      enddo
      HmOvHopt2 = HmOvHopt2 / dble(NumNS)
      if(abs(Hopt2) > rp_Eps) EgHopt2 = HmOvHopt2 * Hopt2
!________________________________________________________________________________________ 	  
!_________________ (2) Expectation Energy of 3NN hopping term --> spin-independent ______
!________________________________________________________________________________________
      HmOvHopt3 = 0.0_rp; EgHopt3 = 0.0_rp
      do I0 = 1, NumNS, +1
!____________________________________________________________________________ 	  
!________________ [0] Lattice sites of NN ___________________________________
!____________________________________________________________________________
         I1 = TNNBond(I0, 1)
         I2 = TNNBond(I0, 2)
!____________________________________________________________________________ 	  
!________________ [1] Spin-up part of NN hopping ____________________________
!____________________________________________________________________________
         HmOvHopt3 = HmOvHopt3 + ( GrFC(I0, I1, 1) + GrFC(I1, I0, 1) ) * dble(TNNStBnd(I0, 1))
         HmOvHopt3 = HmOvHopt3 + ( GrFC(I0, I2, 1) + GrFC(I2, I0, 1) ) * dble(TNNStBnd(I0, 2))
!____________________________________________________________________________ 	  
!________________ [2] Spin-down part of NN hopping __________________________
!____________________________________________________________________________
         HmOvHopt3 = HmOvHopt3 + ( GrFC(I0, I1, 2) + GrFC(I1, I0, 2) ) * dble(TNNStBnd(I0, 1))
         HmOvHopt3 = HmOvHopt3 + ( GrFC(I0, I2, 2) + GrFC(I2, I0, 2) ) * dble(TNNStBnd(I0, 2))
      enddo
      HmOvHopt3 = HmOvHopt3 / dble(NumNS)
      if(abs(Hopt3) > rp_Eps) EgHopt3 = HmOvHopt3 * Hopt3
!________________________________________________________________________________________ 	  
!_________________ (3) Expectation Energy of the Zeeman field in z-direction ____________
!________________________________________________________________________________________
      HmOvZmFld = 0.0_rp; EgZmFld = 0.0_rp   
      do I0 = 1, NumNS, +1
         HmOvZmFld = HmOvZmFld + GrFC(I0, I0, 1) - GrFC(I0, I0, 2)
      enddo
      HmOvZmFld = HmOvZmFld / dble(NumNS)
      if(ZmFdz > rp_Eps) EgZmFld = HmOvZmFld * ZmFdz/2.0_rp
!________________________________________________________________________________________ 	  
!_________________ (4) Expectation Energy of pinning Sz term ____________________________
!____________________________ --> spin-dependent ________________________________________
!________________________________________________________________________________________
      HmOvPinSz = 0.0_rp; EgPinSz = 0.0_rp
      if(abs(PinSz) >= rp_Eps) then
!____________________________________________________________________________ 	  
!________________ [0] Pinning AFM Sz along one edge of ribbon _______________
!____________________________________________________________________________
         if(PinSzType == 0) then
            do Iy = 1, NumL2, +1
               SiteParity = (-1)**(mod(Iy, 2))
               Ix = 1
               I1 = (Iy-1)*NumL1 + Ix
               HmOvPinSz = HmOvPinSz + ( GrFC(I1, I1, 1) - GrFC(I1, I1, 2) ) * dble(SiteParity)
            enddo
!____________________________________________________________________________ 	  
!________________ [1] Pinning AFM Sz along two edges of ribbon ______________
!____________________________________________________________________________
         else if(PinSzType == 1) then
            do Iy = 1, NumL2, +1
               SiteParity = (-1)**(mod(Iy, 2))
               Ix = 1; I1 = (Iy-1)*NumL1 + Ix
               HmOvPinSz = HmOvPinSz + ( GrFC(I1, I1, 1) - GrFC(I1, I1, 2) ) * dble(SiteParity)
               Ix = NumL1; I2 = (Iy-1)*NumL1 + Ix
               HmOvPinSz = HmOvPinSz - ( GrFC(I2, I2, 1) - GrFC(I2, I2, 2) ) * dble(SiteParity)
            enddo
!____________________________________________________________________________ 	  
!________________ [2] Pinning AFM Sz at all lattice sites ___________________
!____________________________________________________________________________
         else if(PinSzType == 2) then
            do Ix = 1, NumL1, +1
               do Iy = 1, NumL2, +1
                  SiteParity = (-1)**(mod(Ix+Iy, 2))
                  I1 = (Iy-1)*NumL1 + Ix
                  HmOvPinSz = HmOvPinSz + ( GrFC(I1, I1, 1) - GrFC(I1, I1, 2) ) * dble(SiteParity)
               enddo
            enddo
         end if
         HmOvPinSz = HmOvPinSz / dble(NumNS)
         EgPinSz   = HmOvPinSz * PinSz/2.0_rp
      end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Energy of the sinusoidal spin pinning fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      tmpHmOvSinusoidalPinSz = 0.0_rp
      EgSinusoidalPinSz = 0.0_rp

      if (ifSinusoidalPinning) then
         do tmp1 = 1, NumL1, +1 
            waveGrid = mod((tmp1 - 1), NumL1)
            do tmp2 = 1, NumL2, +1
               tmp = (tmp2 - 1) * NumL1 + tmp1
               phaseFactor = (-1.0_rp) ** (tmp1 + tmp2)
               tmpHmOvSinusoidalPinSz = tmpHmOvSinusoidalPinSz & 
               & + phaseFactor * ( GrFC(tmp, tmp, 1) - GrFC(tmp, tmp, 2) ) * cos((2.0_rp * rp_pi / LambdaSz) * waveGrid + rp_pi) 
            end do
         end do
         tmpHmOvSinusoidalPinSz = tmpHmOvSinusoidalPinSz / dble(NumNS)
         EgSinusoidalPinSz = tmpHmOvSinusoidalPinSz * SinusoidalPinSz/2.0_rp
      end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!________________________________________________________________________________________ 	  
!_________________ (5) Expectation Energy of dopping chemical potential _________________
!____________________________ --> spin-independent ______________________________________
!________________________________________________________________________________________
      !!!!!!!!!! Filling for spin-up and spin-down channels
      Rtp1 = 0.0_rp; Rtp2 = 0.0_rp
      do I0 = 1, NumNS, +1
         Rtp1 = Rtp1 + GrFC(I0, I0, 1); Rtp2 = Rtp2 + GrFC(I0, I0, 2)
      enddo
      Rtp1 = Rtp1 / dble(NumNS); Rtp2 = Rtp2 / dble(NumNS)
      !!!!!!!!!! The normal part of Chemical potential term energy
      EgDopCh = 0.0_rp
      if(abs(ChemP) > rp_Eps) EgDopCh = ChemP * (Rtp1 + Rtp2)
      !!!!!!!!!! For EkDispType == 1 case, modify EgHopt1 and EgDopCh
      if(EkDispType == 1) then
         Rtp3 = 4.0_rp * ( Hopt1Up*Rtp1 + Hopt1Dn*Rtp2 )
         EgHopt1 = EgHopt1 - Rtp3
         EgDopCh = EgDopCh + Rtp3
      end if
!________________________________________________________________________________________ 	  
!_________________ (6) Expectation Energy of HubbU term --> spin-dependent ______________
!________________________________________________________________________________________
      HmOvHbUCh = 0.0_rp; EgHbUCh = 0.0_rp
      HmOvHubbU = 0.0_rp; EgHubbU = 0.0_rp
      do I0 = 1, NumNS, +1
         Rtp1 = GrFC(I0, I0, 1) * GrFC(I0, I0, 2)
         HmOvHbUCh = HmOvHbUCh + Rtp1
         HmOvHubbU = HmOvHubbU + Rtp1 - ( GrFC(I0, I0, 1) + GrFC(I0, I0, 2) ) / 2.0_rp
      enddo
      HmOvHbUCh = HmOvHbUCh / dble(NumNS)
      HmOvHubbU = HmOvHubbU / dble(NumNS)
      if(abs(HubbU) > rp_Eps) then
         EgHbUCh = HmOvHbUCh * HubbU
         EgHubbU = HmOvHubbU * HubbU
      end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	  
! (7) Revised Version of The Total Energy
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      EgTotCh = EgHopt1 + EgHopt2 + EgHopt3 + EgZmFld + EgPinSz + EgSinusoidalPinSz           + EgHbUCh  
      EgTotal = EgHopt1 + EgHopt2 + EgHopt3 + EgZmFld + EgPinSz + EgSinusoidalPinSz + EgDopCh + EgHubbU
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!**************************************************************************************************	  
!___________________ 1. The fillings, double occupation, and density correlations _________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Total occupation number of electrons _____________________________
!________________________________________________________________________________________
      EgnOccpUp = 0.0_rp; EgnOccpDw = 0.0_rp
      do I0 = 1, NumNS, +1
         EgnOccpUp = EgnOccpUp + GrFC(I0, I0, 1)
         EgnOccpDw = EgnOccpDw + GrFC(I0, I0, 2)
      enddo
      EgnOccpUp = EgnOccpUp / dble(NumNS)
      EgnOccpDw = EgnOccpDw / dble(NumNS)
!________________________________________________________________________________________ 	  
!_________________ (1) The Double occupation of the model system ________________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [0] For PBC, calculate average double occupancy ___________
!____________________________________________________________________________
      if(IfPyObsPBC) then
         EgRDouOcc = 0.0_rp
         do I0 = 1, NumNS, +1
            EgRDouOcc = EgRDouOcc + GrFC(I0, I0, 1)*GrFC(I0, I0, 2)
         enddo
         EgRDouOcc = EgRDouOcc / dble(NumNS)
!____________________________________________________________________________ 	  
!________________ [1] For OBC, calculate double occupancy at site 1 _________
!____________________________________________________________________________
      else
         EgRDouOcc = GrFC(1, 1, 1)*GrFC(1, 1, 2)
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) Expectation value of NN n_{i\sigma}*n_{j,-\sigma} ________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [0] For PBC, Calculate NN n_{i\up}*n_{j,down} _____________
!____________________________________________________________________________
      if(IfPyObsPBC) then
         EgRNNDCrF = 0.0_rp
         do I1 = 1, NumNS, +1
            I2 = FNNBond(I1, 2)
            EgRNNDCrF = EgRNNDCrF + GrFC(I1, I1, 1)*GrFC(I2, I2, 2)
            I2 = FNNBond(I1, 1)
            EgRNNDCrF = EgRNNDCrF + GrFC(I1, I1, 1)*GrFC(I2, I2, 2)
         enddo
         EgRNNDCrF = EgRNNDCrF / dble(2*NumNS)
!____________________________________________________________________________ 	  
!________________ [1] For OBC, Calculate NN n_{1\up}*n_{2,down} _____________
!____________________________________________________________________________
      else
         I1 = 1; I2 = 2
         EgRNNDCrF = GrFC(I1, I1, 1)*GrFC(I2, I2, 2)
      end if
!**************************************************************************************************	  
!___________________ 2. Pick up the constant and Accumulate all the results _______________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Energies, occupations and density correlations ___________________
!________________________________________________________________________________________
      !!!!!!!!!! Expectation energies of all terms in Hamiltonian
		EngOccCrF(01) = EngOccCrF(01) + EgHopt1 * CfgConst
		EngOccCrF(02) = EngOccCrF(02) + EgHopt2 * CfgConst
      EngOccCrF(03) = EngOccCrF(03) + EgHopt3 * CfgConst
      EngOccCrF(04) = EngOccCrF(04) + EgZmFld * CfgConst
      EngOccCrF(05) = EngOccCrF(05) + EgPinSz * CfgConst
      EngOccCrF(06) = EngOccCrF(06) + EgDopCh * CfgConst
      EngOccCrF(07) = EngOccCrF(07) + EgHubbU * CfgConst
      EngOccCrF(08) = EngOccCrF(08) + EgTotal * CfgConst
      !!!!!!!!!! First-order derivative of model parameters
      EngOccCrF(16) = EngOccCrF(16) +   EgHbUCh * CfgConst
      EngOccCrF(17) = EngOccCrF(17) +   EgTotCh * CfgConst
		EngOccCrF(18) = EngOccCrF(18) + HmOvHopt2 * CfgConst
      EngOccCrF(19) = EngOccCrF(19) + HmOvHopt3 * CfgConst
      EngOccCrF(20) = EngOccCrF(20) + HmOvZmFld * CfgConst
      EngOccCrF(21) = EngOccCrF(21) + HmOvHubbU * CfgConst

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      EngOccCrF(22) = EngOccCrF(22) + EgSinusoidalPinSz * CfgConst
      EngOccCrF(23) = EngOccCrF(23) + tmpHmOvSinusoidalPinSz * CfgConst
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !!!!!!!!!! Occupation, DouOcca and density correlations   
      EngOccCrF(31) = EngOccCrF(31) +                        EgnOccpUp  * CfgConst
      EngOccCrF(32) = EngOccCrF(32) +                        EgnOccpDw  * CfgConst
      EngOccCrF(33) = EngOccCrF(33) +             (EgnOccpUp+EgnOccpDw) * CfgConst
      EngOccCrF(34) = EngOccCrF(34) + dble(NumNS)*(EgnOccpUp+EgnOccpDw) * CfgConst
      EngOccCrF(35) = EngOccCrF(35) +                        EgRDouOcc  * CfgConst
      EngOccCrF(36) = EngOccCrF(36) +                        EgRNNDCrF  * CfgConst
      
   end subroutine ObStaEnrgy_SpnDcp
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine ObStaEnrgy_H0FFTW(CfgConst, GFr_MatC, GFk_AllC, EngOccCrF) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ObStaEnrgy_H0FFTW(CfgConst, GFr_MatC, GFk_AllC, EngOccCrF) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate energies using the diagonal elements of GrF matrices in r and k
!                 spaces.
! KEYWORDS: Calculate energies using diagonal elements of GrF matrices.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Calculate the expectation energies for all the terms in the model Hamiltonian.
!
!     Input:  CfgConst --> The configuration/walker related constant (like sign/phase or walker weight);
!             GFr_MatC --> Diagonal elements of r space GrF matrix, GFr_MatC(i, j) = <c_i^+ c_j>;
!             GFk_AllC --> The full matrix of k space GrF matrix,   GFk_AllC(k, q) = <c_k^+ c_q>.
!     
!     Output: EngOccCrF --> Accumulated measuring results of energies and occupations.
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
      real(rp)    GFr_MatC(NumNS, NumNS, 2)
      complex(rp) GFk_AllC(NumNS, NumNS, 2)
      real(rp) EngOccCrF(40)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, I3, I4, Ix, Iy
      real(rp) Rtp0
      real(rp) KxUp, KyUp, KxDn, KyDn
      real(rp) n_iup, n_idn, n_kup, n_kdn
      real(rp) HmOvHopt2, HmOvHopt3, HmOvZmFld, HmOvPinSz, HmOvHubbU, HmOvHbUCh
      real(rp) EgHopt1, EgHopt2, EgHopt3, EgZmFld, EgPinSz, EgDopCh, EgDopChTm, EgHubbU, EgTotal
      real(rp) EgHbUCh, EgTotCh
      real(rp) EgnOccpUp, EgnOccpDw, EgRDouOcc, EgRNNDCrF
      complex(rp) Ztp0
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Measure the expectation energies for the model ____________________________
!**************************************************************************************************              
!________________________________________________________________________________________ 	  
!_________________ (0) Energies for Hopping terms using GFk_AllC ________________________
!________________________________________________________________________________________
      EgHopt1 = 0.0_rp
      HmOvHopt2 = 0.0_rp; EgHopt2 = 0.0_rp
      HmOvHopt3 = 0.0_rp; EgHopt3 = 0.0_rp
      HmOvZmFld = 0.0_rp; EgZmFld = 0.0_rp
      EgDopChTm = 0.0_rp; EgDopCh = 0.0_rp
      do I0 = 1, NumNC, +1
         Ix = mod(I0-1, NumL1); Iy = (I0-1)/NumL1
         KxUp = dble(Ix) * 2.0_rp * rp_pi / dble(NumL1)
         KxDn = dble(Ix) * 2.0_rp * rp_pi / dble(NumL1)
         KyUp = dble(Iy) * 2.0_rp * rp_pi / dble(NumL2)
         KyDn = dble(Iy) * 2.0_rp * rp_pi / dble(NumL2)
         if(KxUp > rp_pi) KxUp = KxUp - 2.0_rp*rp_pi
         if(KxDn > rp_pi) KxDn = KxDn - 2.0_rp*rp_pi
         if(KyUp > rp_pi) KyUp = KyUp - 2.0_rp*rp_pi
         if(KyDn > rp_pi) KyDn = KyDn - 2.0_rp*rp_pi
         n_kup = real(GFk_AllC(I0, I0, 1))
         n_kdn = real(GFk_AllC(I0, I0, 2))
         if(EkDispType == 0) then
            Rtp0 = 2.0_rp * ( (cos(KxUp)+cos(KyUp))*Hopt1Up*n_kup + (cos(KxDn)+cos(KyDn))*Hopt1Dn*n_kdn )
         else if(EkDispType == 1) then
            Rtp0 = 2.0_rp * ( (cos(KxUp)+cos(KyUp)-2.0_rp)*Hopt1Up*n_kup + (cos(KxDn)+cos(KyDn)-2.0_rp)*Hopt1Dn*n_kdn )
         else if(EkDispType == 2) then
            Rtp0 = - Hopt1Up*n_kup*( KxUp*KxUp + KyUp*KyUp ) - Hopt1Dn*n_kdn*( KxDn*KxDn + KyDn*KyDn )
         end if
         EgHopt1   =   EgHopt1 + Rtp0
         HmOvHopt2 = HmOvHopt2 + 4.0_rp * ( cos(KxUp)*cos(KyUp)*n_kup + cos(KxDn)*cos(KyDn)*n_kdn )
         HmOvHopt3 = HmOvHopt3 + 2.0_rp * ( (cos(2.0_rp*KxUp)+cos(2.0_rp*KyUp))*n_kup + &
                                          & (cos(2.0_rp*KxDn)+cos(2.0_rp*KyDn))*n_kdn )
         HmOvZmFld = HmOvZmFld + (n_kup - n_kdn)
         if(EkDispType == 0) then
            Rtp0 = ChemP * (n_kup + n_kdn)
         else if(EkDispType >= 1) then
            Rtp0 = (ChemP + 4.0_rp*Hopt1Up)*n_kup + (ChemP + 4.0_rp*Hopt1Dn)*n_kdn
         end if
         EgDopChTm = EgDopChTm + Rtp0
      enddo
      EgHopt1   =   EgHopt1 / dble(NumNC)
      HmOvHopt2 = HmOvHopt2 / dble(NumNC)
      HmOvHopt3 = HmOvHopt3 / dble(NumNC)
      HmOvZmFld = HmOvZmFld / dble(NumNC)
      EgDopCh   = EgDopChTm / dble(NumNC)
      if(abs(Hopt2) > rp_Eps) EgHopt2 = HmOvHopt2 * Hopt2
      if(abs(Hopt3) > rp_Eps) EgHopt3 = HmOvHopt3 * Hopt3
      if(abs(ZmFdz) > rp_Eps) EgZmFld = HmOvZmFld * ZmFdz/2.0_rp
!________________________________________________________________________________________ 	  
!_________________ (1) Energies for pinning field term with PinSzType==2 ________________
!________________________________________________________________________________________
      HmOvPinSz = 0.0_rp; EgPinSz = 0.0_rp
      if( abs(PinSz) >= rp_Eps .and. PinSzType == 2 ) then
         Ztp0 = rp_Zzero
         do I1 = 1, NumNC, +1
            Ix = mod(I1-1, NumL1) + NumL1/2; Iy = (I1-1)/NumL1 + NumL2/2
            if(Ix >= NumL1) Ix = Ix - NumL1
            if(Iy >= NumL2) Iy = Iy - NumL2
            I2 = Iy*NumL1 + Ix + 1
            Ztp0 = Ztp0 + GFk_AllC(I1, I2, 1) - GFk_AllC(I1, I2, 2)
         enddo
         HmOvPinSz = dble(Ztp0) / dble(NumNS)
         EgPinSz = HmOvPinSz * PinSz/2.0_rp
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) Energies for HubbU terms using GFr_MatC __________________________
!________________________________________________________________________________________
      HmOvHbUCh = 0.0_rp; HmOvHubbU = 0.0_rp
      EgHbUCh   = 0.0_rp; EgHubbU   = 0.0_rp
      do I0 = 1, NumNS, +1
         n_iup = GFr_MatC(I0, I0, 1); n_idn = GFr_MatC(I0, I0, 2)
         Rtp0 = n_iup * n_idn
         HmOvHbUCh = HmOvHbUCh + Rtp0
         HmOvHubbU = HmOvHubbU + Rtp0 - (n_iup + n_idn)/2.0_rp
      enddo
      HmOvHbUCh = HmOvHbUCh / dble(NumNS); HmOvHubbU = HmOvHubbU / dble(NumNS)
      if(abs(HubbU) > rp_Eps) then
         EgHbUCh = HmOvHbUCh * HubbU
         EgHubbU = HmOvHubbU * HubbU
      end if
      EgRDouOcc = HmOvHbUCh
!________________________________________________________________________________________ 	  
!_________________ (3) Total energy _____________________________________________________
!________________________________________________________________________________________
      EgTotCh = EgHopt1 + EgHopt2 + EgHopt3 + EgZmFld + EgPinSz           + EgHbUCh
      EgTotal = EgHopt1 + EgHopt2 + EgHopt3 + EgZmFld + EgPinSz + EgDopCh + EgHubbU
!**************************************************************************************************	  
!___________________ 1. Measure Local densities and density-density Correlation ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Total occupation number of electrons _____________________________
!________________________________________________________________________________________
      EgnOccpUp = 0.0_rp; EgnOccpDw = 0.0_rp
      do I0 = 1, NumNS, +1
         EgnOccpUp = EgnOccpUp + GFr_MatC(I0, I0, 1)
         EgnOccpDw = EgnOccpDw + GFr_MatC(I0, I0, 2)
      enddo
      EgnOccpUp = EgnOccpUp / dble(NumNS)
      EgnOccpDw = EgnOccpDw / dble(NumNS)
!________________________________________________________________________________________ 	  
!_________________ (1) Expectation value of NN n_{i\sigma}*n_{j,-\sigma} ________________
!________________________________________________________________________________________
      EgRNNDCrF = 0.0_rp
      do I0 = 1, NumNS, +1
         I1 = FNNBond(I0, 1); I2 = FNNBond(I0, 2)
         EgRNNDCrF = EgRNNDCrF + GFr_MatC(I0, I0, 1) * ( GFr_MatC(I1, I1, 2) + GFr_MatC(I2, I2, 2) )
      enddo
      EgRNNDCrF = EgRNNDCrF / dble(2*NumNS)
!**************************************************************************************************	  
!___________________ 2. Pick up the constant and Accumulate all the results _______________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Energies, occupations and density correlations ___________________
!________________________________________________________________________________________
      !!!!!!!!!! Expectation energies of all terms in Hamiltonian
		EngOccCrF(01) = EngOccCrF(01) + EgHopt1 * CfgConst
		EngOccCrF(02) = EngOccCrF(02) + EgHopt2 * CfgConst
      EngOccCrF(03) = EngOccCrF(03) + EgHopt3 * CfgConst
      EngOccCrF(04) = EngOccCrF(04) + EgZmFld * CfgConst
      EngOccCrF(05) = EngOccCrF(05) + EgPinSz * CfgConst
      EngOccCrF(06) = EngOccCrF(06) + EgDopCh * CfgConst
      EngOccCrF(07) = EngOccCrF(07) + EgHubbU * CfgConst
      EngOccCrF(08) = EngOccCrF(08) + EgTotal * CfgConst
      !!!!!!!!!! First-order derivative of model parameters
      EngOccCrF(16) = EngOccCrF(16) +   EgHbUCh * CfgConst
      EngOccCrF(17) = EngOccCrF(17) +   EgTotCh * CfgConst
		EngOccCrF(18) = EngOccCrF(18) + HmOvHopt2 * CfgConst
      EngOccCrF(19) = EngOccCrF(19) + HmOvHopt3 * CfgConst
      EngOccCrF(20) = EngOccCrF(20) + HmOvZmFld * CfgConst
      EngOccCrF(21) = EngOccCrF(21) + HmOvHubbU * CfgConst
      !!!!!!!!!! Occupation, DouOcca and density correlations   
      EngOccCrF(31) = EngOccCrF(31) +                        EgnOccpUp  * CfgConst
      EngOccCrF(32) = EngOccCrF(32) +                        EgnOccpDw  * CfgConst
      EngOccCrF(33) = EngOccCrF(33) +             (EgnOccpUp+EgnOccpDw) * CfgConst
      EngOccCrF(34) = EngOccCrF(34) + dble(NumNS)*(EgnOccpUp+EgnOccpDw) * CfgConst
      EngOccCrF(35) = EngOccCrF(35) +                        EgRDouOcc  * CfgConst
      EngOccCrF(36) = EngOccCrF(36) +                        EgRNNDCrF  * CfgConst
      
   end subroutine ObStaEnrgy_H0FFTW
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$