!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the data process for the data of time-displaced measurements.
!             These real space correlation functions are preparing for the calculation of reciprocal space 
!             structure factors.
! COMMENT: Post process for time-displaced physical observables in every BIN simulation.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!    PostDynamc       --> Subroutine to calculate average values for the time-displaced observables;
!    FourierTrans_Tau --> Subroutine to perform Fourier Transformation for correlation functions.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine PostDynamc(NB) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PostDynamc(NB) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the post process for the time-displaced physical observables.
! KEYWORDS: Post Process of time-displaced observables.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Post Process of Dynamic observables.
!
!     Input: NB --> The iteration number of NmBin for the CPMC simulation;  
!
!     Outpt: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use Observable
            use TimeRecord
            use QMCTimeRec
            implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB           ! Number index of the BIN
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for dynamic data process __________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&       
            TimsDatPr = TimsDatPr + 1
            call system_clock(time1)
!______________________________________________________________________________________________________________     
!__________________ Main calculations for Average values of time-displaced observables ________________________
!______________________________________________________________________________________________________________   
!**************************************************************************************************     
!___________________ 0. Process dynamic correlations for every BIN simulation _____________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Results of dynamic correlation functions _________________________
!________________________________________________________________________________________
      call ProcRKDynCrFTau(NB)
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for dynamic data process __________________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&       
      call system_clock(time2)
      TimeDatPr = TimeDatPr + TimeIntrvl(time1, time2)
      
   end subroutine PostDynamc
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!#################### Process dynamic correlations for both OPEN and PERIODIC Boundary conditions #######################
!#################### Process dynamic correlations for both OPEN and PERIODIC Boundary conditions #######################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine ProcRKDynCrFTau(NB) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ProcRKDynCrFTau(NB) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the post process for the time-displaced physical observables, for
!                  PERIODIC boundary conditions.
! KEYWORDS: Post Process of time-displaced observables.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Post Process of observables.
!
!     Input: NB --> The iteration number of NmBin for the CPMC simulation;  
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
            use TimeRecord
            use QMCTimeRec
            use MPISetting
            implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
            integer NB
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
            integer I0, I1, I2, I3, I4, Kd, Hd, NTInd, Idimj, Nk, I0ItrMax
      integer IndZerodist, IndxFNNDist, IndxSNNDist, ITmpA(0:20)
      character(100) FileNameA
      complex(rp) RlCrFTmp(20)
      complex(rp), allocatable :: Collect0(:, :, :)
!______________________________________________________________________________________________________________     
!__________________ Main calculations for Average values of time-displaced observables ________________________
!______________________________________________________________________________________________________________   
!**************************************************************************************************     
!___________________ 0. Process dynamic correlations for PERIODIC boundary conditions _____________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Collect the results of the present BIN ___________________________
!________________________________________________________________________________________
      RealSpCrFTauBIN = RealSpCrFTauBIN / WtMeanSumBIN
!________________________________________________________________________________________         
!_________________ (1) Obtain the local dynamic correlation functions ___________________
!_____________________ Only for the IfPyObsPBC==T case __________________________________
!________________________________________________________________________________________
      if( (amyid == amstr) .and. (IfPyObsPBC) ) then
         !!!!!!!!!! Store local dynamic correlation functions
         IndZerodist = IminusJ(1, 1)
         IndxFNNDist = IminusJ(1, 2); IndxSNNDist = IminusJ(1, SNNBond(1, 1))
         do NTInd = 0, NumTauPnt, +1
            !!!!!!!! Local single-particle Green's function
            RSpCrFLclTauAll(NB, NTInd, 01) = RealSpCrFTauBIN(NTInd, IndZerodist, 01)
            RSpCrFLclTauAll(NB, NTInd, 02) = RealSpCrFTauBIN(NTInd, IndZerodist, 02)
            RSpCrFLclTauAll(NB, NTInd, 03) = RealSpCrFTauBIN(NTInd, IndZerodist, 03)
            RSpCrFLclTauAll(NB, NTInd, 04) = RealSpCrFTauBIN(NTInd, IndZerodist, 04)
            !!!!!!!! Local Sz-Sz and Density-density correlation function
            RSpCrFLclTauAll(NB, NTInd, 05) = RealSpCrFTauBIN(NTInd, IndZerodist, 05)
            RSpCrFLclTauAll(NB, NTInd, 06) = RealSpCrFTauBIN(NTInd, IndZerodist, 07)
            RSpCrFLclTauAll(NB, NTInd, 07) = RSpCrFLclTauAll(NB, NTInd, 06) &
                                                            & + RSpCrFLclTauAll(NB, NTInd, 05)
            RSpCrFLclTauAll(NB, NTInd, 08) = RSpCrFLclTauAll(NB, NTInd, 06) &
                                                            & - RSpCrFLclTauAll(NB, NTInd, 05)
            !!!!!!!! The FNN and SNN Sz-Sz correlations
            RSpCrFLclTauAll(NB, NTInd, 09) = RealSpCrFTauBIN(NTInd, IndxFNNDist, 05)
            RSpCrFLclTauAll(NB, NTInd, 10) = RealSpCrFTauBIN(NTInd, IndxSNNDist, 05)
         enddo
      end if
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%  
!_________________ (2) The Vertex Contribution for pairing correlations _________________
!__________________________ only for IfPyObsPBC == T case _______________________________
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
      if( (amyid == amstr) .and. (IfPyObsPBC) .and. (abs(PinSz) < rp_Eps) ) then
         call VertexContrbDyn()
      end if
!________________________________________________________________________________________         
!_________________ (3) Deduct the background for correlation functions __________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         call DeductCrFBkgDyn()
      end if
!**************************************************************************************************     
!___________________ 1. Postprocess of r- and k-space correlations for both _______________________
!_____________________________ PERIODIC and OPEN boundary conditions ______________________________
!**************************************************************************************************
      if(amyid == amstr) then
!________________________________________________________________________________________         
!_________________ (0) Obtain k-space correlation via Fourier Transformation ____________
!_____________________ for IfPyObsPBC == T case _________________________________________
!_____________________ from RealSpCrFTauBIN to KSpCrFtTauforFT __________________________
!________________________________________________________________________________________ 
         if(IfPyObsPBC) then
            call FourierTransTau()
         end if
!________________________________________________________________________________________         
!_________________ (1) For IfPyObsPBC == T case, store k-space dynamic __________________
!_____________________ correlations at specific k points ________________________________
!________________________________________________________________________________________
         if(IfPyObsPBC) then
!____________________________________________________________________________         
!________________ [0] Define the special k points in BZ _____________________
!____________________________________________________________________________
            ITmpA(0) = InvKpList(0, 0)                      ! The Gamma=(0,0) Point

            ITmpA(1) = InvKpList(+NumL1/2,        0)        ! The X1=(pi,0) point
            ITmpA(2) = InvKpList(       0, +NumL2/2)        ! The X2=(0,pi) point

            ITmpA(3) = InvKpList(+  NumL1/4, +  NumL2/4)    ! The Y1=(+pi/2, +pi/2) point
            ITmpA(4) = InvKpList(+  NumL1/4, +3*NumL2/4)    ! The Y2=(+pi/2, -pi/2) point
            ITmpA(5) = InvKpList(+3*NumL1/4, +  NumL2/4)    ! The Y3=(-pi/2, +pi/2) point
            ITmpA(6) = InvKpList(+3*NumL1/4, +3*NumL2/4)    ! The Y4=(-pi/2, -pi/2) point

            ITmpA(7) = InvKpList(+NumL1/2, +NumL2/2)        ! The M point

            ITmpA(8) = InvKpList(+1,  0)
            ITmpA(9) = InvKpList( 0, +1)
!____________________________________________________________________________         
!________________ [1] Obtain the data at special k points ___________________
!____________________________________________________________________________
            do NTInd = 0, NumTauPnt, +1
!____________________________________________________________________     
!___________________ GFUpUp at (pi,0) and (pi/2,pi/2) points ________
!____________________________________________________________________
               !!!!!!!!!! (pi,0) point
               RlCrFTmp(1:2) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  RlCrFTmp(1) = RlCrFTmp(1) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 01)
                  RlCrFTmp(2) = RlCrFTmp(2) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 02)
               enddo
               RlCrFTmp(1:2) = RlCrFTmp(1:2) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               !!!!!!!!!! (pi/2,pi/2) point
               RlCrFTmp(3:4) = rp_Zzero
               do I0 = 3, 6, +1
                  RlCrFTmp(3) = RlCrFTmp(3) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 01)
                  RlCrFTmp(4) = RlCrFTmp(4) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 02)
               enddo
               RlCrFTmp(3:4) = RlCrFTmp(3:4) / 4.0_rp
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:2, 01) = real( RlCrFTmp(1:3:+2) )
               RorKspCrFTauAll(NB, NTInd, 1:2, 02) = real( RlCrFTmp(2:4:+2) )
!____________________________________________________________________     
!___________________ GFDnDn at (pi,0) and (pi/2,pi/2) points ________
!____________________________________________________________________ 
               !!!!!!!!!! (pi,0) point
               RlCrFTmp(1:2) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  RlCrFTmp(1) = RlCrFTmp(1) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 03)
                  RlCrFTmp(2) = RlCrFTmp(2) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 04)
               enddo
               RlCrFTmp(1:2) = RlCrFTmp(1:2) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               !!!!!!!!!! (pi/2,pi/2) point
               RlCrFTmp(3:4) = rp_Zzero
               do I0 = 3, 6, +1
                  RlCrFTmp(3) = RlCrFTmp(3) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 03)
                  RlCrFTmp(4) = RlCrFTmp(4) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 04)
               enddo
               RlCrFTmp(3:4) = RlCrFTmp(3:4) / 4.0_rp
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:2, 03) = real( RlCrFTmp(1:3:+2) )
               RorKspCrFTauAll(NB, NTInd, 1:2, 04) = real( RlCrFTmp(2:4:+2) )
!____________________________________________________________________     
!___________________ SpinZZ at Gamma, X and M points ________________
!____________________________________________________________________ 
               !!!!!!!!!! The Gamma point
               RlCrFTmp(1) = KSpCrFtTauforFT(NTInd, ITmpA(0), 05)
               !!!!!!!!!! The X point
               RlCrFTmp(2) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  RlCrFTmp(2) = RlCrFTmp(2) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 05)
               enddo
               RlCrFTmp(2) = RlCrFTmp(2) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               !!!!!!!!!! The M point
               RlCrFTmp(3) = KSpCrFtTauforFT(NTInd, ITmpA(7), 05)
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:3, 05) = real( RlCrFTmp(1:3) )
!____________________________________________________________________     
!___________________ SpinPM at Gamma, X and M points ________________
!____________________________________________________________________
               !!!!!!!!!! The Gamma point
               RlCrFTmp(1) =   KSpCrFtTauforFT(NTInd, ITmpA(0), 06)
               !!!!!!!!!! The X point
               RlCrFTmp(2) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  RlCrFTmp(2) = RlCrFTmp(2) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 06)
               enddo
               RlCrFTmp(2) = RlCrFTmp(2) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               !!!!!!!!!! The M point
               RlCrFTmp(3) = KSpCrFtTauforFT(NTInd, ITmpA(7), 06)
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:3, 06) = real( RlCrFTmp(1:3) )
!____________________________________________________________________     
!___________________ DenDen at Gamma, X and M points ________________
!____________________________________________________________________ 
               !!!!!!!!!! The Gamma point
               RlCrFTmp(1) = KSpCrFtTauforFT(NTInd, ITmpA(0), 07)
               !!!!!!!!!! The X point
               RlCrFTmp(2) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  RlCrFTmp(2) = RlCrFTmp(2) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 07)
               enddo
               RlCrFTmp(2) = RlCrFTmp(2) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               !!!!!!!!!! The M point
               RlCrFTmp(3) = KSpCrFtTauforFT(NTInd, ITmpA(7), 07)
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:3, 07) = real( RlCrFTmp(1:3) )
!____________________________________________________________________     
!___________________ Currnt at Gamma, +qx, +qy points _______________
!____________________________________________________________________
               !!!!!!!!!! Gamma, +qx, +qy points
               RlCrFTmp(1) = KSpCrFtTauforFT(NTInd, ITmpA(0), 08)
               RlCrFTmp(2) = KSpCrFtTauforFT(NTInd, ITmpA(8), 08)  
               RlCrFTmp(3) = KSpCrFtTauforFT(NTInd, ITmpA(9), 08)
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:3, 08) = real( RlCrFTmp(1:3) )
!____________________________________________________________________     
!___________________ PairSt at Gamma, X and M points ________________
!____________________________________________________________________ 
               !!!!!!!!!! The Gamma point
               RlCrFTmp(1) = KSpCrFtTauforFT(NTInd, ITmpA(0), 09)
               !!!!!!!!!! The X point
               RlCrFTmp(2) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  RlCrFTmp(2) = RlCrFTmp(2) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 09)
               enddo
               RlCrFTmp(2) = RlCrFTmp(2) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               !!!!!!!!!! The M point 
               RlCrFTmp(3) = KSpCrFtTauforFT(NTInd, ITmpA(7), 09)
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:3, 09) = real( RlCrFTmp(1:3) )
!____________________________________________________________________     
!___________________ EdSPar at Gamma, X and M points ________________
!____________________________________________________________________ 
               !!!!!!!!!! The Gamma point
               RlCrFTmp(1) = KSpCrFtTauforFT(NTInd, ITmpA(0), 10)
               !!!!!!!!!! The X point
               RlCrFTmp(2) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  RlCrFTmp(2) = RlCrFTmp(2) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 10)
               enddo
               RlCrFTmp(2) = RlCrFTmp(2) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               !!!!!!!!!! The M point
               RlCrFTmp(3) = KSpCrFtTauforFT(NTInd, ITmpA(7), 10)
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:3, 10) = real( RlCrFTmp(1:3) )
!____________________________________________________________________     
!___________________ DWvPar at Gamma, X and M points ________________
!____________________________________________________________________ 
               !!!!!!!!!! The Gamma point
               RlCrFTmp(1) = KSpCrFtTauforFT(NTInd, ITmpA(0), 11)
               !!!!!!!!!! The X point
               RlCrFTmp(2) = rp_Zzero
               do I0 = 1, merge(2, 1, NumL1==NumL2), +1
                  RlCrFTmp(2) = RlCrFTmp(2) + KSpCrFtTauforFT(NTInd, ITmpA(I0), 11)
               enddo
               RlCrFTmp(2) = RlCrFTmp(2) / merge(2.0_rp, 1.0_rp, NumL1==NumL2)
               !!!!!!!!!! The M point 
               RlCrFTmp(3) = KSpCrFtTauforFT(NTInd, ITmpA(7), 11)
               !!!!!!!!!! Store the results
               RorKspCrFTauAll(NB, NTInd, 1:3, 11) = real( RlCrFTmp(1:3) )
            enddo
         end if
!________________________________________________________________________________________         
!_________________ (2) For IfPyObsPBC == F case, store r-space dynamic __________________
!_____________________ correlations for all compute site pairs __________________________
!________________________________________________________________________________________
         if(.not. IfPyObsPBC) then
            do I0 = 1, 19, +1 
               do Idimj = 1, NmSitePairTau, +1
                  do NTInd = 0, NumTauPnt, +1
                     RorKspCrFTauAll(NB, NTInd, Idimj, I0) = RealSpCrFTauBIN(NTInd, Idimj, I0)
                  enddo
               enddo
            enddo
         end if
      end if
!**************************************************************************************************     
!___________________ 2. Output the results for the present BIN simulation _________________________
!______________________ for both PERIODIC and OPEN boundary conditions ____________________________
!**************************************************************************************************
      if(amyid == amstr) then
!________________________________________________________________________________________         
!_________________ (0) Output dynamic correlations at all k points ______________________
!_____________________ only for IfPyObsPBC == T case ____________________________________
!________________________________________________________________________________________
         if(IfPyObsPBC) then
            !!!!!!!!!! The dynamic Single particle Green's Function
            if(IfDyGrFOut) then
               do Nk = 1, NumNC, +1
                  write(FileNameA, "('Add_Output/UU_GreenFDynmcData/GrnFct_', I2.2, '_', I2.2, '.txt')") &
                     & KpList(Nk, 1), KpList(Nk, 2)
                  open(291, file = trim(FileNameA), access = "append")
                  do NTInd = 0, NumTauPnt, +1
                     write(291, "(es17.8)", advance = "no") TauPntVal(NTInd)*Dltau
                     write(291, "(A, es17.8, A, es17.8, A, es17.8, A, es17.8)", advance = "no") &
                        & char(9), real(KSpCrFtTauforFT(NTInd, Nk, 01)), char(9), real(KSpCrFtTauforFT(NTInd, Nk, 02)), &
                        & char(9), real(KSpCrFtTauforFT(NTInd, Nk, 03)), char(9), real(KSpCrFtTauforFT(NTInd, Nk, 04))
                     write(291, "()")
                  enddo
                  close(291)
               enddo
            end if
            !!!!!!!!!! The dynamic spin-spin correlation Function
            if(IfDySpnOut) then
               do Nk = 1, NumNC, +1
                  write(FileNameA, "('Add_Output/VV_SpnSpnDynmcData/SpnSpn_', I2.2, '_', I2.2, '.txt')") &
                     & KpList(Nk, 1), KpList(Nk, 2)
                  open(291, file = trim(FileNameA), access = "append")
                  do NTInd = 0, NumTauPnt, +1
                     write(291, "(es17.8)", advance = "no") TauPntVal(NTInd)*Dltau
                     write(291, "(A, es17.8, A, es17.8)", advance = "no") &
                        & char(9), real(KSpCrFtTauforFT(NTInd, Nk, 05)), char(9), real(KSpCrFtTauforFT(NTInd, Nk, 06))
                     write(291, "()")
                  enddo
                  close(291)
               enddo
            end if
            !!!!!!!!!! The dynamic density-density correlation Function
            if(IfDyDenOut) then
               do Nk = 1, NumNC, +1
                  write(FileNameA, "('Add_Output/WW_DenDenDynmcData/DenDen_', I2.2, '_', I2.2, '.txt')") &
                     & KpList(Nk, 1), KpList(Nk, 2)
                  open(291, file = trim(FileNameA), access = "append")
                  do NTInd = 0, NumTauPnt, +1
                     write(291, "(es17.8)", advance = "no") TauPntVal(NTInd)*Dltau
                     write(291, "(A, es17.8)", advance = "no") char(9), real(KSpCrFtTauforFT(NTInd, Nk, 07))
                     write(291, "()")
                  enddo
                  close(291)
               enddo
            end if
            !!!!!!!!!! The dynamic Current-Current correlation Function
            if(IfDyCurrnt) then
               do Nk = 1, NumNC, +1
                  write(FileNameA, "('Add_Output/ZZ_CurrntDynmcData/Currnt_', I2.2, '_', I2.2, '.txt')") &
                     & KpList(Nk, 1), KpList(Nk, 2)
                  open(291, file = trim(FileNameA), access = "append")
                  do NTInd = 0, NumTauPnt, +1
                     write(291, "(es17.8)", advance = "no") TauPntVal(NTInd)*Dltau
                     write(291, "(A, es17.8)", advance = "no") char(9), real(KSpCrFtTauforFT(NTInd, Nk, 08))
                     write(291, "()")
                  enddo
                  close(291)
               enddo
            end if
            !!!!!!!!!! The dynamic PairSt-PairSt correlation Function    
            if(IfDyPstOut) then
               do Nk = 1, NumNC, +1
                  write(FileNameA, "('Add_Output/XX_PairStDynmcData/PairSt_', I2.2, '_', I2.2, '.txt')") &
                     & KpList(Nk, 1), KpList(Nk, 2)
                  open(291, file = trim(FileNameA), access = "append")
                  do NTInd = 0, NumTauPnt, +1
                     write(291, "(es17.8)", advance = "no") TauPntVal(NTInd)*Dltau
                     write(291, "(A, es17.8)", advance = "no") char(9), real(KSpCrFtTauforFT(NTInd, Nk, 09))
                     write(291, "()")
                  enddo
                  close(291)
               enddo
            end if    
            !!!!!!!!!! The dynamic DWvPar-DWvPar correlation Function  
            if(IfDyDWvOut) then
               do Nk = 1, NumNC, +1
                  write(FileNameA, "('Add_Output/YY_DWvParDynmcData/DWvPar_', I2.2, '_', I2.2, '.txt')") &
                     & KpList(Nk, 1), KpList(Nk, 2)
                  open(291, file = trim(FileNameA), access = "append")
                  do NTInd = 0, NumTauPnt, +1
                     write(291, "(es17.8)", advance = "no") TauPntVal(NTInd)*Dltau
                     write(291, "(A, es17.8, A, es17.8)", advance = "no") &
                        & char(9), real(KSpCrFtTauforFT(NTInd, Nk, 10)), char(9), real(KSpCrFtTauforFT(NTInd, Nk, 11))
                     write(291, "()")
                  enddo
                  close(291)
               enddo
            end if
         end if
!________________________________________________________________________________________         
!_________________ (1) Output local dynamic correlations ________________________________
!_____________________ only for IfPyObsPBC == T case ____________________________________
!________________________________________________________________________________________
         if(IfPyObsPBC) then
            !!!!!!!!!! The local correlation functions versus. \tau
            open(342, file = "Add_Output/04_LocalGrnFctcF_Tau.txt", access = "append")
            open(343, file = "Add_Output/04_LocalSpzDencF_Tau.txt", access = "append")
            open(344, file = "Add_Output/04_FNNSNNSpinCrF_Tau.txt", access = "append")
            do NTInd = 0, NumTauPnt, +1
               !!!!!!!! Local single-particle Green's function
               write(342, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
               write(342, "(A, es17.8, A, es17.8, A, es17.8, A, es17.8)", advance = "no") &
                  & char(9), RSpCrFLclTauAll(NB, NTInd, 01), char(9), RSpCrFLclTauAll(NB, NTInd, 02), &
                  & char(9), RSpCrFLclTauAll(NB, NTInd, 03), char(9), RSpCrFLclTauAll(NB, NTInd, 04)
               write(342, "()")
               !!!!!!!! Local Sz-Sz and Density-density correlation function
               write(343, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
               write(343, "(A, es17.8, A, es17.8, A, es17.8, A, es17.8)", advance = "no") &
                  & char(9), RSpCrFLclTauAll(NB, NTInd, 05), char(9), RSpCrFLclTauAll(NB, NTInd, 06), &
                  & char(9), RSpCrFLclTauAll(NB, NTInd, 07), char(9), RSpCrFLclTauAll(NB, NTInd, 08)
               write(343, "()")
               !!!!!!!! The Sz-Sz correlations for FNN and SNN
               write(344, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
               write(344, "(A, es17.8, A, es17.8)", advance = "no") &
                  & char(9), RSpCrFLclTauAll(NB, NTInd, 09), char(9), RSpCrFLclTauAll(NB, NTInd, 10)
               write(344, "()")
            enddo
            close(342); close(343); close(344)
            !!!!!!!!!! \beta * [ S_{13}(\beta/2) + S_{12}(\beta/2) ]
            if( NmTDM >= LTrot/2 ) then
               NTInd = TauHalfBetaT
               Spn13pS12HalfBt(NB) = BetaT * ( RSpCrFLclTauAll(NB, NTInd, 10) + RSpCrFLclTauAll(NB, NTInd, 09) )
               open( 345, file = "Add_Output/S13pS12HfBBINS.txt", access = "append")
               write(345, "(I4.4, A, es17.8)") NB, char(9), Spn13pS12HalfBt(NB)
               close(345)
            end if
         end if
!________________________________________________________________________________________         
!_________________ (2) IfPyObsPBC==T --> k-point Dynamic correlations ___________________
!_____________________ IfPyObsPBC==F --> r-space Dynamic correlations ___________________
!________________________________________________________________________________________
         open(341, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "GFUpUpCorFun_Tau.txt", access = "append")
         open(342, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "GFDnDnCorFun_Tau.txt", access = "append")
         open(343, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "SpinZZCorFun_Tau.txt", access = "append") 
         open(344, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "SpinPMCorFun_Tau.txt", access = "append")
         open(345, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "DenDenCorFun_Tau.txt", access = "append")
         open(346, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "CurrntCorFun_Tau.txt", access = "append")
         open(347, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "PairStCorFun_Tau.txt", access = "append")
         open(348, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "EdSParCorFun_Tau.txt", access = "append")
         open(349, file = "Add_Output/04_" // merge("K", "R", IfPyObsPBC) // "DWvParCorFun_Tau.txt", access = "append")
         do NTInd = 0, NumTauPnt, +1
            !!!!!!!!!! The dynamic single-particle Green's functions
            do I0 = 1, 2, +1
               write(340+I0, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
               do Idimj = 1, merge(2, NmSitePairTau, IfPyObsPBC), +1
                  write(340+I0, "(A, es17.8, A, es17.8)", advance = "no") char(9), &
                     & RorKspCrFTauAll(NB, NTInd, Idimj, 2*I0-1), char(9), RorKspCrFTauAll(NB, NTInd, Idimj, 2*I0)
               enddo
            enddo
            !!!!!!!!!! The dynamic correlation functions in bosonic channels
            do I0 = 3, 9, +1
               write(340+I0, "(es17.8, A)", advance = "no") TauPntVal(NTInd)*Dltau, char(9)
               do Idimj = 1, merge(3, NmSitePairTau, IfPyObsPBC), +1
                  write(340+I0, "(A, es17.8)", advance = "no") char(9), RorKspCrFTauAll(NB, NTInd, Idimj, 2+I0)
               enddo
            enddo
            !!!!!!!!!! Write to a new line
            do I0 = 1, 9, +1
               write(340+I0, "()")
            enddo
         enddo
         close(341); close(342); close(343); close(344); close(345); close(346); close(347); close(348); close(349)
!________________________________________________________________________________________         
!_________________ (3) Store the G(k, Beta/2) for fermi surface _________________________
!_____________________ Store Current(k=0, Beta/2) as dc conductivity ____________________
!_____________________ Only for IfPyObsPBC == T case ____________________________________
!________________________________________________________________________________________
         if( (IfPyObsPBC) .and. (NmTDM >= LTrot/2) ) then
            !!!!!!!!!! Approximated spectral function A(w=e_F, k) = \beta * G(k, Beta/2)
            NTInd = TauHalfBetaT
            open(291, file = "Add_Output/GkHalfBetaBINs.txt", access = "append") 
            do Nk = 1, NumNC, +1
               RlCrFTmp(1) = ( KSpCrFtTauforFT(NTInd, Nk, 1) + KSpCrFtTauforFT(NTInd, Nk, 2) ) / 2.0_rp
               RlCrFTmp(2) = ( KSpCrFtTauforFT(NTInd, Nk, 3) + KSpCrFtTauforFT(NTInd, Nk, 4) ) / 2.0_rp
               GkTauBetaOv2All(NB, Nk, 1) = real(RlCrFTmp(1)) * BetaT
               GkTauBetaOv2All(NB, Nk, 2) = real(RlCrFTmp(2)) * BetaT
               write(291, "(I4, A, I4)", advance = "no") KpList(Nk, 1), char(9), KpList(Nk, 2)
               write(291, "(A, es17.8, A, es17.8)") char(9), GkTauBetaOv2All(NB, Nk, 1), &
                                                  & char(9), GkTauBetaOv2All(NB, Nk, 2)
            enddo
            close(291)
            !!!!!!!!!! Approximated dc conductivity \sigma_{dc} = \beta^2/\pi * Current(k=0, Beta/2)
            NTInd = TauHalfBetaT
            Dc_Conductivity(NB) = RorKspCrFTauAll(NB, NTInd, 1, 08) * BetaT*BetaT/rp_pi
            open( 291, file = "Add_Output/HdcConductBINs.txt", access = "append") 
            write(291, "(I4.4, A, es17.8)") NB, char(9), Dc_Conductivity(NB)
            close(291)
         end if
      end if
!**************************************************************************************************     
!___________________ 3. For IfPyObsPBC==T case, Compute Matsubara Frequency quantities ____________
!**************************************************************************************************
      if( (amyid == amstr) .and. (IfPyObsPBC) .and. (NmTDMType == 0 .or. NmTDMType == 1) ) then
         call ComptGrFCrFkIwn(NB)
      end if
      
   end subroutine ProcRKDynCrFTau
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine VertexContrbDyn()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  VertexContrbDyn() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compute the vertex contribution for various bosonic correlation functions. 
! KEYWORDS: Compute the vertex contribution. For the dynamic correlation functions.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Fourier Transformation.
!
!     Input: (none); Outpt: (none).
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
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
            integer Idimj, I0, I1, I2, I3, I4, I0ItrMax, Kd, Hd, NTInd
            real(rp) RlCrFTmp(15)
      real(rp), allocatable :: Collect0(:, :, :)
!______________________________________________________________________________________________________________     
!_____________________________ Compute the vertex contribution for correlations _______________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Vertex contribution = C_{QMC} - C_{uncorrelated} __________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Allocate the temporary matrice ___________________________________
!________________________________________________________________________________________ 
      allocate(Collect0(0:NumTauPnt, NumNC, 15)); Collect0 = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) Compute the uncorrelated part of the correlations ________________
!________________________________________________________________________________________
      !!!!!!!!!! First compute the uncorrelated pairing correlations
      I0ItrMax = NumNC * NumNC
   !$OMP PARALLEL &
   !$OMP PRIVATE(I0, I1, I2, I3, I4, Kd, Hd, Idimj, NTInd, RlCrFTmp)
   !$OMP DO REDUCTION(+ : Collect0)
      do I0 = 1, I0ItrMax, +1
         !!!!!!!! The integer index for the present term
         I1 = (I0-1)/NumNC + 1; I2 = mod(I0-1, NumNC) + 1
         Idimj = IminusJ(I1, I2)
         !!!!!!!! Iteration for all Tau points
         do NTInd = 0, NumTauPnt, +1
            !!!!!! Spin correlations, <SzSz> and <S+S- + S-S+>/2
            !!!! The <SzSz> correlation
            RlCrFTmp(1) =   + ( RlSpGrnFtBIN(I1, I1, 1, 1) - RlSpGrnFtBIN(I1, I1, 2, 1) ) * &
                                       & ( RlSpGrnFtBIN(I2, I2, 1, 1) - RlSpGrnFtBIN(I2, I2, 2, 1) ) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 02) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 01) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 04) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 03)
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            Collect0(NTInd, Idimj, 01) = Collect0(NTInd, Idimj, 01) + RlCrFTmp(1)
            !!!! The <S+S- + S-S+>/2 correlation
            RlCrFTmp(2) =   + RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 02) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 03) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 04) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 01)
            RlCrFTmp(2) = RlCrFTmp(2) / 2.0_rp
            Collect0(NTInd, Idimj, 02) = Collect0(NTInd, Idimj, 02) + RlCrFTmp(2)
            !!!!!! Density-density correlation function
            RlCrFTmp(1) =   + ( 2.0_rp - RlSpGrnFtBIN(I1, I1, 1, 1) - RlSpGrnFtBIN(I1, I1, 2, 1) ) * &
                                       & ( 2.0_rp - RlSpGrnFtBIN(I2, I2, 1, 1) - RlSpGrnFtBIN(I2, I2, 2, 1) ) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 02) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 01) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 04) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 03)
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            Collect0(NTInd, Idimj, 03) = Collect0(NTInd, Idimj, 03) + RlCrFTmp(1)
            !!!!!! Current-Current correlations in x-direction
            I3 = FNNBond(I1, 2); I4 = FNNBond(I2, 2)
            RlCrFTmp(1) =   ( RlSpGrnFtBIN(I1, I3, 1, 1) - RlSpGrnFtBIN(I3, I1, 1, 1) + &
                                       & RlSpGrnFtBIN(I1, I3, 2, 1) - RlSpGrnFtBIN(I3, I1, 2, 1) ) * &
                          & ( RlSpGrnFtBIN(I2, I4, 1, 1) - RlSpGrnFtBIN(I4, I2, 1, 1) + &
                                       & RlSpGrnFtBIN(I2, I4, 2, 1) - RlSpGrnFtBIN(I4, I2, 2, 1) )
            RlCrFTmp(1) = RlCrFTmp(1) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I4, I1), 02) * RealSpCrFTauBIN(NTInd, IminusJ(I3, I2), 01) &
                          & - RealSpCrFTauBIN(NTInd, IminusJ(I4, I3), 02) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 01) &
                          & - RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 02) * RealSpCrFTauBIN(NTInd, IminusJ(I3, I4), 01) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I3), 02) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I4), 01) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I4, I1), 04) * RealSpCrFTauBIN(NTInd, IminusJ(I3, I2), 03) &
                          & - RealSpCrFTauBIN(NTInd, IminusJ(I4, I3), 04) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 03) &
                          & - RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 04) * RealSpCrFTauBIN(NTInd, IminusJ(I3, I4), 03) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I3), 04) * RealSpCrFTauBIN(NTInd, IminusJ(I1, I4), 03)
            Collect0(NTInd, Idimj, 04) = Collect0(NTInd, Idimj, 04) - RlCrFTmp(1)
            !!!!!! On-site spin-singlet s-wave pairing channel
            RlCrFTmp(1) =   + RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 01)*RealSpCrFTauBIN(NTInd, IminusJ(I1, I2), 03) &
                          & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 02)*RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 04)
            RlCrFTmp(1) = RlCrFTmp(1) / 4.0_rp
            Collect0(NTInd, Idimj, 05) = Collect0(NTInd, Idimj, 05) + RlCrFTmp(1)
            !!!!!! Various off-site pairing correlation functions
            RlCrFTmp = 0.0_rp
            do Kd = 1, 4, +1
               I3 = FNNBond(I1, Kd)
               do Hd = 1, 4, +1
                  I4 = FNNBond(I2, Hd)
                  !!!! The NN spin-singlet pairing channels
                  RlCrFTmp(15) = RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 02)*RealSpCrFTauBIN(NTInd, IminusJ(I4, I3), 04) &
                             & + RealSpCrFTauBIN(NTInd, IminusJ(I2, I3), 02)*RealSpCrFTauBIN(NTInd, IminusJ(I4, I1), 04) &
                             & + RealSpCrFTauBIN(NTInd, IminusJ(I4, I1), 02)*RealSpCrFTauBIN(NTInd, IminusJ(I2, I3), 04) &
                             & + RealSpCrFTauBIN(NTInd, IminusJ(I4, I3), 02)*RealSpCrFTauBIN(NTInd, IminusJ(I2, I1), 04)
                  RlCrFTmp(1) = RlCrFTmp(1) + RlCrFTmp(15)
                  RlCrFTmp(2) = RlCrFTmp(2) + dble(1-2*mod(Kd, 2)) * dble(1-2*mod(Hd, 2)) * RlCrFTmp(15)
               enddo
            enddo
            !!!! The NN spin-singlet pairing channels
            RlCrFTmp(1) = RlCrFTmp(1) / 32.0_rp
            Collect0(NTInd, Idimj, 06) = Collect0(NTInd, Idimj, 06) + RlCrFTmp(1)
            RlCrFTmp(2) = RlCrFTmp(2) / 32.0_rp
            Collect0(NTInd, Idimj, 07) = Collect0(NTInd, Idimj, 07) + RlCrFTmp(2)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
      !!!!!!!!!! Constant related to applying the translational symmetry
      Collect0 = Collect0 / dble(NumNC)
!________________________________________________________________________________________         
!_________________ (2) Compute the vertex contributions _________________________________
!________________________________________________________________________________________
      do Idimj = 1, NumNC, +1
         do NTInd = 0, NumTauPnt, +1
            RealSpCrFTauBIN(NTInd, Idimj, 20:34) = RealSpCrFTauBIN(NTInd, Idimj, 05:19) - Collect0(NTInd, Idimj, 01:15)
         enddo
      enddo
!________________________________________________________________________________________         
!_________________ (3) Deallocate the temporary matrix __________________________________
!________________________________________________________________________________________
      if(allocated(Collect0)) deallocate(Collect0)

   end subroutine VertexContrbDyn
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine DeductCrFBkgDyn()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  DeductCrFBkgDyn() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to deduct the background for the correlation functions in bosonic channels using
!                   RlSpGrnFtBIN and RealSpCrFTauBIN results. For the dynamic correlation functions.
! KEYWORDS: Deduct the background for correlations.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Fourier Transformation.
!
!     Input: (none); Outpt: (none).
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
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
            integer I0, I1, I2, I3, I4, I0ItrMax, Idimj, NTInd
            real(rp) RlCrFTmp(15)
      real(rp), allocatable :: Collect0(:, :, :)
!______________________________________________________________________________________________________________     
!_____________________________ Deduct the background for the correlation functions ____________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Computation for Deducting the background for correlations _________________
!************************************************************************************************** 
!________________________________________________________________________________________         
!_________________ (0) Allocate the temporary matrice ___________________________________
!________________________________________________________________________________________ 
      allocate(Collect0(0:NumTauPnt, NmSitePairTau, 5)); Collect0 = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) For IfCrfDfBkg==F case, reset RealSpCrFTauBIN ____________________
!________________________________________________________________________________________ 
      if(.not. IfCrfDfBkg) then
         RealSpCrFTauBIN(0:NumTauPnt, 1:NmSitePairTau, 35:36) = 0.0_rp
         RealSpCrFTauBIN(0:NumTauPnt, 1:NmSitePairTau, 37:38) = 0.0_rp
         RealSpCrFTauBIN(0:NumTauPnt, 1:NmSitePairTau, 39:40) = 0.0_rp
      end if
!________________________________________________________________________________________         
!_________________ (2) Compute the background for correlations __________________________
!________________________________________________________________________________________
      !!!!!!!!!! Compute the background for correlation functions
      I0ItrMax = merge(NumNC*NumNC, NmSitePairTau, IfPyObsPBC)
   !$OMP PARALLEL &
   !$OMP PRIVATE(I0, I1, I2, Idimj, NTInd, RlCrFTmp)
   !$OMP DO REDUCTION(+ : Collect0)
      do I0 = 1, I0ItrMax, +1
         !!!!!!!! The integer index for the present term
         I1 = merge(   (I0-1)/NumNC  + 1, I0, IfPyObsPBC)
         I2 = merge(mod(I0-1, NumNC) + 1, I0, IfPyObsPBC)
         Idimj = merge(IminusJ(I1, I2), I0, IfPyObsPBC)
         !!!!!!!! Iteration for all Tau points
         do NTInd = 0, NumTauPnt, +1
            RlCrFTmp(1) = RealSpCrFTauBIN(NTInd, I1, 35) * RealSpCrFTauBIN(NTInd, I2, 36) ! Sz-Sz
            RlCrFTmp(2) = RealSpCrFTauBIN(NTInd, I1, 37) * RealSpCrFTauBIN(NTInd, I2, 38) ! density-density
            RlCrFTmp(3) = RealSpCrFTauBIN(NTInd, I1, 39) * RealSpCrFTauBIN(NTInd, I2, 40) ! current-current
            Collect0(NTInd, Idimj, 1) = Collect0(NTInd, Idimj, 1) + RlCrFTmp(1)
            Collect0(NTInd, Idimj, 2) = Collect0(NTInd, Idimj, 2) + RlCrFTmp(2)
            Collect0(NTInd, Idimj, 3) = Collect0(NTInd, Idimj, 3) - RlCrFTmp(3)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
      !!!!!!!!!! The size-related constant for IfPyObsPBC==T case
      if(IfPyObsPBC) Collect0 = Collect0 / dble(NumNC)
!________________________________________________________________________________________         
!_________________ (2) Deduct the background in the correlations ________________________
!________________________________________________________________________________________
      do Idimj = 1, NmSitePairTau, +1
         do NTInd = 0, NumTauPnt, +1
            RealSpCrFTauBIN(NTInd, Idimj, 05) = RealSpCrFTauBIN(NTInd, Idimj, 05) - Collect0(NTInd, Idimj, 1)
            RealSpCrFTauBIN(NTInd, Idimj, 07) = RealSpCrFTauBIN(NTInd, Idimj, 07) - Collect0(NTInd, Idimj, 2)
            RealSpCrFTauBIN(NTInd, Idimj, 08) = RealSpCrFTauBIN(NTInd, Idimj, 08) - Collect0(NTInd, Idimj, 3)
         enddo
      enddo
!________________________________________________________________________________________         
!_________________ (3) Deallocate the temporary matrix __________________________________
!________________________________________________________________________________________
      if(allocated(Collect0)) deallocate(Collect0)

   end subroutine DeductCrFBkgDyn
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine FourierTransTau() 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  FourierTransTau() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to Perform Fourier transformation for GreenT_Tau, SpinZZ_Tau, SpinPM_Tau, 
!                   DenDen_Tau into reciprocal space.
! KEYWORDS: Fourier transformation for time-displaced quantities.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-10
! DESCRIPTION: Fourier Transformation.
!
!     Input: (none).    Outpt: (none).
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
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Nk, NTInd, imj, I1, I2, I0
      real(rp) Rtp1, VecK(2), VecR(2)
!______________________________________________________________________________________________________________     
!_____________________________ Main calculations for Fourier Transformation ___________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Fourier Transformation for all the data ___________________________________
!************************************************************************************************** 
!________________________________________________________________________________________         
!_________________ (0) Initialization of KSpCrFtTauforFT for storing results ____________
!________________________________________________________________________________________          
      KSpCrFtTauforFT = rp_Zzero
!________________________________________________________________________________________         
!_________________ (1) Data process for all NTInd _______________________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(Nk, NTInd, VecK, imj, I1, I2, I0, VecR, Rtp1)
   !$OMP DO
      do Nk = 1, NumNC, +1
         VecK = KpList(Nk, 1) * B1Vec + KpList(Nk, 2) * B2Vec
         do imj = 1, NumNC, +1
            I1 = StList(imj, 1)
            I2 = StList(imj, 2)
            I1 = I1 - NumL1/2 - mod(NumL1, 2)
            I2 = I2 - NumL2/2 - mod(NumL2, 2)
            VecR = I1 * A1Vec + I2 * A2Vec
            Rtp1 = dot_product(VecK, VecR)
            do I0 = 1, 34, +1
               do NTInd = 0, NumTauPnt, +1
                  KSpCrFtTauforFT(NTInd, Nk, I0) = KSpCrFtTauforFT(NTInd, Nk, I0) + &
                                                      & exp(cmplx(0.0_rp, Rtp1, rp)) * RealSpCrFTauBIN(NTInd, imj, I0)
               enddo
                        enddo
                  enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
            
   end subroutine FourierTransTau
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$