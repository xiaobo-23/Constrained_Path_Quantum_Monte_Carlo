!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform Initialization of key matrices used in CPMC before the whole 
!                 CPMC simulations.
! COMMENT: CPMC Initialization process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   InitMCVecMat --> Subroutine to perform Initialization of key matrices used in CPMC in CPMC Simulations.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine InitCoreVcMt()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitMCVecMat()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the initialization for Core CPMC vectors and matrices.
! KEYWORDS: Initialization of key vectors and matrices in CPMC.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-05-24
! DESCRIPTION:
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
      use Observable
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, SpnInd, NT
      real(rp) Rtp0
!______________________________________________________________________________________________________________	  
!_______________________________ Initialization of Core CPMC vectors and matrices _____________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of initialization process _______________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'InitCoreVcMt: Initialization of core vectors and matrices in CPMC!')")
      end if
!**************************************************************************************************	  
!___________________ 3. The Core matrices used in CPMC simulations ________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) UDV matrices at left and right sides for propagation _____________
!________________________________________________________________________________________
      !!!!!!!!!! The UDV matrices for left side
      allocate(ULeft   (NumNS, NumNS, NmSpn)); ULeft    = 0.0_rp
      allocate(DLeftVec(NumNS       , NmSpn)); DLeftVec = 0.0_rp
      allocate(VLeft   (NumNS, NumNS, NmSpn)); VLeft    = 0.0_rp
      !!!!!!!!!! The UDV matrices for right side
      allocate(URght   (NumNS, NumNS, NmSpn, NWalk)); URght    = 0.0_rp
      allocate(DRghtVec(NumNS       , NmSpn, NWalk)); DRghtVec = 0.0_rp
      allocate(VRght   (NumNS, NumNS, NmSpn, NWalk)); VRght    = 0.0_rp
      !!!!!!!!!! The LogScale for left and right
      allocate(LogScaleRght(NmSpn, NWalk)); LogScaleRght = 0.0_rp
      allocate(LogScaleLeft(NmSpn       )); LogScaleLeft = 0.0_rp
      !!!!!!!!!! The LogScale for e^(-\tau*H_T)
      allocate(LogScaleOfHT(NmSpn, 0:LTrot)); LogScaleOfHT = 0.0_rp
      do SpnInd = 1, NmSpn, +1
         Rtp0 = ( HT_EigValu(1, SpnInd) + HT_EigValu(NumNS, SpnInd) ) / 2.0_rp
         do NT = 0, LTrot, +1
            LogScaleOfHT(SpnInd, NT) = - (LTrot-NT)*Dltau * Rtp0
         enddo
      enddo
!________________________________________________________________________________________ 	  
!_________________ (1) The Green's function matrix for propagation ______________________
!________________________________________________________________________________________
      GrFTrcDhld = 32.0_rp
      allocate(GrnFunct(NumNS, NumNS, NmSpn, 0:NWalk)); GrnFunct = 0.0_rp
!________________________________________________________________________________________ 	  
!_________________ (2) Weight and their logs, and all weights vector ____________________
!________________________________________________________________________________________ 
      !!!!!!!!!! Indexes for indepdent random walkers of present process
      allocate(IdptWkIndx(NWalk)); IdptWkIndx = -10000
      !!!!!!!!!! The Weight and their logs of walkers on every process
      allocate(WghtProc(NWalk)); WghtProc = 1.0_rp
      allocate(Log_Wght(NWalk)); Log_Wght = 0.0_rp
      !!!!!!!!!! The weight of all random walkers on all processes, used in population control
#ifdef MPIPROCESS
      if(amyid == amstr) then
         allocate(WghtTotl(NmWalkAllP)); WghtTotl = 0.0_rp
      else
         allocate(WghtTotl(1)); WghtTotl = 0.0_rp
      end if
#endif
!________________________________________________________________________________________ 	  
!_________________ (3) Calculate the growth estimator ___________________________________
!________________________________________________________________________________________
      allocate(GrowthCoefft(1:LTrot)); GrowthCoefft(1:LTrot) = log(1.0_rp)
!________________________________________________________________________________________ 	  
!_________________ (4) The ancestry link and number of ancestry walkers _________________
!________________________________________________________________________________________
      NancestryWalker = NmWalkAllP
      allocate(AncestryLink(NWalk)); AncestryLink = -100
!________________________________________________________________________________________ 	  
!_________________ (5) Store the UDV matrices for measurement in [BetaT, 0] sweep _______
!________________________________________________________________________________________ 
      if(IfM2OneMea) then
         DmUST = LTrot / NvStbM2One
		   if(mod(LTrot, NvStbM2One) /=0) DmUST = DmUST + 1
         allocate(UStMt   (NumNS, NumNS, NmSpn, DmUST-1)); UStMt    = 0.0_rp
         allocate(DVecStMt(NumNS,        NmSpn, DmUST-1)); DVecStMt = 0.0_rp
         allocate(DScalLog(              NmSpn, DmUST-1)); DScalLog = 0.0_rp
         allocate(VStMt   (NumNS, NumNS, NmSpn, DmUST-1)); VStMt    = 0.0_rp
      end if
!________________________________________________________________________________________  
!_________________ (6) Temporary Green's Function matrices for Static and Dynamic _______
!________________________________________________________________________________________  
      allocate(GrnFRTmp00(NumNS, NumNS, 2)); GrnFRTmp00 = 0.0_rp
      allocate(GrnFRTmp11(NumNS, NumNS, 2)); GrnFRTmp11 = 0.0_rp
      allocate(GrnFRTmp22(NumNS, NumNS, 2)); GrnFRTmp22 = 0.0_rp
      if(IfFftEnPar) then
         allocate(GrFKCTmp00(NumNS, NumNS, 2)); GrFKCTmp00 = rp_Zzero
         allocate(GrFKCTmp11(NumNS, NumNS, 2)); GrFKCTmp11 = rp_Zzero
      end if
!________________________________________________________________________________________  
!_________________ (7) Temporary Green's Function matrices for Dynamic __________________
!________________________________________________________________________________________ 
      if(IfTAU) then
         allocate(GrF00(NumNS, NumNS, 2)); GrF00 = 0.0_rp       ! <c_i(\tau=0) c_j^+>
         allocate(GrFT0(NumNS, NumNS, 2)); GrFT0 = 0.0_rp       ! <c_i(tau) c_j^+(0)>
		   allocate(GrF0T(NumNS, NumNS, 2)); GrF0T = 0.0_rp       ! <c_j^+(tau) c_i(0)>
         allocate(GrFTT(NumNS, NumNS, 2)); GrFTT = 0.0_rp       ! <c_i(tau) c_j^+(tau)>
      end if
!________________________________________________________________________________________ 	  
!_________________ (8) Identity matrix, used for initializations ________________________
!________________________________________________________________________________________      
      allocate(IdMtR(NumNS, NumNS, NmSpn)); IdMtR = 0.0_rp
      do SpnInd = 1, NmSpn, +1
         do I0 = 1, NumNS, +1
            IdMtR(I0, I0, SpnInd) = 1.0_rp
         enddo
      enddo
      
   end subroutine InitCoreVcMt
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$