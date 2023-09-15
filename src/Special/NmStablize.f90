!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the numerical stablization process for the CPMC simulations.
! COMMENT: Numerical Stablization process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   NmStablizeOne2M --> Subroutine to perform Numerical stablization for \tau=0 --> \tau=BetaT;
!   NmStablizeM2One --> Subroutine to perform Numerical stablization for \tau=BetaT --> \tau=0;
!
!   Rght_NmStablize --> Subroutine used to perform the numerical stablization for URght;
!   Left_NmStablize --> Subroutine used to perform the numerical stablization for ULeft.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine NmStablizeOne2M(NB, NSW, NT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  NmStablizeOne2M(NB, NSW, NT)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the numerical stablization for in the process of sweeping from 
!                     tau=0 to tau=Beta.
! KEYWORDS: Numerical Stablization of URght for static propagation.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: This subroutine is only used in SweepM2One subroutine, as the main part of numerical stablization.
!
!     Perform the numerical stablization for URght for static propagation, as
!                    URght * DRghtVec * VRght
!        including following steps:
!             (0) Obtain the new URght as URght * DRghtVec;
!             (1) UDV decomposition for new URght = URght * DRghtVec * VRght;
!             (2) Obtain new VRght * onld VRght;
!             (3) Compare the Green's function.
!
!     Input:  NB        --> BIN index;
!             NSW       --> Integer index for sweeps;
!             NT        --> Time slice index.
!
!     Output: (none)
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
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB, NSW, NT
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      integer Iwalk
      real(rp) Xmaxm, Xmean
      complex(rp) LogzDet
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for Numerical stablization process ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsNmStb = TimsNmStb + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of Numerical stablization ______________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Initialization for the numerical stabilization process ____________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Set ChemP_Rght and ChemP_Left for computing GrF __________________
!________________________________________________________________________________________  
      if(.not. ReadFldMea) then
         ChemP_Rght = ChemP
         ChemP_Left = ChemP_BT - HubbU_UHF/2.0_rp
      else
         ChemP_Sgle = ChemP
      end if
!**************************************************************************************************	  
!___________________ 1. Perform the numerical stablization for all random walkers _________________
!______________________________ and compare the GrF difference ____________________________________
!**************************************************************************************************
      do Iwalk = 1, NWalk, +1
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@ First check the weight of present walker @@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         if(WghtProc(Iwalk) <= 0.0_rp) then
            cycle
         end if
!________________________________________________________________________________________ 	  
!_________________ (0) Store the GrnFunct matrix for Error comparison ___________________
!________________________________________________________________________________________
         if(IfNmStbErr .and. mod(NSW, StabOutput) == 0) then
            GrnFRTmp00 = 0.0_rp; LogzDet = rp_Zzero
            if(.not. ReadFldMea) then
               call dMat_Copy_QMC(GrnFunct(1, 1, 1, Iwalk), GrnFRTmp00(1, 1, 1))
            else
               call GrFStatcRB(NT, LogScaleRght(1, Iwalk), URght(1, 1, 1, Iwalk), DRghtVec(1, 1, Iwalk), &
                  & VRght(1, 1, 1, Iwalk), LogzDet, GrnFRTmp00(1, 1, 1))
            end if
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) Perform the QR decomposition for (URght*DRghtVec) ________________
!________________________________________________________________________________________ 
         call Rght_NmStablize(Iwalk)
!________________________________________________________________________________________ 	  
!_________________ (2) Calculate GreenF matrix from UDV decomp __________________________
!________________________________________________________________________________________ 
         GrnFunct(:, :, :, Iwalk) = 0.0_rp; LogzDet = rp_Zzero
         if(.not. ReadFldMea) then
            call GrFStaticR_LR(NT, LogScaleRght(1, Iwalk), LogScaleLeft(1), URght(1, 1, 1, Iwalk), &
               & DRghtVec(1, 1, Iwalk), VRght(1, 1, 1, Iwalk), VLeft(1, 1, 1), DLeftVec(1, 1), ULeft(1, 1, 1), &
               & LogzDet, GrnFunct(1, 1, 1, Iwalk))
         else
            call GrFStatcRB(NT, LogScaleRght(1, Iwalk), URght(1, 1, 1, Iwalk), DRghtVec(1, 1, Iwalk), &
               & VRght(1, 1, 1, Iwalk), LogzDet, GrnFunct(1, 1, 1, Iwalk))
         end if
!________________________________________________________________________________________ 	  
!_________________ (3) Compare GrF difference to check the NumStab ______________________
!________________________________________________________________________________________      
         if(IfNmStbErr .and. mod(NSW, StabOutput) == 0) then      
!____________________________________________________________________________ 	  
!________________ [0] Compare GrnFRTmp00 and GrnFunct matrices ______________
!____________________________________________________________________________
            Xmaxm = 0.0_rp; Xmean = 0.0_rp
            call dMat_Compare_QMC(GrnFRTmp00(1, 1, 1), GrnFunct(1, 1, 1, Iwalk), Xmaxm, Xmean)
            NComp_StaDyn(1) = NComp_StaDyn(1) + 1
            if(Xmaxm > Xmaxm_StaDyn(1)) then
               Xmaxm_StaDyn(1) = Xmaxm
            end if
            Xmean_StaDyn(1) = Xmean_StaDyn(1) + Xmean
!____________________________________________________________________________   
!________________ [1] Output the Green's function differences _______________
!____________________________________________________________________________ 
            if(amyid == amstr) then
               !!!!!!!!!! Print the head of numerical Stabilization process
               if( (NSW == StabOutput) .and. (Iwalk == 1) .and. (NT == NvStb) ) then
                  open(291, file = "Output/00_NmStbGrnFStaDyn.txt", access = "append")
                  if(NB == 0) then
                     write(291, "()")
                     write(291, "()")
                     if(IfMuTqmcNt .or. IfFixnT) then
                        write(291, "()")
                        write(291, "('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')")
                        if(IfMuTqmcNt) then
                           write(291, "('@@@@@@@@@@@ Adjust ChemP_BT to reach the QMC density @@@@@@@@@@@@@@@@')")
                        else if(IfFixnT) then
                           write(291, "('@@@@@@@@@@@ Adjust ChemP to reach fixed electron density @@@@@@@@@@@@')")
                        end if
                        write(291, "('@@@@@@@@@@@@@@@@@@@@@@@@@@ Iteration ', I3.3, ' @@@@@@@@@@@@@@@@@@@@@@@@@@@@')") &
                           & Fix_Iterate
                        write(291, "('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')")
                        write(291, "()")
                        write(291, "('________________________________________________________________')")
                        write(291, "('__________________________ CPMCWarmUp __________________________')") 
                        write(291, "('________________________________________________________________')")
                     else
                        write(291, "()")
                        write(291, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
                        write(291, "('&&&&&&&&&&&&&&&& CPMC Simulation for Data Statistics &&&&&&&&&&&&&&&&')")
                        write(291, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
                        write(291, "()")
                        write(291, "('################################################################')")
                        write(291, "('######################### CPMC WarmUp ##########################')") 
                        write(291, "('################################################################')")
                     end if
                  else if(NB > 0) then
                     if(IfMuTqmcNt) then
                        write(291, "()")
                        write(291, "('________________________________________________________________')")
                        write(291, "('__________________________ CPMCOccMuT __________________________')") 
                        write(291, "('________________________________________________________________')")
                     else if(IfFixnT) then
                        write(291, "()")
                        write(291, "('________________________________________________________________')")
                        write(291, "('__________________________ CPMCMeaFix __________________________')") 
                        write(291, "('________________________________________________________________')")
                     else
                        write(291, "()")
                        write(291, "()")
                        write(291, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
                        write(291, "('++++++++++++++++++++++++++ BIN ', I4.4,' ++++++++++++++++++++++++++++')") NB
                        write(291, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
                     end if
                  end if
                  close(291)
               end if
               !!!!!!!!!! Print Diff of GrF matrices between with and without UDV decomps
               if( (Iwalk == merge(NWalk, NWkBt, .not. ReadFldMea)) .and. (NT == LTrot) ) then
                  open( 291, file = "Output/00_NmStbGrnFStaDyn.txt", access = "append")
                  Xmean_StaDyn(1) = Xmean_StaDyn(1) / dble(NComp_StaDyn(1))
                  write(291, "('One2MStatc: NSW, Xmaxm_One2M, Xmean_One2M = ', I4.4, 2es17.8)") NSW, &
                     & Xmaxm_StaDyn(1), Xmean_StaDyn(1)
                  close(291)
               end if
            end if
         end if
      enddo
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for Numerical stablization process ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeNmStb = TimeNmStb + TimeIntrvl(time1, time2)
      
   end subroutine NmStablizeOne2M
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine NmStablizeM2One(ForBckWard, Iwalk, NB, NSW, NTItr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  NmStablizeM2One(ForBckWard, Iwalk, NB, NSW, NTItr)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the numerical stablization for in the process of sweeping from 
!                     tau=Beta to tau=0.
! KEYWORDS: Numerical Stablization of both ULeft and URght for static propagation.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: This subroutine is only used in SweepOne2M subroutine, as the main part of numerical stablization.
!
!     In SweepOne2M subroutine, there is first a sweep from tau=0 to tau=Beta, used to stablize the matrix
!         multiplication and store the UDV matrices; there is also a sweeo from tau=Beta to tau=0, 
!         in which we perform measurements.
!
!     Input:  ForBckWard --> Character to determine the direction of the sweep;
!             Iwalk      --> Integer index for random walker;
!             NB         --> BIN index;
!             NSW        --> Integer index for sweeps;
!             NTItr      --> Time slice index.
!
!     Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		use RealPrecsn
      use TimeRecord
      use QMCTimeRec
      use Observable
		use CoreParamt
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      character(7) ForBckWard
      integer Iwalk, NB, NSW, NTItr
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer(8) time1, time2
      logical IfOutptNmStb
      integer NT_ST
      real(rp) Xmaxm, Xmean
      complex(rp) LogzDet
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for Numerical stablization process ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsB0Stb = TimsB0Stb + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of Numerical stablization ______________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Initialization for the numerical stabilization process ____________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Determine whether to compare the GrF matrices ____________________
!________________________________________________________________________________________
      IfOutptNmStb = .false.
      if( IfNmStbErr .and. mod(NSW, StabOutput) == 0 ) then
         if( (ForBckWard == "bckward" .and. NTItr == 0) .or. (ForBckWard == "forward") ) then
            IfOutptNmStb = .true.
         else
            IfOutptNmStb = .false.
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (1) Set ChemP_Rght and ChemP_Left for computing GrF __________________
!________________________________________________________________________________________
      if(ForBckWard == "bckward") then
         ChemP_Sgle = ChemP
      else if(ForBckWard == "forward") then  
         ChemP_Rght = ChemP
         ChemP_Left = ChemP
      end if
!**************************************************************************************************	  
!___________________ 1. Perform the numerical stablization for Iwalk walker _______________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Determine time slice index for this UDV decomp ___________________
!________________________________________________________________________________________
      NT_ST = NTItr / NvStbM2One
!________________________________________________________________________________________ 	  
!_________________ (1) For the case of ForBckWard == "bckward" __________________________
!_____________________ Perform the UDV decomposition for (DLeftVec*ULeft)^+ _____________
!________________________________________________________________________________________ 
      if(ForBckWard == "bckward") then
!____________________________________________________________________________ 	  
!________________ [0] Compute the GrF matrix for Error comparison ___________
!____________________________________________________________________________
         if(IfOutptNmStb) then
            GrnFRTmp00 = 0.0_rp; LogzDet = rp_Zzero
            call GrFStatcR0(LogScaleLeft(1), VLeft(1, 1, 1), DLeftVec(1, 1), ULeft(1, 1, 1), LogzDet, &
               & GrnFRTmp00(1, 1, 1))
         end if
!____________________________________________________________________________ 	  
!________________ [1] Perform the UDV decomposition _________________________
!____________________________________________________________________________
         call Left_NmStablize()
!____________________________________________________________________________ 	  
!________________ [2] Record the UDV matrices at left side __________________
!____________________________________________________________________________
         if(NTItr /= 0) then
            call dMat_Copy_QMC( ULeft(1, 1, 1),    UStMt(1, 1, 1, NT_ST)   )
            call dcopy(2*NumNS, DLeftVec(1, 1), 1, DVecStMt(1, 1, NT_ST), 1)
            call dMat_Copy_QMC( VLeft(1, 1, 1),    VStMt(1, 1, 1, NT_ST)   )
            DScalLog(1:NmSpn, NT_ST) = LogScaleLeft(1:NmSpn)
         end if
!____________________________________________________________________________ 	  
!________________ [3] Recompute GrFTau00 matrix at NTItr == 0 point _________
!____________________________________________________________________________
         if(NTItr == 0) then
            GrnFunct(:, :, :, Iwalk) = 0.0_rp; LogzDet = rp_Zzero
            call GrFStatcR0(LogScaleLeft(1), VLeft(1, 1, 1), DLeftVec(1, 1), ULeft(1, 1, 1), LogzDet, &
               & GrnFunct(1, 1, 1, Iwalk))
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) For the case of ForBckWard == "forward" __________________________
!_____________________ Perform the UDV decomposition for (URght*DRghtVec) _______________
!________________________________________________________________________________________ 
      if(ForBckWard == "forward") then
!____________________________________________________________________________ 	  
!________________ [0] Store the GrF matrix for Error comparison _____________
!____________________________________________________________________________
         if(IfOutptNmStb) then
            call dMat_Copy_QMC(GrnFunct(1, 1, 1, Iwalk), GrnFRTmp00(1, 1, 1))
            if(IfTAU .and. NTItr <= NmTDM) then
               call dMat_Copy_QMC(GrFT0(1, 1, 1), GrnFRTmp11(1, 1, 1))
               call dMat_Copy_QMC(GrF0T(1, 1, 1), GrnFRTmp22(1, 1, 1))
            end if
         end if
!____________________________________________________________________________ 	  
!________________ [1] Read ULeft, DLeftVec, VLeft from storage ______________
!____________________________________________________________________________
         if(NTItr == LTrot) then
            VLeft = IdMtR; DLeftVec = 1.0_rp; ULeft = IdMtR; LogScaleLeft = 0.0_rp
         else
            call dMat_Copy_QMC( VStMt(1, 1, 1, NT_ST),    VLeft(1, 1, 1)   )
            call dcopy(2*NumNS, DVecStMt(1, 1, NT_ST), 1, DLeftVec(1, 1), 1)
            call dMat_Copy_QMC( UStMt(1, 1, 1, NT_ST),    ULeft(1, 1, 1)   )
            LogScaleLeft(1:NmSpn) = DScalLog(1:NmSpn, NT_ST)
         end if
!____________________________________________________________________________ 	  
!________________ [2] Perform the UDV decomposition _________________________
!____________________________________________________________________________ 
         call Rght_NmStablize(Iwalk)
!____________________________________________________________________________ 	  
!________________ [3] Recompute all the GrF matrices ________________________
!____________________________________________________________________________
         GrnFunct(:, :, :, Iwalk) = 0.0_rp
         if(IfTAU .and. NTItr <= NmTDM) then
            GrFT0 = 0.0_rp; GrF0T = 0.0_rp
            call GrFDynamcR_LR(NTItr, LogScaleRght(1, Iwalk), LogScaleLeft(1), URght(1, 1, 1, Iwalk), &
               & DRghtVec(1, 1, Iwalk), VRght(1, 1, 1, Iwalk), VLeft(1, 1, 1), DLeftVec(1, 1), ULeft(1, 1, 1), LogzDet, &
               & GrnFunct(1, 1, 1, Iwalk), GrFT0(1, 1, 1), GrF0T(1, 1, 1))
         else
            LogzDet = rp_Zzero
            call GrFStaticR_LR(NTItr, LogScaleRght(1, Iwalk), LogScaleLeft(1), URght(1, 1, 1, Iwalk), &
               & DRghtVec(1, 1, Iwalk), VRght(1, 1, 1, Iwalk), VLeft(1, 1, 1), DLeftVec(1, 1), ULeft(1, 1, 1), &
               & LogzDet, GrnFunct(1, 1, 1, Iwalk))
         end if
      end if
!**************************************************************************************************	  
!___________________ 2. Compare the Green's Functions from wrap and UDV Decomps ___________________
!**************************************************************************************************
      if(IfOutptNmStb) then
!________________________________________________________________________________________ 	  
!_________________ (0) For ForBckWard == "bckward" or "forward" cases ___________________
!_____________________ Compare GrnFRTmp00 and GrnFunct matrices _________________________
!________________________________________________________________________________________
         Xmaxm = 0.0_rp; Xmean = 0.0_rp
         call dMat_Compare_QMC(GrnFRTmp00(1, 1, 1), GrnFunct(1, 1, 1, Iwalk), Xmaxm, Xmean)
         NComp_StaDyn(2) = NComp_StaDyn(2) + 1
         if(Xmaxm > Xmaxm_StaDyn(2)) then
            Xmaxm_StaDyn(2) = Xmaxm
         end if
         Xmean_StaDyn(2) = Xmean_StaDyn(2) + Xmean
!________________________________________________________________________________________ 	  
!_________________ (1) For the ForBckWard == "forward" case _____________________________
!_____________________ Compare GrnFRTmp11 and GrFT0 matrices ____________________________
!_____________________ Compare GrnFRTmp22 and GrF0T matrices ____________________________
!________________________________________________________________________________________
         if( (ForBckWard == "forward") .and. (IfTAU) .and. (NTItr <= NmTDM) ) then
            !!!!!!!!!! Compare GrnFRTmp11 and GrFT0 matrices
            Xmaxm = 0.0_rp; Xmean = 0.0_rp
            call dMat_Compare_QMC(GrnFRTmp11(1, 1, 1), GrFT0(1, 1, 1), Xmaxm, Xmean)
            NComp_StaDyn(3) = NComp_StaDyn(3) + 1
            if(Xmaxm > Xmaxm_StaDyn(3)) then
               Xmaxm_StaDyn(3) = Xmaxm
            end if
            Xmean_StaDyn(3) = Xmean_StaDyn(3) + Xmean
            !!!!!!!!!! Compare GrnFRTmp22 and GrF0T matrices
            Xmaxm = 0.0_rp; Xmean = 0.0_rp
            call dMat_Compare_QMC(GrnFRTmp22(1, 1, 1), GrF0T(1, 1, 1), Xmaxm, Xmean)
            NComp_StaDyn(3) = NComp_StaDyn(3) + 1
            if(Xmaxm > Xmaxm_StaDyn(3)) then
               Xmaxm_StaDyn(3) = Xmaxm
            end if
            Xmean_StaDyn(3) = Xmean_StaDyn(3) + Xmean
         end if
!________________________________________________________________________________________ 	  
!_________________ (2) Output GrFs difference for present sweep process _________________
!________________________________________________________________________________________
         if( (amyid == amstr) .and. (ForBckWard == "forward") .and. (Iwalk == IdptWkIndx(NWkBt)) .and. &
               & (NTItr == LTrot) ) then
            Xmean_StaDyn(2) = Xmean_StaDyn(2) / dble(NComp_StaDyn(2))
            open( 291, file = "Output/00_NmStbGrnFStaDyn.txt", access = "append")
            write(291, "('M2OneStatc: NSW, Xmaxm_M2One, Xmean_M2One = ', I4.4, 2es17.8)") NSW, &
               & Xmaxm_StaDyn(2), Xmean_StaDyn(2)
            if(IfTAU) then
               Xmean_StaDyn(3) = Xmean_StaDyn(3) / dble(NComp_StaDyn(3))
               write(291, "('M2OneDynmc: NSW, Xmaxm_Dynmc, Xmean_Dynmc = ', I4.4, 2es17.8)") NSW, &
                  & Xmaxm_StaDyn(3), Xmean_StaDyn(3)
            end if
            close(291)
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for Numerical stablization process ________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		call system_clock(time2)
      TimeB0Stb = TimeB0Stb + TimeIntrvl(time1, time2)
      
   end subroutine NmStablizeM2One
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine Rght_NmStablize(Iwalk)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  Rght_NmStablize(Iwalk)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the numerical stablization for URght, DRghtVec and VRght, for 
!                     static propagation.
! KEYWORDS: Numerical Stablization of URght for static propagation.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Just as  URght * DRghtVec * VRght = URght' * DRghtVec' * VRght'.
!
!        including following steps:
!             (0) Calculate URght * DRghtVec;
!             (1) UDV decomposition for URght * DRghtVec = UR * DRVec * VR and
!                                       Final URght    = UR
!                                             DRghtVec = DRVec
!                                             VRght    = VR * VRght
!
!     Input:  Iwalk --> Integer index for the random walker;
!
!     Output: (none)
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
      integer Iwalk
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, SpnInd
      real(rp) Rtp0, Rtp1, Rtp2
      real(rp), allocatable :: TmpMatrix1(:, :, :)   ! Temporary real matrix
      real(rp), allocatable :: TmpMatrix2(:, :, :)   ! Temporary real matrix
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of Numerical stablization ______________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Allocate the matrices used in this subroutine _____________________________
!**************************************************************************************************
      allocate(TmpMatrix1(NumNS, NumNS, NmSpn)); TmpMatrix1 = 0.0_rp
		allocate(TmpMatrix2(NumNS, NumNS, NmSpn)); TmpMatrix2 = 0.0_rp
!**************************************************************************************************	  
!___________________ 1. Perform the numerical stablization for (URght*DRghtVec) ___________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Calculate TmpMatrix2 = URght * DRghtVec __________________________
!_____________________ Store     TmpMatrix1 = VRght _____________________________________
!________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               TmpMatrix2(I1, I2, SpnInd) = URght(I1, I2, SpnInd, Iwalk) * DRghtVec(I2, SpnInd, Iwalk)
               TmpMatrix1(I1, I2, SpnInd) = VRght(I1, I2, SpnInd, Iwalk)
            enddo
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!________________________________________________________________________________________ 	  
!_________________ (1) UDV decomposition for (URght*DRghtVec) matrix ____________________
!________________________________________________________________________________________
      URght(:, :, :, Iwalk) = 0.0_rp
      DRghtVec(:, :, Iwalk) = 0.0_rp
      VRght(:, :, :, Iwalk) = 0.0_rp
      call dMatUDV_QMC(NumNS, NumNS, TmpMatrix2(1, 1, 1), URght(1, 1, 1, Iwalk), DRghtVec(1, 1, Iwalk), &
         & VRght(1, 1, 1, Iwalk))
!________________________________________________________________________________________ 	  
!_________________ (2) VRght --> VRght * VRght_{Old} ____________________________________
!________________________________________________________________________________________
      TmpMatrix2 = 0.0_rp
      call dMatPrd_NN_QMC(NumNS, NumNS, NumNS, 1.0_rp, VRght(1, 1, 1, Iwalk), TmpMatrix1(1, 1, 1), &
         & 0.0_rp, TmpMatrix2(1, 1, 1))
      call dMat_Copy_QMC(TmpMatrix2(1, 1, 1), VRght(1, 1, 1, Iwalk))
!________________________________________________________________________________________ 	  
!_________________ (3) Rescale the diagonal values in DRghtVec matrix ___________________
!_____________________ Target: DRghtVec(1, 1)*DRghtVec(NumNS, 1) = 1 ____________________
!________________________________________________________________________________________
      do SpnInd = 1, NmSpn, +1
         Rtp0 = ( log(DRghtVec(1, SpnInd, Iwalk)) + log(DRghtVec(NumNS, SpnInd, Iwalk)) ) / 2.0_rp
         do I1 = 1, NumNS, +1
            DRghtVec(I1, SpnInd, Iwalk) = DRghtVec(I1, SpnInd, Iwalk) * exp(-Rtp0)
         enddo
         LogScaleRght(SpnInd, Iwalk) = LogScaleRght(SpnInd, Iwalk) + Rtp0
      enddo
!**************************************************************************************************	  
!___________________ 2. Deallocate all the used matrices here _____________________________________
!**************************************************************************************************      
      if(allocated(TmpMatrix1)) deallocate(TmpMatrix1)
		if(allocated(TmpMatrix2)) deallocate(TmpMatrix2)
      
   end subroutine Rght_NmStablize
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine Left_NmStablize()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  Left_NmStablize()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the numerical stablization for ULeft, DLeftVec and VLeft, for 
!                    static propagation.
! KEYWORDS: Numerical Stablization of ULeft for static propagation.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Just as  VLeft * DLeftVec * ULeft = VLeft' * DLeftVec' * ULeft'.
!
!     Perform the numerical stablization for ULeft for static propagation, as
!                    ULeft^+ = UL * DLVec * VL
!                    --> ULeft = VL^+ * DLVec * UL^+
!        including following steps:
!             (0) Calculate (DLeftVec*ULeft)^+;
!             (1) UDV decomposition for (DLeftVec*ULeft)^+ = UL * DLVec * VL and
!                                       DLeftVec*ULeft = VL^+ * DLVec * UL^+
!                                       Final ULeft    = UL^+;
!                                             DLeftVec = DLVec
!                                             VLeft = VLeft * VL^+
!
!     Input:  (none). Output: (none)
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
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, SpnInd
      real(rp) Rtp0, Rtp1, Rtp2
      real(rp), allocatable :: TmpMatrix1(:, :, :)   ! Temporary real matrix
      real(rp), allocatable :: TmpMatrix2(:, :, :)   ! Temporary real matrix
      real(rp), allocatable :: TmpMatrix3(:, :, :)   ! Temporary real matrix
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of Numerical stablization ______________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Allocate the matrices used in this subroutine _____________________________
!**************************************************************************************************
      allocate(TmpMatrix1(NumNS, NumNS, NmSpn)); TmpMatrix1 = 0.0_rp
      allocate(TmpMatrix2(NumNS, NumNS, NmSpn)); TmpMatrix2 = 0.0_rp
      allocate(TmpMatrix3(NumNS, NumNS, NmSpn)); TmpMatrix3 = 0.0_rp
!**************************************************************************************************	  
!___________________ 1. Perform the numerical stablization for (DLeftVec*ULeft)^+ _________________
!************************************************************************************************** 
!________________________________________________________________________________________ 	  
!_________________ (0) Calculate DLeftVec * ULeft _______________________________________
!_____________________ Prepare TmpMatrix1 = (DLeftVec*ULeft)^+ __________________________
!________________________________________________________________________________________
      TmpMatrix1 = 0.0_rp
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               TmpMatrix1(I1, I2, SpnInd) = ULeft(I2, I1, SpnInd) * DLeftVec(I2, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!________________________________________________________________________________________ 	  
!_________________ (1) UDV for TmpMatrix1 = U1 * D1 * V1 ________________________________
!_____________________ DLeftVec*ULeft = V1^+ * D1 * U1^+ ________________________________
!________________________________________________________________________________________
      TmpMatrix2 = 0.0_rp; DLeftVec = 0.0_rp; TmpMatrix3 = 0.0_rp
      call dMatUDV_QMC(NumNS, NumNS, TmpMatrix1(1, 1, 1), TmpMatrix2(1, 1, 1), DLeftVec(1, 1), TmpMatrix3(1, 1, 1))
!________________________________________________________________________________________ 	  
!_________________ (2) Obtain the new ULeft = U1^+ ______________________________________
!________________________________________________________________________________________      
      ULeft = 0.0_rp
   !$OMP PARALLEL &
   !$OMP PRIVATE(I2, SpnInd, I1)
   !$OMP DO
      do I2 = 1, NumNS, +1
         do SpnInd = 1, NmSpn, +1
            do I1 = 1, NumNS, +1
               ULeft(I1, I2, SpnInd) = TmpMatrix2(I2, I1, SpnInd)
            enddo
         enddo
      enddo
   !$OMP END DO     
   !$OMP END PARALLEL     
!________________________________________________________________________________________ 	  
!_________________ (3) Obtain the new VLeft = VLeft * V1^+ ______________________________
!________________________________________________________________________________________
      call dMat_Copy_QMC(VLeft(1, 1, 1), TmpMatrix1(1, 1, 1))
      VLeft = 0.0_rp
      call dMatPrd_NT_QMC(NumNS, NumNS, NumNS, 1.0_rp, TmpMatrix1(1, 1, 1), TmpMatrix3(1, 1, 1), &
         & 0.0_rp, VLeft(1, 1, 1))
!________________________________________________________________________________________ 	  
!_________________ (4) Rescale the diagonal values in DLeftVec matrix ___________________
!_____________________ Target: DLeftVec(1, 1)*DLeftVec(NumNS, 1) = 1 ____________________
!________________________________________________________________________________________
      do SpnInd = 1, NmSpn, +1
         Rtp0 = ( log(DLeftVec(1, SpnInd)) + log(DLeftVec(NumNS, SpnInd)) ) / 2.0_rp
         do I1 = 1, NumNS, +1
            DLeftVec(I1, SpnInd) = DLeftVec(I1, SpnInd) * exp(-Rtp0)
         enddo
         LogScaleLeft(SpnInd) = LogScaleLeft(SpnInd) + Rtp0
      enddo
!**************************************************************************************************	  
!___________________ 2. Deallocate all the used matrices here _____________________________________
!**************************************************************************************************      
      if(allocated(TmpMatrix1)) deallocate(TmpMatrix1)
		if(allocated(TmpMatrix2)) deallocate(TmpMatrix2)
      if(allocated(TmpMatrix3)) deallocate(TmpMatrix3)
      
   end subroutine Left_NmStablize
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$