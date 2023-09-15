!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the propagation and updating process from NT=1 to NT=LTrot for URght.
! COMMENT: Sweep and measure for finite-T CPMC simulations.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   SweepOne2M   --> Subroutine to perform the propagation and updating process from NT=1 to NT=LTrot;
!
!   PropgtUpdtBx --> Subroutine to carry out the propagation and updating for URght in B_X;
!   WeightMaxMin --> Subroutine to output the maximum and minimum of weights during propagation.
!           
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine SweepOne2M(IfMeasure, NB, NSW)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SweepOne2M(IfMeasure, NB, NSW) 
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the sweep from \tau=1 to \tau=M time slice, as 0 --> BetaT.
! KEYWORDS: 1-->M (0-->BetaT) sweep.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Perform the sweep from \tau=1 to \tau=M time slice in every Full sweep.
!
!     Input: IfMeasure --> Whether to perform the measurements;
!            NB        --> Integer index for BIN;
!            NSW       --> Integer index for the sweep.
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
      use MPISetting
            implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      logical IfMeasure     ! Whether to measure the physical quantities
      integer NB, NSW
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer NT, I0, Iwalk
!______________________________________________________________________________________________________________     
!_________________________ Main calculations of Sweep from t=M slice to t=0 slice _____________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************     
!___________________ 0. Perform some initializations before the sweep _____________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Initialize walker weights and ancestry links for NT == 0 _________
!________________________________________________________________________________________
      !!!!!!!!!! Weights of random walkers
      WghtProc(1:NWalk) = 1.0_rp; Log_Wght(1:NWalk) = 0.0_rp
      !!!!!!!!!! Number of ancestry walkers and the links
      NancestryWalker = NmWalkAllP
      do Iwalk = 1, NWalk, +1
         AncestryLink(Iwalk) = amyid*NWalk + Iwalk
      enddo  
!________________________________________________________________________________________         
!_________________ (1) Initialize the propagation matrices for NT == 0 __________________
!________________________________________________________________________________________
      !!!!!!!!!! The GrnFunct(:, :, :, Iwalk) matrice
      do Iwalk = 1, NWalk, +1
         call dMat_Copy_QMC(GrnFunct(1, 1, 1, 0), GrnFunct(1, 1, 1, Iwalk))
      enddo
      !!!!!!!!!! The UDV matrices at left and right sides
      call InitUDVMatrx()
!________________________________________________________________________________________         
!_________________ (2) For adjusting GrowthCoefft(:) constants __________________________
!________________________________________________________________________________________
      Sum_WghtOld = dble(NmWalkAllP); Sum_WghtNew = 0.0_rp
      UpdtRtTotMax = -1.0E+64_rp; UpdtRtTotMin = +1.0E+64_rp
!________________________________________________________________________________________         
!_________________ (3) For comparisons of GrFs in NumStab process _______________________
!________________________________________________________________________________________
      NComp_StaDyn = 0; Xmaxm_StaDyn = 0.0_rp; Xmean_StaDyn = 0.0_rp
!________________________________________________________________________________________         
!_________________ (4) Initialize the physical observables ______________________________
!________________________________________________________________________________________       
      EngOccCrFSwp = 0.0_rp
      RlSpGrnFtSwp = 0.0_rp
      if(IfFftEnPar) KSpGreenFSwp = rp_Zzero
      if(IfFftEnPar) then
         NkSgleSwp = rp_Zzero; PairMtSwp = rp_Zzero
      end if
      if(abs(PinSz) < rp_Eps) RealSpCrFSwp = 0.0_rp
!________________________________________________________________________________________         
!_________________ (5) For ReadFldMea == T case, read the path stored before ____________
!________________________________________________________________________________________
      if(ReadFldMea .and. NB >= 1) call SaveReadPathWt("Read")
!**************************************************************************************************     
!___________________ 1. Perform the sweep from  tau == 0  to  tau == BetaT ________________________
!**************************************************************************************************
      do NT = 1, LTrot, +1
!________________________________________________________________________________________         
!_________________ (0) Propagate B_T and B_x operators for NT time slice ________________
!________________________________________________________________________________________
         call PropgtUpdtBx(NB, NSW, NT)
!________________________________________________________________________________________         
!_________________ (1) Perform the Numerical Stablization process _______________________
!________________________________________________________________________________________
         if(mod(NT, NvStb) == 0 .or. NT == LTrot) then
            call NmStablizeOne2M(NB, NSW, NT)
         end if
!________________________________________________________________________________________         
!_________________ (2) Perform the population control process ___________________________
!_____________________ And Adjust the GrowthCoefft constants ____________________________
!________________________________________________________________________________________
         if(.not. ReadFldMea) then
            if( (mod(NT, NvPop) == 0 .or. NT == LTrot) ) then
               call PopControl(NB, NSW, NT)
               if( (NB == 0) .or. (NB > 0 .and. mod(NSW, FrqReCmptGrowth) == 0) ) then
                  call GrwthCoeff(NB, NSW, NT)
               end if
            end if
         end if
      enddo
!**************************************************************************************************     
!___________________ 2. Perform the static measurement at Tau == BetaT point ______________________
!______________________ and stored all the sampled paths for SaveFldMea == T case _________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) Static Measurement at Tau == BetaT point _________________________
!________________________________________________________________________________________      
      if(IfMeasure) then
         NObsStat = NObsStat + 1
         call PhyMeaStatBetaT(NB, NSW)
      end if
!________________________________________________________________________________________         
!_________________ (1) For SaveFldMea == T case, read the path stored before ____________
!________________________________________________________________________________________
      if((SaveFldMea) .and. (.not. IfMuTqmcNt) .and. (.not. IfFixnT) .and. (NB >= 1)) then
         call SaveReadPathWt("Save")
      end if
      
   end subroutine SweepOne2M
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


      
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine PropgtUpdtBx(NB, NSW, NT) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PropgtUpdtBx(NB, NSW, NT)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to propagate and update the B_X matrices at Tau == NT time slices.
! KEYWORDS: Propagate and Update at Tau == NT time slice.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Propagate and update the B_X matrices at Tau == NT time slice.
!
!     Input: NB  --> Integer index for BIN;
!            NSW --> Integer index for the sweep.
!            NT  --> The time slice;
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
      use RandomNumb
            implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NB, NSW, NT
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, Iwalk, SpnInd
!______________________________________________________________________________________________________________     
!_______________________ Main calculations of Propagate and update at Tau == NT slice _________________________
!_______________________________________ For all the random walkers ___________________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************
!___________________ 0. Propagate B_T operator --> Reset VLeft, DLeftVec, ULeft, LogScaleLeft _____
!************************************************************************************************** 
!________________________________________________________________________________________         
!_________________ (0) For NT == LTrot, Set VLeft, DLeftVec, ULeft, LogScaleLeft ________
!________________________________________________________________________________________
      if(NT == LTrot) then
         VLeft = IdMtR; DLeftVec = 1.0_rp; ULeft = IdMtR; LogScaleLeft = 0.0_rp
!________________________________________________________________________________________         
!_________________ (1) For NT /= LTrot, Set DLeftVec and LogScaleLeft ___________________
!________________________________________________________________________________________
      else
         do SpnInd = 1, NmSpn, +1
            do I0 = 1, NumNS, +1
               DLeftVec(I0, SpnInd) = exp( - (LTrot-NT)*Dltau*HT_EigValu(I0, SpnInd) - LogScaleOfHT(SpnInd, NT) )
            enddo
            LogScaleLeft(SpnInd) = LogScaleOfHT(SpnInd, NT)
         enddo
      end if
!**************************************************************************************************
!___________________ 1. Propagate and update B_x for all random walkers ___________________________
!**************************************************************************************************      
      do Iwalk = 1, NWalk, +1
!________________________________________________________________________________________         
!_________________ (0) Choose NT auxiliary fields and update the weights ________________
!________________________________________________________________________________________
         if(.not. ReadFldMea) then
!____________________________________________________________________________         
!________________ [0] Propagate GrF = exp(-dt*H_T) * GrF * exp(+dt*H_T) _____
!____________________________________________________________________________
            call RghtMultExpHT   (.false., GrnFunct(1, 1, 1, Iwalk), NumNS)
            call LeftMultExpHTInv(.false., NumNS, GrnFunct(1, 1, 1, Iwalk)) 
!____________________________________________________________________________         
!________________ [1] Choose or Update the auxiliary fields _________________
!____________________________________________________________________________
            if(abs(HubbU) > rp_Eps) then
               !!!!!!!!!! Upadte the auxiliary fields
               call UpdtIsngbU(Iwalk, NT)
               !!!!!!!!!! Check the weight of present random walker
               if(WghtProc(Iwalk) <= 0.0_rp) cycle
            end if
!____________________________________________________________________________         
!________________ [2] Rescale the weights with growth estimator _____________ 
!____________________________________________________________________________
            Log_Wght(Iwalk) = Log_Wght(Iwalk) +     GrowthCoefft(NT)
            WghtProc(Iwalk) = WghtProc(Iwalk) * exp(GrowthCoefft(NT))
         end if
!________________________________________________________________________________________         
!_________________ (1) Progapate URght for the numerical stablization ___________________
!_____________________ URght = Exp(-dt*H_T) * URght _____________________________________
!_____________________ URght = Exp[-dt*(H_U+H_0-H_T)] * URght ___________________________
!________________________________________________________________________________________
         call RghtMultExpHT(.false., URght(1, 1, 1, Iwalk), NumNS)
         if( abs(HubbU) > rp_Eps ) then
            call RghtMultExpbU_H0T(IsingbU(1, NT, Iwalk), URght(1, 1, 1, Iwalk), NumNS)
         end if
      enddo
!##################################################################################################
!___________________ 2. Output maximum and minimum of Weights of all walkers ______________________
!##################################################################################################
      if(.not. ReadFldMea) then
         if( (.not. IfFixnT) .and. (NB >= 1) .and. (NSW == WghtOutNSW) ) then
            call WeightMaxMin(NT)
         end if
      end if
            
   end subroutine PropgtUpdtBx
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine WeightMaxMin(NT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  WeightMaxMin(NT)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to obtain and output the maximum and minimum of weights of all random walkers.
! KEYWORDS: Maximum and minimum of weights.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Obtain and output the maximum and minimum of weights of all random walkers
!
!     Input:  NT --> Integer index for time slic.
!
!     Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use MPISetting
            use CoreParamt
      implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      integer NT
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________      
      integer Iwalk, I1, I2, SpnInd
      real(rp) ProcWghtInfo(0:4)
      real(rp) FinlWghtInfo(0:4)
      real(rp), allocatable :: GlobWghtInfo(:, :)
!______________________________________________________________________________________________________________     
!_____________________ Main calculations of getting and outputting Max and Min of weights _____________________
!______________________________________________________________________________________________________________
!**************************************************************************************************     
!___________________ 0. Obtain the Max and Min weights of all walkers _____________________________
!**************************************************************************************************
!________________________________________________________________________________________         
!_________________ (0) For randon walkers on every single process _______________________
!________________________________________________________________________________________
      ProcWghtInfo(0) = sum(WghtProc(:))/dble(NWalk)
      ProcWghtInfo(1) = maxval(WghtProc(:), 1)
      ProcWghtInfo(2) = minval(WghtProc(:), 1)
      ProcWghtInfo(3) = UpdtRtTotMax
      ProcWghtInfo(4) = UpdtRtTotMin
!________________________________________________________________________________________         
!_________________ (1) Gather results of all processes and Max, Min _____________________
!________________________________________________________________________________________
#ifdef MPIPROCESS
      allocate(GlobWghtInfo(0:4, 0:anprc-1))
      GlobWghtInfo = 0.0_rp
      call MPI_Barrier(acomm, ierr)
      call MPI_GATHER(ProcWghtInfo(0), 5, rp_MPI_REAL, GlobWghtInfo(0, 0), 5, rp_MPI_REAL, amstr, acomm, ierr)
      call MPI_Barrier(acomm, ierr)
      if(amyid == amstr) then
         FinlWghtInfo(0) = sum(GlobWghtInfo(0, 0:anprc-1))/dble(anprc)
         FinlWghtInfo(1) = maxval(GlobWghtInfo(1, 0:anprc-1), 1)
         FinlWghtInfo(2) = minval(GlobWghtInfo(2, 0:anprc-1), 1)
         FinlWghtInfo(3) = maxval(GlobWghtInfo(3, 0:anprc-1), 1)
         FinlWghtInfo(4) = minval(GlobWghtInfo(4, 0:anprc-1), 1)
      end if
      if(allocated(GlobWghtInfo)) deallocate(GlobWghtInfo)
#else
      FinlWghtInfo = ProcWghtInfo
#endif
!**************************************************************************************************     
!___________________ 1. Output the results of Max and Min weights _________________________________
!**************************************************************************************************
      if(amyid == amstr) then
         open( 293, file = "Output/88_AllWeightMaxMin.txt", access = "append")
         write(293, "(f7.3, A, es17.8, A, es17.8, A, es17.8, A, es17.8, A, es17.8)") dble(NT)*Dltau, char(9), &
            & FinlWghtInfo(0), char(9), FinlWghtInfo(1), char(9), FinlWghtInfo(2), char(9), FinlWghtInfo(3), &
            & char(9), FinlWghtInfo(4)
         close(293)
      end if
      
   end subroutine WeightMaxMin
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$