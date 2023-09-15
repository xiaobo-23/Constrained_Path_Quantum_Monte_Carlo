!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform the population control for all the random walkers during the propagation.
! COMMENT: Population control of random walkers.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-06-15
! PURPOSE: Different subroutines are introduced as following:
!
!   Perform the population control for all the random walkers. 
!             
!   PopControl --> Subroutine used to perform the population control for the CPMC simulation.
!
!   IfPopControl --> Subroutine to determine whether the population control is really needed;
!
!   PopControlNT --> Subroutine to carry out the population control for both parallel and serial cases for NT < LTrot;
!   ReDistrbtWalkNT --> Subroutine to re-distribute all the random walkers according to the weights;
!   PopCtrlPackData --> Subroutine used to   pack the data for communication in population control;
!   PopCtrlUnPkData --> Subroutine used to Unpack the data for communication in population control;
!   CopyWk_Iw_to_Jw --> Subroutine used to copy data from Jw-th walker to Iw-th walker;
!
!   PopCtrlBetaT --> Subroutine to carry out the population control for both parallel and serial cases for NT = LTrot;
!   ReDistrbWkBetaT --> Subroutine to re-distribute all the random walkers according to the weights for NT = LTrot;
!
!   GetNAncestry --> Subroutine to obtain the number ancestry random walkers after population control;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine PopControl(NB, NSW, NT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PopControl(NB, NSW, NT)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the population control for all the random walkers according to 
!                    the weight distribution.
! KEYWORDS: Population control.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Perform the population control for all the random walkers according to the weight distribution. 
!                 This subroutine has some nontrivial application of MPI computations.
!
!     Input:  NB  --> Integer index for the BIN simulation;
!             NSW --> Integer index of the present sweep for both warm-up and measure process.
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
      logical(4) IfPopCtrl
      integer Iwalk, NumZeroWght
      integer NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv
      real(rp) MeanWghtNT
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for population control process ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
		TimsPopCt = TimsPopCt + 1
		call system_clock(time1)
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of population control process __________________________________
!______________________________________________________________________________________________________________      
!**************************************************************************************************	  
!___________________ 0. Perform the population control process ____________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Process the weights to determine whether pop or not ______________
!________________________________________________________________________________________
      NumZeroWght = 0
      call IfPopControl(IfPopCtrl, NumZeroWght, MeanWghtNT)
!________________________________________________________________________________________ 	  
!_________________ (1) If IfPopCtrl == .true., perform population control _______________
!________________________________________________________________________________________
      if(IfPopCtrl) then
         NIndpt = 0; NumChg = 0; MaxNum = 0; NumCom = 0; NumSnd = 0; NumRcv = 0
         if(NT == LTrot) then
            call PopCtrlBetaT(NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv)
         else
            call PopControlNT(NT, NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv)
         end if
         call GetNAncestry()
      else
         NWkBt = NWalk
         do Iwalk = 1, NWkBt, +1
            IdptWkIndx(Iwalk) = Iwalk
         enddo
         NIndpt = NmWalkAllP
      end if
!**************************************************************************************************	  
!___________________ 1. Output information about population control process _______________________
!**************************************************************************************************
      if(amyid == amstr) then
         !!!!!!!!!! Print head of output for the present BIN
         if( (NSW == 1) .and. (NT == NvPop) ) then
            open( 291, file = "Output/00_PopulationCntrl.txt", access = "append")
            write(291, "()")
            write(291, "()")   
            if(NB == 0) then
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
                  write(291, "('________________________________________________________________')")
                  write(291, "('__________________________ CPMCOccMuT __________________________')") 
                  write(291, "('________________________________________________________________')")
               else if(IfFixnT) then
                  write(291, "('________________________________________________________________')")
                  write(291, "('__________________________ CPMCMeaFix __________________________')") 
                  write(291, "('________________________________________________________________')")
               else
                  write(291, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
                  write(291, "('++++++++++++++++++++++++++ BIN ', I4.4,' ++++++++++++++++++++++++++++')") NB
                  write(291, "('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')")
               end if
            end if
            close(291)
         end if
         !!!!!!!!!! Print Information of population control for the PoptOutput sweep
         if( NSW == PoptOutput ) then
            open(291, file = "Output/00_PopulationCntrl.txt", access = "append")
            if(NT == NvPop) then
               write(291, "()")
               write(291, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
               write(291, "('&&&&&&&&&&&&&&&&&&&&&&&& Sweep ', I4.4,' &&&&&&&&&&&&&&&&&&&&&&&')") NSW
               write(291, "('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')")
               write(291, "(1x, 'NT', 8x, 'WeightMean', 4x, 'IfPop', 4x, 'NofW=0', 4x, 'NAnces', 4x, 'NIdpt', 4x, &
                  & 'NmChg', 4x, 'MaxNm', 4x, 'NmCom', 4x, 'NmSnd', 4x, 'NmRcv')")
            end if
            if(IfPopCtrl) then
               write(291, "(I5.5, 4x, es13.6, 5x, l1, 6x, I6, 4x, I6, 6I9)") NT, MeanWghtNT, IfPopCtrl, &
                  & NumZeroWght, NancestryWalker, NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv
            else if(.not. IfPopCtrl) then
               write(291, "(I5.5, 4x, es13.6, 5x, l1, 6x, I6, 4x, I6     )") NT, MeanWghtNT, IfPopCtrl, &
                  & NumZeroWght, NancestryWalker
            end if
            close(291)
         end if
      end if
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&  
!_________________ Counting calling times and time consumed for population control process ____________________
!%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
      call system_clock(time2)
		TimePopCt = TimePopCt + TimeIntrvl(time1, time2)
      
   end subroutine PopControl
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
!########################################################################################################################
!########################################################################################################################
!############################# Determine whether to do the population control or not ####################################
!############################# Determine whether to do the population control or not ####################################
!########################################################################################################################
!########################################################################################################################

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine IfPopControl(IfPopCtrl, NumZeroWght, MeanWghtNT)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  IfPopControl(IfPopCtrl, NumZeroWght, MeanWghtNT)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to determine whether we need to apply the population control for all the 
!                    random walkers.
! KEYWORDS: Determine whether we need to do population control.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Determine whether we need to do population control
!
!     Input:  (none)   
!
!     Output: IfPopCtrl   --> ==T, need to do population control; ==F, don't need population control.
!             NumZeroWght --> Number of random walker with zero weight;
!             MeanWghtNT  --> The average of weights for all random walkers at this NT.
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
      logical(4) IfPopCtrl
      integer NumZeroWght
      real(rp) MeanWghtNT
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0
      real(rp) WghtMaxm, WghtMinm
      real(rp) Weight(NWalk)
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of Determining IfPopCtrl _______________________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Obtain the weight of random walkers _______________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Count number of random walkers killed by propagation _____________
!________________________________________________________________________________________      
      NumZeroWght = 0
      do I0 = 1, NWalk, +1
         if(WghtProc(I0) <= 0.0_rp) NumZeroWght = NumZeroWght + 1
      enddo
!________________________________________________________________________________________ 	  
!_________________ (1) Obtain the weights used for redistributing walkers _______________
!________________________________________________________________________________________       
      do I0 = 1, NWalk, +1
         Weight(I0) = WghtProc(I0)
      enddo
!**************************************************************************************************	  
!___________________ 1. For the MPI parallel computation case _____________________________________
!**************************************************************************************************
#ifdef MPIPROCESS
!________________________________________________________________________________________ 	  
!_________________ (0) Gather the weights from all random walkers _______________________
!________________________________________________________________________________________
      call MPI_Barrier(acomm, ierr)
      WghtTotl = 0.0_rp
      call MPI_GATHER(Weight(1), NWalk, rp_MPI_REAL, WghtTotl(1), NWalk, rp_MPI_REAL, amstr, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (1) Count number of random walkers killed by propagation _____________
!________________________________________________________________________________________      
      I0 = 0
      call MPI_ALLREDUCE(NumZeroWght, I0, 1, MPI_INTEGER, MPI_SUM, acomm, ierr)
      NumZeroWght = I0
!________________________________________________________________________________________ 	  
!_________________ (2) Max, Min and Mean weights of all random walkers __________________
!_____________________ In this step, the weight average is tuned to 1 ___________________
!________________________________________________________________________________________
      if(amyid == amstr) then
!____________________________________________________________________________ 	  
!________________ [0] Get Mean, Max, Min of weights for all walkers _________
!____________________________________________________________________________
         MeanWghtNT = sum(WghtTotl) / dble(NmWalkAllP)
         WghtMaxm = maxval(WghtTotl(:), 1)
         WghtMinm = minval(WghtTotl(:), 1)
!____________________________________________________________________________ 	  
!________________ [1] Determine IfPopCtrl and MeanWghtPop ___________________
!____________________________________________________________________________
         if( (WghtMaxm/MeanWghtNT > PopWghtMax) .or. (WghtMinm/MeanWghtNT < PopWghtMin) ) then
            IfPopCtrl   = .true.
            MeanWghtPop = MeanWghtNT
         else
            IfPopCtrl   = .false.
            MeanWghtPop = 1.0_rp
         end if
      end if      
!________________________________________________________________________________________ 	  
!_________________ (3) BCast these values to all processes ______________________________
!________________________________________________________________________________________
      call MPI_Barrier(acomm, ierr)
      call MPI_BCAST(IfPopCtrl  , 1, MPI_LOGICAL, amstr, acomm, ierr)
      call MPI_BCAST(MeanWghtPop, 1, rp_MPI_REAL, amstr, acomm, ierr)
!**************************************************************************************************	  
!___________________ 2. For the serial computation case ___________________________________________
!**************************************************************************************************
#else
!________________________________________________________________________________________ 	  
!_________________ (0) Max, Min and Mean weights of all random walkers __________________
!________________________________________________________________________________________
      MeanWghtNT = sum(Weight) / dble(NWalk)
      WghtMaxm = maxval(Weight(:), 1)
      WghtMinm = minval(Weight(:), 1)
!________________________________________________________________________________________ 	  
!_________________ (1) Determine IfPopCtrl and MeanWghtPop ______________________________
!________________________________________________________________________________________
   if( (WghtMaxm/MeanWghtNT > PopWghtMax) .or. (WghtMinm/MeanWghtNT < PopWghtMin) ) then
      IfPopCtrl   = .true.
      MeanWghtPop = MeanWghtNT
   else
      IfPopCtrl   = .false.
      MeanWghtPop = 1.0_rp
   end if
#endif
      
   end subroutine IfPopControl
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
!########################################################################################################################
!########################################################################################################################
!############################# Perform the population control for NT<LTrot case #########################################
!############################# Perform the population control for NT<LTrot case #########################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine PopControlNT(NT, NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PopControlNT(NT, NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the population control process for all the random walkers on 
!                 all processes, within the MPI parallel computations. For NT < LTrot case.
! KEYWORDS: Population control for MPI parallel case.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Perform the MPI parallel or serial calculation of population control process.
!
!     Input:  NT --> The integer index of imaginary-time slice;
!  
!     Output: NIndpt --> Number of independent random walkers after the population control;
!             NumChg --> Number of walkers whose information needs to be changed;
!             MaxNum --> Number of copies of the walker that has the largest weight;
!             NumCom --> Number of walkers needing the sending and receiving processes;
!             NumSnd --> Number of walkers coming from different processes;
!             NumRcv --> Number of walkers receiving from different processes.       
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
      integer NT, NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, I3, I4, Itp0
      integer WalkerIndx, OriginWalk, ProcsIndex, NSendTotal, NRecvTotal
      integer Isw_Cnt, Nsw_Cnt, Isw_tmp, Irp_tmp, Irw_tmp, Irw_Cnt, Nrw_Cnt
      !!!!!!!!!! Arrays of Preparing for the sending-receiving
      integer, allocatable :: WkTableAll(:), WkTabProcs(:)        ! Record the walker origin and destination
      integer, allocatable :: SndRcvProc(:, :), SndRcvAllP(:, :)  ! The send out table on every process
      integer, allocatable :: SendNmCont(:)                       ! Number of walkers sent out from each process
      integer, allocatable :: NumScatter(:)                       ! NSend(i)*2 
      integer, allocatable :: DspScatter(:)                       ! Used in MPI_ScatterV
      !!!!!!!!!! The MPI_ISEND and MPI_RECV related
      integer NmBuf
      character, allocatable :: SendDataBf(:, :)
      character, allocatable :: RecvDataBf(:, :)
      integer, allocatable :: SendRequst(:   )
      integer, allocatable :: SendStatus(:, :)
      integer, allocatable :: RecvStatus(:   )
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of population control for parallel case ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Obtain the walker table according to weights ______________________________
!**************************************************************************************************
#ifdef MPIPROCESS
!________________________________________________________________________________________ 	  
!_________________ (0) Obtain the walker table for all random walkers ___________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         allocate(WkTableAll(NmWalkAllP)); WkTableAll = 0
         call ReDistrbtWalkNT(WghtTotl(1), WkTableAll(1), NIndpt, NumChg, MaxNum)
      else
         allocate(WkTableAll(1)); WkTableAll = 0
      end if
!________________________________________________________________________________________ 	  
!_________________ (1) Scatter the walkers table to all processes _______________________
!________________________________________________________________________________________
      allocate(WkTabProcs(NWalk)); WkTabProcs = 0
      call MPI_Barrier(acomm, ierr)
      call MPI_SCATTER(WkTableAll(1), NWalk, MPI_INTEGER, WkTabProcs(1), NWalk, MPI_INTEGER, amstr, acomm, ierr)
!**************************************************************************************************	  
!___________________ 1. Construct the send out table for following use ____________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Count the number of walkers every process sends out ______________
!________________________________________________________________________________________
      !!!!!!!!!! Allocate array for number of sending walkers on every process
      allocate(SendNmCont(0:anprc-1)); SendNmCont = 0
      !!!!!!!!!! Compute SendNmCont from WkTableAll
      if(amyid == amstr) then
         SendNmCont = 0
         do I0 = 0, anprc-1, +1
            do I1 = 1, NWalk, +1
               WalkerIndx = I1 + I0*NWalk
               OriginWalk = WkTableAll(WalkerIndx)
               ProcsIndex = (OriginWalk - 1) / NWalk
               if(ProcsIndex /= I0) then
                  SendNmCont(ProcsIndex) = SendNmCont(ProcsIndex) + 1
               end if
            enddo
         enddo
      end if
      !!!!!!!!!! Broadcast SendNmCont to all processes
      call MPI_Barrier(acomm, ierr)
      call MPI_BCAST(SendNmCont(0), anprc, MPI_INTEGER, amstr, acomm, ierr)
      !!!!!!!!!! Total number of sending walkers
      NSendTotal = sum(SendNmCont(0:anprc-1))
      !!!!!!!!!! Number of communication during population control == NSendTotal
      NumCom = NSendTotal
!________________________________________________________________________________________ 	  
!_________________ (1) Construct the SndRcvAllP matrix and BCast _________________________
!________________________________________________________________________________________
      !!!!!!!!!! Allocate the array for the send-receive table
      if(NSendTotal > 0) then
         allocate(SndRcvAllP(2, NSendTotal))
         SndRcvAllP = 0
      else
         allocate(SndRcvAllP(1, 1))
         SndRcvAllP = 0
      end if
      !!!!!!!!!! Constructing the send-receive table, SndRcvAllP(1, I2) --> SndRcvAllP(2, I2)
      if(amyid == amstr .and. NSendTotal > 0) then
         I2 = 0
         do I0 = 0, anprc-1, +1
            do I1 = 1, NWalk, +1
               WalkerIndx = I1 + I0*NWalk
               OriginWalk = WkTableAll(WalkerIndx)
               ProcsIndex = (OriginWalk - 1) / NWalk
               if(ProcsIndex /= I0) then
                  I2 = I2 + 1
                  SndRcvAllP(1, I2) = OriginWalk
                  SndRcvAllP(2, I2) = WalkerIndx
               end if
            enddo
         enddo
         if(I2 /= NSendTotal) then
            write(*, "(A)") "=================================================================================="
            write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
            write(*, "(3x, 'PopControlNT: I2 /= NSendTotal!!! I2, NSendTotal = ', 2I6)") I2, NSendTotal
            write(*, "(A)") "=================================================================================="
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) ScatterV the SndRcvAllP matrix to all processes __________________
!________________________________________________________________________________________
      !!!!!!!!!! Prepare for the MPI_SCATTERV process
      allocate(NumScatter(0:anprc-1)); NumScatter = 0
      do I0 = 0, anprc-1, +1
         NumScatter(I0) = SendNmCont(I0) * 2
      enddo
      allocate(DspScatter(0:anprc-1)); DspScatter = 0
      do I0 = 1, anprc-1, +1
         DspScatter(I0) = DspScatter(I0-1) + NumScatter(I0-1)
      enddo
      if(SendNmCont(amyid) > 0) then
         allocate(SndRcvProc(2, SendNmCont(amyid)))
      else
         allocate(SndRcvProc(1, 1))
      end if
      !!!!!!!!!! Carry out the MPI_SCATTERV process
      call MPI_Barrier(acomm, ierr)
      call MPI_SCATTERV(SndRcvAllP(1, 1), NumScatter(0:anprc-1), DspScatter(0:anprc-1), MPI_INTEGER, &
         & SndRcvProc(1, 1), NumScatter(amyid), MPI_INTEGER, amstr, acomm, ierr)
!**************************************************************************************************	  
!___________________ 2. Send and receive the data after the population control ____________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Count the number of complex number for sending ___________________
!________________________________________________________________________________________
      NmBuf = 0
!____________________________________________________________________________ 	  
!________________ [0] For the URght, DRghtVec, VRght matrices _______________
!____________________ And for LogScaleRght coefficients _____________________
!____________________________________________________________________________
      NmBuf = NmBuf + 2*NumNS*NumNS * rp
      NmBuf = NmBuf + 2*NumNS       * rp
      NmBuf = NmBuf + 2*NumNS*NumNS * rp
      NmBuf = NmBuf + 2             * rp
!____________________________________________________________________________ 	  
!________________ [1] For GrnFunct(:, :, :, Iwalk) matrix ___________________
!____________________________________________________________________________
      NmBuf = NmBuf + 2*NumNS*NumNS * rp
!____________________________________________________________________________ 	  
!________________ [2] For auxiliary fields of interactions __________________
!____________________________________________________________________________      
      NmBuf = NmBuf + NumNS*NT * 4
!____________________________________________________________________________ 	  
!________________ [3] For the ancestry link of random walkers _______________
!____________________________________________________________________________
      NmBuf = NmBuf + 1*4
!________________________________________________________________________________________ 	  
!_________________ (1) Pack and send message on every process ___________________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [0] Count number of indepdent sending _____________________
!____________________________________________________________________________ 
      !!!!!!!!!! Obtain Nsw_Cnt (number of non-repeated sending) of every process from SndRcvProc
      Isw_Cnt = 0; Isw_tmp = 0; Irp_tmp = -1
      do I0 = 1, SendNmCont(amyid), +1
         WalkerIndx = SndRcvProc(2, I0)
         ProcsIndex = (WalkerIndx - 1) / NWalk
         if( (SndRcvProc(1, I0) /= Isw_tmp) .or. (ProcsIndex /= Irp_tmp) ) then
            Isw_Cnt = Isw_Cnt + 1
            Isw_tmp = SndRcvProc(1, I0)
            Irp_tmp = ProcsIndex
         end if
      enddo
      Nsw_Cnt = Isw_Cnt
      !!!!!!!!!! Obtain the total number of non-repeated sending for all processes as NumSnd
      NumSnd = Nsw_Cnt; Itp0 = 0
      call MPI_Barrier(acomm, ierr)
      call MPI_ALLREDUCE(NumSnd, Itp0, 1, MPI_INTEGER, MPI_SUM, acomm, ierr)
      NumSnd = Itp0
!____________________________________________________________________________ 	  
!________________ [1] Allocate the arrays for sending _______________________
!____________________________________________________________________________       
      if(Nsw_Cnt > 0) then
         allocate(SendDataBf(NmBuf, Nsw_Cnt))
         allocate(SendRequst(Nsw_Cnt))
         allocate(SendStatus(MPI_STATUS_SIZE, Nsw_Cnt))
      end if
!____________________________________________________________________________ 	  
!________________ [2] Pack and send the information _________________________
!____________________________________________________________________________ 
      Isw_Cnt = 0; Isw_tmp = 0; Irp_tmp = -1
      do I0 = 1, SendNmCont(amyid), +1
         !!!!!!!! The original location of the sending data, the present process
         OriginWalk = SndRcvProc(1, I0)
         I1 = (OriginWalk - 1) / NWalk
         OriginWalk = OriginWalk - I1*NWalk
         if(I1 /= amyid) then
            write(*, "(A)") "=================================================================================="
            write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
            write(*, "(3x, 'PopControlNT: I1 /= amyid!!! I1, amyid = ', 2I6)") I1, amyid
            write(*, "(A)") "=================================================================================="
         end if
         !!!!!!!! The destination of the data sending
         WalkerIndx = SndRcvProc(2, I0)
         ProcsIndex = (WalkerIndx - 1) / NWalk
         WalkerIndx = WalkerIndx - ProcsIndex*NWalk
         !!!!!!!! Pack the data and then send it, avoid twice sending the same data on an single process
         if( (SndRcvProc(1, I0) /= Isw_tmp) .or. (ProcsIndex /= Irp_tmp) ) then
            Isw_Cnt = Isw_Cnt + 1
            Isw_tmp = SndRcvProc(1, I0)
            Irp_tmp = ProcsIndex
            call PopCtrlPackData(NT, NmBuf, OriginWalk, SendDataBf(1, Isw_Cnt))
            call MPI_ISEND(SendDataBf(1, Isw_Cnt), NmBuf, MPI_BYTE, ProcsIndex, WalkerIndx, acomm, &
                  & SendRequst(Isw_Cnt), ierr)
         end if
      enddo
!________________________________________________________________________________________ 	  
!_________________ (2) Receive and unpack message on every process ______________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [0] Count number of indepdent receiving ___________________
!____________________________________________________________________________ 
      !!!!!!!!!! Obtain Nrw_Cnt (number of non-repeated receiving) of every process from WkTabProcs
      Irw_Cnt = 0; Isw_tmp = 0
      do I0 = 1, NWalk, +1
         ProcsIndex = (WkTabProcs(I0) - 1) / NWalk
         if( (ProcsIndex /= amyid) .and. (WkTabProcs(I0) /= Isw_tmp) ) then
            Irw_Cnt = Irw_Cnt + 1
            Isw_tmp = WkTabProcs(I0)
         end if
      enddo
      Nrw_Cnt = Irw_Cnt
      !!!!!!!!!! Obtain the total number of non-repeated receiving for all processes as NumRcv
      NumRcv = Nrw_Cnt; Itp0 = 0
      call MPI_Barrier(acomm, ierr)
      call MPI_ALLREDUCE(NumRcv, Itp0, 1, MPI_INTEGER, MPI_SUM, acomm, ierr)
      NumRcv = Itp0
      !!!!!!!!!! Check whether NumSnd == NumRcv
      if(NumSnd /= NumRcv) then
         write(*, "(A)") "=================================================================================="
         write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
         write(*, "(3x, 'PopControlNT: NumSnd /= NumRcv!!! NumSnd, NumRcv = ', 2I6)") NumSnd, NumRcv
         write(*, "(A)") "=================================================================================="
      end if
!____________________________________________________________________________ 	  
!________________ [1] Allocate the arrays for receiving _____________________
!____________________________________________________________________________       
      if(Nrw_Cnt > 0) then
         allocate(RecvDataBf(NmBuf, Nrw_Cnt))
         allocate(RecvStatus(MPI_STATUS_SIZE))
      end if
!____________________________________________________________________________ 	  
!________________ [2] Receive and unpack the information ____________________
!____________________________________________________________________________ 
      Irw_Cnt = 0; Isw_tmp = 0; Irw_tmp = 0
      do I0 = 1, NWalk, +1
         WalkerIndx = I0 + amyid*NWalk
         OriginWalk = WkTabProcs(I0)
         ProcsIndex = (OriginWalk - 1) / NWalk
         if(WalkerIndx /= OriginWalk) then
            if(amyid == ProcsIndex) then
               OriginWalk = OriginWalk - ProcsIndex*NWalk
               call CopyWk_Iw_to_Jw(NT, OriginWalk, I0)
            else
               if(OriginWalk == Isw_tmp) then
                  call CopyWk_Iw_to_Jw(NT, Irw_tmp, I0)
               else
                  Irw_Cnt = Irw_Cnt + 1
                  Isw_tmp = OriginWalk
                  Irw_tmp = I0
                  call MPI_RECV(RecvDataBf(1, Irw_Cnt), NmBuf, MPI_BYTE, ProcsIndex, I0, acomm, RecvStatus, ierr)
                  call PopCtrlUnPkData(NT, NmBuf, I0, RecvDataBf(1, Irw_Cnt))
               end if
            end if
         end if
         WghtProc(I0) = 1.0_rp
         Log_Wght(I0) = 0.0_rp
      enddo
!________________________________________________________________________________________ 	  
!_________________ (2) The finalization for the sending and receiving procedure _________
!________________________________________________________________________________________
      !!!!!!!!!! First finish the non-blocking sending procedure
      if(Nsw_Cnt > 0) call MPI_WAITALL(Nsw_Cnt, SendRequst, SendStatus, ierr)
      !!!!!!!!!! Synchronize all the processes
      call MPI_Barrier(acomm, ierr)      
!**************************************************************************************************	  
!___________________ 3. For the serial computation case ___________________________________________
!**************************************************************************************************
#else
!________________________________________________________________________________________ 	  
!_________________ (0) Redistribute the random walkers by sampling ______________________
!________________________________________________________________________________________
      allocate(WkTabProcs(NWalk)); WkTabProcs = 0
      call ReDistrbtWalkNT(WghtProc(1), WkTabProcs(1), NIndpt, NumChg, MaxNum)
!________________________________________________________________________________________ 	  
!_________________ (1) Population control by simple copying walkers _____________________
!________________________________________________________________________________________
      do I0 = 1, NWalk, +1
         if(I0 /= WkTabProcs(I0)) call CopyWk_Iw_to_Jw(NT, WkTabProcs(I0), I0)
         WghtProc(I0) = 1.0_rp
         Log_Wght(I0) = 0.0_rp
      enddo
#endif
!**************************************************************************************************	  
!___________________ 4. Deallocate all the used matrices __________________________________________
!************************************************************************************************** 
      if(allocated(WkTableAll)) deallocate(WkTableAll)
      if(allocated(WkTabProcs)) deallocate(WkTabProcs)
      if(allocated(SndRcvProc)) deallocate(SndRcvProc)
      if(allocated(SndRcvAllP)) deallocate(SndRcvAllP)
      if(allocated(SendNmCont)) deallocate(SendNmCont)
      if(allocated(NumScatter)) deallocate(NumScatter)
      if(allocated(DspScatter)) deallocate(DspScatter)
      if(allocated(SendDataBf)) deallocate(SendDataBf)
      if(allocated(RecvDataBf)) deallocate(RecvDataBf)
      if(allocated(SendRequst)) deallocate(SendRequst)
      if(allocated(SendStatus)) deallocate(SendStatus)
      if(allocated(RecvStatus)) deallocate(RecvStatus)
      
   end subroutine PopControlNT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine GetNAncestry()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GetNAncestry()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is obtain the number of ancestry random walkers after the population control process.
! KEYWORDS: Obtain the number of ancestry walkers.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Get the ancestry walkers.
!
!     Input:  (none); Output: (none).
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
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________      
      integer I0, Itp0
      integer, allocatable :: AncestryAll(:)
!______________________________________________________________________________________________________________	  
!________________________ Main calculations of Calculating Number of Ancestry walkers _________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For the case of MPI parallel simulations __________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Allocate vector for the ancestry collecttion _____________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         allocate(AncestryAll(NmWalkAllP))
      else
         allocate(AncestryAll(1))
      end if
      AncestryAll = 0
!________________________________________________________________________________________ 	  
!_________________ (1) Gather the ancestry vectors from all processes ___________________
!________________________________________________________________________________________
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr)
      call MPI_GATHER(AncestryLink(1), NWalk, MPI_INTEGER, AncestryAll(1), NWalk, MPI_INTEGER, amstr, acomm, ierr)
#else
      AncestryAll = AncestryLink
#endif
!________________________________________________________________________________________ 	  
!_________________ (2) Calculate NancestryWalker only on the main process _______________
!________________________________________________________________________________________
      if(amyid == amstr) then
!____________________________________________________________________________ 	  
!________________ [0] Sort the AncestryAll vector ___________________________
!____________________________________________________________________________
         call QuckSortI(1, NmWalkAllP, AncestryAll(1))
!____________________________________________________________________________ 	  
!________________ [1] Calculate NancestryWalker _____________________________
!____________________________________________________________________________
         NancestryWalker = 1; Itp0 = AncestryAll(1)
         do I0 = 2, NmWalkAllP, +1
            if(AncestryAll(I0) /= Itp0) then
               NancestryWalker = NancestryWalker + 1
               Itp0 = AncestryAll(I0)
            end if
         enddo
      end if
!________________________________________________________________________________________ 	  
!_________________ (3) Deallocate the used vector here __________________________________
!________________________________________________________________________________________
      if(allocated(AncestryAll)) deallocate(AncestryAll)
      
   end subroutine GetNAncestry
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ReDistrbtWalkNT(Weight, WkIndTable, NIndpt, NumChg, MaxNum)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ReDistrbtWalkNT(Weight, WkIndTable, NIndpt, NumChg, MaxNum)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to obtain the walker table for all the random walkers, after performing the 
!                    population control process. For NT < LTrot case.
! KEYWORDS: Random Walker table after the population control process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Obtain the random walker table after the population control.
!
!     As for the output, we have the vector WkIndTable(I1) = I2. It means that for the I1-th random walker, it is 
!         now a copy of the original I2-th random walker.      
!
!     Input:  Weight --> The weights of all NmWalkAllP random walkers;
!
!     Output: WkIndTable --> The walker table containing integer indexes; 
!             NIndpt     --> Number of independent random walkers after population control process;
!             NumChg     --> Number of random walkers whose information needs to be changed;
!             MaxNum     --> Biggest number of random walkers.
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
      integer NIndpt, NumChg, MaxNum
      integer WkIndTable(NmWalkAllP)
      real(rp) Weight(NmWalkAllP)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, I3, I4
      integer StartWkInd, WalkerIndx
      real(rp) Rtp0, Rtp1
      integer, allocatable :: NumWalkers(:)
      real(rp), allocatable :: NormWeight(:)
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of walker table from population control ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Allocate arrays used in this subroutine ___________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) The normalized weights of all random walkers _____________________
!________________________________________________________________________________________
      allocate(NormWeight(NmWalkAllP))
      NormWeight = Weight / sum(Weight(1:NmWalkAllP))
!________________________________________________________________________________________ 	  
!_________________ (1) Number of each walker after population control ___________________
!________________________________________________________________________________________
      allocate(NumWalkers(NmWalkAllP))
      NumWalkers = 0
!**************************************************************************************************	  
!___________________ 1. Perform the population control ____________________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Initialization for the population control ________________________
!________________________________________________________________________________________
      !!!!!!!!!! Determine the number of random walkers using the weights  
      I1 = 1; Rtp0 = NormWeight(1)
      do I0 = 1, NmWalkAllP, +1
         Rtp1 = ( dble(I0-1) + spring_sfmt_stream() ) / dble(NmWalkAllP)
         do while(Rtp1 > Rtp0)
            I1 = I1 + 1
            Rtp0 = Rtp0 + NormWeight(I1)
         enddo
         NumWalkers(I1) = NumWalkers(I1) + 1
         WkIndTable(I1) = I1
      enddo
      !!!!!!!!!! The number of independent random walkers
      NIndpt = 0
      do I0 = 1, NmWalkAllP, +1
         if(NumWalkers(I0) > 0) NIndpt = NIndpt + 1
      enddo
      !!!!!!!!!! The number of the walkers which neeeds to be changed
      NumChg = 0
      do I0 = 1, NmWalkAllP, +1
         if(WkIndTable(I0) /= I0) NumChg = NumChg + 1
      enddo
      !!!!!!!!!! The number of the walker with largest weight
      MaxNum = maxval(NumWalkers(:), 1)
!________________________________________________________________________________________ 	  
!_________________ (1) Adjust walkers on the same process _______________________________
!________________________________________________________________________________________
      do I0 = 0, anprc-1, +1
         StartWkInd = 1 + I0*NWalk
         do I1 = 1, NWalk, +1
            WalkerIndx = I1 + I0*NWalk
            if(NumWalkers(WalkerIndx) == 0) then
               do I2 = StartWkInd, (I0+1)*NWalk, +1
                  if(NumWalkers(I2) > 1) then
                     NumWalkers(I2) = NumWalkers(I2) - 1
                     NumWalkers(WalkerIndx) = NumWalkers(WalkerIndx) + 1
                     WkIndTable(WalkerIndx) = I2
                     exit
                  end if
               enddo
               StartWkInd = I2
            end if
         enddo
      enddo
!________________________________________________________________________________________ 	  
!_________________ (2) Adjust walkers on different processes ____________________________
!________________________________________________________________________________________       
      if(anprc > 1) then
         StartWkInd = 1
         do I1 = 1, NmWalkAllP, +1
            if(NumWalkers(I1) == 0) then
               do I2 = StartWkInd, NmWalkAllP, +1
                  if(NumWalkers(I2) > 1) then
                     NumWalkers(I2) = NumWalkers(I2) - 1
                     NumWalkers(I1) = NumWalkers(I1) + 1
                     WkIndTable(I1) = I2
                     exit
                  end if
               enddo
               StartWkInd = I2
            end if
         enddo
      end if
!________________________________________________________________________________________ 	  
!_________________ (3) Perform a check for NumWalkers vector ____________________________
!________________________________________________________________________________________ 
      do I0 = 1, NmWalkAllP, +1
         if(NumWalkers(I0) /= 1) then
            write(*, "(A)") "=================================================================================="
            write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
            write(*, "(3x, 'ReDistrbtWalkNT: NumWalkers(I0) /= 1!!! I0, NumWalkers(I0) = ', 2I6)") I0, NumWalkers(I0)
            write(*, "(A)") "=================================================================================="
         end if
      enddo
!**************************************************************************************************	  
!___________________ 2. Deallocate the matrices used in this subroutine ___________________________
!**************************************************************************************************
      if(allocated(NumWalkers)) deallocate(NumWalkers)
      if(allocated(NormWeight)) deallocate(NormWeight)
      
   end subroutine ReDistrbtWalkNT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine PopCtrlPackData(NT, NmBuf, WalkerIndx, SendDataBf)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PopCtrlPackData(NT, NmBuf, WalkerIndx, SendDataBf)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to pack all the data of the WalkerIndx-th walker, preparing for communicating.
! KEYWORDS: Pack data for the population process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Pack all the data information of WalkerIndx-th walker.
!
!     Input:  NT         --> The integer index of imaginary-time slice;
!             NmBuf      --> Number of complex numbers in the pack;
!             WalkerIndx --> Index integer for the random walker;
!
!     Output: SendDataBf --> The packed data.
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
      integer NT, NmBuf, WalkerIndx
      character SendDataBf(NmBuf)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Posit
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of Data packing process ________________________________________
!______________________________________________________________________________________________________________
#ifdef MPIPROCESS
!**************************************************************************************************	  
!___________________ 0. Pack all the data that needs to be transferred ____________________________
!**************************************************************************************************
      Posit = 0
!________________________________________________________________________________________ 	  
!_________________ (0) The URght, DRghtVec, VRght matrices and LogScaleRght _____________
!________________________________________________________________________________________   
      call MPI_PACK( URght(1, 1, 1, WalkerIndx), 2*NumNS*NumNS, rp_MPI_REAL, SendDataBf(1), NmBuf, posit, acomm, ierr)
      call MPI_PACK( DRghtVec(1, 1, WalkerIndx), 2*NumNS      , rp_MPI_REAL, SendDataBf(1), NmBuf, posit, acomm, ierr)
      call MPI_PACK( VRght(1, 1, 1, WalkerIndx), 2*NumNS*NumNS, rp_MPI_REAL, SendDataBf(1), NmBuf, posit, acomm, ierr)
      call MPI_PACK(LogScaleRght(1, WalkerIndx), 2            , rp_MPI_REAL, SendDataBf(1), NmBuf, posit, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (1) The GrnFunct(:, :, :, Iwalk) matrix ______________________________
!________________________________________________________________________________________
      call MPI_PACK(GrnFunct(1, 1, 1, WalkerIndx), 2*NumNS*NumNS, rp_MPI_REAL, SendDataBf(1), NmBuf, posit, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (2) The auxiliary fields for interactions ____________________________
!________________________________________________________________________________________ 
      call MPI_PACK(IsingbU(1, 1, WalkerIndx), NumNS*NT, MPI_INTEGER, SendDataBf(1), NmBuf, posit, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (3) The ancestry link for all the random walkers _____________________
!________________________________________________________________________________________       
      call MPI_PACK(AncestryLink(WalkerIndx), 1, MPI_INTEGER, SendDataBf(1), NmBuf, posit, acomm, ierr)
!**************************************************************************************************	  
!___________________ 2. Check the amount of data in the buf area __________________________________
!**************************************************************************************************
      if(posit /= NmBuf) then
         write(*, "(A)") "=================================================================================="
         write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
         write(*, "(3x, 'PopCtrlPackData: posit /= NmBuf!!! amyid, posit, NmBuf = ', 3I6)") amyid, posit, NmBuf
         write(*, "(A)") "=================================================================================="
      end if
#endif
      
   end subroutine PopCtrlPackData
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine PopCtrlUnPkData(NT, NmBuf, WalkerIndx, RecvDataBf)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PopCtrlUnPkData(NT, NmBuf, WalkerIndx, RecvDataBf)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to unpack all the data of the WalkerIndx-th walker, preparing for communicating.
! KEYWORDS: Pack data for the population process.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Unpack all the data information of WalkerIndx-th walker.
!
!     Input:  NT         --> The integer index of imaginary-time slice;
!             NmBuf      --> Number of complex numbers in the pack;
!             WalkerIndx --> Index integer for the random walker;
!             RecvDataBf --> The data to be unpacked.
!
!     Output: (none).
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
      integer NT, NmBuf, WalkerIndx
      character RecvDataBf(NmBuf)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Posit
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of Data packing process ________________________________________
!______________________________________________________________________________________________________________
#ifdef MPIPROCESS
!**************************************************************************************************	  
!___________________ 0. UnPack all the data that has been transferred _____________________________
!**************************************************************************************************
      posit = 0
!________________________________________________________________________________________ 	  
!_________________ (0) The URght, DRghtVec, VRght matrices and LogScaleRght _____________
!________________________________________________________________________________________   
      call MPI_UNPACK(RecvDataBf(1), NmBuf, posit,  URght(1, 1, 1, WalkerIndx), 2*NumNS*NumNS, rp_MPI_REAL, acomm, ierr)
      call MPI_UNPACK(RecvDataBf(1), NmBuf, posit,  DRghtVec(1, 1, WalkerIndx), 2*NumNS      , rp_MPI_REAL, acomm, ierr)
      call MPI_UNPACK(RecvDataBf(1), NmBuf, posit,  VRght(1, 1, 1, WalkerIndx), 2*NumNS*NumNS, rp_MPI_REAL, acomm, ierr)
      call MPI_UNPACK(RecvDataBf(1), NmBuf, posit, LogScaleRght(1, WalkerIndx), 2            , rp_MPI_REAL, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (1) The GrnFunct(:, :, :, Iwalk) matrix ______________________________
!________________________________________________________________________________________
      call MPI_UNPACK(RecvDataBf(1), NmBuf, posit, GrnFunct(1, 1, 1, WalkerIndx), 2*NumNS*NumNS, rp_MPI_REAL, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (2) The auxiliary fields for interactions ____________________________
!________________________________________________________________________________________
      call MPI_UNPACK(RecvDataBf(1), NmBuf, posit, IsingbU(1, 1, WalkerIndx), NumNS*NT, MPI_INTEGER, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (3) The ancestry link for all the random walkers _____________________
!________________________________________________________________________________________       
      call MPI_UNPACK(RecvDataBf(1), NmBuf, posit, AncestryLink(WalkerIndx), 1, MPI_INTEGER, acomm, ierr)
!**************************************************************************************************	  
!___________________ 1. Check the amount of data in the buf area __________________________________
!**************************************************************************************************
      if(posit /= NmBuf) then
         write(*, "(A)") "=================================================================================="
         write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
         write(*, "(3x, 'PopCtrlUnPkData: posit /= NmBuf!!! amyid, posit, NmBuf = ', 3I6)") amyid, posit, NmBuf
         write(*, "(A)") "=================================================================================="
      end if
#endif
      
   end subroutine PopCtrlUnPkData
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine CopyWk_Iw_to_Jw(NT, Iw, Jw)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CopyWk_Iw_to_Jw(Iw, Jw)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to copy the data information of Iw-th walker to those of Jw-th walker.
! KEYWORDS: Change data information.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Copy the Iw-th walker to the Jw-th walker.
!
!     Input:  NT --> The integer index of imaginary-time slice;
!             Iw --> Integer index of random walker to be changed;
!             Jw --> Integer index of random walker where the data comes from.
!
!     Output: (none).
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
      integer NT, Iw, Jw
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of Data changing process _______________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. For the spin decoupled case _______________________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) The URght, DRghtVec, VRght matrices and LogScaleRght _____________
!________________________________________________________________________________________
      call dMat_Copy_QMC( URght(1, 1, 1, Iw),    URght(1, 1, 1, Jw)   )
      call dcopy(2*NumNS, DRghtVec(1, 1, Iw), 1, DRghtVec(1, 1, Jw), 1)
      call dMat_Copy_QMC( VRght(1, 1, 1, Iw),    VRght(1, 1, 1, Jw)   )
      LogScaleRght(1:NmSpn, Jw) = LogScaleRght(1:NmSpn, Iw)
!________________________________________________________________________________________ 	  
!_________________ (1) The GrnFunct(:, :, :, Iwalk) matrix ______________________________
!________________________________________________________________________________________
      call dMat_Copy_QMC(GrnFunct(1, 1, 1, Iw), GrnFunct(1, 1, 1, Jw))
!________________________________________________________________________________________ 	  
!_________________ (2) For the auxiliary fields for interactions ________________________
!________________________________________________________________________________________
      IsingbU(1:NumNS, 1:NT, Jw) = IsingbU(1:NumNS, 1:NT, Iw)
!________________________________________________________________________________________ 	  
!_________________ (3) The ancestry link for all the random walkers _____________________
!________________________________________________________________________________________
      AncestryLink(Jw) = AncestryLink(Iw)
      
   end subroutine CopyWk_Iw_to_Jw
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


!########################################################################################################################
!########################################################################################################################
!############################# Perform the population control for NT==LTrot case ########################################
!############################# Perform the population control for NT==LTrot case ########################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine PopCtrlBetaT(NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  PopCtrlBetaT(NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the population control process for all the random walkers on 
!                 all processes, within the MPI parallel computations. For NT == LTrot case.
! KEYWORDS: Population control for MPI parallel case.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Perform the population control process with MPI parallel calculation and serial computation.
!
!     Input:  (none)   
!  
!     Output: NIndpt --> Number of independent random walkers after the population control;
!             NumChg --> Number of walkers whose information needs to be changed;
!             MaxNum --> Number of copies of the walker that has the largest weight;
!             NumCom --> Number of walkers needing the sending and receiving processes;
!             NumSnd --> Number of walkers coming from different processes;
!             NumRcv --> Number of walkers receiving from different processes.        
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
      integer NIndpt, NumChg, MaxNum, NumCom, NumSnd, NumRcv
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, Iwalk, Itp0, NAddWk
      integer StartWkInd, WalkerIndx, OriginWalk, ProcsIndex
      integer NSendTotal, NRecvTotal
      !!!!!!!!!! Arrays of Preparing for the sending-receiving
      integer , allocatable :: WkTableAll(:), WkTabProcs(:)        ! Record the walker origin and destination
      real(rp), allocatable :: WlkWghtAll(:), WlkWghtPrc(:)        ! Record the walker weight
      integer, allocatable :: AllIdptWkInd(:)
      integer, allocatable :: SndRcvProc(:, :), SndRcvAllP(:, :)   ! The send out table on every process
      integer, allocatable :: RecvNmCont(:)                        ! Number of walkers sent out from each process
      integer, allocatable :: SendNmCont(:)                        ! Number of walkers sent out from each process
      integer, allocatable :: NumScatter(:)                        ! NSend(i)*2 
      integer, allocatable :: DspScatter(:)                        ! Used in MPI_ScatterV
      !!!!!!!!!! The MPI_ISEND and MPI_RECV related
      integer NmBuf
      character, allocatable :: SendDataBf(:, :)
      character, allocatable :: RecvDataBf(:, :)
      integer, allocatable :: SendRequst(:   )
      integer, allocatable :: SendStatus(:, :)
      integer, allocatable :: RecvStatus(:   )
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of population control for parallel case ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Obtain the walker table according to weights ______________________________
!**************************************************************************************************
#ifdef MPIPROCESS
!________________________________________________________________________________________ 	  
!_________________ (0) Obtain the walker table for all random walkers ___________________
!________________________________________________________________________________________
      !!!!!!!!!! Obtain information of the population control
      if(amyid == amstr) then
         allocate(WkTableAll(NmWalkAllP)); WkTableAll = 0
         allocate(WlkWghtAll(NmWalkAllP)); WlkWghtAll = 0.0_rp
         call ReDistrbWkBetaT(WghtTotl(1), WkTableAll(1), WlkWghtAll(1), NAddWk, NIndpt, NumChg, MaxNum)
      else
         allocate(WkTableAll(1)); WkTableAll = 0
         allocate(WlkWghtAll(1)); WlkWghtAll = 0.0_rp
      end if
      !!!!!!!!!! BCAST NWkBt and NAddWk to all processes
      call MPI_Barrier(acomm, ierr)
      call MPI_BCAST(NWkBt , 1, MPI_INTEGER, amstr, acomm, ierr)
      call MPI_BCAST(NAddWk, 1, MPI_INTEGER, amstr, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (1) Scatter the walkers table to all processes _______________________
!________________________________________________________________________________________
      allocate(WkTabProcs(NWalk)); WkTabProcs = 0
      allocate(WlkWghtPrc(NWalk)); WlkWghtPrc = 0.0_rp
      call MPI_Barrier(acomm, ierr)
      call MPI_SCATTER(WkTableAll(1), NWalk, MPI_INTEGER, WkTabProcs(1), NWalk, MPI_INTEGER, amstr, acomm, ierr)
      call MPI_SCATTER(WlkWghtAll(1), NWalk, rp_MPI_REAL, WlkWghtPrc(1), NWalk, rp_MPI_REAL, amstr, acomm, ierr)
!**************************************************************************************************	  
!___________________ 1. Construct the send out table for following use ____________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Count the number of walkers every process sends out ______________
!________________________________________________________________________________________
      !!!!!!!!!! Allocate array for number of sending and receiving walkers on every process
      allocate(SendNmCont(0:anprc-1)); SendNmCont = 0
      allocate(RecvNmCont(0:anprc-1)); RecvNmCont = 0
      !!!!!!!!!! Compute SendNmCont from WkTableAll
      if(amyid == amstr) then
         SendNmCont = 0
         do I0 = 0, anprc-1, +1
            Itp0 = NWkBt - merge(1, 0, I0 <= NAddWk-1)
            do I1 = 1, Itp0, +1
               WalkerIndx = I1 + I0*NWalk
               OriginWalk = WkTableAll(WalkerIndx)
               ProcsIndex = (OriginWalk - 1) / NWalk
               if(ProcsIndex /= I0) then
                  SendNmCont(ProcsIndex) = SendNmCont(ProcsIndex) + 1
                  RecvNmCont(I0) = RecvNmCont(I0) + 1
               end if
            enddo
         enddo
      end if
      !!!!!!!!!! Broadcast SendNmCont to all processes
      call MPI_Barrier(acomm, ierr)
      call MPI_BCAST(SendNmCont(0), anprc, MPI_INTEGER, amstr, acomm, ierr)
      call MPI_BCAST(RecvNmCont(0), anprc, MPI_INTEGER, amstr, acomm, ierr)
      !!!!!!!!!! Total number of sending walkers
      NSendTotal = sum(SendNmCont(0:anprc-1))
      NRecvTotal = sum(RecvNmCont(0:anprc-1))
      if(amyid == amstr .and. NSendTotal /= NRecvTotal) then
         write(*, "(A)") "=================================================================================="
         write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
         write(*, "(3x, 'PopCtrlBetaT: NSendTotal /= NRecvTotal!!! NSendTotal, NRecvTotal = ', 2I6)") &
            & NSendTotal, NRecvTotal
         write(*, "(A)") "=================================================================================="
      end if
      !!!!!!!!!! Number of communication, sending and receiving during population control == NSendTotal
      NumCom = NSendTotal
      NumSnd = NSendTotal
      NumRcv = NSendTotal
!________________________________________________________________________________________ 	  
!_________________ (1) Obtain table of the indexes of walkers on present process ________
!________________________________________________________________________________________
      !!!!!!!!!! The indexes of NWkBt independent walkers on present process
      IdptWkIndx = -10000
      Itp0 = NWkBt - merge(1, 0, amyid <= NAddWk-1) - RecvNmCont(amyid)
      do I1 = 1, Itp0, +1
         Iwalk = mod(WkTabProcs(I1)-1, NWalk) + 1
         IdptWkIndx(I1) = Iwalk
      enddo
      !!!!!!!!!! Set up the table recording the unoccupied random walkers, and fill the IdptWkIndx table
      if(Itp0 < NWkBt) then
         I2 = 0; StartWkInd = 1
         do I0 = 1, NWalk, +1
            do I1 = StartWkInd, Itp0, +1
               Iwalk = mod(WkTabProcs(I1)-1, NWalk) + 1
               if(Iwalk == I0) then
                  StartWkInd = I1 + 1
                  go to 299
               end if
            enddo
            I2 = I2 + 1; IdptWkIndx(Itp0+I2) = I0
            if(Itp0+I2 >= NWkBt) exit
299         continue
         enddo
      end if
      !!!!!!!!!! The weights of the NWkBt independent walkers on present process
      WghtProc(1:NWalk) = 0.0_rp; Log_Wght(1:NWalk) = -1.0E+100_rp
      do I1 = 1, NWkBt, +1
         Iwalk = IdptWkIndx(I1)
         if(Iwalk < 0) then
            write(*, "(A)") "=================================================================================="
            write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
            write(*, "(3x, 'PopCtrlBetaT: Iwalk < 0 Invalid!!! NWkBt, NAddWk, amyid, I1, Iwalk = ', 5I6)") &
               & NWkBt, NAddWk, amyid, I1, Iwalk
            write(*, "(A)") "=================================================================================="
         end if
         WghtProc(Iwalk) = WlkWghtPrc(I1)
         Log_Wght(Iwalk) = log(WghtProc(I1))
      enddo
      !!!!!!!!!! Gather all IdptWkIndx to the amstr process for following use
      if(amyid == amstr) then
         allocate(AllIdptWkInd(NmWkBtAllP)); AllIdptWkInd = -10000
      else
         allocate(AllIdptWkInd(1)); AllIdptWkInd = -10000
      end if
      call MPI_Barrier(acomm, ierr)
      call MPI_GATHER(IdptWkIndx(1), NWkBt, MPI_INTEGER, AllIdptWkInd(1), NWkBt, MPI_INTEGER, amstr, acomm, ierr)
!________________________________________________________________________________________ 	  
!_________________ (2) Construct the SndRcvAllP matrix and BCast ________________________
!________________________________________________________________________________________
      !!!!!!!!!! Allocate the array for the send-receive table
      if(NSendTotal > 0) then
         allocate(SndRcvAllP(2, NSendTotal))
         SndRcvAllP = 0
      else
         allocate(SndRcvAllP(1, 1))
         SndRcvAllP = 0
      end if
      !!!!!!!!!! Constructing the send-receive table, SndRcvAllP(1, I2) --> SndRcvAllP(2, I2)
      if(amyid == amstr .and. NSendTotal > 0) then
         I2 = 0
         do I0 = 0, anprc-1, +1
            Itp0 = NWkBt - merge(1, 0, I0 <= NAddWk-1)
            do I1 = 1, Itp0, +1
               WalkerIndx = AllIdptWkInd(I1+I0*NWkBt) + I0*NWalk
               OriginWalk = WkTableAll(I1+I0*NWalk)
               ProcsIndex = (OriginWalk - 1) / NWalk
               if(ProcsIndex /= I0) then
                  I2 = I2 + 1
                  SndRcvAllP(1, I2) = OriginWalk
                  SndRcvAllP(2, I2) = WalkerIndx
               end if
            enddo
         enddo
         if(I2 /= NSendTotal) then
            write(*, "(A)") "=================================================================================="
            write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
            write(*, "(3x, 'PopCtrlBetaT: I2 /= NSendTotal!!! I2, NSendTotal = ', 2I6)") I2, NSendTotal
            write(*, "(A)") "=================================================================================="
         end if
      end if
!________________________________________________________________________________________ 	  
!_________________ (3) ScatterV the SndRcvAllP matrix to all processes __________________
!________________________________________________________________________________________
      !!!!!!!!!! Prepare for the MPI_SCATTERV process
      allocate(NumScatter(0:anprc-1)); NumScatter = 0
      do I0 = 0, anprc-1, +1
         NumScatter(I0) = SendNmCont(I0) * 2
      enddo
      allocate(DspScatter(0:anprc-1)); DspScatter = 0
      do I0 = 1, anprc-1, +1
         DspScatter(I0) = DspScatter(I0-1) + NumScatter(I0-1)
      enddo
      if(SendNmCont(amyid) > 0) then
         allocate(SndRcvProc(2, SendNmCont(amyid)))
      else
         allocate(SndRcvProc(1, 1))
      end if
      !!!!!!!!!! Carry out the MPI_SCATTERV process
      call MPI_Barrier(acomm, ierr)
      call MPI_SCATTERV(SndRcvAllP(1, 1), NumScatter(0:anprc-1), DspScatter(0:anprc-1), MPI_INTEGER, &
         & SndRcvProc(1, 1), NumScatter(amyid), MPI_INTEGER, amstr, acomm, ierr)
!**************************************************************************************************	  
!___________________ 2. Send and receive the data after the population control ____________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Count the number of complex number for sending ___________________
!________________________________________________________________________________________
      NmBuf = 0
!____________________________________________________________________________ 	  
!________________ [0] For the URght, DRghtVec, VRght matrices _______________
!____________________ And for LogScaleRght coefficients _____________________
!____________________________________________________________________________
      NmBuf = NmBuf + 2*NumNS*NumNS * rp
      NmBuf = NmBuf + 2*NumNS       * rp
      NmBuf = NmBuf + 2*NumNS*NumNS * rp
      NmBuf = NmBuf + 2             * rp
!____________________________________________________________________________ 	  
!________________ [1] For GrnFunct(:, :, :, Iwalk) matrix ___________________
!____________________________________________________________________________
      NmBuf = NmBuf + 2*NumNS*NumNS * rp
!____________________________________________________________________________ 	  
!________________ [2] For auxiliary fields of interactions __________________
!____________________________________________________________________________      
      NmBuf = NmBuf + NumNS*LTrot * 4
!____________________________________________________________________________ 	  
!________________ [3] For the ancestry link of random walkers _______________
!____________________________________________________________________________
      NmBuf = NmBuf + 1*4
!________________________________________________________________________________________ 	  
!_________________ (1) Pack and send message on every process ___________________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [0] Allocate the arrays for sending _______________________
!____________________________________________________________________________
      if(SendNmCont(amyid) > 0) then
         allocate(SendDataBf(NmBuf, SendNmCont(amyid)))
         allocate(SendRequst(SendNmCont(amyid)))
         allocate(SendStatus(MPI_STATUS_SIZE, SendNmCont(amyid)))
      end if
!____________________________________________________________________________ 	  
!________________ [1] Pack and send the information _________________________
!____________________________________________________________________________
      do I0 = 1, SendNmCont(amyid), +1
         !!!!!!!! The original location of the sending data, the present process
         OriginWalk = SndRcvProc(1, I0)
         I1 = (OriginWalk - 1) / NWalk
         OriginWalk = OriginWalk - I1*NWalk
         if(I1 /= amyid) then
            write(*, "(A)") "=================================================================================="
            write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
            write(*, "(3x, 'PopCtrlBetaT: I1 /= amyid!!! I1, amyid = ', 2I6)") I1, amyid
            write(*, "(A)") "=================================================================================="
         end if
         !!!!!!!! The destination of the data sending
         WalkerIndx = SndRcvProc(2, I0)
         ProcsIndex = (WalkerIndx - 1) / NWalk
         WalkerIndx = WalkerIndx - ProcsIndex*NWalk
         !!!!!!!! Pack the data and then send it, avoid twice sending the same data on an single process
         call PopCtrlPackData(LTrot, NmBuf, OriginWalk, SendDataBf(1, I0))
         call MPI_ISEND(SendDataBf(1, I0), NmBuf, MPI_BYTE, ProcsIndex, WalkerIndx, acomm, SendRequst(I0), ierr)
      enddo
!________________________________________________________________________________________ 	  
!_________________ (2) Receive and unpack message on every process ______________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [0] Allocate the arrays for receiving _____________________
!____________________________________________________________________________
      if(RecvNmCont(amyid) > 0) then
         allocate(RecvDataBf(NmBuf, RecvNmCont(amyid)))
         allocate(RecvStatus(MPI_STATUS_SIZE))
      end if
!____________________________________________________________________________ 	  
!________________ [1] Receive and unpack the information ____________________
!____________________________________________________________________________
      !!!!!!!!!! For every process, receive and unpack the data
      Itp0 = NWkBt - merge(1, 0, amyid <= NAddWk-1); I1 = 0
      do I0 = 1, Itp0, +1
         WalkerIndx = IdptWkIndx(I0) + amyid*NWalk
         OriginWalk = WkTabProcs(I0)
         ProcsIndex = (OriginWalk - 1) / NWalk
         if(ProcsIndex /= amyid) then
            I1 = I1 + 1
            call MPI_RECV(RecvDataBf(1, I1), NmBuf, MPI_BYTE, ProcsIndex, IdptWkIndx(I0), acomm, RecvStatus, ierr)
            call PopCtrlUnPkData(LTrot, NmBuf, IdptWkIndx(I0), RecvDataBf(1, I1))
         end if
      enddo
      !!!!!!!!!! For amyid<=NAddWk-1 case, copy the last random walker
      if(amyid <= NAddWk-1) then
         call CopyWk_Iw_to_Jw(LTrot, IdptWkIndx(NWkBt-1), IdptWkIndx(NWkBt))
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) The finalization for the sending and receiving procedure _________
!________________________________________________________________________________________
      !!!!!!!!!! First finish the non-blocking sending procedure
      if(SendNmCont(amyid) > 0) call MPI_WAITALL(SendNmCont(amyid), SendRequst, SendStatus, ierr)
      !!!!!!!!!! Synchronize all the processes
      call MPI_Barrier(acomm, ierr)      
!**************************************************************************************************	  
!___________________ 3. For the serial computation case ___________________________________________
!**************************************************************************************************
#else
!________________________________________________________________________________________ 	  
!_________________ (0) Redistribute the random walkers by sampling ______________________
!________________________________________________________________________________________
      allocate(WkTabProcs(NWalk)); WkTabProcs = 0
      allocate(WlkWghtPrc(NWalk)); WlkWghtPrc = 0.0_rp
      call ReDistrbWkBetaT(WghtProc(1), WkTabProcs(1), WlkWghtPrc(1), NAddWk, NIndpt, NumChg, MaxNum)
!________________________________________________________________________________________ 	  
!_________________ (1) Set the indexes list and the walker weights ______________________
!________________________________________________________________________________________
      do I1 = 1, NWkBt, +1
         Iwalk = WkTabProcs(I1)
         IdptWkIndx(I1) = Iwalk
         WghtProc(Iwalk) = WlkWghtPrc(I1); Log_Wght(Iwalk) = log(WghtProc(I1))
      enddo
#endif
!**************************************************************************************************	  
!___________________ 4. Deallocate all the used matrices __________________________________________
!************************************************************************************************** 
      if(allocated(WkTableAll)) deallocate(WkTableAll)
      if(allocated(WkTabProcs)) deallocate(WkTabProcs)
      if(allocated(WlkWghtAll)) deallocate(WlkWghtAll)
      if(allocated(WlkWghtPrc)) deallocate(WlkWghtPrc)
      if(allocated(AllIdptWkInd)) deallocate(AllIdptWkInd)
      if(allocated(SndRcvAllP)) deallocate(SndRcvAllP)
      if(allocated(SndRcvProc)) deallocate(SndRcvProc)
      if(allocated(SendNmCont)) deallocate(SendNmCont)
      if(allocated(NumScatter)) deallocate(NumScatter)
      if(allocated(DspScatter)) deallocate(DspScatter)
      if(allocated(SendDataBf)) deallocate(SendDataBf)
      if(allocated(RecvDataBf)) deallocate(RecvDataBf)
      if(allocated(SendRequst)) deallocate(SendRequst)
      if(allocated(SendStatus)) deallocate(SendStatus)
      if(allocated(RecvStatus)) deallocate(RecvStatus)
      
   end subroutine PopCtrlBetaT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine ReDistrbWkBetaT(Weight, WkIndTable, WalkWghtTb, NAddWk, NIndpt, NumChg, MaxNum)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ReDistrbWkBetaT(Weight, WkIndTable, WalkWghtTb, NAddWk, NIndpt, NumChg, MaxNum)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to obtain the walker table for all the random walkers, after performing the 
!                    population control process. For NT == LTrot case.
! KEYWORDS: Random Walker table after the population control process at \tau=\beta.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Obtain the random walker table after the population control.    
!
!     Input:  Weight --> The weights of all NmWalkAllP random walkers;
!
!     Output: WkIndTable --> The walker table containing integer indexes; 
!             WalkWghtTb --> The walker weight of all the remaining independent walkers;
!             NAddWk     --> Number of additional walkers to make mod(NIndpt, anprc)==0 on (anprc-1) process;
!             NIndpt     --> Number of independent random walkers after population control process;
!             NumChg     --> Number of random walkers whose information needs to be changed;
!             MaxNum     --> Biggest number of random walkers.
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
      integer NAddWk, NIndpt, NumChg, MaxNum
      integer WkIndTable(NmWalkAllP)
      real(rp) Weight(NmWalkAllP), WalkWghtTb(NmWalkAllP)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I0, I1, I2, I3, I4
      integer StartWkInd, WalkerIndx, ProcsIndex, TempIndx
      integer NmWlkPop, NmWkPop1, NmWkPop2
      integer IdptWkProc(0:anprc-1)
      real(rp) Rtp0, Rtp1
      integer, allocatable :: NumWalkers(:)
      real(rp), allocatable :: NormWeight(:)
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of walker table from population control ________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Allocate arrays used in this subroutine ___________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) The normalized weights of all random walkers _____________________
!________________________________________________________________________________________
      allocate(NormWeight(NmWalkAllP))
      NormWeight = Weight / sum(Weight(1:NmWalkAllP))
!________________________________________________________________________________________ 	  
!_________________ (1) Number of each walker after population control ___________________
!________________________________________________________________________________________
      allocate(NumWalkers(NmWalkAllP))
      NumWalkers = 0
!**************************************************************************************************	  
!___________________ 1. Perform the population control ____________________________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Initialization for the population control ________________________
!________________________________________________________________________________________
      !!!!!!!!!! Determine the number of random walkers using the weights  
      I1 = 1; Rtp0 = NormWeight(1)
      do I0 = 1, NmWalkAllP, +1
         Rtp1 = ( dble(I0-1) + spring_sfmt_stream() ) / dble(NmWalkAllP)
         do while(Rtp1 > Rtp0)
            I1 = I1 + 1
            Rtp0 = Rtp0 + NormWeight(I1)
         enddo
         NumWalkers(I1) = NumWalkers(I1) + 1
      enddo
      !!!!!!!!!! Number of independent random walkers on every process
      IdptWkProc = 0
      do I0 = 1, NmWalkAllP, +1
         ProcsIndex = (I0 - 1) / NWalk
         if(NumWalkers(I0) > 0) IdptWkProc(ProcsIndex) = IdptWkProc(ProcsIndex) + 1
      enddo
      NIndpt = sum(IdptWkProc(0:anprc-1))
      !!!!!!!!!! Compute NWkBt, NmWkBtAllP, and NAddWk
      NmWkBtAllP = merge(NIndpt, (NIndpt/anprc+1)*anprc, mod(NIndpt, anprc)==0) 
      NWkBt = NmWkBtAllP / dble(anprc)
      NAddWk = merge(NmWkBtAllP-NIndpt, 0, mod(NIndpt, anprc)/=0)
      !!!!!!!!!! The number of the walker with largest weight
      MaxNum = maxval(NumWalkers(:), 1)
      !!!!!!!!!! The number of the walkers which needs to be changed
      NumChg = 0
      do ProcsIndex = 0, anprc-1, +1
         NumChg = NumChg + merge(IdptWkProc(ProcsIndex)-NWkBt, 0, IdptWkProc(ProcsIndex) > NWkBt)
      enddo
!________________________________________________________________________________________ 	  
!_________________ (1) Construct the walker tables of indexes and weights _______________
!_____________________ for the NIndpt independent random walkers ________________________
!________________________________________________________________________________________
      WkIndTable = -1001; WalkWghtTb = -1.00_rp
      do ProcsIndex = 0, anprc-1, +1
         StartWkInd = ProcsIndex * NWalk
         I2 = 0
         do I1 = 1, NWalk, +1
            WalkerIndx = StartWkInd + I1
            if(NumWalkers(WalkerIndx) > 0) then
               I2 = I2 + 1
               WkIndTable(StartWkInd+I2) = WalkerIndx
               WalkWghtTb(StartWkInd+I2) = dble(NumWalkers(WalkerIndx))
            end if
         enddo
      enddo
!________________________________________________________________________________________ 	  
!_________________ (2) Adjust walkers on different processes ____________________________
!________________________________________________________________________________________       
      if(anprc > 1) then
         !!!!!!!!!! For [0, NAddWk-1] processes, NWkBt-1 walkers; For [NAddWk, NWkBt] processes, NWkBt walkers; 
         NmWkPop1 = NWkBt - 1; NmWkPop2 = NWkBt
         StartWkInd = 1
         do ProcsIndex = 0, anprc-1, +1
            TempIndx = ProcsIndex * NWalk
            NmWlkPop = merge(NmWkPop1, NmWkPop2, ProcsIndex <= NAddWk-1)
            do I1 = 1, NmWlkPop, +1
               if(WkIndTable(TempIndx+I1) < 0) then
                  do I2 = StartWkInd, NmWalkAllP, +1
                     I3 = (I2-1) / NWalk; I4 = mod(I2-1, NWalk) + 1
                     if( (I3 <= NAddWk-1 .and. I4 > NmWkPop1) .or. (I3 > NAddWk-1 .and. I4 > NmWkPop2) ) then
                        if(WkIndTable(I2) > 0) then
                           WkIndTable(TempIndx+I1) = WkIndTable(I2)
                           WalkWghtTb(TempIndx+I1) = WalkWghtTb(I2)
                           WkIndTable(I2) = -1001; WalkWghtTb(I2) = -1.00_rp
                           exit
                        end if
                     end if
                  enddo
                  StartWkInd = I2 + 1
               end if
            enddo
         enddo
         !!!!!!!!!! Put the NAddWk additional walkers on [0, NAddWk-1] processes, each with one walker
         if(NAddWk > 0) then
            do ProcsIndex = 0, NAddWk-1, +1
               I0 = ProcsIndex * NWalk
               WalkWghtTb(I0+NWkBt-1) = WalkWghtTb(I0+NWkBt-1) / dble(1+1)
               WkIndTable(I0+NWkBt) = WkIndTable(I0+NWkBt-1)
               WalkWghtTb(I0+NWkBt) = WalkWghtTb(I0+NWkBt-1)
            enddo
         end if 
      end if
!________________________________________________________________________________________ 	  
!_________________ (3) Perform checks for WkIndTable and WalkWghtTb _____________________
!________________________________________________________________________________________ 
      do ProcsIndex = 0, anprc-1, +1
         TempIndx = ProcsIndex * NWalk
         do I1 = 1, NWkBt, +1
            if(WkIndTable(TempIndx+I1) < 0 .or. WalkWghtTb(TempIndx+I1) <= 0.0_rp) then
               write(*, "(A)") "=================================================================================="
               write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
               write(*, "(3x, 'ReDistrbWkBetaT: I1 <= NWkBt!!! amyid, I1, WkIndTable, WalkWghtTb = ', 3I6, es17.8)") &
                  & ProcsIndex, I1, WkIndTable(TempIndx+I1), WalkWghtTb(TempIndx+I1)
               write(*, "(A)") "=================================================================================="
            end if
         enddo
         do I1 = NWkBt+1, NWalk, +1
            if(WkIndTable(TempIndx+I1) > 0 .or. WalkWghtTb(TempIndx+I1) > 0.0_rp) then
               write(*, "(A)") "=================================================================================="
               write(*, "(A)") "============= Warning!!! Warning!!! Warning!!! Warning!!! Warning!!! ============="
               write(*, "(3x, 'ReDistrbWkBetaT: I1 > NWkBt!!! amyid, I1, WkIndTable, WalkWghtTb = ', 3I6, es17.8)") &
                  & ProcsIndex, I1, WkIndTable(TempIndx+I1), WalkWghtTb(TempIndx+I1)
               write(*, "(A)") "=================================================================================="
            end if
         enddo
      enddo
!**************************************************************************************************	  
!___________________ 2. Deallocate the matrices used in this subroutine ___________________________
!**************************************************************************************************
      if(allocated(NumWalkers)) deallocate(NumWalkers)
      if(allocated(NormWeight)) deallocate(NormWeight)
      
   end subroutine ReDistrbWkBetaT
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$