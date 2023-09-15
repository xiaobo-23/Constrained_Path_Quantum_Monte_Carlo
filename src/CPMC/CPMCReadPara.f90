!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSODIAL SPIN PINNING FIELDS; USING PERIODIC BOUNDARY CONDITION (PBC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine CPMCReadPara()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  CPMCReadPara()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to read the input calculation parameters for CPMC calculations.
! KEYWORDS: Read Input model parameters.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Read the simulation parameters from the input file. 
!  
!     If the input file doesn't exit, we simply turn to initialize the simulation parameters by default.
!
!     Input: (none);  Outpt: (none).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use RealPrecsn
      use CoreParamt
      use Observable
      use StdInOutSt
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IndexRead
      character(200) LengthOfPath
!______________________________________________________________________________________________________________     
!______________________________ Initialization of all parameters for CPMC simulation __________________________
!______________________________________________________________________________________________________________
!================================================================================================== 
!___________________ 0. Some preprocess conditions for the code ___________________________________
!==================================================================================================
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (0) Whether to output the QR check ___________________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      IfCheckUDV = .false.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (1) Whether to Compare the Green's function __________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
      IfNmStbErr = .true.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!_________________ (2) Whether to output Max and Min of weights _________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ! IfWghtMxMn = .false.
      IfWghtMxMn = .true.
!================================================================================================== 
!___________________ 1. Read input parameters or give values by default ___________________________
!==================================================================================================
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (0) Read the input parameters from txt file __________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      open(57, err = 58, file = "Input/ModelParam.txt", status = "old")
      read(57, *) NumL1, NumL2                            ! NumL1, NumL2
      read(57, *) IfPr1, IfPr2                            ! IfPr1, IfPr2
      read(57, *) Hopt1Up, Hopt1Dn, EkDispType            ! Hopt1Up, Hopt1Dn, EkDispType
      read(57, *) Hopt2, Hopt3                            ! Hopt2, Hopt3
      read(57, *) ZmFdz                                   ! ZmFdz
      read(57, *) PinSz, PinSzType                        ! PinSz, PinSzType
      read(57, *) HubbU, ChemP                            ! HubbU, ChemP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      read(57, *) ifSinusoidalPinning                     ! ifSinusoidalPinning
      read(57, *) SinusoidalPinSz, LambdaSz               ! SinusoidalPinSz, LambdaSz
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      read(57, *) SymTrotDcp, HS_Type                     ! SymTrotDcp, HS_Type
      read(57, *) FFTEXPDTH0, FFTEXPDTHT, IfFftEnPar      ! FFTEXPDTH0, FFTEXPDTHT, IfFftEnPar
      read(57, *) MuBTType, ChemP_BT, NItrGetMuT          ! MuBTType, ChemP_BT, NItrGetMuT
      read(57, *) UpdtMethod, NblkUDelay                  ! UpdtMethod, NblkUDelay
      read(57, *) NWalk                                   ! NWalk
      read(57, *) PopWghtMax, PopWghtMin                  ! PopWghtMax, PopWghtMin
      read(57, *) IfSetDltau, FixedDltau                  ! IfSetDltau, FixedDltau
      read(57, *) BetaT, LTrot                            ! BetaT, LTrot
      read(57, *) IfFixnT, NSwepFix                       ! IfFixnT, NSwepFix
      read(57, *) FixNenT, NumNe, Fix_nT                  ! FixNenT, NumNe, Fix_nT
      read(57, *) ChPSmallnT, ChPLargenT                  ! ChPSmallnT, ChPLargenT
      read(57, *) IfSwepReWt                              ! IfSwepReWt
      read(57, *) NmBin, NSwep, nWarm                     ! NmBin, NSwep, nWarm
      read(57, *) NvStb, NvStbOut                         ! NvStb, NvStbOut   
      read(57, *) NvPop, NvPopOut                         ! NvPop, NvPopOut
      read(57, *) IfM2OneMea, NvStbM2One, NvMeaM2One      ! IfM2OneMea, NvStbM2One, NvMeaM2One
      read(57, *) IfCutSmLg, NmBinCut                     ! IfCutSmLg, NmBinCut
      read(57, *) IfTAU, IfTau0Rand, NmTDMType, NmTDM     ! IfTAU, IfTau0Rand, NmTDMType, NmTDM
      read(57, *) IfEqDistDt, NvMeaDyn                    ! IfEqDistDt, NvMeaDyn
      read(57, *) IfDyGrFOut, IfDySpnOut, IfDyDenOut      ! IfDyGrFOut, IfDySpnOut, IfDyDenOut
      read(57, *) IfDyPstOut, IfDyDWvOut, IfDyCurrnt      ! IfDyPstOut, IfDyDWvOut, IfDyCurrnt
      read(57, *) IfSaveFlds, SaveFldMea, ReadFldMea      ! IfSaveFlds, SaveFldMea, ReadFldMea
      read(57, "(A500)") MeaFldPath                       ! MeaFldPath
      close(57) 

      IndexRead = 1
      
      go to 59
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (1) Initialize parameters by default _________________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
58    continue
   
      NumL1 = 4; NumL2 = 4
      IfPr1 = .true.; IfPr2 = .true.
      Hopt1Up = - 1.000_rp; Hopt1Dn = - 0.800_rp; EkDispType = 0
      Hopt2 = + 0.200_rp; Hopt3 = + 0.240_rp
      ZmFdz = + 0.000_rp
      PinSz = + 0.10000_rp; PinSzType = 1
      HubbU = + 2.000_rp; ChemP = - 0.500_rp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ifSinusoidalPinning = .false.; 
      SinusoidalPinSz = 0.0_rp; LambdaSz = 0.0_rp
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      SymTrotDcp = .false.; HS_Type = 1
      FFTEXPDTH0 = .false.; FFTEXPDTHT = .false.; IfFftEnPar = .false.
      MuBTType = 1; ChemP_BT = + 1.000_rp; NItrGetMuT = 2
      UpdtMethod = 1; NblkUDelay = NumL1*NumL2/2
      NWalk = 100
      PopWghtMax = + 5.000_rp; PopWghtMin = + 0.200_rp
      IfSetDltau = .true.; FixedDltau = + 0.02_rp
      BetaT = 12.000_rp; LTrot = 1200
      IfFixnT = .false.; NSwepFix = 100
      FixNenT = 0; NumNe = 058.000_rp; Fix_nT = 0.875_rp
      ChPSmallnT = + 3.000_rp; ChPLargenT = + 0.500_rp
      IfSwepReWt = .true.
      NmBin = 020; NSwep = 020; nWarm = 300
      NvStb = 020; NvStbOut = 020
      NvPop = 010; NvPopOut = 020
      IfM2OneMea = .false.; NvStbM2One = 10; NvMeaM2One = 02
      IfCutSmLg = .false.; NmBinCut = 0
      IfTAU = .false.; IfTau0Rand = .false.; NmTDMType = 0; NmTDM = LTrot/2
      IfEqDistDt = .true. ; NvMeaDyn = 10
      IfDyGrFOut = .false.; IfDySpnOut = .false.; IfDyDenOut = .false.
      IfDyPstOut = .false.; IfDyDWvOut = .false.; IfDyCurrnt = .false.
      IfSaveFlds = .false.; SaveFldMea = .false.; ReadFldMea = .false.
      MeaFldPath = "Output"

      IndexRead = 0
      
59    continue
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (2) Reprocess the input model and simulation parameters ______________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!______________________________________________________________________
!________________ For the logic parameter IfPyObsPBC __________________
!______________________________________________________________________
      if( (IfPr1 .and. IfPr2) .and. (abs(PinSz) < rp_Eps) ) then
      ! if( IfPr1 .and. IfPr2 ) then
         IfPyObsPBC = .true.
      else
         IfPyObsPBC = .false.
      end if
!______________________________________________________________________
!________________ For the logic parameter IfCrfDfBkg __________________
!______________________________________________________________________
      IfCrfDfBkg = .false.
!______________________________________________________________________
!________________ For the NumNC, NumNS and NmSpn parameters ___________
!______________________________________________________________________
      NumNC = NumL1 * NumL2
      NumNS = NumNC
      NmSpn = 2
!______________________________________________________________________
!________________ For the NmWalkAllP __________________________________
!______________________________________________________________________
      NmWalkAllP = anprc * NWalk
!______________________________________________________________________
!________________ For the BetaT, TempT and Dltau ______________________
!________________ Keep LTrot to be an even number _____________________
!______________________________________________________________________  
      TempT = 1.0_rp / BetaT
      if(IfSetDltau) LTrot = ceiling(BetaT/FixedDltau)
      if(mod(LTrot, 2) == 1) LTrot = LTrot + 1
      Dltau = BetaT / dble(LTrot)
!______________________________________________________________________
!________________ For FFTEXPDTH0, FFTEXPDTHT, IfFftEnPar parameters ___
!______________________________________________________________________
      if( (.not. IfPr1) .or. (.not. IfPr2) ) then
         FFTEXPDTH0 = .false.; FFTEXPDTHT = .false.
      end if
      if(abs(PinSz) >= rp_Eps) FFTEXPDTHT = .false.
      if(.not. FFTEXPDTH0) IfFftEnPar = .false.
!______________________________________________________________________
!________________ For NumNe and Fix_nT ________________________________
!______________________________________________________________________
      if(FixNenT == 0) then
         Fix_nT = NumNe / dble(NumNS)
      else
         NumNe = Fix_nT * dble(NumNS)
      end if
!______________________________________________________________________
!________________ For MuBTType and IfMuTqmcNt _________________________
!______________________________________________________________________
      IfMuTqmcNt = merge(.true., .false., MuBTType == 4)
!______________________________________________________________________
!________________ For IfSwepReWt and FrqReCmptGrowth __________________
!______________________________________________________________________
      FrqReCmptGrowth = merge(NSwep, 1, IfSwepReWt)
!______________________________________________________________________
!________________ For the NvMeaM2One parameter ________________________
!______________________________________________________________________
      if(IfM2OneMea) then
         if( NvMeaM2One >= LTrot ) then
            NvMeaM2One = LTrot
         else
            do while( mod(LTrot, NvMeaM2One) /= 0 )
               NvMeaM2One = NvMeaM2One + 1
            enddo
         end if
      end if
!______________________________________________________________________
!________________ For the IfTAU and NmTDM parameters __________________
!______________________________________________________________________
      if(.not. IfM2OneMea) IfTAU = .false.
      if(IfTAU) then
         if(NmTDMType == 0) then
            NmTDM = LTrot/2
         else if(NmTDMType == 1) then
            NmTDM = LTrot
         else if(NmTDM > LTrot) then
            NmTDM = LTrot
         end if
      end if
!______________________________________________________________________
!________________ For the nWarm and ReadFldMea parameter ______________
!______________________________________________________________________ 
      if( abs(HubbU) <= rp_Eps ) then
         IfSaveFlds = .false.
         SaveFldMea = .false.
         ReadFldMea = .false.
      end if
      if(ReadFldMea) then
         nWarm = 0
         IfSaveFlds = .false.
      end if
!______________________________________________________________________
!________________ For the NmBinCut, BinStart and Bin_Data _____________
!______________________________________________________________________        
      BinStart = NmBinCut + 1
      Bin_Data = NmBin - NmBinCut
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCReadPara process _________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(   
      if(amyid == amstr) then
         write(*, "(A)") "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         if(IndexRead == 1) then
            write(*, "(2x, 'CPMCReadPara: Parameters are read from input as')")
         else if(IndexRead == 0) then
            write(*, "(2x, 'CPMCReadPara: Parameters are set by default as')")
         end if
         write(*, "()")
         write(*, "(16x, '                 NumL1, NumL2, NumNS = ', 3I6       )") NumL1, NumL2, NumNS
         write(*, "(16x, '                        IfPr1, IfPr2 = ', 2l6       )") IfPr1, IfPr2
         write(*, "(16x, '        Hopt1Up, Hopt1Dn, EkDispType = ',  2f7.3, I6)") Hopt1Up, Hopt1Dn, EkDispType
         write(*, "(16x, '                        Hopt2, Hopt3 = ',      2f7.3)") Hopt2, Hopt3
         write(*, "(16x, '                               ZmFdz = ',       f7.3)") ZmFdz
         write(*, "(16x, '                    PinSz, PinSzType = ',  f9.5,  I6)") PinSz, PinSzType
         write(*, "(16x, '                        HubbU, ChemP = ', sp,2f14.10)") HubbU, ChemP
         write(*, "(16x, '          SymTrotDcp, HS_Type, NmSpn = ', l6,    2I6)") SymTrotDcp, HS_Type, NmSpn
         write(*, "(16x, '  FFTEXPDTH0, FFTEXPDTHT, IfFftEnPar = ', 3l6       )") FFTEXPDTH0, FFTEXPDTHT, IfFftEnPar
         write(*, "(16x, '      MuBTType, ChemP_BT, NItrGetMuT = ', I6, spf9.5)") MuBTType, ChemP_BT, NItrGetMuT
         write(*, "(16x, '              UpdtMethod, NblkUDelay = ', 2I6       )") UpdtMethod, NblkUDelay
         write(*, "(16x, '                   NWalk, NmWalkAllP = ', 2I6       )") NWalk, NmWalkAllP
         write(*, "(16x, '              PopWghtMax, PopWghtMin = ',      2f7.3)") PopWghtMax, PopWghtMin
         write(*, "(16x, '  IfCheckUDV, IfNmStbErr, IfWghtMxMn = ', 3l6       )") IfCheckUDV, IfNmStbErr, IfWghtMxMn
         write(*, "(16x, '              IfSetDltau, FixedDltau = ',  l6,  f7.3)") IfSetDltau, FixedDltau
         write(*, "(16x, '                        BetaT, TempT = ',    2f14.10)") BetaT, TempT
         write(*, "(16x, '                        Dltau, LTrot = ', f13.10, I6)") Dltau, LTrot
         write(*, "(16x, '              IfPyObsPBC, IfCrfDfBkg = ', 2l6       )") IfPyObsPBC, IfCrfDfBkg
         write(*, "(16x, '                   IfFixnT, NSwepFix = ',  l6,    I6)") IfFixnT, NSwepFix
         write(*, "(16x, '              FixNenT, NumNe, Fix_nT = ',  I6,2f12.6)") FixNenT, NumNe, Fix_nT
         write(*, "(16x, '              ChPSmallnT, ChPLargenT = ', sp, 2f12.8)") ChPSmallnT, ChPLargenT
         write(*, "(16x, '         IfSwepReWt, FrqReCmptGrowth = ',  l6,    I6)") IfSwepReWt, FrqReCmptGrowth
         write(*, "(16x, '                 NmBin, NSwep, nWarm = ', 3I6       )") NmBin, NSwep, nWarm
         write(*, "(16x, '                     NvStb, NvStbOut = ', 2I6       )") NvStb, NvStbOut
         write(*, "(16x, '                     NvPop, NvPopOut = ', 2I6       )") NvPop, NvPopOut
         write(*, "(16x, '  IfM2OneMea, NvStbM2One, NvMeaM2One = ',  l6,   2I6)") IfM2OneMea, NvStbM2One, NvMeaM2One
         write(*, "(16x, '       IfCutSmLg, NmBinCut, Bin_Data = ',  l6,   2I6)") IfCutSmLg, NmBinCut, Bin_Data
         write(*, "(16x, ' IfTAU, IfTau0Rand, NmTDMType, NmTDM = ', 2l6,   2I6)") IfTAU, IfTau0Rand, NmTDMType, NmTDM
         write(*, "(16x, '                IfEqDistDt, NvMeaDyn = ',  l6,    I6)") IfEqDistDt, NvMeaDyn
         write(*, "(16x, '  IfDyGrFOut, IfDySpnOut, IfDyDenOut = ', 3l6       )") IfDyGrFOut, IfDySpnOut, IfDyDenOut
         write(*, "(16x, '  IfDyPstOut, IfDyDWvOut, IfDyCurrnt = ', 3l6       )") IfDyPstOut, IfDyDWvOut, IfDyCurrnt
         write(*, "(16x, '  IfSaveFlds, SaveFldMea, ReadFldMea = ', 3l6       )") IfSaveFlds, SaveFldMea, ReadFldMea
         if(ReadFldMea) then
            write(LengthOfPath, *) len_Trim(MeaFldPath)
            write(*, "(16x, '                           MeaFldPath = ',  A" // Trim(LengthOfPath) // ")") &
               & Trim(MeaFldPath)
         end if
      end if
!================================================================================================== 
!___________________ 2. Some post process of the input data _______________________________________
!==================================================================================================
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (0) Check input of PinSz and PinSzType values ________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      if( (abs(PinSz) >= rp_Eps) .and. (PinSzType /= 0) .and. (PinSzType /= 1) .and. (PinSzType /= 2) ) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'ERROR: Wrong input of PinSzType parameter! PinSzType = ', I3)") PinSzType
            write(*, "(16x, '       Will stop running the program!!! Please check the input file!!!')")
         end if
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
         stop
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (1) Check input HS_Type value for validity ___________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      if( (HS_Type /= 0) .and. (HS_Type /= 1) ) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'WARNING: Invalid input value of HS_Type, HS_Type = ', I4)") HS_Type
            write(*, "(16x, '         Will stop running the program!!! Please check the input file!!!')")
         end if
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
         stop
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (2) Check input HS_Type and sign of HubbU ____________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      if(     ( (HubbU < -rp_Eps) .and. (HS_Type == 1) ) .or. &
            & ( (HubbU > +rp_Eps) .and. (HS_Type == 0) ) ) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'WARNING: Inconsistent HubbU and HS_Type!!! HubbU, HS_Type = ', f12.8, I4)") HubbU, HS_Type
            write(*, "(16x, '         Will stop running the program!!! Please check the input file!!!')")
         end if
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
         stop
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (3) Check whether IfPr1, IfPr2 and PinSz is compatible _______________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      if( (IfPr1 .and. IfPr2) .and. (abs(PinSz) >= rp_Eps .and. PinSzType /= 2) ) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'WARNING: The IfPr1, IfPr2 tags are not compatible with nonzero PinSz!')")
            write(*, "(16x, '         Will stop running the code!!! Please check the input file!!!')")
         end if
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
         stop
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (4) Check input UpdtMethod value for validity ________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&       
      if( (UpdtMethod /= 0) .and. (UpdtMethod /= 1)  .and. (UpdtMethod /= 2) ) then
         if(amyid == amstr) then  
            write(*, "()")
            write(*, "(16x, 'WARNING: Illegal input for UpdtMethod!')")
            write(*, "(16x, '         Will stop running the code!!! Please check the input file!!!')")
         end if
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
         stop
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (5) Check EkDispType, FFTEXPDTH0, FFTEXPDTHT for consistence _________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if( ((.not. FFTEXPDTH0) .or. (.not. FFTEXPDTHT)) .and. (EkDispType == 2) ) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'WARNING: The IfFftEnPar tag is not compatible with EkDispType==2!')")
            write(*, "(16x, '         Will stop running the code!!! Please check the input file!!!')")
         end if
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
         stop         
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (6) Check whether IfMuTqmcNt is compatible with IfFixnT ______________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(IfMuTqmcNt .and. IfFixnT) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'WARNING: IfFixnT == T is incompatible with IfMuTqmcNt == T!!!')") 
            write(*, "(16x, '         Will stop running the code!!! Please check the input file!!!')")
         end if
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (7) Check whether to fix the total electron density __________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(IfFixnT) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'Will adjust ChemP to reach desired n_Occ as input IfFixnT = ', l3)") IfFixnT
         end if
         ChemP_Ref(1) = ChPSmallnT    ! Small density
         ChemP_Ref(2) = ChPLargenT    ! Large density
      end if
      Fix_Iterate = 0
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (8) Check input NmBin value for validity _____________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
      if( (IfCutSmLg) .and. (NmBin < 3+NmBinCut) ) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'WARNING: Input NmBin is too small! NmBin = ', I4)") NmBin
            write(*, "(16x, '         Will increase the value of NmBin to 4+NmBinCut = ', I4)") 4+NmBinCut
         end if 
         NmBin = 4+NmBinCut
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (9) Check input IfFixnT, SaveFldMea and ReadFldMea ___________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      !!!!!!!!!! Check IfFixnT and ReadFldMea
      if(IfFixnT .and. ReadFldMea) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'WARNING: IfFixnT == T and ReadFldMea == T can not happen togather!!!')")
            write(*, "(16x, '         Will stop running the code!!! Please check the input file!!!')")
         end if
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
         stop
      end if
      !!!!!!!!!! Check SaveFldMea and ReadFldMea
      if(SaveFldMea .and. ReadFldMea) then
         if(amyid == amstr) then
            write(*, "()")
            write(*, "(16x, 'WARNING: SaveFldMea == T and ReadFldMea == T can not happen togather!!!')")
            write(*, "(16x, '         Will stop running the code!!! Please check the input file!!!')")
         end if
#ifdef MPIPROCESS
         call MPI_Barrier(acomm, ierr)    
#endif
         stop
      end if
!================================================================================================== 
!___________________ 3. Create folders for storing output results _________________________________
!==================================================================================================
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr) 
#endif
      if(amyid == amstr) then
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (0) The Folder for random number seeds _______________________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         call system("mkdir -p Output/RandNumbSeeds")
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (1) Folders for auxiliary fields of configurations ___________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         if(IfSaveFlds) call system("mkdir -p Output/AuxiliaryFlds")
         if(SaveFldMea) call system("mkdir -p Output/AuxiliaryFlds_Measure")
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (2) The Folder for M2One measure for IfM2OneMea == T _________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         if(IfM2OneMea) call system("mkdir -p Add_Output")
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  
!_________________ (3) Folders for dynamical correlation functions ______________________
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         if(IfTAU .and. IfPyObsPBC) then
            call system("mkdir -p Add_Output/TT_GrFCrFIwnDyData")
            if(IfDyGrFOut) call system("mkdir -p Add_Output/UU_GreenFDynmcData")
            if(IfDySpnOut) call system("mkdir -p Add_Output/VV_SpnSpnDynmcData")
            if(IfDyDenOut) call system("mkdir -p Add_Output/WW_DenDenDynmcData")
            if(IfDyPstOut) call system("mkdir -p Add_Output/XX_PairStDynmcData")
            if(IfDyDWvOut) call system("mkdir -p Add_Output/YY_DWvParDynmcData")
            if(IfDyCurrnt) call system("mkdir -p Add_Output/ZZ_CurrntDynmcData")
         end if
      end if
#ifdef MPIPROCESS
      call MPI_Barrier(acomm, ierr) 
#endif
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of CPMCReadPara process _________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "(A)") "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
         write(*, "( )")
         write(*, "( )")
      end if

   end subroutine CPMCReadPara
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$