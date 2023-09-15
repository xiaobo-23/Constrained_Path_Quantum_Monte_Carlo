!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform initally listing og all constants in the program  before the whole 
!                    CPMC simulations.
! COMMENT: CPMC Initialization process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   InitListCons --> Subroutine to list all constants defined in CPMC Simulations;
!   ListItem     --> Subroutine to output a single constant by some format.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine InitListCons()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitListCons()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the output of all constants defined in modules during the 
!                   CPMC simulations.
! KEYWORDS: List all CPMC calculation constants.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Initialization of constants in CPMC, including:
!                  (0) Constants defined in TimeRecord module;
!                  (1) Constants defined in StdInOutSt module;
!                  (2) Constants defined in RealPrecsn module;
!                  (3) Constants defined in CoreParamt module;
!                  (4) Constants defined in Observable module.
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
            use TimeRecord
            use StdInOutSt
            use RealPrecsn
            use CoreParamt
            use Observable
            use MPISetting
            implicit none
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
            character(512) ChaRep
            character(100) FileName
!______________________________________________________________________________________________________________     
!_______________________________ Main calculations of Listing constants _______________________________________
!______________________________________________________________________________________________________________
      if(amyid == amstr) then
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of initialization process _______________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
            write(*, "()")
            write(*, "(16x, 'InitListCons: List all the important, used quantities!')")
!__________________________________________________________________________________________________     
!________________________________ Open the output file ____________________________________________
!__________________________________________________________________________________________________
            open(112, file = Trim(FLog), access = "append")
            write(FileName, *) DatTimLen
!**************************************************************************************************    
!___________________ 0. List Constants in MPISetting module _______________________________________
!**************************************************************************************************
            write(112, "(" // adjustl(FileName) // "x, 4x, A)") "Module MPISetting"
            write(ChaRep,           *) acomm;              call ListItem(112, ChaRep, "acomm")
            write(ChaRep,           *) anprc;              call ListItem(112, ChaRep, "anprc")
            write(ChaRep,           *) amyid;              call ListItem(112, ChaRep, "amyid")
            write(ChaRep,           *) amstr;              call ListItem(112, ChaRep, "amstr")
            write(ChaRep,           *) OneMegaNum;         call ListItem(112, ChaRep, "OneMegaNum")
!**************************************************************************************************     
!___________________ 1. List Constants in TimeRecord module _______________________________________
!**************************************************************************************************
            write(112, "(" // adjustl(FileName) // "x, 4x, A)") "Module TimeRecord"
            write(ChaRep,           *) TimRat;             call ListItem(112, ChaRep, "TimRat")
            write(ChaRep,           *) TimMax;             call ListItem(112, ChaRep, "TimMax")
            write(ChaRep,           *) TimThh;             call ListItem(112, ChaRep, "TimThh")
            write(ChaRep,           *) DatTimLen;          call ListItem(112, ChaRep, "DatTimLen")
!**************************************************************************************************     
!___________________ 2. List Constants in StdInOutSt module _______________________________________
!************************************************************************************************** 
            write(112, "(" // adjustl(FileName) // "x, 4x, A)") "Module StdInOutSt"
            write(ChaRep,           *) FLog;               call ListItem(112, ChaRep, "FLog")
            write(ChaRep,           *) FWrn;               call ListItem(112, ChaRep, "FWrn")
            write(ChaRep,           *) FErr;               call ListItem(112, ChaRep, "FErr")
            write(ChaRep,           *) FMnt;               call ListItem(112, ChaRep, "FMnt")
!**************************************************************************************************     
!___________________ 3. List Constants in RealPrecsn module _______________________________________
!**************************************************************************************************         
            write(112, "(" // adjustl(FileName) // "x, 4x, A)") "Module RealPrecsn"
            write(ChaRep,           *) rp;                 call ListItem(112, ChaRep, "rp")
            write(ChaRep,           *) rp_ow;              call ListItem(112, ChaRep, "rp_ow")
            write(ChaRep,           *) rp_od;              call ListItem(112, ChaRep, "rp_od")
            write(ChaRep, "(es25.16)") rp_Rone;            call ListItem(112, ChaRep, "rp_Rone")
            write(ChaRep, "(es25.16)") rp_pi;              call ListItem(112, ChaRep, "rp_pi")
            write(ChaRep, "(es25.16)") rp_gold;            call ListItem(112, ChaRep, "rp_gold")
            write(ChaRep,           *) rp_huge;            call ListItem(112, ChaRep, "rp_huge")
            write(ChaRep,           *) rp_tiny;            call ListItem(112, ChaRep, "rp_tiny")
            write(ChaRep, "(es25.16)") rp_prec;            call ListItem(112, ChaRep, "rp_prec")
            write(ChaRep, "(es25.16)") rp_Eps;             call ListItem(112, ChaRep, "rp_Eps")
!**************************************************************************************************     
!___________________ 4. List Constants in CoreParamt module _______________________________________
!**************************************************************************************************
            write(112, "(" // adjustl(FileName) // "x, 4x, A)") "Module CoreParamt"
            write(ChaRep,           *) NumL1;              call ListItem(112, ChaRep, "NumL1")
            write(ChaRep,           *) NumL2;              call ListItem(112, ChaRep, "NumL2")
            write(ChaRep,           *) NumNC;              call ListItem(112, ChaRep, "NumNC")
            write(ChaRep,           *) NumNS;              call ListItem(112, ChaRep, "NumNS")
            write(ChaRep,           *) NmSpn;              call ListItem(112, ChaRep, "NmSpn")
         
            write(ChaRep, "(    l12)") IfPr1;              call ListItem(112, ChaRep, "IfPr1")
            write(ChaRep, "(    l12)") IfPr2;              call ListItem(112, ChaRep, "IfPr2")
         
            write(ChaRep, "(2es25.16)") cmplx(A1Vec(1), A1Vec(2), rp);        call ListItem(112, ChaRep, "A1Vec")
            write(ChaRep, "(2es25.16)") cmplx(A2Vec(1), A2Vec(2), rp);        call ListItem(112, ChaRep, "A2Vec")
            write(ChaRep, "(2es25.16)") cmplx(B1Vec(1), B1Vec(2), rp);        call ListItem(112, ChaRep, "B1Vec")
            write(ChaRep, "(2es25.16)") cmplx(B2Vec(1), B2Vec(2), rp);        call ListItem(112, ChaRep, "B2Vec")
         
            write(ChaRep, "(es23.14)") Hopt1Up;            call ListItem(112, ChaRep, "Hopt1Up")
            write(ChaRep, "(es23.14)") Hopt1Dn;            call ListItem(112, ChaRep, "Hopt1Dn")
            write(ChaRep,           *) EkDispType;         call ListItem(112, ChaRep, "EkDispType")
            write(ChaRep, "(es23.14)") Hopt2;              call ListItem(112, ChaRep, "Hopt2")
            write(ChaRep, "(es23.14)") Hopt3;              call ListItem(112, ChaRep, "Hopt3")
            write(ChaRep, "(es23.14)") ZmFdz;              call ListItem(112, ChaRep, "ZmFdz")
            write(ChaRep, "(es23.14)") PinSz;              call ListItem(112, ChaRep, "PinSz")
            write(ChaRep,           *) PinSzType;          call ListItem(112, ChaRep, "PinSzType")
            write(ChaRep, "(es23.14)") ChemP;              call ListItem(112, ChaRep, "ChemP")
            write(ChaRep, "(es23.14)") HubbU;              call ListItem(112, ChaRep, "HubbU")
         
            write(ChaRep,           *) MuBTType;           call ListItem(112, ChaRep, "MuBTType")
            write(ChaRep, "(es25.16)") ChemP_BT;           call ListItem(112, ChaRep, "ChemP_BT")

            write(ChaRep, "(    l12)") IfCheckUDV;         call ListItem(112, ChaRep, "IfCheckUDV")
            write(ChaRep, "(    l12)") IfNmStbErr;         call ListItem(112, ChaRep, "IfNmStbErr")
            
            write(ChaRep, "(    l12)") SymTrotDcp;         call ListItem(112, ChaRep, "SymTrotDcp")
            write(ChaRep, "(    l12)") FFTEXPDTH0;         call ListItem(112, ChaRep, "FFTEXPDTH0")
            write(ChaRep, "(    l12)") FFTEXPDTHT;         call ListItem(112, ChaRep, "FFTEXPDTHT")

            write(ChaRep,           *) HS_Type;            call ListItem(112, ChaRep, "HS_Type")
         
            write(ChaRep, "(es25.16)") LmdbU;              call ListItem(112, ChaRep, "LmdbU")
         
            write(ChaRep, "(    l12)") IfSetDltau;         call ListItem(112, ChaRep, "IfSetDltau")
            write(ChaRep, "(es25.16)") FixedDltau;         call ListItem(112, ChaRep, "FixedDltau")
         
            write(ChaRep, "(es25.16)") BetaT;              call ListItem(112, ChaRep, "BetaT")
            write(ChaRep, "(es25.16)") TempT;              call ListItem(112, ChaRep, "TempT")
            write(ChaRep, "(es25.16)") Dltau;              call ListItem(112, ChaRep, "Dltau")
            write(ChaRep,           *) LTrot;              call ListItem(112, ChaRep, "LTrot")
            write(ChaRep,           *) nWarm;              call ListItem(112, ChaRep, "nWarm")
         
            write(ChaRep, "(    l12)") IfFixnT;            call ListItem(112, ChaRep, "IfFixnT")
            write(ChaRep,           *) NSwepFix;           call ListItem(112, ChaRep, "NSwepFix")

            write(ChaRep,           *) FixNenT;            call ListItem(112, ChaRep, "FixNenT")
            write(ChaRep, "(es25.16)") Fix_nT;             call ListItem(112, ChaRep, "Fix_nT")
            write(ChaRep, "(es25.16)") NumNe;              call ListItem(112, ChaRep, "NumNe")
            write(ChaRep, "(es25.16)") ChPSmallnT;         call ListItem(112, ChaRep, "ChPSmallnT")
            write(ChaRep, "(es25.16)") ChPLargenT;         call ListItem(112, ChaRep, "ChPLargenT")
         
            write(ChaRep,           *) NmBIN;              call ListItem(112, ChaRep, "NmBIN")
            write(ChaRep,           *) NSwep;              call ListItem(112, ChaRep, "NSwep")
         
            write(ChaRep,           *) NvStb;              call ListItem(112, ChaRep, "NvStb")
            write(ChaRep,           *) NvStbOut;           call ListItem(112, ChaRep, "NvStbOut")
         
            write(ChaRep,           *) NvPop;              call ListItem(112, ChaRep, "NvPop")
            write(ChaRep,           *) NvPopOut;           call ListItem(112, ChaRep, "NvPopOut")
         
            write(ChaRep, "(es25.16)") PopWghtMax;         call ListItem(112, ChaRep, "PopWghtMax")
            write(ChaRep, "(es25.16)") PopWghtMin;         call ListItem(112, ChaRep, "PopWghtMin")
         
            write(ChaRep, "(    l12)") IfM2OneMea;         call ListItem(112, ChaRep, "IfM2OneMea")
            write(ChaRep,           *) NvStbM2One;         call ListItem(112, ChaRep, "NvStbM2One")
            write(ChaRep,           *) NvMeaM2One;         call ListItem(112, ChaRep, "NvMeaM2One")
!__________________________________________________________________________________________________     
!___________________ 5. List Constants in Observable module _______________________________________
!__________________________________________________________________________________________________
            write(112, "(" // adjustl(FileName) // "x, 4x, A)") "Module Observable"
            write(ChaRep, "(    l12)") IfTAU;              call ListItem(112, ChaRep, "IfTAU")
            write(ChaRep,           *) NmTDM;              call ListItem(112, ChaRep, "NmTDM")

            write(ChaRep,           *) NumTauPnt;          call ListItem(112, ChaRep, "NumTauPnt")
            write(ChaRep,           *) TauHalfBetaT;       call ListItem(112, ChaRep, "TauHalfBetaT")
         
            write(ChaRep, "(    l12)") IfEqDistDt;         call ListItem(112, ChaRep, "IfEqDistDt")
            write(ChaRep,           *) NvMeaDyn;           call ListItem(112, ChaRep, "NvMeaDyn")
         
            write(ChaRep, "(    l12)") IfDyGrFOut;         call ListItem(112, ChaRep, "IfDyGrFOut")
            write(ChaRep, "(    l12)") IfDySpnOut;         call ListItem(112, ChaRep, "IfDySpnOut")
            write(ChaRep, "(    l12)") IfDyDenOut;         call ListItem(112, ChaRep, "IfDyDenOut")
            write(ChaRep, "(    l12)") IfDyPstOut;         call ListItem(112, ChaRep, "IfDyPstOut")
            write(ChaRep, "(    l12)") IfDyDWvOut;         call ListItem(112, ChaRep, "IfDyDWvOut")
            write(ChaRep, "(    l12)") IfDyCurrnt;         call ListItem(112, ChaRep, "IfDyCurrnt")
!__________________________________________________________________________________________________     
!________________________________ close the output file ___________________________________________
!__________________________________________________________________________________________________
            close(112)
      end if

      end subroutine InitListCons
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine ListItem(FileUnit, ChaRep, Name)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ListItem(FileUnit, ChaRep, Name)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the output process for every single constant in 
!                        subroutine InitListCons.
! KEYWORDS: Output single constant.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Output single constant
!
!     Input: FileUnit --> The integer file handle for the output 
!            ChaRep   --> The character containing the values of the constant;
!            Name     --> Character for the name of the single constant.
!
!     Outpt: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
            use TimeRecord
            implicit none
!______________________________________________________________________________________________________________     
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
            integer FileUnit      ! The integer file handle for the output 
            character(*) ChaRep   ! The character containing the values of the constant
            character(*) Name     ! Character for the name of the single constant
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
            integer SkipLen       ! Length of skipping part
            character(200) FileNameA, FileNameB
!______________________________________________________________________________________________________________     
!_______________________________ Main calculations of Ouputing constant _______________________________________
!______________________________________________________________________________________________________________
            SkipLen = len("                           ")
            SkipLen = SkipLen - len_Trim(Name)
            write(FileNameA, *) DatTimLen
            write(FileNameB, *) SkipLen
            write(FileUnit, "(" // adjustl(FileNameA) // "x, 4x, " // adjustl(FileNameB) // "x, A, ' = ', A)") &
            & Trim(Name), Trim(ChaRep)
            
   end subroutine ListItem
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$