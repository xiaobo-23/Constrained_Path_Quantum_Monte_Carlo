!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A few subroutines used for storing and reading the path (configuration) fields from the constrained-path
!               constructions.
! COMMENT: Read and store the path (configuration) fields.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-06-15
! PURPOSE: Different subroutines are introduced as following:
!             
!   SaveFewCfg     --> Subroutine used to save few paths on every process from the CP construction;
!   SaveReadPathWt --> Subroutine used to save and read the paths stored for measurements;
!
!   BinaryToDecimal --> Subroutine used to transform binary numbers to -1/+1 auxiliary fields;
!   DecimalToBinary --> Subroutine used to transform -1/+1 auxiliary fields to binary numbers.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine SaveFewCfg()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SaveFewCfg()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to output the auxiliary-field configurations after single BIN simulation.
! KEYWORDS: Save auxiliary fields.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Since there might be ~100 random walkers on an single process, and each of them might have an
!                 independent configuration, here we only store the path with the largest weight.
!
!     Input: (none); Output: (none)
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
      integer It, I0, I1, WalkIndx, WkIdMaxWt
      real(rp) Rtp0
      character(010) ProcessID
      character(200) FileName
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of outputting the auxiliary fields _____________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Store the auxiliary fields for HubbU interaction __________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Find the random walker with the largest weight ___________________
!________________________________________________________________________________________
      WkIdMaxWt = 1; Rtp0 = WghtProc(IdptWkIndx(1))
      do WalkIndx = 2, NWkBt, +1
         if(Rtp0 < WghtProc(IdptWkIndx(WalkIndx))) then
            Rtp0 = WghtProc(IdptWkIndx(WalkIndx))
            WkIdMaxWt = WalkIndx
         end if
      enddo
!________________________________________________________________________________________ 	  
!_________________ (1) Transform the fields to 0(1), store in RdWrtField ________________
!_____________________ Take RdWrtField as binary numbers, convert to decimal ____________
!________________________________________________________________________________________
      do It = 1, LTrot, +1
         do I0 = 1, NumNS, +1
            I1 = (It-1)*NumNS + I0
            RdWrtField(I1) = ( IsingbU(I0, It, IdptWkIndx(WkIdMaxWt)) + 1 ) / 2
         enddo
      enddo
      call BinaryToDecimal()
!________________________________________________________________________________________ 	  
!_________________ (2) Output the fields to binary files for every process ______________
!________________________________________________________________________________________
      write(ProcessID, "(I3.3)") amyid
      FileName = "Output/AuxiliaryFlds/ZZ_HubbUIsingFd_" // Trim(ProcessID) // ".bin"
      open( amyid+291, file = Trim(FileName), access = "direct", form = "unformatted", recl = ReclUnitNm)
      write(amyid+291, rec = 1) RdWrtIntgr(1:NmBinaryFd:+1)
      close(amyid+291)
      
   end subroutine SaveFewCfg
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine SaveReadPathWt(SaveOrRead)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SaveReadPathWt(SaveOrRead)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to output the auxiliary fields between certain time slices.
! KEYWORDS: Save auxiliary fields.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Read or Save the CP sampled paths and their weights.
!
!     Input: SaveOrRead --> == "Save" for saving configurations; == "Read" for reading configurations.
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
      use Observable
		use MPISetting
      implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      character(4) SaveOrRead
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer It, I0, I1, WalkIndx, NWalkTmp
      character(010) ProcessID
      character(200) FileName
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of outputting the auxiliary fields _____________________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Save all the CP sampled (independent) paths and their weights _____________
!**************************************************************************************************
      if(SaveOrRead == "Save") then
!________________________________________________________________________________________ 	  
!_________________ (0) Store NWkBt, NWalk and the weights of the walkers ________________
!________________________________________________________________________________________
         write(FileName, "('WeightOfNWlk_Proc', I3.3, '.txt')") amyid
         open( amyid+272, file = "Output/AuxiliaryFlds_Measure/" // Trim(FileName), access = "append")
         write(amyid+272, "(2I4)", advance = "no") NWalk, NWkBt
         do WalkIndx = 1, NWkBt,  +1
            write(amyid+272, "(es20.12)", advance = "no") WghtProc(IdptWkIndx(WalkIndx))
         enddo
         do WalkIndx = NWkBt+1, NWalk,  +1
            write(amyid+272, "(es20.12)", advance = "no") 0.0_rp
         enddo
         close(amyid+272)
!________________________________________________________________________________________ 	  
!_________________ (1) Store the fields of all NWkBt independent walkers ________________
!________________________________________________________________________________________
         do WalkIndx = 1, NWkBt,  +1
            !!!!!!!!!! Transform the fields to 0(1), store in RdWrtField
            do It = 1, LTrot, +1
               do I0 = 1, NumNS, +1
                  I1 = (It-1)*NumNS + I0
                  RdWrtField(I1) = ( IsingbU(I0, It, IdptWkIndx(WalkIndx)) + 1 ) / 2
               enddo
            enddo
            !!!!!!!!!! Take RdWrtField as binary numbers, convert to decimal
            call BinaryToDecimal()
            !!!!!!!!!! Output the fields to binary files for every process
            write(ProcessID, "(I3.3)") amyid
            FileName = "Output/AuxiliaryFlds_Measure/AuxiliaryFld_Proc" // Trim(ProcessID) // ".bin"
            FldIndRdWt = FldIndRdWt + 1
            open( amyid+292, file = Trim(FileName), access = "direct", form = "unformatted", recl = ReclUnitNm)
            write(amyid+292, rec = FldIndRdWt) RdWrtIntgr(1:NmBinaryFd:+1)
            close(amyid+292)
         enddo
      end if
!**************************************************************************************************	  
!___________________ 1. Read the previously stored paths and their weights ________________________
!**************************************************************************************************
      if(SaveOrRead == "Read") then
!________________________________________________________________________________________ 	  
!_________________ (0) Read NWkBt, NWalk and the weights of the walkers _________________
!________________________________________________________________________________________
         !!!!!!!!!! Read the weights of random walkers 
         WghtIdRead = WghtIdRead + 1
         write(FileName, "('WeightOfNWlk_Proc', I3.3, '.txt')") amyid
         open( amyid+272, file = Trim(MeaFldPath) // "/AuxiliaryFlds_Measure/" // Trim(FileName), &
            & access = "direct", recl = (2*4+NWalk*20)+1, form = "formatted", status = "old")
         write(FileName, *) NWalk
         read( amyid+272, fmt = "(2I4, " // adjustl(FileName) // "es20.12" // ")", rec = WghtIdRead) &
            & NWalkTmp, NWkBt, WghtProc(1:NWalk:+1)
         close(amyid+272)
         !!!!!!!!!! Define the array of IdptWkIndx
         IdptWkIndx = -10000
         do WalkIndx = 1, NWkBt, +1
            IdptWkIndx(WalkIndx) = WalkIndx
         enddo
!________________________________________________________________________________________ 	  
!_________________ (1) Read the fields of all NWkBt independent walkers _________________
!________________________________________________________________________________________
         write(ProcessID, "(I3.3)") amyid
         FileName = Trim(MeaFldPath) // "/AuxiliaryFlds_Measure/AuxiliaryFld_Proc" // Trim(ProcessID) // ".bin"
         open(amyid+292, file = Trim(FileName), access = "direct", form = "unformatted", &
            & recl = ReclUnitNm, status = "old")
         do WalkIndx = 1, NWkBt, +1
            !!!!!!!! Read the integers from binary files
            FldIndRdWt = FldIndRdWt + 1
            read(amyid+292, rec = FldIndRdWt) RdWrtIntgr(1:NmBinaryFd:+1)
            !!!!!!!! Convert decimal numbers into binary numbers
            call DecimalToBinary()
            !!!!!!!! Restore -1 and +1 from 0 and 1
            do It = 1, LTrot, +1
               do I0 = 1, NumNS, +1
                  I1 = (It-1)*NumNS + I0
                  IsingbU(I0, It, WalkIndx) = 2 * RdWrtField(I1) - 1
               enddo
            enddo
         enddo
         close(amyid+292)
      end if
      
   end subroutine SaveReadPathWt
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine BinaryToDecimal()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  BinaryToDecimal()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to onvert the binary number into decimal number. 
! KEYWORDS: Number convertion.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION:
!
!     Convert the binary number into decimal number. RdWrtField --> RdWrtIntgr.
!
!     Input: (none); Output: (none)
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
      integer I0, I1, Base
      integer IndBgn, IndEnd
      character(LengthBite) FileStrng
      character(50) Tmp1Strng, Tmp2Strng, Tmp3Strng
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of converting binary number to decimal number __________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Convert a binary number into decimal number _______________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Initialization for the result RdWrtIntgr _________________________
!________________________________________________________________________________________
      RdWrtIntgr = 0
!________________________________________________________________________________________ 	  
!_________________ (1) Iteration to obtain results of RdWrtIntgr ________________________
!________________________________________________________________________________________
      do I0 = 1, NmBinaryFd, +1
         IndBgn = (I0-1)*LengthBite + 1
         IndEnd =     I0*LengthBite
         if( (I0 == NmBinaryFd) .and. (mod(LTrotNumNS, LengthBite) /= 0) ) then
            IndEnd = LTrotNumNS
            write(Tmp1Strng, *) LengthBite-mod(LTrotNumNS, LengthBite)
            write(Tmp2Strng, *) mod(LTrotNumNS, LengthBite)
            write(FileStrng, "(I" // adjustl(Tmp1Strng) // "." // adjustl(Tmp1Strng) // ", " // &
               & adjustl(Tmp2Strng) // "I1)") 0, RdWrtField(IndBgn:IndEnd:+1)
         else
            write(Tmp3Strng, *) LengthBite
            write(FileStrng, "(" // adjustl(Tmp3Strng) // "I1)") RdWrtField(IndBgn:IndEnd:+1)
         end if
         Base = 1
         do I1 = LengthBite, 1, -1
            if(FileStrng(I1:I1) == "1") then
               RdWrtIntgr(I0) = RdWrtIntgr(I0) + Base
            end if
            Base = Base * 2
         enddo
      enddo
      
   end subroutine BinaryToDecimal
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine DecimalToBinary()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  DecimalToBinary()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to onvert the decimal number into binary number. 
! KEYWORDS: Number convertion.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION:
!
!     Convert the decimal number into binary number. RdWrtIntgr --> RdWrtField.
!
!     Input: (none); Output: (none)
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
      integer I0, I1, SiteInd
      integer IndBgn, IndEnd
!______________________________________________________________________________________________________________	  
!___________________________ Main calculations of converting binary number to decimal number __________________
!______________________________________________________________________________________________________________ 
!**************************************************************************************************	  
!___________________ 0. Convert a binary number into decimal number _______________________________
!**************************************************************************************************
      SiteInd = 0
      do I0 = 1, NmBinaryFd, +1
         if( (I0 == NmBinaryFd) .and. (mod(LTrotNumNS, LengthBite) /= 0) ) then
            IndBgn = - mod(LTrotNumNS, LengthBite) + 1
            IndEnd = 0
         else
            IndBgn = -LengthBite + 1
            IndEnd = 0
         end if
         do I1 = IndBgn, IndEnd, +1
            SiteInd = SiteInd + 1
            RdWrtField(SiteInd) = iand(ishft(RdWrtIntgr(I0), I1), 1)
         enddo
      enddo
      
   end subroutine DecimalToBinary
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
