!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to perform Initialization of exp(-\Delta\tau\hat{V}) related quantities before the
!               whole CPMC simulations.
! COMMENT: CPMC Initialization process.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!      Since we always write the simulated Hamiltonian as H = H_T + (H_U + H_0 - H_T), here in these subroutines, we
!      deal with exp[-dt*(H_U + H_0 - H_T)]. Since [H_U, H_0-H_T] = 0, So we we have
!        exp[-dt*(H_U + H_0 - H_T)] = exp(-dt*H_U) * exp[-dt*(H_0 - H_T)].
!             
!   InitExpHmltV --> Subroutine to perform initialization of HS Tranformation related quantities used in CPMC;
!
!   InitUpdatebU --> Initialization of updating related quantities;
!   ResetDeltabU --> Reset the quantities by the (H_0-H_T) Hamiltonian.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine InitExpHmltV()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitExpHmltV()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the initialization for auxiliary fields and the updating related 
!                  quantities.
! KEYWORDS: Initialization of constants in CPMC.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Initialization of Hubbard-Stratonovich Transformation quantities and auxiliary fields.
!
!     For both B_X and B_T;
!
!     For B_T: (0) Choose RHF or UHF method and determine corresponding HF_Correction value;
!
!     For B_X: (0) Allocate the Field matrice;
!              (1) Generate the HS-related quantities.
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
      use RandomNumb
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!__________________________ Main calculations of quantities in HS Transformation ______________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of initialization process _______________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'InitExpHmltV: Initialization of exp(-Dltau*H_I)!')")
      end if
!**************************************************************************************************	  
!___________________ 0. Set up the V term related quantities in B_X matrix ________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Initializations for the field path and related ___________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         write(*, "(30x, 'Set up the auxiliary fields related quantities!')")
      end if
      call InitFldConfg()
!________________________________________________________________________________________ 	  
!_________________ (1) Calculate Updating related quantities for interactions ___________
!________________________________________________________________________________________
      if(amyid == amstr) then
         write(*, "(30x, 'Set up the Hubbard-Stratonovich transformations for interactions!')")
      end if
      call InitUpdatebU()
!________________________________________________________________________________________ 	  
!_________________ (2) Reset the DeltbU_H0T matrix by (H_0-H_T) factor __________________
!________________________________________________________________________________________
      if(amyid == amstr) then
         write(*, "(30x, 'Set up the \Delta matrix used for updating auxiliary fields!')")
      end if
      call ResetDeltabU()
      
   end subroutine InitExpHmltV
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine InitFldConfg()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitFldConfg()
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to initialize the field configurations and the related reading and storing 
!               of these fields.
! KEYWORDS: Field configuration initializations.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-06-15
! DESCRIPTION: Initialize the field configuration, and reading and storing related.
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
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer  anprcTmp, NumL1Tmp, NumL2Tmp, LTrotTmp, NWalkTmp
      real(rp) MuConst, ExpMuConst
!______________________________________________________________________________________________________________	  
!__________________________ Main calculations of updating related quantities __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Initialization of field configuration and related _________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Allocate the field configuration matrix __________________________
!________________________________________________________________________________________
      allocate(IsingbU(NumNS, LTrot, NWalk))
      IsingbU = +1
!________________________________________________________________________________________ 	  
!_________________ (1) For storing or reading the auxiliary fields ______________________
!________________________________________________________________________________________
      LengthBite = 31
      LTrotNumNS = LTrot * NumNS
      NmBinaryFd = LTrotNumNS / LengthBite
      if( mod(LTrotNumNS, LengthBite) /= 0 ) then
         NmBinaryFd = NmBinaryFd + 1
      end if
      allocate(RdWrtField(LTrotNumNS))
      allocate(RdWrtIntgr(NmBinaryFd))
      inquire(iolength=ReclUnitNm) RdWrtIntgr
      WghtIdRead = 0; FldIndRdWt = 0
!**************************************************************************************************	  
!___________________ 1. Store or check the information of the paths _______________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For IfSaveFlds == T case, store the information __________________
!________________________________________________________________________________________
      if( (IfSaveFlds) .and. (amyid == amstr) ) then
         open( 197, file = "Output/AuxiliaryFlds/AuxiliaryFieldInfo.txt")
         write(197, "(I4, A, I4, A, I4, A, I6)") anprc, char(9), NumL1, char(9), NumL2, char(9), LTrot
         close(197)
      end if
!________________________________________________________________________________________ 	  
!_________________ (1) For SaveFldMea == T case, store the information __________________
!________________________________________________________________________________________
      if( (SaveFldMea) .and. (amyid == amstr) ) then
         open( 197, file = "Output/AuxiliaryFlds_Measure/AuxiliaryFieldInfo.txt")
         write(197, "(I4, A, I4, A, I4, A, I6, A, I4)") anprc, char(9), NumL1, char(9), NumL2, char(9), LTrot, &
            & char(9), NWalk
         close(197)
      end if
!________________________________________________________________________________________ 	  
!_________________ (2) For ReadFldMea == T case, read the information ___________________
!________________________________________________________________________________________
      if(ReadFldMea) then
         !!!!!!!!!! Read information from input file
         open( amyid+197, err = 377, file = Trim(MeaFldPath) // "/AuxiliaryFlds_Measure/AuxiliaryFieldInfo.txt", &
            & status = "old")
         read( amyid+197, *) anprcTmp, NumL1Tmp, NumL2Tmp, LTrotTmp, NWalkTmp
         close(amyid+197)
         if( (anprcTmp /= anprc) .or. (NumL1Tmp /= NumL1) .or. (NumL2Tmp /= NumL2) .or. (LTrotTmp /= LTrot) .or. &
               & (NWalkTmp /= NWalk) ) then
            write(*, "(35x, 'Input auxiliary fields are not compatible, will stop runnning!!! amyid = ', I4)") amyid
            stop
         end if
         go to 399
         !!!!!!!!!! If the input file does not exist
377      continue
         write(*, "(35x, 'Input file .../AuxiliaryFlds_Measure/AuxiliaryFieldInfo.txt does not exist!!! Stop!!!')")
         stop
399      continue
      end if

   end subroutine InitFldConfg
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine InitUpdatebU()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitUpdatebU()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the initialization for HubbU interaction updating related quantities.
! KEYWORDS: Initialization of updating in CPMC.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Initialization of updating related quantities.
!
!     Since we always write the simulated Hamiltonian as H = H_T + (H_U + H_0 - H_T)
!     the additional (H_0 - H_T) term simply enters into the H_U term.
!     For the case of RHF type of H_T, (H_0 - H_T) is simply a chemical potential term;
!     For the case of UHF type of H_T, (H_0 - H_T) is a chemical potential term plus density terms.
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
		use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, ISite, Ix, Iy, SiteParity
      real(rp) TmpMat(2)
      real(rp) ExpobU_Init(   2, -1:+1)
      real(rp) ExpobUInv_Init(2, -1:+1)
!______________________________________________________________________________________________________________	  
!__________________________ Main calculations of updating related quantities __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Initialization of the updating related quantities _________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) For the coupling constant LmdbU __________________________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [0] For HubbU > 0 --> HS_Type = 1 _________________________
!____________________________________________________________________________
      if(HubbU > +rp_Eps) then
         LmdbU = acosh(exp(+Dltau*HubbU/2.0_rp))
!____________________________________________________________________________ 	  
!________________ [1] For HubbU < 0 --> HS_Type = 0 _________________________
!____________________________________________________________________________
      else if(HubbU < -rp_Eps) then
         LmdbU = acosh(exp(-Dltau*HubbU/2.0_rp))
!____________________________________________________________________________ 	  
!________________ [2] For HubbU = 0 --> LmdbU = 0.0_rp ______________________
!____________________________________________________________________________
      else
         LmdbU = 0.0_rp
      end if
!________________________________________________________________________________________ 	  
!_________________ (1) Initialize ExpobU_Init and ExpobUInv_Init ________________________
!________________________________________________________________________________________
      ExpobU_Init = 0.0_rp
      ExpobUInv_Init = 0.0_rp
!________________________________________________________________________________________ 	  
!_________________ (2) ExpobU, ExpobUInv, DeltbU, DeltbUCnst for updates ________________
!________________________________________________________________________________________
      allocate(ExpobU(   2, -1:+1, NumNS)); ExpobU     = 0.0_rp
      allocate(ExpobUInv(2, -1:+1, NumNS)); ExpobUInv  = 0.0_rp
!________________________________________________________________________________________ 	  
!_________________ (3) The ExpobU_H0T and ExpobUInv_H0T constants for updates ___________
!________________________________________________________________________________________
      allocate(ExpobU_H0T   (2, -1:+1, NumNS)); ExpobU_H0T    = 0.0_rp
      allocate(ExpobUInv_H0T(2, -1:+1, NumNS)); ExpobUInv_H0T = 0.0_rp
      allocate(DeltbU_H0T   (2, -1:+1, NumNS)); DeltbU_H0T    = 0.0_rp
!**************************************************************************************************	  
!___________________ 1. Obtain LmdbU, ExpobU, ExpobUInv for HS transformations ____________________
!**************************************************************************************************      
      if( abs(HubbU) > rp_Eps ) then
!________________________________________________________________________________________ 	  
!_________________ (0) For HS Transf into CDW channel --> HubbU < 0 _____________________
!_____________________ For HS Transf into Sz  channel --> HubbU > 0 _____________________
!________________________________________________________________________________________
!____________________________________________________________________________ 	  
!________________ [0] ExpobU_Init and ExpobUInv_Init ________________________
!____________________________________________________________________________
         if(HS_Type == 0) then
            do I1 = -1, +1, +2
               ExpobU_Init   (1, I1) = exp( + dble(I1) * LmdbU ) 
               ExpobU_Init   (2, I1) = exp( + dble(I1) * LmdbU ) 
               ExpobUInv_Init(1, I1) = exp( - dble(I1) * LmdbU ) 
               ExpobUInv_Init(2, I1) = exp( - dble(I1) * LmdbU ) 
            enddo
         else if(HS_Type == 1) then
            do I1 = -1, +1, +2
               ExpobU_Init   (1, I1) = exp( + dble(I1) * LmdbU ) 
               ExpobU_Init   (2, I1) = exp( - dble(I1) * LmdbU ) 
               ExpobUInv_Init(1, I1) = exp( - dble(I1) * LmdbU ) 
               ExpobUInv_Init(2, I1) = exp( + dble(I1) * LmdbU ) 
            enddo
         end if
!____________________________________________________________________________ 	  
!________________ [1] Output the HS transformation information ______________
!____________________________________________________________________________            
         if(amyid == amstr) then
            if(HS_Type == 0) then
               write(*, "(35x, 'HubbU --> Apply HS transformation into CDW channel!')")
               write(*, "(35x, '          With LmdbU = ', sp, es23.16)") LmdbU
            else if(HS_Type == 1) then
               write(*, "(35x, 'HubbU --> Apply HS transformation into SDW of Sz channel!')")
               write(*, "(35x, '          With LmdbU = ', sp, es23.16)") LmdbU
            end if
         end if
      end if
!**************************************************************************************************	  
!___________________ 2. Include the AFM pinning fields for PBC and FFT case _______________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (0) Obtain ExpobU and ExpobUInv within Pinning fields ________________
!________________________________________________________________________________________
      if( (abs(PinSz) >= rp_Eps) .and. (PinSzType == 2) .and. (FFTEXPDTH0) ) then
         if(amyid == amstr) then
            write(*, "(35x, '          Rescale ExpobU and ExpobUInv by the AFM pinning fields!')")
         end if
         do ISite = 1, NumNS, +1
            Ix = mod(ISite-1, NumL1) + 1
            Iy = (ISite-1)/NumL1 + 1
            SiteParity = (-1)**(mod(Ix+Iy, 2))
            do I1 = -1, +1, +2
               ExpobU(1, I1, ISite) = ExpobU_Init(1, I1) * exp(-Dltau*PinSz/2.0_rp*dble(SiteParity))
               ExpobU(2, I1, ISite) = ExpobU_Init(2, I1) * exp(+Dltau*PinSz/2.0_rp*dble(SiteParity))
               ExpobUInv(1, I1, ISite) = exp(+Dltau*PinSz/2.0_rp*dble(SiteParity)) * ExpobUInv_Init(1, I1)
               ExpobUInv(2, I1, ISite) = exp(-Dltau*PinSz/2.0_rp*dble(SiteParity)) * ExpobUInv_Init(2, I1)
            enddo
         enddo
      else
         do ISite = 1, NumNS, +1
            do I1 = -1, +1, +2
               ExpobU(1, I1, ISite) = ExpobU_Init(1, I1)
               ExpobU(2, I1, ISite) = ExpobU_Init(2, I1)
               ExpobUInv(1, I1, ISite) = ExpobUInv_Init(1, I1)
               ExpobUInv(2, I1, ISite) = ExpobUInv_Init(2, I1)
            enddo
         enddo
      end if
!**************************************************************************************************	  
!___________________ 3. Calculate ExpobU_H0T and ExpobUInv_H0T as _________________________________
!______________________ ExpobU_H0T    = ExpobU_Init * exp[-dt*(H_0-H_T)] __________________________
!______________________ ExpobUInv_H0T = exp[+dt*(H_0-H_T)] * ExpobUInv_Init _______________________
!______________________ Here, (H_0-H_T) term doesn't contain diagonal terms _______________________
!**************************************************************************************************
      if( abs(HubbU) > rp_Eps ) then
!________________________________________________________________________________________ 	  
!_________________ (0) Reset ExpobU_H0T --> ExpobU * exp(-dt*(H_0-H_T)) _________________
!________________________________________________________________________________________         
         ExpobU_H0T = 0.0_rp
         do ISite = 1, NumNS, +1
!____________________________________________________________________________ 	  
!________________ [0] First obtain the (H_0-H_T) Diag matrix ________________
!____________________________________________________________________________
            TmpMat(1) = exp( + Dltau * HubbU_UHF * MagMoment(ISite, 2) )
            TmpMat(2) = exp( + Dltau * HubbU_UHF * MagMoment(ISite, 1) )
!____________________________________________________________________________ 	  
!________________ [1] Obtain ExpobU_H0T Diag matrix on ISite-th site ________
!____________________________________________________________________________
            do I1 = -1, +1, 2
               ExpobU_H0T(1, I1, ISite) = ExpobU_Init(1, I1)*TmpMat(1)
               ExpobU_H0T(2, I1, ISite) = ExpobU_Init(2, I1)*TmpMat(2)
            enddo
         enddo
!________________________________________________________________________________________ 	  
!_________________ (1) Reset ExpobUInv_H0T --> exp[+dt*(H_0-H_T)]* ExpobUInv_Init _______
!________________________________________________________________________________________          
         ExpobUInv_H0T = 0.0_rp
         do ISite = 1, NumNS, +1
!____________________________________________________________________________ 	  
!________________ [0] First obtain the (H_0-H_T) diag matrix ________________
!____________________________________________________________________________
            TmpMat = 0.0_rp
            TmpMat(1) = exp( - Dltau * HubbU_UHF * MagMoment(ISite, 2) )
            TmpMat(2) = exp( - Dltau * HubbU_UHF * MagMoment(ISite, 1) )
!____________________________________________________________________________ 	  
!________________ [1] Obtain ExpobUInv_H0T Diag matrix on ISite-th site _____
!____________________________________________________________________________            
            do I1 = -1, +1, 2
               ExpobUInv_H0T(1, I1, ISite) = TmpMat(1)*ExpobUInv_Init(1, I1)
               ExpobUInv_H0T(2, I1, ISite) = TmpMat(2)*ExpobUInv_Init(2, I1)
            enddo
         enddo
      end if
            
   end subroutine InitUpdatebU
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine ResetDeltabU()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  ResetDeltabU()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to reset the DeltbU_H0T matrix in HS transformation by the additional 
!                chemical potential term (ChemP-ChemP_BT+HubbU_UHF/2.0_rp) in (H_0-H_T) Hamiltonian.
! KEYWORDS: Reset update related quantities.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Reset the DeltbU_H0T matrix by the (ChemP-ChemP_BT+HubbU_UHF/2.0_rp) term in (H_0-H_T) Hamiltonian.
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
		implicit none
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, ISite
      real(rp) MuConst, ExpMuConst
!______________________________________________________________________________________________________________	  
!__________________________ Main calculations of updating related quantities __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Reset the updating related quantities for HubbU interaction _______________
!**************************************************************************************************
      if( abs(HubbU) > rp_Eps ) then
!________________________________________________________________________________________ 	  
!_________________ (0) The Chemical potential term in (H_0-H_T) Hamlt ___________________
!_____________________ (ChemP-ChemP_BT+HubbU_UHF/2.0_rp) * Identity matrix ______________
!________________________________________________________________________________________
         MuConst = ChemP - ChemP_BT + HubbU_UHF/2.0_rp
         ExpMuConst = exp( - Dltau * MuConst )
!________________________________________________________________________________________ 	  
!_________________ (1) Reset DeltbU_H0T --> DeltbU_H0T = ExpobU_H0T * ExpMuConst - I ____
!________________________________________________________________________________________         
         DeltbU_H0T = 0.0_rp
         do ISite = 1, NumNS, +1
            do I1 = -1, 1, +2
               DeltbU_H0T(1, I1, ISite) = ExpobU_H0T(1, I1, ISite)*ExpMuConst - 1.0_rp
               DeltbU_H0T(2, I1, ISite) = ExpobU_H0T(2, I1, ISite)*ExpMuConst - 1.0_rp
            enddo
         enddo
      end if

   end subroutine ResetDeltabU
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine Obtain_Gamma(x, y)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  Obtain_Gamma(x, y)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the equation of cosh(x) = exp(y), to get the coupling constant.
! KEYWORDS: Solve cosh(x) = exp(y).
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Solve the equation of cosh(x) = exp(y), to get the coupling constant.
!
!     Input:  y --> Input real number;
!
!     Output: x --> Output real number.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		use RealPrecsn
      use MPISetting
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
      real(rp) x, y
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      real(rp) Rtp0
!______________________________________________________________________________________________________________	  
!__________________________ Main calculations of Solving cosh(x) = exp(y) _____________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Solve the equation ________________________________________________________
!************************************************************************************************** 
      if( abs(y) < rp_Eps ) then
         x = 0.0_rp
      else
         Rtp0 = exp(y)
         x = log( Rtp0 + sqrt(Rtp0**2-1.0_rp) )
      end if
!**************************************************************************************************	  
!___________________ 1. Check the solution ________________________________________________________
!************************************************************************************************** 
      Rtp0 = log( (exp(+x) + exp(-x))/2.0_rp )
      if( abs(Rtp0-y) > 1.0E-10_rp ) then
         write(*, "('Subroutine Obtain_Gamma: For process ', I3.3, ', Wrong solution!!!')") amyid
      end if
      
   end subroutine Obtain_Gamma
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$