!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to construct the lattice geometry for the finite-size lattice we applied to
!              simulate with CPMC method.
! COMMENT: Construct Lattice Geometry.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   InitLattStru --> Subroutine used to allocate all the necessary quantities related with the lattice geometry;
!
!   GenerateList --> Subroutine used to generate all the lists used to construct the lattice geometry.
!   Bound_Matrix --> Subroutine used to generate the matrices used for boundary conditions.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine InitLattStru()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  InitLattStru()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the initialization for the lattice geometry we applied in the CPMC 
!                     calculations.
! KEYWORDS: Initialization of CPMC finite lattice.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Initialization of CPMC finite lattice, including:
!             (0) Some constants for the lattice;
!             (1) Allocate the corresponding arrays;
!             (2) Calculate the arrays for using in CPMC.
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
!_______________________________ Main calculations of CPMC finite lattice _____________________________________
!______________________________________________________________________________________________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
!_________________________ Monitor output of initialization process _______________________________
!&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(&(
      if(amyid == amstr) then
         write(*, "()")
         write(*, "(16x, 'InitLattStru: Initialization of the lattice geometry!')")
      end if
!**************************************************************************************************	  
!___________________ 0. Some constants for the lattice ____________________________________________
!**************************************************************************************************
		A1Vec(1) = 1.0_rp; A1Vec(2) = 0.0_rp
      A2Vec(1) = 0.0_rp; A2Vec(2) = 1.0_rp
      
      B1Vec(1) = 1.0_rp; B1Vec(2) = 0.0_rp
      B2Vec(1) = 0.0_rp; B2Vec(2) = 1.0_rp
		B1Vec = 2.0_rp * rp_pi * B1Vec / NumL1
      B2Vec = 2.0_rp * rp_pi * B2Vec / NumL2
!**************************************************************************************************	  
!___________________ 1. Allocate and the arrays for the finite lattice ____________________________
!**************************************************************************************************
!________________________________________________________________________________________ 	  
!_________________ (1). Basic lattice information _______________________________________
!________________________________________________________________________________________ 
      allocate(   StList(NumNS,     2))         ! Return the unit cell integer index for site and 1 or 2 for sublattice
		allocate(InvStList(NumL1, NumL2))         ! From UC index and Orbit to get the lattice site index
		allocate(   KpList(NumNC, 2))             ! The two integer indexes for a single k point
		allocate(InvKpList(0:NumL1-1, 0:NumL2-1)) ! From two integer indexes to get k point number
		StList    = 0
		InvStList = 0
		KpList    = 0
		InvKpList = 0
!______________________________________________________________________________________	  
!_________________ (2). Distance indexes for unit cells in the finite lattice _________
!______________________________________________________________________________________
		allocate(IminusJ(NumNC, NumNC))         ! The integer index for (i-j) and i,j are unit cell indexes
		IminusJ = 0
!______________________________________________________________________________________  
!_________________ (3). The Nearest-neighbor and Next-nearest-neighbor bonds __________
!______________________________________________________________________________________
		allocate(FNNBond(NumNS, 4))
		allocate(SNNBond(NumNS, 4))  
      	allocate(TNNBond(NumNS, 4))  
      	FNNBond = 0
     	SNNBond = 0
      	TNNBond = 0
!______________________________________________________________________________________ 
!_________________ (4). Open or periodic boundary conditions __________________________
!______________________________________________________________________________________
		allocate(FNNStBnd(NumNS, 4))
		allocate(SNNStBnd(NumNS, 4))
		allocate(TNNStBnd(NumNS, 4))
      	FNNStBnd = 1
      	SNNStBnd = 1
      	TNNStBnd = 1
!**************************************************************************************************	  
!___________________ 2. Generate the Boundary Condition matrix for the calculations _______________
!**************************************************************************************************
		call Bound_Matrix()
!**************************************************************************************************	  
!___________________ 3. Calculate all the above defined lists for the finite lattice ______________
!**************************************************************************************************
		call GenerateList()
      
	end subroutine InitLattStru
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine GenerateList() 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  GenerateList()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to generate the lists matrices defined in the CoreParamt module for square lattice.
! KEYWORDS: List matrices.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: In this subroutine, we need to evaluate all the following matrices: 
!            (4) StList(NumNS, 2)
!            (5) InvStList(NumL1, NumL2)
!            (6) NNBond(NumNS, 0:3)
!            (7) TNNBond(NumNS, 0:3)
!            (8) KpList(NumNC, 2)
!            (9) InvKpList(-NumL1 : NumL1, -NumL2 : NumL2)
!
!     And we adopt the periodic boundary conditions.
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
		use StdInOutSt
		implicit none
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1    ! Loop integer
		integer I2    ! Loop integer
		integer Id    ! Loop integer
		integer Jd    ! Loop integer
		integer Itp1  ! Temporary integer used in calculations
		integer Itp2  ! Temporary integer used in calculations
		integer Itp3  ! Temporary integer used in calculations
		integer Itp4  ! Temporary integer used in calculations
      
      	integer NT
      	integer NTMnus
      	integer NTPlus
		
		integer Idx       ! Index integer for Id-th unit cell in a1 direction
		integer Idy       ! Index integer for Id-th unit cell in a2 direction
		integer Jdx       ! Index integer for Jd-th unit cell in a1 direction
		integer Jdy       ! Index integer for Jd-th unit cell in a2 direction
      
      	integer Nu
      	integer No
		
		integer ImJ       ! The integer index for the distance between two unit cells
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the List matrices _______________________________________
!______________________________________________________________________________________________________________
!____________________________________________________________________________________________________	  
!___________________ 0. Calculate StList(NumNS, 2) and InvStList(NumL1, NumL2) ______________________
!____________________________________________________________________________________________________
		Itp1 = 0
		do I2 = 1, NumL2, +1                        ! Loop over a2 direction unit cells
			do I1 = 1, NumL1, +1                    ! Loop over a1 direction unit cells
				Itp1 = Itp1 + 1                     ! Count the number of unit cells
				StList(Itp1, 1)   = I1              ! The number of a_1 in position vector of Itp1 unit cell 
				StList(Itp1, 2)   = I2              ! The number of a_2 in position vector of Itp1 unit cell
				InvStList(I1, I2) = Itp1            ! determine the number of unit cell according to the (I1, I2) integer
			enddo
      	enddo
!____________________________________________________________________________________________________	  
!___________________ 1. Calculate FNNBond(NumNS, 1:4) _______________________________________________
!____________________________________________________________________________________________________
		do Id = 1, NumNS
			I1 = StList(Id, 1)
			I2 = StList(Id, 2)
			
			Itp2 = I2 + 1
			if(Itp2 > NumL2) Itp2 = Itp2 - NumL2
			FNNBond(Id, 1) = InvStList(I1, Itp2)
			
			Itp1 = I1 + 1
			if(Itp1 > NumL1) Itp1 = Itp1 - NumL1
			FNNBond(Id, 2) = InvStList(Itp1, I2)
			
			Itp2 = I2 - 1
			if(Itp2 <= 0) Itp2 = Itp2 + NumL2
			FNNBond(Id, 3) = InvStList(I1, Itp2)
			
			Itp1 = I1 - 1
			if(Itp1 <= 0) Itp1 = Itp1 + NumL1
			FNNBond(Id, 4) = InvStList(Itp1, I2)
      enddo
!____________________________________________________________________________________________________	  
!___________________ 2. Calculate SNNBond(NumNS, 1:4) _______________________________________________
!____________________________________________________________________________________________________
		do Id = 1, NumNS
			I1 = StList(Id, 1)
			I2 = StList(Id, 2)
			
			Itp1 = I1 + 1
			Itp2 = I2 + 1
			if(Itp1 > NumL1) Itp1 = Itp1 - NumL1
			if(Itp2 > NumL2) Itp2 = Itp2 - NumL2
			SNNBond(Id, 1) = InvStList(Itp1, Itp2)
			
			Itp1 = I1 + 1
			Itp2 = I2 - 1
			if(Itp1 > NumL1) Itp1 = Itp1 - NumL1
			if(Itp2 <= 0   ) Itp2 = Itp2 + NumL2
			SNNBond(Id, 2) = InvStList(Itp1, Itp2)
			
			Itp1 = I1 - 1
			Itp2 = I2 - 1
			if(Itp1 <= 0) Itp1 = Itp1 + NumL1
			if(Itp2 <= 0) Itp2 = Itp2 + NumL2
			SNNBond(Id, 3) = InvStList(Itp1, Itp2)
			
			Itp1 = I1 - 1
			Itp2 = I2 + 1
			if(Itp1 <= 0   ) Itp1 = Itp1 + NumL1
			if(Itp2 > NumL2) Itp2 = Itp2 - NumL2
			SNNBond(Id, 4) = InvStList(Itp1, Itp2)
      enddo
!____________________________________________________________________________________________________	  
!___________________ 3. Calculate TNNBond(NumNS, 1:4) _______________________________________________
!____________________________________________________________________________________________________
      do Id = 1, NumNS, +1
			I1 = StList(Id, 1)
			I2 = StList(Id, 2)
			
			Itp2 = I2 + 2
			if(Itp2 > NumL2) Itp2 = Itp2 - NumL2
			TNNBond(Id, 1) = InvStList(I1, Itp2)
			
			Itp1 = I1 + 2
			if(Itp1 > NumL1) Itp1 = Itp1 - NumL1
			TNNBond(Id, 2) = InvStList(Itp1, I2)
			
			Itp2 = I2 - 2
			if(Itp2 <= 0) Itp2 = Itp2 + NumL2
			TNNBond(Id, 3) = InvStList(I1, Itp2)
			
			Itp1 = I1 - 2
			if(Itp1 <= 0) Itp1 = Itp1 + NumL1
			TNNBond(Id, 4) = InvStList(Itp1, I2)
      enddo
!____________________________________________________________________________________________________	  
!___________________ 4. Calculate KpList(NumNC, 2) and InvKpList(0 : NumL1, 0 : NumL2) ______________
!____________________________________________________________________________________________________      
      Itp1 = 0
      do I2 = 0, NumL2-1, +1
		   do I1 = 0, NumL1-1, +1
				Itp1 = I2*NumL1 + I1 + 1
				KpList(Itp1, 1) = I1
				KpList(Itp1, 2) = I2
				InvKpList(I1, I2) = Itp1
			enddo
		enddo
!____________________________________________________________________________________________________	  
!___________________ 5. Calculate IminusJ(NumNC, NumNC) matrix ______________________________________
!____________________________________________________________________________________________________
		do Id = 1, NumNC
			do Jd = 1, NumNC
				Idx = StList(Id, 1)
				Idy = StList(Id, 2)
				Jdx = StList(Jd, 1)
				Jdy = StList(Jd, 2)
				I1 = Idx - Jdx
				I2 = Idy - Jdy
				
				if( I1*1.0_rp < -NumL1/2.0_rp+1.0E-8_rp ) then
					I1 = I1 + NumL1
				else if( I1*1.0_rp > NumL1/2.0_rp+1.0E-8_rp ) then
					I1 = I1 - NumL1
				end if
				
				if( I2*1.0_rp < -NumL2/2.0_rp+1.0E-8_rp ) then
					I2 = I2 + NumL2
				else if( I2*1.0_rp > NumL2/2.0_rp+1.0E-8_rp ) then
					I2 = I2 - NumL2
				end if
				
				I1 = I1 + NumL1/2 + mod(NumL1, 2)
				I2 = I2 + NumL2/2 + mod(NumL2, 2)
				
				ImJ = InvStList(I1, I2)
				
				IminusJ(Id, Jd) = ImJ
			enddo
      enddo   
      
   end subroutine GenerateList
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine Bound_Matrix() 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  Bound_Matrix()
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to construct the NNStBnd, NNUCBnd and TNNsBnd integer matrices.
! KEYWORDS: Calculate the switch integers.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: In this subroutine, we calculate the NNStBnd, NNUCBnd and TNNsBnd integer matrices.
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
		integer I1                       ! Loop integer for all lattice sites
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the List matrices _______________________________________
!______________________________________________________________________________________________________________
		FNNStBnd = 0                      ! Initialization
		SNNStBnd = 0                      ! Initialization
		TNNStBnd = 0                      ! Initialization
!*****************************************************************************************************	  
!___________________ 0. Construct BndMat, HbUMat, ExJMat, HbVMat matrices for both PBC _______________
!*****************************************************************************************************
		FNNStBnd = 1                      ! For PBC condition
		SNNStBnd = 1                      ! For PBC condition 
		TNNStBnd = 1                      ! For PBC condition
!*****************************************************************************************************	  
!___________________ 1. Construct the BndMat matrix for Open boundary condition case _________________
!*****************************************************************************************************
!__________________________________________________________________________________________	  
!________________________ (0) OBC in a1 direction for the system __________________________
!__________________________________________________________________________________________
		if(.not. IfPr1) then
!____________________________________________________________________________  
!____________________ <0> For 1NN lattice sites ______________________________
!____________________________________________________________________________
			do I1 = NumL1, NumNC, NumL1
				FNNStBnd(        I1, 2) = 0
				FNNStBnd(I1-NumL1+1, 4) = 0
			enddo
!____________________________________________________________________________  
!____________________ <1> For 2NN lattice sites _____________________________
!____________________________________________________________________________
			do I1 = NumL1, NumNC, NumL1
				SNNStBnd(        I1, 1) = 0
				SNNStBnd(        I1, 2) = 0
				SNNStBnd(I1-NumL1+1, 3) = 0
				SNNStBnd(I1-NumL1+1, 4) = 0
			enddo
!____________________________________________________________________________  
!____________________ <2> For 3NN lattice sites _____________________________
!____________________________________________________________________________
			do I1 = NumL1, NumNC, NumL1
				TNNStBnd(        I1, 2) = 0
				TNNStBnd(      I1-1, 2) = 0
				TNNStBnd(I1-NumL1+1, 4) = 0
				TNNStBnd(I1-NumL1+2, 4) = 0
			enddo
		end if
!__________________________________________________________________________________________	  
!________________________ (1) OBC in a2 direction for the system __________________________
!__________________________________________________________________________________________
		if(.not. IfPr2) then
!____________________________________________________________________________  
!____________________ <0> For NN lattice sites ______________________________
!____________________________________________________________________________	
			do I1 = 1, NumL1, 1
				FNNStBnd(            I1, 3) = 0
				FNNStBnd(NumNC-NumL1+I1, 1) = 0
			enddo
!____________________________________________________________________________  
!____________________ <1> For NNN lattice sites _____________________________
!____________________________________________________________________________
			do I1 = 1, NumL1, 1
				SNNStBnd(            I1, 2) = 0
				SNNStBnd(            I1, 3) = 0
				SNNStBnd(NumNC-NumL1+I1, 1) = 0
				SNNStBnd(NumNC-NumL1+I1, 4) = 0
			enddo
!____________________________________________________________________________  
!____________________ <2> For 3NN lattice sites _____________________________
!____________________________________________________________________________
			do I1 = 1, NumL1, 1
				TNNStBnd(              I1, 3) = 0
				TNNStBnd(        NumL1+I1, 3) = 0
				TNNStBnd(NumNC-  NumL1+I1, 1) = 0
				TNNStBnd(NumNC-2*NumL1+I1, 1) = 0
			enddo
		end if
		
   end subroutine Bound_Matrix
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
