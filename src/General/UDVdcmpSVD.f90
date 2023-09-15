!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to computes the UDV decomposition of general RDim-by-CDim matrix, with U(RDim, CDim), 
!              D(CDim, CDim) and V(CDim, CDim). D(CDim, CDim) is REAL diagonal matrix.
!          These subroutines calculates the UDV decomposition by the singular-value decomposition (SVD) method. So both U
!              and V matrices are unitary matrix. U(RDim, CDim) is part of unitary matrix (N columns). 
!
!          Based on LAPACK package.
!
!          Only works for RDim >= CDim, with AMat = UMat * DVec * VMat
!          If RDim < CDim, you can first calculate A^+ = UMat * DVec * VMat  --> A = (UMat * DVec * VMat)^+
!
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!    UDVdcmpsZ_SVD --> Subroutine to calculate the UDV decomposition for complex RDim-by-CDim matrix;
!    UDVdcmpsR_SVD --> Subroutine to calculate the UDV decomposition for real    RDim-by-CDim matrix.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine UDVdcmpsZ_SVD(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  UDVdcmpsZ_SVD(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the UDV decomposition of RDim*CDim complex matrix of A by 
!                  Singular-Value Decomposition.
! KEYWORDS: Calculate the SVD decomposition, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the SVD decomposition, complex version.
!     By the method of Singular-Value Decomposition (SVD) from Lapack.
!
!     Only works for RDim >= CDim, with AMat = UMat * DVec * VMat
!     If RDim < CDim, you can first calculate A^+ = UMat * DVec * VMat  --> A = (UMat * DVec * VMat)^+
!
!     Input: RDim       --> Dimension of AMat matrix;
!            CDim       --> Dimension of AMat matrix;
!            AMat       --> Input real RDim*CDim matrix (N>M);
!            LDA        --> Leading dimension of AMat matrix;
!            LDU        --> Leading dimension of UMat matrix;
!            LDV        --> Leading dimension of VMat matrix;
!            IfUDVCheck --> Whether to check the UDV decomposition.
!
!     Outpt: UMat       --> The U matrix;             --> By output, UMat is a complex unitary matrix;
!            DVec       --> The D diagonal matrix;    --> By output, DVec is a real, positive vector;
!            VMat       --> The V matrix.             --> By output, VMat is a complex upper triangular matrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!________________________________________ Modules used in this subroutine _____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!________________________________________ All Input and Output Quantities _____________________________________
!______________________________________________________________________________________________________________
      logical IfUDVCheck
		integer RDim, CDim                           ! Dimension of A matrix A(RDim, CDim)
      integer LDA, LDU, LDV                        ! Leading dimensions
		complex(kind=kind(0.d0)) AMat(LDA, *)        ! Input A(RDim, CDim) matrix, Result: A=UDV
		complex(kind=kind(0.d0)) UMat(LDU, *)        ! Result Output U(RDim, CDim) matrix
		real(kind=kind(0.d0))    DVec(CDim)          ! Result Output D(CDim, CDim) matrix --> Only output the diagonal matrix elements
		complex(kind=kind(0.d0)) VMat(LDV, *)        ! Result Output V(CDim, CDim) matrix --> An upper triangular matrix
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info
      integer lwork
      integer I0, I1, I2
      
      real(kind=kind(0.d0)) XMaxm, XMean
      
      complex(kind=kind(0.d0)) Ztp0
      
      complex(kind=kind(0.d0)), allocatable :: work(:)
      real(kind=kind(0.d0)), allocatable :: rwork(:)
      
      complex(kind=kind(0.d0)), allocatable :: AMatTmp(:, :)
      complex(kind=kind(0.d0)), allocatable :: CheckMt(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!_______________________________________________________________________________________________	  
!____________________________ 0. Initialization and allocation _________________________________
!_______________________________________________________________________________________________
      allocate(AMatTmp(RDim, CDim))
      AMatTmp(1:RDim, 1:CDim) = AMat(1:RDim, 1:CDim)
      
      allocate(rwork(5*CDim))
!_______________________________________________________________________________________________	  
!____________________________ 1. Query optimal amount of memory ________________________________
!_______________________________________________________________________________________________
      allocate(work(10))
      info = 0
      call ZGESVD('S', 'S', RDim, CDim, AMatTmp, RDim, DVec, UMat, LDU, VMat, LDV, work, -1, rwork, info)
      lwork = int(dble(work(1)))
      if(allocated(work)) deallocate(work)
!_______________________________________________________________________________________________	  
!____________________________ 2. Perform the SVD decomposition _________________________________
!_______________________________________________________________________________________________
      allocate(work(lwork))
      call ZGESVD('S', 'S', RDim, CDim, AMatTmp, RDim, DVec, UMat, LDU, VMat, LDV, work, lwork, rwork, info)
      if(allocated(work)) deallocate(work)
      
      if( info /= 0 ) then
         write(*, "('UDVdcmpsZ_SVD: error in lapack subroutine ZGESVD! ierror = ', I4)") info
         stop
      end if
!_______________________________________________________________________________________________	  
!____________________________ 3. Check the result of AMat = UMat * DVec * VMat _________________
!_______________________________________________________________________________________________
      if(IfUDVCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = dcmplx(0.0d0, 0.0d0)
   !$OMP PARALLEL &
   !$OMP PRIVATE(I0, I1, I2, Ztp0)
   !$OMP DO
         do I1 = 1, RDim
            do I2 = 1, CDim
               Ztp0 = dcmplx(0.0d0, 0.0d0)
               do I0 = 1, CDim
                  Ztp0 = Ztp0 + UMat(I1, I0) * dcmplx(DVec(I0), 0.0d0) * VMat(I0, I2)
               enddo
               CheckMt(I1, I2) = Ztp0
            enddo
         enddo
   !$OMP END DO
   !$OMP END PARALLEL
         XMaxm = 0.0d0
			XMean = 0.0d0
			call MatrCmpr_C(RDim, CDim, CheckMt(1, 1), RDim, AMat(1, 1), LDA, XMaxm, XMean)
			write(*, "('UDVdcmpsZ_SVD: UDV Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(rwork  )) deallocate(rwork  )
      if(allocated( work  )) deallocate( work  )
				
   end subroutine UDVdcmpsZ_SVD
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine UDVdcmpsR_SVD(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  UDVdcmpsR_SVD(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the UDV decomposition of RDim*CDim real matrix of A by 
!                  Singular-Value Decomposition.
! KEYWORDS: Calculate the UDV decomposition, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the UDV decomposition, real version.
!     By the method of Singular-Value Decomposition (SVD) from Lapack.
!
!     Only works for RDim >= CDim, with AMat = UMat * DVec * VMat
!     If RDim < CDim, you can first calculate A^+ = UMat * DVec * VMat  --> A = (UMat * DVec * VMat)^+
!
!     Input: RDim       --> Dimension of AMat matrix;
!            CDim       --> Dimension of AMat matrix;
!            AMat       --> Input real RDim*CDim matrix (N>M);
!            LDA        --> Leading dimension of AMat matrix;
!            LDU        --> Leading dimension of UMat matrix;
!            LDV        --> Leading dimension of VMat matrix;
!            IfUDVCheck --> Whether to check the UDV decomposition.
!
!     Outpt: UMat       --> The U matrix;             --> By output, UMat is a complex unitary matrix;
!            DVec       --> The D diagonal matrix;    --> By output, DVec is a real, positive vector;
!            VMat       --> The V matrix.             --> By output, VMat is a complex upper triangular matrix;
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!________________________________________ Modules used in this subroutine _____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!________________________________________ All Input and Output Quantities _____________________________________
!______________________________________________________________________________________________________________
      logical IfUDVCheck
		integer RDim, CDim                           ! Dimension of A matrix A(RDim, CDim)
      integer LDA, LDU, LDV                        ! Leading dimensions
		real(kind=kind(0.d0)) AMat(LDA, *)           ! Input A(RDim, CDim) matrix, Result: A=UDV
		real(kind=kind(0.d0)) UMat(LDU, *)           ! Result Output U(RDim, CDim) matrix
		real(kind=kind(0.d0)) DVec(CDim)             ! Result Output D(CDim, CDim) matrix --> Only output the diagonal matrix elements
		real(kind=kind(0.d0)) VMat(LDV, *)           ! Result Output V(CDim, CDim) matrix --> An upper triangular matrix
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info
      integer lwork
      integer I0, I1, I2
      
      real(kind=kind(0.d0)) XMaxm
		real(kind=kind(0.d0)) XMean
      real(kind=kind(0.d0)) Rtp0
      
      real(kind=kind(0.d0)), allocatable :: work(:)
      real(kind=kind(0.d0)), allocatable :: rwork(:)
      
      real(kind=kind(0.d0)), allocatable :: AMatTmp(:, :)
      real(kind=kind(0.d0)), allocatable :: CheckMt(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!_______________________________________________________________________________________________	  
!____________________________ 0. Initialization and allocation _________________________________
!_______________________________________________________________________________________________
      allocate(AMatTmp(RDim, CDim))
      AMatTmp(1:RDim, 1:CDim) = AMat(1:RDim, 1:CDim)
      
      allocate(rwork(5*CDim))
!_______________________________________________________________________________________________	  
!____________________________ 1. Query optimal amount of memory ________________________________
!_______________________________________________________________________________________________
      allocate(work(10))
      info = 0
      call DGESVD('S', 'S', RDim, CDim, AMatTmp, RDim, DVec, UMat, LDU, VMat, LDV, work, -1, rwork, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!_______________________________________________________________________________________________	  
!____________________________ 2. Perform the SVD decomposition _________________________________
!_______________________________________________________________________________________________
      allocate(work(lwork))
      call DGESVD('S', 'S', RDim, CDim, AMatTmp, RDim, DVec, UMat, LDU, VMat, LDV, work, lwork, rwork, info)
      if(allocated( work  )) deallocate( work  )
      
      if( info /= 0 ) then
         write(*, "('UDVdcmpsR_SVD: error in lapack subroutine DGESVD! ierror = ', I4)") info
         stop
      end if
!_______________________________________________________________________________________________	  
!____________________________ 3. Check the result of AMat = UMat * DVec * VMat _________________
!_______________________________________________________________________________________________
      if(IfUDVCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = 0.0d0
   !$OMP PARALLEL &
   !$OMP PRIVATE(I0, I1, I2, Rtp0)
   !$OMP DO
         do I1 = 1, RDim
            do I2 = 1, CDim
               Rtp0 = 0.0d0
               do I0 = 1, CDim
                  Rtp0 = Rtp0 + UMat(I1, I0) * DVec(I0) * VMat(I0, I2)
               enddo
               CheckMt(I1, I2) = Rtp0
            enddo
         enddo
   !$OMP END DO
   !$OMP END PARALLEL
         XMaxm = 0.0d0
			XMean = 0.0d0
			call MatrCmpr_R(RDim, CDim, CheckMt(1, 1), RDim, AMat(1, 1), LDA, XMaxm, XMean)
			write(*, "('UDVdcmpsR_SVD: UDV Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(rwork  )) deallocate(rwork  )
      if(allocated( work  )) deallocate( work  )
				
   end subroutine UDVdcmpsR_SVD
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   