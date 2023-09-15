!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to computes the UDV decomposition of general RDim-by-CDim matrix, with U(RDim, CDim), 
!              D(CDim, CDim) and V(CDim, CDim). D(CDim, CDim) is REAL diagonal matrix.
!          These subroutines calculates the UDV decomposition by the QR decomposition method. 
!              With U(RDim, CDim) as part of the unitary matrix, D(CDim, CDim) as real-nonnegative diagonal matrix,
!              and V(CDim, CDim) as upper-triangular matrix.
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
!******************* UDV decomposition as AMat = UMat * DVec * VMat ***********************************
!    UDVdcmpsZ_QR      --> Obtain UDV decomp for complex RDim-by-CDim matrix, by QR algorithm;
!    UDVdcmpsZ_QRPivot --> Obtain UDV decomp for complex RDim-by-CDim matrix, by QR algorithm with column pivoting;
!
!    UDVdcmpsR_QR      --> Obtain UDV decomp for real    RDim-by-CDim matrix, by QR algorithm;
!    UDVdcmpsR_QRPivot --> Obtain UDV decomp for real    RDim-by-CDim matrix, by QR algorithm with column pivoting;
!
!******************* QR decomposition as AMat = QMat * RMat *******************************************
!    QRdecmpsZ_QR      --> Obtain QR decomp for complex RDim-by-CDim matrix, by QR algorithm;
!    QRdecmpsZ_QRPivot --> Obtain QR decomp for complex RDim-by-CDim matrix, by QR algorithm with column pivoting;
!
!    QRdecmpsR_QR      --> Obtain QR decomp for real    RDim-by-CDim matrix, by QR algorithm;
!    QRdecmpsR_QRPivot --> Obtain QR decomp for real    RDim-by-CDim matrix, by QR algorithm with column pivoting;
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   
!########################################################################################################################
!########################################################################################################################
!####################################### UDV decomposition --> For Complex Version ######################################
!####################################### UDV decomposition --> For Complex Version ######################################
!####################################### UDV decomposition --> For Complex Version ######################################
!########################################################################################################################
!########################################################################################################################
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine UDVdcmpsZ_QR(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  UDVdcmpsZ_QR(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the UDV decomposition of N*M real matrix of A by QR decomposition.
! KEYWORDS: Calculate the UDV decomposition, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the UDV decomposition, complex version.
!     By the method of QR decomposition from Lapack.
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
!            VMat       --> The V matrix.             --> By output, VMat is a complex upper triangular matrix.
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
		integer RDim, CDim                           ! Dimension of A matrix A(RDim, CDim) for UDV decomp
      integer LDA, LDU, LDV
		complex(kind=kind(0.d0)) AMat(LDA, *)        ! Input A(RDim, CDim) matrix, Result: A=UDV
		complex(kind=kind(0.d0)) UMat(LDU, *)        ! Result Output U(RDim, CDim) matrix of UDV
		real(kind=kind(0.d0))    DVec(CDim)          ! Result Output D(CDim, CDim) matrix of UDV
		complex(kind=kind(0.d0)) VMat(LDV, *)        ! Result Output V(CDim, CDim) matrix of UDV
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info, lwork
      integer I0, I1, I2
      integer IVPT(CDim)
      real(kind=kind(0.d0)) XMaxm, XMean
      real(kind=kind(0.d0)), external :: dznrm2
      complex(kind=kind(0.d0)) Ztp0
      complex(kind=kind(0.d0)), allocatable :: work(:)
      complex(kind=kind(0.d0)), allocatable ::  tau(:)
      complex(kind=kind(0.d0)), allocatable :: AMatTmp(:, :)
      complex(kind=kind(0.d0)), allocatable :: VMatTmp(:, :)
      complex(kind=kind(0.d0)), allocatable :: CheckMt(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!_______________________________________________________________________________________________	  
!____________________________ 0. Initialization and allocation _________________________________
!_______________________________________________________________________________________________
      allocate(AMatTmp(RDim, CDim))
      AMatTmp(1:RDim, 1:CDim) = AMat(1:RDim, 1:CDim)
      
      allocate(VMatTmp(RDim, CDim))
      VMatTmp = cmplx(0.0d0, 0.0d0)
      
      allocate(tau(CDim))
      tau = dcmplx(0.d0, 0.d0)
!_______________________________________________________________________________________________	  
!____________________________ 1. Perform the QR decomposition to obtain UDV ____________________
!_______________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________________ (0) Calculate AMatTmp = Q * R ___________________________
!____________________________________________________________________________________  
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0
      call zgeqrf(RDim, CDim, AMatTmp, RDim, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Perform QR decomposition for AMatTmp ____________
!________________________________________________________________________
      allocate(work(lwork))
      call zgeqrf(RDim, CDim, AMatTmp, RDim, tau, work, lwork, info)
      deallocate(work)
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the upper triangular R matrix ________________
!____________________________________________________________________________________
      VMat(1:CDim, 1:CDim) = dcmplx(0.0d0, 0.0d0)
      call zlacpy("U", CDim, CDim, AMatTmp, RDim, VMat, LDV)
!____________________________________________________________________________________	  
!__________________________ (2) Obtain the U matrix _________________________________
!____________________________________________________________________________________ 
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0
      call zungqr(RDim, CDim, CDim, AMatTmp, RDim, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Use zungqr to obtain Q matrix (UMat) ____________
!________________________________________________________________________
      allocate(work(lwork))
      call zungqr(RDim, CDim, CDim, AMatTmp, RDim, tau, work, lwork, info)
      deallocate(work)
      
      UMat(1:RDim, 1:CDim) = dcmplx(0.0d0, 0.0d0)
      call zlacpy("A", RDim, CDim, AMatTmp, RDim, UMat, LDU)
!____________________________________________________________________________________	  
!__________________________ (3) Separate DVec and Rescale R to obtain VMat __________
!____________________________________________________________________________________  
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, CDim, +1
         !!! DVec(I1) = dabs(dble(VMat(I1, I1)))               ! Diagonal element of R matrix
         DVec(I1) = dznrm2(CDim-I1+1, VMat(I1, I1), CDim)      ! Norm of present row of R matrix
         do I2 = I1, CDim, +1
            VMat(I1, I2) = VMat(I1, I2) / dcmplx(DVec(I1), 0.0d0)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&	  
!____________________________ 2. Sort the DVec applying descending order _______________________
!_____________________ AMat = U * D * V = (U*P) * (P^T*D*P) * (P^T*V) __________________________
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
!____________________________________________________________________________________	  
!__________________________ (0) First sort DVec with descending order _______________
!____________________________________________________________________________________
      do I1 = 1, CDim, +1
         IVPT(I1) = I1
      enddo
      call QuckSortRI_Descend(1, CDim, DVec, IVPT) 
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the correct UMat and VMat ____________________
!____________________________________________________________________________________
      do I2 = 1, CDim, +1
         do I1 = 1, RDim, +1
            AMatTmp(I1, I2) = UMat(I1, IVPT(I2))
         enddo
         do I1 = 1, CDim, +1
            VMatTmp(I1, I2) = VMat(IVPT(I1), I2)
         enddo
      enddo
      UMat(1:RDim, 1:CDim) = AMatTmp(1:RDim, 1:CDim)
      VMat(1:CDim, 1:CDim) = VMatTmp(1:CDim, 1:CDim) 
!_______________________________________________________________________________________________	  
!____________________________ 3. Check the result of AMat = UMat * DVec * VMat _________________
!_______________________________________________________________________________________________
      if(IfUDVCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = dcmplx(0.0d0, 0.0d0)
   !$OMP PARALLEL &
   !$OMP PRIVATE(I0, I1, I2, Ztp0)
   !$OMP DO
         do I1 = 1, RDim, +1
            do I2 = 1, CDim, +1
               Ztp0 = dcmplx(0.0d0, 0.0d0)
               do I0 = 1, CDim, +1
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
			write(*, "('UDVdcmpsQR_Z: UDV Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(work   )) deallocate(work   )
      if(allocated(tau    )) deallocate(tau    )
      if(allocated(CheckMt)) deallocate(CheckMt)
		if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(VMatTmp)) deallocate(VMatTmp)
				
   end subroutine UDVdcmpsZ_QR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine UDVdcmpsZ_QRPivot(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  UDVdcmpsZ_QRPivot(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the UDV decomposition of N*M real matrix of A by QR decomposition.
! KEYWORDS: Calculate the UDV decomposition, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the UDV decomposition, complex version.
!     By the method of QR decomposition from Lapack with artificial column pivoting.
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
!            VMat       --> The V matrix.             --> By output, VMat is a complex matrix;
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
		real(kind=kind(0.d0))    DVec(CDim)          ! Result Output D(CDim, CDim) matrix
		complex(kind=kind(0.d0)) VMat(LDV, *)        ! Result Output V(CDim, CDim) matrix
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info, lwork
      integer I0, I1, I2, IMaxm
      integer IVPT(CDim), IVPTM1(CDim)
      real(kind=kind(0.d0)) XMaxm, XMean
      real(kind=kind(0.d0)), external :: dznrm2
      complex(kind=kind(0.d0)) Ztp0
      real(kind=kind(0.d0)), allocatable :: Rwork(:)
      complex(kind=kind(0.d0)), allocatable :: work(:)
      complex(kind=kind(0.d0)), allocatable ::  tau(:)
      complex(kind=kind(0.d0)), allocatable :: AMatTmp(:, :)
      complex(kind=kind(0.d0)), allocatable :: VMatTmp(:, :)
      complex(kind=kind(0.d0)), allocatable :: TmpMat1(:, :)
      complex(kind=kind(0.d0)), allocatable :: TmpMat2(:, :)
      complex(kind=kind(0.d0)), allocatable :: CheckMt(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!_______________________________________________________________________________________________	  
!____________________________ 0. Initialization and allocation _________________________________
!_______________________________________________________________________________________________
      allocate(AMatTmp(RDim, CDim))
      AMatTmp(1:RDim, 1:CDim) = AMat(1:RDim, 1:CDim)
      
      allocate(VMatTmp(CDim, CDim))
      VMatTmp = dcmplx(0.0d0, 0.0d0)
      
      allocate(TmpMat1(CDim, CDim))
      allocate(TmpMat2(RDim, CDim))
      TmpMat1 = dcmplx(0.0d0, 0.0d0)
      TmpMat2 = dcmplx(0.0d0, 0.0d0)
      
      allocate(Rwork(2*CDim))
      Rwork = 0.0d0
      
      allocate(tau(CDim))
      tau = dcmplx(0.0d0, 0.0d0)
!_______________________________________________________________________________________________	  
!____________________________ 1. Perform the column-pivoted QR decomps for AMatTmp _____________
!_______________________________ AMatTmp * P = UMat * VMat _____________________________________
!_______________________________________________________________________________________________ 
!____________________________________________________________________________________	  
!______________________ (0) AMatTmp * P = UMat * VMat _______________________________
!____________________________________________________________________________________
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0; IVPT = 0
      call zgeqp3(RDim, CDim, AMatTmp, RDim, IVPT, tau, work, -1, Rwork, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Perform colum-pivoted QR for AMatTmp ____________
!________________________________________________________________________
      allocate(work(lwork))
      IVPT = 0
      call zgeqp3(RDim, CDim, AMatTmp, RDim, IVPT, tau, work, lwork, Rwork, info)
      deallocate(work) 
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the upper triangular R matrix ________________
!____________________________________________________________________________________
      VMatTmp = dcmplx(0.d0, 0.d0)
      call zlacpy("U", CDim, CDim, AMatTmp, RDim, VMatTmp, CDim)
!____________________________________________________________________________________	  
!__________________________ (2) R * P^T = D * R' ____________________________________
!__________________________ R * P^T = D * (D^{-1}*R) * P^T __________________________
!__________________________ R' = (D^{-1}*R) * P^T ___________________________________
!____________________________________________________________________________________
      do I1 = 1, CDim, +1
         IVPTM1(IVPT(I1)) = I1
         !!! DVec(I1) = dabs(dble(VMatTmp(I1, I1)))               ! Diagonal element of R matrix
         DVec(I1) = dznrm2(CDim-I1+1, VMatTmp(I1, I1), CDim)      ! Norm of present row of R matrix
      enddo
            
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, CDim, +1
         do I2 = 1, CDim, +1
            TmpMat1(I1, I2) = VMatTmp(I1, IVPTM1(I2)) / DVec(I1)
         enddo
      enddo 
   !$OMP END DO
   !$OMP END PARALLEL
!____________________________________________________________________________________	  
!__________________________ (3) Obtain the U matrix _________________________________
!____________________________________________________________________________________ 
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0
      call zungqr(RDim, CDim, CDim, AMatTmp, RDim, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Use zungqr to obtain Q matrix (UMat) ____________
!________________________________________________________________________
      allocate(work(lwork))
      call zungqr(RDim, CDim, CDim, AMatTmp, RDim, tau, work, lwork, info)
      deallocate(work)
      deallocate(tau)
      
      TmpMat2 = dcmplx(0.0d0, 0.0d0)
      call zlacpy("A", RDim, CDim, AMatTmp, RDim, TmpMat2, RDim)
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&	  
!____________________________ 2. Sort the DVec applying descending order _______________________
!_____________________ AMat = U * D * V = (U*P) * (P^T*D*P) * (P^T*V) __________________________
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
!____________________________________________________________________________________	  
!__________________________ (0) First sort DVec with descending order _______________
!____________________________________________________________________________________
      do I1 = 1, CDim, +1
         IVPT(I1) = I1
      enddo
      call QuckSortRI_Descend(1, CDim, DVec, IVPT) 
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the correct UMat and VMat ____________________
!____________________________________________________________________________________
      do I2 = 1, CDim, +1
         do I1 = 1, RDim, +1
            UMat(I1, I2) = TmpMat2(I1, IVPT(I2))
         enddo
         do I1 = 1, CDim, +1
            VMat(I1, I2) = TmpMat1(IVPT(I1), I2)
         enddo
      enddo
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
			write(*, "('UDVdcmpsZ_QRPivot: UDV Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
      if(allocated(work   )) deallocate(work   )
      if(allocated(Rwork  )) deallocate(Rwork  )
      if(allocated(tau    )) deallocate(tau    )
      if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(VMatTmp)) deallocate(VMatTmp)
      if(allocated(TmpMat1)) deallocate(TmpMat1)
      if(allocated(TmpMat2)) deallocate(TmpMat2)
      if(allocated(CheckMt)) deallocate(CheckMt)
				
   end subroutine UDVdcmpsZ_QRPivot
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
!########################################################################################################################
!########################################################################################################################
!####################################### UDV decomposition --> For real Version #########################################
!####################################### UDV decomposition --> For real Version #########################################
!####################################### UDV decomposition --> For real Version #########################################
!########################################################################################################################
!########################################################################################################################


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine UDVdcmpsR_QR(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  UDVdcmpsR_QR(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the UDV decomposition of N*M real matrix of A by QR decomposition.
! KEYWORDS: Calculate the UDV decomposition, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the UDV decomposition, real version.
!     By the method of QR decomposition from Lapack.
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
!            VMat       --> The V matrix.             --> By output, VMat is a complex upper triangular matrix.
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
		integer RDim, CDim                        ! Dimension of A matrix A(RDim, CDim) for UDV decomp
      integer LDA, LDU, LDV
		real(kind=kind(0.d0)) AMat(LDA, *)        ! Input A(RDim, CDim) matrix, Result: A=UDV
		real(kind=kind(0.d0)) UMat(LDU, *)        ! Result Output U(RDim, CDim) matrix of UDV
		real(kind=kind(0.d0)) DVec(CDim)          ! Result Output D(CDim, CDim) matrix of UDV
		real(kind=kind(0.d0)) VMat(LDV, *)        ! Result Output V(CDim, CDim) matrix of UDV
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info, lwork
      integer I0, I1, I2
      integer IVPT(CDim)
      real(kind=kind(0.d0)) XMaxm, XMean, Rtp0
      real(kind=kind(0.d0)), external :: dnrm2
      real(kind=kind(0.d0)), allocatable :: work(:)
      real(kind=kind(0.d0)), allocatable ::  tau(:)
      real(kind=kind(0.d0)), allocatable :: AMatTmp(:, :)
      real(kind=kind(0.d0)), allocatable :: VMatTmp(:, :)
      real(kind=kind(0.d0)), allocatable :: CheckMt(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!_______________________________________________________________________________________________	  
!____________________________ 0. Initialization and allocation _________________________________
!_______________________________________________________________________________________________
      allocate(AMatTmp(RDim, CDim))
      AMatTmp(1:RDim, 1:CDim) = AMat(1:RDim, 1:CDim)
      
      allocate(VMatTmp(CDim, CDim))
      VMatTmp = 0.0d0
   
      allocate(tau(CDim))
      tau = 0.0d0
!_______________________________________________________________________________________________	  
!____________________________ 1. Perform the QR decomposition to obtain UDV ____________________
!_______________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________________ (0) Calculate AMatTmp = Q * R ___________________________
!____________________________________________________________________________________  
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0
      call dgeqrf(RDim, CDim, AMatTmp, RDim, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Perform QR decomposition for AMatTmp ____________
!________________________________________________________________________
      allocate(work(lwork))
      call dgeqrf(RDim, CDim, AMatTmp, RDim, tau, work, lwork, info)
      deallocate(work)
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the upper triangular R matrix ________________
!____________________________________________________________________________________
      VMat(1:CDim, 1:CDim) = 0.0d0
      call dlacpy("U", CDim, CDim, AMatTmp, RDim, VMat, LDV)
!____________________________________________________________________________________	  
!__________________________ (2) Obtain the U matrix _________________________________
!____________________________________________________________________________________ 
      allocate(work(CDim))
      call dorg2r(RDim, CDim, CDim, AMatTmp, RDim, tau, work, info)
      deallocate(work)
      deallocate(tau )
      
      UMat(1:RDim, 1:CDim) = 0.0d0
      call dlacpy("A", RDim, CDim, AMatTmp, RDim, UMat, LDU)
!____________________________________________________________________________________	  
!__________________________ (3) Separate DVec and Rescale R to obtain VMat __________
!____________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, CDim, +1
         !!! DVec(I1) = dabs(VMat(I1, I1))                   ! Diagonal element
         DVec(I1) = dnrm2(CDim-I1+1, VMat(I1, I1), CDim)     ! Norm of the present row
         do I2 = I1, CDim, +1
            VMat(I1, I2) = VMat(I1, I2) / DVec(I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&	  
!____________________________ 4. Sort the DVec applying descending order _______________________
!_____________________ AMat = U * D * V = (U*P) * (P^T*D*P) * (P^T*V) __________________________
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
!____________________________________________________________________________________	  
!__________________________ (0) First sort DVec with descending order _______________
!____________________________________________________________________________________
      do I1 = 1, CDim, +1
         IVPT(I1) = I1
      enddo
      call QuckSortRI_Descend(1, CDim, DVec, IVPT) 
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the correct UMat and VMat ____________________
!____________________________________________________________________________________
      do I2 = 1, CDim, +1
         do I1 = 1, RDim, +1
            AMatTmp(I1, I2) = UMat(I1, IVPT(I2))
         enddo
         do I1 = 1, CDim, +1
            VMatTmp(I1, I2) = VMat(IVPT(I1), I2)
         enddo
      enddo
      UMat(1:RDim, 1:CDim) = AMatTmp(1:RDim, 1:CDim)
      VMat(1:CDim, 1:CDim) = VMatTmp(1:CDim, 1:CDim) 
!_______________________________________________________________________________________________	  
!____________________________ 2. Check the result of AMat = UMat * DVec * VMat _________________
!_______________________________________________________________________________________________
      if(IfUDVCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = 0.0d0
   !$OMP PARALLEL &
   !$OMP PRIVATE(I0, I1, I2, Rtp0)
   !$OMP DO
         do I1 = 1, RDim, +1
            do I2 = 1, CDim, +1
               Rtp0 = 0.0d0
               do I0 = 1, CDim, +1
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
			write(*, "('UDVdcmpsR_QR: UDV Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
      if(allocated(work   )) deallocate(work   )
      if(allocated(tau    )) deallocate(tau    )
      if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(VMatTmp)) deallocate(VMatTmp)
      if(allocated(CheckMt)) deallocate(CheckMt)
				
   end subroutine UDVdcmpsR_QR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine UDVdcmpsR_QRPivot(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  UDVdcmpsR_QRPivot(RDim, CDim, AMat, LDA, UMat, LDU, DVec, VMat, LDV, IfUDVCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the UDV decomposition of N*M real matrix of A by QR decomposition.
! KEYWORDS: Calculate the UDV decomposition, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the UDV decomposition, real version.
!     By the method of QR decomposition from Lapack with artificial column pivoting.
!
!     The column-pivoted QR decomposition is done by the DGEQP3 subroutine from MKL.
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
!     Outpt: UMat       --> The U matrix;             --> By output, UMat is a real unitary matrix;
!            DVec       --> The D diagonal matrix;    --> By output, DVec is a real, positive vector;
!            VMat       --> The V matrix.             --> By output, VMat is a real matrix.
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
		integer RDim, CDim                        ! Dimension of A matrix A(RDim, CDim) for UDV decomp
      integer LDA, LDU, LDV
		real(kind=kind(0.d0)) AMat(LDA, *)        ! Input A(RDim, CDim) matrix, Result: A=UDV
		real(kind=kind(0.d0)) UMat(LDU, *)        ! Result Output U(RDim, CDim) matrix of UDV
		real(kind=kind(0.d0)) DVec(CDim)          ! Result Output D(CDim, CDim) matrix of UDV
		real(kind=kind(0.d0)) VMat(LDV, *)        ! Result Output V(CDim, CDim) matrix of UDV
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info, lwork
      integer I0, I1, I2
      integer IVPT(CDim), IVPTM1(CDim)
      real(kind=kind(0.d0)) XMaxm, XMean, Rtp0
      real(kind=kind(0.d0)), external :: dnrm2
      real(kind=kind(0.d0)), allocatable :: work(:)
      real(kind=kind(0.d0)), allocatable ::  tau(:)
      real(kind=kind(0.d0)), allocatable :: AMatTmp(:, :)
      real(kind=kind(0.d0)), allocatable :: VMatTmp(:, :)
      real(kind=kind(0.d0)), allocatable :: TmpMat1(:, :)
      real(kind=kind(0.d0)), allocatable :: TmpMat2(:, :)
      real(kind=kind(0.d0)), allocatable :: CheckMt(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!_______________________________________________________________________________________________	  
!____________________________ 0. Initialization and allocation _________________________________
!_______________________________________________________________________________________________
      allocate(AMatTmp(RDim, CDim))
      AMatTmp(1:RDim, 1:CDim) = AMat(1:RDim, 1:CDim)
      
      allocate(VMatTmp(CDim, CDim))
      VMatTmp = 0.0d0
      
      allocate(TmpMat1(CDim, CDim))
      allocate(TmpMat2(RDim, CDim))
      TmpMat1 = 0.0d0
      TmpMat2 = 0.0d0
      
      allocate(tau(CDim))
      tau = 0.0d0
!_______________________________________________________________________________________________	  
!____________________________ 1. Perform the column-pivoted QR decomps for AMatTmp _____________
!_______________________________ AMatTmp * P = UMat * VMat _____________________________________
!_______________________________________________________________________________________________ 
!____________________________________________________________________________________	  
!______________________ (0) AMatTmp * P = UMat * VMat _______________________________
!____________________________________________________________________________________
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0; IVPT = 0
      call dgeqp3(RDim, CDim, AMatTmp, RDim, IVPT, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Perform colum-pivoted QR for AMatTmp ____________
!________________________________________________________________________
      allocate(work(lwork))
      IVPT = 0
      call dgeqp3(RDim, CDim, AMatTmp, RDim, IVPT, tau, work, lwork, info)
      deallocate(work) 
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the upper triangular R matrix ________________
!____________________________________________________________________________________
      VMatTmp = 0.0d0
      call dlacpy("U", CDim, CDim, AMatTmp, RDim, VMatTmp, CDim)
!____________________________________________________________________________________	  
!__________________________ (2) R * P^T = D * R' ____________________________________
!__________________________ R * P^T = D * (D^{-1}*R) * P^T __________________________
!__________________________ R' = (D^{-1}*R) * P^T ___________________________________
!____________________________________________________________________________________
      do I1 = 1, CDim, +1
         IVPTM1(IVPT(I1)) = I1
         !!! DVec(I1) = dabs(VMatTmp(I1, I1))                   ! Diagonal element
         DVec(I1) = dnrm2(CDim-I1+1, VMatTmp(I1, I1), CDim)     ! Norm of the present row
      enddo
            
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, CDim, +1
         do I2 = 1, CDim, +1
            TmpMat1(I1, I2) = VMatTmp(I1, IVPTM1(I2)) / DVec(I1)
         enddo
      enddo 
   !$OMP END DO
   !$OMP END PARALLEL
!____________________________________________________________________________________	  
!__________________________ (3) Obtain the U matrix _________________________________
!____________________________________________________________________________________ 
      allocate(work(CDim))
      call dorg2r(RDim, CDim, CDim, AMatTmp, RDim, tau, work, info)
      deallocate(work)
      deallocate(tau)
      
      TmpMat2 = 0.0d0
      call dlacpy("A", RDim, CDim, AMatTmp, RDim, TmpMat2, RDim)
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&	  
!____________________________ 2. Sort the DVec applying descending order _______________________
!_____________________ AMat = U * D * V = (U*P) * (P^T*D*P) * (P^T*V) __________________________
!&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
!____________________________________________________________________________________	  
!__________________________ (0) First sort DVec with descending order _______________
!____________________________________________________________________________________
      do I1 = 1, CDim, +1
         IVPT(I1) = I1
      enddo
      call QuckSortRI_Descend(1, CDim, DVec, IVPT) 
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the correct UMat and VMat ____________________
!____________________________________________________________________________________
      do I2 = 1, CDim, +1
         do I1 = 1, RDim, +1
            UMat(I1, I2) = TmpMat2(I1, IVPT(I2))
         enddo
         do I1 = 1, CDim, +1
            VMat(I1, I2) = TmpMat1(IVPT(I1), I2)
         enddo
      enddo
!_______________________________________________________________________________________________	  
!____________________________ 3. Check the result of AMat = UMat * DVec * VMat _________________
!_______________________________________________________________________________________________
      if(IfUDVCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = 0.0d0
   !$OMP PARALLEL &
   !$OMP PRIVATE(I0, I1, I2, Rtp0)
   !$OMP DO
         do I1 = 1, RDim, +1
            do I2 = 1, CDim, +1
               Rtp0 = 0.0d0
               do I0 = 1, CDim, +1
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
			write(*, "('UDVdcmpsR_QRPivot: UDV Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
      if(allocated(work   )) deallocate(work   )
      if(allocated(tau    )) deallocate(tau    )
      if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(VMatTmp)) deallocate(VMatTmp)
      if(allocated(TmpMat1)) deallocate(TmpMat1)
      if(allocated(TmpMat2)) deallocate(TmpMat2)
      if(allocated(CheckMt)) deallocate(CheckMt)
				
   end subroutine UDVdcmpsR_QRPivot
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
!########################################################################################################################
!########################################################################################################################
!####################################### QR decomposition --> For Complex Version #######################################
!####################################### QR decomposition --> For Complex Version #######################################
!####################################### QR decomposition --> For Complex Version #######################################
!########################################################################################################################
!########################################################################################################################
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine QRdecmpsZ_QR(RDim, CDim, AMat, LDA, QMat, LDQ, RMat, LDR, IfQRCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  QRdecmpsZ_QR(RDim, CDim, AMat, LDA, QMat, LDQ, RMat, LDR, IfQRCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create AMat(RDim, CDim) = QMat(RDim, CDim) * RMat(CDim, CDim) by 
!                 simple QR decomposition.
! KEYWORDS: Calculate the column-pivoted QR decomposition, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the QR decomposition, complex version.
!     By the method of QR decomposition from Lapack.
!
!     Only works for RDim >= CDim, with AMat = QMat * RMat
!     If RDim < CDim, you can first calculate A^+ = QMat * RMat  --> A = (QMat * RMat)^+
!
!     Input: RDim      --> Dimension of AMat matrix;
!            CDim      --> Dimension of AMat matrix;
!            AMat      --> Input real RDim*CDim matrix (N>M);
!            LDA       --> Leading dimension of AMat matrix;
!            LDQ       --> Leading dimension of QMat matrix;
!            LDR       --> Leading dimension of RMat matrix;
!            IfQRCheck --> Whether to check the QR decomposition.
!
!     Outpt: QMat      --> The Q matrix;  --> By output, QMat is a complex unitary matrix;
!            RMat      --> The R matrix.  --> By output, RMat is a complex upper triangular matrix.
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
      logical IfQRCheck
		integer RDim, CDim                           ! Dimension of A matrix A(RDim, CDim) for QR decomp
      integer LDA, LDQ, LDR
		complex(kind=kind(0.d0)) AMat(LDA, *)        ! Input A(RDim, CDim) matrix, Result: A=QR
		complex(kind=kind(0.d0)) QMat(LDQ, *)        ! Result Output Q(RDim, CDim) matrix of QR
		complex(kind=kind(0.d0)) RMat(LDR, *)        ! Result Output R(CDim, CDim) matrix of QR
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info, lwork
      integer I1, I2
      real(kind=kind(0.d0)) XMaxm, XMean
      real(kind=kind(0.d0)) LogdDetR
      complex(kind=kind(0.d0)) LogzDet
      
      complex(kind=kind(0.d0)), allocatable :: work(:)
      complex(kind=kind(0.d0)), allocatable ::  tau(:)
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
      
      allocate(tau(CDim))
      tau = dcmplx(0.d0, 0.d0)
!_______________________________________________________________________________________________	  
!____________________________ 1. Perform the QR decomposition to obtain QR _____________________
!_______________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________________ (0) Calculate AMatTmp = Q * R ___________________________
!____________________________________________________________________________________  
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0
      call zgeqrf(RDim, CDim, AMatTmp, RDim, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Perform QR decomposition for AMatTmp ____________
!________________________________________________________________________
      allocate(work(lwork))
      call zgeqrf(RDim, CDim, AMatTmp, RDim, tau, work, lwork, info)
      deallocate(work)
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the upper triangular R matrix ________________
!______________________________ and the log of its determinant ______________________
!____________________________________________________________________________________
      RMat(1:CDim, 1:CDim) = dcmplx(0.0d0, 0.0d0)
      call zlacpy("U", CDim, CDim, AMatTmp, RDim, RMat, LDR)
      
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, CDim, +1
         LogzDet = LogzDet + log( RMat(I1, I1) )
      enddo
!____________________________________________________________________________________	  
!__________________________ (2) Obtain the Q matrix _________________________________
!____________________________________________________________________________________ 
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0
      call zungqr(RDim, CDim, CDim, AMatTmp, RDim, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Use zungqr to obtain Q matrix (QMat) ____________
!________________________________________________________________________
      allocate(work(lwork))
      call zungqr(RDim, CDim, CDim, AMatTmp, RDim, tau, work, lwork, info)
      deallocate(work)
      
      QMat(1:RDim, 1:CDim) = dcmplx(0.0d0, 0.0d0)
      call zlacpy("A", RDim, CDim, AMatTmp, RDim, QMat, LDQ)
!____________________________________________________________________________________	  
!__________________________ (3) Choose positive dDetR and Rescale Q and R ___________
!____________________________________________________________________________________
      LogdDetR = dreal(LogzDet)      
      if( dcos(dimag(LogzDet)) < 0.0d0 ) then
         do I1 = 1, RDim, +1
            QMat(I1, 1) = - QMat(I1, 1)
         enddo
         do I2 = 1, CDim, +1   
            RMat(1, I2) = - RMat(1, I2)
         enddo
      end if
!_______________________________________________________________________________________________	  
!____________________________ 2. Check the result of AMat = QMat * RMat ________________________
!_______________________________________________________________________________________________
      if(IfQRCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = dcmplx(0.0d0, 0.0d0)
         call ZGEMM("N", "N", RDim, CDim, CDim, dcmplx(1.0d0, 0.0d0), QMat, LDQ, RMat, LDR, dcmplx(0.0d0, 0.0d0), &
            & CheckMt, RDim)
         XMaxm = 0.0d0; XMean = 0.0d0
			call MatrCmpr_C(RDim, CDim, CheckMt(1, 1), RDim, AMat(1, 1), LDA, XMaxm, XMean)
			write(*, "('QRdecmpsZ_QR: QR Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(work   )) deallocate(work   )
      if(allocated(tau    )) deallocate(tau    )
      if(allocated(CheckMt)) deallocate(CheckMt)
				
   end subroutine QRdecmpsZ_QR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine QRdecmpsZ_QRPivot(RDim, CDim, AMat, LDA, QMat, LDQ, RMat, LDR, IfQRCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  QRdecmpsZ_QRPivot(RDim, CDim, AMat, LDA, QMat, LDQ, RMat, LDR, IfQRCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create AMat(RDim, CDim) = QMat(RDim, CDim) * RMat(CDim, CDim) by 
!                 column-pivoted QR decomposition.
! KEYWORDS: Calculate the column-pivoted QR decomposition, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the QR decomposition, complex version.
!     By the method of QR decomposition from Lapack with artificial column pivoting.
!
!     Only works for RDim >= CDim, with AMat = QMat * RMat;
!     If RDim < CDim, you can first calculate A^+ = QMat * RMat  --> A = (QMat * RMat)^+
!
!     Input: RDim      --> Dimension of AMat matrix;
!            CDim      --> Dimension of AMat matrix;
!            AMat      --> Input real RDim*CDim matrix (N>M);
!            LDA       --> Leading dimension of AMat matrix;
!            LDQ       --> Leading dimension of QMat matrix;
!            LDR       --> Leading dimension of RMat matrix;
!            IfQRCheck --> Whether to check the QR decomposition.
!
!     Outpt: QMat      --> The Q matrix;  --> By output, QMat is a complex unitary matrix;
!            RMat      --> The R matrix.  --> By output, RMat is a complex matrix;
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
      logical IfQRCheck
		integer RDim, CDim                           ! Dimension of A matrix A(RDim, CDim)
      integer LDA, LDQ, LDR                        ! Leading dimensions
		complex(kind=kind(0.d0)) AMat(LDA, *)        ! Input A(RDim, CDim) matrix, Result: A=QR
		complex(kind=kind(0.d0)) QMat(LDQ, *)        ! Result Output Q(RDim, CDim) matrix
		complex(kind=kind(0.d0)) RMat(LDR, *)        ! Result Output R(CDim, CDim) matrix
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info, lwork
      integer I1, I2
      integer IVPT(CDim), IVPTM1(CDim)
      real(kind=kind(0.d0)) LogdDetR
      real(kind=kind(0.d0)) XMaxm, XMean
      complex(kind=kind(0.d0)) LogzDet
      real(kind=kind(0.d0)), allocatable :: Rwork(:)
      complex(kind=kind(0.d0)), allocatable :: work(:)
      complex(kind=kind(0.d0)), allocatable ::  tau(:)
      complex(kind=kind(0.d0)), allocatable :: AMatTmp(:, :)
      complex(kind=kind(0.d0)), allocatable :: RMatTmp(:, :)
      complex(kind=kind(0.d0)), allocatable :: CheckMt(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!_______________________________________________________________________________________________	  
!____________________________ 0. Initialization and allocation _________________________________
!_______________________________________________________________________________________________
      allocate(AMatTmp(RDim, CDim))
      AMatTmp(1:RDim, 1:CDim) = AMat(1:RDim, 1:CDim)
      
      allocate(RMatTmp(RDim, CDim))
      RMatTmp = dcmplx(0.0d0, 0.0d0)
      
      allocate(Rwork(2*CDim))
      Rwork = 0.0d0
      
      allocate(tau(CDim))
      tau = dcmplx(0.d0, 0.d0)
!_______________________________________________________________________________________________	  
!____________________________ 1. Perform the column-pivoted QR decomps for AMatTmp _____________
!_______________________________ AMatTmp * P = QMat * RMat _____________________________________
!_______________________________________________________________________________________________ 
!____________________________________________________________________________________	  
!__________________________ (0) AMatTmp * P = QMat * RMat ___________________________
!____________________________________________________________________________________
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0; IVPT = 0
      call zgeqp3(RDim, CDim, AMatTmp, RDim, IVPT, tau, work, -1, Rwork, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Perform colum-pivoted QR for AMatTmp ____________
!________________________________________________________________________
      allocate(work(lwork))
      IVPT = 0
      call zgeqp3(RDim, CDim, AMatTmp, RDim, IVPT, tau, work, lwork, Rwork, info)
      deallocate(work) 
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the RMatTmp matrix ___________________________
!______________________________ and the log of its determinant ______________________
!____________________________________________________________________________________
      RMatTmp = dcmplx(0.0d0, 0.0d0)
      call zlacpy("U", CDim, CDim, AMatTmp, RDim, RMatTmp, CDim)
      
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, CDim, +1
         LogzDet = LogzDet + log( RMatTmp(I1, I1) )
      enddo
!____________________________________________________________________________________	  
!__________________________ (2) Calculate R * P^T ___________________________________
!____________________________________________________________________________________
      do I1 = 1, CDim, +1
         IVPTM1(IVPT(I1)) = I1
      enddo

      do I1 = 1, CDim, +1
         do I2 = 1, CDim, +1
            RMat(I1, I2) = RMatTmp(I1, IVPTM1(I2))
         enddo
      enddo
!____________________________________________________________________________________	  
!__________________________ (3) Obtain the QMat matrix ______________________________
!____________________________________________________________________________________ 
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0
      call zungqr(RDim, CDim, CDim, AMatTmp, RDim, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Use zungqr to obtain Q matrix (QMat) ____________
!________________________________________________________________________
      allocate(work(lwork))
      call zungqr(RDim, CDim, CDim, AMatTmp, RDim, tau, work, lwork, info)
      deallocate(work)
      deallocate(tau)
      
      QMat(1:RDim, 1:CDim) = dcmplx(0.0d0, 0.0d0)
      call zlacpy("A", RDim, CDim, AMatTmp, RDim, QMat, LDQ)
!_______________________________________________________________________________________________	  
!____________________________ 2. Obtain the final QMat, RMat matrices as _______________________
!_______________________________ and get the logdDet of RMat matrix ____________________________
!_______________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________________ (0) Accumulate LogzDet of P^T matrix ____________________
!____________________________________________________________________________________
      I2 = 0
      do I1 = 1, CDim, +1
         if(IVPTM1(I1) /= I1) then
            I2 = I2 + 1
         end if
      enddo
      LogzDet = LogzDet + log( dcmplx((-1.0d0)**(mod(I2, 2)), 0.d0) )
!____________________________________________________________________________________	  
!__________________________ (1) Choose positive dDetR and Rescale Q and R ___________
!____________________________________________________________________________________
      LogdDetR = dreal(LogzDet)      
      if( dcos(dimag(LogzDet)) < 0.0d0 ) then
         do I1 = 1, RDim, +1
            QMat(I1, 1) = - QMat(I1, 1)
         enddo
         do I2 = 1, CDim, +1   
            RMat(1, I2) = - RMat(1, I2)
         enddo
      end if
!_______________________________________________________________________________________________	  
!____________________________ 3. Check the result of AMat = QMat * RMat ________________________
!_______________________________________________________________________________________________
      if(IfQRCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = dcmplx(0.0d0, 0.0d0)
         call ZGEMM("N", "N", RDim, CDim, CDim, dcmplx(1.0d0, 0.0d0), QMat, LDQ, RMat, LDR, dcmplx(0.0d0, 0.0d0), &
            & CheckMt, RDim)
         XMaxm = 0.0d0; XMean = 0.0d0
			call MatrCmpr_C(RDim, CDim, CheckMt(1, 1), RDim, AMat(1, 1), LDA, XMaxm, XMean)
			write(*, "('QRdecmpsZ_QRPivot: QR Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
      if(allocated(work   )) deallocate(work   )
      if(allocated(Rwork  )) deallocate(Rwork  )
      if(allocated(tau    )) deallocate(tau    )
      if(allocated(CheckMt)) deallocate(CheckMt)
		if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(RMatTmp)) deallocate(RMatTmp)
				
   end subroutine QRdecmpsZ_QRPivot
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
!########################################################################################################################
!########################################################################################################################
!####################################### QR decomposition --> For Real Version ##########################################
!####################################### QR decomposition --> For Real Version ##########################################
!####################################### QR decomposition --> For Real Version ##########################################
!########################################################################################################################
!########################################################################################################################
  
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine QRdecmpsR_QR(RDim, CDim, AMat, LDA, QMat, LDQ, RMat, LDR, IfQRCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  QRdecmpsR_QR(RDim, CDim, AMat, LDA, QMat, LDQ, RMat, LDR, IfQRCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the QR decomposition of N*M real matrix of A by QR decomposition.
! KEYWORDS: Calculate the QR decomposition, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the QR decomposition, real version.
!     By the method of QR decomposition from Lapack.
!
!     Only works for RDim >= CDim, with AMat = QMat * RMat
!     If RDim < CDim, you can first calculate A^+ = QMat * RMat  --> A = (QMat * RMat)^+
!
!     Input: RDim       --> Dimension of AMat matrix;
!            CDim       --> Dimension of AMat matrix;
!            AMat       --> Input real RDim*CDim matrix (N>M);
!            LDA        --> Leading dimension of AMat matrix;
!            LDQ        --> Leading dimension of QMat matrix;
!            LDR        --> Leading dimension of RMat matrix;
!            IfQRCheck --> Whether to check the QR decomposition.
!
!     Outpt: QMat       --> The Q matrix;  --> By output, QMat is a real symmetric matrix;
!            RMat       --> The R matrix.  --> By output, RMat is a real upper triangular matrix.
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
      logical IfQRCheck
		integer RDim, CDim                        ! Dimension of A matrix A(RDim, CDim) for QR decomp
      integer LDA, LDQ, LDR
		real(kind=kind(0.d0)) AMat(LDA, *)        ! Input A(RDim, CDim) matrix, Result: A=QR
		real(kind=kind(0.d0)) QMat(LDQ, *)        ! Result Output Q(RDim, CDim) matrix of QR
		real(kind=kind(0.d0)) RMat(LDR, *)        ! Result Output R(CDim, CDim) matrix of QR
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info, lwork
      integer I1, I2
      real(kind=kind(0.d0)) XMaxm, XMean, Rtp0
      real(kind=kind(0.d0)) LogdDetR
      complex(kind=kind(0.d0)) LogzDet
      real(kind=kind(0.d0)), allocatable :: work(:)
      real(kind=kind(0.d0)), allocatable ::  tau(:)
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
   
      allocate(tau(CDim))
      tau = 0.0d0
!_______________________________________________________________________________________________	  
!____________________________ 1. Perform the QR decomposition to obtain QR ____________________
!_______________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________________ (0) Calculate AMatTmp = Q * R ___________________________
!____________________________________________________________________________________  
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0
      call dgeqrf(RDim, CDim, AMatTmp, RDim, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Perform QR decomposition for AMatTmp ____________
!________________________________________________________________________
      allocate(work(lwork))
      call dgeqrf(RDim, CDim, AMatTmp, RDim, tau, work, lwork, info)
      deallocate(work)
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the upper triangular R matrix ________________
!______________________________ and the log of its determinant ______________________
!____________________________________________________________________________________
      RMat(1:CDim, 1:CDim) = 0.0d0
      call dlacpy("U", CDim, CDim, AMatTmp, RDim, RMat, LDR)
      
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, CDim, +1
         LogzDet = LogzDet + log( dcmplx(RMat(I1, I1), 0.0d0) )
      enddo  
!____________________________________________________________________________________	  
!__________________________ (2) Obtain the Q matrix _________________________________
!____________________________________________________________________________________ 
      allocate(work(CDim))
      call dorg2r(RDim, CDim, CDim, AMatTmp, RDim, tau, work, info)
      deallocate(work)
      deallocate(tau )
      
      QMat(1:RDim, 1:CDim) = 0.0d0
      call dlacpy("A", RDim, CDim, AMatTmp, RDim, QMat, LDQ)
!____________________________________________________________________________________	  
!__________________________ (3) Choose positive dDetR and Rescale Q and R ___________
!____________________________________________________________________________________
      LogdDetR = dreal(LogzDet)      
      if( dcos(dimag(LogzDet)) < 0.0d0 ) then
         do I1 = 1, RDim, +1
            QMat(I1, 1) = - QMat(I1, 1)
         enddo
         do I2 = 1, CDim, +1   
            RMat(1, I2) = - RMat(1, I2)
         enddo
      end if
!_______________________________________________________________________________________________	  
!____________________________ 2. Check the result of AMat = QMat * RMat ________________________
!_______________________________________________________________________________________________
      if(IfQRCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = 0.0d0
         call DGEMM("N", "N", RDim, CDim, CDim, dcmplx(1.0d0, 0.0d0), QMat, LDQ, RMat, LDR, dcmplx(0.0d0, 0.0d0), &
            & CheckMt, RDim)
         XMaxm = 0.0d0; XMean = 0.0d0
			call MatrCmpr_R(RDim, CDim, CheckMt(1, 1), RDim, AMat(1, 1), LDA, XMaxm, XMean)
			write(*, "('QRdecmpsR_QR: QR Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
      if(allocated(work   )) deallocate(work   )
      if(allocated(tau    )) deallocate(tau    )
      if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(CheckMt)) deallocate(CheckMt)
				
   end subroutine QRdecmpsR_QR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine QRdecmpsR_QRPivot(RDim, CDim, AMat, LDA, QMat, LDQ, RMat, LDR, IfQRCheck)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  QRdecmpsR_QRPivot(RDim, CDim, AMat, LDA, QMat, LDQ, RMat, LDR, IfQRCheck)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create AMat(RDim, CDim) = QMat(RDim, CDim) * RMat(CDim, CDim) by 
!                 column-pivoted QR decomposition.
! KEYWORDS: Calculate the column-pivoted QR decomposition, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate the QR decomposition, complex version.
!     By the method of QR decomposition from Lapack with artificial column pivoting.
!
!     Only works for RDim >= CDim, with AMat = QMat * RMat;
!     If RDim < CDim, you can first calculate A^+ = QMat * RMat  --> A = (QMat * RMat)^+
!
!     Input: RDim      --> Dimension of AMat matrix;
!            CDim      --> Dimension of AMat matrix;
!            AMat      --> Input real RDim*CDim matrix (N>M);
!            LDA       --> Leading dimension of AMat matrix;
!            LDQ       --> Leading dimension of QMat matrix;
!            LDR       --> Leading dimension of RMat matrix;
!            IfQRCheck --> Whether to check the QR decomposition.
!
!     Outpt: QMat      --> The Q matrix;  --> By output, QMat is a real symmetric matrix;
!            RMat      --> The R matrix.  --> By output, RMat is a real matrix;
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
      logical IfQRCheck
		integer RDim, CDim                        ! Dimension of A matrix A(RDim, CDim) for QR decomp
      integer LDA, LDQ, LDR
		real(kind=kind(0.d0)) AMat(LDA, *)        ! Input A(RDim, CDim) matrix, Result: A=QR
		real(kind=kind(0.d0)) QMat(LDQ, *)        ! Result Output Q(RDim, CDim) matrix of QR
		real(kind=kind(0.d0)) RMat(LDR, *)        ! Result Output R(CDim, CDim) matrix of QR
!______________________________________________________________________________________________________________	  
!________________________________________ Temperory Quantities used ___________________________________________
!______________________________________________________________________________________________________________
      integer info, lwork, IMaxm
      integer I1, I2
      integer IVPT(CDim), IVPTM1(CDim)
      real(kind=kind(0.d0)) LogdDetR
      real(kind=kind(0.d0)) XMaxm, XMean
      real(kind=kind(0.d0)) Rtp0, Rtp1
      real(kind=kind(0.d0)) ColNorm(CDim), VecHelp(CDim)
      complex(kind=kind(0.d0)) LogzDet
      real(kind=kind(0.d0)), allocatable :: work(:)
      real(kind=kind(0.d0)), allocatable ::  tau(:)
      real(kind=kind(0.d0)), allocatable :: AMatTmp(:, :)
      real(kind=kind(0.d0)), allocatable :: RMatTmp(:, :)
      real(kind=kind(0.d0)), allocatable :: CheckMt(:, :)
!______________________________________________________________________________________________________________	  
!________________________________________ Main calculations ___________________________________________________
!______________________________________________________________________________________________________________
!_______________________________________________________________________________________________	  
!____________________________ 0. Initialization and allocation _________________________________
!_______________________________________________________________________________________________
      allocate(AMatTmp(RDim, CDim))
      AMatTmp(1:RDim, 1:CDim) = AMat(1:RDim, 1:CDim)
      
      allocate(RMatTmp(CDim, CDim))
      RMatTmp = 0.0d0
      
      allocate(tau(CDim))
      tau = 0.0d0
!_______________________________________________________________________________________________	  
!____________________________ 1. Perform the column-pivoted QR decomps for AMatTmp _____________
!_______________________________ AMatTmp * P = QMat * RMat _____________________________________
!_______________________________________________________________________________________________ 
!____________________________________________________________________________________	  
!__________________________ (0) AMatTmp * P = QMat * RMat ___________________________
!____________________________________________________________________________________
!________________________________________________________________________  
!__________________ [0] Query optimal amount of memory __________________
!________________________________________________________________________
      allocate(work(1))
      info = 0; IVPT = 0
      call dgeqp3(RDim, CDim, AMatTmp, RDim, IVPT, tau, work, -1, info)
      lwork = int(dble(work(1)))
      deallocate(work)
!________________________________________________________________________  
!__________________ [1] Perform colum-pivoted QR for AMatTmp ____________
!________________________________________________________________________
      allocate(work(lwork))
      IVPT = 0
      call dgeqp3(RDim, CDim, AMatTmp, RDim, IVPT, tau, work, lwork, info)
      deallocate(work) 
!____________________________________________________________________________________	  
!__________________________ (1) Obtain the RMatTmp matrix ___________________________
!______________________________ and the log of its determinant ______________________
!____________________________________________________________________________________
      RMatTmp = dcmplx(0.0d0, 0.0d0)
      call dlacpy("U", CDim, CDim, AMatTmp, RDim, RMatTmp, CDim)
      
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, CDim, +1
         LogzDet = LogzDet + log( cmplx(RMatTmp(I1, I1), 0.0d0) )
      enddo
!____________________________________________________________________________________	  
!__________________________ (2) Calculate R * P^T ___________________________________
!____________________________________________________________________________________
      do I1 = 1, CDim, +1
         IVPTM1(IVPT(I1)) = I1
      enddo

      do I1 = 1, CDim, +1
         do I2 = 1, CDim, +1
            RMat(I1, I2) = RMatTmp(I1, IVPTM1(I2))
         enddo
      enddo
!____________________________________________________________________________________	  
!__________________________ (3) Obtain the QMat matrix ______________________________
!____________________________________________________________________________________ 
      allocate(work(CDim))
      call dorg2r(RDim, CDim, CDim, AMatTmp, RDim, tau, work, info)
      deallocate(work)
      deallocate(tau)
      
      QMat(1:RDim, 1:CDim) = dcmplx(0.0d0, 0.0d0)
      call dlacpy("A", RDim, CDim, AMatTmp, RDim, QMat, LDQ)
!_______________________________________________________________________________________________	  
!____________________________ 2. Obtain the final QMat, RMat matrices as _______________________
!_______________________________ and get the logdDet of RMat matrix ____________________________
!_______________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________________ (0) Accumulate LogzDet of P^T matrix ____________________
!____________________________________________________________________________________
      I2 = 0
      do I1 = 1, CDim, +1
         if(IVPTM1(I1) /= I1) then
            I2 = I2 + 1
         end if
      enddo
      LogzDet = LogzDet + log( dcmplx((-1.0d0)**(mod(I2, 2)), 0.d0) )
!____________________________________________________________________________________	  
!__________________________ (1) Choose positive dDetR and Rescale Q and R ___________
!____________________________________________________________________________________
      LogdDetR = dreal(LogzDet)      
      if( dcos(dimag(LogzDet)) < 0.0d0 ) then
         do I1 = 1, RDim, +1
            QMat(I1, 1) = - QMat(I1, 1)
         enddo
         do I2 = 1, CDim, +1   
            RMat(1, I2) = - RMat(1, I2)
         enddo
      end if
!_______________________________________________________________________________________________	  
!____________________________ 3. Check the result of AMat = QMat * RMat ________________________
!_______________________________________________________________________________________________
      if(IfQRCheck) then
         allocate(CheckMt(RDim, CDim))
         CheckMt = 0.0d0
         call DGEMM("N", "N", RDim, CDim, CDim, dcmplx(1.0d0, 0.0d0), QMat, LDQ, RMat, LDR, dcmplx(0.0d0, 0.0d0), &
            & CheckMt, RDim)
         XMaxm = 0.0d0; XMean = 0.0d0
			call MatrCmpr_R(RDim, CDim, CheckMt(1, 1), RDim, AMat(1, 1), LDA, XMaxm, XMean)
			write(*, "('QRdecmpsR_QRPivot: QR Accurancy --> XMaxm, XMean = ', 2es25.16)") XMaxm, XMean
			deallocate(CheckMt)
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
      if(allocated(work   )) deallocate(work   )
      if(allocated(tau    )) deallocate(tau    )
      if(allocated(AMatTmp)) deallocate(AMatTmp)
      if(allocated(RMatTmp)) deallocate(RMatTmp)
      if(allocated(CheckMt)) deallocate(CheckMt)
				
   end subroutine QRdecmpsR_QRPivot
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
