!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Two subroutines used to construct compare two different matrices and calculate the difference of these two
!              matrices. Calculate the largest difference element and calculate the average difference for all elements.
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!    MatrCmpr_R --> Subroutine to compare and calculate difference of two real matrices;
!    MatrCmpr_C --> Subroutine to compare and calculate difference of two complex matrices;
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrCmpr_R(N, M, A, LDA, B, LDB, XMax, XMean)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrCmpr_R(N, M, A, LDA, B, LDB, XMax, XMean)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compare two real different matrices and calculate 
!                  their differences in elements.
! KEYWORDS: Compare two matrices, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Compare two matrices, real version. 
!
!     Input: N   --> Dimension of A, B matrix as A(N, M) and B(N, M);
!            M   --> Dimension of A, B matrix as A(N, M) and B(N, M);
!            A   --> Input matrix A(N, M);
!            LDA --> Leading dimension of A matrix;
!            B   --> Input matrix A(N, M);
!            LDB --> Leading dimension of B matrix;
!     
!     Outpt: XMax  --> The largest absolute difference;
!            XMean --> The average difference in all elements.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer N           ! Dimension of A, B matrices
		integer M           ! Dimension of A, B matrices
      integer LDA, LDB
		real(kind=kind(0.d0)) XMax       ! Maximum difference for matrix elements in A, B matrices
		real(kind=kind(0.d0)) XMean      ! Average difference for matrix elements in A, B matrices
		real(kind=kind(0.d0)) A(LDA, *)  ! Input A matrix
		real(kind=kind(0.d0)) B(LDB, *)  ! Input B matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1      ! Loop integer 
		integer I2      ! Loop integer
		real(kind=kind(0.d0)) Diff   ! Difference between the elements in A and B
      
      real(kind=kind(0.d0)), allocatable :: XMaxVec(:)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Creating matrix _________________________________________
!______________________________________________________________________________________________________________
		XMax  = 0.d0
		XMean = 0.d0
      
      allocate(XMaxVec(N))
      XMaxVec = 0.d0
      
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2, Diff)
   !$OMP DO REDUCTION(+ : XMean)
		do I1 = 1, N
			do I2 = 1, M
				Diff = dabs(B(I1, I2) - A(I1, I2))
				if(Diff .gt. XMaxVec(I1)) then
					XMaxVec(I1) = Diff
				end if
				XMean = XMean + Diff
			enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
				
		XMean = XMean / dble(M * N)
      XMax = maxval(XMaxVec)
      
      deallocate(XMaxVec)
		
	end subroutine MatrCmpr_R
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrCmpr_C(N, M, A, LDA, B, LDB, XMax, XMean)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrCmpr_C(N, M, A, LDA, B, LDB, XMax, XMean)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compare two different complex matrices and calculate 
!                 their differences in elements.
! KEYWORDS: Compare two matrices, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Compare two matrices, complex version. 
!
!     Input: N   --> Dimension of A, B matrix as A(N, M) and B(N, M);
!            M   --> Dimension of A, B matrix as A(N, M) and B(N, M);
!            A   --> Input matrix A(N, M);
!            LDA --> Leading dimension of A matrix;
!            B   --> Input matrix A(N, M);
!            LDB --> Leading dimension of B matrix;
!     
!     Outpt: XMax  --> The largest absolute difference;
!            XMean --> The average difference in all elements.
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer N              ! Dimension of A, B matrices
		integer M              ! Dimension of A, B matrices
      integer LDA, LDB
		real(kind=kind(0.d0)) XMax          ! Maximum difference for matrix elements in A, B matrices
		real(kind=kind(0.d0)) XMean         ! Average difference for matrix elements in A, B matrices
		complex(kind=kind(0.d0)) A(LDA, *)  ! Input A matrix
		complex(kind=kind(0.d0)) B(LDB, *)  ! Input B matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1      ! Loop integer 
		integer I2      ! Loop integer
		real(kind=kind(0.d0)) Diff   ! Difference between the elements in A and B
      
      real(kind=kind(0.d0)), allocatable :: XMaxVec(:)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Creating matrix _________________________________________
!______________________________________________________________________________________________________________			
		XMax  = 0.d0
		XMean = 0.d0
      
      allocate(XMaxVec(N))
      XMaxVec = 0.d0
      
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2, Diff)
   !$OMP DO REDUCTION(+ : XMean)
		do I1 = 1, N
			do I2 = 1, M
				Diff = dsqrt( dreal( (B(I1, I2) - A(I1, I2)) * dconjg(B(I1, I2) - A(I1, I2)) ) )
				if(Diff .gt. XMaxVec(I1)) then
					XMaxVec(I1) = Diff
				end if
				XMean = XMean + Diff
			enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL

		XMean = XMean / dble(M * N)
      XMax = maxval(XMaxVec)
      
      deallocate(XMaxVec)
		
   end subroutine MatrCmpr_C
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
