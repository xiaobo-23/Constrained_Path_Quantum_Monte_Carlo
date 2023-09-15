!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Two subroutines used to construct the interface for matrix multiplication, including the real matrices
!              multiplication and complex matrices multiplication. The kernel subroutine used to perform the 
!              matrix multiplication is subroutines from the LAPACK mathematical library. 
!
!          Need LAPACK mathematical Library.
!
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!    MatrxMultR --> Subroutine to perform the real matrices multiplication as C=A*B;
!    MatrxMultZ --> Subroutine to perform the complex matrices multiplication as C=A*B.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrxMultR(N, M, K, A, B, C)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrxMultR(N, M, K, A, B, C)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the real matrices multiplization  as: C = A * B
! KEYWORDS: Matrice multiplication, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrice multiplication, real version. 
!
!     Input: N --> Dimension of A matrix as A(N, M);
!            N --> Dimension of A matrix as A(N, M), and B matrix as B(M, K);
!            K --> Dimension of B matrix as B(M, K);
!            A --> Input matrix A(N, M);
!            B --> Input matrix B(M, K).
!             
!     Outpt: C --> Result output matrix as C(N, K).
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
		integer N            ! Dimension of A matrix A(N, M)
		integer M            ! Dimension of A matrix A(N, M)
		integer K            ! Dimension of B matrix B(M, K) 
		real(kind=kind(0.d0)) A(N, M)     ! Input A matrix
		real(kind=kind(0.d0)) B(M, K)     ! Input B matrix
		real(kind=kind(0.d0)) C(N, K)     ! Ouput C matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer LDA   ! leading dimension of A matrix
		integer LDB   ! leading dimension of B matrix
		integer LDC   ! leading dimension of C matrix
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LDA = N
		LDB = M
		LDC = N
				
		call DGEMM("N", "N", N, K, M, 1.0d0, A, LDA, B, LDB, 0.0d0, C, LDC) ! Call BLAS subroutine
		
	end subroutine MatrxMultR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrxMultZ(N, M, K, A, B, C)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrxMultZ(N, M, K, A, B, C)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the complex matrices multiplization  as: C = A * B
! KEYWORDS: Matrice multiplication, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrice multiplication, complex version. 
!
!     Input: N --> Dimension of A matrix as A(N, M);
!            N --> Dimension of A matrix as A(N, M), and B matrix as B(M, K);
!            K --> Dimension of B matrix as B(M, K);
!            A --> Input matrix A(N, M);
!            B --> Input matrix B(M, K).
!             
!     Outpt: C --> Result output matrix as C(N, K).
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
		integer N            ! Dimension of A matrix A(N, M)
		integer M            ! Dimension of A matrix A(N, M)
		integer K            ! Dimension of B matrix B(M, K) 
		complex(kind=kind(0.d0)) A(N, M)  ! Input A matrix
		complex(kind=kind(0.d0)) B(M, K)  ! Input B matrix
		complex(kind=kind(0.d0)) C(N, K)  ! Ouput C matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer LDA   ! leading dimension of A matrix
		integer LDB   ! leading dimension of B matrix
		integer LDC   ! leading dimension of C matrix
		
		complex(kind=kind(0.d0)) Ztp1  ! Complex temporary number used in calculations
		complex(kind=kind(0.d0)) Ztp2  ! Complex temporary number used in calculations
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LDA = N
		LDB = M
		LDC = N
		Ztp1 = dcmplx(1.0d0, 0.0d0)
		Ztp2 = dcmplx(0.0d0, 0.0d0)
				
		call ZGEMM("N", "N", N, K, M, Ztp1, A, LDA, B, LDB, Ztp2, C, LDC) ! Call BLAS subroutine
		
	end subroutine MatrxMultZ
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
