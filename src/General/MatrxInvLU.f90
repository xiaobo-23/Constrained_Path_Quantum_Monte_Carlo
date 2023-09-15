!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to calculate the matrix inverse for real and complex square matrices separately.
!              There are also alternative subroutines used to calculate both the inverse and determinant for matrices.
!          Need the LAPACK mathematical library.
!          Note: A real matrix must have real matrix determinant.
!        
!          These subroutines calculate matrix inverse by partial row-pivoted LU decomposition method.
!          A = P*L*U with L as lower-triangular matrix with U as upper-triangular matrix, and diagonal elements 
!                    of L are one.
!          det(A) = det(P*L*U) = det(L) * det(P*U) = det(P*U), here P is the partial pivoting matrix with 
!                    row interchanges
!  
!          (1) The subroutine names with "Lq" means Applying the method of solving linear equation set to 
!                 calculate inv(A) as:
!                    A = P*L*U   -->   inv(A) * P * L = inv(U)
!                    In subroutine ZGETRI: First solve inv(A)*P = BMat  -->  inv(A) = BMat * P^T
!          (2) The subroutine names with "Iv" means Applying the method of simple inverse to calculate inv(A) as:
!                    A = P*L*U   -->   inv(A) = inv(P*L*U) = inv(L*U) * inv(P) = inv(U) * inv(L) * P^T
!
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!    MatrxInvLULqZ_NoDet  --> Subroutine to calculate the inverse of complex square matrix without  determinant , by Linear Equation method;  
!    MatrxInvLULqZ_Det    --> Subroutine to calculate the inverse of complex square matrix with     determinant , by Linear Equation method;
!    MatrxInvLULqZ_LogDet --> Subroutine to calculate the inverse of complex square matrix with log(determinant), by Linear Equation method;
!    MatrxInvLUIvZ_NoDet  --> Subroutine to calculate the inverse of complex square matrix without  determinant , by Direct Inverse method; 
!    MatrxInvLUIvZ_Det    --> Subroutine to calculate the inverse of complex square matrix with     determinant , by Direct Inverse method; 
!    MatrxInvLUIvZ_LogDet --> Subroutine to calculate the inverse of complex square matrix with log(determinant), by Direct Inverse method;
!  
!    MatrxInvLULqR_NoDet  --> Subroutine to calculate the inverse of real square matrix without  determinant , by Linear Equation method;  
!    MatrxInvLULqR_Det    --> Subroutine to calculate the inverse of real square matrix with     determinant , by Linear Equation method;
!    MatrxInvLULqR_LogDet --> Subroutine to calculate the inverse of real square matrix with log(determinant), by Linear Equation method;
!    MatrxInvLUIvR_NoDet  --> Subroutine to calculate the inverse of real square matrix without  determinant , by Direct Inverse method; 
!    MatrxInvLUIvR_Det    --> Subroutine to calculate the inverse of real square matrix with     determinant , by Direct Inverse method; 
!    MatrxInvLUIvR_LogDet --> Subroutine to calculate the inverse of real square matrix with log(determinant), by Direct Inverse method;
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   

!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
!################################################# For Complex Version ##################################################
!################################################# For Complex Version ##################################################
!################################################# For Complex Version ##################################################
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLULqZ_NoDet(NDim, zMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLULqZ_NoDet(NDim, zMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for complex square matrix, and the input A 
!                  is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, complex version. 
!
!     Input: NDim --> Dimension of input A matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: zMat --> Result output complex inv matrix of A matrix.
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
		integer NDim                     ! Dimension of A square matrix
		complex(kind=kind(0.d0)) zMat(NDim, NDim)     ! The input A matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: work(:)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_NoDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for zMat matrix ________________________
!________________________________________________________________________________________________
		call ZGETRF(NDim, NDim, zMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_NoDet: error in lapack subroutine ZGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Use inv(A) * P * L = inv(U) to slove inv(A) ____________________
!________________________________________________________________________________________________      
      call ZGETRI(NDim, zMat, NDim, ipiv, work, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_NoDet: error in lapack subroutine ZGETRI! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
		
   end subroutine MatrxInvLULqZ_NoDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLULqZ_Det(NDim, zMat, zDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLULqZ_Det(NDim, zMat, zDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for complex square matrix and its determinant, and 
!                  the input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, complex version. 
!
!     Input: NDim --> Dimension of input zMat matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: zMat --> Result output complex inv matrix of A matrix.
!            zDet --> Complex determinant.     
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
		integer NDim                     ! Dimension of zMat square matrix
      complex(kind=kind(0.d0)) zDet                 ! Determinant of zMat matrix
		complex(kind=kind(0.d0)) zMat(NDim, NDim)     ! The input zMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1
      
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: work(:)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_Det: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for zMat matrix ________________________
!________________________________________________________________________________________________
		call ZGETRF(NDim, NDim, zMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_Det: error in lapack subroutine ZGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Get the complex determinant of zMat matrix _____________________
!________________________________________________________________________________________________          
      zDet = dcmplx(1.0d0, 0.0d0)
      do I1 = 1, ndim
         if( ipiv(I1) == I1 ) then
            zDet = zDet * ( +zMat(I1, I1) )
         else
            zDet = zDet * ( -zMat(I1, I1) )
         end if  
     enddo
!________________________________________________________________________________________________	  
!____________________________ 2. Use inv(A) * P * L = inv(U) to slove inv(A) ____________________
!________________________________________________________________________________________________         
      call ZGETRI(NDim, zMat, NDim, ipiv, work, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_Det: error in lapack subroutine ZGETRI! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
		
   end subroutine MatrxInvLULqZ_Det
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLULqZ_LogDet(NDim, zMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLULqZ_LogDet(NDim, zMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for complex square matrix and its determinant, and 
!                  the input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, complex version. 
!
!     Input: NDim --> Dimension of input zMat matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: zMat --> Result output complex inv matrix of A matrix.
!            LogzDet --> Complex log(determinant).     
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
		integer NDim                     ! Dimension of zMat square matrix
      complex(kind=kind(0.d0)) LogzDet              ! Determinant of zMat matrix
		complex(kind=kind(0.d0)) zMat(NDim, NDim)     ! The input zMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1
      integer Itp0
      
      real(kind=kind(0.d0)) rp_pi
      real(kind=kind(0.d0)) ReLogzDet
      real(kind=kind(0.d0)) ImLogzDet
      
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: work(:)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_LogDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for zMat matrix ________________________
!________________________________________________________________________________________________
		call ZGETRF(NDim, NDim, zMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_LogDet: error in lapack subroutine ZGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Get the complex determinant of zMat matrix _____________________
!________________________________________________________________________________________________          
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, ndim
         if( ipiv(I1) == I1 ) then
            LogzDet = LogzDet + log(+zMat(I1, I1))
         else
            LogzDet = LogzDet + log(-zMat(I1, I1))
         end if  
      enddo
      
      rp_pi = dacos( -1.0d0 )
      ReLogzDet = dreal(LogzDet)
      ImLogzDet = dimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = dcmplx(ReLogzDet, ImLogzDet)
!________________________________________________________________________________________________	  
!____________________________ 2. Use inv(A) * P * L = inv(U) to slove inv(A) ____________________
!________________________________________________________________________________________________           
      call ZGETRI(NDim, zMat, NDim, ipiv, work, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqZ_LogDet: error in lapack subroutine ZGETRI! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
		
   end subroutine MatrxInvLULqZ_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLUIvZ_NoDet(NDim, zMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLUIvZ_NoDet(NDim, zMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for complex square matrix, and the 
!                  input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, complex version. 
!
!     Input: NDim --> Dimension of input A matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: zMat --> Result output complex inv matrix of A matrix.
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
		integer NDim                     ! Dimension of A square matrix
		complex(kind=kind(0.d0)) zMat(NDim, NDim)     ! The input A matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1, I2
      integer pv
      
      complex(kind=kind(0.d0)) Ztp1, Ztp2
      
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: work(:)
      
      complex(kind=kind(0.d0)), allocatable :: LMat(:, :)
      complex(kind=kind(0.d0)), allocatable :: UMat(:, :)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
      allocate(LMat(NDim, NDim), stat=ierror)
      allocate(UMat(NDim, NDim), stat=ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_NoDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
      
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_NoDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for zMat matrix ________________________
!________________________________________________________________________________________________
		call ZGETRF(NDim, NDim, zMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_NoDet: error in lapack subroutine ZGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Calculate inv(U) matrix ________________________________________
!________________________________________________________________________________________________      
      call ZTRTRI("U", "N", NDim, zMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_NoDet: error in lapack subroutine ZTRTRI for inv(U)! ierror = ', I4)") ierror
         stop
      end if
      
      UMat = dcmplx(0.0d0, 0.0d0)
      do I2 = 1, NDim
         do I1 = 1, I2
            UMat(I1, I2) = zMat(I1, I2)
         enddo
      enddo
!________________________________________________________________________________________________	  
!____________________________ 3. Calculate inv(L) matrix ________________________________________
!________________________________________________________________________________________________
      call ZTRTRI("L", "U", NDim, zMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_NoDet: error in lapack subroutine ZTRTRI for inv(L)! ierror = ', I4)") ierror
         stop
      end if
      
      LMat = dcmplx(0.0d0, 0.0d0)
      do I2 = 1, NDim-1
         LMat(I2, I2) = dcmplx(1.0d0, 0.0d0)
         do I1 = I2+1, NDim
            LMat(I1, I2) = zMat(I1, I2)
         enddo
      enddo
      LMat(NDim, NDim) = dcmplx(1.0d0, 0.0d0)
!________________________________________________________________________________________________	  
!____________________________ 4. Calculate inv(A) = inv(U) * inv(L) * P^T _______________________
!________________________________________________________________________________________________      
!_____________________________________________________________________________________  
!__________________________ (0) Calculate inv(L) * P^T _______________________________
!_____________________________________________________________________________________
      do I2 = NDim-1, 1, -1
         pv = ipiv(I2)
         if( pv .ne. I2 ) then
            call zswap(NDim, LMat(1, I2), 1, LMat(1, pv), 1)
         end if
      enddo
!_____________________________________________________________________________________  
!__________________________ (1) Calculate inv(U) * inv(L)*P^T ________________________
!_____________________________________________________________________________________
      zMat = dcmplx(0.0d0, 0.0d0)
      
      Ztp1 = dcmplx(1.0d0, 0.0d0)
      Ztp2 = dcmplx(0.0d0, 0.0d0)
      call zgemm("N", "N", NDim, NDim, NDim, Ztp1, UMat, NDim, LMat, NDim, Ztp2, zMat, NDim)
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
      if(allocated(LMat)) deallocate(LMat)
      if(allocated(UMat)) deallocate(UMat)
		
   end subroutine MatrxInvLUIvZ_NoDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLUIvZ_Det(NDim, zMat, zDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLUIvZ_Det(NDim, zMat, zDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for complex square matrix and its determinant, and 
!                  the input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, complex version. 
!
!     Input: NDim --> Dimension of input zMat matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: zMat --> Result output complex inv matrix of A matrix.
!            zDet --> Complex determinant.     
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
		integer NDim                     ! Dimension of zMat square matrix
      complex(kind=kind(0.d0)) zDet                 ! Determinant of zMat matrix
		complex(kind=kind(0.d0)) zMat(NDim, NDim)     ! The input zMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1, I2
      integer pv
      
      complex(kind=kind(0.d0)) Ztp1, Ztp2
      
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: work(:)
      
      complex(kind=kind(0.d0)), allocatable :: LMat(:, :)
      complex(kind=kind(0.d0)), allocatable :: UMat(:, :)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
      allocate(LMat(NDim, NDim), stat=ierror)
      allocate(UMat(NDim, NDim), stat=ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_Det: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
      
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_Det: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for zMat matrix ________________________
!________________________________________________________________________________________________
		call ZGETRF(NDim, NDim, zMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_Det: error in lapack subroutine ZGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Get the complex determinant of zMat matrix _____________________
!________________________________________________________________________________________________
      zDet = dcmplx(1.0d0, 0.0d0)
      do I1 = 1, ndim
         if( ipiv(I1) == I1 ) then
            zDet = zDet * ( +zMat(I1, I1) )
         else
            zDet = zDet * ( -zMat(I1, I1) )
         end if  
     enddo
!________________________________________________________________________________________________	  
!____________________________ 3. Calculate inv(U) matrix ________________________________________
!________________________________________________________________________________________________      
      call ZTRTRI("U", "N", NDim, zMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_Det: error in lapack subroutine ZTRTRI for inv(U)! ierror = ', I4)") ierror
         stop
      end if
      
      UMat = dcmplx(0.0d0, 0.0d0)
      do I2 = 1, NDim
         do I1 = 1, I2
            UMat(I1, I2) = zMat(I1, I2)
         enddo
      enddo
!________________________________________________________________________________________________	  
!____________________________ 4. Calculate inv(L) matrix ________________________________________
!________________________________________________________________________________________________
      call ZTRTRI("L", "U", NDim, zMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_Det: error in lapack subroutine ZTRTRI for inv(L)! ierror = ', I4)") ierror
         stop
      end if
      
      LMat = dcmplx(0.0d0, 0.0d0)
      do I2 = 1, NDim-1
         LMat(I2, I2) = dcmplx(1.0d0, 0.0d0)
         do I1 = I2+1, NDim
            LMat(I1, I2) = zMat(I1, I2)
         enddo
      enddo
      LMat(NDim, NDim) = dcmplx(1.0d0, 0.0d0)
!________________________________________________________________________________________________	  
!____________________________ 5. Calculate inv(A) = inv(U) * inv(L) * P^T _______________________
!________________________________________________________________________________________________      
!_____________________________________________________________________________________  
!__________________________ (0) Calculate inv(L) * P^T _______________________________
!_____________________________________________________________________________________
      do I2 = NDim-1, 1, -1
         pv = ipiv(I2)
         if( pv .ne. I2 ) then
            call zswap(NDim, LMat(1, I2), 1, LMat(1, pv), 1)
         end if
      enddo
!_____________________________________________________________________________________  
!__________________________ (1) Calculate inv(U) * inv(L)*P^T ________________________
!_____________________________________________________________________________________
      zMat = dcmplx(0.0d0, 0.0d0)
      
      Ztp1 = dcmplx(1.0d0, 0.0d0)
      Ztp2 = dcmplx(0.0d0, 0.0d0)
      call zgemm("N", "N", NDim, NDim, NDim, Ztp1, UMat, NDim, LMat, NDim, Ztp2, zMat, NDim)
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
      if(allocated(LMat)) deallocate(LMat)
      if(allocated(UMat)) deallocate(UMat)
		
   end subroutine MatrxInvLUIvZ_Det
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLUIvZ_LogDet(NDim, zMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLUIvZ_LogDet(NDim, zMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for complex square matrix and its determinant, and 
!                  the input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, complex version. 
!
!     Input: NDim --> Dimension of input zMat matrix;
!            zMat --> Input complex square matrix;
!
!     Outpt: zMat --> Result output complex inv matrix of A matrix.
!            LogzDet --> Complex log(determinant).    
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
		integer NDim                     ! Dimension of zMat square matrix
      complex(kind=kind(0.d0)) LogzDet              ! Determinant of zMat matrix
		complex(kind=kind(0.d0)) zMat(NDim, NDim)     ! The input zMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1, I2
      integer pv
      integer Itp0
      
      real(kind=kind(0.d0)) rp_pi
      real(kind=kind(0.d0)) ReLogzDet
      real(kind=kind(0.d0)) ImLogzDet
      
      complex(kind=kind(0.d0)) Ztp1, Ztp2
      
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: work(:)
      
      complex(kind=kind(0.d0)), allocatable :: LMat(:, :)
      complex(kind=kind(0.d0)), allocatable :: UMat(:, :)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
      allocate(LMat(NDim, NDim), stat=ierror)
      allocate(UMat(NDim, NDim), stat=ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_LogDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
      
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_LogDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for zMat matrix ________________________
!________________________________________________________________________________________________
		call ZGETRF(NDim, NDim, zMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_LogDet: error in lapack subroutine ZGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Get the complex determinant of zMat matrix _____________________
!________________________________________________________________________________________________
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, NDim
         if( ipiv(I1) == I1 ) then
            LogzDet = LogzDet + log(+zMat(I1, I1))
         else
            LogzDet = LogzDet + log(-zMat(I1, I1))
         end if  
      enddo
      
      rp_pi = dacos( -1.0d0 )
      ReLogzDet = dreal(LogzDet)
      ImLogzDet = dimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = dcmplx(ReLogzDet, ImLogzDet)
!________________________________________________________________________________________________	  
!____________________________ 3. Calculate inv(U) matrix ________________________________________
!________________________________________________________________________________________________      
      call ZTRTRI("U", "N", NDim, zMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_LogDet: error in lapack subroutine ZTRTRI for inv(U)! ierror = ', I4)") ierror
         stop
      end if
      
      UMat = dcmplx(0.0d0, 0.0d0)
      do I2 = 1, NDim
         do I1 = 1, I2
            UMat(I1, I2) = zMat(I1, I2)
         enddo
      enddo
!________________________________________________________________________________________________	  
!____________________________ 4. Calculate inv(L) matrix ________________________________________
!________________________________________________________________________________________________
      call ZTRTRI("L", "U", NDim, zMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvZ_LogDet: error in lapack subroutine ZTRTRI for inv(L)! ierror = ', I4)") ierror
         stop
      end if
      
      LMat = dcmplx(0.0d0, 0.0d0)
      do I2 = 1, NDim-1
         LMat(I2, I2) = dcmplx(1.0d0, 0.0d0)
         do I1 = I2+1, NDim
            LMat(I1, I2) = zMat(I1, I2)
         enddo
      enddo
      LMat(NDim, NDim) = dcmplx(1.0d0, 0.0d0)
!________________________________________________________________________________________________	  
!____________________________ 5. Calculate inv(A) = inv(U) * inv(L) * P^T _______________________
!________________________________________________________________________________________________      
!_____________________________________________________________________________________  
!__________________________ (0) Calculate inv(L) * P^T _______________________________
!_____________________________________________________________________________________
      do I2 = NDim-1, 1, -1
         pv = ipiv(I2)
         if( pv .ne. I2 ) then
            call zswap(NDim, LMat(1, I2), 1, LMat(1, pv), 1)
         end if
      enddo
!_____________________________________________________________________________________  
!__________________________ (1) Calculate inv(U) * inv(L)*P^T ________________________
!_____________________________________________________________________________________
      zMat = dcmplx(0.0d0, 0.0d0)
      
      Ztp1 = dcmplx(1.0d0, 0.0d0)
      Ztp2 = dcmplx(0.0d0, 0.0d0)
      call zgemm("N", "N", NDim, NDim, NDim, Ztp1, UMat, NDim, LMat, NDim, Ztp2, zMat, NDim)
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
      if(allocated(LMat)) deallocate(LMat)
      if(allocated(UMat)) deallocate(UMat)
		
   end subroutine MatrxInvLUIvZ_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
!################################################# For Real Version #####################################################
!################################################# For Real Version #####################################################
!################################################# For Real Version #####################################################
!########################################################################################################################
!########################################################################################################################
!######################################################################################################################## 
   
   
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLULqR_NoDet(NDim, dMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLULqR_NoDet(NDim, dMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for real square matrix, and the input A is 
!                    overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, real version. 
!
!     Input: NDim --> Dimension of input A matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: dMat --> Result output real inv matrix of A matrix.
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
		integer NDim                     ! Dimension of A square matrix
		real(kind=kind(0.d0)) dMat(NDim, NDim)        ! The input A matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: work(:)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_NoDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for dMat matrix ________________________
!________________________________________________________________________________________________
		call DGETRF(NDim, NDim, dMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_NoDet: error in lapack subroutine DGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Use inv(A) * P * L = inv(U) to slove inv(A) ____________________
!________________________________________________________________________________________________      
      call DGETRI(NDim, dMat, NDim, ipiv, work, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_NoDet: error in lapack subroutine DGETRI! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
		
   end subroutine MatrxInvLULqR_NoDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLULqR_Det(NDim, dMat, dDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLULqR_Det(NDim, dMat, dDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for real square matrix and its determinant, and 
!                  the input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, real version. 
!
!     Input: NDim --> Dimension of input dMat matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: dMat --> Result output real inv matrix of A matrix.
!            dDet --> real determinant.     
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
		integer NDim                     ! Dimension of dMat square matrix
      real(kind=kind(0.d0)) dDet                    ! Determinant of dMat matrix
		real(kind=kind(0.d0)) dMat(NDim, NDim)        ! The input dMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1
      
      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: work(:)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_Det: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for dMat matrix ________________________
!________________________________________________________________________________________________
		call DGETRF(NDim, NDim, dMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_Det: error in lapack subroutine DGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Get the real determinant of dMat matrix ________________________
!________________________________________________________________________________________________          
      dDet = dcmplx(1.0d0, 0.0d0)
      do I1 = 1, ndim
         if( ipiv(I1) == I1 ) then
            dDet = dDet * ( +dMat(I1, I1) )
         else
            dDet = dDet * ( -dMat(I1, I1) )
         end if  
     enddo
!________________________________________________________________________________________________	  
!____________________________ 2. Use inv(A) * P * L = inv(U) to slove inv(A) ____________________
!________________________________________________________________________________________________         
      call DGETRI(NDim, dMat, NDim, ipiv, work, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_Det: error in lapack subroutine DGETRI! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
		
   end subroutine MatrxInvLULqR_Det
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLULqR_LogDet(NDim, dMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLULqR_LogDet(NDim, dMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for real square matrix and its determinant, and 
!                  the input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, real version. 
!     The matrix must have real determinant, but we want to incorporate the -1 sign, so we use complex LogzDet here.
!
!     Input: NDim --> Dimension of input dMat matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: dMat --> Result output real inv matrix of A matrix.
!            LogzDet --> Complex log(determinant).     
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
		integer NDim                     ! Dimension of dMat square matrix
      complex(kind=kind(0.d0)) LogzDet              ! Determinant of dMat matrix
		real(kind=kind(0.d0)) dMat(NDim, NDim)        ! The input dMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1
      integer Itp0
      
      real(kind=kind(0.d0)) rp_pi
      real(kind=kind(0.d0)) ReLogzDet
      real(kind=kind(0.d0)) ImLogzDet
      
      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: work(:)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_LogDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for dMat matrix ________________________
!________________________________________________________________________________________________
		call DGETRF(NDim, NDim, dMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_LogDet: error in lapack subroutine DGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Get the complex determinant of dMat matrix _____________________
!________________________________________________________________________________________________          
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, ndim
         if( ipiv(I1) == I1 ) then
            LogzDet = LogzDet + log(dcmplx(+dMat(I1, I1), 0.d0))
         else
            LogzDet = LogzDet + log(dcmplx(-dMat(I1, I1), 0.d0))
         end if  
      enddo
      
      rp_pi = dacos( -1.0d0 )
      ReLogzDet = dreal(LogzDet)
      ImLogzDet = dimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = dcmplx(ReLogzDet, ImLogzDet)
!________________________________________________________________________________________________	  
!____________________________ 2. Use inv(A) * P * L = inv(U) to slove inv(A) ____________________
!________________________________________________________________________________________________           
      call DGETRI(NDim, dMat, NDim, ipiv, work, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLULqR_LogDet: error in lapack subroutine DGETRI! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
		
   end subroutine MatrxInvLULqR_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLUIvR_NoDet(NDim, dMat)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLUIvR_NoDet(NDim, dMat)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for real square matrix, and the input 
!                   A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, real version. 
!
!     Input: NDim --> Dimension of input A matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: dMat --> Result output real inv matrix of A matrix.
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
		integer NDim                     ! Dimension of A square matrix
		real(kind=kind(0.d0)) dMat(NDim, NDim)         ! The input A matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1, I2
      integer pv
      
      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: work(:)
      
      real(kind=kind(0.d0)), allocatable :: LMat(:, :)
      real(kind=kind(0.d0)), allocatable :: UMat(:, :)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
      allocate(LMat(NDim, NDim), stat=ierror)
      allocate(UMat(NDim, NDim), stat=ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_NoDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
      
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_NoDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for dMat matrix ________________________
!________________________________________________________________________________________________
		call DGETRF(NDim, NDim, dMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_NoDet: error in lapack subroutine DGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Calculate inv(U) matrix ________________________________________
!________________________________________________________________________________________________      
      call DTRTRI("U", "N", NDim, dMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_NoDet: error in lapack subroutine DTRTRI for inv(U)! ierror = ', I4)") ierror
         stop
      end if
      
      UMat = 0.0d0
      do I2 = 1, NDim
         do I1 = 1, I2
            UMat(I1, I2) = dMat(I1, I2)
         enddo
      enddo
!________________________________________________________________________________________________	  
!____________________________ 3. Calculate inv(L) matrix ________________________________________
!________________________________________________________________________________________________
      call DTRTRI("L", "U", NDim, dMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_NoDet: error in lapack subroutine DTRTRI for inv(L)! ierror = ', I4)") ierror
         stop
      end if
      
      LMat = 0.0d0
      do I2 = 1, NDim-1
         LMat(I2, I2) = 1.0d0
         do I1 = I2+1, NDim
            LMat(I1, I2) = dMat(I1, I2)
         enddo
      enddo
      LMat(NDim, NDim) = 1.0d0
!________________________________________________________________________________________________	  
!____________________________ 4. Calculate inv(A) = inv(U) * inv(L) * P^T _______________________
!________________________________________________________________________________________________      
!_____________________________________________________________________________________  
!__________________________ (0) Calculate inv(L) * P^T _______________________________
!_____________________________________________________________________________________
      do I2 = NDim-1, 1, -1
         pv = ipiv(I2)
         if( pv .ne. I2 ) then
            call dswap(NDim, LMat(1, I2), 1, LMat(1, pv), 1)
         end if
      enddo
!_____________________________________________________________________________________  
!__________________________ (1) Calculate inv(U) * inv(L)*P^T ________________________
!_____________________________________________________________________________________
      dMat = 0.0d0
      call dgemm("N", "N", NDim, NDim, NDim, 1.0d0, UMat, NDim, LMat, NDim, 0.0d0, dMat, NDim)
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
      if(allocated(LMat)) deallocate(LMat)
      if(allocated(UMat)) deallocate(UMat)
		
   end subroutine MatrxInvLUIvR_NoDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLUIvR_Det(NDim, dMat, dDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLUIvR_Det(NDim, dMat, dDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for real square matrix and its determinant, and 
!                  the input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, real version. 
!
!     Input: NDim --> Dimension of input dMat matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: dMat --> Result output real inv matrix of A matrix.
!            dDet --> real determinant.     
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
		integer NDim                     ! Dimension of dMat square matrix
      real(kind=kind(0.d0)) dDet                    ! Determinant of dMat matrix
		real(kind=kind(0.d0)) dMat(NDim, NDim)        ! The input dMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1, I2
      integer pv
      
      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: work(:)
      
      real(kind=kind(0.d0)), allocatable :: LMat(:, :)
      real(kind=kind(0.d0)), allocatable :: UMat(:, :)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
      allocate(LMat(NDim, NDim), stat=ierror)
      allocate(UMat(NDim, NDim), stat=ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_Det: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
      
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_Det: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for dMat matrix ________________________
!________________________________________________________________________________________________
		call DGETRF(NDim, NDim, dMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_Det: error in lapack subroutine DGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Get the real determinant of dMat matrix _________________________
!________________________________________________________________________________________________
      dDet = 1.0d0
      do I1 = 1, ndim
         if( ipiv(I1) == I1 ) then
            dDet = dDet * ( +dMat(I1, I1) )
         else
            dDet = dDet * ( -dMat(I1, I1) )
         end if  
     enddo
!________________________________________________________________________________________________	  
!____________________________ 3. Calculate inv(U) matrix ________________________________________
!________________________________________________________________________________________________      
      call DTRTRI("U", "N", NDim, dMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_Det: error in lapack subroutine DTRTRI for inv(U)! ierror = ', I4)") ierror
         stop
      end if
      
      UMat = 0.0d0
      do I2 = 1, NDim
         do I1 = 1, I2
            UMat(I1, I2) = dMat(I1, I2)
         enddo
      enddo
!________________________________________________________________________________________________	  
!____________________________ 4. Calculate inv(L) matrix ________________________________________
!________________________________________________________________________________________________
      call DTRTRI("L", "U", NDim, dMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_Det: error in lapack subroutine DTRTRI for inv(L)! ierror = ', I4)") ierror
         stop
      end if
      
      LMat = 0.0d0
      do I2 = 1, NDim-1
         LMat(I2, I2) = 1.0d0
         do I1 = I2+1, NDim
            LMat(I1, I2) = dMat(I1, I2)
         enddo
      enddo
      LMat(NDim, NDim) = 1.0d0
!________________________________________________________________________________________________	  
!____________________________ 5. Calculate inv(A) = inv(U) * inv(L) * P^T _______________________
!________________________________________________________________________________________________      
!_____________________________________________________________________________________  
!__________________________ (0) Calculate inv(L) * P^T _______________________________
!_____________________________________________________________________________________
      do I2 = NDim-1, 1, -1
         pv = ipiv(I2)
         if( pv .ne. I2 ) then
            call dswap(NDim, LMat(1, I2), 1, LMat(1, pv), 1)
         end if
      enddo
!_____________________________________________________________________________________  
!__________________________ (1) Calculate inv(U) * inv(L)*P^T ________________________
!_____________________________________________________________________________________
      dMat = 0.0d0
      call dgemm("N", "N", NDim, NDim, NDim, 1.0d0, UMat, NDim, LMat, NDim, 0.0d0, dMat, NDim)
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
      if(allocated(LMat)) deallocate(LMat)
      if(allocated(UMat)) deallocate(UMat)
		
   end subroutine MatrxInvLUIvR_Det
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine MatrxInvLUIvR_LogDet(NDim, dMat, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  MatrxInvLUIvR_LogDet(NDim, dMat, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate the matrix inverse for real square matrix and its determinant, and 
!                  the input A is overwrite by Inv(A).
! KEYWORDS: Calculate Matrix inverse, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Calculate Matrix inverse, real version. 
!     The matrix must have real determinant, but we want to incorporate the -1 sign, so we use complex LogzDet here.
!
!     Input: NDim --> Dimension of input dMat matrix;
!            dMat --> Input real square matrix;
!
!     Outpt: dMat --> Result output real inv matrix of A matrix.
!            LogzDet --> complex log(determinant).    
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
		integer NDim                     ! Dimension of dMat square matrix
      complex(kind=kind(0.d0)) LogzDet              ! Determinant of dMat matrix
		real(kind=kind(0.d0)) dMat(NDim, NDim)        ! The input dMat matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer IError
      
      integer I1, I2
      integer pv
      integer Itp0
      
      real(kind=kind(0.d0)) rp_pi
      real(kind=kind(0.d0)) ReLogzDet
      real(kind=kind(0.d0)) ImLogzDet

      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: work(:)
      
      real(kind=kind(0.d0)), allocatable :: LMat(:, :)
      real(kind=kind(0.d0)), allocatable :: UMat(:, :)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
      allocate(LMat(NDim, NDim), stat=ierror)
      allocate(UMat(NDim, NDim), stat=ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_LogDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
      
		allocate(ipiv(NDim), stat=ierror)
		allocate(work(NDim), stat=ierror)
		if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_LogDet: can not allocate enough memory! ierror = ', I4)") ierror
         stop
      end if
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Matrix Inverse __________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!____________________________ 1. A = P*L*U decomposition for dMat matrix ________________________
!________________________________________________________________________________________________
		call DGETRF(NDim, NDim, dMat, NDim, ipiv, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_LogDet: error in lapack subroutine DGETRF! ierror = ', I4)") ierror
         stop
      end if
!________________________________________________________________________________________________	  
!____________________________ 2. Get the complex determinant of dMat matrix _____________________
!________________________________________________________________________________________________
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, ndim
         if( ipiv(I1) == I1 ) then
            LogzDet = LogzDet + log(dcmplx(+dMat(I1, I1), 0.d0))
         else
            LogzDet = LogzDet + log(dcmplx(-dMat(I1, I1), 0.d0))
         end if 
      enddo
      
      rp_pi = dacos( -1.0d0 )
      ReLogzDet = dreal(LogzDet)
      ImLogzDet = dimag(LogzDet)
      Itp0 = nint(ImLogzDet/2.0d0/rp_pi)
      ImLogzDet = ImLogzDet - 2.0d0*rp_pi*dble(Itp0)
      LogzDet = dcmplx(ReLogzDet, ImLogzDet)
!________________________________________________________________________________________________	  
!____________________________ 3. Calculate inv(U) matrix ________________________________________
!________________________________________________________________________________________________      
      call DTRTRI("U", "N", NDim, dMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_LogDet: error in lapack subroutine DTRTRI for inv(U)! ierror = ', I4)") ierror
         stop
      end if
      
      UMat = 0.0d0
      do I2 = 1, NDim
         do I1 = 1, I2
            UMat(I1, I2) = dMat(I1, I2)
         enddo
      enddo
!________________________________________________________________________________________________	  
!____________________________ 4. Calculate inv(L) matrix ________________________________________
!________________________________________________________________________________________________
      call DTRTRI("L", "U", NDim, dMat, NDim, Ierror)
      if ( ierror /= 0 ) then
         write(*, "('MatrxInvLUIvR_LogDet: error in lapack subroutine DTRTRI for inv(L)! ierror = ', I4)") ierror
         stop
      end if
      
      LMat = 0.0d0
      do I2 = 1, NDim-1
         LMat(I2, I2) = 1.0d0
         do I1 = I2+1, NDim
            LMat(I1, I2) = dMat(I1, I2)
         enddo
      enddo
      LMat(NDim, NDim) = 1.0d0
!________________________________________________________________________________________________	  
!____________________________ 5. Calculate inv(A) = inv(U) * inv(L) * P^T _______________________
!________________________________________________________________________________________________      
!_____________________________________________________________________________________  
!__________________________ (0) Calculate inv(L) * P^T _______________________________
!_____________________________________________________________________________________
      do I2 = NDim-1, 1, -1
         pv = ipiv(I2)
         if( pv .ne. I2 ) then
            call dswap(NDim, LMat(1, I2), 1, LMat(1, pv), 1)
         end if
      enddo
!_____________________________________________________________________________________  
!__________________________ (1) Calculate inv(U) * inv(L)*P^T ________________________
!_____________________________________________________________________________________
      dMat = 0.0d0
      call dgemm("N", "N", NDim, NDim, NDim, 1.0d0, UMat, NDim, LMat, NDim, 0.0d0, dMat, NDim)
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(work)) deallocate(work)
      if(allocated(LMat)) deallocate(LMat)
      if(allocated(UMat)) deallocate(UMat)
		
   end subroutine MatrxInvLUIvR_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   