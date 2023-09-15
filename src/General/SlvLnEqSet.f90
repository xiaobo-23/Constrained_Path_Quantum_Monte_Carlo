!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Several subroutines used to solve the linear equation set XMat*AMat = BMat or AMat*XMat = BMat, where AMat 
!              is always a square matrix, and both XMat and BMat can be rectangular matrices. Apply the LU 
!              decomposition method.
!
!          First perform the partial-pivoted LU decomposition for AMat as AMat = P * L * U, then solve XMat matrix.
!
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!
!    In all the following subroutines, the output result XMat matrix is written into BMat matrix.
!
!    SlvLnEqSetZ_Left_NoDet  --> Subroutine to solve XMat*AMat = BMat, without     determinant  of AMat, complex version;  
!    SlvLnEqSetZ_Left_Det    --> Subroutine to solve XMat*AMat = BMat, with        determinant  of AMat, complex version;  
!    SlvLnEqSetZ_Left_LogDet --> Subroutine to solve XMat*AMat = BMat, with    log(determinant) of AMat, complex version;
!  
!    SlvLnEqSetZ_Rght_NoDet  --> Subroutine to solve AMat*XMat = BMat, without     determinant  of AMat, complex version;  
!    SlvLnEqSetZ_Rght_Det    --> Subroutine to solve AMat*XMat = BMat, with        determinant  of AMat, complex version;  
!    SlvLnEqSetZ_Rght_LogDet --> Subroutine to solve AMat*XMat = BMat, with    log(determinant) of AMat, complex version.
!  
!
!    SlvLnEqSetR_Left_NoDet  --> Subroutine to solve XMat*AMat = BMat, without     determinant  of AMat, real version;  
!    SlvLnEqSetR_Left_Det    --> Subroutine to solve XMat*AMat = BMat, with        determinant  of AMat, real version;  
!    SlvLnEqSetR_Left_LogDet --> Subroutine to solve XMat*AMat = BMat, with    log(determinant) of AMat, real version;
!  
!    SlvLnEqSetR_Rght_NoDet  --> Subroutine to solve AMat*XMat = BMat, without     determinant  of AMat, real version;  
!    SlvLnEqSetR_Rght_Det    --> Subroutine to solve AMat*XMat = BMat, with        determinant  of AMat, real version;  
!    SlvLnEqSetR_Rght_LogDet --> Subroutine to solve AMat*XMat = BMat, with    log(determinant) of AMat, real version.
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
	subroutine SlvLnEqSetZ_Left_NoDet(ND1, ND2, AMat, LDA, BMat, LDB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetZ_Left_NoDet(ND1, ND2, AMat, LDA, BMat, LDB)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set XMat * AMat = BMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2). 
!     First , take a transpose (or conjugatetranspose) as AMat^T * XMat^T = BMat^T (or AMat^+ * XMat^+ = BMat^+);
!     Second, partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Third , Solve AMat^T from AMat^T * XMat^T = BMat^T (or XMat^+ from AMat^+ * XMat^+ = BMat^+);
!     Forth , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
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
		integer ND1, ND2
      integer LDA, LDB
      complex(kind=kind(0.d0)) AMat(LDA, *)
      complex(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, I2
      integer Info
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: BMatT(:, :)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND2))
      ipiv = 0
      
		allocate(BMatT(ND2, ND1))
      BMatT = dcmplx(0.0d0, 0.0d0)
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
      call ZGETRF(ND2, ND2, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetZ_Left_NoDet: ZGETRF error with Info = ', I4)") Info
         write(*, "('                        The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                        Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetZ_Left_NoDet: ZGETRF error with Info = ', I4)") Info
         write(*, "('                        LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                        Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat^T * XMat^T = BMat^T equation set ________________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) Obtain the transpose of BMat matrix _________________________
!____________________________________________________________________________________     
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND2, +1
         do I2 = 1, ND1, +1
            BMatT(I1, I2) = BMat(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!____________________________________________________________________________________	  
!__________________ (1) Solve AMat^T * XMat^T = BMat^T ______________________________
!____________________________________________________________________________________ 
      call ZGETRS("T", ND2, ND1, AMat, LDA, ipiv, BMatT, ND2, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetZ_Left_NoDet: ZGETRS error with Info = ', I4)") Info
         write(*, "('                        Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Output the result XMat matrx to BMat matrix ________________________
!________________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND1, +1
         do I2 = 1, ND2, +1
            BMat(I1, I2) = BMatT(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL     
!________________________________________________________________________________________________	  
!________________________ 4. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(BMatT)) deallocate(BMatT)
		
   end subroutine SlvLnEqSetZ_Left_NoDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetZ_Left_Det(ND1, ND2, AMat, LDA, BMat, LDB, zDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetZ_Left_Det(ND1, ND2, AMat, LDA, BMat, LDB, zDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set XMat * AMat = BMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2). 
!     First , take a transpose (or conjugatetranspose) as AMat^T * XMat^T = BMat^T (or AMat^+ * XMat^+ = BMat^+);
!     Second, partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Third , Solve AMat^T from AMat^T * XMat^T = BMat^T (or XMat^+ from AMat^+ * XMat^+ = BMat^+);
!     Forth , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
!            zDet --> Determinant of AMat(ND2, ND2) matrix.
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
		integer ND1, ND2
      integer LDA, LDB
      complex(kind=kind(0.d0)) zDet
      complex(kind=kind(0.d0)) AMat(LDA, *)
      complex(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, I2
      integer Info
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: BMatT(:, :)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND2))
      ipiv = 0
      
		allocate(BMatT(ND2, ND1))
      BMatT = dcmplx(0.0d0, 0.0d0)
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) LU decomposition for AMat matrix ____________________________
!____________________________________________________________________________________ 
      call ZGETRF(ND2, ND2, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetZ_Left_Det: ZGETRF error with Info = ', I4)") Info
         write(*, "('                      The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                      Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetZ_Left_Det: ZGETRF error with Info = ', I4)") Info
         write(*, "('                      LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                      Will stop running here!!!')")
         stop
      end if
!____________________________________________________________________________________	  
!__________________ (1) Determinant of AMat(ND2, ND2) matrix ________________________
!____________________________________________________________________________________
      zDet = dcmplx(1.0d0, 0.0d0)
      do I1 = 1, ND2, +1
         if(ipiv(I1) == I1) then
            zDet = zDet * ( + AMat(I1, I1) )
         else
            zDet = zDet * ( - AMat(I1, I1) )
         end if
      enddo
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat^T * XMat^T = BMat^T equation set ________________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) Obtain the transpose of BMat matrix _________________________
!____________________________________________________________________________________     
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND2, +1
         do I2 = 1, ND1, +1
            BMatT(I1, I2) = BMat(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!____________________________________________________________________________________	  
!__________________ (1) Solve AMat^T * XMat^T = BMat^T ______________________________
!____________________________________________________________________________________ 
      call ZGETRS("T", ND2, ND1, AMat, LDA, ipiv, BMatT, ND2, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetZ_Left_Det: ZGETRS error with Info = ', I4)") Info
         write(*, "('                      Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Output the result XMat matrx to BMat matrix ________________________
!________________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND1, +1
         do I2 = 1, ND2, +1
            BMat(I1, I2) = BMatT(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL     
!________________________________________________________________________________________________	  
!________________________ 4. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(BMatT)) deallocate(BMatT)
		
   end subroutine SlvLnEqSetZ_Left_Det
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetZ_Left_LogDet(ND1, ND2, AMat, LDA, BMat, LDB, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetZ_Left_LogDet(ND1, ND2, AMat, LDA, BMat, LDB, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set XMat * AMat = BMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2). 
!     First , take a transpose (or conjugatetranspose) as AMat^T * XMat^T = BMat^T (or AMat^+ * XMat^+ = BMat^+);
!     Second, partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Third , Solve AMat^T from AMat^T * XMat^T = BMat^T (or XMat^+ from AMat^+ * XMat^+ = BMat^+);
!     Forth , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
!            LogzDet --> Log(Determinant) of AMat(ND2, ND2) matrix.
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
		integer ND1, ND2
      integer LDA, LDB
      complex(kind=kind(0.d0)) LogzDet
      complex(kind=kind(0.d0)) AMat(LDA, *)
      complex(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, I2
      integer Info
      integer, allocatable :: ipiv(:)
      complex(kind=kind(0.d0)), allocatable :: BMatT(:, :)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND2))
      ipiv = 0
      
		allocate(BMatT(ND2, ND1))
      BMatT = dcmplx(0.0d0, 0.0d0)
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) LU decomposition for AMat matrix ____________________________
!____________________________________________________________________________________ 
      call ZGETRF(ND2, ND2, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetZ_Left_LogDet: ZGETRF error with Info = ', I4)") Info
         write(*, "('                         The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                         Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetZ_Left_LogDet: ZGETRF error with Info = ', I4)") Info
         write(*, "('                         LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                         Will stop running here!!!')")
         stop
      end if
!____________________________________________________________________________________	  
!__________________ (1) Log(Determinant) of AMat(ND2, ND2) matrix ___________________
!____________________________________________________________________________________
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, ND2, +1
         if( ipiv(I1) == I1 ) then
            LogzDet = LogzDet + log( + AMat(I1, I1) )
         else
            LogzDet = LogzDet + log( - AMat(I1, I1) )
         end if  
      enddo
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat^T * XMat^T = BMat^T equation set ________________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) Obtain the transpose of BMat matrix _________________________
!____________________________________________________________________________________     
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND2, +1
         do I2 = 1, ND1, +1
            BMatT(I1, I2) = BMat(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!____________________________________________________________________________________	  
!__________________ (1) Solve AMat^T * XMat^T = BMat^T ______________________________
!____________________________________________________________________________________ 
      call ZGETRS("T", ND2, ND1, AMat, LDA, ipiv, BMatT, ND2, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetZ_Left_LogDet: ZGETRS error with Info = ', I4)") Info
         write(*, "('                         Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Output the result XMat matrx to BMat matrix ________________________
!________________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND1, +1
         do I2 = 1, ND2, +1
            BMat(I1, I2) = BMatT(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL     
!________________________________________________________________________________________________	  
!________________________ 4. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(BMatT)) deallocate(BMatT)
		
   end subroutine SlvLnEqSetZ_Left_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetZ_Rght_NoDet(ND1, ND2, AMat, LDA, BMat, LDB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetZ_Rght_NoDet(ND1, ND2, AMat, LDA, BMat, LDB)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set AMat * XMat = BMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2). 
!     First , partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Second, Solve AMat * XMat = BMat;
!     Third , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
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
		integer ND1, ND2
      integer LDA, LDB
      complex(kind=kind(0.d0)) AMat(LDA, *)
      complex(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Info
      integer, allocatable :: ipiv(:)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND1))
      ipiv = 0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
      call ZGETRF(ND1, ND1, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetZ_Rght_NoDet: ZGETRF error with Info = ', I4)") Info
         write(*, "('                        The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                        Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetZ_Rght_NoDet: ZGETRF error with Info = ', I4)") Info
         write(*, "('                        LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                        Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat * XMat = BMat equation set ______________________________
!________________________________________________________________________________________________
      call ZGETRS("N", ND1, ND2, AMat, LDA, ipiv, BMat, LDB, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetZ_Rght_NoDet: ZGETRS error with Info = ', I4)") Info
         write(*, "('                        Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		
   end subroutine SlvLnEqSetZ_Rght_NoDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetZ_Rght_Det(ND1, ND2, AMat, LDA, BMat, LDB, zDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetZ_Rght_Det(ND1, ND2, AMat, LDA, BMat, LDB, zDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set AMat * XMat = BMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2). 
!     First , partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Second, Solve AMat * XMat = BMat;
!     Third , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
!            zMat --> Determinant of AMat(ND2, ND2) matrix.
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
		integer ND1, ND2
      integer LDA, LDB
      complex(kind=kind(0.d0)) zDet
      complex(kind=kind(0.d0)) AMat(LDA, *)
      complex(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, Info
      integer, allocatable :: ipiv(:)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND1))
      ipiv = 0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) LU decomposition for AMat matrix ____________________________
!____________________________________________________________________________________ 
      call ZGETRF(ND1, ND1, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetZ_Rght_Det: ZGETRF error with Info = ', I4)") Info
         write(*, "('                      The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                      Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetZ_Rght_Det: ZGETRF error with Info = ', I4)") Info
         write(*, "('                      LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                      Will stop running here!!!')")
         stop
      end if
!____________________________________________________________________________________	  
!__________________ (1) Determinant of AMat(ND1, ND1) matrix ________________________
!____________________________________________________________________________________
      zDet = dcmplx(1.0d0, 0.0d0)
      do I1 = 1, ND1, +1
         if(ipiv(I1) == I1) then
            zDet = zDet * ( + AMat(I1, I1) )
         else
            zDet = zDet * ( - AMat(I1, I1) )
         end if
      enddo
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat * XMat = BMat equation set ______________________________
!________________________________________________________________________________________________
      call ZGETRS("N", ND1, ND2, AMat, LDA, ipiv, BMat, LDB, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetZ_Rght_Det: ZGETRS error with Info = ', I4)") Info
         write(*, "('                      Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		
   end subroutine SlvLnEqSetZ_Rght_Det
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetZ_Rght_LogDet(ND1, ND2, AMat, LDA, BMat, LDB, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetZ_Rght_LogDet(ND1, ND2, AMat, LDA, BMat, LDB, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set AMat * XMat = BMat, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2). 
!     First , partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Second, Solve AMat * XMat = BMat;
!     Third , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
!            LogzDet --> Log(Determinant) of AMat(ND1, ND1) matrix.
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
		integer ND1, ND2
      integer LDA, LDB
      complex(kind=kind(0.d0)) LogzDet
      complex(kind=kind(0.d0)) AMat(LDA, *)
      complex(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, Info
      integer, allocatable :: ipiv(:)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND1))
      ipiv = 0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) LU decomposition for AMat matrix ____________________________
!____________________________________________________________________________________ 
      call ZGETRF(ND1, ND1, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetZ_Rght_LogDet: ZGETRF error with Info = ', I4)") Info
         write(*, "('                         The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                         Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetZ_Rght_LogDet: ZGETRF error with Info = ', I4)") Info
         write(*, "('                         LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                         Will stop running here!!!')")
         stop
      end if
!____________________________________________________________________________________	  
!__________________ (1) Log(Determinant) of AMat(ND1, ND1) matrix ________________________
!____________________________________________________________________________________
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, ND1, +1
         if(ipiv(I1) == I1) then
            LogzDet = LogzDet + log( + AMat(I1, I1) )
         else
            LogzDet = LogzDet + log( - AMat(I1, I1) )
         end if
      enddo
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat * XMat = BMat equation set ______________________________
!________________________________________________________________________________________________
      call ZGETRS("N", ND1, ND2, AMat, LDA, ipiv, BMat, LDB, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetZ_Rght_LogDet: ZGETRS error with Info = ', I4)") Info
         write(*, "('                         Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		
   end subroutine SlvLnEqSetZ_Rght_LogDet
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
	subroutine SlvLnEqSetR_Left_NoDet(ND1, ND2, AMat, LDA, BMat, LDB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetR_Left_NoDet(ND1, ND2, AMat, LDA, BMat, LDB)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set XMat * AMat = BMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2). 
!     First , take a transpose (or conjugatetranspose) as AMat^T * XMat^T = BMat^T (or AMat^+ * XMat^+ = BMat^+);
!     Second, partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Third , Solve AMat^T from AMat^T * XMat^T = BMat^T (or XMat^+ from AMat^+ * XMat^+ = BMat^+);
!     Forth , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
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
		integer ND1, ND2
      integer LDA, LDB
      real(kind=kind(0.d0)) AMat(LDA, *)
      real(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, I2
      integer Info
      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: BMatT(:, :)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND2))
      ipiv = 0
      
		allocate(BMatT(ND2, ND1))
      BMatT = 0.0d0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
      call DGETRF(ND2, ND2, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetR_Left_NoDet: DGETRF error with Info = ', I4)") Info
         write(*, "('                        The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                        Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetR_Left_NoDet: DGETRF error with Info = ', I4)") Info
         write(*, "('                        LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                        Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat^T * XMat^T = BMat^T equation set ________________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) Obtain the transpose of BMat matrix _________________________
!____________________________________________________________________________________     
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND2, +1
         do I2 = 1, ND1, +1
            BMatT(I1, I2) = BMat(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!____________________________________________________________________________________	  
!__________________ (1) Solve AMat^T * XMat^T = BMat^T ______________________________
!____________________________________________________________________________________ 
      call DGETRS("T", ND2, ND1, AMat, LDA, ipiv, BMatT, ND2, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetR_Left_NoDet: DGETRS error with Info = ', I4)") Info
         write(*, "('                        Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Output the result XMat matrx to BMat matrix ________________________
!________________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND1, +1
         do I2 = 1, ND2, +1
            BMat(I1, I2) = BMatT(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL     
!________________________________________________________________________________________________	  
!________________________ 4. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(BMatT)) deallocate(BMatT)
		
   end subroutine SlvLnEqSetR_Left_NoDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetR_Left_Det(ND1, ND2, AMat, LDA, BMat, LDB, dDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetR_Left_Det(ND1, ND2, AMat, LDA, BMat, LDB, dDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set XMat * AMat = BMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2). 
!     First , take a transpose (or conjugatetranspose) as AMat^T * XMat^T = BMat^T (or AMat^+ * XMat^+ = BMat^+);
!     Second, partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Third , Solve AMat^T from AMat^T * XMat^T = BMat^T (or XMat^+ from AMat^+ * XMat^+ = BMat^+);
!     Forth , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
!            dDet --> Determinant of AMat(ND2, ND2) matrix.
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
		integer ND1, ND2
      integer LDA, LDB
      real(kind=kind(0.d0)) dDet
      real(kind=kind(0.d0)) AMat(LDA, *)
      real(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, I2
      integer Info
      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: BMatT(:, :)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND2))
      ipiv = 0
      
		allocate(BMatT(ND2, ND1))
      BMatT = 0.0d0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) LU decomposition for AMat matrix ____________________________
!____________________________________________________________________________________ 
      call DGETRF(ND2, ND2, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetR_Left_Det: DGETRF error with Info = ', I4)") Info
         write(*, "('                      The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                      Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetR_Left_Det: DGETRF error with Info = ', I4)") Info
         write(*, "('                      LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                      Will stop running here!!!')")
         stop
      end if
!____________________________________________________________________________________	  
!__________________ (1) Determinant of AMat(ND2, ND2) matrix ________________________
!____________________________________________________________________________________
      dDet = 1.0d0
      do I1 = 1, ND2, +1
         if(ipiv(I1) == I1) then
            dDet = dDet * ( + AMat(I1, I1) )
         else
            dDet = dDet * ( - AMat(I1, I1) )
         end if
      enddo
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat^T * XMat^T = BMat^T equation set ________________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) Obtain the transpose of BMat matrix _________________________
!____________________________________________________________________________________     
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND2, +1
         do I2 = 1, ND1, +1
            BMatT(I1, I2) = BMat(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!____________________________________________________________________________________	  
!__________________ (1) Solve AMat^T * XMat^T = BMat^T ______________________________
!____________________________________________________________________________________ 
      call DGETRS("T", ND2, ND1, AMat, LDA, ipiv, BMatT, ND2, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetR_Left_Det: DGETRS error with Info = ', I4)") Info
         write(*, "('                      Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Output the result XMat matrx to BMat matrix ________________________
!________________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND1, +1
         do I2 = 1, ND2, +1
            BMat(I1, I2) = BMatT(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL    
!________________________________________________________________________________________________	  
!________________________ 4. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(BMatT)) deallocate(BMatT)
		
   end subroutine SlvLnEqSetR_Left_Det
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetR_Left_LogDet(ND1, ND2, AMat, LDA, BMat, LDB, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetR_Left_LogDet(ND1, ND2, AMat, LDA, BMat, LDB, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set XMat * AMat = BMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve XMat(ND1, ND2) * AMat(ND2, ND2) = BMat(ND1, ND2). 
!     First , take a transpose (or conjugatetranspose) as AMat^T * XMat^T = BMat^T (or AMat^+ * XMat^+ = BMat^+);
!     Second, partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Third , Solve AMat^T from AMat^T * XMat^T = BMat^T (or XMat^+ from AMat^+ * XMat^+ = BMat^+);
!     Forth , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND2, ND2), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
!            LogzDet --> Log(Determinant) of AMat(ND2, ND2) matrix.
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
		integer ND1, ND2
      integer LDA, LDB
      complex(kind=kind(0.d0)) LogzDet
      real(kind=kind(0.d0)) AMat(LDA, *)
      real(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1, I2
      integer Info
      integer, allocatable :: ipiv(:)
      real(kind=kind(0.d0)), allocatable :: BMatT(:, :)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND2))
      ipiv = 0
      
		allocate(BMatT(ND2, ND1))
      BMatT = 0.0d0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) LU decomposition for AMat matrix ____________________________
!____________________________________________________________________________________ 
      call DGETRF(ND2, ND2, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetR_Left_LogDet: DGETRF error with Info = ', I4)") Info
         write(*, "('                         The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                         Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetR_Left_LogDet: DGETRF error with Info = ', I4)") Info
         write(*, "('                         LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                         Will stop running here!!!')")
         stop
      end if
!____________________________________________________________________________________	  
!__________________ (1) Log(Determinant) of AMat(ND2, ND2) matrix ___________________
!____________________________________________________________________________________
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, ND2, +1
         if( ipiv(I1) == I1 ) then
            LogzDet = LogzDet + log( + dcmplx(AMat(I1, I1), 0.0d0) )
         else
            LogzDet = LogzDet + log( - dcmplx(AMat(I1, I1), 0.0d0) )
         end if  
      enddo
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat^T * XMat^T = BMat^T equation set ________________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) Obtain the transpose of BMat matrix _________________________
!____________________________________________________________________________________     
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND2, +1
         do I2 = 1, ND1, +1
            BMatT(I1, I2) = BMat(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL
!____________________________________________________________________________________	  
!__________________ (1) Solve AMat^T * XMat^T = BMat^T ______________________________
!____________________________________________________________________________________ 
      call DGETRS("T", ND2, ND1, AMat, LDA, ipiv, BMatT, ND2, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetR_Left_LogDet: DGETRS error with Info = ', I4)") Info
         write(*, "('                         Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Output the result XMat matrx to BMat matrix ________________________
!________________________________________________________________________________________________
   !$OMP PARALLEL &
   !$OMP PRIVATE(I1, I2)
   !$OMP DO
      do I1 = 1, ND1, +1
         do I2 = 1, ND2, +1
            BMat(I1, I2) = BMatT(I2, I1)
         enddo
      enddo
   !$OMP END DO
   !$OMP END PARALLEL     
!________________________________________________________________________________________________	  
!________________________ 4. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		if(allocated(BMatT)) deallocate(BMatT)
		
   end subroutine SlvLnEqSetR_Left_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   
   
   

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetR_Rght_NoDet(ND1, ND2, AMat, LDA, BMat, LDB)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetR_Rght_NoDet(ND1, ND2, AMat, LDA, BMat, LDB)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set AMat * XMat = BMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2). 
!     First , partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Second, Solve AMat * XMat = BMat;
!     Third , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
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
		integer ND1, ND2
      integer LDA, LDB
      real(kind=kind(0.d0)) AMat(LDA, *)
      real(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer Info
      integer, allocatable :: ipiv(:)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND1))
      ipiv = 0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
      call DGETRF(ND1, ND1, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetR_Rght_NoDet: DGETRF error with Info = ', I4)") Info
         write(*, "('                        The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                        Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetR_Rght_NoDet: DGETRF error with Info = ', I4)") Info
         write(*, "('                        LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                        Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat * XMat = BMat equation set ______________________________
!________________________________________________________________________________________________
      call DGETRS("N", ND1, ND2, AMat, LDA, ipiv, BMat, LDB, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetR_Rght_NoDet: DGETRS error with Info = ', I4)") Info
         write(*, "('                        Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		
   end subroutine SlvLnEqSetR_Rght_NoDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetR_Rght_Det(ND1, ND2, AMat, LDA, BMat, LDB, dDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetR_Rght_Det(ND1, ND2, AMat, LDA, BMat, LDB, dDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set AMat * XMat = BMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2). 
!     First , partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Second, Solve AMat * XMat = BMat;
!     Third , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
!            zMat --> Determinant of AMat(ND2, ND2) matrix.
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
		integer ND1, ND2
      integer LDA, LDB
      real(kind=kind(0.d0)) dDet
      real(kind=kind(0.d0)) AMat(LDA, *)
      real(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1
      integer Info
      integer, allocatable :: ipiv(:)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND1))
      ipiv = 0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) LU decomposition for AMat matrix ____________________________
!____________________________________________________________________________________ 
      call DGETRF(ND1, ND1, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetR_Rght_Det: DGETRF error with Info = ', I4)") Info
         write(*, "('                      The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                      Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetR_Rght_Det: DGETRF error with Info = ', I4)") Info
         write(*, "('                      LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                      Will stop running here!!!')")
         stop
      end if
!____________________________________________________________________________________	  
!__________________ (1) Determinant of AMat(ND1, ND1) matrix ________________________
!____________________________________________________________________________________
      dDet = 1.0d0
      do I1 = 1, ND1, +1
         if(ipiv(I1) == I1) then
            dDet = dDet * ( + AMat(I1, I1) )
         else
            dDet = dDet * ( - AMat(I1, I1) )
         end if
      enddo
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat * XMat = BMat equation set ______________________________
!________________________________________________________________________________________________
      call DGETRS("N", ND1, ND2, AMat, LDA, ipiv, BMat, LDB, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetR_Rght_Det: DGETRS error with Info = ', I4)") Info
         write(*, "('                      Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		
   end subroutine SlvLnEqSetR_Rght_Det
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	subroutine SlvLnEqSetR_Rght_LogDet(ND1, ND2, AMat, LDA, BMat, LDB, LogzDet)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! PROGRAM:  SlvLnEqSetR_Rght_LogDet(ND1, ND2, AMat, LDA, BMat, LDB, LogzDet)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to solve the linear equation set AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2), and
!                 finally write XMat into BMat matrix.
! KEYWORDS: Solve linear equation set AMat * XMat = BMat, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Solve AMat(ND1, ND1) * XMat(ND1, ND2) = BMat(ND1, ND2). 
!     First , partial-pivoted LU decomposition for AMat as AMat = P * L * U;
!     Second, Solve AMat * XMat = BMat;
!     Third , BMat = XMat.
!
!     Input: ND1  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            ND2  --> Dimension as AMat(ND1, ND1), BMat(ND1, ND2);
!            AMat --> Input AMat square matrix;
!            LDA  --> Leading dimension of input AMat matrix;
!            BMat --> Input BMat matrix and output XMat matrix; 
!            LDB  --> Leading dimension of input BMat matrix.
!
!     Outpt: BMat --> Input BMat matrix and output XMat matrix; 
!            LogzDet --> Log(Determinant) of AMat(ND2, ND2) matrix.
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
		integer ND1, ND2
      integer LDA, LDB
      complex(kind=kind(0.d0)) LogzDet
      real(kind=kind(0.d0)) AMat(LDA, *)
      real(kind=kind(0.d0)) BMat(LDB, *)
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer I1
      integer Info
      integer, allocatable :: ipiv(:)
!______________________________________________________________________________________________________________	  
!____________________________ Main calculations of solving XMat * AMat = BMat _________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________________ 0. Allocate necessary matrices and vectors ____________________________
!________________________________________________________________________________________________
      allocate(ipiv(ND1))
      ipiv = 0
!________________________________________________________________________________________________	  
!________________________ 1. Perform AMat = P*L*U decomposition for AMat matrix _________________
!________________________________________________________________________________________________
!____________________________________________________________________________________	  
!__________________ (0) LU decomposition for AMat matrix ____________________________
!____________________________________________________________________________________ 
      call DGETRF(ND1, ND1, AMat, LDA, ipiv, Info)
      if(Info < 0) then
         write(*, "('SlvLnEqSetR_Rght_LogDet: DGETRF error with Info = ', I4)") Info
         write(*, "('                         The ', I4, '-th input parameter is illegal!')") Info
         write(*, "('                         Will stop running here!!!')") 
         stop
      else if(Info > 0) then
         write(*, "('SlvLnEqSetR_Rght_LogDet: DGETRF error with Info = ', I4)") Info
         write(*, "('                         LU decomposition is done. But the ', I4, '-th diagonal element is zero!')") Info
         write(*, "('                         Will stop running here!!!')")
         stop
      end if
!____________________________________________________________________________________	  
!__________________ (1) Log(Determinant) of AMat(ND1, ND1) matrix ________________________
!____________________________________________________________________________________
      LogzDet = dcmplx(0.0d0, 0.0d0)
      do I1 = 1, ND1, +1
         if(ipiv(I1) == I1) then
            LogzDet = LogzDet + log( + dcmplx(AMat(I1, I1), 0.0d0) )
         else
            LogzDet = LogzDet + log( - dcmplx(AMat(I1, I1), 0.0d0) )
         end if
      enddo
!________________________________________________________________________________________________	  
!________________________ 2. Solve AMat * XMat = BMat equation set ______________________________
!________________________________________________________________________________________________
      call DGETRS("N", ND1, ND2, AMat, LDA, ipiv, BMat, LDB, Info)
      if(Info /= 0) then
         write(*, "('SlvLnEqSetR_Rght_LogDet: DGETRS error with Info = ', I4)") Info
         write(*, "('                         Will stop running here!!!')")
         stop
      end if
!________________________________________________________________________________________________	  
!________________________ 3. Deallocate the allocated matrices and vectors ______________________
!________________________________________________________________________________________________
		if(allocated(ipiv)) deallocate(ipiv)
		
   end subroutine SlvLnEqSetR_Rght_LogDet
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   