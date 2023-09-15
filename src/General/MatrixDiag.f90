!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A few subroutines used to diagonalize matrices and obtain the eigenvalues and eigenvaectors based on LAPACK. 
!              Both real and complex versions are presented.
!           
! COMMENT: Based on subroutines in LAPACK.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!   MatrDiagR1(LDim, NDim, AMat, EVal)                   --> Eigenvalues                  for real symmetric matrix   
!   MatrDiagR2(LDim, NDim, AMat, EVal, EVec)             --> Eigenvalues and eigenvectors for real symmetric matrix
!   MatrDiagR3(LDim, NDim, AMat, EVlr, EVli)             --> Eigenvalues                  for general real matrix
!   MatrDiagR4(LDim, NDim, AMat, EVlr, EVli, EVcr, EVci) --> Eigenvalues and eigenvectors for general real matrix
!   MatrDiagZ1(LDim, NDim, AMat, EVal)                   --> Eigenvalues                  for complex Hermitian matrix   
!   MatrDiagZ2(LDim, NDim, AMat, EVal, EVec)             --> Eigenvalues and eigenvectors for complex Hermitian matrix   
!   MatrDiagZ3(LDim, NDim, AMat, EVal)                   --> Eigenvalues                  for general complex matrix 
!   MatrDiagZ4(LDim, NDim, AMat, EVal, EVec)             --> Eigenvalues and eigenvectors for general complex matrix
!  
! PURPOSE: Get the Eigenvalues and Eigenvectors for general Matrices.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   
   
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
!################################################# For Complex Version ##################################################
!################################################# For Complex Version ##################################################
!################################################# For Complex Version ##################################################
!################################################# For Complex Version ##################################################
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrDiagZ1(LDim, NDim, AMat, EVal)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrDiagZ1(LDim, NDim, AMat, EVal)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate only the eigenvalues of complex Hermitian matrix. 
! KEYWORDS: Eigenvalues of complex Hermitian matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrix Diaogonalization of complex Hermitian matrix. 
!
!     Input: LDim --> Leading dimension of AMat matrix as AMat(LDim, NDim);
!            NDim --> Other   dimension of AMat matrix as AMat(LDim, NDim);
!            AMat --> Input complex Hermitian matrix;
!             
!     Outpt: EVal --> The eigenvalues of AMat matrix.
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
		integer LDim                    ! Dimension of A matrix A(N, M)
		integer NDim                    ! Dimension of A matrix A(N, M)
		complex(kind=kind(0.d0)) AMat(LDim, NDim)    ! Input AMat matrix
		real(kind=kind(0.d0)) EVal(NDim)             ! Output eigenvalues
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer Istat
      integer Info
      integer LWork
      integer LrWork
      
      real(kind=kind(0.d0)), allocatable :: rWork(:)
      complex(kind=kind(0.d0)), allocatable :: Work(:)  
      complex(kind=kind(0.d0)), allocatable :: TMat(:, :)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LWork  = 2 * NDim - 1
      LrWork = 3 * NDim - 2
      
      allocate( Work(LWork     ), stat=istat)
      allocate(rWork(LrWork    ), stat=istat)
      allocate( TMat(LDim, NDim), stat=istat)
      if ( istat /= 0 ) then
         write(*, "('MatrDiagZ1: can not allocate enough memory! ierror = ', I4)") istat
      end if
      
      Eval = 0.0d0
      TMat = AMat
				
		call ZHEEV('N', 'U', Ndim, TMat, Ldim, Eval, Work, LWork, rWork, Info)
		if ( info /= 0 ) then
         write(*, "('MatrDiagZ1: error in lapack subroutine ZHEEV! ierror = ', I4)") info
      end if
      
      if(allocated( Work)) deallocate( Work)
      if(allocated(rWork)) deallocate(rWork)
      if(allocated( TMat)) deallocate( TMat)
      
   end subroutine MatrDiagZ1
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrDiagZ2(LDim, NDim, AMat, EVal, EVec)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrDiagZ2(LDim, NDim, AMat, EVal, EVec)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate only the eigenvalues and eigenvectors of complex Hermitian matrix. 
! KEYWORDS: Eigenvalues and eigenvectors of complex Hermitian matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrix Diaogonalization of complex Hermitian matrix. 
!
!     Input: LDim --> Leading dimension of AMat matrix as AMat(LDim, NDim);
!            NDim --> Other   dimension of AMat matrix as AMat(LDim, NDim);
!            AMat --> Input complex Hermitian matrix;
!             
!     Outpt: EVal --> The eigenvalues  of AMat matrix;
!            EVec --> The eigenvectors of AMat matrix.
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
		integer LDim                    ! Dimension of A matrix A(N, M)
		integer NDim                    ! Dimension of A matrix A(N, M)
		complex(kind=kind(0.d0)) AMat(LDim, NDim)    ! Input AMat matrix
		real(kind=kind(0.d0)) EVal(NDim)             ! Output eigenvalues
      complex(kind=kind(0.d0)) EVec(LDim, NDim)    ! Output real      part of the eigenvector matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer Istat
      integer Info
      integer LWork
      integer LrWork
      
      real(kind=kind(0.d0)), allocatable :: rWork(:)
      complex(kind=kind(0.d0)), allocatable :: Work(:)  
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LWork  = 2 * NDim - 1
      LrWork = 3 * NDim - 2
      
      allocate( Work(LWork     ), stat=istat)
      allocate(rWork(LrWork    ), stat=istat)
      if ( istat /= 0 ) then
         write(*, "('MatrDiagZ2: can not allocate enough memory! ierror = ', I4)") istat
      end if
      
      Eval = 0.0d0
      EVec = AMat
				
		call ZHEEV('V', 'U', Ndim, EVec, Ldim, Eval, Work, LWork, rWork, Info)
		if ( info /= 0 ) then
         write(*, "('MatrDiagZ2: error in lapack subroutine ZHEEV! ierror = ', I4)") info
      end if
      
      if(allocated( Work)) deallocate( Work)
      if(allocated(rWork)) deallocate(rWork)
      
   end subroutine MatrDiagZ2
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrDiagZ3(LDim, NDim, AMat, EVal)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrDiagZ3(LDim, NDim, AMat, EVal)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate only the eigenvalues of general complex matrix. 
! KEYWORDS: Eigenvalues of general complex matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrix Diaogonalization of general complex matrix. 
!
!     Input: LDim --> Leading dimension of AMat matrix as AMat(LDim, NDim);
!            NDim --> Other   dimension of AMat matrix as AMat(LDim, NDim);
!            AMat --> Input general complex matrix;
!             
!     Outpt: EVal --> The eigenvalues of AMat matrix.
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
		integer LDim                    ! Dimension of A matrix A(N, M)
		integer NDim                    ! Dimension of A matrix A(N, M)
		complex(kind=kind(0.d0)) AMat(LDim, NDim)    ! Input AMat matrix
		complex(kind=kind(0.d0)) EVal(NDim)          ! Output eigenvalues
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer Istat
      integer Info
      integer LWork
      
      complex(kind=kind(0.d0)), allocatable :: rWork(:   )
      complex(kind=kind(0.d0)), allocatable ::  Work(:   )  
      complex(kind=kind(0.d0)), allocatable ::  TMat(:, :)
      complex(kind=kind(0.d0)), allocatable ::    vr(:, :)
      complex(kind=kind(0.d0)), allocatable ::    vl(:, :)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LWork  = 2 * NDim
      
      allocate( Work(LWork     ), stat=istat)
      allocate(rWork(LWork     ), stat=istat)
      allocate( TMat(LDim, NDim), stat=istat)
      allocate(   vr(NDim, NDim), stat=istat)
      allocate(   vl(NDim, NDim), stat=istat)
      if ( istat /= 0 ) then
         write(*, "('MatrDiagZ3: can not allocate enough memory! ierror = ', I4)") istat
      end if
      
      Eval = dcmplx(0.0d0, 0.0d0)
      TMat = AMat
				
		call ZGEEV('N', 'N', Ndim, TMat, Ldim, Eval, vl, NDim, vr, NDim, Work, LWork, rWork, Info)
		if ( info /= 0 ) then
         write(*, "('MatrDiagZ3: error in lapack subroutine ZGEEV! ierror = ', I4)") info
      end if
      
      if(allocated( Work)) deallocate( Work)
      if(allocated(rWork)) deallocate(rWork)
      if(allocated( TMat)) deallocate( TMat)
      if(allocated(   vr)) deallocate(   vr)
      if(allocated(   vl)) deallocate(   vl)
      
   end subroutine MatrDiagZ3
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrDiagZ4(LDim, NDim, AMat, EVal, EVec)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrDiagZ4(LDim, NDim, AMat, EVal, EVec)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate only the eigenvalues and eigenvectors of general complex matrix. 
! KEYWORDS: Eigenvalues and eigenvectors of general complex matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrix Diaogonalization of general complex matrix. 
!
!     Input: LDim --> Leading dimension of AMat matrix as AMat(LDim, NDim);
!            NDim --> Other   dimension of AMat matrix as AMat(LDim, NDim);
!            AMat --> Input general complex matrix;
!             
!     Outpt: EVal --> The eigenvalues  of AMat matrix;
!            EVec --> The eigenvectors of AMat matrix.
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
		integer LDim                    ! Dimension of A matrix A(N, M)
		integer NDim                    ! Dimension of A matrix A(N, M)
		complex(kind=kind(0.d0)) AMat(LDim, NDim)    ! Input AMat matrix
		complex(kind=kind(0.d0)) EVal(NDim)          ! Output eigenvalues
      complex(kind=kind(0.d0)) EVec(LDim, NDim)    ! Output eigenvectors
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer Istat
      integer Info
      integer LWork
      
      complex(kind=kind(0.d0)), allocatable :: rWork(:   )
      complex(kind=kind(0.d0)), allocatable ::  Work(:   )  
      complex(kind=kind(0.d0)), allocatable ::    vr(:, :)
      complex(kind=kind(0.d0)), allocatable ::    vl(:, :)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LWork  = 2 * NDim
      
      allocate( Work(LWork     ), stat=istat)
      allocate(rWork(LWork     ), stat=istat)
      allocate(   vr(NDim, NDim), stat=istat)
      allocate(   vl(NDim, NDim), stat=istat)
      if ( istat /= 0 ) then
         write(*, "('MatrDiagZ4: can not allocate enough memory! ierror = ', I4)") istat
      end if
      
      Eval = dcmplx(0.0d0, 0.0d0)
      EVec = AMat
				
		call ZGEEV('N', 'V', Ndim, EVec, Ldim, Eval, vl, NDim, vr, NDim, Work, LWork, rWork, Info)
		if ( info /= 0 ) then
         write(*, "('MatrDiagZ4: error in lapack subroutine ZGEEV! ierror = ', I4)") info
      end if
      
      EVec = vr
      
      if(allocated( Work)) deallocate( Work)
      if(allocated(rWork)) deallocate(rWork)
      if(allocated(   vr)) deallocate(   vr)
      if(allocated(   vl)) deallocate(   vl)
      
   end subroutine MatrDiagZ4
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################
!################################################# For Real Version #####################################################
!################################################# For Real Version #####################################################
!################################################# For Real Version #####################################################
!################################################# For Real Version #####################################################
!########################################################################################################################
!########################################################################################################################
!########################################################################################################################   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrDiagR1(LDim, NDim, AMat, EVal)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrDiagR1(LDim, NDim, AMat, EVal)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate only the eigenvalues of real symmetric matrix. 
! KEYWORDS: Eigenvalues of real symmetric matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrix Diaogonalization of real symmetric matrix. 
!
!     Input: LDim --> Leading dimension of AMat matrix as AMat(LDim, NDim);
!            NDim --> Other   dimension of AMat matrix as AMat(LDim, NDim);
!            AMat --> Input real symmetric matrix;
!             
!     Outpt: EVal --> The eigenvalues of AMat matrix.
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
		integer LDim                 ! Dimension of A matrix A(N, M)
		integer NDim                 ! Dimension of A matrix A(N, M)
		real(kind=kind(0.d0)) AMat(LDim, NDim)    ! Input A matrix
		real(kind=kind(0.d0)) EVal(NDim)          ! Output eigenvalues
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer Istat
      integer Info
      integer LWork
      
      real(kind=kind(0.d0)), allocatable :: Work(:)  
      real(kind=kind(0.d0)), allocatable :: TMat(:, :)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LWork = 3 * NDim - 1
      
      allocate(Work(LWork     ), stat=istat)
      allocate(TMat(LDim, NDim), stat=istat)
      if ( istat /= 0 ) then
         write(*, "('MatrDiagR1: can not allocate enough memory! ierror = ', I4)") istat
      end if
      
      Eval = 0.0d0
      TMat = AMat
				
		call DSYEV('N', 'U', Ndim, TMat, Ldim, Eval, Work, LWork, Info)
		if ( info /= 0 ) then
         write(*, "('MatrDiagR1: error in lapack subroutine DSYEV! ierror = ', I4)") info
      end if
      
      if(allocated(Work)) deallocate(Work)
      if(allocated(TMat)) deallocate(TMat)
      
   end subroutine MatrDiagR1
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrDiagR2(LDim, NDim, AMat, EVal, EVec)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrDiagR2(LDim, NDim, AMat, EVal, EVec)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate all the eigenvalues and eigenvectors of real symmetric matrix. 
! KEYWORDS: Eigenvalues and eigenvectors of real symmetric matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrix Diaogonalization of real symmetric matrix. 
!
!     Input: LDim --> Leading dimension of AMat matrix as AMat(LDim, NDim);
!            NDim --> Other   dimension of AMat matrix as AMat(LDim, NDim);
!            AMat --> Input real symmetric matrix;
!             
!     Outpt: EVal --> The eigenvalues of AMat matrix.
!            EVec --> The eigenvalues of AMat matrix.
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
		integer LDim                 ! Dimension of A matrix A(N, M)
		integer NDim                 ! Dimension of A matrix A(N, M)
		real(kind=kind(0.d0)) AMat(LDim, NDim)    ! Input AMat matrix
		real(kind=kind(0.d0)) EVal(NDim)          ! Output eigenvalues
      real(kind=kind(0.d0)) EVec(LDim, NDim)    ! Eigenvector matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer Istat
      integer Info
      integer LWork
      
      real(kind=kind(0.d0)), allocatable :: Work(:)  
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LWork = 3 * NDim - 1
      
      allocate(Work(LWork     ), stat=istat)
      if ( istat /= 0 ) then
         write(*, "('MatrDiagR2: can not allocate enough memory! ierror = ', I4)") istat
      end if
      
      Eval = 0.0d0
      Evec = AMat
				
		call DSYEV('V', 'U', Ndim, EVec, Ldim, Eval, Work, LWork, Info)
		if ( info /= 0 ) then
         write(*, "('MatrDiagR2: error in lapack subroutine DSYEV! ierror = ', I4)") Info
      end if
      
      if(allocated(Work)) deallocate(Work)
      
   end subroutine MatrDiagR2
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrDiagR3(LDim, NDim, AMat, EVlr, EVli)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrDiagR3(LDim, NDim, AMat, EVlr, EVli)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate only the eigenvalues of general real matrix. 
! KEYWORDS: Eigenvalues of general real matrix. 
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrix Diaogonalization of general real matrix. 
!
!     Input: LDim --> Leading dimension of AMat matrix as AMat(LDim, NDim);
!            NDim --> Other   dimension of AMat matrix as AMat(LDim, NDim);
!            AMat --> Input real symmetric matrix;
!             
!     Outpt: EVlr --> The real      part of eigenvalues of AMat matrix;
!            EVli --> The imaginary part of eigenvalues of AMat matrix.
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
		integer LDim                 ! Dimension of A matrix A(N, M)
		integer NDim                 ! Dimension of A matrix A(N, M)
		real(kind=kind(0.d0)) AMat(LDim, NDim)    ! Input A matrix
		real(kind=kind(0.d0)) EVlr(NDim)          ! Output real      part of the eigenvalues
      real(kind=kind(0.d0)) EVli(NDim)          ! Output imajinary part of the eigenvalues
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer Istat
      integer Info
      integer LWork
      
      real(kind=kind(0.d0)), allocatable :: Work(:)      ! workspace array
      real(kind=kind(0.d0)), allocatable :: TMat(:, :)   ! workspace array, used to store AMat matrix
      real(kind=kind(0.d0)), allocatable :: wr(:)        ! real      part of eigenvalues
      real(kind=kind(0.d0)), allocatable :: wi(:)        ! imaginary part of eigenvalues
      real(kind=kind(0.d0)), allocatable :: vr(:, :)     ! real      part of eigenvectors
      real(kind=kind(0.d0)), allocatable :: vl(:, :)     ! imaginary part of eigenvectors
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LWork = 4 * NDim
      
      allocate(Work(LWork     ), stat=istat)
      allocate(TMat(LDim, NDim), stat=istat)
      allocate(  wr(NDim      ), stat=istat)
      allocate(  wi(NDim      ), stat=istat)
      allocate(  vr(NDim, NDim), stat=istat)
      allocate(  vl(NDim, NDim), stat=istat)
      if ( istat /= 0 ) then
         write(*, "('MatrDiagR3: can not allocate enough memory! ierror = ', I4)") istat
      end if
      
      EVlr = 0.0d0
      EVli = 0.0d0
      TMat = AMat
				
		call DGEEV('N', 'N', Ndim, TMat, Ldim, wr, wi, vl, Ndim, vr, Ndim, Work, LWork, Info)
		if ( Info /= 0 ) then
         write(*, "('MatrDiagR3: error in lapack subroutine DGEEV! ierror = ', I4)") Info
      end if
      
      EVlr(1:NDim) = wr(1:NDim)
      EVli(1:NDim) = wi(1:NDim)
      
      if(allocated(Work)) deallocate(Work)
      if(allocated(TMat)) deallocate(TMat)
      if(allocated(wr  )) deallocate(wr  )
      if(allocated(wi  )) deallocate(wi  )
      if(allocated(vr  )) deallocate(vr  )
      if(allocated(vl  )) deallocate(vl  )
      
   end subroutine MatrDiagR3
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine MatrDiagR4(LDim, NDim, AMat, EVlr, EVli, EVcr, EVci)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrDiagR4(LDim, NDim, AMat, EVlr, EVli, EVcr, EVci)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to calculate only the eigenvalues and eigenvectors of general real matrix. 
! KEYWORDS: Eigenvalues and eigenvectors of general real matrix. 
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     Matrix Diaogonalization of general real matrix. 
!
!     Input: LDim --> Leading dimension of AMat matrix as AMat(LDim, NDim);
!            NDim --> Other   dimension of AMat matrix as AMat(LDim, NDim);
!            AMat --> Input real symmetric matrix;
!             
!     Outpt: EVlr --> The real      part of eigenvalues  of AMat matrix;
!            EVli --> The imaginary part of eigenvalues  of AMat matrix;
!            EVcr --> The real      part of eigenvectors of AMat matrix;
!            EVci --> The imaginary part of eigenvectors of AMat matrix.
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
		integer LDim                 ! Dimension of A matrix A(N, M)
		integer NDim                 ! Dimension of A matrix A(N, M)
		real(kind=kind(0.d0)) AMat(LDim, NDim)    ! Input A matrix
		real(kind=kind(0.d0)) EVlr(NDim)          ! Output real      part of the eigenvalues
      real(kind=kind(0.d0)) EVli(NDim)          ! Output imajinary part of the eigenvalues
      real(kind=kind(0.d0)) EVcr(LDim, NDim)    ! Output real      part of the eigenvector matrix
      real(kind=kind(0.d0)) EVci(LDim, NDim)    ! Output imajinary part of the eigenvector matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer Istat
      integer Info
      integer LWork
      
      real(kind=kind(0.d0)), allocatable :: Work(:)      ! workspace array
      real(kind=kind(0.d0)), allocatable :: TMat(:, :)   ! workspace array, used to store AMat matrix
      real(kind=kind(0.d0)), allocatable :: wr(:)        ! real      part of eigenvalues
      real(kind=kind(0.d0)), allocatable :: wi(:)        ! imaginary part of eigenvalues
      real(kind=kind(0.d0)), allocatable :: vr(:, :)     ! real      part of eigenvectors
      real(kind=kind(0.d0)), allocatable :: vl(:, :)     ! imaginary part of eigenvectors
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LWork = 4 * NDim
      
      allocate(Work(LWork     ), stat=istat)
      allocate(TMat(LDim, NDim), stat=istat)
      allocate(  wr(NDim      ), stat=istat)
      allocate(  wi(NDim      ), stat=istat)
      allocate(  vr(NDim, NDim), stat=istat)
      allocate(  vl(NDim, NDim), stat=istat)
      if ( istat /= 0 ) then
         write(*, "('MatrDiagR4: can not allocate enough memory! ierror = ', I4)") istat
      end if
      
      EVlr = 0.0d0
      EVli = 0.0d0
      TMat = AMat
				
		call DGEEV('N', 'V', Ndim, TMat, Ldim, wr, wi, vl, Ndim, vr, Ndim, Work, LWork, Info)
		if ( Info /= 0 ) then
         write(*, "('MatrDiagR4: error in lapack subroutine DGEEV! ierror = ', I4)") Info
      end if
      
      EVlr(1:NDim) = wr(1:NDim)
      EVli(1:NDim) = wi(1:NDim)
      EVcr(1:NDim, 1:NDim) = vr(1:NDim, 1:NDim)
      EVci(1:NDim, 1:NDim) = vl(1:NDim, 1:NDim)
      
      if(allocated(Work)) deallocate(Work)
      if(allocated(TMat)) deallocate(TMat)
      if(allocated(wr  )) deallocate(wr  )
      if(allocated(wi  )) deallocate(wi  )
      if(allocated(vr  )) deallocate(vr  )
      if(allocated(vl  )) deallocate(vl  )
      
   end subroutine MatrDiagR4
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$