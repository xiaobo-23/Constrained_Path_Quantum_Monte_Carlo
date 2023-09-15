!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A few common subroutines used for sorting of arrays applying the QuckSort Algorithm.
! COMMENT: Used in the energies sorting for different subspaces.  
! AUTHOR:  Yuan-Yao He  -->  770
! DATE:    2020-02-20
! PURPOSE: Different subroutines are introduced as following:
!             
!    QuckSortR(IBgn, IEnd, Arr)              --> Ascending order of Real Array Arr, Real Version;
!    QuckSortI(IBgn, IEnd, Arr)              --> Ascending order of Integer Array Arr, Integer Version;
!    QuckSortIR(IBgn, IEnd, Arr, Brr)        --> Ascending order of Integer Array Arr, Integer-Real Version;
!    QuckSortIC(IBgn, IEnd, Arr, Brr)        --> Ascending order of Integer Array Arr, Integer-Real Version;
!    QuckSortII(IBgn, IEnd, Arr, Brr)        --> Ascending order of Integer Array Arr, Integer-Integer Version;
!    QuckSortRR(IBgn, IEnd, Arr, Brr)        --> Ascending order of Real Array Arr, Real-Real Version;
!    QuckSortRI(IBgn, IEnd, Arr, Brr)        --> Ascending order of Real Array Arr, Real-Integer Version;
!    QuckSortRII(IBgn, IEnd, Arr, Brr, Crr)  --> Ascending order of Real Array Arr, Real-Integer-Integer Version;
!    QuckSortRIR(IBgn, IEnd, Arr, Brr, Crr)  --> Ascending order of Real Array Arr, Real-Integer-Real Version;
!    QuckSortRIIR(IBgn, IEnd, Arr, Brr, Crr) --> Ascending order of Real Array Arr, Real-Integer-Real Version;
!    QuckSortIII(IBgn, IEnd, Arr, Brr, Crr)  --> Ascending order of Integer Array Arr, Integer-Integer-Integer Version;
!    SwapReal(a, b)                          --> Exachange the values of Real a and b;
!    SwapCmplx(a, b)                         --> Exachange the values of complex a and b;
!    SwapInt(a, b)                           --> Exachange the values of Integer a and b;
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortR(IBgn, IEnd, Arr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortR(IBgn, IEnd, Arr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Real Version.
! I/O:     Input:  Arr--The Real array needing to be sorted; IBgn--the initial left number; IEnd--the initial
!                  right number
!          Output: Arr--the array after sorting in ascending order.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find 
!                      Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have 
!                        Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		real(kind=kind(0.d0)) K                        ! The Key value and the temporary quantity
		real(kind=kind(0.d0)) Arr(IBgn : IEnd)         ! The Array needing to be sorted
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if(Arr(L) > K) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                           ! If (R<=L), exit for the do while(.true.)
			if(L /= R) call SwapReal(Arr(L), Arr(R))  ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
		enddo
		if(IBgn /= R) call SwapReal(Arr(IBgn), Arr(R))
		if(R - 1 > IBgn) call QuckSortR(IBgn, R - 1, Arr(IBgn))   ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortR(R + 1, IEnd, Arr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortI(IBgn, IEnd, Arr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortI(IBgn, IEnd, Arr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Integer Version.
! I/O:     Input:  Arr--The Integer array needing to be sorted; IBgn--the initial left number; IEnd--the initial
!                  right number
!          Output: Arr--the array after sorting in ascending order.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we 
!                   find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have 
!                        Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		integer K                        ! The Key value and the temporary quantity
		integer Arr(IBgn : IEnd)         ! The Array needing to be sorted
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                           ! If (R<=L), exit for the do while(.true.)
			if(L /= R) call SwapInt(Arr(L), Arr(R))  ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
		enddo
		if(IBgn /= R) call SwapInt(Arr(IBgn), Arr(R))
		if(R - 1 > IBgn) call QuckSortI(IBgn, R - 1, Arr(IBgn))   ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortI(R + 1, IEnd, Arr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortI
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortIR(IBgn, IEnd, Arr, Brr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortIR(IBgn, IEnd, Arr, Brr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Integer and Real Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Integer array needing to be sorted; 
!                  Brr--The Real array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Integer array after sorting in ascending order;
!                  Brr--The Real array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we 
!                   find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have 
!                        Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		integer K                        ! The Key value and the temporary quantity
		integer Arr(IBgn : IEnd)         ! The Array needing to be sorted
		real(kind=kind(0.d0)) Brr(IBgn : IEnd)         ! The Array needing to be Rearranged
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapInt(Arr(L), Arr(R))
				call SwapReal(Brr(L), Brr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapInt(Arr(IBgn), Arr(R))
			call SwapReal(Brr(IBgn), Brr(R))
		end if
		if(R - 1 > IBgn) call QuckSortIR(IBgn, R - 1, Arr(IBgn), Brr(IBgn))   ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortIR(R + 1, IEnd, Arr(R + 1), Brr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortIR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortIC(IBgn, IEnd, Arr, Brr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortIC(IBgn, IEnd, Arr, Brr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Integer-complex Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Integer array needing to be sorted; 
!                  Brr--The complex array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Integer array after sorting in ascending order;
!                  Brr--The complex array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		integer K                        ! The Key value and the temporary quantity
		integer Arr(IBgn : IEnd)         ! The Array needing to be sorted
		complex(kind=kind(0.d0)) Brr(IBgn : IEnd)     ! The Array needing to be Rearranged
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapInt(Arr(L), Arr(R))
				call SwapCmplx(Brr(L), Brr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapInt(Arr(IBgn), Arr(R))
			call SwapCmplx(Brr(IBgn), Brr(R))
		end if
		if(R - 1 > IBgn) call QuckSortIC( IBgn, R - 1, Arr(IBgn), Brr(IBgn))    ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortIC(R + 1,  IEnd, Arr(R + 1), Brr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortIC
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortRR(IBgn, IEnd, Arr, Brr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortRR(IBgn, IEnd, Arr, Brr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Real and Real Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Real array needing to be sorted; 
!                  Brr--The Real array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Real array after sorting in ascending order;
!                  Brr--The Real array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		real(kind=kind(0.d0)) K                        ! The Key value and the temporary quantity
		real(kind=kind(0.d0)) Arr(IBgn : IEnd)         ! The Array needing to be sorted
		real(kind=kind(0.d0)) Brr(IBgn : IEnd)         ! The Array needing to be Rearranged
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapReal(Arr(L), Arr(R))
				call SwapReal(Brr(L), Brr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapReal(Arr(IBgn), Arr(R))
			call SwapReal(Brr(IBgn), Brr(R))
		end if
		if(R - 1 > IBgn) call QuckSortRR(IBgn, R - 1, Arr(IBgn), Brr(IBgn))   ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortRR(R + 1, IEnd, Arr(R + 1), Brr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortRR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortRI(IBgn, IEnd, Arr, Brr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortRI(IBgn, IEnd, Arr, Brr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Real and Integer Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Real array needing to be sorted; 
!                  Brr--The Integer array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Real array after sorting in ascending order;
!                  Brr--The Integer array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		real(kind=kind(0.d0)) K                        ! The Key value and the temporary quantity
		real(kind=kind(0.d0)) Arr(IBgn : IEnd)         ! The Array needing to be sorted
		integer Brr(IBgn : IEnd)         ! The Array needing to be Rearranged
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapReal(Arr(L), Arr(R))
				call SwapInt(Brr(L), Brr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapReal(Arr(IBgn), Arr(R))
			call SwapInt(Brr(IBgn), Brr(R))
		end if
		if(R - 1 > IBgn) call QuckSortRI(IBgn, R - 1, Arr(IBgn), Brr(IBgn))   ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortRI(R + 1, IEnd, Arr(R + 1), Brr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortRI
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortII(IBgn, IEnd, Arr, Brr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortII(IBgn, IEnd, Arr, Brr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Integer and Integer Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Integer array needing to be sorted; 
!                  Brr--The Integer array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Integer array after sorting in ascending order;
!                  Brr--The Integer array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		integer K                        ! The Key value and the temporary quantity
		integer Arr(IBgn : IEnd)         ! The Array needing to be sorted
		integer Brr(IBgn : IEnd)         ! The Array needing to be Rearranged
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapInt(Arr(L), Arr(R))
				call SwapInt(Brr(L), Brr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapInt(Arr(IBgn), Arr(R))
			call SwapInt(Brr(IBgn), Brr(R))
		end if
		if(R - 1 > IBgn) call QuckSortII(IBgn, R - 1, Arr(IBgn), Brr(IBgn))   ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortII(R + 1, IEnd, Arr(R + 1), Brr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortII
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortRII(IBgn, IEnd, Arr, Brr, Crr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortRII(IBgn, IEnd, Arr, Brr, Crr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Real-Integer-Integer Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Real array needing to be sorted; 
!                  Brr--The Integer array to be rearranged according to the sorting of Arr;
!                  Crr--The Integer array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Real array after sorting in ascending order;
!                  Brr--The Integer array rearranged according to the sorting of Arr.
!                  Crr--The Integer array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		real(kind=kind(0.d0)) K                        ! The Key value and the temporary quantity
		real(kind=kind(0.d0)) Arr(IBgn : IEnd)         ! The Array needing to be sorted
		integer Brr(IBgn : IEnd)         ! The Array needing to be Rearranged
		integer Crr(IBgn : IEnd)         ! The Array needing to be Rearranged
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapReal(Arr(L), Arr(R))
				call SwapInt(Brr(L), Brr(R))
				call SwapInt(Crr(L), Crr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapReal(Arr(IBgn), Arr(R))
			call SwapInt(Brr(IBgn), Brr(R))
			call SwapInt(Crr(IBgn), Crr(R))
		end if
		if(R - 1 > IBgn) call QuckSortRII(IBgn, R - 1, Arr(IBgn), Brr(IBgn), Crr(IBgn))     ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortRII(R + 1, IEnd, Arr(R + 1), Brr(R + 1), Crr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortRII
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortRIR(IBgn, IEnd, Arr, Brr, Crr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortRIR(IBgn, IEnd, Arr, Brr, Crr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Real-Integer-Integer Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Real array needing to be sorted; 
!                  Brr--The Integer array to be rearranged according to the sorting of Arr;
!                  Crr--The Integer array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Real array after sorting in ascending order;
!                  Brr--The Integer array rearranged according to the sorting of Arr.
!                  Crr--The Integer array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		real(kind=kind(0.d0)) K                        ! The Key value and the temporary quantity
		real(kind=kind(0.d0)) Arr(IBgn : IEnd)         ! The Array needing to be sorted
		real(kind=kind(0.d0)) Crr(IBgn : IEnd)         ! The Array needing to be Rearranged
		integer Brr(IBgn : IEnd)         ! The Array needing to be Rearranged	  
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapReal(Arr(L), Arr(R))
				call SwapInt(Brr(L), Brr(R))
				call SwapReal(Crr(L), Crr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapReal(Arr(IBgn), Arr(R))
			call SwapInt(Brr(IBgn), Brr(R))
			call SwapReal(Crr(IBgn), Crr(R))
		end if
		if(R - 1 > IBgn) call QuckSortRIR(IBgn, R - 1, Arr(IBgn), Brr(IBgn), Crr(IBgn))     ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortRIR(R + 1, IEnd, Arr(R + 1), Brr(R + 1), Crr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortRIR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortRIIR(IBgn, IEnd, Arr, Brr, Crr, Drr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortRIIR(IBgn, IEnd, Arr, Brr, Crr, Drr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Real-Integer-Integer Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Real array needing to be sorted; 
!                  Brr--The Integer array to be rearranged according to the sorting of Arr;
!                  Crr--The Integer array to be rearranged according to the sorting of Arr;
!                  Drr--The Real array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Real array after sorting in ascending order;
!                  Brr--The Integer array rearranged according to the sorting of Arr.
!                  Crr--The Integer array rearranged according to the sorting of Arr.
!                  Drr--The Real array to be rearranged according to the sorting of Arr;
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		real(kind=kind(0.d0)) K                        ! The Key value and the temporary quantity
		real(kind=kind(0.d0)) Arr(IBgn : IEnd)         ! The Array needing to be sorted
		real(kind=kind(0.d0)) Drr(IBgn : IEnd)         ! The Array needing to be Rearranged
		integer Brr(IBgn : IEnd)         ! The Array needing to be Rearranged	  
		integer Crr(IBgn : IEnd)         ! The Array needing to be Rearranged	  
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapReal(Arr(L), Arr(R))
				call SwapInt(Brr(L), Brr(R))
				call SwapInt(Crr(L), Crr(R))
				call SwapReal(Drr(L), Drr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapReal(Arr(IBgn), Arr(R))
			call SwapInt(Brr(IBgn), Brr(R))
			call SwapInt(Crr(IBgn), Crr(R))
			call SwapReal(Drr(IBgn), Drr(R))
		end if
		if(R - 1 > IBgn) call QuckSortRIIR(IBgn, R - 1, Arr(IBgn), Brr(IBgn), Crr(IBgn), Drr(IBgn))     ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortRIIR(R + 1, IEnd, Arr(R + 1), Brr(R + 1), Crr(R + 1), Drr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortRIIR
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortIII(IBgn, IEnd, Arr, Brr, Crr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortIII(IBgn, IEnd, Arr, Brr, Crr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in ascending order, Integer-Integer-Integer Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Integer array needing to be sorted; 
!                  Brr--The Integer array to be rearranged according to the sorting of Arr;
!                  Crr--The Integer array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Integer array after sorting in ascending order;
!                  Brr--The Integer array rearranged according to the sorting of Arr.
!                  Crr--The Integer array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		integer K                        ! The Key value and the temporary quantity
		integer Arr(IBgn : IEnd)         ! The Array needing to be sorted
		integer Brr(IBgn : IEnd)         ! The Array needing to be Rearranged
		integer Crr(IBgn : IEnd)         ! The Array needing to be Rearranged
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) > K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) < K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapInt(Arr(L), Arr(R))
				call SwapInt(Brr(L), Brr(R))
				call SwapInt(Crr(L), Crr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapInt(Arr(IBgn), Arr(R))
			call SwapInt(Brr(IBgn), Brr(R))
			call SwapInt(Crr(IBgn), Crr(R))
		end if
		if(R - 1 > IBgn) call QuckSortIII(IBgn, R - 1, Arr(IBgn), Brr(IBgn), Crr(IBgn))     ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortIII(R + 1, IEnd, Arr(R + 1), Brr(R + 1), Crr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortIII
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine SwapReal(a, b)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SwapReal(a, b)
! TYPE:     subroutine
! PURPOSE:  This Subroutine exchange the value of a and b, Real Version.
! KEYWORDS: Exchange Value 
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Small Algorithm.
!
!     Input:  (a,b)--Two Input Real numbers.
!
!     Output: (a,b)--Two Output Real numbers. 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
      implicit none
	  
		real(kind=kind(0.d0)) a, b, c
	  
		c = a
		a = b
		b = c
	  
	end subroutine SwapReal 
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine SwapCmplx(a, b)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SwapCmplx(a, b)
! TYPE:     subroutine
! PURPOSE:  This Subroutine exchange the value of a and b, complex Version.
! KEYWORDS: Exchange Value 
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Small Algorithm.
!
!     Input:  (a,b)--Two Input complex numbers.
!
!     Output: (a,b)--Two Output complex numbers. 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
      implicit none
	  
		complex(kind=kind(0.d0)) a, b, c
	  
		c = a
		a = b
		b = c
	  
	end subroutine SwapCmplx 
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine SwapInt(a, b)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SwapInt(a, b)
! TYPE:     subroutine
! PURPOSE:  This Subroutine exchange the value of a and b, Integer Version.
! KEYWORDS: Exchange Value 
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-20
! DESCRIPTION: Small Algorithm.
!
!     Input:  (a,b)--Two Input Integer numbers.
!
!     Output: (a,b)--Two Output Integer numbers. 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

		implicit none
	  
		integer a, b, c
	  
		c = a
		a = b
		b = c
	  
   end subroutine SwapInt	 
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________  
	recursive subroutine QuckSortRI_Descend(IBgn, IEnd, Arr, Brr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: QuckSortRI_Descend(IBgn, IEnd, Arr, Brr)
! TYPE:    recursive subroutine
! PURPOSE: Using the Quick Sort method to sort an array in descending order, Real and Integer Version.
! I/O:     Input:  (IBgn, IEnd)--The Sorting Range. 
!                  Arr--The Real array needing to be sorted; 
!                  Brr--The Integer array to be rearranged according to the sorting of Arr;
!          Output: Arr--The Real array after sorting in ascending order;
!                  Brr--The Integer array rearranged according to the sorting of Arr.
! VERSION: 2020-02-20
! AUTHER:  Yuan-Yao He
! COMMENT: Complexity--O(n*log_{2}n),when Dim is very big, this method doesn't work well.
! THEORY:  (1)Set a key value first as K=Arr(L),using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1. 
!          (2)Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
!          (3)If p1<p2, interexchange the value of Arr(p1) and Arr(p2). Then go the first step in (1) again.
!          (4)Until we find p1>p2,then we can determine the position of Arr(R), we have Arr(1--(R-1))=<Arr(R),Arr((R+1)--Dim)>=Arr(R)
!          (5)Using the above steps to sort the to array as Arr(1--(R-1)) and Arr((R+1)--Dim).
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
		implicit none
		
		integer L, R, IBgn, IEnd         ! L,R the temporary Left and right number of sorting
		real(kind=kind(0.d0)) K                        ! The Key value and the temporary quantity
		real(kind=kind(0.d0)) Arr(IBgn : IEnd)         ! The Array needing to be sorted
		integer Brr(IBgn : IEnd)         ! The Array needing to be Rearranged
		
		if(IBgn >= IEnd) return          ! if IBgn >= IEnd, there is no need for sorting
		L = IBgn                          ! The initial value for the left number
		R = IEnd + 1                      ! The initial value for the right number
		K = Arr(IBgn)                     ! the initial Key value 
		do while(.true.)
			do while(L < IEnd)             ! Using the values after IBgn to compare with K, until we find Arr(p1)>=K, setting L=p1.
				L = L + 1
				if((Arr(L) < K)) exit
			enddo
			do while(R > IBgn)              ! Using the values before IEnd to compare with K, until we find Arr(p2)<=K, setting IEnd=p2.
				R = R - 1
				if(Arr(R) > K) exit
			enddo
			if(R <= L) exit                ! If (R<=L), exit for the do while(.true.)
			if(L /= R) then                ! If p1<p2, interexchange the value of Arr(p1) and Arr(p2).
				call SwapReal(Arr(L), Arr(R))
				call SwapInt(Brr(L), Brr(R))
			end if
		enddo
		if(IBgn /= R) then
			call SwapReal(Arr(IBgn), Arr(R))
			call SwapInt(Brr(IBgn), Brr(R))
		end if
		if(R - 1 > IBgn) call QuckSortRI_Descend(IBgn, R - 1, Arr(IBgn), Brr(IBgn))   ! sort the to array as Arr(1--(R-1)).
		if(R + 1 < IEnd) call QuckSortRI_Descend(R + 1, IEnd, Arr(R + 1), Brr(R + 1))  ! sort the to array as Arr((R+1)--Dim).

	end subroutine QuckSortRI_Descend
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________ 
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$