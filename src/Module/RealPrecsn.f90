!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A module and a few subroutines used for setting the precision control for all the calculations of the programs.
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!          RealPrecsn--Module to define the related integers and real numbers for all the following calculation 
!                          precision control;
!          MachineRealPrecsn--A small subroutine to calculate the machine precision for the rp kind of real number;
!          RealPrecsnInit--Subroutine used for setting the numbers of calculation precision before the formal calculation.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Module ______________________________________________________________________
!________________________________________________________________________________________________________________________
   module RealPrecsn
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________    
      implicit none
!______________________________________________________________________________________________________________     
!_______________________ Setting to some Accurancy Control and Output Setting _________________________________
!______________________________________________________________________________________________________________      
      integer,parameter :: rp = 8   ! The integer to control the digits of all the calculations
      integer rp_ow                 ! The number of digits for the storing before the decimal point   es(rp_ow,rp_od) --> es25.16
      integer rp_od                 ! The number of digits for the storing before the decimal point
!______________________________________________________________________________________________________________     
!__________________________________ Setting to some Constants In Calculation __________________________________
!______________________________________________________________________________________________________________      
      real(rp) rp_pi                 ! The quantity pi in the precision control rp
      real(rp) rp_one                ! The number 1.0 in the precision control rp
      real(rp) rp_gold               ! The golden section ratio in the precision control rp
      real(rp) rp_huge               ! The largest number in the precision control rp
      real(rp) rp_tiny               ! The smallest number in the precision control rp
      real(rp) rp_prec               ! The machine accurancy in the precision control rp
      real(rp) rp_Eps                ! The small tolerance for some calculations
!______________________________________________________________________________________________________________     
!__________________________________ Some constants ____________________________________________________________
!______________________________________________________________________________________________________________       
      real(rp) rp_Rone
      real(rp) rp_Rzero
      complex(rp) rp_Z_one
      complex(rp) rp_Zzero
!______________________________________________________________________________________________________________     
!__________________________________ MPI related integers for parallel computing _______________________________
!______________________________________________________________________________________________________________
      integer rp_MPI_REAL            ! The real data type used in MPI reduce and other subroutines
      integer rp_MPI_COMPLEX         ! The complex data type used in MPI reduce and other subroutines
      integer CmplxOverReal          ! The complex over real data type
   
   end module RealPrecsn
!________________________________________________________________________________________________________________________  
!____________________________________ End Module ________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine RealPrecsnInit()  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  RealPrecsnInit() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine Calculates the Hamiltonian Matrix for subspace (Nu,Nd).
! KEYWORDS: Initialization of the precision control module
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: We give the values to the quanties defined in the module RealPrecsn. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use MPISetting
      use RealPrecsn
      implicit none
!______________________________________________________________________________________________________________     
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________       
      real(rp) x
!______________________________________________________________________________________________________________     
!_________________________________ The Real quantities Output Setting _________________________________________
!______________________________________________________________________________________________________________     
      if(rp == 8) then            ! If rp=8, we store the real number as es25.16
         rp_ow = 25       
         rp_od = 16
      else if(rp == 16)  then    ! If rp=8, we store the real number as es44.34
         rp_ow = 44
         rp_od = 34
      end if
!______________________________________________________________________________________________________________     
!__________________________________ Setting to some Constants In Calculation __________________________________
!______________________________________________________________________________________________________________     
      rp_one = 1.0_rp                          ! number 1
      rp_pi = acos( -rp_one )                  ! quantity Pi
      rp_gold = (sqrt(5.0_rp) - 1.0_rp) / 2    ! golden section ratio
      rp_huge = huge(x)                        ! The largest number for this kind of rp variable
      rp_tiny = tiny(x)                        ! The smallest number for this kind of rp variable
      rp_Eps  = 1.0E-12_rp                     ! Small number for some calculations
!______________________________________________________________________________________________________________     
!__________________________________ Some constants ____________________________________________________________
!______________________________________________________________________________________________________________      
      rp_Rone  = 1.0_rp
      rp_Rzero = 0.0_rp
      rp_Z_One = cmplx(1.0_rp, 0.0_rp, rp)
      rp_Zzero = cmplx(0.0_rp, 0.0_rp, rp)
!______________________________________________________________________________________________________________     
!__________________________________ Calculate the Machine Real Precision ______________________________________
!______________________________________________________________________________________________________________     
      call MachineRealPrecsn()             ! Subroutine to calculate the machine accurancy
!______________________________________________________________________________________________________________     
!__________________________________ MPI related integers for parallel computing _______________________________
!______________________________________________________________________________________________________________
#ifdef MPIPROCESS
      if(rp == 4) then
         rp_MPI_REAL    = MPI_REAL4
         rp_MPI_COMPLEX = MPI_COMPLEX8
         CmplxOverREal  = 2
      else if(rp == 8) then
         rp_MPI_REAL    = MPI_REAL8
         rp_MPI_COMPLEX = MPI_COMPLEX16
         CmplxOverREal  = 2
      else if(rp == 16) then
         rp_MPI_REAL    = MPI_REAL16
         rp_MPI_COMPLEX = MPI_COMPLEX32
         CmplxOverREal  = 2
      end if
#endif

   end subroutine RealPrecsnInit
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
   subroutine MachineRealPrecsn()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  MachineRealPrecsn()
! TYPE:     subroutine
! PURPOSE:  This Subroutine Calculates the machine precision which will be used in the following.
! KEYWORDS: Machine Real Precision
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: We give the values to the quanties defined in the module RealPrecsn. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
      use RealPrecsn
      implicit none
      
      integer i
      real(rp) x, y, d
     
      d = 2._rp
      x = rp_one
      y = x + rp_one

      do i = 1, 100
         do while (y > x .and. rp_one + x > rp_one)
            y = x
            x = x / d
         enddo
         x = y
         y = x + rp_one
         d = sqrt(d)
      enddo
   
      rp_prec = x

   end subroutine MachineRealPrecsn
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
