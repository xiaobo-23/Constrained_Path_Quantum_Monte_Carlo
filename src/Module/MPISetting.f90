!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A module used to define some integer quantities for the MPI parallel calculations for the CPMC program.
! COMMENT: MPI parallel calculations.  
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!    MPISetting     --> Subroutine to define the related quantities for MPI parallel computations;
!    MPISettingInit --> Subroutine to perform the initializations of MPI parallel computations;
!    MPISettingFinl --> Subroutine to perform finalization of MPI parallel computations.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Module ______________________________________________________________________
!________________________________________________________________________________________________________________________
      module MPISetting
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      implicit none

#ifdef MPIPROCESS
      include "mpif.h"
#endif
      
!**************************************************************************************************     
!___________________ 0. Integer numbers for the cores during the simulations ______________________
!**************************************************************************************************
      integer acomm                      ! The integer index for the communication region
      integer anprc                      ! Number of cores used in simulations
      integer amyid                      ! Integer index for the current used core during the simulation
      integer amstr                      ! Integer index for the core used to process and output data
      integer ierr                       ! The Error setting for the MPI subroutines
#ifdef MPIPROCESS
      integer stts(MPI_STATUS_SIZE)      ! The MPI status
#endif
!**************************************************************************************************     
!___________________ 1. Integer parameters used in time-displayed measurements ____________________
!**************************************************************************************************
      integer OneMegaNum                 ! How manny double precision numbers 1 MB contains

end module MPISetting
!________________________________________________________________________________________________________________________  
!____________________________________ End Module ________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine MPISettingInit()  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  MPISettingInit() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the initialization for the MPISetting module.
! KEYWORDS: Initialization of the MPISetting module. 
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: We give the values to the quanties defined in the module MPISetting. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!______________________________________ Initialization process ________________________________________________
!______________________________________________________________________________________________________________   
#ifdef MPIPROCESS
            call MPI_INIT(ierr)                       ! Initialization of MPI enviaroment
            acomm = MPI_COMM_WORLD                    ! Get the communication region
            call MPI_COMM_SIZE(acomm, anprc, ierr)    ! Get the number of cores as the size of communication region
            call MPI_COMM_RANK(acomm, amyid, ierr)    ! Get the current core integer number
            amstr = 0                                 ! The core used for read in parameter, data process and output
#else
      acomm = 0
      anprc = 1
      amstr = 0
      amyid = 0
      ierr  = 0
#endif      
            OneMegaNum = ISHFT(1, 17)                 ! How manny double precision numbers 1 MB contains
      
      end subroutine MPISettingInit
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine MPISettingFinl()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  MPISettingFinl() 
! TYPE:     subroutine
! PURPOSE:  This Subroutine performs the finalization for the MPISetting module.
! KEYWORDS: Finalization of the MPISetting module. 
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: Finalization of the MPISetting module. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
      use MPISetting
      implicit none
!______________________________________________________________________________________________________________     
!______________________________________ Finalization process __________________________________________________
!______________________________________________________________________________________________________________
#ifdef MPIPROCESS
      call MPI_FINALIZE(ierr)                       ! Finalization of MPI enviaroment
#endif

end subroutine MPISettingFinl
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$