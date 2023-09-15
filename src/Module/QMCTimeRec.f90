!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A module and a few subroutines used for count the time consumed during the CPMC simulations for every 
!             BIN calculations. The consumed time mainly includes two parts: 
!               (0) All the main parts in single BIN simulation;
!               (1) Time consumed by all the matrix operations.
! COMMENT: This file is used to check Calculation part for every part in CPMC code.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!     QMCTimeRec     --> Module to define the time recording quantities used in CPMC simulations;
!     QMCTimeRecInit --> Subroutine to perform the initialization for QMCTimeRec module.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Module ______________________________________________________________________
!________________________________________________________________________________________________________________________
      module QMCTimeRec
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
            implicit none
!______________________________________________________________________________________________________________     
!_________________________ Integer calculation times in QMC --> Count times ___________________________________
!______________________________________________________________________________________________________________
      !!!!!!!!!! Most time-consuming portions in QMC simulations, summation=total
      integer TimsH0Prp       ! Number of H0 propagation
      integer TimsbUPrp       ! Number of HU propagation
            integer TimsUptbU       ! Number of calling HubbU interaction updating
      integer TimsNmStb       ! Number of Numerical stablization
      integer TimsPopCt       ! Number of population control
      integer TimsConst       ! Number of calculating the growth estimator
      integer TimsStaMe       ! Number of static measurements in a single BIN simulation
      integer TimsDatPr       ! Number of data post process in a single BIN simulation
      !!!!!!!!!! Operations in BetaT --> 0 sweep simulation
      integer TimsM2One       ! Number of times of Beta --> 0 sweep
      integer TimsB0Pgt       ! Number of dynamic propagations
      integer TimsB0Stb       ! Number of Numerical stabilizations in Beta --> 0 sweep
      integer TimsB0Mea       ! Number of static measurements in Beta --> 0 sweep
      integer TimsDyMea       ! Number of dynamic measurements
      !!!!!!!!!! The mathematical opearations for matrices using MKL
      integer TimsMtFFT       ! Number of FFT calculations for matrix product
            integer TimsUDVOt       ! Number of QR decompositions in the CPMC simulations in a single BIN simulation
      integer TimsMtprd       ! Number of matrix product beyond the kinetic energy part matrix multiplication 
            integer TimsMtInv       ! Number of matrix inverse in the CPMC simulations in a single B\IN simulation
      integer TimsMtDet       ! Number of matrix determinants
      integer TimsEqSet       ! Number of applying solving linear equation set for AMat*BMat^{-1}
      !!!!!!!!!! Computations for single-particle Green's functions
      integer TimsGFSta       ! Number of calculating Static Green's Function Matrices
      integer TimsGFDyn       ! Number of calculating Dynamic Green's Function Matrices
!______________________________________________________________________________________________________________     
!_________________________ Real calculation time in QMC --> Count time ________________________________________
!______________________________________________________________________________________________________________
      !!!!!!!!!! Most time-consuming portions in QMC simulations, summation=total
            real(8) TimeH0Prp       ! Time consumed for H0 propagation
      real(8) TimebUPrp       ! Time consumed for HU propagation
            real(8) TimeUptbU       ! Time consumed for calling HubbU interaction updating
      real(8) TimeNmStb       ! Time consumed for Numerical stablization
      real(8) TimePopCt       ! Time consumed for population control
      real(8) TimeConst       ! Time consumed for calculating the growth estimator
      real(8) TimeStaMe       ! Time consumed for static measurements in a single BIN simulation
      real(8) TimeDatPr       ! Time consumed for data post process in a single BIN simulation
      !!!!!!!!!! Operations in BetaT --> 0 sweep simulation
      real(8) TimeM2One       ! Time consumed for Beta --> 0 sweep
      real(8) TimeB0Pgt       ! Time consumed for propagations as Beta --> 0
      real(8) TimeB0Stb       ! Time consumed for numerical stabilizations in Beta --> 0
      real(8) TimeB0Mea       ! Time consumed for static  measurements in Beta --> 0
      real(8) TimeDyMea       ! Time consumed for dynamic measurements in Beta --> 0
      !!!!!!!!!! The mathematical opearations for matrices using MKL
      real(8) TimeMtFFT       ! Time consumed for FFT for kinetic term in Hamiltonian
            real(8) TimeUDVOt       ! Time consumed for QR decompositions in the CPMC simulations in a single BIN simulation
      real(8) TimeMtprd       ! Time consumed for matrix product beyond the kinetic energy part matrix multiplication 
            real(8) TimeMtInv       ! Time consumed for matrix inverse in the CPMC simulations in a single B\IN simulation
      real(8) TimeMtDet       ! Time consumed for matrix determinants
      real(8) TimeEqSet       ! Time consumed for applying solving linear equation set for AMat*BMat^{-1}
      real(8) TimeMtSum       ! Time consumed for matrix operation by LaPack
      !!!!!!!!!! Computations for single-particle Green's functions
      real(8) TimeGFSta       ! Time consumed for calculating static  Green's function matrices
      real(8) TimeGFDyn       ! Time consumed for calculating Dynamic Green's function matrices
      !!!!!!!!!! Count the total time consumed for a single BIN
            real(8) TimeTotal       ! Time consumed for all the above calculations
            real(8) TimeSgBIN       ! Time consumed for a single BIN simulation
            
      end module QMCTimeRec
!________________________________________________________________________________________________________________________  
!____________________________________ End Module ________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
      subroutine QMCTimeRecInit()  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  QMCTimeRecInit()  
! TYPE:     subroutine
! PURPOSE:  This Subroutine initialize the quantities defined in QMCTimeRec module.
! KEYWORDS: Initialization of time recording quantities.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION: We give the values to the quanties defined in the module QMCTimeRec. 
!
!     Input:  (none)   Output: (none)
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________
            use QMCTimeRec 
            implicit none
!______________________________________________________________________________________________________________     
!_________________________ Integer calculation times in QMC --> Count times ___________________________________
!______________________________________________________________________________________________________________
      !!!!!!!!!! Most time-consuming portions in QMC simulations, summation=total
            TimsH0Prp = 0       ! Number of H0 propagation
      TimsbUPrp = 0       ! Number of HU propagation
            TimsUptbU = 0       ! Number of calling HubbU interaction updating
      TimsNmStb = 0       ! Number of Numerical stablization
      TimsPopCt = 0       ! Number of population control
      TimsConst = 0       ! Number of calculating the growth estimator
      TimsStaMe = 0       ! Number of static measurements in a single BIN simulation
      TimsDatPr = 0       ! Number of data post process in a single BIN simulation
      !!!!!!!!!! Operations in BetaT --> 0 sweep simulation
      TimsM2One = 0       ! Number of times of Beta --> 0 sweep
      TimsB0Pgt = 0       ! Number of dynamic propagations
      TimsB0Stb = 0       ! Number of numerical stabilizations in Beta --> 0
      TimsB0Mea = 0       ! Number of static  measurements in Beta --> 0
      TimsDyMea = 0       ! Number of dynamic measurements
      !!!!!!!!!! The mathematical opearations for matrices using MKL
      TimsMtFFT = 0       ! Number of FFT calculations for matrix product
            TimsUDVOt = 0       ! Number of QR decompositions in the CPMC simulations in a single BIN simulation
      TimsMtprd = 0       ! Number of matrix product beyond the kinetic energy part matrix multiplication 
            TimsMtInv = 0       ! Number of matrix inverse in the CPMC simulations in a single B\IN simulation
      TimsMtDet = 0       ! Number of matrix determinants
      TimsEqSet = 0       ! Number of applying solving linear equation set for AMat*BMat^{-1}
      !!!!!!!!!! Computations for single-particle Green's functions
      TimsGFSta = 0       ! Number of calculating Static Green's Function Matrices
      TimsGFDyn = 0       ! Number of calculating Dynamic Green's Function Matrices
!______________________________________________________________________________________________________________     
!_________________________ Real calculation time in QMC --> Count time ________________________________________
!______________________________________________________________________________________________________________
      !!!!!!!!!! Most time-consuming portions in QMC simulations, summation=total
            TimeH0Prp = 0.0_8       ! Time consumed for H0 propagation
      TimebUPrp = 0.0_8       ! Time consumed for HU propagation
            TimeUptbU = 0.0_8       ! Time consumed for calling HubbU interaction updating
      TimeNmStb = 0.0_8       ! Time consumed for Numerical stablization
      TimePopCt = 0.0_8       ! Time consumed for population control
      TimeConst = 0.0_8       ! Time consumed for calculating the growth estimator
      TimeStaMe = 0.0_8       ! Time consumed for static measurements in a single BIN simulation
      TimeDatPr = 0.0_8       ! Time consumed for data post process in a single BIN simulation
      !!!!!!!!!! Operations in BetaT --> 0 sweep simulation
      TimeM2One = 0.0_8       ! Time consumed for Beta --> 0 sweep
      TimeB0Pgt = 0.0_8       ! Time consumed for propagations as Beta --> 0
      TimeB0Stb = 0.0_8       ! Time consumed for numerical stabilizations in Beta --> 0
      TimeB0Mea = 0.0_8       ! Time consumed for static  measurements in Beta --> 0
      TimeDyMea = 0.0_8       ! Time consumed for dynamic measurements in Beta --> 0
      !!!!!!!!!! The mathematical opearations for matrices using MKL
      TimeMtFFT = 0.0_8       ! Time consumed for FFT for kinetic term in Hamiltonian
            TimeUDVOt = 0.0_8       ! Time consumed for QR decompositions in the CPMC simulations in a single BIN simulation
      TimeMtprd = 0.0_8       ! Time consumed for matrix product beyond the kinetic energy part matrix multiplication 
            TimeMtInv = 0.0_8       ! Time consumed for matrix inverse in the CPMC simulations in a single B\IN simulation
      TimeMtDet = 0.0_8       ! Time consumed for matrix determinants
      TimeEqSet = 0.0_8       ! Time consumed for applying solving linear equation set for AMat*BMat^{-1}
      TimeMtSum = 0.0_8       ! Time consumed for matrix operation by LaPack
      !!!!!!!!!! Computations for single-particle Green's functions
      TimeGFSta = 0.0_8       ! Time consumed for calculating static  Green's function matrices
      TimeGFDyn = 0.0_8       ! Time consumed for calculating Dynamic Green's function matrices
      !!!!!!!!!! Count the total time consumed for a single BIN
            TimeTotal = 0.0_8       ! Time consumed for all the above calculations (not including matrix multiplication in UDV)
            TimeSgBIN = 0.0_8       ! Time consumed for a single BIN simulation
            
   end subroutine QMCTimeRecInit
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$