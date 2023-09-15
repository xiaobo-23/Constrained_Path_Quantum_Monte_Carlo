!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 09/25/2022
! ADD SINUSOIDAL SPIN PINNING FIELDS WITH PERIODIC BOUNDARY CONDITION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: A module used for defining the related quantities for CPMC calculations. All the quantities include
!              the following different parts:
!               0. The lattice geometry for the simulation;
!               1. Model and simulation parameters set;
!               2. Exp(-\Delta\tau\hat{V}) and Exp(-\Delta\tau\hat{K}) related quantities;
!               3. Core matrices used in CPMC simulations;
!               4. Some other quantities used in CPMC simulations.
! COMMENT: For extended Hubbard model with sign problem on square lattice.
! AUTHOR:  Yuan-Yao He
! DATE:    2020-02-27
! PURPOSE: Different subroutines are introduced as following:
!             
!   CoreParamt --> Define the core variables for the calculation; 
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Module ______________________________________________________________________
!________________________________________________________________________________________________________________________
      module CoreParamt
!______________________________________________________________________________________________________________     
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________      
            use RealPrecsn
            implicit none
!**************************************************************************************************  
!___________________ 0. Define the finite lattice size ____________________________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) Lattice size of the system _______________________________________
!________________________________________________________________________________________       
            integer NumL1       ! Number of unit cells in a1-direction               
            integer NumL2       ! Number of unit cells in a2-direction             
            integer NumNC       ! Total number of unit cells
            integer NumNS       ! Total number of lattice sites for the finite lattice system = NumNC
!________________________________________________________________________________________  
!_________________ (1) The basic lattice vectors in real and reciprocal _________________
!________________________________________________________________________________________ 
            real(rp) A1Vec(2)   ! The a1 lattice vector in real space                 
            real(rp) A2Vec(2)   ! The a2 lattice vector in real space 
            real(rp) B1Vec(2)   ! The b1 basic vector in reciprocal space            
            real(rp) B2Vec(2)   ! The b2 basic vector in reciprocal space
!________________________________________________________________________________________  
!_________________ (2) Nearest-Neighbor unit cells of all unit cells ____________________
!________________________________________________________________________________________
            integer, allocatable ::    StList(:, :)    ! I = Ix * a1 + Iy * a2; StList(I, 1) = Ix, StList(I, 2) = Iy
            integer, allocatable :: InvStList(:, :)    ! InvStList(Ix, Iy) = I, Ix \in [1, L1],  Iy \in [1, L2]
            integer, allocatable ::    KpList(:, :)    ! wavevector points
            integer, allocatable :: InvKpList(:, :)    ! 
!________________________________________________________________________________________  
!_________________ (3) Distance indexes for two unit cells in the finite lattice ________
!________________________________________________________________________________________
            integer, allocatable :: IminusJ(:, :)
!________________________________________________________________________________________  
!_________________ (4) NN and 2NN lattice sites information _____________________________
!________________________________________________________________________________________
            integer, allocatable :: FNNBond(:, :)
            integer, allocatable :: SNNBond(:, :)
            integer, allocatable :: TNNBond(:, :)
!________________________________________________________________________________________  
!_________________ (5) Open or periodic boundary conditions in a1 and a2 directions _____
!________________________________________________________________________________________
            logical IfPr1                               ! Whether to apply the periodic boundary condition in a1 direction
            logical IfPr2                               ! Whether to apply the periodic boundary condition in a2 direction
            integer, allocatable :: FNNStBnd(:, :)     ! Boundary condition for First  NN sites
            integer, allocatable :: SNNStBnd(:, :)     ! Boundary condition for Second NN sites
            integer, allocatable :: TNNStBnd(:, :)     ! Boundary condition for Third  NN sites
!**************************************************************************************************     
!___________________ 1. Model Parameters and simulation parameters ________________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) The model parameters for the t-t2-U Model ________________________
!________________________________________________________________________________________
!____________________________________________________________________________         
!________________ [0] Model Hamiltonian parameters __________________________
!____________________________________________________________________________
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            logical ifSinusoidalPinning ! Whether to add sinusoidal spin pinnign fields 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            real(rp) Hopt1Up         ! The Nearest-Neighbor hopping parameter for spin-up electrons
            real(rp) Hopt1Dn         ! The Nearest-Neighbor hopping parameter for spin-down electrons
            real(rp) Hopt2           ! The Next-Nearest-Neighbor hopping parameter
            real(rp) Hopt3           ! The third-Nearest-Neighbor hopping parameter
            real(rp) ZmFdz           ! The Zeeman field in z-direction
            real(rp) PinSz           ! The pinning field for Sz order (only along the edges)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            real(rp) SinusoidalPinSz ! The amplitude of sinusoidal pinning fields
            real(rp) LambdaSz        ! The wavelength of sinusoidal pinning fields
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            real(rp) ChemP           ! The chemical potential for dopping
            real(rp) HubbU           ! The on-site Hubbard Repulsion interaction strength
            integer PinSzType        ! Type of pinning field of Sz. ==0 --> Only along one edge; == 1 --> both edges; == 2 --> All sites.
!________________________________________________________________________________________  
!_________________ (1) Parameters used in CPMC simulations ______________________________
!________________________________________________________________________________________
!____________________________________________________________________________         
!________________ [0] Inverse temperature, dltau, LTrot _____________________
!____________________________________________________________________________
            integer LTrot            ! Number of time slices. LTrot = BetaT / Dltau
            logical IfSetDltau       ! ==F, don't fix dtau, use BetaT and LTrot to get dtau; ==T, fix dtau. 
            real(rp) BetaT           ! The inverse temperature in CPMC
            real(rp) TempT           ! The temperature for the simulation
            real(rp) Dltau           ! The imaginary time slice, Dltau = BetaT / LTrot      
            real(rp) FixedDltau      ! The value of the fixed dtau
!____________________________________________________________________________         
!________________ [1] Number of random walkers ______________________________
!____________________________________________________________________________
            integer NWalk            ! Number of random walkers on single process
            integer NmWalkAllP       ! Number of random walkers on all processes
            integer NWkBt            ! Number of independent random walkers at \tau=\beta on single process
            integer NmWkBtAllP       ! Number of independent random walkers at \tau=\beta on all processes
!____________________________________________________________________________         
!________________ [2] Warmup, BIN, Sweep and measure ________________________
!____________________________________________________________________________
            integer nWarm            ! Number of sweeps of the warm up for CPMC simulation
            integer NmBin            ! Number of measurements independently
            integer NSwep            ! Number of sweep in a single Bin
!____________________________________________________________________________         
!________________ [3] Numerical stablization process related ________________
!____________________________________________________________________________    
            integer NvStb            ! Interval number of time slices for the re-othogonalization
            integer NvStbOut         ! Number of sweep intervals for comparing Green's function error
            integer StabOutput       ! Interval number of sweeps for the re-othogonalization --> warmup and measure sweep
      
            integer NvPop            ! Interval number of time slices for the population control
            integer NvPopOut         ! Number of sweep intervals for Output the population information
            integer PoptOutput       ! Interval number of sweeps for population control      --> warmup and measure sweep
      
            real(rp) PopWghtMax      ! Allowed maximum weight for population control
            real(rp) PopWghtMin      ! Allowed minimum weight for population control  
!________________________________________________________________________________________  
!_________________ (2) Whether to fix the total electron density ________________________
!________________________________________________________________________________________
            integer Fix_Iterate      ! Number of iterations to get the desired total electron density
      
            logical IfFixnT          ! ==0, don't fix nOccpy and use the input chemical potential; ==1, fix n_occupy
            integer NSwepFix         ! Number of sweeps for the measurement of n_occupy for fixed nOccpy
      
            integer FixNenT          ! ==0, fix Ne; == 1, fix nT
            real(rp) NumNe           ! Desired number of electrons, up + down
            real(rp) Fix_nT          ! The desired total electron density
            real(rp) ChPSmallnT      ! Input chemical potential corresponding to small electron density
            real(rp) ChPLargenT      ! Input chemical potential corresponding to large electron density
            
            real(rp) ChemP_Ref(2)    ! The small step size we change ChemP
      
            real(rp) Ne_Now          ! The present number of total electrons
            real(rp) Ne_Now_Err      ! The errorbar of present number of total electrons
            real(rp) nT_Now          ! The present total electron density
            real(rp) nT_Now_Err      ! The errorbar of present total electron density
            real(rp) WgtNow          ! Sign average in this fixed n_Occ measurements
            real(rp) WgtNow_Err      ! The errorbar of Sign average in this fixed n_Occ measurements
!**************************************************************************************************  
!___________________ 2. Exp(-\Delta\tau\hat{V}) and Exp(-\Delta\tau\hat{K}) _______________________
!**************************************************************************************************  
!________________________________________________________________________________________  
!_________________ (0) Number of spin channels we need to simulate ______________________
!_____________________ Whether to use symmetric Trotter Decomposition ___________________
!________________________________________________________________________________________
            integer NmSpn
            logical SymTrotDcp                 ! == .true. --> asymmetric; == .false. --> symmetric.
!________________________________________________________________________________________  
!_________________ (1) For the kinetic propagator of H_0 part in many-body H ____________
!________________________________________________________________________________________
!____________________________________________________________________________   
!________________ [0] Eigenvalues and eigenvectors for H_0 Hamiltonian ______
!____________________________________________________________________________
            real(rp), allocatable :: H0_EgValue(:   , :)        ! The eigenvalues of K(B_T)
            real(rp), allocatable :: H0_EgVctor(:, :, :)        ! The eigenvector matrix of K(B_T)
!____________________________________________________________________________   
!________________ [1] Use FFT method to calculate ___________________________
!____________________________________________________________________________ 
            logical FFTEXPDTH0                                  ! Whether to use FFT method
            integer EkDispType                                  ! ==0 for usual Hubbard dispersion, -2*(cos(kx)+cos(ky))
                                                          ! ==1 for Hubbard   dispersion, 2*(2-cos(kx)-cos(ky))
                                                          ! ==2 for Quadratic dispersion, kx^2 + ky^2
            real(rp), allocatable :: ExpdtEofH0(:, :, :)        ! Exp(-/+ Dltau*E_0(k)) and Exp(-/+ Dltau*E_0(k)/2.0)        
!____________________________________________________________________________   
!________________ [2] Use matrix product of full matrices ___________________
!____________________________________________________________________________
            real(rp), allocatable :: ExpdtH0Mat(:, :, :, :)     ! Exp(-/+Dltau*H_0) and Exp(-/+Dltau*H_0/2)
!________________________________________________________________________________________  
!_________________ (2) Exp(-\Delta\tau\hat{V}) related quantities _______________________
!________________________________________________________________________________________
!____________________________________________________________________________   
!________________ [0] The type of the HS transformation _____________________
!____________________________________________________________________________    
            integer HS_Type                                  ! HS_Type==0 for CDW; HS_Type==1 for Sz
!____________________________________________________________________________         
!________________ [1] Update method and the size of block ___________________
!____________________________________________________________________________       
            integer UpdtMethod                ! ==0, local update; ==1, delayed update; ==2, force-bias update.
            integer NblkUDelay                     ! Size of the block used in delayed updates and force-bias updates.
!____________________________________________________________________________ 
!________________ [2] Auxiliary fields for HubbU interaction ________________
!____________________________________________________________________________ 
            integer, allocatable :: IsingbU(:, :, :)         ! Auxiliary fields for HubbU interaction
!____________________________________________________________________________   
!________________ [3] Two-component auxiliary fields and ____________________
!____________________ Corresonding \Delta_{ii} and phase ____________________
!____________________________________________________________________________  
            real(rp) LmdbU                                   ! The constant in HS transformation for U interaction term
            real(rp), allocatable :: ExpobU   (:, :, :)      ! exp(+\lambda*x_{\ell,i})
            real(rp), allocatable :: ExpobUInv(:, :, :)      ! exp(-\lambda*x_{\ell,i})
            real(rp), allocatable :: ExpobU_H0T   (:, :, :)  ! Constant used for the (H_U + H_0 - H_T) term
            real(rp), allocatable :: ExpobUInv_H0T(:, :, :)  ! Constant used for the (H_U + H_0 - H_T) term
            real(rp), allocatable :: DeltbU_H0T   (:, :, :)  ! Constant used for the (H_U + H_0 - H_T) term
!________________________________________________________________________________________  
!_________________ (3) For the trial Hamiltonian H_T used in B_T propagator _____________
!________________________________________________________________________________________
!____________________________________________________________________________   
!________________ [0] Eigenvalues and eigenvectors for H_T Hamiltonian ______
!____________________________________________________________________________       
            real(rp), allocatable :: HT_EigValu(:   , :)      ! The eigenvalues of K(B_T)
            real(rp), allocatable :: HT_EigVect(:, :, :)      ! The eigenvector matrix of K(B_T)
!____________________________________________________________________________   
!________________ [1] Deal with the exp(+/-dt*H_T) matrices _________________
!____________________________________________________________________________ 
            !!!!!!!!!! Use FFT method to calculate
            logical FFTEXPDTHT                               ! Whether to use FFT method
            real(rp), allocatable :: ExpdtOfHTe(:, :, :)     ! Exp(-/+ Dltau*E_T(k)) and Exp(-/+ Dltau*E_T(k)/2.0) 
            !!!!!!!!!! Use matrix product of full matrices
            real(rp), allocatable :: ExpdtTryHT(:, :, :, :)  ! Exp(-/+Dltau*H_T) and Exp(-/+Dltau*H_T/2)
!____________________________________________________________________________   
!________________ [2] Quantities used to determine ChemP_BT parameter _______
!____________________________________________________________________________
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@ MuBTType == 0 --> Simply choose ChemP_BT = ChemP      @@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@ MuBTType == 1 --> Read ChemP_BT from input parameter  @@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@ MuBTType == 2 --> Obtain ChemP_BT from input Dnsty_BX @@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@ MuBTType == 3 --> Obtain ChemP_BT from Fix_nT         @@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@ MuBTType == 4 --> Obtain ChemP_BT by n_BT = n_QMC     @@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
            integer MuBTType      ! Method to set ChemP_BT, ==0, just read; ==1 calculate from Dnsty_BT or just use ChemP
            logical IfMuTqmcNt    ! Whether to determine ChemP_BT by n_T = n from QMC simulations
            integer NItrGetMuT    ! Number of iterations to get the Mu_T for n_T = n
            real(rp) ChemP_BT     ! The chemical potential used in H_T Hamiltonian
!____________________________________________________________________________   
!________________ [3] For the UHF kind of trial Hamiltonian H_T _____________
!____________________________________________________________________________      
            real(rp) HubbU_UHF
            real(rp), allocatable :: UHF_Const(:   )
            real(rp), allocatable :: MagMoment(:, :)
!************************************************************************************************** 
!___________________ 3. The Core matrices used in CPMC simulations ________________________________
!**************************************************************************************************
!________________________________________________________________________________________  
!_________________ (0) UDV matrices at left and right sides for propagation _____________
!________________________________________________________________________________________
            real(rp), allocatable :: ULeft(:, :, :)         ! UDV decomposition for Left part
            real(rp), allocatable :: DLeftVec(:, :)         ! UDV decomposition for Left part
            real(rp), allocatable :: VLeft(:, :, :)         ! UDV decomposition for Left part
      
            real(rp), allocatable :: URght(:, :, :, :)      ! UDV decomposition for Right part
            real(rp), allocatable :: DRghtVec(:, :, :)      ! UDV decomposition for Right part
            real(rp), allocatable :: VRght(:, :, :, :)      ! UDV decomposition for Right part

            real(rp), allocatable :: LogScaleRght(:, :)     ! Log of the scale for the right side
            real(rp), allocatable :: LogScaleLeft(:   )     ! Log of the scale for the left  side
            real(rp), allocatable :: LogScaleOfHT(:, :)     ! Log of the scale for the H_T Hamiltonian
!________________________________________________________________________________________  
!_________________ (1) The Green's function matrix for propagation ______________________
!________________________________________________________________________________________
            real(rp) GrFTrcDhld                             ! Numbers <= 10^(-GrFTrcDhld), take it as zero
            real(rp) ChemP_Left, ChemP_Rght, ChemP_Sgle     ! ChemP terms for left, right sides and Single
            real(rp), allocatable :: GrnFunct(:, :, :, :)   ! The equal-time single-particle Green's function
!________________________________________________________________________________________  
!_________________ (2) Weight and their logs, and all weights vector ____________________
!________________________________________________________________________________________     
            integer,  allocatable :: IdptWkIndx(:)          ! Indexes for indepdent random walkers of present process
            real(rp), allocatable :: WghtProc(:)            ! Weights of walkers in every process  
            real(rp), allocatable :: Log_Wght(:)            ! Log of weights in every process
            real(rp), allocatable :: WghtTotl(:)            ! Weights of all walkers on all processes
!________________________________________________________________________________________  
!_________________ (3) Calculate the growth estimator ___________________________________
!________________________________________________________________________________________
            integer FrqReCmptGrowth                         ! Frequecy for recomputing growth estimator
            real(rp) Sum_WghtOld, Sum_WghtNew               ! Used to calculate GrowthCoefft
            real(rp) MeanWghtPop                            ! Used to calculate GrowthCoefft
            real(rp), allocatable :: GrowthCoefft(:)        ! The accumulated constants during propagation
!________________________________________________________________________________________  
!_________________ (4) The ancestry link and number of ancestry walkers _________________
!________________________________________________________________________________________
            integer NancestryWalker
            integer, allocatable :: AncestryLink(:)
!________________________________________________________________________________________   
!_________________ (5) Store the UDV matrices for measurement in [BetaT, 0] sweep _______
!________________________________________________________________________________________
            integer NvStbM2One                          ! Step size for the numerical stablization in [BetaT, 0] sweep
            integer NvMeaM2One                          ! Step size for static measurements in [BetaT, 0] sweep
            integer DmUST                               ! Third dimension of storing matrices
            real(rp), allocatable :: UStMt(:, :, :, :)  ! Store all the UDV matrices along the sweeps
            real(rp), allocatable :: DVecStMt(:, :, :)  ! Store all the UDV matrices along the sweeps
            real(rp), allocatable :: DScalLog(   :, :)  ! Store the log scales of values in DVec matrices
            real(rp), allocatable :: VStMt(:, :, :, :)  ! Store all the UDV matrices along the sweeps
!________________________________________________________________________________________  
!_________________ (6) Temporary Green's Function matrices for Static and Dynamic _______
!________________________________________________________________________________________  
            real(rp)   , allocatable :: GrnFRTmp00(:, :, :)    ! Temporary Single-particle Green's function, real
            real(rp)   , allocatable :: GrnFRTmp11(:, :, :)    ! Temporary Single-particle Green's function, real
            real(rp)   , allocatable :: GrnFRTmp22(:, :, :)    ! Temporary Single-particle Green's function, real
            complex(rp), allocatable :: GrFKCTmp00(:, :, :)    ! Temporary Single-particle Green's function, complex
            complex(rp), allocatable :: GrFKCTmp11(:, :, :)    ! Temporary Single-particle Green's function, complex
!________________________________________________________________________________________  
!_________________ (7) Temporary Green's Function matrices for Dynamic __________________
!________________________________________________________________________________________ 
            real(rp), allocatable :: GrF00(:, :, :)     ! The equal-time Green Function <c_i(\tau=0) c_j^+>
            real(rp), allocatable :: GrF0T(:, :, :)     ! The time-displaced Green Function <c_i(0) c_j^+(tau)>
            real(rp), allocatable :: GrFT0(:, :, :)     ! The time-displaced Green Function <c_i(tau) c_j^+(0)>
            real(rp), allocatable :: GrFTT(:, :, :)     ! The equal-time Green Function <c_i(tau) c_j^+(tau)>
!________________________________________________________________________________________  
!_________________ (8) Identity matrix, used for initializations ________________________
!________________________________________________________________________________________      
            real(rp), allocatable :: IdMtR(:, :, :)
!**************************************************************************************************     
!___________________ 4. Some other quantities used in CPMC simulations ____________________________
!**************************************************************************************************    
!________________________________________________________________________________________         
!_________________ (0) For the final data process --> Average and errorbar ______________
!________________________________________________________________________________________
            logical IfCutSmLg                     ! Whether we cut the smallest and largest value
            integer NmBinCut                      ! Number of Bins that is cut (not count in data statistics)
            integer BinStart                      ! The starting BIN number for data statistics
            integer Bin_Data                      ! Total number of BIN for data statistics     
!________________________________________________________________________________________         
!_________________ (1) Compare error of static and dynamic Green's Function _____________
!________________________________________________________________________________________
            logical IfCheckUDV                    ! Whether to check the UDV decomposition
            logical IfNmStbErr                    ! Whether to compare and output the static Green's function difference
            integer NComp_StaDyn(3)               ! Number of times we compare the Green's function
            real(rp) Xmaxm_StaDyn(3)              ! The maximum difference of the comparison
            real(rp) Xmean_StaDyn(3)              ! The mean difference of the comparison
!________________________________________________________________________________________         
!_________________ (2) Whether to output Max and Min of weights _________________________
!________________________________________________________________________________________
            logical IfWghtMxMn                    ! Whether to output Max and Min of weights                        
            integer WghtOutNSW                    ! The sweep index for output in every BIN
            real(rp) UpdtRtTotl                   ! The normalization factor for every time slice
            real(rp) UpdtRtTotMin                 ! The Min of the normalization factor of all walkers
            real(rp) UpdtRtTotMax                 ! The Max of the normalization factor of all walkers
!________________________________________________________________________________________  
!_________________ (3) Whether to save the auxiliary-field configurations _______________
!________________________________________________________________________________________
            logical IfSaveFlds                    ! Whether to save the auxiliary fields
            integer LengthBite                    ! Combine LengthBite-digits of 0(1) as a binary number
            integer NWlkFldOut                    ! Number of walker configuration to output
            integer LTrotNumNS                    ! LTrotNumNS = LTrot * NumNS
            integer NmBinaryFd                    ! Number of binary numbers for all the fields in IsingbU
            integer ReclUnitNm                    ! ReclUnitNm == LTrot*NumNS for ifort and == LTrot*NumNS*4 for gfortran
            logical SaveFldMea                    ! Whether to save fields for all the measurements
            logical ReadFldMea                    ! Whether to read fields for all the measurements
            integer(8) FldIndRdWt                 ! Number of rows (time slices) of auxiliary field when reading
            integer(8) WghtIdRead                 ! Number of rows (time slices) of measurements    when reading
            character(500) MeaFldPath             ! The path of measuring configurations
            integer, allocatable :: RdWrtField(:) ! Temperory for writing or reading fields
            integer, allocatable :: RdWrtIntgr(:) ! Temperory for writing or reading fields
      
   end module CoreParamt
!________________________________________________________________________________________________________________________  
!____________________________________ End Module ________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$