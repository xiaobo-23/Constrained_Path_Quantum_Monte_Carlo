!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine __________________________________________________________________
!________________________________________________________________________________________________________________________
	subroutine SwpDatProc() 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  SwpDatProc()
! TYPE:     subroutine
! PURPOSE:  This subroutine is used to process the static measurements of physical observables for single sweep. 
! KEYWORDS: Data process for a single sweep.
! AUTHOR:   Yuan-Yao He
! TIME:     2020-02-27
! DESCRIPTION:
!
!     The static measurements of single sweep includes the measure at Tau=0 point as well as [BetaT, 0] measure.
!          So we process the static results for them individually here. 
!
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		use RealPrecsn
      use CoreParamt
      use MPISetting
      use Observable
      implicit none
!______________________________________________________________________________________________________________
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
      integer MeaInd
!______________________________________________________________________________________________________________	  
!__________________________ Process the static measurements for single sweep __________________________________
!______________________________________________________________________________________________________________
!**************************************************************************************************	  
!___________________ 0. Postprocess of the computation results ____________________________________
!************************************************************************************************** 
      do MeaInd = 0, merge(1, 0, IfM2OneMea), +1
!________________________________________________________________________________________ 	  
!_________________ (0) Accumulate results of energies and the GrF matrices ______________
!________________________________________________________________________________________
         !!!!!!!!!! The energies, fillings and density correlations
         EngOccCrFBIN(:, MeaInd) = EngOccCrFBIN(:, MeaInd) + EngOccCrFSwp(:, MeaInd)*WeightAvgSwp
         !!!!!!!!!! The r-space and k-space single-particle GrFs
         RlSpGrnFtBIN(:, :, :, MeaInd) = RlSpGrnFtBIN(:, :, :, MeaInd) + RlSpGrnFtSwp(:, :, :, MeaInd)*WeightAvgSwp
         if(IfFftEnPar) then
            KSpGreenFBIN(:, :, :, MeaInd) = KSpGreenFBIN(:, :, :, MeaInd) + KSpGreenFSwp(:, :, :, MeaInd)*WeightAvgSwp
         end if
!________________________________________________________________________________________ 	  
!_________________ (1) Accumulate results of n(k) and pairing matrices __________________
!________________________________________________________________________________________
         if(IfFftEnPar) then
            NkSgleBIN(:, :, MeaInd) = NkSgleBIN(:, :, MeaInd) + NkSgleSwp(:, :, MeaInd)*WeightAvgSwp
            PairMtBIN(:, :, MeaInd) = PairMtBIN(:, :, MeaInd) + PairMtSwp(:, :, MeaInd)*WeightAvgSwp
         end if
!________________________________________________________________________________________ 	  
!_________________ (2) Accumulate results of r-space correlations (PBC or OBC) __________
!________________________________________________________________________________________  
         if( abs(PinSz) < rp_Eps ) then
            RealSpCrFBIN(:, :, MeaInd) = RealSpCrFBIN(:, :, MeaInd) + RealSpCrFSwp(:, :, MeaInd)*WeightAvgSwp 
         end if
!________________________________________________________________________________________ 	  
!_________________ (3) Accumulate results of dynamic correlations _______________________
!________________________________________________________________________________________  
         if( (MeaInd == 1) .and. (IfTAU) ) then
            RealSpCrFTauBIN(:, :, :) = RealSpCrFTauBIN(:, :, :) + RealSpCrFTauSwp(:, :, :)*WeightAvgSwp
         end if
      enddo
      
   end subroutine SwpDatProc
!________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ____________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$