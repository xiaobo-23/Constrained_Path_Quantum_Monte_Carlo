######################################################################################################################
########################################## Information of this file ##################################################
######################################################################################################################
##
## Introduction
## ============
##
## This Makefile is used to compile the .f90 files in the "Special" folder.
##
## Usage
## =====
##
## 'make' or 'make clean'
##
## Dependence
## ==========
##
## This file depends on the "make.sys" file in the upper-level folder.
##
## Author
## ======
##
## This building system is designed, created, implemented, and maintained by
##
## Author: Yuan-Yao He
## email : heyuanyao@163.com
##
## History
## =======
##
## 04/27/2017 by Yuan-Yao He
##
######################################################################################################################
######################################################################################################################
######################################################################################################################



######################################################################################################################
############################################ Commands of this file ###################################################
######################################################################################################################
.SUFFIXES: .f90

include ../make.sys

QMCCore = ExpH0MultM.o ExpbUMultM.o MatOperatA.o MatOperatB.o NmStablize.o PopControl.o GrFCompute.o UpdtIsngbU.o \
          SubroutnTm.o SaveCalTim.o GrwthCoeff.o SaveCfgFld.o PhyMeaStat.o ObStaEnrgy.o PairMatrix.o ObStaCrFct.o \
			 PostStatic.o HTStatcObs.o SwpDatProc.o
QMCDynm = ObStaDynmc.o ProcCrFIwn.o PostDynamc.o
QMCInit = A0RandPrep.o A1LattStru.o A2TrialHmt.o A3ExpHmltK.o A4ExpHmltV.o A5CoreVcMt.o A6PhyObsev.o A7ListCons.o \
          A8MCOutput.o
QMCMain = CPMCSmuBgn.o CPMCWarmUp.o CPMCSwpMea.o CPMCGetMuT.o CPMCFixdnT.o
QMCSwep = SweepM2One.o SweepOne2M.o 

default: all

all: cpmod build-core build-dynm build-init build-main build-swep

build-core: ${QMCCore} 
build-dynm: $(QMCDynm)
build-init: $(QMCInit)
build-main: $(QMCMain)
build-swep: $(QMCSwep)

cpmod:
	cp ../Module/*.mod .

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o

clean-dat:
	rm -f *.dat
	rm -f *.bin.*
	rm -f *.out

clean-all: clean clean-dat
######################################################################################################################
######################################################################################################################
######################################################################################################################
