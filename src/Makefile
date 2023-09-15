######################################################################################################################
########################################## Information of this file ##################################################
######################################################################################################################
##
## Introduction
## ============
##
## Makefile for the finite-T CPMC program of t1-t2-U-V model (generalized Hubbard model) on square lattice.
## This Makefile is used to compile the whole program consisted of .f90 files from five different folders:
## General --> Some subroutines as connection to MKL math library, without dependence on modules;
## Module  --> Some modules to define all the related quantities for ftdqmc simulations;
## Special --> All the core subroutines performing the simulations;
## CPMC    --> Some subroutines to divide the simulation into five main steps;
## Main    --> The main program controlling the whole simulation.
##
## Usage
## =====
##
## 'make' or 'make clean'
##
## Dependence
## ==========
##
## This file depends on the "make.sys" file in the same folder as this file.
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
include ./make.sys

all: Compile TotalExe ## RmfModO

Compile:
	(cd General; $(MAKE))
	(cd Module;  $(MAKE))
	(cd Special; $(MAKE))
	(cd CPMC;    $(MAKE))
	(cd Main;    $(MAKE))

AllOFile = General/*.o  Module/*.o  Special/*.o  CPMC/*.o  Main/*.o

TotalExe:
	$(LINKER) $(AllOFile) -o ftcpmc $(LFLAGS) $(LIBS)

RmfModO:
	(cd General; $(MAKE) clean)
	(cd Module;  $(MAKE) clean)
	(cd Special; $(MAKE) clean)
	(cd CPMC;    $(MAKE) clean)
	(cd Main;    $(MAKE) clean)

clean:
	(cd General; $(MAKE) clean)
	(cd Module;  $(MAKE) clean)
	(cd Special; $(MAKE) clean)
	(cd CPMC;    $(MAKE) clean)
	(cd Main;    $(MAKE) clean)
	rm -rf ftcpmc
	rm -rf Add_Output/
	rm -rf Output/*

clean-dat:
	rm -rf *.dat
	rm -rf *.bin.*
	rm -rf *.out

clean-all: clean clean-dat
######################################################################################################################
######################################################################################################################
######################################################################################################################
