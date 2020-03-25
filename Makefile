#FC = ifort
FC = gfortran
FFLAGS = -O3 -g

DAM_DIR= $(CURDIR)/

GLOBAL  = TDAMGLOBAL16.F90
BUILD   = TDAMBUILD16.F90
COMMON  = TDAMCOM16.F90
DAMPOT  = TDAMPOT16.F90
LBFGSB  = lbfgsb.f
LINPACK = linpack.f
BLAS    = blas.f
TIMER   = timer.f
HESSIAN = HESSIAN.F90

all : MESPIMIZER

MESPIMIZER : $(GLOBAL) $(BUILD) $(COMMON) $(DAMPOT) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(HESSIAN) 
	$(FC) $(FFLAGS) $(GLOBAL) $(BUILD) $(COMMON) $(DAMPOT) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(HESSIAN) -o DAMESPIMIZER.exe

