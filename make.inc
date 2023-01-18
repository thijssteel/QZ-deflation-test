# Project directory
INSTALLDIR := .

# fortran compiler and flags
FC := gfortran
# FFLAGS := -O0 -ggdb3 -fcheck=bounds -ffpe-trap=zero,invalid,overflow,underflow
FFLAGS := -O2 -march=native
INCFLAG := -J
#INCFLAG := -module

# libraries required by compiler (LAPACK AND BLAS)
LIBS := -lblas -llapack -lstdc++ -lm -lpthread

# Dependent directories
SRCDIR := $(INSTALLDIR)/src
OBJDIR := $(INSTALLDIR)/obj
INCDIR := $(INSTALLDIR)/include
BINDIR := $(INSTALLDIR)/bin