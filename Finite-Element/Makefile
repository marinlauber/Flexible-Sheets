# Makefile

# Fortran compiler
FC = gfortran

# Suffixes rules to control how .o and .mod files are treated, inference rules
.SUFFIXES:
.SUFFIXES: .f90 .o .mod 

.f90.o:
	$(FC) $(FFLAGS) -c $<
.f90.mod:
	$(FC) $(FFLAGS) -c $<
#
# Declare all the separate files in one as OBJ
# These are the files created by fsplit and are what the executable depends on
#
OBJ = finiteelement.o
#
##########################################################
#
# Rules for making executable and cleaning directory
#
# This rule says that the executable md depends on the object files in OBJ
# and that they need to be complied then linked to produce the executable.
#	$@ is the target of the rule
#	$^ is all the dependencies of each rule
#	$< only the first element of the dependencies
# 
.PHONY: main
main: $(OBJ)
	$(FC) $(OBJ) -o $@
	./$@
	python *.py
#
# This is a rule for cleaning out object files and the executable.
# If the comand "make clean" is typed it will work
#
.PHONY: clean
clean:
	rm *.o *.mod *.csv main
#