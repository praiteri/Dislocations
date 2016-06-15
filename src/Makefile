#####################################################################
# Copyright (C) 2016 by Paolo Raiteri                               #
#                                                                   #
# p.raiteri@curtin.edu.au                                           #
#                                                                   #
# This program is free software; you can redistribute it and/or     #
# modify it under the terms of the GNU General Public License       #
# as published by the Free Software Foundation; either version 2    #
# of the License, or (at your option) any later version.            #
#                                                                   #
# This program is distributed in the hope that it will be useful,   #
# but WITHOUT ANY WARRANTY; without even the implied warranty of    #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      #
# GNU General Public License for more details.                      #
#                                                                   #
# For the full text of the GNU General Public License,              #
# write to: Free Software Foundation, Inc.,                         #
#           675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                   #
# The GNU GPL can also be found at http://www.gnu.org               #
#                                                                   #
# No claim is made that this program is free from errors and        #
# no liability will be accepted for any loss or damage that         #
# may result. The user is responsible for checking the validity     #
# of their results.                                                 #
#####################################################################

# Linux-specific options
VERSION=1.0
EXEDIR = ../bin
EXE = disloc.x

# Fortran source code
SRCF90 = $(wildcard *F90)
OBJ=$(SRCF90:.F90=.o)

CPPFLAGS=

# Compiling options
help:
	@echo '--------------------------------------------------'
	@echo 'gf         : gfortran compiler - standard build'
	@echo 'ifort      : intel fortran compiler - standard build'
	@echo 'gf-dbg     : gfortran compiler - debugging build'
	@echo 'ifort-dbg  : intel fortran compiler - debugging build'
	@echo 'tar        : create tarball of the sources'
	@echo 'clean      : clean build'
	@echo 'clean-all  : clean everything'
	@echo '--------------------------------------------------'

# Start up
start: 
	@if [ ! -d ../bin ]; then mkdir ../bin; fi 
	@if [ ! -d ../obj ]; then mkdir ../obj; fi 
	-mv -f ../obj/*.mod ../obj/*.o ./ 

# Executables
gf: start
	$(MAKE) FC="gfortran" CC="gcc" \
	FFLAGS="-funroll-all-loops -O3 -ffast-math -fno-second-underscore" \
	CFLAGS="-funroll-all-loops -O3 -ffast-math" \
	$(EXE) 
	mv $(EXE) $(EXEDIR)/

gf-dbg: start
	$(MAKE) FC="gfortran" CC="gcc" \
	CPPFLAGS="-DDEBUG" \
	FFLAGS="-fbounds-check -fbacktrace -fno-second-underscore -ffpe-trap=zero,overflow,invalid" \
	CFLAGS="-fbounds-check" \
	$(EXE)
	mv $(EXE) $(EXEDIR)/

ifort: start
	$(MAKE) FC="ifort" CC="icc" \
	FFLAGS="-O3 -ip -w" \
	CFLAGS="-O3 -ip -w" \
	$(EXE)
	mv $(EXE) $(EXEDIR)/

ifort-dbg: start
	$(MAKE) FC="ifort" CC="icc" \
	CPPFLAGS="-DDEBUG" \
	FFLAGS="-debug all -check all -traceback" \
	CFLAGS="-O3 -ip -w" \
	$(EXE)
	mv $(EXE) $(EXEDIR)/

clean:
	rm -fr ../obj ./*.o ./*.mod 

clean-all:
	rm -fr ../obj ./*.o ./*.mod ../bin/disloc.x 

tar:    clean-all
	cd ../../ ; tar -zcf disloc.tgz disloc/src disloc/bin disloc/test ; cd - ; mv ../../disloc.tgz ../
#
# Executable generation
#
$(EXE): $(OBJ) $(OBJC)
	$(FC) $(CPPFLAGS) $(FFLAGS) $(OBJ) -o $(EXE)
	-mv *.mod *.o ../obj/.
#
# Fortran Dependencies
#
disloc.o: ./disloc.F90 variables.o
	$(FC) $(CPPFLAGS) $(FFLAGS) -c	./disloc.F90
rw_pdb.o: ./rw_pdb.F90 variables.o
	$(FC) $(CPPFLAGS) $(FFLAGS) -c	./rw_pdb.F90
tools.o: ./tools.F90 variables.o
	$(FC) $(CPPFLAGS) $(FFLAGS) -c	./tools.F90
variables.o: ./variables.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c	./variables.F90
