#Directories
PWD    = $(shell pwd)
OBJDIR = ../../obj
SRCDIR = $(PWD)
BINDIR = ../../bin
#LIBDIR = $(PWD)/UTILS/lib
#IDIR   = $(PWD)/UTILS/include
#EVALLIB=/home/lei/ASDF/evalresp/.libs

## set ADIOS_DIR here or before doing make
#override ADIOS_DIR:=/home/lei/bin/adios-1.5.0
#override ADIOS_INC:=` ${ADIOS_DIR}/bin/adios_config -c -f`
#override ADIOS_FLIB:=`${ADIOS_DIR}/bin/adios_config -l -f`

############################
#compiler option
#OPT = -I${SHARED}
#OPT = -std03
FC = ifort
CC = icc
MPIFC = mpif90
MPICC = mpicc
CFLAGS= -g -O0 -check all -free

#############################
#link TT and SAC library
#TAULIBDIR=$(PWD)/ttimes_mod
#SACLIBDIR=${SACHOME}/lib
#LIBS= -lsacio -lsac -ltau -lm

############################
#user defined library
#LIB_FLEXWIN = libflexwin.a
#LIB_IO = libio.a
#LIB_VAR = libvar.a

#src file
#SRC= asdf_subs.f90 asdf_data.f90 copy_float.c distaz.f user_parameters.f90  seismo_subs.f90 user_functions.f90 maxima.f90 travel_times.f90 io_subs.f90 select_windows_stalta2.f90 flexwin.f90 FLEXWIN_ADIOS.f90

_OBJ = constants.o lagrange_poly.o gll_library.o exit_mpi.o \
			 wave2d_variables.o wave2d_define_der_matrices.o \
			 read_mesh.o assemble_MPI.o wave2d_simul_type.o\
			 wave2d_io_subs.o source_time_function.o\
			 set_src_and_rec.o wave2d_solver.o main.o\
	
OBJ = $(patsubst %, $(OBJDIR)/%,$(_OBJ))

##########################################################
PROG = $(BINDIR)/xspecfem2d 
default: ${PROG}

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)

$(OBJDIR)/%.o: $(SRCDIR)/%.f
	$(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)
.F.o:
$(OBJDIR)/%.o: $(SRCDIR)/%.F
	$(MPIFC) ${CFLAGS} -c -o $@ $< -module $(OBJDIR) $(ADIOS_INC)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(MPICC) -c -o $@ $< 

${PROG}: $(OBJ)
	${MPIFC} ${CFLAGS} -o $@ $^ -I$(OBJDIR)

.PHONY:clean print_var cleanall

print_var:
	echo $(OBJ)

clean:
	rm -f  ${LIB_ALL} ${PROGS} *.o *.mod *.a $(OBJDIR)/*

cleanall:
	rm -f  iasp91.*
	cd ${TAULIBDUR} ; make -f make_gfortran clean

