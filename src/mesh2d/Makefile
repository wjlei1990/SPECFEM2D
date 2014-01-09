CC=gcc
FC=gfortran
CFLAGS=-g -lm
FFLAGS=-g -lm

all_obj=mesh_par.o build_mesh_interface_subs.o calcu_xigll_zigll.o init_variables.o  mesh2d_subs.o read_parameters.o save_mesh.o gll_library.o main.o

PROG=mesh2d
default: ${PROG}

mesh_par.o :mesh_par.f90
	${FC} ${FFLAGS} -c -o $@ $<

build_mesh_interface_subs.o :build_mesh_interface_subs.f90
	${FC} ${FFLAGS} -c -o $@ $<

calcu_xigll_zigll.o : calcu_xigll_zigll.f90
	${FC} ${FFLAGS} -c -o $@ $<

init_variables.o :init_variables.f90
	${FC} ${FFLAGS} -c -o $@ $<

mesh2d_subs.o :mesh2d_subs.f90
	${FC} ${FFLAGS} -c -o $@ $<

read_parameters.o :read_parameters.f90
	${FC} ${FFLAGS} -c -o $@ $<

save_mesh.o :save_mesh.f90
	${FC} ${FFLAGS} -c -o $@ $<

gll_library.o :gll_library.f90
	${FC} ${FFLAGS} -c -o $@ $<

main.o :main.f90
	${FC} ${FFLAGS} -c -o $@ $<

${PROG} :${all_obj}
	${FC} ${FFLAGS} -o $@ ${all_obj}

.PHONY :clean

clean:
	rm -f *.o ${PROG}
