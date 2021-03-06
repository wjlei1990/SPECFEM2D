F90 = gfortran   # /opt/intel/fc/10.1.018/bin/ifort 
F90_FLAGS = -O3  # -O3 -$(MOD_FLAG) $(MOD_DIR)

MOD = wave2d_constants wave2d_variables wave2d_simul_type wave2d_define_der_matrices solve_eg wave2d_mesher wave2d_solver wave2d_sub
SUB = source_time_function  gll_library lagrange_poly  numerical_recipes

SRC_DIR = src
MOD_DIR = mod
OBJ_DIR = obj
BIN_DIR = bin
MAIN = wave2d
MOD_FLAG= module


MOD_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(MOD))
F90_OBJ = $(patsubst %,$(OBJ_DIR)/%.o,$(SUB))
OBJ = $(F90_OBJ) $(MOD_OBJ)

PROG:mesh specfem

mesh:
	cd src/mesh2d; make

specfem:
	cd src/specfem2d; make


.PHONY : clean

clean:
	rm -f $(OBJ_DIR)/* $(BIN_DIR)/* 

