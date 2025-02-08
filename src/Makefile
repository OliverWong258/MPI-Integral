cc = g++

mpi_cc = mpicxx

prom = final_openmp

ompFLAG = -fopenmp

LFLAGS = -llapacke -llapack -lcblas -lrefblas -lm -lgfortran

SFLAGS = -lscalapack -lblas

LAPACK_DIR = lapack-3.11
LIB_DIR = $(LAPACK_DIR)/lib
CBLASH_DIR = $(LAPACK_DIR)/CBLAS/include
LAPACKEH_DIR = $(LAPACK_DIR)/LAPACKE/include

SCALAPACK_DIR = scalapack-2.2.0

omp_mpi: *.cpp *.h $(LIB_DIR)/libcblas.a $(LIB_DIR)/liblapack.a $(LIB_DIR)/liblapacke.a $(LIB_DIR)/librefblas.a $(SCALAPACK_DIR)/libscalapack.a
	$(mpi_cc) -D __MPI -o $(prom) *.cpp -I$(CBLASH_DIR) -I$(LAPACKEH_DIR) -L$(LIB_DIR) -L$(SCALAPACK_DIR) $(SFLAGS) $(LFLAGS) $(ompFLAG)
	 
clean:
	rm -f final_openmp
