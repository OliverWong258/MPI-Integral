#include "Read.h"
#include "Calc.h"
#include"Matrix.h"
#include "openmpi-4.1.4/ompi/include/mpi.h"
#include <cmath>
#include <iostream>
#include<iomanip>
#include<fstream>
#include <stdlib.h>
#include <cstdio>
#include<omp.h>
using namespace std;

extern "C" {
    void pdsyevd_(char* jobz, char* uplo, int* n, double* a, int* ia, int* ja, int* desca, double* w, double* z, int* iz, int* jz, int* descz, double* work, int* lwork, int* iwork, int* liwork, int* info);
    void descinit_(int* desc, int* m, int* n, int* mb, int* nb, int* irsrc, int* icsrc, int* ictxt, int* lld, int* info);
    void Cblacs_pinfo(int* mypnum, int* nprocs);
    void Cblacs_get(int ConTxt, int what, int* val);
    void Cblacs_gridinit(int* ConTxt, char* order, int nprow, int npcol);
    void Cblacs_gridinfo(int ConTxt, int* nprow, int* npcol, int* myrow, int* mycol);
    int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
}

const int tag = 99;

int main(){
    // Initialize MPI
    MPI_Init(NULL, NULL);

    double start = 0.0, end = 0.0;

    int np, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int mpi_scalapack = 0;
    int n = 0;
    double* M = nullptr;

    if(rank == 0){
        MyRead read_in;
        cout << "Enter the input file path" << endl;
        string input_path;
        cin >> input_path;
        if(input_path == "default")
            input_path = "../input/INPUT.txt";
        read_in.ReadInput(input_path);
        read_in.ReadPoints();
        read_in.ReadV();
        read_in.ReadDistr();

        cout << "Enter the num of threads" << endl;
        int nthreads;
        cin >> nthreads;
        omp_set_num_threads(nthreads);
        start = MPI_Wtime();//timer begins 

        long long V_num = read_in.nx * read_in.ny * read_in.nz;
        double sum = 0.0;
        Matrix mat(read_in.point_num, read_in.point_num);

        double *m = new double[read_in.mesh];
        Calc::capway(m, read_in.mesh, read_in.dr, read_in.distr);

        double dv = (read_in.lx/read_in.nx) * (read_in.ly/read_in.ny) * (read_in.lz/read_in.nz);

        int vpoint_num = 200000000 / read_in.point_num;
        int total_round = read_in.nx*read_in.ny*read_in.nz / (vpoint_num);
        int rest_num = (read_in.nx*read_in.ny*read_in.nz) % (vpoint_num);

        double **ff = new double*[read_in.point_num];
        for(int i = 0;i < read_in.point_num;++i){
            ff[i] = new double[vpoint_num];
        }

        for(int ii = 0;ii < total_round;++ii){
            //init f
                for(int i = 0;i < read_in.point_num;++i){
                    #pragma omp parallel for schedule(static)
                    for(long long k = 0;k < vpoint_num;++k){
                        double f = 0.0, t = 0.0;
                        long long index = ii*vpoint_num+k;

                        long long x = 0, y = 0, z = 0;
		                x = index / (read_in.ny*read_in.nz);
		                y = (index-read_in.ny*read_in.nz*x) / read_in.nz;
		                z = index - read_in.ny*read_in.nz*x - read_in.nz*y;

                        t = Calc::dist(read_in.point_vec[i], double(x) * double(read_in.lx / read_in.nx), double(y) * double(read_in.ly / read_in.ny), 
                            double(z) * double(read_in.lz/read_in.nz));

                        f = Calc::spline(read_in.mesh, read_in.dr, m, read_in.distr, t);

                        ff[i][k] = f;
                    }
                }


            for(int i = 0;i < read_in.point_num;++i){
                for(int j = i;j < read_in.point_num;++j){
                    sum = 0;
                    if(Calc::dist(read_in.point_vec[i], read_in.point_vec[j]) <= read_in.cutoff*2){
                        #pragma omp parallel for schedule(static) reduction(+:sum)
                        for(long long k = 0;k < vpoint_num;++k){
                                double fa=0.0, fb=0.0, t = 0.0;
                                fa = ff[i][k];
                                fb = ff[j][k];
                                sum += fa * read_in.V[ii*vpoint_num+k] * fb * dv;
                        }
                    }
                
                    mat(i,j) += sum;
                }
            }
        }

        //deal with the rest
        //init f
        for(int i = 0;i < read_in.point_num;++i){
            #pragma omp parallel for schedule(static)
            for(long long k = 0;k < rest_num;++k){
                double f = 0.0, t = 0.0;
                long long index = total_round*vpoint_num+k;

                long long x, y, z;
		        x = index / (read_in.ny*read_in.nz);
                y = (index-read_in.ny*read_in.nz*x) / read_in.nz;
                z = index - read_in.ny*read_in.nz*x - read_in.nz*y;

                t = Calc::dist(read_in.point_vec[i], double(x) * double(read_in.lx / read_in.nx), double(y) * double(read_in.ly / read_in.ny), 
                    double(z) * double(read_in.lz/read_in.nz));

                f = Calc::spline(read_in.mesh, read_in.dr, m, read_in.distr, t);

                ff[i][k] = f;
            }
        }
        for(int i = 0;i < read_in.point_num;++i){
            for(int j = i;j < read_in.point_num;++j){
                sum = 0;

                if(Calc::dist(read_in.point_vec[i], read_in.point_vec[j]) <= read_in.cutoff*2){
                    #pragma omp parallel for schedule(static) reduction(+:sum)
                    for(long long k = 0;k < rest_num;++k){
                        double fa=0.0, fb=0.0, t = 0.0;
                        fa = ff[i][k];
                        fb = ff[j][k];
                        sum += fa * read_in.V[total_round*vpoint_num+k] * fb * dv;
                    }
                }
            
                mat(i,j) += sum;
            }
        }

        mat.set_diago();
        /*
        cout << "The Matrix is:" << endl;
        cout << mat;
        */

        if(mat.isRSM()){
            if(read_in.arg2meth["diago_lib"] == "lapack"){
                mat.lapack_diago();
                MPI_Bcast(&mpi_scalapack, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
            else if(read_in.arg2meth["diago_lib"] == "scalapack"){
                mpi_scalapack = 1;
                n = read_in.point_num;
                M = new double[n*n];
                for(int ii = 0;ii < n;++ii)
                    for(int jj = 0;jj < n;++jj)
                        M[ii*n+jj] = mat.data[ii*n+jj];

                MPI_Bcast(&mpi_scalapack, 1, MPI_INT, 0, MPI_COMM_WORLD);
            }
            else{
                cout << "ERROR, I do not use mylib" << endl;
                exit(1);
            }
        }
        else{
            cout << "Error, the Matrix is not diago" << endl;
            exit(1);
        }
    }
    else{
        MPI_Bcast(&mpi_scalapack, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if(mpi_scalapack){
        //init
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(rank != 0)
            M = new double[n*n];
        MPI_Bcast(M, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Set the grid parameters
            int nprow = sqrt(np); // Number of process rows
            int npcol = 0;
            for(;nprow >= 1;--nprow){
                npcol = int(np/nprow);
                if(np == nprow*npcol)
                    break;
            }

            // Get the BLACS context
            int mypnum, nprocs;
            int myrow, mycol;
            int ictxt;
            Cblacs_pinfo(&mypnum, &nprocs);
            Cblacs_get(-1, 0, &ictxt);
            Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
            Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

            // Set the grid parameters
            int nb = 1; // The block size
            int zero = 0; // Base index in C
            int one = 1; // Base index in Fortran
            int descA[9], descZ[9]; // The array descriptors

            // Initialize the array descriptors
            int info;
            // Allocate and initialize the matrix
            // Calculate the local dimensions
            int mA = numroc_(&n, &nb, &mypnum, &zero, &nprow);
            int nA = numroc_(&n, &nb, &mypnum, &zero, &npcol);

            descinit_(descA, &n, &n, &nb, &nb, &zero, &zero, &ictxt, &mA, &info);
            descinit_(descZ, &n, &n, &nb, &nb, &zero, &zero, &ictxt, &mA, &info);

            // Allocate local arrays
            double* A = new double[mA * nA];

            // Distribute the matrix
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int prow = (i / nb) % nprow;
                    int pcol = (j / nb) % npcol;
                    if (myrow == prow && mycol == pcol)
                    {
                        int local_i = (i / nb) / nprow * nb + i % nb;
                        int local_j = (j / nb) / npcol * nb + j % nb;
                        A[local_i + local_j * mA] = M[i * n + j];
                    }
                }
            }

            // Allocate the array for the eigenvalues
            double* W = new double[n];
            double* V = new double[n*n];

            // Allocate and initialize the matrix for the eigenvectors
            double* Z = new double[mA * nA];
            for (int i = 0; i < mA * nA; i++)
                Z[i] = 0.0;

            // Allocate the work arrays
            int lwork = -1;
            int liwork = -1;
            double workSize;
            int iworkSize;
            pdsyevd_("V", "U", &n, A, &one, &one, descA, W, Z, &one, &one, descZ, &workSize, &lwork, &iworkSize, &liwork, &info);
            lwork = (int)workSize;
            liwork = iworkSize;
            double* work = new double[lwork];
            int* iwork = new int[liwork];

            // Compute the eigenvalues and eigenvectors
            pdsyevd_("V", "U", &n, A, &one, &one, descA, W, Z, &one, &one, descZ, work, &lwork, iwork, &liwork, &info);

            // Handle errors
            if (info != 0) {
                std::cout << "pdsyevd failed with exit code: " << info << std::endl;
            }

            if(rank != 0){
                MPI_Send(Z, mA*nA, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                MPI_Send(&myrow, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
                MPI_Send(&mycol, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
            }
            else{
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        int prow = (i / nb) % nprow;
                        int pcol = (j / nb) % npcol;

                        if (myrow == prow && mycol == pcol){
                            int local_i = (i / nb) / nprow * nb + i % nb;
                            int local_j = (j / nb) / npcol * nb + j % nb;
                            V[i*n+j] = Z[local_i + local_j * mA];
                        }
                    }
                }

                MPI_Status status;
                int myrow, mycol;
                for(int k = 1;k < np;++k){
                    MPI_Recv(Z, mA*nA, MPI_DOUBLE, k, tag, MPI_COMM_WORLD, &status);
                    MPI_Recv(&myrow, 1, MPI_INT, k, tag, MPI_COMM_WORLD, &status);
                    MPI_Recv(&mycol, 1, MPI_INT, k, tag, MPI_COMM_WORLD, &status);

                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                            int prow = (i / nb) % nprow;
                            int pcol = (j / nb) % npcol;

                            if (myrow == prow && mycol == pcol)
                            {
                                int local_i = (i / nb) / nprow * nb + i % nb;
                                int local_j = (j / nb) / npcol * nb + j % nb;
                                V[i*n+j] = Z[local_i + local_j * mA];
                            }
                        }
                    }
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);

            // Only process with rank 0 outputs the eigenvalues
            if (rank == 0) {
                ofstream ofs2("../output/eigenvectors.log");
	            ofstream ofs1("../output/eigenvalues.log");
                ofs1 << "eigenvalues: " << endl;
                for(int i = 0;i < n;++i){
                    ofs1 << fixed << setprecision(6) << W[i] << std::endl;
                }
                ofs1.close();
                ofs2 << "eigenvectors(rowwise): " << endl;
                for(int i = 0;i < n;++i){
                    for(int j = 0;j < n;++j){
                        ofs2 << fixed << setprecision(6) << V[i*n+j] << ", ";
                    }
                    ofs2 << std::endl;
                }
                ofs2.close();
            }

            // Clean up
            delete[] A;
            delete[] W;
            delete[] V;
            delete[] Z;
            delete[] work;
            delete[] iwork;
    }

    if(rank == 0){
        end = MPI_Wtime();//timer ends
        double tmp = end - start;
        std::cout << "complete" << std::endl;
        std::cout << "Total time: " << tmp << 's' << endl;
    }
    MPI_Finalize();

    return 0;
}