#pragma once
#include<iomanip>
#include<iostream>
#include<algorithm>
#include<fstream>

#include "lapack-3.11/CBLAS/include/cblas.h"
#include "lapack-3.11/LAPACKE/include/lapacke.h"

using namespace std;

class Matrix {
public:
	double* data = nullptr;//store elements of the matrix
	int nrow = 0, ncol = 0;//the row number nad column number
	Matrix(int row, int col);
	~Matrix();
	double& operator()(int i, int j);
    void set_diago();
	bool isRSM();
	void lapack_diago();
};

ostream& operator<<(ostream& os, Matrix& mat);
