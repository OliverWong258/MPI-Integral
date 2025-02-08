#include "Matrix.h"
using namespace std;


Matrix::Matrix(int row, int col) {//create a matrix without initialized elements
	nrow = row, ncol = col;
    int num = row*col;
	data = new double [num];
    for(int i = 0;i < num;++i){
        data[i] = 0.0;
    }
}

Matrix::~Matrix() {
	delete[] data;
}

double& Matrix::operator()(int i, int j){
	return this->data[i*this->ncol+j];
}

ostream& operator<<(ostream& os, Matrix& mat) {//output the elements of the matrix
	for (int i = 0;i < mat.nrow;++i) {
		for (int j = 0;j < mat.ncol-1;++j) {
			os << fixed << setprecision(4) << mat(i,j) << ", ";
		}
		os << mat(i,mat.ncol-1) << endl;
	}
	return os;
}

void Matrix::set_diago(){
    for(int i = 1;i < nrow;++i){
        for(int j = 0;j < i;++j){
            data[i*ncol+j] = data[j*ncol+i];
        }
    }
}

bool Matrix::isRSM(){
	//Mytimer::tick("Matrix", "isRSM");
	if(nrow != ncol)
		return 0;
	for(int i = 0;i < nrow-1;++i){
		for(int j = 1;j < ncol;++j){
			if(abs(data[i*ncol+j]-data[j*ncol+i]) > 1e-6)
				return 0;
		}
	}
	//Mytimer::tick("Matrix", "isRSM");
	return 1;
}

void Matrix::lapack_diago(){
	ofstream ofs1("../output/eigenvalues.log");
	ofstream ofs2("../output/eigenvectors.log");
	int num = nrow*ncol;
	double* tmp = new double[num];
	for(int i = 0;i < num;++i){
		tmp[i] = data[i];
	}
	double eigval[nrow];
	int info;
	info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nrow, tmp, nrow, eigval);
	if(info > 0){
		cout << "fail to compute eigenvalues" << endl;
		exit(1);
	}
	ofs1 << "eigenvalues: " << endl;
	for(int i = 0;i < nrow;++i){
		ofs1 << fixed << setprecision(6) << eigval[i] << endl;
	}
	ofs1.close();
	ofs2 << "eigenvectors(rowwise): " << endl;
	for(int j = 0;j < ncol;++j){
		for(int i = 0;i < nrow;++i){
			ofs2 << fixed << setprecision(6) << tmp[i*nrow+j] << ", ";
		}
		ofs2 << endl;
	}
	ofs2 << endl;
	ofs2.close();
	delete tmp;
}
