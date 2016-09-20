#ifndef SPECTRAL_H
#define SPECTRAL_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include "utils.h"
#include <functional>
#include "MersenneTwister.h"
#include <iomanip>


class Spectral{
public:
	Spectral();
	void load(double h);
	void make_states(int charge);
	void check(Eigen::VectorXd eigvals);
	
private:
	double h;
	int N,rows;
	Eigen::MatrixXd eigvecs,Hnn;
	Eigen::VectorXd eigvals;
	bool loading;
	vector<double> alphax, alphaz;
	vector<int> states;

	void HeisenbergXZ(double *v, double *w, double J);
	void HeisenbergZ(double *v, double *w, double J);
	void Sz(double *v, double *w, int site);
	void makeDense(function<void(double  *v, double *w)> matvec, Eigen::MatrixXd &EigenDense );
	void makeSparse(function<void(double  *v, double *w)> matvec, Eigen::SparseMatrix<double> &EigenSparse );
	
	int next(int);
};

#endif
