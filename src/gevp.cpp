/*
 * gevp.cpp
 *
 *  Created on: Apr 24, 2014
 *      Author: Thibaut Metivet
 */

 #include "gevp.hpp"
 #include "Z001.hpp"

 #include <Eigen/Eigenvalues>
 #include <Eigen/SVD> 
 #include <complex>

 using namespace LQCDA;

 // put a line of reals into the gevp form (square complex matrix)
 template<typename Derived>
 Matrix<std::complex<double>> gevp_form(const MatrixExpr<Derived>& mat)
 {
 	int rows = sqrt(mat.cols() / 2);
 	Matrix<std::complex<double>> res(rows, rows);
 	FOR_MAT(res, i, j)
 	{
 		index_t n = i*rows + j;
 		res(i, j) = std::complex<double>(mat(0, 2*n), mat(0, 2*n+1));
 	}

 	return res;
 }

 double real(const std::complex<double>& c)
 {
 	return c.real();
 }

 Matrix<std::complex<double>> gev(const Matrix<double>& mat, unsigned int t0)
 {
 	int Nt = mat.rows();
 	unsigned int ngev = sqrt(mat.cols() / 2);

 	Matrix<std::complex<double>> result(Nt, ngev);

 	Matrix<std::complex<double>> Ct0_inv = gevp_form(mat.row(t0)).inverse();
 	Matrix<std::complex<double>> Ct;

 	for(int t = 0; t < Nt; ++t)
 	{
 		Ct = gevp_form(mat.row(t));
 		Eigen::ComplexEigenSolver<Matrix<std::complex<double>>> eigsolv(Ct0_inv * Ct, false);
 		result.row(t) = eigsolv.eigenvalues().transpose();
 		// std::cout << Ct << '\n';
 	}

 	return result;
 }

 Matrix<double> gev_svd(const Matrix<double>& mat, unsigned int t0)
 {
 	int Nt = mat.rows();
 	unsigned int ngev = sqrt(mat.cols() / 2);

 	Matrix<double> result(Nt, ngev);

 	Matrix<std::complex<double>> Ct0_inv = gevp_form(mat.row(t0)).inverse();
 	Matrix<std::complex<double>> Ct;

 	for(int t = 0; t < Nt; ++t)
 	{
 		Ct = gevp_form(mat.row(t));
 		// Eigen::ComplexEigenSolver<Matrix<std::complex<double>>> eigsolv(Ct0_inv * Ct, false);
 		Eigen::JacobiSVD<Matrix<std::complex<double>>> svd(Ct0_inv * Ct);
 		// result.row(t) = eigsolv.eigenvalues().transpose();
 		result.row(t) = svd.singularValues().transpose();
 		// std::cout << Ct << '\n';
 	}

 	return result;
 }

 double deltaM2_COM(double W, double MPi, double L)
 {
 	double E_COM = W;
    double q2 = (E_COM*E_COM/4. - MPi*MPi)*(L*L)/(4.*M_PI*M_PI);
    double delta;
    if(q2 < 0.) {
	delta = 0.;
	std::cout << "q2 < 0 !" << std::endl;
    }
    else
	delta = -8.*M_PI*M_SQRTPI/(L*L*L)*q2/E_COM*Z001::z001q2(q2);
    std::cout << "q2 = " << q2 << std::endl;
    return delta;
 }