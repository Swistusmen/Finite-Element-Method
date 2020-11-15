#include "Solver.h"
#include <memory>
#include <Eigen/Dense>
#include <utility>
#include <cmath>
#include <vector>

Matrix4d& slv::Solver::getEtaMatrixG2()
{
	double ksi = 1 / sqrt(3);
	Matrix4d* matrix = new Matrix4d();

	for (int i = 0; i < 4; i++)
	{
		ksi *= -1;
		matrix->operator()(i, 0) = (ksi - 1) / 4;
		matrix->operator()(i, 1) = (-ksi - 1) / 4;
		matrix->operator()(i, 2) = (ksi + 1) / 4;
		matrix->operator()(i, 3) = (-ksi + 1) / 4;
	}

	return *matrix;
}

Matrix4d& slv::Solver::getKsiMatrixG2()
{
	double eta = 1 / sqrt(3);
	auto *matrix = new Matrix4d();
	for (int i = 0; i < 4; i++)
	{
		eta *= -1;
		matrix->operator()(i, 0) = (eta - 1) / 4;
		matrix->operator()(i, 1) = (-eta + 1) / 4;
		matrix->operator()(i, 2) = (eta + 1) / 4;
		matrix->operator()(i, 3) = (-eta - 1) / 4;
	}
	return *matrix;
}

Matrix2d& slv::Solver::getJacobyMatrix2(Matrix4d& eta, Matrix4d& ksi, double* x, double* y, int point)
{
	auto matrix = new Matrix2d();
	for (int i = 0; i < 4; i++)
		matrix->operator()(0, 0) += x[i] * ksi(point, i);
	for (int i = 0; i < 4; i++)
		matrix->operator()(0, 1) += y[i] * ksi(point, i);
	for (int i = 0; i < 4; i++)
		matrix->operator()(1, 0) += x[i] * eta(point, i);
	for (int i = 0; i < 4; i++)
		matrix->operator()(1, 1) += y[i] * eta(point, i);
	return *matrix;
}

Vector2d& slv::Solver::getDerivativeOfNByCoordinate_XY(Matrix2d& inversedJacoby, double detJ, Vector2d& derivatives) 
{
	auto vec = new Vector2d();
	vec->operator()(0) = (inversedJacoby(0, 0)*derivatives(0) + inversedJacoby(0, 1)*derivatives(1)) ; ///zgodnie ze wzorem powinno byc dzielenie przez detJ
	vec->operator()(1) = (inversedJacoby(1, 0)*derivatives(0) + inversedJacoby(1, 1)*derivatives(1)) ;
	return *vec; 
}

Vector2d& slv::Solver::getVectorOfDerivatives(Matrix4d& ksi, Matrix4d& eta, int fShape, int point)
{
	auto vec = new Vector2d();
	vec->operator()(0) = ksi(point, fShape);
	vec->operator()(1) = eta(point, fShape);
	return *vec;
}

Vector4d* slv::Solver::getXYDerivativesForPoint(Matrix4d& eta, Matrix4d& ksi, Matrix2d& jacoby, int point)
{
	Vector4d* xy = new Vector4d[2];
	for (int i = 0; i < 4; i++)
	{
		auto buffer = this->getDerivativeOfNByCoordinate_XY(jacoby.inverse(), jacoby.determinant(),
			this->getVectorOfDerivatives(ksi, eta, i, point));
		xy[0].operator()(i) = (buffer)(0);
		xy[1].operator()(i) = (buffer)(1);
	}
	return xy;
}


Matrix4d* slv::Solver::getHSumbatricies(Matrix4d& eta, Matrix4d& ksi, Matrix2d& jacoby, int point)
{
	auto result = getXYDerivativesForPoint(eta, ksi, jacoby, point);
	Matrix4d* H1 = new Matrix4d(vecAndvecTMultiplication(result[0]));
	Matrix4d* H2 = new Matrix4d(vecAndvecTMultiplication(result[1]));
	delete result;
	auto a = new Matrix4d(*H1 + *H2);
	return a;
}

Matrix4d* slv::Solver::getHMatrix(Matrix4d& eta, Matrix4d& ksi, Matrix2d& jacoby, double k)
{
	double detJ = jacoby.determinant();
	Matrix4d* H = new Matrix4d();
	for (int i = 0; i < 4; i++)
	{
		auto buffer= this->getHSumbatricies(eta, ksi, jacoby, i);
		*buffer *= k*detJ;
		(*H) += buffer;
	}
	return H;
}