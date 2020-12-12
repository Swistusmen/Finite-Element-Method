#include "Solver.h"
#include <memory>
#include <Eigen/Dense>
#include <utility>
#include <cmath>
#include <vector>
#include <array>
#include "InputData.h"

namespace {
	std::vector<double> generateGaussIntegrationalPoints(size_t scheme) {
		std::vector<double> points;
		switch (scheme) {
		case 2: {
			points.push_back(-1 / sqrt(3.0));
			points.push_back(1 / sqrt(3.0));
		}break;
		case 3: {
			points.push_back(-std::sqrt(3.0 / 5.0));
			points.push_back(0.0);
			points.push_back(std::sqrt(3.0 / 5.0));
		}break;
		case 4: {
			points.push_back(-1.0*std::sqrt(3.0 / 7.0 + 2.0 / 7.0*std::sqrt(1.2)));
			points.push_back(-1.0*std::sqrt(3.0 / 7.0 - 2.0 / 7.0*std::sqrt(1.2)));
			points.push_back(std::sqrt(3.0 / 7.0 - 2.0 / 7.0*std::sqrt(1.2)));
			points.push_back(std::sqrt(3.0 / 7.0 + 2.0 / 7.0*std::sqrt(1.2)));
		}break;
		}
		return points;
	}
}

namespace slv {
	LocalOperations::LocalOperations()
	{
		Ksi.push_back([](double a) {return -0.25*(1 - a); });
		Ksi.push_back([](double a) {return 0.25*(-a + 1); });
		Ksi.push_back([](double a) {return 0.25*(a + 1); });
		Ksi.push_back([](double a) {return -0.25*(1 + a); });
		Eta.push_back([](double a) {return -0.25*(1 - a); });
		Eta.push_back([](double a) {return -0.25*(1 + a); });
		Eta.push_back([](double a) {return 0.25*(a + 1); });
		Eta.push_back([](double a) {return 0.25*(-a + 1); });	
		Shape.push_back([](double e, double n) {return 0.25*(1 - e)*(1 - n); });
		Shape.push_back([](double e, double n) {return 0.25*(1 + e)*(1 - n); });
		Shape.push_back([](double e, double n) {return 0.25*(1 + e)*(1 + n); });
		Shape.push_back([](double e, double n) {return 0.25*(1 - e)*(1 + n); });
	}

	Solver::Solver(int points)
	{
		this->gaussIntegralScheme = points;
		double *tempData = new double[points];
		if (points == 2) {
			tempData[0] = 1;
			tempData[1] = 1;
		}
		else if (points == 3) {
			tempData[0] = 5.0 / 9.0;// 0.55;
			tempData[1] = 8.0 / 9.0;// 0.88;
			tempData[2] = 5.0 / 9.0;// 0.55;
		}
		else {
			tempData[0] = (18.0 - std::sqrt(30.0)) / 36.0;// 0.347;
			tempData[1] = (18.0 + std::sqrt(30.0)) / 36.0;//0.652;
			tempData[2] = (18.0 + std::sqrt(30.0)) / 36.0;//0.652;
			tempData[3] = (18.0 - std::sqrt(30.0)) / 36.0;// 0.347;
		}
		for (size_t i = 0; i < points; ++i)
		{
			for (size_t j = 0; j < points; ++j)
			{
				this->wages.push_back(tempData[i] * tempData[j]);
			}
		}
		delete tempData;

		Eta = (getMatrixOfLocalTransformation(LocalType::ETA));
		Ksi = (getMatrixOfLocalTransformation(LocalType::KSI));
		Shape = (getMatrixOfLocalTransformation(LocalType::SHAPE));
	}

	MatSPtr Solver::getMatrixOfLocalTransformation(LocalType type)
	{
		int noIntPoints = this->gaussIntegralScheme;
		size_t size = std::pow(this->gaussIntegralScheme, 2);
		auto matrix = std::make_shared<Matrix>(4, size);
		std::vector<double> points= std::move(generateGaussIntegrationalPoints(noIntPoints)) ;
		switch (type) {
		case LocalType::ETA: case LocalType::KSI: {
			auto functions = localOperations.getFunctions(type);
			for (size_t i = 0; i < size; i++)
			{
				for (size_t j = 0; j < 4; j++)
					matrix->operator()(j, i) = functions.at(j)(points.at(i % noIntPoints));
			}
		}break;
		case LocalType::SHAPE: {
			for (size_t i = 0; i < size; i++)
			{
				for (size_t j = 0; j < 4; j++)
					matrix->operator()(j, i) = (localOperations.Shape.at(j))(points.at(i % noIntPoints),
						points.at(static_cast<int>(i/noIntPoints)));
			}
		}break;
		}
		return matrix;
	}

	MatSPtr Solver::getJacobyMatrix2( double* x, double* y, int point)
	{
		auto matrix = std::make_shared<Matrix>(2);
		const size_t size = 4;
		for (size_t i = 0; i < size; i++)
			matrix->operator()(0, 0) += x[i] * (*Ksi).operator()(i,point);
		for (size_t i = 0; i < size; i++)
			matrix->operator()(0, 1) += y[i] * (*Ksi).operator()(i,point);
		for (size_t i = 0; i < size; i++)
			matrix->operator()(1, 0) += x[i] * (*Eta).operator()(i,point);
		for (size_t i = 0; i < size; i++)
			matrix->operator()(1, 1) += y[i] * (*Eta).operator()(i,point);

		return matrix;
	}

	VecSPtr Solver::getDerivativeOfNByCoordinate_XY(Matrix& inversedJacoby, double detJ, VecSPtr& derivatives)
	{
		auto vec = std::make_shared<Vector>(2);
		vec->operator()(0) = (inversedJacoby(0, 0)*derivatives->operator()(0) + inversedJacoby(0, 1)*derivatives->operator()(1)); ///zgodnie ze wzorem powinno byc dzielenie przez detJ
		vec->operator()(1) = (inversedJacoby(1, 0)*derivatives->operator()(0) + inversedJacoby(1, 1)*derivatives->operator()(1));
		return vec;
	}

	VecSPtr Solver::getVectorOfDerivatives(int fShape, int point)
	{
		auto vec = std::make_shared<Vector>(2);
		vec->operator()(0) = Ksi->operator()(fShape, point);
		vec->operator()(1) = Eta->operator()(fShape, point);
		return vec;
	}

	std::vector<VecSPtr> Solver::getXYDerivativesForPoint( MatSPtr& jacoby, int point)
	{
		std::vector<VecSPtr> xy;
		xy.push_back(std::make_shared<Vector>(4));
		xy.push_back(std::make_shared<Vector>(4));
		const size_t size = std::pow(this->gaussIntegralScheme, 2);
		for (int i = 0; i < 4; i++)
		{
			auto mat=std::make_unique<Matrix>(*jacoby);
			inverseMat2(*mat);
			auto vec = this->getVectorOfDerivatives( i, point);
			auto buffer = this->getDerivativeOfNByCoordinate_XY(*mat, determinantMat2(*jacoby),vec);
			xy[0]->operator()(i) = (*buffer)(0);
			xy[1]->operator()(i) = (*buffer)(1);
			mat.reset();
		}
		return xy;
	}


	MatUPtr Solver::getMatrixForPoint(MatrixType type, double* X, double*Y, int point)
	{
		auto jacoby = this->getJacobyMatrix2(X, Y, point);
		double detJ = determinantMat2(*jacoby);
		auto mat = std::make_unique<Matrix>(4);
		if (type == MatrixType::H) {
			auto result = getXYDerivativesForPoint(jacoby, point);
			auto H1 = vecAndvecTMultiplication(*result[0]);
			auto H2 = vecAndvecTMultiplication(*result[0]); //?????????????????
			//delete result;
			*mat = H1->operator+= (*H2);
		}
		else if (type == MatrixType::C) {
			auto vec = std::make_shared<Vector>(4);
			for (size_t i = 0; i < 4; i++) {
				vec->operator()(i)=Shape->operator()(i, point);
			}
			mat = std::move(vecAndvecTMultiplication(*vec));// dostan wektor funkcji ksztaltu dla punktu
		}
		*mat *= detJ;
		return mat;
	}

	MatUPtr Solver::getMatrixForElement(MatrixType type,double*X, double*Y,std::vector<double>& multipliers)
	{//matricies, coordinates, multipliers
		auto H = std::make_unique<Matrix>(4);
		const int size = std::pow(this->gaussIntegralScheme, 2);
		const size_t noMultipliers = multipliers.size();

		for (int i = 0; i < size; i++)
		{
			auto buffer = std::move(this->getMatrixForPoint(type, X, Y, i));
			for (size_t j = 0; j < noMultipliers; j++) {
				*buffer *= multipliers.at(j);
			}
			*buffer *= this->wages[i];
			H->operator+=(*buffer);
		}
		return H;
	}

	void Solver::aggregateGlobalMatrix(Matrix& mat, std::vector<MatUPtr>& locals, std::vector<std::array<int, 4>>& nodes)
	{
		const size_t numberOfLocals = locals.size();
		for (size_t i = 0; i < numberOfLocals; i++)
		{
			auto& cur = *locals[i];
			std::array<int, 4>& curNo = nodes[i];
			for (int j = 0; j < 4; j++)
			{
				for (int c = 0; c < 4; c++)
				{
					mat(curNo.at(j) - 1, curNo.at(c) - 1) += cur(j, c);
				}
			}
		}
	}
};