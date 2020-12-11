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
	}

	MatSPtr Solver::getLocalMatrixOfLocalDerivatives(LocalType type)
	{
		int noIntPoints = this->gaussIntegralScheme;
		size_t size = std::pow(this->gaussIntegralScheme, 2);
		auto matrix = std::make_shared<Matrix>(4, size);
		std::vector<double> points= std::move(generateGaussIntegrationalPoints(noIntPoints)) ;
		switch (type) {
		case LocalType::ETA: {
			for (size_t i = 0; i < size; i++)
			{
				matrix->operator()(0, i) = -0.25*(1 - points.at(i % noIntPoints));
				matrix->operator()(1, i) = -0.25*(1 + points.at(i % noIntPoints));
				matrix->operator()(2, i) = 0.25*(points.at(i % noIntPoints) + 1);
				matrix->operator()(3, i) = 0.25*(-points.at(i % noIntPoints) + 1);
			}
			}break;
		case LocalType::KSI: {
			for (size_t i = 0; i < size; i++)
			{
				matrix->operator()(0, i) = -0.25*(1 - points.at(i % noIntPoints));
				matrix->operator()(1, i) = 0.25*(-points.at(i % noIntPoints) + 1);
				matrix->operator()(2, i) = 0.25*(points.at(i % noIntPoints) + 1);
				matrix->operator()(3, i) = -0.25*(points.at(i % noIntPoints) + 1);
			}
		}
		}
		return matrix;
	}

	VecSPtr Solver::getDerivativeOfNByCoordinate_XY(Matrix& inversedJacoby, double detJ, VecSPtr& derivatives)
	{
		auto vec = std::make_shared<Vector>(2);
		vec->operator()(0) = (inversedJacoby(0, 0)*derivatives->operator()(0) + inversedJacoby(0, 1)*derivatives->operator()(1)); ///zgodnie ze wzorem powinno byc dzielenie przez detJ
		vec->operator()(1) = (inversedJacoby(1, 0)*derivatives->operator()(0) + inversedJacoby(1, 1)*derivatives->operator()(1));
		return vec;
	}

	VecSPtr Solver::getVectorOfDerivatives(MatSPtr& ksi, MatSPtr& eta, int fShape, int point)
	{
		auto vec = std::make_shared<Vector>(2);
		vec->operator()(0) = ksi->operator()(fShape, point);
		vec->operator()(1) = eta->operator()(fShape, point);
		return vec;
	}

	std::vector<VecSPtr> Solver::getXYDerivativesForPoint(MatSPtr& eta, MatSPtr& ksi, MatSPtr& jacoby, int point)
	{
		std::vector<VecSPtr> xy;
		xy.push_back(std::make_shared<Vector>(4));
		xy.push_back(std::make_shared<Vector>(4));
		const size_t size = std::pow(this->gaussIntegralScheme, 2);
		for (int i = 0; i < 4; i++)
		{
			auto mat=std::make_unique<Matrix>(*jacoby);
			inverseMat2(*mat);
			auto vec = this->getVectorOfDerivatives(ksi, eta, i, point);
			auto buffer = this->getDerivativeOfNByCoordinate_XY(*mat, determinantMat2(*jacoby),vec);
			xy[0]->operator()(i) = (*buffer)(0);
			xy[1]->operator()(i) = (*buffer)(1);
			mat.reset();
		}
		return xy;
	}


	MatUPtr Solver::getHSumbatricies(MatSPtr& eta, MatSPtr& ksi, MatSPtr& jacoby, int point)
	{
		auto result = getXYDerivativesForPoint(eta, ksi, jacoby, point);
		auto H1 = vecAndvecTMultiplication(*result[0]);
		auto H2 = vecAndvecTMultiplication(*result[0]); //?????????????????
		//delete result;
		auto mat = std::make_unique<Matrix>(4);
		*mat = H1->operator+= (*H2);
		return mat;
	}

	MatUPtr Solver::getHMatrix(MatSPtr& eta, MatSPtr& ksi, double* X, double*Y,double k)
	{
		auto H = std::make_unique<Matrix>(4);
		const int size = std::pow(this->gaussIntegralScheme, 2);

		for (int i = 0; i < size; i++)
		{
			auto jacoby = this->getJacobyMatrix2(eta, ksi, X, Y, i);
			double detJ = determinantMat2(*jacoby);
			auto buffer = std::move(this->getHSumbatricies(eta, ksi, jacoby, i));
			*buffer *= k * detJ;
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
	/*
	MatUPtr Solver::getCLocalMatrix(double Eta, double Ksi)
	{
		Vector* vec = new Vector(4);
		vec->operator()(0) = 0.25*(1 - Eta)*(1 - Ksi);
		vec->operator()(1) = 0.25*(1 - Eta)*(1 + Ksi);
		vec->operator()(2) = 0.25*(1 + Eta)*(1 + Ksi);
		vec->operator()(3) = 0.25*(1 + Eta)*(1 - Ksi);
		auto mat = std::make_unique<Matrix>( std::move(vecAndvecTMultiplication(*vec)));
		return mat;
	}

	MatUPtr Solver::getCMatrix(data::Elem4& data, MatSPtr& jacoby, double ro, double temp)
	{
		size_t size = this->gaussIntegralScheme;
		double det = determinantMat2(*jacoby);

		auto mat = std::make_unique<Matrix>();
		for (size_t i = 0; i < size; ++i)
		{
			for (size_t j = 0; j < size; j++)
			{
				auto buffer = this->getCLocalMatrix(data.tabEta[i], data.tabKsi[j]);
				*buffer *= ro;
				*buffer *= temp;
				*buffer *= det;
				mat->operator +=(*buffer);
			}
		}
		return mat;
	}
	*/
};