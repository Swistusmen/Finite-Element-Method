#include "Solver.h"
#include <memory>
#include <Eigen/Dense>
#include <utility>
#include <cmath>
#include <vector>
#include <array>
#include "InputData.h"
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
			tempData[0] = 0.55;
			tempData[1] = 0.88;
			tempData[2] = 0.55;
		}
		else {
			tempData[0] = 0.347;
			tempData[1] = 0.652;
			tempData[2] = 0.652;
			tempData[3] = 0.347;
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

	MatSPtr Solver::getEtaMatrix()
	{
		int size = std::pow(this->gaussIntegralScheme, 2);
		auto matrix = std::make_shared<Matrix>(4, size);
		if (size == 4)
		{
			double ksi = 1 / sqrt(3);
			for (int i = 0; i < 4; i++)
			{
				ksi *= -1;
				matrix->operator()(i, 0) = (ksi - 1) / 4;
				matrix->operator()(i, 1) = (-ksi - 1) / 4;
				matrix->operator()(i, 2) = (ksi + 1) / 4;
				matrix->operator()(i, 3) = (-ksi + 1) / 4;
			}
			return matrix;
		}
		else if (size == 9)
		{
			double values[12] = { -0.04365,-0.25,-0.05635 ,-0.05635,-0.25,-0.44365,0.05635,0.25,0.44365,0.44365,0.25,0.05635 };
			for (int i = 0, c = 0; i < 4; i++, c++)
			{
				matrix->operator()(i, 0) = values[c * 3];
				matrix->operator()(i, 1) = matrix->operator()(i, 0);
				matrix->operator()(i, 2) = matrix->operator()(i, 0);
				matrix->operator()(i, 3) = values[c * 3 + 1];
				matrix->operator()(i, 4) = matrix->operator()(i, 3);
				matrix->operator()(i, 5) = matrix->operator()(i, 3);
				matrix->operator()(i, 6) = values[c * 3 + 2];
				matrix->operator()(i, 7) = matrix->operator()(i, 6);
				matrix->operator()(i, 8) = matrix->operator()(i, 6);
			}
			return matrix;
		}
		else if (size == 16) //need to be checked if is correct
		{
			double values[16] = { -0.46528,-0.335,-0.165,-0.3472,0.465284,0.334995,0.165005,0.034716,0.034716 ,0.165005,0.334995,0.465284,-0.03472,-0.165,-0.335,-0.46528 };
			for (int i = 0, c = 0; i < 4; i++, c++)
			{
				matrix->operator()(i, 0) = values[c * 4];
				matrix->operator()(i, 1) = values[c * 4];
				matrix->operator()(i, 2) = values[c * 4];
				matrix->operator()(i, 3) = values[c * 4];
				matrix->operator()(i, 4) = values[c * 4 + 1];
				matrix->operator()(i, 5) = values[c * 4 + 1];
				matrix->operator()(i, 6) = values[c * 4 + 1];
				matrix->operator()(i, 7) = values[c * 4 + 1];
				matrix->operator()(i, 8) = values[c * 4 + 2];
				matrix->operator()(i, 9) = values[c * 4 + 2];
				matrix->operator()(i, 10) = values[c * 4 + 2];
				matrix->operator()(i, 11) = values[c * 4 + 2];
				matrix->operator()(i, 12) = values[c * 4 + 3];
				matrix->operator()(i, 13) = values[c * 4 + 3];
				matrix->operator()(i, 14) = values[c * 4 + 3];
				matrix->operator()(i, 15) = values[c * 4 + 3];
			}
			return matrix;
		}
	}

	MatSPtr Solver::getKsiMatrix()
	{
		int size = std::pow(this->gaussIntegralScheme, 2);
		auto matrix = std::make_shared<Matrix>(4, size);
		if (size == 4)
		{
			double eta = 1 / sqrt(3);
			for (int i = 0; i < 4; i++)
			{
				eta *= -1;
				matrix->operator()(i, 0) = (eta - 1) / 4;
				matrix->operator()(i, 1) = (-eta + 1) / 4;
				matrix->operator()(i, 2) = (eta + 1) / 4;
				matrix->operator()(i, 3) = (-eta - 1) / 4;
			}
			return matrix;
		}
		else if (size == 9)
		{
			double values[12] = { -0.04365,-0.25,-0.05635 ,0.04365,0.25,0.05635,0.05635,0.25,0.44365,-0.05635,-0.25,-0.44365 };
			for (int i = 0, c = 0; i < 4; i++, c++)
			{
				matrix->operator()(i, 0) = values[c * 3 + 0];
				matrix->operator()(i, 1) = values[c * 3 + 1];
				matrix->operator()(i, 2) = values[c * 3 + 2];
				matrix->operator()(i, 3) = matrix->operator()(i, 0);
				matrix->operator()(i, 4) = matrix->operator()(i, 1);
				matrix->operator()(i, 5) = matrix->operator()(i, 2);
				matrix->operator()(i, 6) = matrix->operator()(i, 0);
				matrix->operator()(i, 7) = matrix->operator()(i, 1);
				matrix->operator()(i, 8) = matrix->operator()(i, 2);
			}
			return matrix;
		}
		else if (size == 16)
		{
			double values[16] = { -0.46528,-0.335,-0.165,-0.3472,0.465284,0.334995,0.165005,0.034716,0.034716 ,0.165005,0.334995,0.465284,-0.03472,-0.165,-0.335,-0.46528 };
			for (int i = 0, c = 0; i < 4; i++, c++)
			{
				matrix->operator()(i, 0) = values[c * 4];
				matrix->operator()(i, 1) = values[c * 4 + 1];
				matrix->operator()(i, 2) = values[c * 4 + 2];
				matrix->operator()(i, 3) = values[c * 4 + 3];
				matrix->operator()(i, 4) = values[c * 4];
				matrix->operator()(i, 5) = values[c * 4 + 1];
				matrix->operator()(i, 6) = values[c * 4 + 2];
				matrix->operator()(i, 7) = values[c * 4 + 3];
				matrix->operator()(i, 8) = values[c * 4];
				matrix->operator()(i, 9) = values[c * 4 + 1];
				matrix->operator()(i, 10) = values[c * 4 + 2];
				matrix->operator()(i, 11) = values[c * 4 + 3];
				matrix->operator()(i, 12) = values[c * 4];
				matrix->operator()(i, 13) = values[c * 4 + 1];
				matrix->operator()(i, 14) = values[c * 4 + 2];
				matrix->operator()(i, 15) = values[c * 4 + 3];
			}
			return matrix;
		}
	}

	MatSPtr Solver::getJacobyMatrix2(MatSPtr& eta, MatSPtr& ksi, double* x, double* y, int point)
	{
		auto matrix = std::make_shared<Matrix>(2);
		const int size = std::pow(this->gaussIntegralScheme, 2);
		for (int i = 0; i < size; i++)
			matrix->operator()(0, 0) += x[i] * (*ksi).operator()(point, i);
		for (int i = 0; i < size; i++)
			matrix->operator()(0, 1) += y[i] * (*ksi).operator()(point, i);
		for (int i = 0; i < size; i++)
			matrix->operator()(1, 0) += x[i] * (*eta).operator()(point, i);
		for (int i = 0; i < size; i++)
			matrix->operator()(1, 1) += y[i] * (*eta).operator()(point, i);

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
		//vec->operator()(0) = ksi(point, fShape);
		//vec->operator()(1) = eta(point, fShape);
		vec->operator()(0) = ksi->operator()(fShape, point);
		vec->operator()(1) = eta->operator()(fShape, point); //check!!!
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
		auto H2 = vecAndvecTMultiplication(*result[0]); //?????????????????????????????????????????
		//delete result;
		auto mat = std::make_unique<Matrix>(4);
		*mat = H1->operator+= (*H2);
		return mat;
	}

	MatUPtr Solver::getHMatrix(MatSPtr& eta, MatSPtr& ksi, MatSPtr& jacoby, double k)
	{
		double detJ = determinantMat2(*jacoby);
		auto H = std::make_unique<Matrix>(4);
		const int size = std::pow(this->gaussIntegralScheme, 2);
		for (int i = 0; i < size; i++)
		{
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