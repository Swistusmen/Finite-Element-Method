#include "Solver.h"
#include <memory>
#include <Eigen/Dense>
#include <utility>
#include <cmath>
#include <vector>
#include <array>
#include "InputData.h"
#include "AlgebraicOperations.h"

#define Epsilon 0.000001

namespace {
	std::vector<double> getWagesOneDimension(int scheme)
	{
		std::vector<double> wagesOneDimension;
		if (scheme == 2) {
			wagesOneDimension.push_back(1);
			wagesOneDimension.push_back(1);
		}
		else if (scheme == 3) {
			wagesOneDimension.push_back(5.0 / 9.0);
			wagesOneDimension.push_back(8.0 / 9.0);
			wagesOneDimension.push_back(5.0 / 9.0);
		}
		else {
			wagesOneDimension.push_back((18.0 - std::sqrt(30.0)) / 36.0);
			wagesOneDimension.push_back((18.0 + std::sqrt(30.0)) / 36.0);
			wagesOneDimension.push_back((18.0 + std::sqrt(30.0)) / 36.0);
			wagesOneDimension.push_back((18.0 - std::sqrt(30.0)) / 36.0);
		}
		return wagesOneDimension;
	}

	std::vector<double> getWagesTwoDimensions(int scheme)
	{
		std::vector<double> wagesOne = std::move(getWagesOneDimension(scheme));
		std::vector<double> wagesTwoDimension;
			for (size_t i = 0; i < scheme; ++i)
			{
				for (size_t j = 0; j < scheme; ++j)
				{
					wagesTwoDimension.push_back(wagesOne[i] * wagesOne[j]);
				}
			}
			return wagesTwoDimension;
	}
}

namespace slv {
	Solver::Solver(LocalOperations& localOperations, IntegrationalPoints& iPoints,int scheme,data::RectangleMeshInput& data)
	{
		this->gaussIntegralScheme = scheme;
		this->wagesOneDim = getWagesOneDimension(scheme);
		this->wagesTwoDim = getWagesTwoDimensions(scheme);
		this->iPoints = iPoints;
		this->localOperations = localOperations;

		//Resolution only for rectangle
		this->HBound = data.H;
		this->WBound = data.W;

		this->Eta = (getMatrixOfLocalTransformation(LocalType::ETA));
		this->Ksi = (getMatrixOfLocalTransformation(LocalType::KSI));
		this->Shape = (getMatrixOfLocalTransformation(LocalType::SHAPE2D));

		std::cout <<"Eta matrix\n"<< *Eta << std::endl;
		std::cout << "Ksi matrix\n" << *Ksi << std::endl;
	}
	
	MatSPtr Solver::getMatrixOfLocalTransformation(LocalType type)
	{
		int noIntPoints = this->gaussIntegralScheme;
		size_t size = std::pow(this->gaussIntegralScheme, 2);
		auto matrix = std::make_shared<Matrix>(4, size);
		std::vector<double> points = std::move(this->iPoints.getDependingonScheme(this->gaussIntegralScheme));
		switch (type) {
		case LocalType::ETA: case LocalType::KSI: {
			auto functions = localOperations.getFunctions(type);
			for (size_t i = 0; i < size; i++)
			{
				for (size_t j = 0; j < 4; j++)
					matrix->operator()(j, i) = functions.at(j)(points.at(i % noIntPoints));
			}
		}break;
		case LocalType::SHAPE2D: {
			for (size_t i = 0; i < size; i++)
			{
				for (size_t j = 0; j < 4; j++)
					matrix->operator()(j, i) = (localOperations.Shape2D.at(j))(points.at(i % noIntPoints),
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
			gjInverseMatrix(mat);
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
			auto H2 = vecAndvecTMultiplication(*result[1]); 
			*mat = H1->operator+= (*H2);
		}
		else if (type == MatrixType::C) {
			auto vec = std::make_shared<Vector>(4);
			for (size_t i = 0; i < 4; i++) {
				vec->operator()(i)=Shape->operator()(i, point);
			}
			mat = std::move(vecAndvecTMultiplication(*vec));
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
			*buffer *= this->wagesTwoDim[i];
			H->operator+=(*buffer);
		}
		return H;
	}

	MatUPtr Solver::getBoundaryMatrixForSide(std::array<int, 4> config, double det, std::vector<double>& multipliers)
	{
		const size_t noOneDShapeFunctions = this->gaussIntegralScheme;
		const size_t noMultipliers = multipliers.size();
		auto Hbc= std::make_unique<Matrix>(4);
		auto points = iPoints.getDependingonScheme(noOneDShapeFunctions);
		auto wages = getWagesOneDimension(noOneDShapeFunctions);
		for (size_t i = 0; i < noOneDShapeFunctions; i++)
		{
			auto vec = std::make_shared<Vector>(4);
			vec->operator()(config.at(0)) +=localOperations.Shape1D.at(config.at(2))(points.at(i));
			vec->operator()(config.at(1)) +=localOperations.Shape1D.at(config.at(3))(points.at(i));
			auto temp = vecAndvecTMultiplication(*vec);
			temp->operator*=(wages[i]);
			Hbc->operator+=(*temp);
		}
		for (size_t j = 0; j < noMultipliers; j++)
		{
			Hbc->operator*=  (multipliers.at(j));
		}
		Hbc->operator*=(det);
		return Hbc;
	}

	MatUPtr Solver::getBoundaryMatrixForElement(double* X, double* Y, std::vector<double>& multipliers)
	{
		std::vector<MatUPtr> Hbc;
		double detW = (X[1] - X[0]) / 2.0;
		double detH = (Y[3] - Y[0]) / 2.0; //tu
		if (X[0] == 0.0)
			Hbc.push_back(std::move(this->getBoundaryMatrixForSide(bCondition.Left, detW, multipliers)));
		if (std::abs(X[1] - this->WBound) < Epsilon)
			Hbc.push_back(std::move(this->getBoundaryMatrixForSide(bCondition.Right, detW, multipliers)));
		if (Y[0] == 0.0)
			Hbc.push_back(std::move(this->getBoundaryMatrixForSide(bCondition.Bottom, detH, multipliers)));
		if (std::abs(Y[3] - HBound) < Epsilon)
			Hbc.push_back(std::move(this->getBoundaryMatrixForSide(bCondition.Top, detH, multipliers)));
		const size_t size = Hbc.size();
		auto result = std::make_unique<Matrix>(4);
		for (size_t i = 0; i < size; i++)
		{
			result->operator += (*Hbc[i]);
		}
		return result;
	}

	VecUPtr Solver::getPreassureVectorForSide(std::array<int, 4> config, double det, std::vector<double>& multipliers)
	{
		const size_t noOneDShapeFunctions = this->gaussIntegralScheme;
		const size_t noMultipliers = multipliers.size();
		auto Vec = std::make_unique<Vector>(4);
		auto points = iPoints.getDependingonScheme(noOneDShapeFunctions);
		for (size_t i = 0; i < noOneDShapeFunctions; i++)
		{
			Vec->operator()(config.at(0)) += localOperations.Shape1D.at(config.at(2))(points.at(i));
			Vec->operator()(config.at(1)) += localOperations.Shape1D.at(config.at(3))(points.at(i));
		}
		for (size_t j = 0; j < noMultipliers; j++)
		{
			Vec->operator*=(multipliers.at(j));
		}
		return Vec;
	}

	VecUPtr Solver::getPreassureVectorForElement(double* X, double* Y, std::vector<double>& multipliers)
	{
		auto Hbc = std::make_unique<Vector>(4);
		double detW = ((X[1] - X[0]) / 2.0);//check correctness of this
		double detH = ((Y[3] - Y[0]) / 2.0);

		if (X[0] == 0.0) {
			auto temp = *(getPreassureVectorForSide(bCondition.Left, detW, multipliers));
			temp *= detW;
			Hbc->operator+=(temp);
		}
		if (std::abs(X[1] - this->WBound) < Epsilon) {
			auto temp = *(getPreassureVectorForSide(bCondition.Right, detW, multipliers));
			temp *= detW;
			Hbc->operator+=(temp);
		}
		if (Y[0] == 0.0) {
			auto temp = *(getPreassureVectorForSide(bCondition.Bottom, detH, multipliers));
			temp *= detH;
			Hbc->operator+=(temp);
		}
		if (std::abs(Y[3] - HBound) < Epsilon) {
			auto temp = *(getPreassureVectorForSide(bCondition.Top, detH, multipliers));
			temp *= detH;
			Hbc->operator+=(temp);
		}
		//Hbc->operator*=(-1.0);//it's from the formula specific to the fourier equation (q=alfa(t-t*alfa))
		return Hbc;
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

	void Solver::aggregateGlobalVector(Vector& vec, std::vector<VecUPtr>& locals, std::vector<std::array<int, 4>>& nodes)
	{
		const size_t numberOfLocals = locals.size();
		for (size_t i = 0; i < numberOfLocals; i++)
		{
			auto& cur = *locals[i];
			std::array<int, 4>& curNo = nodes[i];
			for (int j = 0; j < 4; j++)
			{
					vec(curNo.at(j) - 1) += cur(j);
			}
		}
	}
};