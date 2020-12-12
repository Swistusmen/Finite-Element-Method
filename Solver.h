#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include <utility>
#include <vector>
#include "Algebra.h"
#include "InputData.h"

namespace slv {
	using VecSPtr = std::shared_ptr<Vector>;
	using VecUPtr = std::unique_ptr<Vector>;
	using MatSPtr = std::shared_ptr<Matrix>;
	using MatUPtr = std::unique_ptr<Matrix>;

	enum LocalType{ETA, KSI,SHAPE};
	enum MatrixType {H , C};

	struct LocalOperations {
		LocalOperations();
		std::vector < std::function<double(double)>> Ksi;
		std::vector < std::function<double(double)>> Eta;
		std::vector < std::function<double(double, double)>> Shape;

		std::vector <std::function<double(double)>>& getFunctions(LocalType type) {
			return type == ETA ? Eta : Ksi;}
		};

	class Solver {
	public:
		Solver(int points);

		MatSPtr getMatrixOfLocalTransformation(LocalType type);

		MatSPtr getJacobyMatrix2(double* x, double* y, int point);
		VecSPtr getVectorOfDerivatives( int fShape, int point); //zwaraca 1 pochodn ksi i eta w wektorze dla odpowiedniej funkcji kszta³tu i pc

		VecSPtr getDerivativeOfNByCoordinate_XY(Matrix& inversedJacoby, double detJ, VecSPtr& derivatives); //derivate of shape funciton and cooridnate
		std::vector<VecSPtr> getXYDerivativesForPoint(MatSPtr& jacoby, int point);
		
		MatUPtr getMatrixForPoint(MatrixType type,  double* X, double*Y, int point);
		MatUPtr getMatrixForElement(MatrixType type, double*X, double*Y, 
			std::vector<double>& multipliers);
		
		void aggregateGlobalMatrix(Matrix& mat, std::vector<MatUPtr>& locals, std::vector<std::array<int, 4>>& nodes);
	private:
		int gaussIntegralScheme;
		std::vector<double> wages;
		LocalOperations localOperations;

		MatSPtr Eta=nullptr;
		MatSPtr Ksi=nullptr;
		MatSPtr Shape=nullptr;

	};

}
