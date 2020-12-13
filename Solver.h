#pragma once

#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include "Algebra.h"
#include "InputData.h"
#include "NumericalStaticData.h"

namespace slv {
	using VecSPtr = std::shared_ptr<Vector>;
	using VecUPtr = std::unique_ptr<Vector>;
	using MatSPtr = std::shared_ptr<Matrix>;
	using MatUPtr = std::unique_ptr<Matrix>;

	class Solver {
	public:
		Solver(LocalOperations& localOperations, IntegrationalPoints& iPoints,int scheme, data::RectangleMeshInput& data);

		MatSPtr getMatrixOfLocalTransformation(LocalType type);

		MatSPtr getJacobyMatrix2(double* x, double* y, int point);
		VecSPtr getVectorOfDerivatives( int fShape, int point);

		VecSPtr getDerivativeOfNByCoordinate_XY(Matrix& inversedJacoby, double detJ, VecSPtr& derivatives);
		std::vector<VecSPtr> getXYDerivativesForPoint(MatSPtr& jacoby, int point);
		
		MatUPtr getMatrixForPoint(MatrixType type,  double* X, double*Y, int point);
		MatUPtr getMatrixForElement(MatrixType type, double*X, double*Y, 
			std::vector<double>& multipliers);
		MatUPtr getBoundaryMatrixForElement(double* X, double* Y, std::vector<double>& multipliers);
		
		void aggregateGlobalMatrix(Matrix& mat, std::vector<MatUPtr>& locals, std::vector<std::array<int, 4>>& nodes);
	private:
		int gaussIntegralScheme;
		LocalOperations localOperations;

		std::vector<double> wagesOneDim;
		std::vector<double> wagesTwoDim;
		IntegrationalPoints iPoints;

		double HBound = -1.0;
		double WBound = -1.0;

		MatSPtr Eta=nullptr;
		MatSPtr Ksi=nullptr;
		MatSPtr Shape=nullptr;
	};

}
