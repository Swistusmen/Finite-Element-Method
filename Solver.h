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

		MatSPtr getLocalMatrixOfLocalTransformation(LocalType type);

		MatSPtr getJacobyMatrix2(MatSPtr& eta, MatSPtr& ksi, double* x, double* y, int point);
		VecSPtr getVectorOfDerivatives(MatSPtr& ksi, MatSPtr& eta, int fShape, int point); //zwaraca 1 pochodn ksi i eta w wektorze dla odpowiedniej funkcji kszta³tu i pc

		VecSPtr getDerivativeOfNByCoordinate_XY(Matrix& inversedJacoby, double detJ, VecSPtr& derivatives); //derivate of shape funciton and cooridnate
		std::vector<VecSPtr> getXYDerivativesForPoint(MatSPtr& eta, MatSPtr& ksi, MatSPtr& jacoby, int point);
		
		MatUPtr getHSumbatricies(MatSPtr& eta, MatSPtr& ksi, MatSPtr& jacoby, int point);
		MatUPtr getHMatrix(MatSPtr& eta, MatSPtr& ksi, double*X, double*Y, double k);
		
		void aggregateGlobalMatrix(Matrix& mat, std::vector<MatUPtr>& locals, std::vector<std::array<int, 4>>& nodes);
	private:
		int gaussIntegralScheme;
		std::vector<double> wages;
		LocalOperations localOperations;
	};

}
