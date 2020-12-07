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

	class Solver {
	public:
		Solver(int points);

		MatSPtr getEtaMatrix();
		MatSPtr getKsiMatrix();
		Matrix* getNMatrix();

		MatSPtr getJacobyMatrix2(MatSPtr& eta, MatSPtr& ksi, double* x, double* y, int point);
		VecSPtr getVectorOfDerivatives(MatSPtr& ksi, MatSPtr& eta, int fShape, int point); //zwaraca 1 pochodn ksi i eta w wektorze dla odpowiedniej funkcji kszta³tu i pc

		VecSPtr getDerivativeOfNByCoordinate_XY(Matrix& inversedJacoby, double detJ, VecSPtr& derivatives); //derivate of shape funciton and cooridnate
		std::vector<VecSPtr> getXYDerivativesForPoint(MatSPtr& eta, MatSPtr& ksi, MatSPtr& jacoby, int point);
		
		MatUPtr getHSumbatricies(MatSPtr& eta, MatSPtr& ksi, MatSPtr& jacoby, int point);
		MatUPtr getHMatrix(MatSPtr& eta, MatSPtr& ksi, MatSPtr& jacoby, double k);
		/*
		MatUPtr getCLocalMatrix(double Eta, double Ksi);
		MatUPtr getCMatrix(data::Elem4& data, MatSPtr& jacoby, double ro, double temp);
		*/
		void aggregateGlobalMatrix(Matrix& mat, std::vector<MatUPtr>& locals, std::vector<std::array<int, 4>>& nodes);
	private:
		int gaussIntegralScheme;
		std::vector<double> wages;
	};

}
