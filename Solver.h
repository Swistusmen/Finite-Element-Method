#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include <utility>
#include <vector>
#include "Algebra.h"
#include "InputData.h"

namespace slv {
	class Solver {
	public:
		Solver(int points);

		Matrix& getEtaMatrix();
		Matrix& getKsiMatrix();
		Matrix* getNMatrix();

		Matrix& getJacobyMatrix2(Matrix& eta, Matrix& ksi, double* x, double* y, int point);
		Vector& getVectorOfDerivatives(Matrix& ksi, Matrix& eta, int fShape, int point); //zwaraca 1 pochodn ksi i eta w wektorze dla odpowiedniej funkcji kszta³tu i pc

		Vector& getDerivativeOfNByCoordinate_XY(Matrix& inversedJacoby, double detJ, Vector& derivatives); //derivate of shape funciton and cooridnate
		Vector* getXYDerivativesForPoint(Matrix& eta, Matrix& ksi, Matrix& jacoby, int point);
		
		Matrix* getHSumbatricies(Matrix& eta, Matrix& ksi, Matrix& jacoby, int point);
		Matrix* getHMatrix(Matrix& eta, Matrix& ksi, Matrix& jacoby, double k);

		Matrix* getCLocalMatrix(double Eta, double Ksi);
		Matrix* getCMatrix(data::Elem4& data, Matrix& jacoby, double ro, double temp);
		
		void aggregateGlobalMatrix(Matrix& mat, std::vector<Matrix*>& locals, std::vector<std::array<int, 4>>& nodes);
	private:
		int gaussIntegralScheme;
		std::vector<double> wages;
	};

}
