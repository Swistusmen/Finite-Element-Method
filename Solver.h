#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include <utility>
#include <vector>
#include "Algebra.h"

namespace slv {
	class Solver {
	public:
		Solver(int points);

		Matrix& getEtaMatrix();
		Matrix& getKsiMatrix();

		Matrix2d& getJacobyMatrix2(Matrix& eta, Matrix& ksi, double* x, double* y, int point);
		Vector2d& getVectorOfDerivatives(Matrix& ksi, Matrix& eta, int fShape, int point); //zwaraca 1 pochodn ksi i eta w wektorze dla odpowiedniej funkcji kszta³tu i pc

		Vector2d& getDerivativeOfNByCoordinate_XY(Matrix2d& inversedJacoby, double detJ, Vector2d& derivatives); //derivate of shape funciton and cooridnate
		Vector4d* getXYDerivativesForPoint(Matrix& eta, Matrix& ksi, Matrix2d& jacoby, int point);
		
		Matrix4d* getHSumbatricies(Matrix& eta, Matrix& ksi, Matrix2d& jacoby, int point);
		Matrix4d* getHMatrix(Matrix& eta, Matrix& ksi, Matrix2d& jacoby, double k);
		
		void aggregateGlobalMatrix(Matrix& mat, std::vector<Matrix4d*>& locals, std::vector<std::array<int, 4>>& nodes);
	private:
		int gaussIntegralScheme;
		std::vector<double> wages;
	};

}
