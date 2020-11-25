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
		MatrixXd& getEtaMatrixG2(int size);
		MatrixXd& getKsiMatrixG2(int size);

		Matrix2d& getJacobyMatrix2(MatrixXd& eta, MatrixXd& ksi, double* x, double* y, int point);
		Vector2d& getVectorOfDerivatives(MatrixXd& ksi, MatrixXd& eta, int fShape, int point); //zwaraca 1 pochodn ksi i eta w wektorze dla odpowiedniej funkcji kszta³tu i pc

		Vector2d& getDerivativeOfNByCoordinate_XY(Matrix2d& inversedJacoby, double detJ, Vector2d& derivatives); //derivate of shape funciton and cooridnate
		Vector4d* getXYDerivativesForPoint(MatrixXd& eta, MatrixXd& ksi, Matrix2d& jacoby, int point);
		
		Matrix4d* getHSumbatricies(MatrixXd& eta, MatrixXd& ksi, Matrix2d& jacoby, int point);
		Matrix4d* getHMatrix(MatrixXd& eta, MatrixXd& ksi, Matrix2d& jacoby, double k);
		
		void aggregateGlobalMatrix(MatrixXd& mat, std::vector<Matrix4d*>& locals, std::vector<std::array<int, 4>>& nodes);
	};

}
