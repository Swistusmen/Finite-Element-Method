#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include <utility>
#include <vector>
#include "Algebra.h"

namespace slv {
	using mat4U = std::unique_ptr<Matrix4d>;
	using mat4S = std::shared_ptr<Matrix4d>;
	using mat2U = std::unique_ptr<Matrix2d>;
	using mat2S = std::shared_ptr<Matrix2d>;
	using vec4U = std::unique_ptr<Vector4d>;
	using vec4S = std::shared_ptr<Vector4d>;
	using vec2U = std::unique_ptr<Vector2d>;
	using vec2S = std::shared_ptr<Vector2d>;

	class Solver {
	public:
		Matrix4d& getEtaMatrixG2();
		Matrix4d& getKsiMatrixG2();

		Matrix2d& getJacobyMatrix2(Matrix4d& eta, Matrix4d& ksi, double* x, double* y, int point);
		Vector2d& getVectorOfDerivatives(Matrix4d& ksi, Matrix4d& eta, int fShape, int point); //zwaraca 1 pochodn ksi i eta w wektorze dla odpowiedniej funkcji kszta³tu i pc

		Vector2d& getDerivativeOfNByCoordinate_XY(Matrix2d& inversedJacoby, double detJ, Vector2d& derivatives); //derivate of shape funciton and cooridnate
		Vector4d* getXYDerivativesForPoint(Matrix4d& eta, Matrix4d& ksi, Matrix2d& jacoby, int point);
		
		Matrix4d* getHSumbatricies(Matrix4d& eta, Matrix4d& ksi, Matrix2d& jacoby, int point);
		Matrix4d* getHMatrix(Matrix4d& eta, Matrix4d& ksi, Matrix2d& jacoby, double k);
	};

}
