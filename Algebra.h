#pragma once
#include <iostream>
#include <memory>
#include <utility>
#include <vector>


struct Vector;

	class Matrix {
	public:
		Matrix() = default;
		Matrix(int size, int size2 = 0); explicit
		Matrix(const Matrix& mat);
		Matrix(Matrix&& mat);
		Matrix(double* table, const std::pair<int, int> size);

		std::pair<int, int> getSize() const;

		double& operator()(int x, int y);
		double operator()(int x, int y) const;
		Matrix& operator=(const Matrix& mat);
		Matrix& operator=(Matrix&& mat);
		Matrix& operator+=(Matrix& mat);
		Matrix& operator*= (double mul);

		friend std::ostream& operator<< (std::ostream& os, const Matrix& mat);
		friend Matrix& operator/ (Matrix& mat, double scalar);
		friend Matrix& operator+ (Matrix& mat, double scalar);
		friend std::shared_ptr<Matrix> operator+ (Matrix& mat, Matrix& tab); 
		friend std::unique_ptr<Matrix> operator* (Matrix& mat, Matrix& mat1); 
		std::unique_ptr<Vector> multiply (Matrix& mat, Vector& vec);

		void abs();

	protected:
		size_t sizeH, sizeW;
		std::vector<std::vector<double>> tab;
	};

	struct Vector {
		Vector() = default;
		Vector(size_t size);
		Vector(double* data, const size_t size);
		Vector(const Vector& vec);
		Vector(std::unique_ptr<Vector> vec);

		Vector& operator=(Vector& vec);
		Vector& operator=(std::unique_ptr<Vector> vec);
		Vector& operator+=(Vector& mat);
		Vector& operator*=(double scalar);
		
		friend std::unique_ptr<Vector> operator*(Vector& vec, Matrix& mat);
		friend std::unique_ptr<Vector> operator+ (Vector& a, Vector& b);
		friend std::ostream& operator<< (std::ostream& os, const Vector& vec);
		friend std::unique_ptr<Matrix> vecAndvecTMultiplication(Vector& a);
		double& operator()(int);
		double operator()(int) const;

		void  fullFill(double value);

		std::vector<double> tab;
	};

	double determinantMat2(Matrix& mat);
	void transposeMat2(Matrix& mat);
	void inverseMat2(Matrix& mat);
