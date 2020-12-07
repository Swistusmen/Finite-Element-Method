#pragma once
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

class Matrix {
public:
	Matrix()=default;
	Matrix(int size,int size2=0);
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
	friend std::shared_ptr<Matrix> operator+ (Matrix& mat, Matrix& tab); //TODO
	friend std::shared_ptr<Matrix> operator* (Matrix& mat, Matrix& mat1); //TODO

protected:
	unsigned int sizeH, sizeW;
	std::vector<std::vector<double>> tab;
};

class Matrix2d:public Matrix {
public:
	Matrix2d() :Matrix(2) {};
	Matrix2d(double* tab1) :Matrix(tab1, std::make_pair(2, 2)) {};
	Matrix2d(const Matrix& mat) :Matrix(mat) {};

	friend Matrix2d& operator+ (Matrix2d& mat, Matrix2d& tab);
	friend Matrix2d& operator* (Matrix2d& mat, Matrix2d& mat1);

	Matrix2d& inverse();
	Matrix2d& transpose();
	double determinant();
};

class Matrix4d : public Matrix {
public:
	Matrix4d() :Matrix(4) {};
	Matrix4d(double* tab1) :Matrix(tab1, std::make_pair(4, 4)) {};
	Matrix4d(const Matrix4d& mat) :Matrix(mat) {};

	friend Matrix4d& operator+ (Matrix4d& mat, Matrix4d& tab);
	Matrix4d* operator+=(Matrix4d* mat);
	friend Matrix4d& operator* (Matrix4d& mat, Matrix4d& mat1);
	Matrix4d* operator*= ( double mul);
};

class Vector2d {
public:
	Vector2d();
	Vector2d(double* data);
	Vector2d(Vector2d& vec);

	Vector2d& operator=(Vector2d& vec);
	friend Vector2d& operator+ (Vector2d& a, Vector2d& b);
	friend std::ostream& operator <<(std::ostream& os, Vector2d& vec);
	Matrix2d& vecAndvecTMultiplication(Vector2d& a);
	double& operator()(int);
private:
	double* tab;
};

class Vector4d {
public:
	Vector4d();
	Vector4d(double* data);
	Vector4d(Vector4d& vec);

	Vector4d& operator=(Vector4d& vec);
	friend Vector4d& operator+ (Vector4d& a, Vector4d& b);
	friend std::ostream& operator<< (std::ostream& os, Vector4d& vec);
	friend Matrix4d& vecAndvecTMultiplication(Vector4d& a);
	double& operator()(int);

private:
	double* tab;
};

//Matrix* multiplyVecAndTransposedVec(std::vector<double>& vec1, std::vector<double>& vec2);
