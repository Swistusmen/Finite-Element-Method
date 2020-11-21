#pragma once
#include <iostream>
#include <memory>
#include <utility>

class Matrix2d {
public:
	Matrix2d();
	Matrix2d(double* tab1);
	Matrix2d(Matrix2d& mat);

	double& operator() (int x, int y);
	friend Matrix2d& operator+ (Matrix2d& mat, Matrix2d& tab);
	Matrix2d& operator= (Matrix2d& mat);
	friend Matrix2d& operator* (Matrix2d& mat, Matrix2d& mat1);
	friend std::ostream& operator<< (std::ostream& os, Matrix2d& mat);
	friend Matrix2d& operator/ (Matrix2d& mat, double scalar);

	Matrix2d& inverse();
	Matrix2d& transpose();
	double determinant();
	
private:
	double **tab;
};

class Matrix4d {
public:
	Matrix4d();
	Matrix4d(double* tab1);
	Matrix4d(Matrix4d& mat);

	double& operator() (int x, int y);
	friend Matrix4d& operator+ (Matrix4d& mat, Matrix4d& tab);
	Matrix4d& operator= (Matrix4d& mat);
	Matrix4d* operator+=(Matrix4d* mat);
	friend Matrix4d& operator* (Matrix4d& mat, Matrix4d& mat1);
	Matrix4d* operator*= ( double mul);
	friend std::ostream& operator<< (std::ostream& os, Matrix4d& mat);
	friend Matrix4d& operator/ (Matrix4d& mat, double scalar);

private:
	double **tab;
};

class MatrixXd
{
public:
	MatrixXd(int size);
	friend std::ostream& operator<< (std::ostream& os, MatrixXd& mat);
	double& operator() (int x, int y);

	int getSize();
private:
	double**tab;
	int size;
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