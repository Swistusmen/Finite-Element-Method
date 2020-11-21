#include <iostream>
#include "Algebra.h"

//Matrix2d
Matrix2d::Matrix2d()
{
	this->tab = new double*[2];
	tab[0] = new double;
	tab[1] = new double;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			tab[i][j] = 0;
}

Matrix2d::Matrix2d(Matrix2d& mat)
{
	this->tab = new double*[2];
	tab[0] = new double;
	tab[1] = new double;
	for(int i=0;i<2;i++)
		for (int j = 0; j < 2; j++)
		{
			tab[i][j] = mat(i, j);
		}
}

Matrix2d::Matrix2d(double* tab1)
{
	this->tab = new double*[2];
	tab[0] = new double;
	tab[1] = new double;
	int c = 0;
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++,c++)
			tab[i][j] = tab1[c];
}

double& Matrix2d::operator() (int x, int y)
{
	return tab[x][y];
}

Matrix2d& operator+ (Matrix2d& mat, Matrix2d& tab)
{
	Matrix2d* matrix = new Matrix2d();
	for(int i=0;i<2;i++)
		for (int j = 0; j < 2; j++)
		{
			matrix->operator()(i, j) = tab.operator()(i, j) + mat.operator()(i, j);
		}
	return *matrix;
}

Matrix2d& Matrix2d::operator= (Matrix2d& mat)
{
	Matrix2d *matrix = new Matrix2d();
	for (int i = 0; i<2; i++)
		for (int j = 0; j < 2; j++)
		{
			matrix->operator()(i, j) =  mat.operator()(i, j);
		}
	return *matrix;
}

Matrix2d& operator* (Matrix2d& mat, Matrix2d& mat1)
{
	Matrix2d* result = new Matrix2d();
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int c = 0; c < 2; c++) {
				result->operator()(i, j) += mat(i, c)*mat1(c, j);
			}
		}
	}
	return *result;
}

std::ostream& operator<< (std::ostream& os, Matrix2d& mat)
{
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			os << mat(i, j) << " ";
		}
		os << "\n";
	}
	return os;
}

Matrix2d& operator/ (Matrix2d& mat, double scalar)
{
	Matrix2d * matrix = new Matrix2d;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			matrix->operator()(i, j) = mat(i, j) / scalar;
		}
	}
	return *matrix;
}

Matrix2d& Matrix2d::inverse()
{
	Matrix2d* matrix = new Matrix2d;
	double det = this->determinant();
	matrix->operator()(0, 0) = tab[1][1]/det;
	matrix->operator()(1, 1) = tab[0][0]/det;
	matrix->operator()(0, 1) = tab[1][0]/det;
	matrix->operator()(1, 0) = tab[0][1]/det;
	return *matrix;
}

double Matrix2d::determinant()
{
	return tab[0][0] * tab[1][1] - tab[0][1] * tab[1][0];
}

Matrix2d& Matrix2d::transpose()
{
	Matrix2d* matrix = new Matrix2d();
	matrix->operator()(0, 0) = tab[0][0];
	matrix->operator()(1, 1) = tab[1][1];
	matrix->operator()(0, 1) = tab[1][0];
	matrix->operator()(1, 0) = tab[0][1];
	return *matrix;
}

//////////////////////////////////////////////////////////////// Matrix4d
Matrix4d::Matrix4d(Matrix4d& mat)
{
	this->tab = new double*[4];
	tab[0] = new double[4];
	tab[1] = new double[4];
	tab[2] = new double[4];
	tab[3] = new double[4];
	for (int i = 0; i<4; i++)
		for (int j = 0; j < 4; j++)
		{
			tab[i][j] = mat(i, j);
		}
}

Matrix4d::Matrix4d()
{
	this->tab = new double*[4];
	tab[0] = new double[4];
	tab[1] = new double[4];
	tab[2] = new double[4];
	tab[3] = new double[4];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			tab[i][j] = 0;
}

Matrix4d::Matrix4d(double* tab1)
{
	this->tab = new double*[4];
	tab[0] = new double[4];
	tab[1] = new double[4];
	tab[2] = new double[4];
	tab[3] = new double[4];
	int c = 0;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++, c++)
			tab[i][j] = tab1[c];
}

double& Matrix4d::operator() (int x, int y)
{
	return tab[x][y];
}

Matrix4d& operator+ (Matrix4d& mat, Matrix4d& tab)
{
	Matrix4d* matrix = new Matrix4d();
	for (int i = 0; i<4; i++)
		for (int j = 0; j < 4; j++)
		{
			matrix->operator()(i, j) = tab.operator()(i, j) + mat.operator()(i, j);
		}
	return *matrix;
}

Matrix4d& Matrix4d::operator= (Matrix4d& mat)
{
	Matrix4d *matrix = new Matrix4d();
	for (int i = 0; i<4; i++)
		for (int j = 0; j < 4; j++)
		{
			matrix->operator()(i, j) = mat.operator()(i, j);
		}
	return *matrix;
}

Matrix4d& operator* (Matrix4d& mat, Matrix4d& mat1)
{
	Matrix4d* result = new Matrix4d();
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int c = 0; c < 4; c++) {
				result->operator()(i, j) += mat(i, c)*mat1(c, j);
			}
		}
	}
	return *result;
}

std::ostream& operator<< (std::ostream& os, Matrix4d& mat)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			os << mat(i, j) << " ";
		}
		os << "\n";
	}
	return os;
}

Matrix4d& operator/ (Matrix4d& mat, double scalar)
{
	Matrix4d * matrix = new Matrix4d;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrix->operator()(i, j) = mat(i, j) / scalar;
		}
	}
	return *matrix;
}

Matrix4d* Matrix4d::operator+=(Matrix4d* mat)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			this->tab[i][j] += (*mat).operator()(i, j);
		}
	}
	return this;
}

Matrix4d* Matrix4d::operator*= ( double mul)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			tab[i][j] *= mul;
		}
	}

	return this;
}

///////////////////////////////MatrixXd
MatrixXd::MatrixXd(int size)
{
	if (size <= 0)
		throw new std::exception("Dimension of matrix has to be positive number");
	this->size = size;
	this->tab = new double*[this->size];
	for (size_t i = 0; i < size; i++)
	{
		this->tab[i] = new double[size] {0};
	}
}

std::ostream& operator<< (std::ostream& os, MatrixXd& mat)
{
	const size_t size = mat.getSize();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			os << mat(i, j) << " ";
		}
		os << "\n";
	}
	return os;
}

double& MatrixXd::operator() (int x, int y)
{
	if ((x < 0) || (y < 0) || (x >= this->size) || (y >= this->size))
	{
		throw new std::exception("There is no such a cell in this matrix, coordinates are wrong");
	}
	//if (this->tab == nullptr)
	//	throw new std::exception("Uninitialized");
	return this->tab[x][y];
}

int MatrixXd::getSize()
{
	return this->size;
}

/////////////////////////////// Vector2d
Vector2d::Vector2d()
{
	this->tab = new double[2]{ 0 };
}

Vector2d::Vector2d(double* data)
{
	this->tab = new double[2]{ data[0],data[1] };
}

Vector2d::Vector2d(Vector2d& vec)
{
	this->tab = new double[2];
	tab[0] = vec(0);
	tab[1] = vec(1);
}

double& Vector2d::operator()(int i)
{
	return this->tab[i];
}

Vector2d& Vector2d::operator=(Vector2d& vec)
{
	Vector2d* tab = new Vector2d();
	tab->operator()(0) = vec(0);
	tab->operator()(1) = vec(1);
	return *tab;
}

Vector2d& operator+ (Vector2d& a, Vector2d& b)
{
	Vector2d * vec = new Vector2d();
	vec->operator()(0) = a(0) + b(0);
	vec->operator()(1) = a(1) + b(1);
	return *vec;
}

std::ostream& operator <<(std::ostream& os, Vector2d& vec)
{
	os << vec(0) << " " << vec(1);
	return os;
}

Matrix2d& Vector2d::vecAndvecTMultiplication(Vector2d& a)
{
	Matrix2d *mat = new Matrix2d();
	mat->operator()(0, 0) = a(0)*a(0);
	mat->operator()(0, 1) = a(0)*a(1);
	mat->operator()(1, 0) = a(1)*a(0);
	mat->operator()(1, 1) = a(0)*a(0);
	return *mat;
}

/////////////////////////////// Vector4d
Vector4d::Vector4d()
{
	this->tab = new double[4]{ 0 };
}

Vector4d::Vector4d(double* data)
{
	this->tab = new double[4]{ data[0],data[1],data[2],data[3]};
}

Vector4d::Vector4d(Vector4d& vec)
{
	this->tab = new double[4];
	tab[0] = vec(0);
	tab[1] = vec(1);
	tab[2] = vec(2);
	tab[3] = vec(3);
}

double& Vector4d::operator()(int i)
{
	return this->tab[i];
}

Vector4d& Vector4d::operator=( Vector4d& vec)
{
	Vector4d* tab = new Vector4d();
	tab->operator()(0) = vec(0);
	tab->operator()(1) = vec(1);
	tab->operator()(2) = vec(2);
	tab->operator()(3) = vec(3);
	return *tab;
}

Vector4d& operator+ (Vector4d& a, Vector4d& b)
{
	Vector4d * vec = new Vector4d();
	vec->operator()(0) = a(0) + b(0);
	vec->operator()(1) = a(1) + b(1);
	vec->operator()(2) = a(2) + b(2);
	vec->operator()(3) = a(3) + b(3);
	return *vec;
}

std::ostream& operator <<(std::ostream& os, Vector4d& vec)
{
	os << vec(0) << " " << vec(1)<<" "<<vec(2)<<" "<<vec(3);
	return os;
}

Matrix4d& vecAndvecTMultiplication(Vector4d& a)
{
	Matrix4d *mat = new Matrix4d();
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat->operator()(i, j) = a(i)*a(j);
		}
	}
	return *mat;
}