#include <iostream>
#include "Algebra.h"

Matrix::Matrix(int size, int size2)
{
	if (size <= 0 || size2 < 0)
		throw new std::exception("Matrix constructor error");
	if (size2 == 0)
		size2 = size;
	this->sizeH = size;
	this->sizeW = size2 == 0 ? size : size2;
	this->tab = new double*[this->sizeH];
	for (int i = 0; i < this->sizeH; i++)
	{
		this->tab[i] = new double[this->sizeW];
		for (int j = 0; j < this->sizeW; j++)
		{
			this->tab[i][j] = 0.0;
		}
	}
}

Matrix::Matrix(Matrix& mat)
{
	auto size = mat.getSize();
	this->sizeH = size.first;
	this->sizeW = size.second;
	this->tab = new double*[this->sizeH];
	for (int i = 0; i < this->sizeH; i++)
	{
		this->tab[i] = new double[this->sizeW];
		for (int j = 0; j < this->sizeW; j++)
		{
			this->tab[i][j] = mat(i,j);
		}
	}
}

Matrix::Matrix(double* table, const std::pair<int, int> size)
{
	if (size.first <= 0 || size.second <= 0)
		throw std::exception("Wrong size of table");
	if (table == nullptr)
		throw std::exception("Initializaing table is empty");
	this->sizeH = size.first;
	this->sizeW = size.second;
	this->tab = new double*[this->sizeH];
	int iterator = 0;
	for (size_t i = 0; i < this->sizeH; i++)
	{
		this->tab[i] = new double[this->sizeW];
		for (size_t j = 0; j < this->sizeW; j++)
		{
			this->tab[i][j] = std::move(table[iterator]);
			iterator++;
		}
	}
}

std::pair<int, int> Matrix::getSize()
{
	return std::make_pair(this->sizeH,this->sizeW);
}

double& Matrix::operator()(int x, int y)
{
	if (x < 0 || x > (this->sizeH - 1) || y < 0 || y > (this->sizeW - 1))
		throw new std::exception("Wrong index in () matrix operator");
	return this->tab[x][y];
}

std::ostream& operator<< (std::ostream& os, Matrix& mat)
{
	auto size = mat.getSize();
	for (int i = 0; i < size.first; i++)
	{
		for (int j = 0; j < size.second; j++)
		{
			os << mat(i,j) << " ";
		}
		os << std::endl;
	}
	return os;
}

Matrix& operator/ (Matrix& mat, double scalar)
{
	if (scalar == 0)
		throw std::exception("Dividing matix by 0\n");
	auto size = mat.getSize();
	const size_t h = size.first;
	const size_t w = size.second;
	for (size_t i = 0; i < h; i++)
	{
		for (size_t j = 0; j < w; j++)
		{
			mat(i, j) /= scalar;
		}
	}
	return mat;
}

Matrix& operator+ (Matrix& mat, double scalar)
{
	auto size = mat.getSize();
	const size_t h = size.first;
	const size_t w = size.second;
	for (size_t i = 0; i < h; i++)
	{
		for (size_t j = 0; j < w; j++)
		{
			mat(i, j) += scalar;
		}
	}
	return mat;
}

//Matrix2d
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

//Fuctions
/*
Matrix* multiplyVecAndTransposedVec(std::vector<double>& vec1, std::vector<double>& vec2)
{
	const int size1 = vec1.size();
	int size2 = vec2.size();
	if (size1 != size2)
		throw new std::exception("Sizes of 2 vectors are not equal, not possible to multiply them");
		
	Matrix* mat1=new Matrix(size1,size1); //compiler chooses wrong constructor, don't know why

	for (int i = 0; i < size1; i++)
	{
		for (int j = 0; j < size1; j++)
		{
			mat1->operator()(i, j) = vec1[i]*vec2[j];
		}
	}
	return mat1;
}
*/
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
