#include <iostream>
#include "Algebra.h"

////////////////////// Constructors
Matrix::Matrix(int size, int size2)
{
	if (size <= 0 || size2 < 0)
		throw new std::exception("Matrix constructor error");
	if (size2 == 0)
		size2 = size;
	this->sizeH = size;
	this->sizeW = size2 == 0 ? size : size2;
	for (int i = 0; i < this->sizeH; i++)
	{
		this->tab.push_back(std::vector<double>());
		for (int j = 0; j < this->sizeW; j++)
		{
			this->tab[i].push_back(0);
		}
	}
}

Matrix::Matrix(const Matrix& mat)
{
	auto size= mat.getSize();
	this->sizeH = size.first;
	this->sizeW = size.second;

	for (int i = 0; i < this->sizeH; i++)
	{
		this->tab.push_back(std::vector<double>());
		for (int j = 0; j < this->sizeW; j++)
		{
			this->tab[i].push_back(mat(i,j));
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
	int iterator = 0;

	for (int i = 0; i < this->sizeH; i++)
	{
		this->tab.push_back(std::vector<double>());
		for (int j = 0; j < this->sizeW; j++)
		{
			this->tab[i].push_back(std::move(table[iterator]));
			iterator++;
		}
	}
}

Matrix::Matrix(Matrix&& mat)
{
	auto size = mat.getSize();
	this->sizeH = size.first;
	this->sizeW = size.second;

	for (int i = 0; i < this->sizeH; i++)
	{
		this->tab.push_back(std::vector<double>());
		for (int j = 0; j < this->sizeW; j++)
		{
			this->tab[i].push_back(std::move(mat(i, j)));
		}
	}
}

/////////////////////////// The most used functions
std::pair<int, int> Matrix::getSize() const
{
	return std::make_pair(this->sizeH,this->sizeW);
}

double& Matrix::operator()(int x, int y)
{
	if (x < 0 || x > (this->sizeH - 1) || y < 0 || y > (this->sizeW - 1))
		throw new std::exception("Wrong index in () matrix operator");
	return this->tab[x][y];
}

double Matrix::operator()(int x, int y) const
{
	if (x < 0 || x >(this->sizeH - 1) || y < 0 || y >(this->sizeW - 1))
		throw new std::exception("Wrong index in () matrix operator");
	return this->tab[x][y];
}

std::ostream& operator<< (std::ostream& os, const Matrix& mat)
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
///////////////////////////////////// Operators 
Matrix& Matrix::operator=(const Matrix& mat)
{
	auto size = this->getSize();
	const size_t sizeH = size.first;
	const size_t sizeW = size.first;
	if ((this->sizeH != sizeH) && (this->sizeH != sizeH))
		throw new std::exception("Copy operations is inpossible due to differnet size");
	for (size_t i = 0; i < sizeH; i++)
	{
		for (size_t j = 0; j < sizeW; j++)
		{
			this->operator()(i, j) = mat(i, j);
		}
	}
	return *this;
}

Matrix& Matrix::operator=( Matrix&& mat)
{
	auto size = this->getSize();
	const size_t sizeH = size.first;
	const size_t sizeW = size.first;
	if ((this->sizeH != sizeH) && (this->sizeH != sizeH))
		throw new std::exception("Copy operations is inpossible due to differnet size");
	for (size_t i = 0; i < sizeH; i++)
	{
		for (size_t j = 0; j < sizeW; j++)
		{
			this->operator()(i, j) = std::move(mat(i, j));
		}
	}
	return *this;
}

Matrix& Matrix::operator+=(Matrix& mat)
{
	auto size = mat.getSize();
	if ((this->sizeH != size.first) || (this->sizeW != size.second))
	{
		throw new std::exception("Wrong sizes of operations");
	}
	for (size_t i = 0; i < size.first; i++)
	{
		for (size_t j = 0; j < size.second; j++)
		{
			this->tab[i][j] += mat(i, j);
		}
	}
	return *this;
}

Matrix& Matrix::operator*=(double scalar)
{
	const size_t h = this->sizeH, w = this->sizeW;
	for (size_t i = 0; i < h; i++)
	{
		for (size_t j = 0; j < w; j++)
		{
			this->tab[i][j] *= scalar;
		}
	}
	return *this;
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

std::shared_ptr<Matrix> operator+ (Matrix& mat, Matrix& tab)
{
	const auto size1 = mat.getSize();
	const auto size2 = tab.getSize();
	if (size1 != size2)
	{
		throw new std::exception("Wrong sizes of matricies, cannot add them");
	}
	auto newMat = std::make_shared<Matrix>(size1.first, size1.second);
	for (size_t i = 0; i < size1.first; i++)
	{
		for (size_t j = 0; j < size1.second; j++)
		{
			newMat->operator()(i, j) = mat(i, j) + tab(i, j);
		}
	}
	return newMat;
}

std::shared_ptr<Matrix> operator* (Matrix& mat, Matrix& mat1)
{
	const auto size1 = mat.getSize();
	const auto size2 = mat1.getSize();
	if (size1.first != size2.first!= size1.second!=size2.second)
	{
		throw new std::exception("Wrong sizes of matricies, cannot add them. You can multiply only quadratic matrisies");
	}
	auto newMat = std::make_shared<Matrix>(size1.first, size1.second);
	const size_t size = size2.first;
	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			for (int c = 0; c < size; c++) {
				newMat->operator()(i, j) += mat(i, c)*mat1(c, j);
			}
		}
	}
	return newMat;
}

double determinantMat2(Matrix& mat)
{
	auto size = mat.getSize();
	if ((size.first == size.second)&&(size.first != 2))
		throw new std::exception("Matrix has to be 2x2\n");
	return mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
}

void transposeMat2(Matrix& mat)
{
	auto size = mat.getSize();
	if ((size.first == size.second) && (size.first != 2))
		throw new std::exception("Matrix has to be 2x2\n");
	mat(0, 1) = mat(1,0);
	mat(1, 0) = mat(0,1);
}

void inverseMat2(Matrix& mat)
{
	auto size = mat.getSize();
	if ((size.first == size.second) && (size.first != 2))
		throw new std::exception("Matrix has to be 2x2\n");
	double det = determinantMat2(mat);
	mat(0, 0) = mat(1,1) / det;
	mat(1, 1) = mat(0,0) / det;
	mat(0, 1) = mat(1,0) / det;
	mat(1, 0) = mat(0,1) / det;

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

std::unique_ptr<Matrix> vecAndvecTMultiplication(Vector2d& a)
{
	const size_t size = 2;
	auto mat = std::make_unique<Matrix>(size);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			mat->operator()(i, j) = a(i)*a(j);
		}
	}
	return mat;
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

std::unique_ptr<Matrix> vecAndvecTMultiplication(Vector4d& a)
{
	const size_t size = 4;
	auto mat = std::make_unique<Matrix>(size);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			mat->operator()(i, j) = a(i)*a(j);
		}
	}
	return mat;
}
