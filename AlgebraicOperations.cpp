#include "Algebra.h"
#include "AlgebraicOperations.h"
#include <algorithm>

void gjInverseMatrix(std::unique_ptr<Matrix>& mat)
{
	const auto size = mat->getSize();
	if (size.first != size.second)
		throw new std::exception("Matrix is not square matrix, cannot be inversed");
	auto I = std::make_unique<Matrix>(size.first);
	for (size_t i = 0; i < size.first; i++)
	{
		I->operator()(i, i) = 1;
	}
	//making triangle up matrix
	for (size_t i = 1; i < size.first; i++)
	{
		for (size_t c = i; c < size.first; c++) 
		{
			double multiplier = mat->operator()(c, i - 1) / mat->operator()(i - 1, i - 1);
			for (size_t j = 0; j < size.first; j++)
			{
				mat->operator()(c, j) = mat->operator()(c, j) - multiplier * mat->operator()(i - 1, j);
				I->operator()(c, j) = I->operator()(c, j) - multiplier * I->operator()(i - 1, j);
			}
		}
	}
	//making diagonal matrix
	for (size_t i = size.first-2; i >=0 &&i<size.first; i--)
	{
		for (size_t c = i; c >= 0 && c < size.first; c--) 
		{
			double multiplier = mat->operator()(c, i + 1) / mat->operator()(i + 1, i + 1); //for now i assume it's not zero
			for (size_t j = size.first - 1; j >= 0 && j < size.first; j--)
			{
				mat->operator()(c, j) = mat->operator()(c, j) - multiplier * mat->operator()(i + 1, j);
				I->operator()(c, j) = I->operator()(c, j) - multiplier * I->operator()(i + 1, j);
			}
		}
	}
	//making unit matrix
	for (size_t i = 0; i < size.first; i++)
	{
		double divider = mat->operator()(i, i);
		for (size_t j = 0; j < size.first; j++)
		{
			mat->operator()(i, j) /= divider;
			I->operator()(i, j) /= divider;
		}
	}
	std::swap(*mat, *I);
	I.reset(nullptr);
}

void solveJacobiEquation(Vector& preassure, Matrix& H, Vector& temperatures)
{
	const auto size = H.getSize();
	auto D = std::make_unique<Matrix>(size.first);
	for (size_t i = 0; i < size.first; i++)
	{
		D->operator()(i, i) = H(i, i);
		H(i, i) = 0;
	}
	gjInverseMatrix(D);

	H *= -1.0;
	for (int i = 0; i < 100; i++)
	{
		auto LX = temperatures * H; 
		LX->operator+=(preassure);
		auto result = (*LX)*(*D);
		temperatures = std::move(result);
	}
}