#pragma once
#include "Algebra.h"
#include <iostream>
#include <vector>
#include <memory>

#define EQUATION_ITERATIONS 10

void gjInverseMatrix(std::unique_ptr<Matrix>& mat);

void solveJacobiEquation(Vector& presassure, Matrix& H, Vector& temperatures);



