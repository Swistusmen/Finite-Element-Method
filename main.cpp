#include <iostream>
#include <string>

#include "Element.h"
#include "RectangleMesh.h"
#include "Reader.h"
#include "Solver.h"

#include "Algebra.h"
#include <memory>
#include <array>
#include <cmath>
#include <algorithm>
#include "AlgebraicOperations.h"


int main()
{
	auto InputData = fs::readProgramData("data.txt");

	int scheme = 2;
	slv::LocalOperations localOperations;
	slv::IntegrationalPoints integrationalPoints;

	data::RectangleMesh mesh(InputData.meshData);

	slv::Solver solver(localOperations, integrationalPoints, scheme, InputData.meshData);

	const int size = mesh.maxIndexOfElement();

	std::vector<slv::MatUPtr> localHMatricies;
	std::vector<slv::MatUPtr> localCMatricies;
	std::vector<slv::VecUPtr> preassures;
	std::vector<std::array<int, 4>> nodes;
	std::cout << size << std::endl;

	int noNodes = mesh.getNH()*mesh.getNW();
	Vector t0(noNodes);
	t0.fullFill(InputData.t0);

	Matrix globalHMat(noNodes);
	Matrix globalCMat(noNodes);
	Vector globalVecOfPreassure(noNodes);
	Matrix HPrim(noNodes);
	Vector PPrim(noNodes);

	int iterations = InputData.experimentParameters[0] / InputData.experimentParameters[1];

	for (int i = 1; i <= size; i++)
	{
		double *X = mesh.getX(i); //need to be refactiored, can be std::array<double,4>
		double *Y = mesh.getY(i);
		localHMatricies.push_back(std::move(solver.getMatrixForElement(slv::MatrixType::H, X, Y, InputData.HMultipliers)));
		localHMatricies.back()->operator+=(*(solver.getBoundaryMatrixForElement(X, Y, InputData.HBCMultipliers)));
		localCMatricies.push_back(std::move(solver.getMatrixForElement(slv::MatrixType::C, X, Y, InputData.CMultipliers)));
		preassures.push_back(std::move(solver.getPreassureVectorForElement(X, Y, InputData.Preassure)));
		nodes.push_back(mesh.getElementNodesIndexes(i));
	}

	solver.aggregateGlobalMatrix(globalHMat, localHMatricies, nodes);
	std::cout << "H global\n" << globalHMat << std::endl;

	solver.aggregateGlobalMatrix(globalCMat, localCMatricies, nodes);
	std::cout << "C global\n" << globalCMat << std::endl;

	solver.aggregateGlobalVector(globalVecOfPreassure, preassures, nodes);

	localCMatricies.clear();
	localHMatricies.clear();
	nodes.clear();
	preassures.clear();

	std::vector<std::pair<double, double>> MinMaxTemperatures;
	for (int c = 0; c < iterations; c++)
	{
		std::cout << "Iteration: " << c << std::endl;
		Matrix Ccurrent(globalCMat);
		Matrix HCurrent(globalHMat);
		Ccurrent.abs();

		HPrim = HCurrent;
		Matrix CCopy(Ccurrent);
		HPrim += Ccurrent / InputData.experimentParameters[1];

		CCopy *= 1.0 / InputData.experimentParameters[1];
		PPrim = CCopy.multiply(CCopy, t0);
		PPrim += globalVecOfPreassure;
		//std::cout << "PPrim\n" << PPrim << std::endl;

		//solve set of equations
		solveJacobiEquation(PPrim, HPrim, t0);
		auto MinMax = std::minmax_element(t0.tab.begin(), t0.tab.end());
		std::cout << *MinMax.first << " " << *MinMax.second << std::endl;
		MinMaxTemperatures.push_back({ *MinMax.first ,*MinMax.second });
	}
	for (int c = 0; c < iterations; c++)
	{
		std::cout << MinMaxTemperatures[c].first << " " << MinMaxTemperatures[c].second << std::endl;
	}

	system("pause");
	return 0;
}