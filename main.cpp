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


int main()
{
	int scheme = 4;
	slv::LocalOperations localOperations;
	slv::IntegrationalPoints integrationalPoints;

	auto data = fs::readRectangleMeshInput("a.txt");
	data::RectangleMesh mesh(data);

	slv::Solver solver(localOperations, integrationalPoints,scheme,data);

	const int size = mesh.maxIndexOfElement();
	
	std::vector<slv::MatUPtr> localHMatricies;
	std::vector<slv::MatUPtr> localCMatricies;
	std::vector<std::array<int, 4>> nodes;

	std::cout << size << std::endl;
	std::vector<double> HMultipliers;
	HMultipliers.push_back(25.0);
	std::vector<double> CMultipliers;
	CMultipliers.push_back(700.0);
	CMultipliers.push_back(7800.0);
	
	for (int i = 1; i <= size; i++)
	{
		double *X = mesh.getX(i); //need to be refactiored, can be std::array<double,4>
		double *Y = mesh.getY(i);
		localHMatricies.push_back(std::move(solver.getMatrixForElement(slv::MatrixType::H, X,Y, HMultipliers)));
		localHMatricies.back()->operator+=(*(solver.getBoundaryMatrixForElement(X, Y, HMultipliers)));
		localCMatricies.push_back(std::move(solver.getMatrixForElement(slv::MatrixType::C, X, Y, CMultipliers)));
		nodes.push_back(mesh.getElementNodesIndexes(i));
	}
	int sizeOfGlobalMatrix = mesh.getNH()*mesh.getNW();
	Matrix globalHMat(sizeOfGlobalMatrix);
	solver.aggregateGlobalMatrix(globalHMat, localHMatricies, nodes);
	std::cout << globalHMat << std::endl;
	
	
	Matrix globalCMat(sizeOfGlobalMatrix);
	solver.aggregateGlobalMatrix(globalCMat, localCMatricies, nodes);
	std::cout << globalCMat << std::endl;
	
	system("pause");
	return 0;
}