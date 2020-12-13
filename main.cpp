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

//TODO make copy and move assignment for matrix and vector class

int main()
{
	int scheme = 2;
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
	std::vector<double> HbcMultipliers;
	HbcMultipliers.push_back(25.0);
	for (int i = 1; i <= size; i++)
	{
		localHMatricies.push_back(std::move(solver.getMatrixForElement(slv::MatrixType::H, mesh.getX(i),mesh.getY(i), HMultipliers)));
		localCMatricies.push_back(std::move(solver.getMatrixForElement(slv::MatrixType::C, mesh.getX(i), mesh.getY(i), CMultipliers)));
		std::cout << "Element " << i << std::endl;
		auto d = solver.getBoundaryMatrixForElement(mesh.getX(i), mesh.getY(i), HbcMultipliers);
		//std::cout << *d << std::endl;
		nodes.push_back(mesh.getElementNodesIndexes(i));
		//std::cout << *localCMatricies.back() << std::endl;
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