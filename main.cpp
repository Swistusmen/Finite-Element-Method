#include <iostream>
#include <string>
#include "Element.h"
#include "RectangleMesh.h"
#include <Eigen/Core>
#include "Reader.h"
#include "Solver.h"
#include <memory>
#include "Algebra.h"
#include <array>
#include <cmath>

//TODO: implement RAII
//TODO: make code more elegant
//TODO: make nodes semi-separate from the Elements- for saving memory- 2 neigbour elements duplicate 2 nodes
//TODO: refactor RectangleMesh class	

int main()
{
	
	auto data = fs::readRectangleMeshInput("a.txt");
	data::RectangleMesh mesh(data);
	
	int numberOfIntegrationalPoints = 3;

	slv::Solver solver;
	MatrixXd eta = solver.getEtaMatrixG2(std::pow(numberOfIntegrationalPoints,2));
	MatrixXd ksi = solver.getKsiMatrixG2(std::pow(numberOfIntegrationalPoints, 2));
	std::cout << eta << std::endl;
	std::cout << ksi << std::endl;

	/*
	MatrixXd ksi = solver.getKsiMatrixG2(std::pow(numberOfIntegrationalPoints,2));
	const int size = mesh.maxIndexOfElement();
	std::cout << eta << std::endl;
	std::cout << ksi << std::endl;
	std::vector<Matrix4d*> localMatricies;
	std::vector<std::array<int, 4>> nodes;
	
	for (int i = 1; i <= size; i++)
	{
		Matrix2d jacoby = solver.getJacobyMatrix2(eta, ksi, mesh.getX(i), mesh.getY(i), 0);
		localMatricies.push_back(solver.getHMatrix(eta, ksi, jacoby, 25));
		nodes.push_back(mesh.getElementNodesIndexes(i));
	}
	
	int sizeOfGlobalMAtrix = mesh.getNH()*mesh.getNW();
	MatrixXd globalMat(sizeOfGlobalMAtrix);
	solver.aggregateGlobalMatrix(globalMat, localMatricies, nodes);
	std::cout << globalMat << std::endl;

	
	*/
	system("pause");
	return 0;
}