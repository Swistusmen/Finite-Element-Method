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

//TODO: implement RAII
//TODO: make code more elegant
//TODO: make nodes semi-separate from the Elements- for saving memory- 2 neigbour elements duplicate 2 nodes
//TODO: refactor RectangleMesh class	

int main()
{
	
	auto data = fs::readRectangleMeshInput("a.txt");
	data::RectangleMesh mesh(data);
	
	slv::Solver solver;
	Matrix4d eta = solver.getEtaMatrixG2();
	Matrix4d ksi = solver.getKsiMatrixG2();
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
	for (int i = 0; i < 9; i++)
	{
		//std::cout << *localMatricies[i] << std::endl;
	}
	int sizeOfGlobalMAtrix = mesh.getNH()*mesh.getNW();
	MatrixXd globalMat(sizeOfGlobalMAtrix);
	solver.aggregateGlobalMatrix(globalMat, localMatricies, nodes);
	std::cout << globalMat << std::endl;

	
	/*
	MatrixXd mat(9);
	std::cout << mat << std::endl;
	
	mat(3, 4) = 1.5;
	mat(5, 1) = 3.4;
	mat(7, 5) = 3.14;
	std::cout << mat << std::endl;
	*/
	system("pause");
	return 0;
}