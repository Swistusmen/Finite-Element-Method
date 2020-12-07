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
	/*
	Matrix mat(3);
	std::cout << mat << std::endl;
	*/
	
	auto data = fs::readRectangleMeshInput("a.txt");
	data::RectangleMesh mesh(data);
	
	int numberOfIntegrationalPoints = 3;

	slv::Solver solver(2);
	Matrix eta (solver.getEtaMatrix());
	Matrix ksi (solver.getKsiMatrix());
	const int size = mesh.maxIndexOfElement();
	std::cout << eta << std::endl;
	std::cout << ksi << std::endl;
	
	std::vector<Matrix*> localHMatricies;
	std::vector<Matrix*> localCMatricies;
	std::vector<std::array<int, 4>> nodes;
	
	data::Elem4 elem;
	double ro = 5;
	double temp = 700.0;


	std::cout << size << std::endl;
	for (int i = 1; i <= size; i++)
	{
		Matrix jacoby (solver.getJacobyMatrix2(eta, ksi, mesh.getX(i), mesh.getY(i), 0));
		std::cout << jacoby << std::endl << std::endl;
		localHMatricies.push_back(solver.getHMatrix(eta, ksi, jacoby, 25));
		//localCMatricies.push_back(solver.getCMatrix(elem,jacoby,ro,temp));
		nodes.push_back(mesh.getElementNodesIndexes(i));
	}
	int sizeOfGlobalMatrix = mesh.getNH()*mesh.getNW();
	Matrix globalHMat(sizeOfGlobalMatrix);
	solver.aggregateGlobalMatrix(globalHMat, localHMatricies, nodes);
	std::cout << globalHMat << std::endl;
	
	/*
	Matrix globalCMat(sizeOfGlobalMatrix);
	solver.aggregateGlobalMatrix(globalCMat, localCMatricies, nodes);
	std::cout << globalCMat << std::endl;
	*/
	system("pause");
	return 0;
}