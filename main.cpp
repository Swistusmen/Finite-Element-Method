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
	auto data = fs::readRectangleMeshInput("a.txt");
	data::RectangleMesh mesh(data);
	
	int numberOfIntegrationalPoints = 3;

	slv::Solver solver(4);
	auto eta (solver.getLocalMatrixOfLocalTransformation(slv::LocalType::ETA));
	auto ksi(solver.getLocalMatrixOfLocalTransformation(slv::LocalType::KSI));
	const int size = mesh.maxIndexOfElement();
	std::cout << "Eta matrix\n"<<*eta << std::endl;
	std::cout << "Ksi matrix\n"<< *ksi << std::endl;
	
	std::vector<slv::MatUPtr> localHMatricies;
	//std::vector<slv::MatUPtr> localCMatricies;
	std::vector<std::array<int, 4>> nodes;
	
	/*
	data::Elem4 elem;
	double ro = 5;
	double temp = 700.0;
	*/

	std::cout << size << std::endl;
	for (int i = 1; i <= size; i++)
	{
		//auto jacoby (solver.getJacobyMatrix2(eta, ksi, mesh.getX(i), mesh.getY(i), 0));
		//std::cout << *jacoby << std::endl << std::endl;
		
		localHMatricies.push_back(std::move(solver.getHMatrix(eta, ksi, mesh.getX(i),mesh.getY(i), 25)));
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