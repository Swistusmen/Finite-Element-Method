#include <iostream>
#include <string>
#include "Element.h"
#include "RectangleMesh.h"
#include "Reader.h"
#include "Solver.h"
#include <memory>
#include "Algebra.h"
#include <array>
#include <cmath>

int main()
{
	
	auto data = fs::readRectangleMeshInput("a.txt");
	data::RectangleMesh mesh(data);
	
	int numberOfIntegrationalPoints = 3;

	slv::Solver solver(4);
	Matrix eta = solver.getEtaMatrix();
	Matrix ksi = solver.getKsiMatrix();
	const int size = mesh.maxIndexOfElement();
	std::cout << eta << std::endl;
	std::cout << ksi << std::endl;
	
	std::vector<Matrix4d*> localMatricies;
	std::vector<std::array<int, 4>> nodes;
	
	std::cout << size << std::endl;
	for (int i = 1; i <= size; i++)
	{
		Matrix2d jacoby = solver.getJacobyMatrix2(eta, ksi, mesh.getX(i), mesh.getY(i), 0);
		std::cout << jacoby << std::endl << std::endl;
		localMatricies.push_back(solver.getHMatrix(eta, ksi, jacoby, 25));
		nodes.push_back(mesh.getElementNodesIndexes(i));
	}
	std::cout << "AAAAA\n";
	int sizeOfGlobalMAtrix = mesh.getNH()*mesh.getNW();
	Matrix globalMat(sizeOfGlobalMAtrix);
	solver.aggregateGlobalMatrix(globalMat, localMatricies, nodes);
	std::cout << globalMat << std::endl;
	
	system("pause");
	return 0;
}