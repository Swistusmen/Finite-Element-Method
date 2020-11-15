#include <iostream>
#include <string>
#include "Element.h"
#include "RectangleMesh.h"
#include <Eigen/Core>
#include "Reader.h"
#include "Solver.h"
#include <memory>
#include "Algebra.h"

//TODO: implement RAII
//TODO: make code more elegant
//TODO: make nodes semi-separate from the Elements- for saving memory- 2 neigbour elements duplicate 2 nodes
//TODO: refactor RectangleMesh class

int main()
{
	data::Element* tabOfElements = fs::readRectangleElementsFromFile("");
	data::RectangleMesh mesh(1.0, 1.0, 3, 3, tabOfElements);
	slv::Solver solver;
	Matrix4d eta = solver.getEtaMatrixG2();
	Matrix4d ksi = solver.getKsiMatrixG2();
	const int size = mesh.maxIndexOfElement();
	
	for (int i = 1; i <= size; i++)
	{
		auto currentElement = mesh.getElement(i);
		Matrix2d jacoby = solver.getJacobyMatrix2(eta, ksi, currentElement.getX(), currentElement.getY(), 0);
		auto H = solver.getHMatrix(eta, ksi, jacoby, 30);
		std::cout << *H << std::endl;
	}
	
	system("pause");
	return 0;
}