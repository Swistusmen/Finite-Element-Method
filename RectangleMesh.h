#pragma once
#include <iostream>
#include "Element.h"
#include "InputData.h"

namespace data {
	class RectangleMesh {
	public:
		RectangleMesh(data::RectangleMeshInput& data);
		int maxIndexOfNode(); 
		int maxIndexOfElement();
		Node* getElementNodes(int index);
		std::array<int, 4> getElementNodesIndexes(int index);
		double* getX(int index);
		double* getY(int index);

		int getNH();
		int getNW();

	private:
		float h, w; 
		Element** elements;
		Node** nodes;
		int nH; 
		int nW;
	};
}