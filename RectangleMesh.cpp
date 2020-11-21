#include <iostream>
#include "Element.h"
#include "RectangleMesh.h"
#include <utility>
#include <array>
#include <memory>
using namespace data;

std::array<int,4> RectangleMesh::getElementNodesIndexes(int index)
{
	auto a = this->getElementNodes(index);

	std::array<int, 4> myTab = { a[0].index,a[2].index,a[3].index,a[1].index };
	return myTab;
}

int RectangleMesh::getNH()
{
	return this->nH;
}

int RectangleMesh::getNW()
{
	return this->nW;
}

Node* RectangleMesh::getElementNodes(int index)
{
	if ((index > this->maxIndexOfElement()) || (index < 0))
		throw new std::exception("There is no such an index of element in rectangle Mesh");
	Node*myTab = new Node[4];
	int index2 = index;
	int i = 0;
	while (index2 > (nH - 1))
	{
		index2 -= (nH - 1);
		i++;
	}

	myTab[0] = this->nodes[i][index2-1];
	myTab[1] = this->nodes[i][index2];
	myTab[2] = this->nodes[i+1][index2-1];
	myTab[3] = this->nodes[i+1][index2];
	return myTab;
}

int RectangleMesh::maxIndexOfNode()
{
	return this->nH *this->nW ;
}

int RectangleMesh::maxIndexOfElement()
{
	return (this->nH - 1)*(this->nW - 1);
}

RectangleMesh::RectangleMesh(data::RectangleMeshInput& mesh)
{
	this->nW = mesh.nW;
	this->nH = mesh.nH;
	this->h = mesh.H;
	this->w = mesh.W;
	
	double wStep = double(w / (nW-1));
	double hStep = double(h / (nH-1));
	this->nodes = new Node*[nW];
	double wPoint = 0.0;
	double hPoint = 0.0;
	
	for (size_t i = 0; i < nW; i++,wPoint+=wStep)
	{
		this->nodes[i] = new Node[nH];
		for (size_t j = 0; j < nH; ++j, hPoint+=hStep)
		{
			this->nodes[i][j].X = wPoint;
			this->nodes[i][j].Y = hPoint;
		}
		hPoint = 0;
	}
	
	this->elements= new data::Element*[nW - 1];
	for (size_t i = 0; i < nW - 1; i++)
	{
		this->elements[i] = new data::Element[nH - 1];
	}
}

double* RectangleMesh::getX(int index)
{
	auto temp = this->getElementNodes(index);
	double* myTab = new double[4];

	myTab[0] = temp[0].X;
	myTab[1] = temp[2].X;
	myTab[2] = temp[3].X;
	myTab[3] = temp[1].X; //it's because of differecness between id number of node, and order of numerating node in mesh

	return myTab;
}

double* RectangleMesh::getY(int index)
{
	auto temp = this->getElementNodes(index);
	double* myTab = new double[4];
	myTab[0] = temp[0].Y;
	myTab[1] = temp[2].Y;
	myTab[2] = temp[3].Y;
	myTab[3] = temp[1].Y;
	return myTab;
}