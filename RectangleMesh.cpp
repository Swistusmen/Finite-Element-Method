#include <iostream>
#include "Element.h"
#include "RectangleMesh.h"
using namespace data;

Element& RectangleMesh::getElement(int ID)
{
	return elements[ID-1];
}

Element* RectangleMesh::getNodesID(int ID)
{
		Element* temp = new Element[4];
		int i = 0;
		int ID2 = ID;
		while (ID2 > (nH - 1))
		{
			ID2 -= (nH - 1);
			i++;
		}
		int first = i * nH + ID2-1; //gdyby chodzilo o numer indexu +1, ale chodzi o miejsce w tablicy
		
		temp[0] = this->elements[first];
		temp[1] = this->elements[first + 1];
		temp[2] = this->elements[first + nH];
		temp[3] = this->elements[first + 1 + nH];
		
		return temp;
}

int RectangleMesh::maxIndexOfElement()
{
	return (this->nH - 1)*(this->nW - 1);
}

RectangleMesh::RectangleMesh(float H, float W, int nH, int nW, int* tab, int sizeofTab)
{
	this->h = H;
	this->w = W;
	this->nH = nH;
	this->nW = nW;
	//this->elements = new Element[sizeofTab / 2];
	/*
	for (int i = 0, j=0; i < sizeofTab/2 ; i++, j+=2)
	{
		this->elements[i] =Element(tab[j], tab[j + 1]);
	}
	*/
	
}


RectangleMesh::RectangleMesh(float H, float W, int nH, int nW, Element* elements)
{
	this->h = H;
	this->w = W;
	this->nH = nH;
	this->nW = nW;
	this->elements = std::move(elements);
}