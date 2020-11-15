#include <iostream>
#include "Element.h"

using namespace data;

int Element::numberOfElements = 1;

int Element::setIndex()
{
	int a = numberOfElements;
	numberOfElements += 1;
	return a;
}

int Element::getIndex()
{
	return this->index;
}

Element::Element(double* tabX, double* tabY)
{
	this->nodes = new Node[4];
	for (int i = 0; i < 4; i++)
	{
		this->nodes[i].X = tabX[i];
		this->nodes[i].Y = tabY[i];
	}
	this->index = setIndex();
}

double* Element::getX()
{
	double* ret = new double[4];
	double a;
	for (int i = 0; i < 4; i++)
	{
		this->nodes;
		a = this->nodes[i].X;
		ret[i] = a;
	}
	return ret;
}

double* Element::getY()
{
	double* ret = new double[4];
	for (int i = 0; i < 4; i++)
	{
		ret[i] = this->nodes[i].Y;
	}
	return ret;
}

