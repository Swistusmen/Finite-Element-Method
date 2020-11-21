#include <iostream>
#include "Element.h"

using namespace data;

int Element::numberOfElements = 1;
int Node::numberOfNodes = 1;

int Node::setIndex()
{
	int a = numberOfNodes;
	numberOfNodes += 1;
	return a;
}

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

