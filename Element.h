#pragma once
#include <iostream>

namespace data {

	struct Node {
		double X;
		double Y;
	};

	class Element {
	public:
		Element(double* tabX,double* tabY);
		Element() = default; //by moc alokowac miejsce na stercie

		int getIndex();
		double* getX();
		double* getY();

		static int setIndex();
	private:
		static int numberOfElements;
		Node* nodes = nullptr;
		int index;
	};
}
