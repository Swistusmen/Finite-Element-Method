#pragma once
#include <iostream>

namespace data {

	struct Node {
		double X;
		double Y;
		int index;
		static int numberOfNodes;

		static int setIndex();
		Node(double X, double Y)
		{
			this->X = X;
			this->Y = Y;
			index = setIndex();
		}

		Node()
		{
			index = setIndex();
		}
	};

	class Element {
	public:
		Element() = default; //by moc alokowac miejsce na stercie

		int getIndex();

		static int setIndex();
	private:
		static int numberOfElements;
		int index;
	};
}