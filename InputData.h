#pragma once
#include <iostream>

namespace data {

	struct Points {
		double* x=nullptr;
		double* y=nullptr;
		int size = 0;;
	};

	struct Elem4 {
		double tabKsi[4] = { -1 / sqrt(3),1 / sqrt(3),1 / sqrt(3),-1 / sqrt(3) };
		double tabEta[4] = { -1 / sqrt(3),-1 / sqrt(3),1 / sqrt(3),1 / sqrt(3) };
		double tabWage[4] = { 1,1,1,1 };
	};

}
