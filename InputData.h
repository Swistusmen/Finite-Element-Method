#pragma once
#include <iostream>
#include <vector>

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

	struct RectangleMeshInput {
		int nH, nW;
		double H, W;
		RectangleMeshInput()= default;

		void init(std::vector<double> vec)
		{
			nH = (int) vec[0];
			nW = (int)vec[1];
			H = vec[2];
			W = vec[3];
		}
	};

	//awful code, to improve and rewrite some day- enum and algorithm
	struct InputArguments {
		data::RectangleMeshInput meshData;
		std::vector<double> HMultipliers;
		std::vector<double> HBCMultipliers;
		std::vector<double> CMultipliers;
		std::vector<double> Preassure;
		double t0;
		std::vector<double> experimentParameters;

		void fullFill(std::vector<double> vec , int n)
		{
			switch (n)
			{
			case 0: {
				meshData.init(vec);
			}break;
			case 1: {
				HMultipliers = vec;
			}break;
			case 2: {
				HBCMultipliers = vec;
			}break;
			case 3: {
				CMultipliers = vec;
			}break;
			case 4: {
				Preassure = vec;
			}break;
			case 5: {
				t0 = vec[0];
			}break;
			case 6: {
				experimentParameters = vec;
			}break;
			}
		}
	};
}
