#pragma once
#include <vector>
#include <iostream>
#include <utility>
#include <memory>
#include <functional>
#include <array>

namespace slv {
	enum LocalType { ETA, KSI, SHAPE2D, SHAPE1D };
	enum MatrixType { H, C};

	struct LocalOperations{
		std::vector < std::function<double(double)>> Ksi;
		std::vector < std::function<double(double)>> Eta;
		std::vector < std::function<double(double)>> Shape1D;
		std::vector < std::function<double(double, double)>> Shape2D;

		LocalOperations() {
			Ksi.push_back([](double a) {return -0.25*(1 - a); });
			Ksi.push_back([](double a) {return 0.25*(-a + 1); });
			Ksi.push_back([](double a) {return 0.25*(a + 1); });
			Ksi.push_back([](double a) {return -0.25*(1 + a); });
			Eta.push_back([](double a) {return -0.25*(1 - a); });
			Eta.push_back([](double a) {return -0.25*(1 + a); });
			Eta.push_back([](double a) {return 0.25*(a + 1); });
			Eta.push_back([](double a) {return 0.25*(-a + 1); });
			Shape1D.push_back([](double e) {return 0.5*(1 - e); });
			Shape1D.push_back([](double e) {return 0.5*(1 + e); });
			Shape2D.push_back([](double e, double n) {return 0.25*(1 - e)*(1 - n); });
			Shape2D.push_back([](double e, double n) {return 0.25*(1 + e)*(1 - n); });
			Shape2D.push_back([](double e, double n) {return 0.25*(1 + e)*(1 + n); });
			Shape2D.push_back([](double e, double n) {return 0.25*(1 - e)*(1 + n); });
		};

		std::vector <std::function<double(double)>>& getFunctions(LocalType type) {
			return type == ETA ? Eta : Ksi;
		}
	};

	struct IntegrationalPoints {
		std::vector<double> twoPointScheme;
		std::vector<double> threePointScheme;
		std::vector<double> fourPointScheme;

		IntegrationalPoints() {
				twoPointScheme.push_back(-1 / sqrt(3.0));
				twoPointScheme.push_back(1 / sqrt(3.0));
				threePointScheme.push_back(-std::sqrt(3.0 / 5.0));
				threePointScheme.push_back(0.0);
				threePointScheme.push_back(std::sqrt(3.0 / 5.0));
				fourPointScheme.push_back(-1.0*std::sqrt(3.0 / 7.0 + 2.0 / 7.0*std::sqrt(1.2)));
				fourPointScheme.push_back(-1.0*std::sqrt(3.0 / 7.0 - 2.0 / 7.0*std::sqrt(1.2)));
				fourPointScheme.push_back(std::sqrt(3.0 / 7.0 - 2.0 / 7.0*std::sqrt(1.2)));
				fourPointScheme.push_back(std::sqrt(3.0 / 7.0 + 2.0 / 7.0*std::sqrt(1.2)));
		}

		const std::vector<double>& getDependingonScheme(int scheme) {
			if (scheme == 2)return twoPointScheme;
			if (scheme == 3)return threePointScheme;
			if (scheme == 4)return fourPointScheme;
		}
	};

	struct BoundaryConditionConfuguration
	{
		std::array<int, 4> Left= { 0,3,0,1 } ;
		std::array<int, 4> Right=  {1,2,1,0 } ;
		std::array<int, 4> Bottom= { 0,1,0,1 } ;
		std::array<int, 4> Top= { 2,3,1,0 } ;
	};
}
