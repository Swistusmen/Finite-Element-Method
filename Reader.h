#include <fstream>
#include <string>
#include <iostream>
#include "InputData.h"


namespace fs {
	data::Points* readTable(std::string filepath) {
		data::Points* points=new data::Points;
		std::ifstream reader("points.txt");
		if (reader.is_open())
		{
			reader >> points->size;

			points->x = new double[points->size];
			points->y = new double[points->size];
			const int size=points->size;
			for (int i = 0; i < size; i++)
			{
				reader >> points->x[i];
				reader >> points->y[i];
			}
		}
			return points;
	}

	data::Element* readRectangleElementsFromFile(std::string filepath)
	{
		auto points = readTable(filepath);
		const int size = points->size / 4;
		data::Element* elements = new data::Element[size];
		int pIterator = 0;
		double tabX[4], tabY[4];
		for (int i = 0; i < size; i++)
		{
			int a = 0;
			while (a < 4)
			{
				tabX[a] = points->x[pIterator];
				tabY[a] = points->y[pIterator];
				a++;
				pIterator++;
			}
			elements[i] = data::Element(tabX, tabY);
		}
		std::cout << std::endl;
		return elements;
	}

}
	
