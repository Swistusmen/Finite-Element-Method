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

	data::RectangleMeshInput& readRectangleMeshInput(std::string&& file)
	{
		std::ifstream reader(std::move(file));
		data::RectangleMeshInput mesh;
		if (reader.is_open())
		{
			reader >> mesh.nH;
			reader >> mesh.nW;
			reader >> mesh.H;
			reader >> mesh.W;

			return mesh;
		}
		throw new std::exception("Error in reading from file");
	}

}
	
