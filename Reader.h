#include <fstream>
#include <string>
#include <iostream>
#include <tuple>
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

	data::InputArguments readProgramData(std::string filename)
	{
		std::ifstream file(filename);
		data::InputArguments inputArguments;
		std::string currentLine = "";
		std::string::iterator it;
		std::vector<double> numbers;
		if (file.is_open())
		{
			int i = 0;
			while (std::getline(file, currentLine))
			{
				it = currentLine.begin();
				std::string value = "";
				bool comment = false;
				while (it != currentLine.end())
				{
					if ((*it) == '#')
					{
						if (comment == false)
							comment = true;
						else
							comment = false;
					}
					else if ((*it) != '#'&& (*it) != '\n'&& (*it) != ' ')
					{
						if (comment != true)
						{
							value += *it;
						}
					}
					else if (((*it) == ' ')&&(comment==false))
					{
						numbers.push_back(std::stod(value));
						value = "";
					}
					it++;
				}
				inputArguments.fullFill(numbers, i);
				numbers.clear();
				i++;

			}
		}
		else {
			std::cerr << "Error with reading initializing data, abort\n";
			exit(0);
		}
		return inputArguments;
	}

}
	
