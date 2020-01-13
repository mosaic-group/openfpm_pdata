/*
 * domain_vtk_gen.cpp
 *
 *  Created on: Aug 30, 2015
 *      Author: Pietro Incardona
 */

#include "config.h"
#include <iostream>
#include "Space/Shape/Box.hpp"
#include "Vector/map_vector.hpp"
#include "VTKWriter/VTKWriter.hpp"

int main(int argc, char ** argv)
{
	//! [Output a vector of boxes]

	// Physical domain
	Box<2,float> box({0.0,0.0},{1.0,1.0});

	// Cell
	Box<2,float> cell = box;

	// division on each direction
	size_t div[2] = {20,20};
	Point<2,float> p_div({20.0,20.0});
	cell /= p_div;

	// create 20 cell on each direction

	openfpm::vector<Box<2,float>> v_box;

	for (size_t i = 0; i <= div[0] ; i++)
	{
		for (size_t j = 0 ; j <= div[1] ; j++)
		{
			Point<2,float> p({(float)i,(float)j});
			Box<2,float> box = cell * p;

			v_box.add(box);
		}
	}

	// write the vector of boxes
	VTKWriter<openfpm::vector<Box<2,float>>,VECTOR_BOX> vtk_box1;
	vtk_box1.add(v_box);
	vtk_box1.write("CartDecomposition/dom_box.vtk");

	//! [Output a vector of boxes]
}


