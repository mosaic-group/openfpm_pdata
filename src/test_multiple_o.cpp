/*
 * test_multiple_o.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: i-bird
 *
 *
 *  It just test that the compilation with multiple translation unit (*.o) does not
 *  produce error, if we have duplicated symbol in the translation unit we will get error
 *
 */

#include "Vector/vector_dist.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"

void f()
{
	// Ghost
	Ghost<3,float> g(0.01);

	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t sz[3];
	sz[0] = 100;
	sz[1] = 100;
	sz[2] = 100;

	vector_dist<3,float, aggregate<float>, CartDecomposition<3,float> > vd(4096,domain,bc,g);
	grid_dist_id<3, float, aggregate<float[3]>, CartDecomposition<3,float>> g_dist(sz,domain,g);
}


