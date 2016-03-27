/*
 * CartDecomposition_gen_vtk.hpp
 *
 *  Created on: Aug 28, 2015
 *      Author: i-bird
 */
#include "Decomposition/CartDecomposition.hpp"


int main(int argc, char ** argv)
{
	// Initialize the global VCluster
	init_global_v_cluster(&argc,&argv);

	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	//! [Create CartDecomposition vtk gen]
	CartDecomposition<2,float> dec(vcl);

	// Physical domain
	Box<2,float> box({0.0,0.0},{1.0,1.0});

	// division on each direction
	size_t div[2] = {20,20};

	// Define ghost
	Ghost<2,float> g(0.01);

	// boundary conditions
	size_t bc[2] = {PERIODIC,PERIODIC};

	// Decompose and write the decomposed graph
	dec.setParameters(div,box,bc,g);
	dec.decompose();

	// create a ghost border
	dec.calculateGhostBoxes();

	// Write the decomposition
	dec.write("CartDecomposition/out_");

	//! [Create CartDecomposition]

	delete &vcl;
}

