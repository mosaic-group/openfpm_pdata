/*
 * Metis_gen_vtk.hpp
 *
 *  Created on: Aug 29, 2015
 *      Author: i-bird
 */

#ifndef VTK_METIS_GEN_VTK_CPP_
#define VTK_METIS_GEN_VTK_CPP_

#include <iostream>
#include "Graph/CartesianGraphFactory.hpp"
#include "Graph/map_graph.hpp"
#include "Decomposition/Distribution/metis_util.hpp"
#include "SubdomainGraphNodes.hpp"

int main(int argc, char ** argv)
{
	CartesianGraphFactory<2,Graph_CSR<nm_v<2>,nm_e>> g_factory;

	// Cartesian grid
	size_t sz[2] = {20,20};

	// Box
	Box<2,float> box({0.0,0.0},{1.0,1.0});

	const size_t bc[] = {NON_PERIODIC,NON_PERIODIC};

	// Graph to decompose

	Graph_CSR<nm_v<2>,nm_e> g = g_factory.construct<nm_e::communication,NO_VERTEX_ID,float,1,0,1>(sz,box,bc);

	// Convert the graph to metis

	Metis<Graph_CSR<nm_v<2>,nm_e>> met(g,4);

	// decompose

	met.decompose<nm_v_id>();

	// Write the decomposition

	VTKWriter<Graph_CSR<nm_v<2>,nm_e>,VTK_GRAPH> vtk(g);
	vtk.write("Metis/vtk_partition.vtk");
}


#endif /* VTK_METIS_GEN_VTK_CPP_ */
