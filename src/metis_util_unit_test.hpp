/*
 * metis_util_unit_test.hpp
 *
 *  Created on: Dec 21, 2014
 *      Author: i-bird
 */

#ifndef METIS_UTIL_UNIT_TEST_HPP_
#define METIS_UTIL_UNIT_TEST_HPP_

#include "Graph/CartesianGraphFactory.hpp"
#include "map_graph.hpp"
#include "metis_util.hpp"

#undef GS_SIZE
#define GS_SIZE 8

/*!
 *
 * Test node
 *
 */

BOOST_AUTO_TEST_SUITE( Metis_test )

BOOST_AUTO_TEST_CASE( Metis_test_use)
{
	CartesianGraphFactory<3,Graph_CSR<nm_v,nm_e>> g_factory;
	CartesianGraphFactory<3,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

	// Cartesian grid
	size_t sz[3] = {GS_SIZE,GS_SIZE,GS_SIZE};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Graph to decompose

	Graph_CSR<nm_v,nm_e> g = g_factory.construct<nm_e::communication,float,2,0,1,2>(sz,box);

	// Processor graph

	Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.construct<NO_EDGE,float,2>(sz,box);

	// Convert the graph to metis

	Metis<Graph_CSR<nm_v,nm_e>> met(g,16);

	// decompose

	met.decompose<nm_part_v::id>(gp);
	met.decompose<nm_v::id>();

	// Write the VTK file

	VTKWriter<Graph_CSR<nm_part_v,nm_part_e>> vtk(gp);
	vtk.write("vtk_partition.vtk");
	VTKWriter<Graph_CSR<nm_v,nm_e>> vtk2(g);
	vtk2.write("vtk_partition2.vtk");
}

BOOST_AUTO_TEST_SUITE_END()

#endif
