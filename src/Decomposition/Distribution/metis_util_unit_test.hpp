/*
 * metis_util_unit_test.hpp
 *
 *  Created on: Dec 21, 2014
 *      Author: i-bird
 */

#ifndef METIS_UTIL_UNIT_TEST_HPP_
#define METIS_UTIL_UNIT_TEST_HPP_

#include "Graph/CartesianGraphFactory.hpp"
#include "Graph/map_graph.hpp"
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
	Vcluster & v_cl = *global_v_cluster;

	if (v_cl.getProcessingUnits() != 3)
		return;

	if (v_cl.getProcessUnitID() != 0)
		return;

	CartesianGraphFactory<3,Graph_CSR<nm_v,nm_e>> g_factory;
	CartesianGraphFactory<3,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

	// Cartesian grid
	size_t sz[3] = {GS_SIZE,GS_SIZE,GS_SIZE};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions, non periodic
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// Graph to decompose
	Graph_CSR<nm_v,nm_e> g = g_factory.construct<nm_e::communication,NO_VERTEX_ID,float,2,0,1,2>(sz,box,bc);

	// Processor graph
	Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.construct<NO_EDGE,NO_VERTEX_ID,float,2>(sz,box,bc);

	// Convert the graph to metis

	Metis<Graph_CSR<nm_v,nm_e>> met(g,8);

	// decompose

	met.decompose<nm_part_v::id>(gp);
	met.decompose<nm_v::proc_id>();

	// Write the VTK file

	VTKWriter<Graph_CSR<nm_part_v,nm_part_e>,VTK_GRAPH> vtk(gp);
	vtk.write("vtk_metis_util_gp.vtk");

	VTKWriter<Graph_CSR<nm_v,nm_e>,VTK_GRAPH> vtk2(g);
	vtk2.write("vtk_metis_util_g.vtk");

	// check that match

	bool test = compare("vtk_metis_util_gp.vtk","src/Decomposition/Distribution/test_data/vtk_metis_util_gp_test.vtk");
	bool test2 = compare("vtk_metis_util_g.vtk","src/Decomposition/Distribution/test_data/vtk_metis_util_g_test.vtk");
	BOOST_REQUIRE_EQUAL(true,test);
	BOOST_REQUIRE_EQUAL(true,test2);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
