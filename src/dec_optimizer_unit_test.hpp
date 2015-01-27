/*
 * dec_optimize.hpp
 *
 *  Created on: Jan 16, 2015
 *      Author: i-bird
 */

#ifndef DEC_OPTIMIZE_HPP_
#define DEC_OPTIMIZE_HPP_


#include "Graph/CartesianGraphFactory.hpp"
#include "map_graph.hpp"
#include "metis_util.hpp"
#include "dec_optimizer.hpp"

#undef GS_SIZE
#define GS_SIZE 8

BOOST_AUTO_TEST_SUITE( dec_optimizer_test )

BOOST_AUTO_TEST_CASE( dec_optimizer_test_use)
{
	CartesianGraphFactory<3,Graph_CSR<nm_v,nm_e>> g_factory;
	CartesianGraphFactory<3,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

	// Cartesian grid
	std::vector<size_t> sz;
	sz.push_back(4);
	sz.push_back(4);
	sz.push_back(1);

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

	// optimize

	dec_optimizer<3,Graph_CSR<nm_v,nm_e>> d_o(g,sz);

	grid_key_dx<3> keyZero(0,0,0);
	d_o.optimize<nm_v::sub_id,nm_v::id>(keyZero,g);

	// Write the VTK file

	VTKWriter<Graph_CSR<nm_part_v,nm_part_e>> vtk(gp);
	vtk.write("vtk_partition.vtk");
	VTKWriter<Graph_CSR<nm_v,nm_e>> vtk2(g);
	vtk2.write("vtk_partition2.vtk");
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* DEC_OPTIMIZE_HPP_ */
