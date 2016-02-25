/*
 * dec_optimize.hpp
 *
 *  Created on: Jan 16, 2015
 *      Author: Pietro Incardona
 */

#ifndef DEC_OPTIMIZE_HPP_
#define DEC_OPTIMIZE_HPP_


#include "Graph/CartesianGraphFactory.hpp"
#include "Graph/map_graph.hpp"
#include "metis_util.hpp"
#include "dec_optimizer.hpp"


#undef GS_SIZE
#define GS_SIZE 8

BOOST_AUTO_TEST_SUITE( dec_optimizer_test )

BOOST_AUTO_TEST_CASE( dec_optimizer_test_use_np)
{
	CartesianGraphFactory<3,Graph_CSR<nm_v,nm_e>> g_factory;
	CartesianGraphFactory<3,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

	// Cartesian grid
	size_t sz[3] = {GS_SIZE,GS_SIZE,GS_SIZE};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions, non periodic
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// Graph to decompose
	Graph_CSR<nm_v,nm_e> g = g_factory.construct<nm_e::communication,float,2,0,1,2>(sz,box,bc);

	// Processor graph
	Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.construct<NO_EDGE,float,2>(sz,box,bc);

	// Convert the graph to metis
	Metis<Graph_CSR<nm_v,nm_e>> met(g,16);

	// decompose
	met.decompose<nm_part_v::id>(gp);
	met.decompose<nm_v::id>();

	// optimize
	dec_optimizer<3,Graph_CSR<nm_v,nm_e>> d_o(g,sz);

	grid_key_dx<3> keyZero(0,0,0);
	d_o.optimize<nm_v::sub_id,nm_v::id>(keyZero,g,bc);
}

BOOST_AUTO_TEST_CASE( dec_optimizer_test_use_p)
{
	CartesianGraphFactory<3,Graph_CSR<nm_v,nm_e>> g_factory;
	CartesianGraphFactory<3,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

	// Cartesian grid
	size_t sz[3] = {GS_SIZE,GS_SIZE,GS_SIZE};

	//! Grid info
	grid_sm<3,void> gs(sz);

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions, non periodic
	size_t bc[] = {PERIODIC,PERIODIC,PERIODIC};

	// Graph to decompose
	Graph_CSR<nm_v,nm_e> g = g_factory.construct<nm_e::communication,float,2,0,1,2>(sz,box,bc);

	// Processor graph
	Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.construct<NO_EDGE,float,2>(sz,box,bc);

	bool p[3];

	// Divide in 8 parts your graph

	// decompose
	for (size_t i = 0 ; i < GS_SIZE ; i++)
	{
		p[0] = (i < GS_SIZE/2)?false:true;
		for (size_t j = 0 ; j < GS_SIZE ; j++)
		{
			p[1] = (j < GS_SIZE/2)?false:true;
			for (size_t k = 0 ; k < GS_SIZE ; k++)
			{
				p[2] = (k < GS_SIZE/2)?false:true;
				size_t id = 4*p[2] + 2*p[1] + p[0];

				grid_key_dx<3> key(i,j,k);
				gp.vertex(gs.LinId(key)).template get<nm_part_v::id>() = id;
				g.vertex(gs.LinId(key)).template get<nm_v::id>() = id;
			}
		}
	}

	// optimize
	dec_optimizer<3,Graph_CSR<nm_v,nm_e>> d_o(g,sz);

	grid_key_dx<3> keyZero(0,0,0);

	// Set of sub-domain produced by dec-optimizer
	openfpm::vector<Box<3,size_t>> dec_o;

	// For each sub-domain check the neighborhood processors
	openfpm::vector< openfpm::vector<size_t> > box_nn_processor;

	// key
	grid_key_dx<3> zero;
	zero.zero();

	// gp,p_id,loc_box,box_nn_processor,bc
	d_o.optimize<nm_v::sub_id,nm_v::id>(zero,g,-1,dec_o,box_nn_processor,bc);

	BOOST_REQUIRE_EQUAL(box_nn_processor.size(),8ul);

	for(size_t i = 0 ; i < box_nn_processor.size() ; i++)
	{
		bool nn[] = {false,false,false,false,false,false,false,false};
		BOOST_REQUIRE_EQUAL(box_nn_processor.get(i).size(),7ul);
		for (size_t j = 0 ; j < box_nn_processor.get(i).size(); j++)
		{
			BOOST_REQUIRE(box_nn_processor.get(i).get(j) < 8);
			nn[box_nn_processor.get(i).get(j)] = true;
		}

		// search the neighborhood

		size_t cnt = 0;
		for(size_t i = 0 ; i < 8 ; i++)
		{
			if (nn[i] == false)
				cnt++;
		}

		BOOST_REQUIRE_EQUAL(cnt,1ul);
	}

	// check
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* DEC_OPTIMIZE_HPP_ */
