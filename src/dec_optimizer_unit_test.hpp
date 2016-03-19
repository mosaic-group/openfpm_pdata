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
#include "Decomposition/Distribution/metis_util.hpp"
#include "dec_optimizer.hpp"
#include "util/SimpleRNG.hpp"

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
	Graph_CSR<nm_v,nm_e> g = g_factory.construct<nm_e::communication,NO_VERTEX_ID,float,2,0>(sz,box,bc);

	// Processor graph
	Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.construct<NO_EDGE,NO_VERTEX_ID,float,2>(sz,box,bc);

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
	Graph_CSR<nm_v,nm_e> g = g_factory.construct<nm_e::communication,NO_VERTEX_ID,float,2,0>(sz,box,bc);

	// Processor graph
	Graph_CSR<nm_part_v,nm_part_e> gp = g_factory_part.construct<NO_EDGE,NO_VERTEX_ID,float,2>(sz,box,bc);

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

	// gp,p_id,loc_box,box_nn_processor,bc
	d_o.optimize<nm_v::sub_id,nm_v::id>(g,-1,dec_o,box_nn_processor,bc);

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

BOOST_AUTO_TEST_CASE( dec_optimizer_disconnected_subdomains_np)
{
	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Test for only 3 processors
	if (vcl.getProcessingUnits() > 3)
		return;

	CartesianGraphFactory<2,Graph_CSR<nm_v,nm_e>> g_factory;

	// Cartesian grid
	size_t sz[2] = {GS_SIZE,GS_SIZE};

	// Box
	Box<2,float> box({0.0,0.0},{1.0,1.0});

	// Boundary conditions, non periodic
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC};

	// Graph to decompose
	Graph_CSR<nm_v,nm_e> g = g_factory.construct<nm_e::communication,NO_VERTEX_ID,float,1,0>(sz,box,bc);

	SimpleRNG rng;

	auto vit = g.getVertexIterator();

	while (vit.isNext())
	{
		auto vk = vit.get();

		g.template vertex_p<nm_v::proc_id>(vk) = rng.GetUniform() * 2.9999;
		g.template vertex_p<nm_v::sub_id>(vk) = 100;

		++vit;
	}

	// optimize
	dec_optimizer<2,Graph_CSR<nm_v,nm_e>> d_o(g,sz);

	// set of Boxes produced by the decomposition optimizer
	openfpm::vector<::Box<2, size_t>> loc_box;

	//! for each sub-domain, contain the list of the neighborhood processors
	openfpm::vector<openfpm::vector<long unsigned int> > box_nn_processor;

	d_o.optimize<nm_v::sub_id, nm_v::proc_id>(g, vcl.getProcessUnitID(), loc_box, box_nn_processor,bc);

	std::stringstream str_g;
	str_g << "dec_optimizer_disc_graph" << vcl.getProcessUnitID() << ".vtk";
	std::stringstream str_gt;
	str_gt << "src/Decomposition/Distribution/test_data/dec_optimizer_disc_graph" << vcl.getProcessUnitID() << "_test.vtk";

	std::stringstream str_s;
	str_s << "dec_optimizer_disc_sub" << vcl.getProcessUnitID() << ".vtk";
	std::stringstream str_st;
	str_st << "src/Decomposition/Distribution/test_data/dec_optimizer_disc_sub" << vcl.getProcessUnitID() << "_test.vtk";

	VTKWriter<Graph_CSR<nm_v,nm_e>,VTK_GRAPH> wrt(g);
	wrt.write("dec_optimizer_disc_graph" + std::to_string(vcl.getProcessUnitID()) + ".vtk");

	VTKWriter<openfpm::vector<::Box<2, size_t>>, VECTOR_BOX> vtk_box1;
	vtk_box1.add(loc_box);
	vtk_box1.write("dec_optimizer_disc_sub" + std::to_string(vcl.getProcessUnitID()) + std::string(".vtk"));

	bool test = compare(str_g.str(), str_gt.str());
	BOOST_REQUIRE_EQUAL(true,test);

	test = compare(str_s.str(),str_st.str());
	BOOST_REQUIRE_EQUAL(true,test);
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* DEC_OPTIMIZE_HPP_ */
