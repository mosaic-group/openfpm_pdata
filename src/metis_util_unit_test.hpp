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

struct nm_v
{
	//! The node contain 3 unsigned long integer for communication computation memory and id
	typedef boost::fusion::vector<float,float,float,size_t,size_t,size_t,size_t,size_t> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! type of the positional field
	typedef float s_type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! computation property id in boost::fusion::vector
	static const unsigned int x = 0;
	//! computation property id in boost::fusion::vector
	static const unsigned int y = 1;
	//! memory property id in boost::fusion::vector
	static const unsigned int z = 2;
	//! computation property id in boost::fusion::vector
	static const unsigned int communication = 3;
	//! computation property id in boost::fusion::vector
	static const unsigned int computation = 4;
	//! memory property id in boost::fusion::vector
	static const unsigned int memory = 5;
	//! memory property id in boost::fusion::vector
	static const unsigned int id = 6;
	//! memory property sub_id in boost::fusion::vector
	static const unsigned int sub_id = 7;

	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 8;
};

const std::string nm_v::attributes::name[] = {"x","y","z","communication","computation","memory","id","sub_id"};

/*!
 *
 * Test node
 *
 */

struct nm_e
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<size_t> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! computation property id in boost::fusion::vector
	static const unsigned int communication = 0;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 1;
};

const std::string nm_e::attributes::name[] = {"communication"};

struct nm_part_v
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<size_t> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	typedef float s_type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data

	type data;

	//! partition id in the boost::fusion::vector
	static const unsigned int id = 0;

	//! total number of properties
	static const unsigned int max_prop = 1;
};

const std::string nm_part_v::attributes::name[] = {"id"};

struct nm_part_e
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! The data

	type data;

	//! total number of properties
	static const unsigned int max_prop = 0;
};

BOOST_AUTO_TEST_SUITE( Metis_test )

BOOST_AUTO_TEST_CASE( Metis_test_use)
{
	CartesianGraphFactory<3,Graph_CSR<nm_v,nm_e>> g_factory;
	CartesianGraphFactory<3,Graph_CSR<nm_part_v,nm_part_e>> g_factory_part;

	// Cartesian grid
	std::vector<size_t> sz;
	sz.push_back(GS_SIZE);
	sz.push_back(GS_SIZE);
	sz.push_back(GS_SIZE);

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
